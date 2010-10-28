/// \file TestPoissonPar.cpp
/// \brief Testing parallel solving of the Poisson problem with a large variety of linear solvers
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

 // include parallel computing!
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/exchange.h"          // transfer of numerical datas

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construuct the initial multigrid
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid

 // include numeric computing!
#include "num/parsolver.h"              // various parallel solvers
#include "num/parprecond.h"             // various parallel preconditioners
#include "num/parlanczos.h"
#include "num/fe.h"

 // include problem class
#include "poisson/poisson.h"      // setting up the Poisson problem

 // include misc
#include "out/output.h"                 // output in geomview or other formats
#include "out/ensightOut.h"             // output in ensight format
#include "partests/params.h"            // reading parameter-files

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;
/****************************************************************************
* G L O B A L  V A R I B L E S                                              *
****************************************************************************/
DROPS::TimeStoreCL Times(10);   // time table all times are stored within this class
enum TimePart{                  // numbers for accesing the time table
    T_Refine,
    T_HandleUnkAfterRef,
    T_Migration,
    T_Interpolate,
    T_HandleNewIdx,
    T_Discretize,
    T_CreateNumbering,
    T_Solve,
    T_Check,
    T_SolCheck
};

/****************************************************************************
* S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_Refine, "Refine");
    Times.SetDescriber(T_HandleUnkAfterRef, "Handle unknowns after refine");
    Times.SetDescriber(T_Migration, "Migration");
    Times.SetDescriber(T_HandleNewIdx, "Handle new index");
    Times.SetDescriber(T_Discretize, "Discretize");
    Times.SetDescriber(T_CreateNumbering, "Create Numbering");
    Times.SetDescriber(T_Solve, "Solve");
    Times.SetDescriber(T_Check, "Check MG");
    Times.SetDescriber(T_SolCheck, "Check Solution");
    Times.SetDescriber(T_Interpolate, "Interpolate to new MG");
    Times.SetCounterDescriber("Moved MultiNodes");
}

DROPS::ParamParPoissonCL C;     // parameter class
const char line[] ="----------------------------------------------------------------------------------";
using DROPS::ProcCL;

/****************************************************************************
* P R I N T  M U L T I G R I D  I N  G E O M V I E W  F O R M A T           *
****************************************************************************/
void PrintGEO(const DROPS::ParMultiGridCL& pmg)
{
    const int me=ProcCL::MyRank();
    static int num=0;
    char filename[30];
    sprintf (filename, "output/geo_%i_GEOM_%i.geo",me,num++);
    ofstream file(filename);
    file << DROPS::GeomMGOutCL(pmg.GetMG(),-1,false,0.25,0.25) << std::endl;
    file.close();
}
/****************************************************************************
* P R I N T  M U L T I G R I D  I N  A S C I I  F I L E                     *
****************************************************************************/
// void PrintMG(const DROPS::ParMultiGridCL& pmg, int type=MIG)
// {
//     const int me=DROPS::ProcCL::MyRank();
//     static int REFnum=0;
//     static int MIGnum=0;
//     char filename[30];
//
//     if (type==REF)
//         sprintf(filename, "output/%i_MG_REF_%i.mg",me,REFnum++);
//     else if (type == MIG)
//         sprintf(filename, "output/%i_MG_MIG_%i.mg",me,MIGnum++);
//     if (me==0)
//         std::cout << " - Writing multigrid into: " << filename<< " (for proc 0)"<<std::endl;
//
//     ofstream file(filename);
//     pmg.DebugInfo(file);
//     file.close();
// }

void PrintExchange(const DROPS::ExchangeCL& ex)
{
    const int me =ProcCL::MyRank();
    static int ExNum=0;
    char filename[30];
    sprintf(filename, "output/%i_Exchange_%i.mg",me,ExNum++);
    ofstream file(filename);
    ex.DebugInfo(file);
    file.close();
    if (me==0)
        cout << "    --> Schreibe in Datei: " << filename << "  \n";
}

/****************************************************************************
* M A R K I N G  P R O C E D U R E S                                        *
*****************************************************************************
*  Mark special tetras for refinement or removement                         *
****************************************************************************/
void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// Mark around droplet
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void MarkLeftDownCorner (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// mark around corner with coordinates (0,0,0)
{
    DROPS::Point3DCL Corner; Corner[0]=0.; Corner[1]=0.; Corner[2]=0.;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Corner).norm()<=0.3)
            It->SetRegRefMark();
    }
}

void MarkRightUpCorner (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// mark around corner with coordinates (1,1,1)
{
    DROPS::Point3DCL Corner; Corner[0]=1.; Corner[1]=1.; Corner[2]=1.;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Corner).norm()<=0.3)
            It->SetRegRefMark();
    }
}

void UnMarkRightUpCorner (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// unmark around corner with coordinates (1,1,1)
{
    DROPS::Point3DCL Corner; Corner[0]=1.; Corner[1]=1.; Corner[2]=1.;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Corner).norm()<=0.3)
            It->SetRemoveMark();
    }
}

void UnMarkForGhostKill (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// search for a ghost tetra and unmark all children
{
    for (DROPS::MultiGridCL::const_TetraIterator  It(mg.GetTetrasBegin(maxLevel-1)),
         ItEnd(mg.GetTetrasEnd(maxLevel-1)); It!=ItEnd; ++It)
    {
        if (It->IsGhost() && It->IsRegularlyRef()){
            for (DROPS::TetraCL::const_ChildPIterator ch(It->GetChildBegin()),
                 chEnd(It->GetChildEnd()); ch!=chEnd; ++ch)
                (*ch)->SetRemoveMark();
            return;
        }
    }
}

bool MarkIrregular (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// mark all irregular refined tetras
{
    bool newmarks=false;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if (!It->IsRegular()){
            It->SetRegRefMark();
            newmarks=true;
        }
    }
    return newmarks;
}

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
****************************************************************************/
bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    ofstream check(dat);

    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);

    check.close();
    bool sane=DROPS::ProcCL::Check(pmg_sane && mg_sane);
    if (!sane)
        throw DROPS::DROPSErrCL("CheckParMultiGrid: Error within multigrid!");
    return sane;
}

/****************************************************************************
* P O I S S O N  C O E F F  C L                                             *
*****************************************************************************
*   Coeffs of the Poisson problem:                                          *
*       - Laplace u + Vel.(grad u) +q*u = f                                 *
*       with q=0 and Vel=0                                                  *
****************************************************************************/
class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL, double){
        DROPS::Point3DCL ret; //ret[0]=1.; ret[1]=1.; ret[2]=1.;
        return ret;
    }
    static double f(const DROPS::Point3DCL& p, double= 0.0){
        return 128*(p[1]*p[2]*(1.-p[1])*(1.-p[2])
                + p[0]*p[2]*(1.-p[0])*(1.-p[2])
                + p[0]*p[1]*(1.-p[0])*(1.-p[1]));
    }
    static double alpha(const DROPS::Point3DCL&, double=0.0){
        return 1.;
    }
};

/****************************************************************************
* S O L U T I O N                                                           *
*****************************************************************************
*   Solution of the poisson problem to the above coeffs with Dirichlet      *
*   boundary (zero)                                                         *
****************************************************************************/
inline double Lsg( const DROPS::Point3DCL& p, double=0.0){
    return 64*p[0]*p[1]*p[2]*(1.-p[0])*(1.-p[1])*(1.-p[2]);
}

// boundary function
double Null (const DROPS::Point3DCL&, double) {return 0.;}

namespace DROPS
{
/****************************************************************************
* S O L V E                                                                 *
*****************************************************************************
*   Embrace all functions solve a linear system. This function constructs   *
*   all solvers, that can be used, but uses only the wanted one             *
****************************************************************************/
template <typename Mat, typename Vec>
void Solve(const Mat &A, Vec &x, const Vec &b, const IdxDescCL& idx)
{
    ParTimerCL time; double duration;

    // Definition and initialisation of preconditioner
    ParDummyPcCL PCDummy(idx);
    ParJac0CL    PCJac(idx, C.relax);
    ParSSOR0CL   PCSSOR(idx, C.relax);

    // BiCGStab as preconditioner
    typedef SolverAsPreCL< ParBiCGSTABSolverCL<ParDummyPcCL> > PCBiCGStabT;
    ParBiCGSTABSolverCL<ParDummyPcCL> PC_BiCGStab(C.pciter, C.pctol, idx, PCDummy, /*relative*/true);
    PCBiCGStabT PCBiCGStab(PC_BiCGStab);

    // Definition of solvers
    // CG-Type Solvers
    ParCGSolverCL                CGSolver(C.iter, C.tol, idx);
    ParPCGSolverCL<ParDummyPcCL> P_D_CGSolver(C.iter, C.tol, idx, PCDummy, C.relative, C.accur);
    ParPCGSolverCL<ParJac0CL>    P_J_CGSolver(C.iter, C.tol, idx, PCJac, C.relative, C.accur);
    ParPCGSolverCL<ParSSOR0CL>   P_S_CGSolver(C.iter, C.tol, idx, PCSSOR, C.relative, C.accur);
    ParPCGSolverCL<PCBiCGStabT>  P_B_CGSolver(C.iter, C.tol, idx, PCBiCGStab, C.relative, C.accur);

    // GMRES-Type Solvers
    PreMethGMRES method=(C.preCondMeth==0 ? LeftPreconditioning : RightPreconditioning);
    ParPreGMResSolverCL<ParDummyPcCL>
            P_D_GMRESSolver(C.restart, C.iter, C.tol, idx, PCDummy, C.relative, C.accur, C.useMGS, method, C.modified);
    ParPreGMResSolverCL<ParJac0CL>
            P_J_GMRESSolver(C.restart, C.iter, C.tol, idx, PCJac, C.relative, C.accur, C.useMGS, method, C.modified);
    ParPreGMResSolverCL<ParSSOR0CL>
            P_S_GMRESSolver(C.restart, C.iter, C.tol, idx, PCSSOR, C.relative, C.accur, C.useMGS, method, C.modified);
    ParPreGMResSolverCL<PCBiCGStabT>
            P_B_GMRESSolver(C.restart, C.iter, C.tol, idx, PCBiCGStab, C.relative, C.accur, C.useMGS, method, C.modified);


    // GCR-Type Solvers
    ParPreGCRSolverCL<ParDummyPcCL>
            P_D_GCRSolver(C.restart,C.iter,C.tol, idx, PCDummy, C.modified, C.relative, C.accur);
    ParPreGCRSolverCL<ParJac0CL>
            P_J_GCRSolver(C.restart,C.iter,C.tol, idx, PCJac, C.modified, C.relative, C.accur);
    ParPreGCRSolverCL<ParSSOR0CL>
            P_S_GCRSolver(C.restart,C.iter,C.tol, idx, PCSSOR, C.modified, C.relative, C.accur);
    ParPreGCRSolverCL<PCBiCGStabT>
            P_B_GCRSolver(C.restart,C.iter,C.tol, idx, PCBiCGStab, C.modified, C.relative, C.accur);

    // Lanczos-Algorithms
    typedef ParLanczos2CL<MatrixCL,VectorCL,ExchangeCL> LanczosT;
    LanczosT Lan(1e-32);
    typedef ParPreLanczos2CL<MatrixCL,VectorCL,ParDummyPcCL,ExchangeCL> P_D_LanczosT;
    P_D_LanczosT LanDummy(PCDummy,1e-32);
    typedef ParPreLanczos2CL<MatrixCL,VectorCL,ParJac0CL,ExchangeCL> P_J_LanczosT;
    P_J_LanczosT LanJacobi(PCJac,1e-32);

    // QMR-Method
    ParQMRSolverCL<LanczosT>           ParQMR(C.iter,  C.tol, idx, Lan);
    ParQMRSolverCL<P_D_LanczosT>       P_D_QMR(C.iter, C.tol, idx, LanDummy);
    ParQMRSolverCL<P_J_LanczosT>       P_J_QMR(C.iter, C.tol, idx, LanJacobi);

    // BiCGStab-Solver
    ParBiCGSTABSolverCL<ParDummyPcCL> P_D_BiCGStab(C.iter, C.tol, idx, PCDummy, C.relative);
    ParBiCGSTABSolverCL<ParJac0CL>    P_J_BiCGStab(C.iter, C.tol, idx, PCJac, C.relative);
    ParBiCGSTABSolverCL<ParSSOR0CL>   P_S_BiCGStab(C.iter, C.tol, idx, PCSSOR, C.relative);
    ParBiCGSTABSolverCL<PCBiCGStabT>  P_B_BiCGStab(C.iter, C.tol, idx, PCBiCGStab, C.relative);


    if (ProcCL::IamMaster()){
        std::cout <<" Solve the linear system with ";
        switch(C.solver){
            case  0: std::cout<<"CG method: iter=" <<C.iter<< ", tol="<<C.tol<<std::endl; break;
            case  1: std::cout<<"P-CG method: iter=" <<C.iter<< ", tol="<<C.tol; break;
            case  2: std::cout<<"P-GMRES method: iter=" <<C.iter<< ", tol="<<C.tol<<", restart="<<C.restart; break;
            case  3: std::cout<<"P-GCR method: iter=" <<C.iter<< ", tol="<<C.tol<<", truncate="<<C.restart; break;
            case  4: std::cout<<"QMR method: iter=" <<C.iter<< ", tol="<<C.tol<<", BreakDownEps="<<Lan.GetBreakDownEps()<<std::endl; break;
            case  5: std::cout<<"P-QMR method: iter=" <<C.iter<< ", tol="<<C.tol<<", BreakDownEps="<<Lan.GetBreakDownEps(); break;
            case  6: std::cout <<"P-BiCGStab: iter=" <<C.iter<< ", tol="<<C.tol; break;
            default: std::cout<<" \"I do not know this method, breaking\"\n"; return;
        }
        if ((C.solver!=0 && C.solver!=4) && C.precond==0)
            std::cout<<", PreCond=Dummy,\n";
        if ((C.solver!=0 && C.solver!=4) && C.precond==1)
            std::cout<<", PreCond=Jac("<<C.relax<<"),\n";
        if ((C.solver!=0 && C.solver!=4) && C.precond==2)
            std::cout<<", inexact PreCond=SSOR0("<<C.relax<<"),\n";
        if ((C.solver!=0 && C.solver!=4) && C.precond==4)
            std::cout<<", BiCGStab as Preconditioner(iter "<<C.pciter<<", tol "<<C.pctol<<", relative true),\n";
        std::cout << " using "<<(C.accur ? "accur " : "not accur ")<<"version for inner products\n";
        if (C.solver!=0 && C.solver!=4)
            std::cout << " using "<<(C.relative?"relative":"no relative")<<" tolerance meassurement"<<std::endl;
        if (C.modified)
            std::cout << " using modificated version"<< std::endl;
        else
            std::cout << " using serial version"<< std::endl;

        if (C.useMGS==0)
            std::cout << " using standard Gramm-Schmidt for orthogonalization"<<std::endl;
        else
            std::cout << " using modified Gramm-Schmidt for orthogonalization"<<std::endl;

        if (C.solver==2)
            std::cout << " using "<<(C.preCondMeth==0 ? "left " : "right ")<<"preconditioning\n";

    }
    if (C.solver==5 && (C.precond==2 || C.precond==3)){
        if (ProcCL::IamMaster())
            std::cout << ">>> No SSOR or BiCGStab for Lanczos as preconditioner impemented, abort solving system ..."<<std::endl;
        return;
    }

    time.Reset();
    int steps; double res;
    switch(C.solver*10+C.precond){      // Choose solver and try solving the system
        case  0: case  1: case  2: case  3:
            CGSolver.Solve(A, x, b); steps=CGSolver.GetIter(); res=CGSolver.GetResid(); break;

        case 10: P_D_CGSolver.Solve(A, x, b);    steps=P_D_CGSolver.GetIter(); res=P_D_CGSolver.GetResid(); break;
        case 11: P_J_CGSolver.Solve(A, x, b);    steps=P_J_CGSolver.GetIter(); res=P_J_CGSolver.GetResid(); break;
        case 12: P_S_CGSolver.Solve(A, x, b);    steps=P_S_CGSolver.GetIter(); res=P_S_CGSolver.GetResid(); break;
        case 13: P_B_CGSolver.Solve(A, x, b);    steps=P_S_CGSolver.GetIter(); res=P_S_CGSolver.GetResid(); break;

        case 20: P_D_GMRESSolver.Solve(A, x, b); steps=P_D_GMRESSolver.GetIter(); res=P_D_GMRESSolver.GetResid(); break;
        case 21: P_J_GMRESSolver.Solve(A, x, b); steps=P_J_GMRESSolver.GetIter(); res=P_J_GMRESSolver.GetResid(); break;
        case 22: P_S_GMRESSolver.Solve(A, x, b); steps=P_S_GMRESSolver.GetIter(); res=P_S_GMRESSolver.GetResid(); break;
        case 23: P_B_GMRESSolver.Solve(A, x, b); steps=P_S_GMRESSolver.GetIter(); res=P_S_GMRESSolver.GetResid(); break;

        case 30: P_D_GCRSolver.Solve(A, x, b);   steps=P_D_GCRSolver.GetIter(); res=P_D_GCRSolver.GetResid(); break;
        case 31: P_J_GCRSolver.Solve(A, x, b);   steps=P_J_GCRSolver.GetIter(); res=P_J_GCRSolver.GetResid(); break;
        case 32: P_S_GCRSolver.Solve(A, x, b);   steps=P_S_GCRSolver.GetIter(); res=P_S_GCRSolver.GetResid(); break;
        case 33: P_B_GCRSolver.Solve(A, x, b);   steps=P_S_GCRSolver.GetIter(); res=P_S_GCRSolver.GetResid(); break;

        case 40: case 41: case 42: case 43:
            ParQMR.Solve(A,x,b);            steps=ParQMR.GetIter(); res=ParQMR.GetResid(); break;

        case 50: P_D_QMR.Solve(A,x,b);           steps=P_D_QMR.GetIter(); res=P_D_QMR.GetResid(); break;
        case 51: P_J_QMR.Solve(A,x,b);           steps=P_J_QMR.GetIter(); res=P_J_QMR.GetResid(); break;

        case 60: P_D_BiCGStab.Solve(A,x,b);      steps=P_D_BiCGStab.GetIter(); res=P_D_BiCGStab.GetResid(); break;
        case 61: P_J_BiCGStab.Solve(A,x,b);      steps=P_J_BiCGStab.GetIter(); res=P_J_BiCGStab.GetResid(); break;
        case 62: P_S_BiCGStab.Solve(A,x,b);      steps=P_S_BiCGStab.GetIter(); res=P_S_BiCGStab.GetResid(); break;
        case 63: P_B_BiCGStab.Solve(A,x,b);      steps=P_S_BiCGStab.GetIter(); res=P_S_BiCGStab.GetResid(); break;
        default: std::cout << "No solver found\n"; return;
    }
    time.Stop(); duration = time.GetMaxTime(); Times.AddTime(T_Solve, duration);
    double realresid = idx.GetEx().Norm(VectorCL(A*x-b),false);
    double normb     = idx.GetEx().Norm(b,false);
    if (ProcCL::IamMaster()){
        std::cout << " - Iterations:  " << steps << std::endl;
        std::cout << " - Residuum:    " << res  << std::endl;
        std::cout << " - real Resid:  " << realresid << std::endl;
        if (C.poisson_printTime)
            std::cout << " - Time:        " << duration << std::endl;
        std::cout << " - Norm(b):     " << normb << std::endl;
    }
}


/****************************************************************************
* S T R A T E G Y                                                           *
*****************************************************************************
* The Poisson problem is solved by this non adaptive procedure. These are   *
* the main steps:                                                           *
*   1. Numbering the unknwons (map unknowns to geometry) and create         *
*      ExchangeCL                                                           *
*   2. Discretize                                                           *
*   3. Solve the linear equation system                                     *
*   4. Check the solution                                                   *
****************************************************************************/
template<class PoissonCoeffCL>
void Strategy(PoissonP1CL<PoissonCoeffCL>& Poisson)
{
    ParTimerCL time; double duration;

    typedef PoissonP1CL<PoissonCoeffCL> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();

    MLIdxDescCL* idx= &Poisson.idx;
    VecDescCL* x= &Poisson.x; VecDescCL* b= &Poisson.b;
    VecDescCL vU, vA, vM, vf;
    MLMatDescCL* A= &Poisson.A; MLMatDescCL* M= &Poisson.M;
    MLMatDescCL* U= &Poisson.U;
    MLMatrixCL AU, AUM;

    // Point 1
    //--------
    idx->SetFE( P1_FE); // linear finite elements ...

    // Create numbering of the unknowns according to the index. Also the
    // ExchangeCL within the PoissonCL is initialised
    time.Reset();
    Poisson.CreateNumbering(MG.GetLastLevel(), idx);
    time.Stop(); Times.AddTime(T_CreateNumbering, time.GetMaxTime());

    // tell vectors and matrices about the unknowns
    b->SetIdx( idx); x->SetIdx( idx); vU.SetIdx( idx);
    vA.SetIdx( idx); vM.SetIdx( idx); vf.SetIdx( idx);
    A->SetIdx(idx, idx); M->SetIdx(idx, idx); U->SetIdx(idx, idx);

    if (ProcCL::IamMaster() && C.ensight){
        std::cout << line << std::endl;
        std::cout << " Write ensight case and geometry files ... " << std::endl;
    }

    Ensight6OutCL *ensight= 0;
    if (C.ensight)
    {
        // Create ensight case and geometry file
        std::string ensf( C.EnsDir + "/" + C.EnsCase);
        ensight = new Ensight6OutCL( C.EnsCase + ".case", 0, false, true);
        ensight->Register( make_Ensight6Geom      ( MG, idx->GetFinest().TriangLevel(), C.geomName, ensf + ".geo"));
        ensight->Register( make_Ensight6Scalar    ( Poisson.GetSolution(),              C.varName,  ensf + ".tmp"));
        ensight->SetMasterOut();
    }

    size_t *Sizes = new size_t[ProcCL::Size()];
    ProcCL::Gather(x->Data.size(), Sizes,0);
    int Gsize=0;
    for (int i=0; i<ProcCL::Size(); ++i) Gsize+=Sizes[i];

    IdxT numUnknowns=idx->GetGlobalNumUnknowns(MG);
    if (ProcCL::IamMaster()){
        std::cout << " Number of unknowns (accumulated): " << Gsize << std::endl;
        std::cout << " Number of real unknowns:          " << numUnknowns << std::endl;
        for (int i=0; i<ProcCL::Size(); ++i) std::cout << " - Proc " <<i<< ": " <<Sizes[i]<<" unknowns\n";
    }

    // Point 2
    //--------
    if (ProcCL::IamMaster()){
        std::cout << line << std::endl;
        std::cout << " Discretize ... " << std::endl;
    }

    time.Reset();
    Poisson.SetupInstatSystem(*A, *M,0.);
    Poisson.SetupConvection(*U,vU,0.);
    Poisson.SetupInstatRhs(vA,vM,0.,vf,0.);
    AU.LinComb( C.nu, A->Data, 1., U->Data);
    AUM.LinComb( 1., AU, 1., M->Data);
    b->Data = vf.Data + vU.Data + vA.Data  + vM.Data;
    time.Stop(); ::Times.AddTime(T_Discretize, time.GetMaxTime());

    // Point 3
    //--------
    // Solve
    if (ProcCL::IamMaster())
        std::cout << line << std::endl;
    Solve(AUM.GetFinest(), x->Data, b->Data, idx->GetFinest());

    // Point 4
    //--------
    if (ProcCL::IamMaster()){
        std::cout << line << std::endl;
        std::cout << " Check the solution ... " << std::endl;
    }

    time.Reset();
    Poisson.CheckSolution(*x, Lsg,0.);
    time.Stop(); duration = time.GetMaxTime(); Times.AddTime(T_SolCheck,duration);

    if (C.ensight){
        if (ProcCL::IamMaster()){
            std::cout << line << std::endl;
            std::cout << " Write ensight-Data files" << std::endl;
        }

        ensight->Write();
    }
}

/****************************************************************************
* S T R A T E G Y  A D A P T I V E                                          *
*****************************************************************************
* The Poisson problem is solved by this adaptive procedure. These are the   *
* main steps:                                                               *
*   while (there are new marks and step<max_step) do                        *
*       1. Refinement                                                       *
*       2. Handle unknowns after refine                                     *
*       3. Migration                                                        *
*       4. Creation of Numbering and ExchangeCL                             *
*       5. Resize vectors and matrices                                      *
*       6. Handle transfered unknowns                                       *
*       7. Repair and Interpolation of the P1-Function                      *
*       8. Discretization                                                   *
*       9. Solve of the linear equation system                              *
*      10. Error estimation                                                 *
*   end do                                                                  *
****************************************************************************/
/// \brief Solve the Poisson equation by an adaptive strategy
template<class PoissonCoeffCL>
void Strategy_Adaptive(PoissonP1CL<PoissonCoeffCL>& Poisson, ParMultiGridCL &pmg, LoadBalHandlerCL &lb)
{
    ParTimerCL time;                                                        // time stamping

    typedef PoissonP1CL<PoissonCoeffCL> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();                                       // MultiGrid
    const typename MyPoissonCL::BndDataCL& BndData= Poisson.GetBndData();   // Boundary-Values

    MLIdxDescCL  loc_idx;
    VecDescCL  loc_x;

    // Vector x is "adaptive"
    MLIdxDescCL* old_idx= &loc_idx;       // old index
    MLIdxDescCL* new_idx= &Poisson.idx;   // new index after load balancing
    VecDescCL* old_x= &loc_x;             // old solution
    VecDescCL* new_x= &Poisson.x;         // new solution

    VecDescCL* b= &Poisson.b;             // rhs of the system
    MLMatDescCL* A= &Poisson.A;           // coeffs of the linear system
    VecDescCL vU, vA, vM, vf;
    MLMatDescCL* M= &Poisson.M; MLMatDescCL* U= &Poisson.U;
    MLMatrixCL AU, AUM;                   // linear combinations of matrices

    old_idx->SetFE( P1_FE);               // linear finite elements ( 1 unknowns of each vertex, that does not lie on the dirichlet boundary)
    new_idx->SetFE( P1_FE);

    int step= 0;                            // number of the actual step
    bool newmarks=false;                    // set by the estimator, if it markes tetras
    std::vector<VecDescCL*> DescWrapper(1); // for ParMultiGridCL we need a vector of (VecDescCL*)

    // Output of the solution in ensight format
    Ensight6OutCL *ensight=0;

    old_x->SetIdx(old_idx);

    do{                                     // start the adaptive cycle
        if (step==8)
            std::cout << "In step "<<step<<std::endl;
        ParTimerCL steptimer; steptimer.Start();
        if (ProcCL::IamMaster()){
            std::cout << line <<"\n Step " <<step<< " of the adaptive cycle\n";
            if (step>0){
                std::cout << "   Indices:  old_idx: "<<old_idx->GetIdx()<< " ("<<old_idx->TriangLevel()<<") "<<std::endl
                          << "             new_idx: "<<new_idx->GetIdx()<< " ("<<new_idx->TriangLevel()<<") "<<std::endl;
            }
            std::cout << line << std::endl;
        }

        // Point 1: Refinement
        // -----------------------------------------------------------------
        if (ProcCL::IamMaster())
            std::cout << " - Refine ... " << std::endl;

        pmg.AttachTo(old_x, &BndData);  // unknowns and boundary conditions

        time.Reset();
        pmg.Refine();
        time.Stop(); Times.AddTime(T_Refine, time.GetMaxTime());
        PrintMG(pmg, REF);



        // Point 2: Handle unknowns after refine
        // -----------------------------------------------------------------
        if (ProcCL::IamMaster())
            std::cout << " - Refine ... " << std::endl;

        if (MG.KilledGhosts() && ProcCL::IamMaster())
                std::cout << "   + Refine with killed ghosts!\n";
        time.Start();
        pmg.HandleUnknownsAfterRefine();
        time.Stop(); Times.AddTime(T_HandleUnkAfterRef, time.GetMaxTime());


        // Point 3: Migration
        // -----------------------------------------------------------------
        if (ProcCL::IamMaster())
            std::cout << " - Load balancing ... \n";

        time.Reset();
        lb.DoMigration();
        time.Stop(); Times.AddTime(T_Migration, time.GetMaxTime());
        Times.IncCounter(lb.GetMovedMultiNodes());
        PrintMG(pmg, MIG);

        // Point 4: Creation of Numbering and ExchangeCL
        // -----------------------------------------------------------------
        if (ProcCL::IamMaster())
            std::cout << " - Create Numbering with ExchangeCL of Index "<<new_idx->GetIdx()<<" for Discretization and Solving ..." << std::endl;

        time.Reset();
        Poisson.CreateNumbering(MG.GetLastLevel(), new_idx);
        time.Stop(); Times.AddTime(T_CreateNumbering, time.GetMaxTime());

        // Point 5: Resize vectors and matrices
        // -----------------------------------------------------------------
        new_x->SetIdx( new_idx); b->SetIdx( new_idx);
        vU.SetIdx( new_idx);     vA.SetIdx( new_idx);
        vM.SetIdx( new_idx);     vf.SetIdx( new_idx);
        A->Reset();              A->SetIdx(new_idx, new_idx);
        M->Reset();              M->SetIdx(new_idx, new_idx);
        U->Reset();              U->SetIdx(new_idx, new_idx);
        AU.clear();              AUM.clear();

        // Point 6: Handle transfered unknowns
        // -----------------------------------------------------------------
        if (ProcCL::IamMaster())
            std::cout << " - Handle transfered unknowns because of transfering index "
                      <<new_x->RowIdx->GetIdx()<<" ... \n";

        time.Reset();
        pmg.HandleNewIdx(old_idx, new_x);
        time.Stop(); Times.AddTime(T_HandleNewIdx, time.GetMaxTime());
        ParMultiGridCL::DeleteUnksOnGhosts();

        if (C.printMG)                   // output of the parallel multi grid in ascii file
            PrintMG(pmg,MIG);

        if (C.printUnknowns)                    // output of the problem size
        {
            IdxT old_realUnknowns = old_idx->GetGlobalNumUnknowns(MG);
            IdxT new_realUnknowns = new_idx->GetGlobalNumUnknowns(MG);
            IdxT old_accUnknwons = ProcCL::GlobalSum(old_x->Data.size());
            IdxT new_accUnknwons = ProcCL::GlobalSum(new_x->Data.size());
            if (ProcCL::IamMaster())
                std::cout << " - Number of Unknowns before after balancing\n"
                          << "   + Number of accumulated unknowns: "
                          << " old level " << old_accUnknwons << ", new level "<<new_accUnknwons<<std::endl
                          << "   + Number of real unknowns:         "
                          << "old level " << old_realUnknowns<< ", new level "<< new_realUnknowns << std::endl;
        }

        if (step>0)
        {
            // Point 7: Repair and Interpolation of the P1-Function
            // -----------------------------------------------------------------
            if (ProcCL::IamMaster())
                std::cout << " - Interpolating old index ("<<old_idx->GetIdx()<<") to the new one ("<<new_idx->GetIdx()<<") ...\n";

            P1EvalCL<double, const PoissonBndDataCL, const VecDescCL>  oldDesc(new_x, &BndData, &MG);

            time.Reset();
            RepairAfterRefineP1 (oldDesc, *new_x);
            pmg.CompleteRepair(new_x);
            pmg.DelAllUnkRecv();
            pmg.DeleteRecvBuffer();
            pmg.DeleteVecDesc();
            time.Stop(); Times.AddTime(T_Interpolate, time.GetMaxTime());

            Poisson.DeleteNumbering( old_idx );
            old_x->Clear(new_x->t);

            if (C.ensight)
                ensight->Write( step+1);
        }
        else    // init ensight
        {
            if (C.ensight)
            {
                // Create ensight files
                std::string ensf( C.EnsDir + "/" + C.EnsCase);
                ensight = new Ensight6OutCL( C.EnsCase + ".case", 7, false);
                ensight->Register( make_Ensight6Geom  ( MG, new_idx->GetFinest().TriangLevel(), C.geomName, ensf + ".geo", true));
                ensight->Register( make_Ensight6Scalar( Poisson.GetSolution(),              C.varName,  ensf + ".tmp", true));
                P1EvalCL<double, const PoissonBndDataCL, const VecDescCL>  ensSol(new_x, &BndData, &MG);
                ensight->Register( make_Ensight6Scalar( ensSol, "Interpolation",  ensf + ".interpol", true));
                ensight->SetMasterOut();

                ensight->Write( step+1);
            }
        }

//         if (step==8){
//             pmg.AttachTo(0, new_x);
//             PrintMG(pmg, REF);
//         }


        if (C.check)                    // checking of the parallel multigrid and solution
        {
            if (ProcCL::IamMaster())
                std::cout << " - Checking\n   + the parallel multigrid ... " << std::flush;
            time.Reset();
            if (CheckParMultiGrid(pmg)){
                if (ProcCL::IamMaster()) std::cout << " OK\n";
            }
            else
                if (ProcCL::IamMaster()) std::cout << "  !!! not OK !!!\n";

            if (ProcCL::IamMaster())
                std::cout << "   + the interpolated solution ... " << std::endl;

            Poisson.CheckSolution(*new_x, &::Lsg, 0);
            time.Stop(); Times.AddTime(T_Check, time.GetMaxTime());
        }

        if (ProcCL::IamMaster())
            std::cout << "   + Check if solution is still accumulated ... ";

        const bool isacc= ProcCL::Check(Poisson.idx.GetEx().IsAcc(new_x->Data));
        if (isacc){
            if (ProcCL::IamMaster())
                std::cout << "OK\n";
        }
        else{
            if (ProcCL::IamMaster())
                std::cout << "not OK!!!\n";
            throw DROPSErrCL("Solution is not accumulated any more!");
        }

        if (C.printSize){
            if (ProcCL::IamMaster())
                std::cout << " - Distribution of elemente:\n";
            MG.SizeInfo(cout);
        }

        // Point 8: Discretization
        // -----------------------------------------------------------------
        if (ProcCL::IamMaster()) std::cout << " - Discretize ... \n";

        time.Reset();
        Poisson.SetupInstatSystem(*A, *M,0.);
        Poisson.SetupConvection(*U,vU,0.);
        Poisson.SetupInstatRhs(vA,vM,0.,vf,0.);
        AU.LinComb( C.nu, A->Data, 1., U->Data);
        AUM.LinComb( 1., AU, 1., M->Data);
        b->Data = vf.Data + vU.Data + vA.Data  + vM.Data;
        time.Stop(); ::Times.AddTime(T_Discretize, time.GetMaxTime());


        // Point 9: Solve of the linear equation system
        // -----------------------------------------------------------------
        time.Reset();
        Solve(A->Data.GetFinest(), new_x->Data, b->Data, Poisson.idx.GetFinest());
        time.Stop(); ::Times.AddTime(T_Solve, time.GetMaxTime());

        A->Reset();
        b->Reset();

        if (C.ensight)
            ensight->Write( step+1);

        // Point 10: Error estimation
        // -----------------------------------------------------------------
        //  Test different refinements. These are not necessarily reasonable
        //  for solving the Poisson-equation
        newmarks=false;
        if (step==0)
        {
            if (ProcCL::IamMaster()) std::cout << " - Mark all ... \n";
            MarkAll(MG);
            newmarks=true;
        }
        if (step==1)
        {
            if (ProcCL::IamMaster()) std::cout << " - Mark tetras around drop ... \n";
            MarkDrop(MG,MG.GetLastLevel());
            newmarks=true;
        }
        if (step==2)
        {
            if (ProcCL::IamMaster()) std::cout << " - Mark left down corner ... \n";
            ::MarkLeftDownCorner(MG,MG.GetLastLevel());
            newmarks=true;
        }
        if (step==3)
        {
            if (ProcCL::IamMaster()) std::cout << " - Mark right up corner ... \n";
            ::MarkRightUpCorner(MG,MG.GetLastLevel());
            newmarks=true;
        }
        if (step==4)
        {
            if (ProcCL::IamMaster()) std::cout << " - UnMark right up corner ... \n";
            ::UnMarkRightUpCorner(MG,MG.GetLastLevel());
            newmarks=true;
        }
        if (step==5)
        {
            if (ProcCL::IamMaster()) std::cout << " - UnMark for ghost kill ... \n";
            ::UnMarkForGhostKill(MG,MG.GetLastLevel());
            newmarks=true;
        }
        if (step>5)
        {
            if (step%2)
                MarkAll(MG);
            else
                UnMarkAll(MG);
            newmarks=true;
        }

        if (ProcCL::IamMaster()) std::cout << " - Check solution:\n";

        time.Reset();
        Poisson.CheckSolution(*new_x, &::Lsg, 0);
        time.Stop(); Times.AddTime(T_SolCheck, time.GetMaxTime());

        // swap old and new solution
        std::swap(old_x, new_x);
        std::swap(old_idx, new_idx);
        step++;

        steptimer.Stop();
        double duration=steptimer.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cout << "==> Step "<<step<<" took "<<duration<<" sec.\n";
    } while(newmarks && step<=C.refall);

}
} // end of namespace DROPS

int main (int argc, char** argv)
{
    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
    try
    {
        SetDescriber();
        //DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE/*|XFER_SHOW_MSGSALL*/);

        if (argc!=2 && ProcCL::IamMaster()){
            std::cout << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param && ProcCL::IamMaster()){
            std::cout << "error while opening parameter file\n"; return 1;
        }
        param >> C; param.close();
        if (ProcCL::IamMaster())
            std::cout << C << std::endl;

        DROPS::ParTimerCL time, alltime;

        typedef DROPS::PoissonP1CL<PoissonCoeffCL> PoissonOnBCL;
        typedef PoissonOnBCL                       MyPoissonCL;

        // Init of the parallel structurs.
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

        DROPS::Point3DCL orig(0.);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= e2[1]= e3[2]= 1.0;

        // just Dirichlet boundarys
        const bool IsNeumann[6]=
            {false, false, false, false, false, false};
        const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
            { &Null, &Null, &Null, &Null, &Null, &Null};
        DROPS::PoissonBndDataCL bdata(6, IsNeumann, bnd_fun);

        if (ProcCL::IamMaster()){
            std::cout << line << std::endl
                 << " Create initial grid and distribution ... \n";
        }

        DROPS::MGBuilderCL * mgb;
        if (ProcCL::IamMaster())
            mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, C.brk_BasicRefX, C.brk_BasicRefY, C.brk_BasicRefZ);
        else
            mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);

        // Setup the problem
        PoissonOnBCL prob(*mgb, PoissonCoeffCL(), bdata);
        DROPS::MultiGridCL &mg = prob.GetMG();
        pmg.AttachTo(mg);

        // Init the LoadBalHandler (create an MultiGridCL and distribute the multigrid)
        DROPS::LoadBalHandlerCL lb(mg, DROPS::metis);
        lb.DoInitDistribution(ProcCL::Master());
        lb.SetXferUnknowns(true);
        switch (C.refineStrategy){
            case 0 : lb.SetStrategy(DROPS::NoMig);     break;
            case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
            case 2 : lb.SetStrategy(DROPS::Recursive); break;
        }
        if (C.printSize) lb.SetDebugMode(true);

        int steps = ( C.adaptiv ? 0 :  C.refall + C.markdrop + C.markcorner);

        for (int ref=0; ref<steps; ++ref)
        {
            // Markieren und verfeinern
            if (ProcCL::IamMaster()) cout << " Refine (" << ref << ") ";
            if (ref<C.refall)
            {
                if (ProcCL::IamMaster()) cout << "regular ...\n";
                DROPS::MarkAll(mg);
            }
            else if (ref<C.refall+C.markdrop)
            {
                if (ProcCL::IamMaster()) cout << "drop ...\n";
                MarkDrop(mg,mg.GetLastLevel());
            }
            else
            {
                if (ProcCL::IamMaster()) cout << "corner ...\n";
                MarkLeftDownCorner(mg,mg.GetLastLevel());
            }
            time.Reset();
            pmg.Refine();
            time.Stop(); Times.AddTime(T_Refine,time.GetMaxTime());
            // Last-Balancierung berechnen und neu Verteilen
            time.Reset();
            lb.DoMigration();
            time.Stop(); Times.AddTime(T_Migration, time.GetMaxTime());
            Times.IncCounter(lb.GetMovedMultiNodes());
        }

        if (C.check)
        {
            if (ProcCL::IamMaster()) cout << " Check the parallel multigrid ... " << flush;
            time.Reset();
            if (CheckParMultiGrid(pmg)){
                if (ProcCL::IamMaster()) cout << " OK\n";
            }
            else
                if (ProcCL::IamMaster()) cout << "  !!! not OK !!!\n";
            time.Stop();
            Times.AddTime(T_Check, time.GetMaxTime());
        }
        if (C.printSize){
            if (ProcCL::IamMaster())
                cout << " Distribution of elemente:\n";
            mg.SizeInfo(cout);
        }

        if (C.printGEO)
            PrintGEO(pmg);

        if (C.printMG)
            PrintMG(pmg,DROPS::REF);

        if (!C.adaptiv)
            DROPS::Strategy(prob);
        else
            DROPS::Strategy_Adaptive(prob, pmg, lb);

        if (ProcCL::IamMaster()){
            cout << line << std::endl;
        }

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        if (C.poisson_printTime)
            Times.Print(cout);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}

