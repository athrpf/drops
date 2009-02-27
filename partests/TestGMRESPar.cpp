//**************************************************************************
// File:    TestGMRESPar.cpp                                             *
// Content: parallel solver for Poisson problems                           *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
// test different setting and solvers for different multigrids for         *
// parallel solving of Poisson problems                                    *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestGMRESPar.cpp
/// \brief Testing parallel solving of the Poisson problem

 // include parallel computing!
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/exchange.h"          // transfer of numerical datas
#include <ddd.h>

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
#include "poisson/instatpoisson.h"      // setting up the Poisson problem

 // include misc
#include "out/output.h"                 // output in geomview or other formats
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
//         std::cerr << " - Writing multigrid into: " << filename<< " (for proc 0)"<<std::endl;
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
        cerr << "    --> Schreibe in Datei: " << filename << "  \n";
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
    bool sane= DROPS::ProcCL::Check(pmg_sane && mg_sane);
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
void Solve(const Mat &A, Vec &x, const Vec &b, IdxDescCL &idx)
{
    ParTimerCL time;
    double duration[3];

    // Create preconditioner class
    typedef ParDummyPcCL PreCondT;
    PreCondT precond(idx);

    // Create solver class
    bool relative=false, accure=true, mod4par=false, modGS=true;
    std::ostream* output=0;
    typedef ParPreGMResSolverCL<PreCondT> SolverT;
    SolverT solver_v1(C.restart, C.iter, C.tol, idx, precond, relative, accure, modGS, RightPreconditioning, mod4par, output);
    mod4par=true, modGS=true;
    SolverT solver_v2(C.restart, C.iter, C.tol, idx, precond, relative, accure, modGS, RightPreconditioning, mod4par, output);
    mod4par=true, modGS=false;
    SolverT solver_v3(C.restart, C.iter, C.tol, idx, precond, relative, accure, modGS, RightPreconditioning, mod4par, output);

    x=0.;
    time.Reset();
    solver_v1.Solve(A, x, b);
    time.Stop();
    duration[0]= time.GetMaxTime();

    x=0.;
    time.Reset();
    solver_v2.Solve(A, x, b);
    time.Stop(); duration[1]= time.GetMaxTime();

    x=0.;
    time.Reset();
    solver_v3.Solve(A, x, b);
    time.Stop(); duration[2]= time.GetMaxTime();

    if (ProcCL::IamMaster())
        std::cerr << "LES solved with three different versions of GMRES:\n"
                  << "           VERSION                         : ITER\t"<<std::setw(11)<<"RESID\t"<<"TIME\n"
                  << " - GMRES Version 1 (no mod for par, mod GS): "
                  << std::setw(5)<<solver_v1.GetIter()<<'\t'<<solver_v1.GetResid()<<'\t'<<duration[0]<<'\n'
                  << " - GMRES Version 2 (mod for par,  mod GS)  : "
                  << std::setw(5)<<solver_v2.GetIter()<<'\t'<<solver_v2.GetResid()<<'\t'<<duration[1]<<'\n'
                  << " - GMRES Version 3 (mod for par,  std GS)  : "
                  << std::setw(5)<<solver_v3.GetIter()<<'\t'<<solver_v3.GetResid()<<'\t'<<duration[2]
                  << std::endl;
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
void Strategy(InstatPoissonP1CL<PoissonCoeffCL>& Poisson)
{
    ParTimerCL time; double duration;

    typedef InstatPoissonP1CL<PoissonCoeffCL> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();

    MLIdxDescCL* idx= &Poisson.idx;
    VecDescCL* x= &Poisson.x; VecDescCL* b= &Poisson.b;
    VecDescCL vU, vA, vM, vf;
    MLMatDescCL* A= &Poisson.A; MLMatDescCL* M= &Poisson.M;
    MLMatDescCL* U= &Poisson.U;
    MLMatrixCL AU, AUM;

    // Point 1
    //--------
    idx->SetFE( P1_FE); // 1 unknown on vertex

    // Create numbering of the unknowns according to the index. Also the
    // ExchangeCL within the PoissonCL is initialised
    time.Reset();
    Poisson.CreateNumbering(MG.GetLastLevel(), idx);
    time.Stop(); Times.AddTime(T_CreateNumbering, time.GetMaxTime());

    // tell vectors and matrices about the unknowns
    b->SetIdx( idx); x->SetIdx( idx); vU.SetIdx( idx);
    vA.SetIdx( idx); vM.SetIdx( idx); vf.SetIdx( idx);
    A->SetIdx(idx, idx); M->SetIdx(idx, idx); U->SetIdx(idx, idx);

    Ulint *Sizes = new Ulint[ProcCL::Size()];
    ProcCL::Gather(x->Data.size(), (size_t*)Sizes,0);
    int Gsize=0;
    for (int i=0; i<ProcCL::Size(); ++i) Gsize+=Sizes[i];

    IdxT numUnknowns=idx->GetGlobalNumUnknowns(MG);
    if (ProcCL::IamMaster()){
        std::cerr << line << std::endl;
        std::cerr << " Number of unknowns (accumulated): " << Gsize << std::endl;
        std::cerr << " Number of real unknowns:          " << numUnknowns << std::endl;
        for (int i=0; i<ProcCL::Size(); ++i) std::cerr << " - Proc " <<i<< ": " <<Sizes[i]<<" unknowns\n";
    }

    // Point 2
    //--------
    if (ProcCL::IamMaster()){
        std::cerr << line << std::endl;
        std::cerr << " Discretize ... " << std::endl;
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
        std::cerr << line << std::endl;
    Solve(AUM, x->Data, b->Data, idx->GetFinest());

    // Point 4
    //--------
    if (ProcCL::IamMaster()){
        std::cerr << line << std::endl;
        std::cerr << " Check the solution ... " << std::endl;
    }

    time.Reset();
    Poisson.CheckSolution(*x, Lsg,0.);
    time.Stop(); duration = time.GetMaxTime(); Times.AddTime(T_SolCheck,duration);
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
            std::cerr << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param && ProcCL::IamMaster()){
            std::cerr << "error while opening parameter file\n"; return 1;
        }
        param >> C; param.close();
        if (ProcCL::IamMaster())
            std::cerr << C << std::endl;

        DROPS::ParTimerCL time, alltime;

        typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL> PoissonOnBCL;
        typedef PoissonOnBCL                             MyPoissonCL;

        // Init of the parallel structurs.
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

        DROPS::Point3DCL orig(0.);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= e2[1]= e3[2]= 1.0;

        // just Dirichlet boundarys
        const bool IsNeumann[6]=
            {false, false, false, false, false, false};
        const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
            { &Null, &Null, &Null, &Null, &Null, &Null};
        DROPS::InstatPoissonBndDataCL bdata(6, IsNeumann, bnd_fun);

        if (ProcCL::IamMaster()){
            std::cerr << line << std::endl
                 << " Create initial grid and distribution ... \n";
        }

        DROPS::MGBuilderCL * mgb;
        if (ProcCL::IamMaster())
            mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
        else
            mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);

        // Setup the problem
        PoissonOnBCL prob(*mgb, PoissonCoeffCL(), bdata);
        DROPS::MultiGridCL &mg = prob.GetMG();
        pmg.AttachTo(mg);

        // Init the LoadBalHandler (create an MultiGridCL and distribute the multigrid)
        DROPS::LoadBalHandlerCL lb(mg);
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
            if (ProcCL::IamMaster()) cerr << " Refine (" << ref << ") ";
            if (ref<C.refall)
            {
                if (ProcCL::IamMaster()) cerr << "regular ...\n";
                DROPS::MarkAll(mg);
            }
            else if (ref<C.refall+C.markdrop)
            {
                if (ProcCL::IamMaster()) cerr << "drop ...\n";
                MarkDrop(mg,mg.GetLastLevel());
            }
            else
            {
                if (ProcCL::IamMaster()) cerr << "corner ...\n";
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
            if (ProcCL::IamMaster()) cerr << " Check the parallel multigrid ... " << flush;
            time.Reset();
            if (CheckParMultiGrid(pmg)){
                if (ProcCL::IamMaster()) cerr << " OK\n";
            }
            else
                if (ProcCL::IamMaster()) cerr << "  !!! not OK !!!\n";
            time.Stop();
            Times.AddTime(T_Check, time.GetMaxTime());
        }
        if (ProcCL::IamMaster())
            cerr << line << endl << " Distribution of elemente:\n";
        mg.SizeInfo(cerr);

        if (C.printGEO)
            PrintGEO(pmg);

        if (C.printMG)
            PrintMG(pmg,DROPS::REF);

        DROPS::Strategy(prob);

        if (ProcCL::IamMaster()){
            cerr << line << std::endl;
        }

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        Times.Print(cerr);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
