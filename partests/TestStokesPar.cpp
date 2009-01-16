//**************************************************************************
// File:    TestStokesPar.cpp                                              *
// Content: parallel solver for Stokes problems                            *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestStokesPar.cpp
/// \brief Testing parallel solvers for the stat. Stokes problem

// include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include "parallel/exchange.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "num/spmat.h"
#include "num/parsolver.h"
#include "num/parprecond.h"
#include "num/stokessolver.h"
#include "num/parstokessolver.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"

 // include problem class
#include "stokes/stokes.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#endif

using namespace std;

enum TimePart{
    T_init,
    T_ref,
    T_ex,
    T_disc,
    T_solve
};
DROPS::ParamParStokesCL C;
DROPS::TimeStoreCL Times(5);
const char line[] ="------------------------------------------------------------";
using DROPS::ProcCL;

/****************************************************************************
    * S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_init, "Initialization");
    Times.SetDescriber(T_ref, "Refinement");
    Times.SetDescriber(T_ex, "Create ExchangeCL");
    Times.SetDescriber(T_disc, "Discretize");
    Times.SetDescriber(T_solve, "Solve");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
*****************************************************************************
*   Checkt, ob die parallele Verteilung und die MultiGrid-Struktur gesund   *
*   zu sein scheint.                                                        *
****************************************************************************/
void CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    sprintf(dat,"output/sane%i.chk",ProcCL::MyRank());
    ofstream check(dat);
    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);
    check.close();
    if( DROPS::Check(pmg_sane && mg_sane) ){
        IF_MASTER
          std::cerr << " As far as I can tell, the multigrid is sane\n";
    }
    else
        throw DROPS::DROPSErrCL("Found error in multigrid!");
}

/****************************************************************************
* S O L U T I O N                                                           *
****************************************************************************/
inline DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    sin(p[0])*sin(p[1])*sin(p[2])/3.;
    ret[1]=  - cos(p[0])*cos(p[1])*sin(p[2])/3.;
    ret[2]= 2.*cos(p[0])*sin(p[1])*cos(p[2])/3.;
    return ret;
}

// Jacobi-matrix od exact solution
inline DROPS::SMatrixCL<3, 3> DLsgVel(const DROPS::Point3DCL& p)
{
    DROPS::SMatrixCL<3, 3> ret;
    ret(0,0)= cos(p[0])*sin(p[1])*sin(p[2])/3.;
    ret(0,1)= sin(p[0])*cos(p[1])*sin(p[2])/3.;
    ret(0,2)= sin(p[0])*sin(p[1])*cos(p[2])/3.;

    ret(1,0)=   sin(p[0])*cos(p[1])*sin(p[2])/3.;
    ret(1,1)=   cos(p[0])*sin(p[1])*sin(p[2])/3.;
    ret(1,2)= - cos(p[0])*cos(p[1])*cos(p[2])/3.;

    ret(2,0)= -2.*sin(p[0])*sin(p[1])*cos(p[2])/3.;
    ret(2,1)=  2.*cos(p[0])*cos(p[1])*cos(p[2])/3.;
    ret(2,2)= -2.*cos(p[0])*sin(p[1])*sin(p[2])/3.;
    return ret;
}

// Volume of the box: 0.484473073129685
// int(p)/vol = -0.125208551608365
inline double LsgPr(const DROPS::Point3DCL& p, double){
    return cos(p[0])*sin(p[1])*sin(p[2]) - 0.125208551608365;
}
inline double LsgPr(const DROPS::Point3DCL& p){
    return cos(p[0])*sin(p[1])*sin(p[2]) - 0.125208551608365;
}

class StokesCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double)
        { DROPS::SVectorCL<3> ret(0.0); ret[2]= 3.*cos(p[0])*sin(p[1])*cos(p[2]); return ret; }
    const double nu;
    StokesCoeffCL(double Nu) : nu(Nu) {}
};

typedef DROPS::StokesP2P1CL<StokesCoeffCL> StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

namespace DROPS // for Strategy
{
/****************************************************************************
* S T R A T E G Y                                                           *
*****************************************************************************
*   Das ist die Strategy, um die Poisson-Gleichung, die durch die           *
*   Koeffizienten von oben gegeben sind, zu loesen. Es werden lineare FE    *
*   und das CG-Verfahren für das Loesen des linearen Gleichungssystems      *
*   verwendet. Es wird nicht adaptiv vorgegangen                            *
****************************************************************************/
template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, const ParMultiGridCL& /*pmg*/)
{
    ParTimerCL time;
    double duration;

    MultiGridCL& MG= Stokes.GetMG();

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    VelVecDescCL* v= &Stokes.v;
    VecDescCL*    p= &Stokes.p;
    VelVecDescCL* b= &Stokes.b;
    VecDescCL*    c= &Stokes.c;

    MLMatDescCL* A= &Stokes.A;
    MLMatDescCL* B= &Stokes.B;

//     ExchangeBlockCL& Ex=Stokes.GetEx();
    const Uint pressure=Stokes.pressure, velocity=Stokes.velocity;
    ExchangeCL& ExP= Stokes.GetEx(pressure);
    ExchangeCL& ExV= Stokes.GetEx(velocity);

    vidx->SetFE( vecP2_FE);
    pidx->SetFE( P1_FE);

    if (ProcCL::IamMaster()){
        std::cerr << line << std::endl;
        std::cerr << " - Numbering DOFs ... \n";
    }

    // erzeuge Nummerierung zu diesem Index und fuelle auch die ExchangeCL
    Stokes.CreateNumberingVel(MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(MG.GetLastLevel(), pidx);

    if (C.printInfo){
        if (ProcCL::IamMaster())
            std::cerr << "   + ExchangeCL size for velocity:\n";
        Stokes.GetEx(velocity).SizeInfo(std::cerr);
        if (ProcCL::IamMaster())
            std::cerr << "\n   + ExchangeCL size for pressure:\n";
        Stokes.GetEx(pressure).SizeInfo(std::cerr);
    }

    // Teile den numerischen Daten diese Nummerierung mit
    b->SetIdx(vidx); v->SetIdx(vidx);
    c->SetIdx(pidx); p->SetIdx(pidx);
    A->SetIdx(vidx, vidx); B->SetIdx(pidx, vidx);

    Ulint GPsize_acc = GlobalSum(p->Data.size());
    Ulint GVsize_acc = GlobalSum(v->Data.size());
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);

    if (ProcCL::IamMaster()){
        std::cerr << "  + Number of pressure DOF (accumulated/global):  " <<GPsize_acc<< "/" <<GPsize<< std::endl;
        std::cerr << "  + Number of velocity DOF (accumulated/global):  " <<GVsize_acc<< "/" <<GVsize<< std::endl;
    }

    if (ProcCL::IamMaster())
        std::cerr << line << std::endl << " - Setup matrices and right hand sides ... " << std::endl;

    time.Reset();
    Stokes.SetupSystem(A, b, B, c);
    time.Stop(); duration=time.GetMaxTime(); ::Times.AddTime(T_disc, time.GetMaxTime());

    size_t A_nonzeros = A->Data.GetFinest().num_acc_nonzeros(),
           B_nonzeros = B->Data.GetFinest().num_acc_nonzeros();

    if (ProcCL::IamMaster())
        std::cerr << "  + "<<A_nonzeros<<" nonzeros (accumulated) in A"
                  << "\n  + "<<B_nonzeros<<" nonzeros (accumulated) in B" << std::endl;

    // Solver for Oseen-Problem
    // Preconditioner for A (must be quite exact, so use PCG)
    typedef ParJac0CL<ExchangeCL> APcT;
    APcT APc(ExV, C.relax);
    ParPCGSolverCL<APcT, ExchangeCL> cgsolver(C.pc_iter, C.pc_rel_tol, ExV, APc, true, true);
    typedef SolverAsPreCL<ParPCGSolverCL<APcT, ExchangeCL> > APcSolverT;
    APcSolverT APcSolver(cgsolver);
    // Preconditioner for Schur-Complement
    typedef ParDummyPcCL<ExchangeCL> SPcT; SPcT Spc(ExP);
    ParInexactUzawaCL<APcSolverT, SPcT, APC_SYM>
            symmSchurSolver(APcSolver, Spc, ExV, ExP, C.outer_iter, C.outer_tol, C.reduction, C.inner_iter);

    if (ProcCL::IamMaster())
        std::cerr << line << std::endl<<" - Solve system with InexactUzawa (PCG for approximate Schur-Complement-Matrix) ..." <<std::endl;

    time.Reset();
    symmSchurSolver.Solve(A->Data.GetFinest(), B->Data.GetFinest(), v->Data, p->Data, b->Data, c->Data);
    time.Stop(); duration=time.GetMaxTime(); ::Times.AddTime(T_solve,time.GetMaxTime());

    VectorCL diff_v(A->Data*v->Data +transp_mul(B->Data,p->Data) - b->Data);
    VectorCL diff_p(B->Data*v->Data - c->Data);
    double real_resid = std::sqrt( ExV.Norm_sq(diff_v, false, true) + ExP.Norm_sq(diff_p, false, true) );

    if (ProcCL::IamMaster())
        std::cerr << "   + Time:       " << duration << std::endl
                  << "   + Steps:      " << symmSchurSolver.GetIter() << std::endl
                  << "   + Resid:      " << symmSchurSolver.GetResid() << std::endl
                  << "   + real Resid: " << real_resid << std::endl;

    if (C.ensight)
    {
        if (ProcCL::IamMaster())
            std::cerr << line << std::endl << " - Write solution out in ensight format ... " << std::endl;

        // Erzeuge ensight case File und geom-File
        EnsightP2SolOutCL ensight( MG, pidx->GetFinestPtr(), C.binary, C.masterOut);
        const string filename= C.ensDir + "/" + C.ensCase;
        const string datgeo= filename+".geo",
        datvel = filename+".vel",
        datpr  = filename+".pr";

        ensight.CaseBegin( string(C.ensCase+".case").c_str());
        ensight.DescribeGeom(   C.geomName.c_str(), datgeo);
        ensight.DescribeScalar( "Pressure", datpr);
        ensight.DescribeVector( "Velocity", datvel);
        ensight.Commit();
        ensight.CaseEnd();
        ensight.putGeom( datgeo);

        ensight.putScalar(datpr, Stokes.GetPrSolution());
        ensight.putVector(datvel, Stokes.GetVelSolution());
    }

    if (C.vtk){
        const std::string filenames=C.vtkDir + "/" + C.vtkName;
        VTKOutCL vtkwriter(MG, "StokesLoesung", 2, filenames.c_str(), C.vtkBinary);
        vtkwriter.PutGeom(0.0);
        vtkwriter.PutScalar("pressure", Stokes.GetPrSolution());
        vtkwriter.PutVector("velocity", Stokes.GetVelSolution());
        vtkwriter.Commit();
        vtkwriter.PutGeom(0.1);
        vtkwriter.PutScalar("pressure", Stokes.GetPrSolution());
        vtkwriter.PutVector("velocity", Stokes.GetVelSolution());
        vtkwriter.Commit();
    }

    if (ProcCL::IamMaster())
        std::cerr << line << std::endl << " - Check solution"<<std::endl;

    Stokes.CheckSolution(v, p, &LsgVel, &DLsgVel, &LsgPr);
}
} // end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

// boundary functions (neumann, dirichlet type)
// used for BndSegCL-object of a UnitCube
inline double neu_val(const DROPS::Point2DCL& p) { return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]); }
inline double dir_val(const DROPS::Point2DCL&) { return 0.0; }

// dirichlet value for planes of cube, that has been cut out
inline double dir_val0(const DROPS::Point2DCL& p) { return (1. - p[0]*p[0])*(1. - p[1]*p[1]); }

int main (int argc, char** argv)
{
    ProcCL Proc(&argc, &argv);
    try
    {
        SetDescriber();
        //DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE/*|XFER_SHOW_MSGSALL*/);

        if (argc<2 && ProcCL::IamMaster()){
            std::cerr << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param && ProcCL::IamMaster()){
            std::cerr << "error while opening parameter file: "<<argv[1]<<"\n";
            return 1;
        }

        param >> C;
        param.close();
        if (ProcCL::IamMaster())
            std::cerr << C << std::endl;

        DROPS::ParTimerCL time, alltime;

        // Initialisierung der parallelen Strukturen und Mitteilung, dass es zwei Unbekannten-Index (pressure und velocity) geben wird
        DROPS::ParMultiGridCL pmg;

        DROPS::Point3DCL orig(0.0);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= e2[1]= e3[2]= M_PI/4.;

//      DROPS::BrickBuilderCL brick(null, e1, e2, e3, 3, 3, 3);
        const bool IsNeumann[6]=
            {false, false, false, false, false, false};
        const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
            { &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel};

        if (ProcCL::IamMaster())
        {
            std::cerr << line << std::endl;
            std::cerr << " - Create init grid and distribute ... \n";
        }
        time.Reset();
        DROPS::MGBuilderCL * mgb;
        if (ProcCL::IamMaster())
            mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
        else
            mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);

        // Setup the problem
        StokesOnBrickCL prob(*mgb, StokesCoeffCL(C.nu), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
        DROPS::MultiGridCL &mg = prob.GetMG();
        pmg.AttachTo(mg);

        // Init the LoadBalHandler (create an MultiGridCL and distribute the multigrid)
        DROPS::LoadBalHandlerCL lb(mg);
        lb.SetDebugMode(C.printInfo);

        lb.DoInitDistribution(ProcCL::Master());
        time.Stop(); Times.AddTime(T_init,time.GetMaxTime());
        Times.IncCounter(lb.GetMovedMultiNodes());

        switch (C.refineStrategy){
            case 0 : lb.SetStrategy(DROPS::NoMig);     break;
            case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
            case 2 : lb.SetStrategy(DROPS::Recursive); break;
        }

        if (C.printInfo){
            if (ProcCL::IamMaster())
                std::cerr << " - Distribution of elements:\n";
            mg.SizeInfo(cerr);
        }

        if (ProcCL::IamMaster()){
            std::cerr << line << std::endl;
            std::cerr << " - Refine the grid "<<C.refall<<" regulary\n   and use the following load balancing strategy: ";
            switch (C.refineStrategy){
                case 0: std::cerr << "No LoadBalancing\n"; break;
                case 1: std::cerr << "adaptive\n"; break;
                case 2: std::cerr << "PartKWay\n"; break;
                default: std::cerr << "unknown strategy\nusing no strategy"; C.refineStrategy=0;
            }
        }

        time.Reset();
        for (int ref=0; ref<C.refall; ++ref)
        {
            // Markieren und verfeinern
            if (ProcCL::IamMaster())
                std::cerr << "   + Refine all ("<<ref<<") ";
            DROPS::MarkAll(mg);
            pmg.Refine();
            if (ProcCL::IamMaster())
                std::cerr << "and migrate ...\n";
            lb.DoMigration();
            Times.IncCounter(lb.GetMovedMultiNodes());
        }
        time.Stop(); Times.AddTime(T_ref,time.GetMaxTime());

        if (C.printInfo){
            if (ProcCL::IamMaster())
                std::cerr << " - Distribution of elements:\n";
            mg.SizeInfo(cerr);
        }

        DROPS::Strategy(prob, pmg);

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        if (ProcCL::IamMaster())
            std::cerr << line << std::endl;
        Times.Print(cerr);

        if (ProcCL::IamMaster())
            std::cerr << line<<std::endl<<" - Check parallel multigrid ... ";
        CheckParMultiGrid(pmg);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
