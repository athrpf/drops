//**************************************************************************
// File:    TestInstatStokesPar.cpp                                        *
// Content: parallel solver for instat. Stokes problems                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestInstatStokesPar.cpp
/// \brief Testing the parallel solving of instat. Stokes problems

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

 // include problem class
#include "stokes/stokes.h"
#include "stokes/integrTime.h"

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
DROPS::ParamParInstatStokesCL C;
DROPS::TimeStoreCL Times(5);
const char line[] ="------------------------------------------------------------";

/****************************************************************************
    * S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_init,  "Initialization");
    Times.SetDescriber(T_ref,   "Refinement");
    Times.SetDescriber(T_ex,    "Create ExchangeCL and numbering");
    Times.SetDescriber(T_disc,  "Discretize");
    Times.SetDescriber(T_solve, "Solve");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
*****************************************************************************
*   Checkt, ob die parallele Verteilung und die MultiGrid-Struktur gesund   *
*   zu sein scheint.                                                        *
****************************************************************************/
bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    ofstream check(dat);
    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);
    check.close();
    return DROPS::ProcCL::Check(pmg_sane && mg_sane);
}

/****************************************************************************
* S O L U T I O N                                                           *
****************************************************************************/
inline DROPS::SVectorCL<3> helper( const DROPS::Point3DCL& p)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=     std::cos(p[0])*std::sin(p[1])*std::sin(p[2]);
    ret[1]=     std::sin(p[0])*std::cos(p[1])*std::sin(p[2]);
    ret[2]= -2.*std::sin(p[0])*std::sin(p[1])*std::cos(p[2]);
    return ret;
}

DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double t)
{
    return t*helper(p);
}

double LsgPr(const DROPS::Point3DCL&, double)
{
    return 0;
}

class StokesCoeffCL
{
    public:
        static double q(const DROPS::Point3DCL&) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double t)
            { return (1 + 3*t)*helper(p); }
        const double nu;

    StokesCoeffCL(const double Nu=1e-6/*, const DROPS::Point3DCL& G=DROPS::Point3DCL()*/) : nu(Nu)/*, g(G)*/ {}
};

typedef DROPS::StokesP2P1CL<StokesCoeffCL> StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{
    return DROPS::SVectorCL<3>(0.);
}

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x   = p[0]*(C.dx-p[0]) / (C.dx/2.*C.dx/2.),
                 z   = p[2]*(C.dx-p[2]) / (C.dx/2.*C.dx/2.),
                 freq= C.frequence,
                 ampl= C.ampl;

    ret[1]= x * z * C.inflowVel * (1-ampl*cos(2*M_PI*freq*t));  // Rohr
    return ret;
}

namespace DROPS // for Strategy
{
/****************************************************************************
* S T R A T E G Y                                                           *
*****************************************************************************
*   Das ist die Strategy, um die Poisson-Gleichung, die durch die           *
*   Koeffizienten von oben gegeben sind, zu lösen. Es werden lineare FE     *
*   und das CG-Verfahren für das Lösen des linearen Gleichungssystems       *
*   verwendet. Es wird nicht adaptiv vorgegangen                            *
****************************************************************************/
template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, const ParMultiGridCL& /*pmg*/)
{
    typedef StokesP2P1CL<Coeff> StokesProblemT;
    ParTimerCL timer;
    double duration;
    double t= 0.;
    const double dt= C.stepsize;
    Stokes.t= 0;

    MultiGridCL& MG= Stokes.GetMG();

    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    VelVecDescCL* v= &Stokes.v;
    VecDescCL*    p= &Stokes.p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;

    MLMatDescCL* A= &Stokes.A;
    MLMatDescCL* B= &Stokes.B;
    MLMatDescCL* M= &Stokes.M;

    ExchangeCL&      ExV= Stokes.vel_idx.GetEx();
    ExchangeCL&      ExP= Stokes.pr_idx.GetEx();

    vidx->SetFE( vecP2_FE);
    pidx->SetFE( P1_FE);

    if (ProcCL::IamMaster()){
        std::cerr << line << std::endl;
        std::cerr << " - Numbering DOFs ... \n";
    }

    // erzeuge Nummerierung zu diesem Index und fuelle die ExchangeCL
    timer.Reset();
    Stokes.CreateNumberingVel(MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(MG.GetLastLevel(), pidx);
    timer.Stop(); ::Times.AddTime(T_ex, timer.GetMaxTime());

    // Teile den numerischen Daten diese Nummerierung mit
    b->SetIdx(vidx); v->SetIdx(vidx);
    c->SetIdx(pidx); p->SetIdx(pidx);
    A->SetIdx(vidx, vidx);
    B->SetIdx(pidx, vidx);
    M->SetIdx(vidx, vidx);

    Ulint GPsize_acc = ProcCL::GlobalSum(p->Data.size());
    Ulint GVsize_acc = ProcCL::GlobalSum(v->Data.size());
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);

    if (ProcCL::IamMaster()){
        std::cerr << "  + Number of pressure DOF (accumulated/global):  " << GPsize_acc <<"/"<<GPsize << std::endl;
        std::cerr << "  + Number of velocity DOF (accumulated):  " << GVsize_acc <<"/"<<GVsize << std::endl;
    }

        // Erzeuge ensight case File und geom-File
    EnsightP2SolOutCL *ensight=0;
    const string EnsCase = C.EnsCase;
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datvel= filename+".vel",
                 datpr = filename+".pr";

    if (C.ensight){
        ensight = new EnsightP2SolOutCL( MG, pidx->GetFinestPtr(), false);
        ensight->CaseBegin( string(EnsCase+".case").c_str(), C.timesteps);
        ensight->DescribeGeom(   (C.geomName).c_str(), datgeo);
        ensight->DescribeScalar( "Pressure", datpr, true);
        ensight->DescribeVector( "Velocity", datvel, true);
        ensight->putGeom( datgeo);
    }

    if (C.printInfo){
        if (ProcCL::IamMaster())
            std::cerr << '\n' << "   + ExchangeCL size for velocity:\n";
        ExV.SizeInfo(std::cerr);
        if (ProcCL::IamMaster())
            std::cerr << "\n   + ExchangeCL size for pressure:\n";
        ExP.SizeInfo(std::cerr);
    }

    if (ProcCL::IamMaster())
        std::cerr << line << std::endl << " Seting up matrices ...\n";
    timer.Reset();
    Stokes.SetupInstatSystem(A,B,M);
    timer.Stop(); duration = timer.GetMaxTime();
    ::Times.AddTime(T_disc,duration);
    if (ProcCL::IamMaster())
        std::cerr << "     ---> " << duration << " sec \n";

    Stokes.InitVel(v,&Null,0);

      // Preconditioner
//     typedef ParJac0CL<ExchangeCL> APcT;      APcT APc(ExV, C.relax);
//     typedef ParDummyPcCL<ExchangeCL> BPcT;   BPcT BPc(ExP);
//     typedef BlockPreCL<APcT, BPcT> OseenPCT; OseenPCT OseenPC(APc, BPc);

      // Base Solver for the Oseen problem
//     typedef ParPreGCRSolverCL<OseenPCT, ExchangeBlockCL> OseenBaseSolT;
//     OseenBaseSolT OseenBaseSol(C.restart,C.iter,C.tol, &Ex, OseenPC, false);

//     typedef ParPCGSolverCL<OseenPCT,ExchangeBlockCL> OseenBaseSolT;
//     OseenBaseSolT OseenBaseSol(C.outer_iter,C.outer_tol, Ex, OseenPC, C.relative, C.accur);

     // Solver for the stat. Stokes problem
//     typedef BlockMatrixSolverCL<OseenBaseSolT> SolverT;
//     SolverT Solver(OseenBaseSol);
    typedef ParDummyPcCL SPcT;
    SPcT ispc( Stokes.pr_idx.GetFinest());
    typedef ParJac0CL  APcPcT;
    APcPcT Apcpc( Stokes.vel_idx.GetFinest());
    typedef ParPCGSolverCL<APcPcT> ASolverT;
    ASolverT Asolver( 500, 0.02, Stokes.vel_idx.GetFinest(), Apcpc, /*relative=*/ true, /*accur*/ true);
    typedef SolverAsPreCL<ASolverT> APcT;
    APcT Apc( Asolver/*, &std::cerr*/);
    typedef ParInexactUzawaCL<APcT, SPcT, APC_SYM> OseenSolverT;
    OseenSolverT Solver( Apc, ispc, Stokes.vel_idx.GetFinest(), Stokes.pr_idx.GetFinest(), C.outer_iter, C.outer_tol, 0.1);

      // Solver for the time integration
    typedef InstatStokesThetaSchemeCL<StokesP2P1CL<Coeff>, OseenSolverT> InstatSolverT;
    InstatSolverT instatsolver(Stokes, Solver, C.theta);

    instatsolver.SetTimeStep(dt);
    VectorCL diff_v(A->Data.num_rows());
    VectorCL diff_p(B->Data.num_rows());

    for (int time=0; time<C.timesteps; ++time, t+=dt)
    {
        if (ProcCL::IamMaster())
            std::cerr << line<< '\n' << "  Solving timestep " << time << ", time " << t << '\n'
                      << "  ---------------- " << '\n' << '\n'
                      << " - Solving stat. Stokes problem" << std::endl;
        timer.Reset();
        instatsolver.DoStep(v->Data, p->Data);
        timer.Stop(); duration = timer.GetMaxTime();
        ::Times.AddTime(T_solve, duration);
        if (ProcCL::IamMaster())
            std::cerr << "   + Time     " << duration << '\n'
                      << "   + Steps    " << Solver.GetIter() << '\n'
                      << "   + Resid    " << Solver.GetResid() <<'\n'<< std::endl;

        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,DROPS::VecDescCL> pr(&Stokes.p, &Stokes.GetBndData().Pr, &MG);
        if (C.ensight){
            ensight->putScalar(datpr, Stokes.GetPrSolution(),t);
            ensight->putVector(datvel, Stokes.GetVelSolution(),t);
            ensight->Commit();
        }

        if (ProcCL::IamMaster())
            std::cerr << " - Checking solution ...\n";

        Stokes.CheckSolution(v,p,&LsgVel,&LsgPr,t);
    }

    if (C.ensight){
        ensight->CaseEnd();
        delete ensight;
    }
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

void MarkInflow (DROPS::MultiGridCL& mg, double e2, double percent_of_ref)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin()),
         ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
    {
        if ( e2-GetBaryCenter(*It)[1]<(e2*percent_of_ref/100.) )
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
    DROPS::ProcCL::Instance(&argc, &argv);
    try
    {
        SetDescriber();
        //DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE/*|XFER_SHOW_MSGSALL*/);

        if (argc<2 && DROPS::ProcCL::IamMaster()){
            std::cerr << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param && DROPS::ProcCL::IamMaster()){
            std::cerr << "error while opening parameter file: "<<argv[1]<<"\n";
            return 1;
        }

        param >> C;
        param.close();
        if (DROPS::ProcCL::IamMaster())
            std::cerr << C << std::endl;

        DROPS::ParTimerCL time, alltime;
        double duration;

        // Initialisierung der parallelen Strukturen
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

        DROPS::Point3DCL orig(0.0);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0] = C.dx; e2[1] = C.dy; e3[2]=C.dz;

        const bool IsNeumann[6]=
                {false, false, true, false, false, false};
//              {false, false, true, false, true, true};
        const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
                { &Null, &Null, &Null, &Inflow, &Null, &Null};

        StokesOnBrickCL *prob=0;

        if (DROPS::ProcCL::IamMaster())
        {
            std::cerr << line << std::endl;
            std::cerr << " - Create init grid and distribute ... \n";
        }
        time.Reset();
        if (DROPS::ProcCL::IamMaster())
        {
            DROPS::BrickBuilderCL brick(orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
            prob= new StokesOnBrickCL(brick, StokesCoeffCL(C.nu), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
            prob->GetMG().MakeConsistentNumbering();
        }
        else
        {
            DROPS::EmptyBrickBuilderCL emptyBrick( orig, e1, e2, e3, C.basicref_x);
            prob= new StokesOnBrickCL(emptyBrick, StokesCoeffCL(C.nu), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
        }
        // Weise dem parallelen MultiGrid das MultiGrid zu
        DROPS::MultiGridCL& mg = prob->GetMG();
        pmg.AttachTo(mg);

        // Init the LoadBalHandler (create an MultiGridCL and distribute the multigrid)
        DROPS::LoadBalHandlerCL lb(mg);
        lb.DoInitDistribution(DROPS::ProcCL::Master());
        lb.SetXferUnknowns(true);
        switch (C.refineStrategy){
            case 0 : lb.SetStrategy(DROPS::NoMig);     break;
            case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
            case 2 : lb.SetStrategy(DROPS::Recursive); break;
        }
        time.Stop(); duration = time.GetMaxTime(); Times.AddTime(T_init,duration);
        if (DROPS::ProcCL::IamMaster())
            std::cerr << "     ---> " << duration << " sec \n";
        if (C.printInfo)
            mg.SizeInfo(cerr);


        if (DROPS::ProcCL::IamMaster()){
            std::cerr << line << std::endl;
            std::cerr << " - Refine the grid "<<C.refall<<" at inflow\n   and use the following load balancing strategy: ";
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
            if (DROPS::ProcCL::IamMaster()){
                std::cerr << "   + Refine inflow ("<<ref<<") and migrate ...\n";
            }
            MarkInflow(mg,C.dy,10.0);
            pmg.Refine();
            lb.DoMigration();
        }

        time.Stop(); duration=time.GetMaxTime(); Times.AddTime(T_ref,duration);
        if (C.printInfo){
            if (DROPS::ProcCL::IamMaster())
                std::cerr << "     ---> " << duration << " sec \n" << " - Verteilung der Elemente:\n";
            mg.SizeInfo(cerr);
        }

        DROPS::Strategy(*prob, pmg);

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        if (DROPS::ProcCL::IamMaster())
            std::cerr << line << std::endl;
        Times.Print(cerr);

        if (DROPS::ProcCL::IamMaster())
            std::cerr << line<<std::endl<<" - Check parallel multigrid ... ";
        if (CheckParMultiGrid(pmg)){
            if (DROPS::ProcCL::IamMaster()) cerr << " OK\n";
        }
        else
            if (DROPS::ProcCL::IamMaster()) cerr << "  !!! not OK !!!\n";
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
