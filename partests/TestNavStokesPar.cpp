//**************************************************************************
// File:    TestNavStokesPar.cpp                                           *
// Content: parallel solver for instat Navier-Stokes problem               *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestNavStokesPar.cpp
/// \brief Testing parallel solvers for the instat. Navier-Stokes problem

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
#include "num/nssolver.h"
#include "levelset/surfacetension.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"

 // include problem class
#include "navstokes/navstokes.h"
#include "stokes/instatstokes2phase.h"  // for repair classes
#include "navstokes/integrTime.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"


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

DROPS::ParamParNavStokesCL C;
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
          std::cout << " As far as I can tell, the multigrid is sane\n";
    }
    else
        throw DROPS::DROPSErrCL("Found error in multigrid!");
}

/****************************************************************************
* S T O K E S  C O E F F  C L                                               *
*****************************************************************************
*   Koeffizienten der der Navier-Stokes-Gleichung. Es soll eine Driven-     *
*   Cavity simuliert werden. Daher sind f und q identisch 0. Die Navier-    *
*   Stokes-Gleichung lautet:                                                *
*****************************************************************************
*   ( u * grad ) u - div(grad u) + grad p = 0    im Brick                   *
*                                   div u = 0    im Brick                   *
*                                       u = flow auf Deckel vom Brick       *
*                                       u = 0    auf restl. Rand            *
****************************************************************************/
struct InstatNSCL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]= (2.*t - 1.)*p[0];
        ret[1]= (2.*t - 1.)*p[1];
        ret[2]= (2. - 4.*t)*p[2];
        return ret;
    }

    // int_{x=0..1, y=0..1,z=0..1} p(x,y,z,t) dx dy dz = 0 for all t.
    static double LsgPr(const DROPS::Point3DCL& p, double t)
    {
        return (0.5 - t)*p.norm_sq() + t - 0.5;
    }

    // Jacobi-matrix of exact solution (only in the spatial variables)
    static inline DROPS::SMatrixCL<3, 3> DxLsgVel(const DROPS::Point3DCL&, double t)
    {
        DROPS::SMatrixCL<3, 3> ret(0.0);
        ret(0,0)= 2.*t - 1.;
        ret(1,1)= 2.*t - 1.;
        ret(2,2)= 2. - 4.*t;
        return ret;
    }

    // Time-derivative of exact solution
    static inline DROPS::SVectorCL<3> DtLsgVel(const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.0);
        ret[0]= 2.*p[0];
        ret[1]= 2.*p[1];
        ret[2]= -4.*p[2];
        return ret;
    }
    // u_t + q*u - nu*laplace u + (u*D)u + Dp = f
    //                                 -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&, double) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double t)
        {
            DROPS::SVectorCL<3> ret;
            ret[0]= (4.*t*t - 6.*t + 4.)*p[0];
            ret[1]= (4.*t*t - 6.*t + 4.)*p[1];
            ret[2]= (16.*t*t -18.*t + 1.)*p[2];
            return ret;
        }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };

    static StokesCoeffCL Coeff;
};

InstatNSCL::StokesCoeffCL InstatNSCL::Coeff;
typedef InstatNSCL MyPdeCL;
typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&, double);

const double radiusorbit= 0.3; // Radius of the drops' orbit.
const double radiusdrop= 0.15; // Initial radius of the drop.

// positive outside the drop, negative inside the drop.
double
SignedDistToInterface(const DROPS::Point3DCL& p, double t)
{
   DROPS::Point3DCL c;
   c[0]= 0.5 + radiusorbit*std::cos( 2.*M_PI*t);
   c[1]= 0.5 + radiusorbit*std::sin( 2.*M_PI*t);
   c[2]= 0.5;
   return (p-c).norm() - radiusdrop;
}

double actual_time= 0.;
double SignedDistToInterfaceWOTime(const DROPS::Point3DCL& p)
{
    return SignedDistToInterface(p, actual_time);
}

typedef double (*signed_dist_fun)(const DROPS::Point3DCL& p, double t);


inline DROPS::SVectorCL<3> uD(const DROPS::Point3DCL&, double)
{
    DROPS::SVectorCL<3> ret(0.);
    ret[0]= 1.;
    return ret;
}

inline DROPS::SVectorCL<3> Null(const DROPS::Point3DCL&, double)
{
    DROPS::SVectorCL<3> ret(0.);
    return ret;
}

template<class Coeff>
void
SetMatVecIndices(DROPS::NavierStokesP2P1CL<Coeff>& NS,
                 DROPS::MLIdxDescCL* const vidx,
                 DROPS::MLIdxDescCL* const pidx)
{
    std::cout << "#Druck-Unbekannte: " << pidx->NumUnknowns() << std::endl;
    std::cout << "#Geschwindigkeitsunbekannte: " << vidx->NumUnknowns() << std::endl;
    NS.b.SetIdx( vidx);
    NS.c.SetIdx( pidx);
    NS.cplM.SetIdx( vidx);
    NS.cplN.SetIdx( vidx);
    NS.A.SetIdx( vidx, vidx);
    NS.B.SetIdx( pidx, vidx);
    NS.M.SetIdx( vidx, vidx);
    NS.N.SetIdx( vidx, vidx);
}

template<class Coeff>
void
ResetSystem(DROPS::NavierStokesP2P1CL<Coeff>& NS)
{
    NS.A.Reset(); NS.B.Reset();
    NS.M.Reset(); NS.N.Reset();
    NS.b.Reset(); NS.c.Reset();
    NS.cplM.Reset(); NS.cplN.Reset();
}

namespace DROPS // for Strategy
{
/****************************************************************************
* S T R A T E G Y                                                           *
*****************************************************************************
*   This is the strategy to solve the instat. Navier-Stokes equation.       *
****************************************************************************/
template<class Coeff>
void Strategy(NavierStokesP2P1CL<Coeff> & NavStokes, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    typedef NavierStokesP2P1CL<Coeff> StokesProblemT;
    ParTimerCL time;
    double duration;

    MultiGridCL& MG= NavStokes.GetMG();
    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    LevelsetP2CL lset( MG, sf);                 // levelset just for refinement
    AdapTriangCL adapt(pmg, lb, 0.1, 1, 2);
    actual_time=0.;
    adapt.MakeInitialTriang(::SignedDistToInterfaceWOTime);

    double dt= 1./C.timesteps;

    LevelsetRepairCL lsetrepair( lset, pmg);
    adapt.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( NavStokes, pmg);
    adapt.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( NavStokes, lset, pmg);
    adapt.push_back( &prrepair);

    MLIdxDescCL* vidx= &NavStokes.vel_idx;
    MLIdxDescCL* pidx= &NavStokes.pr_idx;
    IdxDescCL*   lidx= &lset.idx;

    VelVecDescCL* v   = &NavStokes.v;
    VecDescCL*    p   = &NavStokes.p;
    VelVecDescCL* b   = &NavStokes.b;
    VecDescCL*    c   = &NavStokes.c;
    VelVecDescCL* cplN= &NavStokes.cplN;
    VelVecDescCL* cplM= &NavStokes.cplM;

    MLMatDescCL* A= &NavStokes.A;
    MLMatDescCL* B= &NavStokes.B;
    MLMatDescCL* M= &NavStokes.M;
    MLMatDescCL* N= &NavStokes.N;

    const Uint pressure=NavStokes.pressure,
               velocity=NavStokes.velocity;
    ExchangeBlockCL& Ex   =NavStokes.GetEx();
//     ExchangeCL&      ExVel=NavStokes.GetEx(velocity);
//     ExchangeCL&      ExPr =NavStokes.GetEx(pressure);

    vidx->SetFE( vecP2_FE);
    pidx->SetFE( P1_FE);

    if (ProcCL::IamMaster()){
        std::cout << line << std::endl;
        std::cout << " - Numbering DOFs ... \n";
    }

    // erzeuge Nummerierung zu diesem Index und fülle auch die ExchangeCL
    NavStokes.CreateNumberingVel(MG.GetLastLevel(), vidx);
    NavStokes.CreateNumberingPr(MG.GetLastLevel(), pidx);
    lset.CreateNumbering(MG.GetLastLevel(), lidx);

    if (C.printInfo){
        if (ProcCL::IamMaster())
            std::cout << "   + ExchangeCL size for velocity:\n";
        NavStokes.GetEx(velocity).SizeInfo(std::cout);
        if (ProcCL::IamMaster())
            std::cout << "\n   + ExchangeCL size for pressure:\n";
        NavStokes.GetEx(pressure).SizeInfo(std::cout);
    }

    // Ensight
    EnsightP2SolOutCL* ensight=0;
    const string EnsCase=C.ensCase;
    const string filename= C.ensDir +"/"+ EnsCase;
    const string datgeo= filename+".geo",
          datvel = filename+".vel",
          datpr  = filename+".pr";

    if (C.ensight)
    {
        if (ProcCL::IamMaster())
            std::cout << line << std::endl << " - Create ensight case ... " << std::endl;

        // Erzeuge ensight case File und geom-File
        ensight = new EnsightP2SolOutCL( MG, pidx->GetFinestPtr(), false);
        ensight->CaseBegin( string(EnsCase+".case").c_str(), C.timesteps);
        ensight->DescribeGeom(   "Messzelle", datgeo, true);
        ensight->DescribeScalar( "Pressure", datpr, true);
        ensight->DescribeVector( "Velocity", datvel, true);
        ensight->putGeom( datgeo);
    }

    // Teile den numerischen Daten diese Nummerierung mit
    b->SetIdx(vidx);    v->SetIdx(vidx);
    c->SetIdx(pidx);    p->SetIdx(pidx);
    cplN->SetIdx(vidx); cplM->SetIdx(vidx);
    A->SetIdx(vidx, vidx); B->SetIdx(pidx, vidx);
    N->SetIdx(vidx, vidx); M->SetIdx(vidx, vidx);
    lset.Phi.SetIdx(lidx);

    Ulint GPsize_acc = GlobalSum(p->Data.size());
    Ulint GVsize_acc = GlobalSum(v->Data.size());
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);

    if (ProcCL::IamMaster()){
        std::cout << "  + Number of pressure DOF (accumulated/global):  " <<GPsize_acc<< "/" <<GPsize<< std::endl;
        std::cout << "  + Number of velocity DOF (accumulated/global):  " <<GVsize_acc<< "/" <<GVsize<< std::endl;
    }

    NavStokes.InitVel( v, &MyPdeCL::LsgVel);

    // Setup stat. part of matrices
    if (ProcCL::IamMaster())
        std::cout << line << std::endl << " - Setup matrices and right hand sides ... " << std::endl;
    time.Reset();
    NavStokes.SetupInstatSystem(A,B,M);
    NavStokes.SetupNonlinear( N, v, cplN);
    time.Stop(); duration=time.GetMaxTime(); ::Times.AddTime(T_disc, time.GetMaxTime());

    size_t A_nonzeros = A->Data.GetFinest().num_acc_nonzeros(),
           B_nonzeros = B->Data.GetFinest().num_acc_nonzeros(),
           N_nonzeros = N->Data.GetFinest().num_acc_nonzeros(),
           M_nonzeros = M->Data.GetFinest().num_acc_nonzeros();

    if (ProcCL::IamMaster())
        std::cout << "  + "<<A_nonzeros<<" nonzeros (accumulated) in A"<<'\n'
                  << "  + "<<B_nonzeros<<" nonzeros (accumulated) in B"<<'\n'
                  << "  + "<<N_nonzeros<<" nonzeros (accumulated) in N"<<'\n'
                  << "  + "<<M_nonzeros<<" nonzeros (accumulated) in M"<<std::endl;

    // Solver for Oseen-Problem
      // Preconditioner for A
    typedef ParJac0CL<ExchangeCL> APcT; APcT APc(NavStokes.GetEx(velocity), C.relax);
//  typedef ParDummyPcCL APcT; APcT APc;
      // Preconditioner for B
    typedef ParDummyPcCL<ExchangeCL> BPcT; BPcT BPc(NavStokes.GetEx(pressure));
      // Preconditioner for Oseen-Matrix
    typedef BlockPreCL<APcT, BPcT> OseenPCT; OseenPCT OseenPC(APc, BPc);
      // Base Solver for the Oseen-Problem1
    typedef ParPreGMResSolverCL<OseenPCT,ExchangeBlockCL> OseenBaseSolT;
    OseenBaseSolT OseenBaseSol(C.restart, C.outer_iter, C.outer_tol, Ex, OseenPC, /*relative*/C.relative, C.accur);
//     typedef ParPCGSolverCL<OseenPCT,ExchangeBlockCL> OseenBaseSolT;
//     OseenBaseSolT OseenBaseSol(/*C.restart, */C.outer_iter, C.outer_tol, Ex, OseenPC,/*relative*/C.relative, C.accur);
      // Solver for the Oseen-Problem
    typedef BlockMatrixSolverCL<OseenBaseSolT> SolverT;
    SolverT Solver(OseenBaseSol);

//     typedef ParPreGMResSolverCL<APcT, ExchangeCL> innerSolT;
//     innerSolT innerSolver(C.restart, C.inner_iter, C.inner_tol, ExVel, APc, /*ModGS*/false, /*relative*/C.relative);
//     typedef ParPreGMResSolverCL<BPcT, ExchangeCL> outerSolT;
//     outerSolT outerSolver(C.restart, C.outer_iter, C.outer_tol, ExPr, BPc, /*ModGS*/false, /*relative*/C.relative);
//     typedef PSchurSolver2CL<innerSolT, outerSolT> SchurSolT;
//     SchurSolT SchurSolver(innerSolver, outerSolver, C.outer_iter, C.outer_tol);

    typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
    NSSolverT *statsolver=0; // (NavStokes, Solver, C.nav_iter, C.nav_tol, C.reduction);

    typedef InstatNavStokesThetaSchemeCL<StokesProblemT, NSSolverT> NSThetaCL;
    NSThetaCL *instatsolver= 0; // thetaScheme(NavStokes, nssolver, C.theta);

    for (int step=0; step<C.timesteps; ++step, actual_time+=dt, NavStokes.t+= dt)
    {
        if (ProcCL::IamMaster())
            std::cout << line<< '\n' << "  Solving timestep " << step << ", time " << NavStokes.t << '\n'
                      << "  ---------------- " << '\n' << '\n';

        if (step%(C.timesteps/10) == 0) { // modify the grid
            if (step>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                ResetSystem( NavStokes);
                lset.Init(::SignedDistToInterfaceWOTime);
                adapt.UpdateTriang(lset);
            }
            SetMatVecIndices( NavStokes, vidx, pidx);
            time.Reset(); time.Start();
            NavStokes.SetupInstatSystem( &NavStokes.A, &NavStokes.B, &NavStokes.M);
            time.Stop();
            std::cout << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();

            statsolver= new NSSolverT( NavStokes, Solver, C.nav_iter, C.nav_tol, C.reduction);

            NavStokes.SetupNonlinear( &NavStokes.N, v, &NavStokes.cplN, actual_time, actual_time);
            time.Stop();
            std::cout << "SetupNonlinear: " << time.GetTime() << " seconds" << std::endl;
            NavStokes.SetupInstatRhs( &NavStokes.b, &NavStokes.c, &NavStokes.cplM, actual_time, &NavStokes.b, actual_time);
            instatsolver= new NSThetaCL( NavStokes, *statsolver, C.theta);
        }
        NavStokes.SetTime( actual_time+dt); // We have to set the new time!
        instatsolver->SetTimeStep( dt);
        std::cout << "Before timestep." << std::endl;
        instatsolver->DoStep( *v, p->Data);
        std::cout << "After timestep." << std::endl;
        NavStokes.CheckSolution( v, vidx, p, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, actual_time+dt);
        if (C.ensight){
            ensight->putGeom(   datgeo, actual_time);
            ensight->putScalar( datpr, NavStokes.GetPrSolution(), actual_time);
            ensight->putVector( datvel, NavStokes.GetVelSolution(), actual_time);
        }
    }
    delete statsolver; statsolver= 0;
    delete instatsolver; instatsolver= 0;
    ResetSystem( NavStokes);

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

void MarkTop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel, double range)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        DROPS::Point3DCL TetraMitte= GetBaryCenter(*It);
        if ( (C.dz-TetraMitte[2])/C.dz < range ){
            It->SetRegRefMark();
        }
    }
}

int main (int argc, char** argv)
{
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
  try
  {
    SetDescriber();
    //DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE/*|XFER_SHOW_MSGSALL*/);

    if (argc<2 && ProcCL::IamMaster()){
        std::cout << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
    }
    std::ifstream param( argv[1]);
    if (!param && ProcCL::IamMaster()){
        std::cout << "error while opening parameter file: "<<argv[1]<<"\n";
        return 1;
    }

    param >> C;
    param.close();
    if (ProcCL::IamMaster())
        std::cout << C << std::endl;

    DROPS::ParTimerCL time, alltime;

    // Initialisierung der parallelen Strukturen
    DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

    DROPS::Point3DCL orig(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]=C.dx; e2[1]=C.dy; e3[2]=C.dz;

    //      DROPS::BrickBuilderCL brick(null, e1, e2, e3, 3, 3, 3);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
          &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };

    if (ProcCL::IamMaster())
    {
        std::cout << line << std::endl;
        std::cout << " - Create init grid and distribute ... \n";
    }

    time.Reset();
    DROPS::MGBuilderCL * mgb;
    if (ProcCL::IamMaster())
        mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
    else
        mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);

    // Setup the problem
    typedef DROPS::NavierStokesP2P1CL<MyPdeCL::StokesCoeffCL> NSOnBrickCL;
    typedef NSOnBrickCL MyNavierStokesCL;
    MyNavierStokesCL prob(*mgb, MyPdeCL::StokesCoeffCL(),
                          DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));

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
            std::cout << " - Distribution of elements:\n";
        mg.SizeInfo(cout);
    }

    if (ProcCL::IamMaster()){
        std::cout << line << std::endl;
        std::cout << " - Refine the grid "<<C.refall<<" regulary\n   and use the following load balancing strategy: ";
        switch (C.refineStrategy){
            case 0: std::cout << "No LoadBalancing\n"; break;
            case 1: std::cout << "adaptive\n"; break;
            case 2: std::cout << "PartKWay\n"; break;
            default: std::cout << "unknown strategy\nusing no strategy"; C.refineStrategy=0;
        }
    }

    time.Reset();
    for (int ref=0; ref<C.refall+C.markTop; ++ref)
    {
        // Markieren und verfeinern
        if (ref<C.refall){
            if (ProcCL::IamMaster())
                std::cout << "   + Refine all ("<<ref<<") ";
            DROPS::MarkAll(mg);
        }
        else{
            if (ProcCL::IamMaster())
                std::cout << "   + Refine top ("<<ref<<") ";
            MarkTop(mg, mg.GetLastLevel(), 0.4);
        }
        pmg.Refine();
        if (ProcCL::IamMaster())
            std::cout << "and migrate ...\n";
        lb.DoMigration();
        Times.IncCounter(lb.GetMovedMultiNodes());
    }
    time.Stop(); Times.AddTime(T_ref,time.GetMaxTime());

    if (C.printInfo){
        if (ProcCL::IamMaster())
            std::cout << " - Distribution of elements:\n";
        mg.SizeInfo(cout);
    }

    DROPS::Strategy(prob, pmg, lb);

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    if (ProcCL::IamMaster())
        std::cout << line << std::endl;
    Times.Print(cout);

    if (ProcCL::IamMaster())
        std::cout << line<<std::endl<<" - Check parallel multigrid ... ";
    CheckParMultiGrid(pmg);

    return 0;
}
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
