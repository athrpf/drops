//**************************************************************************
// File:    TestMzellePar.cpp                                              *
// Content: parallel solver for instat Navier-Stokes problem               *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. September 2006                                             *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestMzellePar.cpp
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
#include "num/parstokessolver.h"
#include "num/nssolver.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/gridOut.h"

 // include problem class
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "levelset/coupling.h"
#include "parallel/parfastmarch.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#endif

// **************************************************************
// *  C H E C K  P A R  M U L T I G R I D                       *
// **************************************************************
// * Check the parallel multigrid and write sanity-results in   *
// * a file.                                                    *
// **************************************************************
void CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream check(dat);
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


// **************************************************************
// *  D I S P L A Y  U N K S                                    *
// **************************************************************
// * Display statistics about unknowns                          *
// **************************************************************
using DROPS::Ulint;
using DROPS::GlobalSum;
using DROPS::GlobalMax;
using DROPS::GlobalMin;
void DisplayUnks(const DROPS::MLIdxDescCL* vidx, const DROPS::MLIdxDescCL* pidx, const DROPS::IdxDescCL* lidx,
                 const DROPS::ExchangeCL& ExV, const DROPS::ExchangeCL& ExP, const DROPS::ExchangeCL& ExL,
                 const DROPS::MultiGridCL& MG)
{
    // local number on unknowns
    Ulint Psize      = pidx->NumUnknowns();
    Ulint Vsize      = vidx->NumUnknowns();
    Ulint Lsize      = lidx->NumUnknowns();

    // global number of unknowns
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);
    Ulint GLsize     = lidx->GetGlobalNumUnknowns(MG);

    // accumulated size of unknwons
    Ulint Psize_acc = GlobalSum(Psize);
    Ulint Vsize_acc = GlobalSum(Vsize);
    Ulint Lsize_acc = GlobalSum(Lsize);

    // maximal and minimal number of unknowns
    Ulint P_min= GlobalMin(Psize); Ulint P_max= GlobalMax(Psize);
    Ulint V_min= GlobalMin(Vsize); Ulint V_max= GlobalMax(Vsize);
    Ulint L_min= GlobalMin(Lsize); Ulint L_max= GlobalMax(Lsize);

    // ratios between maximal number of unknowns/proc and minimal number
    double P_ratio   = (double)P_max/(double)P_min;
    double V_ratio   = (double)V_max/(double)V_min;
    double L_ratio   = (double)L_max/(double)L_min;

    // number on boundaries
    Ulint P_accmax=GlobalMax(ExP.AccDistIndex.size()), P_accmin=GlobalMin(ExP.AccDistIndex.size());
    Ulint V_accmax=GlobalMax(ExV.AccDistIndex.size()), V_accmin=GlobalMin(ExV.AccDistIndex.size());
    Ulint L_accmax=GlobalMax(ExL.AccDistIndex.size()), L_accmin=GlobalMin(ExL.AccDistIndex.size());

    // ratio of these unknowns
    double P_accratio= (double)P_accmax / (double)P_accmin;
    double V_accratio= (double)V_accmax / (double)V_accmin;
    double L_accratio= (double)L_accmax / (double)L_accmin;

    // output on screen
    if (DROPS::ProcCL::IamMaster()){
        std::cerr << "  + Number of DOF\n        "
                  << std::setw(10)<<"global"<<std::setw(10)<<"accum"<<std::setw(10)
                  << "max"<<std::setw(10)<<"min"<<std::setw(10)<<"ratio"<<"  |  "
                  << std::setw(10)<<"max_acc" <<std::setw(10)<<"min_acc"<<std::setw(10)<<"ratio_acc"<<std::endl;

        std::cerr << "    "<<"pr  "
                  << std::setw(10)<<GPsize<<std::setw(10)<<Psize_acc<<std::setw(10)<<P_max
                  << std::setw(10)<<P_min<< std::setw(10)<<P_ratio<<"  |  "
                  << std::setw(10)<<P_accmax<<std::setw(10)<<P_accmin<<std::setw(10)<<P_accratio<<std::endl;

        std::cerr << "    "<<"vel "
                  << std::setw(10)<<GVsize<<std::setw(10)<<Vsize_acc<<std::setw(10)<<V_max
                  << std::setw(10)<<V_min<< std::setw(10)<<V_ratio<<"  |  "
                  << std::setw(10)<<V_accmax<<std::setw(10)<<V_accmin<<std::setw(10)<<V_accratio<<std::endl;

        std::cerr << "    "<<"scl "
                  << std::setw(10)<<GLsize<<std::setw(10)<<Lsize_acc<<std::setw(10)<<L_max
                  << std::setw(10)<<L_min<< std::setw(10)<<L_ratio<<"  |  "
                  << std::setw(10)<<L_accmax<<std::setw(10)<<L_accmin<<std::setw(10)<<L_accratio<<std::endl;

        std::cerr << std::endl;
    }
}


// **************************************************************
// * Parameterfile where all parameters can be found            *
// **************************************************************
DROPS::ParParamMesszelleNsCL C;

// **************************************************************
// * Z E R O  F L O W  C L                                      *
// **************************************************************
// * Coefficients of the PDE                                    *
// * rho*du/dt - mu*laplace u + Dp = f + rho*g - okn            *
// *                        -div u = 0                          *
// *                             u = u0, t=t0                   *
// **************************************************************
class ZeroFlowCL
{
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( C.rhoD, C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muD,  C.muF),   DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

// Not used a.t.m
class DimLessCoeffCL
{
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( 1., C.rhoF/C.rhoD ), DROPS::H_sm, C.sm_eps),
        mu ( DROPS::JumpCL( 1., C.muF/C.muD),    DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma/C.rhoD), g( C.g)    {}
};


// **************************************************************
// * Zero function                                              *
// **************************************************************
DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{
    return DROPS::SVectorCL<3>(0.);
}

// **************************************************************
// Inflow at top of the cell                                    *
// **************************************************************
DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

// **************************************************************
// * D I S T A N C E  F C T                                     *
// **************************************************************
// * distance function for initial droplet                      *
// **************************************************************
double DistanceFct( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.Mitte-p;
    const double avgRad = cbrt(C.Radius[0]*C.Radius[1]* C.Radius[2]);
    d/= C.Radius/avgRad;
    return d.norm()-avgRad;
}

// **************************************************************
// * S U R F A C E   T E N S I O N                              *
// **************************************************************
// * surface tension at droplet interphase                      *
// **************************************************************
double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }

namespace DROPS
{

// **************************************************************
// * S T R A T E G Y                                            *
// **************************************************************
// * flow control                                               *
// **************************************************************
template<class Coeff>
void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes)
{
    ParTimerCL time;        // time measurement
    double duration;        // used time

    // Typedefinition of the problem
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    MultiGridCL& MG= Stokes.GetMG();

    sigma= Stokes.GetCoeff().SurfTens;

    // Create the levelset class
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0,
                       C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff, C.NarrowBand);

    // Index describer for levelset, velocity and pressure
    IdxDescCL*   lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    // Matrices and vectors outside of problem classes
    VecDescCL cplN;

    // Create the numbering of unknowns
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);
    lset.CreateNumbering(         MG.GetLastLevel(), lidx);

    QuadOutCL *brickout=0;
    // If wished print out the geometry on a quadrilateral grid
    if (C.quad){
        if (ProcCL::IamMaster())
            std::cerr << "Write out geometry of the quadrilateral grid" << std::endl;
        brickout = new QuadOutCL( MG, lidx);
        brickout->Init( C.gridX, C.gridY, C.gridZ, C.stepsize, C.barycenter, C.rotation);
        brickout->putGeom(C.quadFileName + std::string(".geo"));
    }

    // References of Exchange classes
    ExchangeCL& ExV = Stokes.GetEx(Stokes.velocity);
    ExchangeCL& ExP = Stokes.GetEx(Stokes.pressure);
//     ExchangeBlockCL& Ex = Stokes.GetEx();
    ExchangeCL& ExL = lset.GetEx();

    // Get information about problem size and out them onto std::cerr
    DisplayUnks(vidx, pidx, lidx, ExV, ExP, ExL, MG);

    // Init ensight-output class
    EnsightP2SolOutCL ensight( MG, lidx, C.binary, C.masterOut);
    const string filename= C.ensDir + "/" + C.ensCase;
    const string datgeo= filename+".geo",
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.ensCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr,  true);
    ensight.DescribeVector( "Velocity", datvec, true);
    ensight.putGeom( datgeo);


    // Tell matrices and vectors about the numbering
    lset.Phi.SetIdx( lidx);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    cplN.SetIdx( vidx);
    Stokes.A.SetIdx(vidx, vidx);
    Stokes.B.SetIdx(pidx, vidx);
    Stokes.M.SetIdx(vidx, vidx);
    Stokes.N.SetIdx(vidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx);
    Stokes.prA.SetIdx( pidx, pidx);

    // Init velocity with zero and setup pressure matrices
    Stokes.InitVel( &Stokes.v, Null);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    // Compute intitial state
    switch (C.IniCond)
    {
      // stationary flow with/without drop
      case 1: case 2:
      {
//         typedef ParDummyPcCL SPcPcT;
//         SPcPcT SPcPc;
//
//         typedef ParPCGSolverCL<SPcPcT, ExchangeCL> SPcSolverT;
//         SPcSolverT CGsolver( 50, 0.02, ExP, SPcPc, /*relative*/ true, /*accur*/ true);
//
//         typedef ISNonlinearPreCL<SPcSolverT> SPcT;
//         SPcT ispc( CGsolver, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 0., /*kM*/ 1.);

        // PC for Schur-Complement
        typedef ParDummyPcCL<ExchangeCL> SPcT;
        SPcT ispc(Stokes.GetEx(Stokes.pressure));

        // PC for A-Block-PC
        typedef ParJac0CL<ExchangeCL>  APcPcT;
        APcPcT Apcpc(Stokes.GetEx(Stokes.velocity));

        // PC for A-block
        typedef ParPCGSolverCL<APcPcT,ExchangeCL> ASolverT;        // CG-based APcT
        ASolverT Asolver( 500, 0.02, ExV, Apcpc, /*relative=*/ true, /*accur*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver/*, &std::cerr*/);

        // Oseen solver
        typedef ParInexactUzawaCL<APcT, SPcT, APC_SYM, ExchangeCL, ExchangeCL> OseenSolverT;
        OseenSolverT schurSolver( Apc, ispc, ExV, ExP, C.outer_iter, C.outer_tol, 0.1);

        const Point3DCL old_Radius= C.Radius;
        if (C.IniCond==2) // stationary flow without drop
            for (int i=0; i<3; ++i)
                C.Radius[i]= -10.;
        lset.Init( DistanceFct);

        // Setup initial problem
        VelVecDescCL curv( vidx);
        if (ProcCL::IamMaster())
            std::cerr << "=================================================================================== Init:\n"
                      << "==> Initialize Problem\n";
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        curv.Data=0.;
        lset.AccumulateBndIntegral( curv);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cerr << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";

        // Solve initial problem
        double theta= C.theta;
        time.Reset();
        int step=0;
        int iters=0;
        do
        {
            Stokes.SetupNonlinear( &Stokes.N, &Stokes.v, &cplN, lset, Stokes.t);
            cplN.Data-= (1-theta) * (Stokes.N.Data * Stokes.v.Data);
            schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
                               Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
            if (ProcCL::IamMaster())
                std::cerr << "- Solving lin. Stokes ("<<step<<"): iter "<<schurSolver.GetIter()
                          <<", resid "<<schurSolver.GetResid()<<std::endl;
            ++step; iters+= schurSolver.GetIter();
        } while (schurSolver.GetIter() > 0);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cerr << "- Solving Stokes for initialization took "<<duration<<" sec, "
                      << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<schurSolver.GetResid()<<'\n';

        if (C.IniCond==2) // stationary flow without drop
        {
            C.Radius= old_Radius;
            lset.Init( DistanceFct);
        }
      } break;

      // read from file
      case 3:
      {
        ReadEnsightP2SolCL reader( MG, false);
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel, ExV);
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr,  ExP);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData()     ,  ExL);
      } break;

      // Use distance function
      default:
        lset.Init( DistanceFct);
    }

    const double Vol= 4./3.*M_PI*C.Radius[0]*C.Radius[1]*C.Radius[2];
    double relVol = lset.GetVolume()/Vol;
    if (ProcCL::IamMaster())
        std::cerr << "- Relative Volume is " << relVol << std::endl;

    // Write solution out in ensight format
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    // Use fractional step for solving coupled Navier-Stokes problem
    if (C.scheme)
    {
        // Linear equation system solver
        // Preconditioner
        // PC for instat. Schur complement
//         typedef ParAccPcCL SPcT;
//         SPcT ispc(ExP);
        typedef ParDummyPcCL<ExchangeCL>  SPcT;
        SPcT ispc(Stokes.GetEx(Stokes.pressure));

        // PC for A-Block-PC
        typedef ParJac0CL<ExchangeCL>  APcPcT;
        APcPcT Apcpc(Stokes.GetEx(Stokes.velocity));

        // PC for A-block
        typedef ParPreGMResSolverCL<APcPcT, ExchangeCL>    ASolverT;        // GMRes-based APcT
        ASolverT Asolver( /*restart*/    100, /*iter*/ 500,  /*tol*/   2e-5, ExV, Apcpc,
                          /*relative=*/ true, /*acc*/ true, /*modGS*/ false, RightPreconditioning, /*mod4par*/true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver/*,&std::cerr*/);

        // PC for Oseen solver
//         typedef DiagBlockPreCL<APcT, SPcT> OseenPcT;
//         OseenPcT oseenpc( Apc, ispc);

        // Oseen Solver
//         typedef ParPreGCRSolverCL<OseenPcT,ExchangeBlockCL> OseenBaseSolverT;
//         OseenBaseSolverT oseensolver0( /*trunc*/C.outer_iter, C.outer_iter, C.outer_tol, Ex, oseenpc,
//                                        /*mod*/false, /*rel*/ false, /*acc*/true);
//         typedef ParPreGMResSolverCL<OseenPcT,ExchangeBlockCL> OseenBaseSolverT;
//         OseenBaseSolverT oseensolver0( /*trunc*/C.outer_iter, C.outer_iter, C.outer_tol, Ex, oseenpc,
//                                        /*rel*/ false, /*acc*/true, /*modGS*/false,  LeftPreconditioning, /*mod4par*/true);

//         typedef BlockMatrixSolverCL<OseenBaseSolverT> OseenSolverT;
//         OseenSolverT oseensolver( oseensolver0);

        typedef ParInexactUzawaCL<APcT, SPcT, APC_OTHER, ExchangeCL, ExchangeCL> OseenSolverT;
        OseenSolverT oseensolver( Apc, ispc, ExV, ExP, C.outer_iter, C.outer_tol, 0.2, 500, &std::cerr);

        typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
        NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red);

        time.Reset();
        // Coupling Navier-Stokes with Levelset
//         typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//         CouplingT cpl( Stokes, lset, nssolver, C.nonlinear);

        typedef RecThetaScheme2PhaseCL <StokesProblemT, NSSolverT> CouplingT;
        CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear);

        time.Stop();
        duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cerr << "- Updating discretization took "<<duration<<" sec.\n";

        // Set time step and create matrices
        cpl.SetTimeStep( C.dt);

        for (int step= 1; step<=C.num_steps; ++step)
        {
            ParTimerCL step_time;
            step_time.Reset();
            if (ProcCL::IamMaster())
                std::cerr << "=================================================================================== Schritt " << step << ":\n"
                          << "==> Solving coupled Levelset-Navier-Stokes problem ....\n";
            time.Reset();
            cpl.DoStep( C.cpl_iter);
            time.Stop(); duration=time.GetMaxTime();
            if (ProcCL::IamMaster())
                std::cerr << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec."<<std::endl;
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster())
                std::cerr << "- rel. Volume: " << relVol << std::endl;
            if (C.VolCorr)
            {
                if (ProcCL::IamMaster())
                    std::cerr << "\n==> Adjust volume ...\n";
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                lset.Phi.Data+= dphi;
                relVol = lset.GetVolume()/Vol;
                if (ProcCL::IamMaster())
                    std::cerr << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< std::endl;
            }
            if ((C.ensight && step%C.ensight==0) || step==C.num_steps)
            {
                step_time.Stop();
                ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
                ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
                ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
                ensight.Commit();
                step_time.Start();
            }

            // Reparametrization of levelset function
            if (C.RepFreq && step%C.RepFreq==0)
            {
                if (ProcCL::IamMaster())
                    std::cerr << "\n==> Reparametrization with FastMarching algorithm"<<std::endl;
                time.Reset();
                lset.ReparamFastMarching( ExL, C.RepMethod, C.RepMethod==3);
                time.Stop(); duration=time.GetMaxTime();
                relVol = lset.GetVolume()/Vol;
                if (ProcCL::IamMaster())
                    std::cerr << "- FastMarching took "<<duration<<" sec."<<std::endl;

                if (ProcCL::IamMaster())
                    std::cerr << "- rel. Volume: " << relVol << std::endl;
                if (C.VolCorr)
                {
                    if (ProcCL::IamMaster())
                        std::cerr << "\n==> Adjust volume ...\n";

                    double dphi= lset.AdjustVolume( Vol, 1e-9);
                    lset.Phi.Data+= dphi;
                    relVol = lset.GetVolume()/Vol;
                    if (ProcCL::IamMaster())
                        std::cerr << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< std::endl;
                }
                step_time.Stop();
                // Write solution out after reparametrization
                if ((C.ensight && step%C.ensight==0) || step==C.num_steps)
                {
                    ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.1)*C.dt);
                    ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.1)*C.dt);
                    ensight.putScalar( datscl, lset.GetSolution(), (step+0.1)*C.dt);
                    ensight.Commit();
                }
                step_time.Start();
            }
            step_time.Stop(); duration=step_time.GetMaxTime();
            if (ProcCL::IamMaster())
			{
                std::cerr <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
			}
        }
        //ExP.CheckSends();
        //ExV.CheckSends();
    }
    else {} // No other method yet implemented
    ensight.CaseEnd();

    double min= GlobalMin(Stokes.p.Data.min()),
           max= GlobalMax(Stokes.p.Data.max());
    if (ProcCL::IamMaster())
        std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    if (C.quad){
        LevelsetP2CL::const_DiscSolCL lset_sol=lset.GetSolution();
        if (ProcCL::IamMaster())
            std::cerr << "Writing out velocity on quadrilateral grid"<<std::endl;
        brickout->putVector(C.quadFileName + std::string(".vel_norm"),
                            C.quadFileName + std::string(".velY"),
                            C.quadFileName + std::string(".velZ"),
                            Stokes.GetVelSolution(), &lset_sol);
        delete brickout;
    }
}
} // end of namespace DROPS


// **************************************************************
// * M A R K  D R O P                                           *
// **************************************************************
// * Mark tetras for refining in the near of phase-border       *
// **************************************************************
void MarkDrop (DROPS::MultiGridCL& mg, int maxLevel= -1)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( DistanceFct(GetBaryCenter(*It)) <= C.ref_width)
            It->SetRegRefMark();
    }
}


// **************************************************************
// * M A I N                                                    *
// **************************************************************
int main (int argc, char** argv)
{
  // Init parallel enviroment before using try, so error handling works correct
  DROPS::ProcCL Proc(&argc, &argv);
  try
  {
    DROPS::ParTimerCL alltime, time;
    double duration;
    if (DROPS::ProcCL::IamMaster())
        std::cerr << "TestMzellePar: Running on "<<DROPS::ProcCL::Size()<<" processors!\n";

    // Read parameter
    if (argc!=2)
    {
        std::cerr << "You have to specify one parameter:\n\t"
                  << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();

    if (DROPS::ProcCL::IamMaster())
        std::cerr << C << std::endl;

    // Typedefinition of the problem
    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyNavStokesCL;

    // Create geometry on proc 0
    DROPS::ParMultiGridCL pmg;
    std::ifstream meshfile( C.meshfile.c_str());
    if (!meshfile)
    {
        std::cerr << "error while opening mesh file " << C.meshfile << "\n";
        return 1;
    }
    DROPS::ReadMeshBuilderCL * mgb;
    if (DROPS::ProcCL::IamMaster())
        mgb = new DROPS::ReadMeshBuilderCL( meshfile );
    else
        mgb = new DROPS::EmptyReadMeshBuilderCL( meshfile );
    DROPS::MultiGridCL mg( *mgb );
    pmg.AttachTo(mg);

    // Init load balancing class
    DROPS::LoadBalHandlerCL lb(mg, C.quality);
    switch (C.refineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
    }

    // Create boundary-information
    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();
    if (num_bnd>10) {
        std::cerr << "Increase size of BndSegs in main() for proper use!\n";
        return 1;
    }
    DROPS::BndCondT bc[10];
    DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[10];
    for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
    {
        bnd_fun[i]= (bc[i]= mgb->GetBC( i))==DROPS::DirBC ? &Inflow : &Null;
        if (DROPS::ProcCL::IamMaster())
            std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }

    // Setup problem
    MyNavStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun));

    time.Reset();
    // Distribute and refine multigrid
    lb.DoInitDistribution(DROPS::ProcCL::Master());
    for (int i=0; i<C.ref_flevel; ++i)
    {
        if (DROPS::ProcCL::IamMaster())
            std::cerr << "+ Refine drop "<<i<<std::endl;
        MarkDrop( mg);
        mg.Refine();
    }
    lb.DoMigration();
    time.Stop(); duration=time.GetMaxTime();
    if (DROPS::ProcCL::IamMaster())
        std::cerr << " Creating and distributing of multigrid took "<<duration<<" sec."<<std::endl;

    // Check the parallel multigrid
    CheckParMultiGrid(pmg);
    if (DROPS::ProcCL::IamMaster())
        std::cerr << " Number of simplices over procs:\n";
    mg.SizeInfo(std::cerr);
    double tetra_ratio=lb.GetTetraBalance();
    if (DROPS::ProcCL::IamMaster())
        std::cerr << " Maximal number of tetras over minimal number of tetras: "<<tetra_ratio<<std::endl;

    DROPS::Uint *numTetrasAllProc=0;
    DROPS::Uint *numFacesAllProc=0;
    DROPS::Uint *numDistFaceAllProc=0;
    if (DROPS::ProcCL::IamMaster()){
        numTetrasAllProc  = new DROPS::Uint[DROPS::ProcCL::Size()];
        numFacesAllProc   = new DROPS::Uint[DROPS::ProcCL::Size()];
        numDistFaceAllProc= new DROPS::Uint[DROPS::ProcCL::Size()];
    }
    DROPS::Gather(mg.GetNumTriangTetra(),      numTetrasAllProc,   DROPS::ProcCL::Master());
    DROPS::Gather(mg.GetNumTriangFace(),       numFacesAllProc,    DROPS::ProcCL::Master());
    DROPS::Gather(mg.GetNumDistributedFaces(), numDistFaceAllProc, DROPS::ProcCL::Master());

    if (DROPS::ProcCL::IamMaster()){
        double ratioTetra       =  (double)*std::max_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size())
                                  /(double)*std::min_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size());
        DROPS::Uint allTetra    =  std::accumulate(numTetrasAllProc, numTetrasAllProc+DROPS::ProcCL::Size(), 0),
                    allFace     =  std::accumulate(numFacesAllProc, numFacesAllProc+DROPS::ProcCL::Size(), 0),
                    allDistFace =  std::accumulate(numDistFaceAllProc, numDistFaceAllProc+DROPS::ProcCL::Size(), 0);
        double      *ratioDistFace=new double[DROPS::ProcCL::Size()];
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            ratioDistFace[i]= ((double)numDistFaceAllProc[i]/(double)numFacesAllProc[i]*100.);
        double maxRatio= *std::max_element(ratioDistFace, ratioDistFace+DROPS::ProcCL::Size());
        std::cerr << "#(master tetras in finest level): "<<allTetra<<", #(distributed Faces): "<<allDistFace<<", #(all Faces): "<<allFace<<"\n"
                  << "Proc\tTetra\tFace\t(%Distributed Faces)\n";
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            std::cerr << i << '\t'<< numTetrasAllProc[i]<<'\t'<<numDistFaceAllProc[i]
                      <<'\t'<<ratioDistFace[i]<<std::endl;
        std::cerr << "Ratio between max/min Tetra: "<<ratioTetra<<" max Ratio DistFace/AllFace: "<<maxRatio<<std::endl;
        delete[] numTetrasAllProc;
        delete[] numFacesAllProc;
        delete[] numDistFaceAllProc;
        delete[] ratioDistFace;
    }
    // Solve the coupled Navier-Stokes equation
    Strategy( prob);
    alltime.Stop(); duration=alltime.GetMaxTime();
    if (DROPS::ProcCL::IamMaster())
        std::cerr << " The whole programm took "<<duration<<" sec."<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }   // error handling
}
