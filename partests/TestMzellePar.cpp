/// \file TestMzellePar.cpp
/// \brief Testing parallel solvers for the two-phase flow in the measurement cell
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
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include "parallel/exchange.h"

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
#include "num/stokessolverfactory.h"
#include "levelset/surfacetension.h"
#include "levelset/surfacetension.h"

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
    if( DROPS::ProcCL::Check(pmg_sane && mg_sane) ){
        IF_MASTER
                std::cout << " As far as I can tell, the multigrid is sane\n";
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
using DROPS::ProcCL;
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
    Ulint Psize_acc = ProcCL::GlobalSum(Psize);
    Ulint Vsize_acc = ProcCL::GlobalSum(Vsize);
    Ulint Lsize_acc = ProcCL::GlobalSum(Lsize);

    // maximal and minimal number of unknowns
    Ulint P_min= ProcCL::GlobalMin(Psize); Ulint P_max= ProcCL::GlobalMax(Psize);
    Ulint V_min= ProcCL::GlobalMin(Vsize); Ulint V_max= ProcCL::GlobalMax(Vsize);
    Ulint L_min= ProcCL::GlobalMin(Lsize); Ulint L_max= ProcCL::GlobalMax(Lsize);

    // ratios between maximal number of unknowns/proc and minimal number
    double P_ratio   = (double)P_max/(double)P_min;
    double V_ratio   = (double)V_max/(double)V_min;
    double L_ratio   = (double)L_max/(double)L_min;

    // number on boundaries
    Ulint P_accmax= ProcCL::GlobalMax(ExP.AccDistIndex.size()), P_accmin= ProcCL::GlobalMin(ExP.AccDistIndex.size());
    Ulint V_accmax= ProcCL::GlobalMax(ExV.AccDistIndex.size()), V_accmin= ProcCL::GlobalMin(ExV.AccDistIndex.size());
    Ulint L_accmax= ProcCL::GlobalMax(ExL.AccDistIndex.size()), L_accmin= ProcCL::GlobalMin(ExL.AccDistIndex.size());

    // ratio of these unknowns
    double P_accratio= (double)P_accmax / (double)P_accmin;
    double V_accratio= (double)V_accmax / (double)V_accmin;
    double L_accratio= (double)L_accmax / (double)L_accmin;

    // output on screen
    if (DROPS::ProcCL::IamMaster()){
        std::cout << "  + Number of DOF\n        "
                  << std::setw(10)<<"global"<<std::setw(10)<<"accum"<<std::setw(10)
                  << "max"<<std::setw(10)<<"min"<<std::setw(10)<<"ratio"<<"  |  "
                  << std::setw(10)<<"max_acc" <<std::setw(10)<<"min_acc"<<std::setw(10)<<"ratio_acc"<<std::endl;

        std::cout << "    "<<"pr  "
                  << std::setw(10)<<GPsize<<std::setw(10)<<Psize_acc<<std::setw(10)<<P_max
                  << std::setw(10)<<P_min<< std::setw(10)<<P_ratio<<"  |  "
                  << std::setw(10)<<P_accmax<<std::setw(10)<<P_accmin<<std::setw(10)<<P_accratio<<std::endl;

        std::cout << "    "<<"vel "
                  << std::setw(10)<<GVsize<<std::setw(10)<<Vsize_acc<<std::setw(10)<<V_max
                  << std::setw(10)<<V_min<< std::setw(10)<<V_ratio<<"  |  "
                  << std::setw(10)<<V_accmax<<std::setw(10)<<V_accmin<<std::setw(10)<<V_accratio<<std::endl;

        std::cout << "    "<<"scl "
                  << std::setw(10)<<GLsize<<std::setw(10)<<Lsize_acc<<std::setw(10)<<L_max
                  << std::setw(10)<<L_min<< std::setw(10)<<L_ratio<<"  |  "
                  << std::setw(10)<<L_accmax<<std::setw(10)<<L_accmin<<std::setw(10)<<L_accratio<<std::endl;

        std::cout << std::endl;
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

    ZeroFlowCL( const DROPS::ParamMesszelleNsCL& C)
      : rho( DROPS::JumpCL( C.mat_DensDrop, C.mat_DensFluid ), DROPS::H_sm, C.mat_SmoothZone),
        mu(  DROPS::JumpCL( C.mat_ViscDrop,  C.mat_ViscFluid),   DROPS::H_sm, C.mat_SmoothZone),
        SurfTens( C.sft_SurfTension), g( C.exp_Gravity)    {}
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

    DimLessCoeffCL( const DROPS::ParamMesszelleNsCL& C)
      : rho( DROPS::JumpCL( 1., C.mat_DensFluid/C.mat_DensDrop ), DROPS::H_sm, C.mat_SmoothZone),
        mu ( DROPS::JumpCL( 1., C.mat_ViscFluid/C.mat_ViscDrop),    DROPS::H_sm, C.mat_SmoothZone),
        SurfTens( C.sft_SurfTension/C.mat_DensDrop), g( C.exp_Gravity)    {}
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
    const double s2= C.exp_RadInlet*C.exp_RadInlet,
                 r2= p.norm_sq() - p[C.exp_FlowDir]*p[C.exp_FlowDir];
    ret[C.exp_FlowDir]= -(r2-s2)/s2*C.exp_InflowVel;
    return ret;
}

// **************************************************************
// * D I S T A N C E  F C T                                     *
// **************************************************************
// * distance function for initial droplet                      *
// **************************************************************
double DistanceFct( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.exp_PosDrop-p;
    const double avgRad = cbrt(C.exp_RadDrop[0]*C.exp_RadDrop[1]* C.exp_RadDrop[2]);
    d/= C.exp_RadDrop/avgRad;
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
    SurfaceTensionCL sf( sigmaf, 0);
    LevelsetP2CL lset( MG, sf, C.lvs_SD, C.lvs_CurvDiff, C.rpm_NarrowBand);

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
    if (C.qlg_Quad){
        if (ProcCL::IamMaster())
            std::cout << "Write out geometry of the quadrilateral grid" << std::endl;
        brickout = new QuadOutCL( MG, lidx);
        brickout->Init( C.qlg_GridX, C.qlg_GridY, C.qlg_GridZ, C.qlg_Stepsize, C.qlg_Barycenter, C.qlg_Rotation);
        brickout->putGeom(C.qlg_FileName + std::string(".geo"));
    }

    // References of Exchange classes
    ExchangeCL& ExV = Stokes.vel_idx.GetEx();
    ExchangeCL& ExP = Stokes.pr_idx.GetEx();
//     ExchangeBlockCL& Ex = Stokes.GetEx();
    ExchangeCL& ExL = lset.idx.GetEx();

    // Get information about problem size and out them onto std::cout
    DisplayUnks(vidx, pidx, lidx, ExV, ExP, ExL, MG);

    // Init ensight-output class
    std::string ensf( C.ens_EnsDir + "/" + C.ens_EnsCase);
    Ensight6OutCL ensight( C.ens_EnsCase + ".case", C.tm_NumSteps+1);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   "Messzelle",     ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));

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
    switch (C.dmc_InitialCond)
    {
      // stationary flow with/without drop
      case 1: case 2:
      {
        StokesSolverParamST statStokesParam(C);
        statStokesParam.stk_StokesMethod= 20301;
        statStokesParam.tm_NumSteps   = 0;
        statStokesParam.tm_StepSize   = 0.;
        statStokesParam.stk_PcATol     = 0.02;
        StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> statStokesSolverFactory(Stokes, statStokesParam);
        StokesSolverBaseCL* statStokesSolver= statStokesSolverFactory.CreateStokesSolver();

        const Point3DCL old_Radius= C.exp_RadDrop;
        if (C.dmc_InitialCond==2) // stationary flow without drop
            for (int i=0; i<3; ++i)
                C.exp_RadDrop[i]= -10.;
        lset.Init( DistanceFct);

        // Setup initial problem
        VelVecDescCL curv( vidx);
        if (ProcCL::IamMaster())
            std::cout << "=================================================================================== Init:\n"
                      << "==> Initialize Problem\n";
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        curv.Data=0.;
        lset.AccumulateBndIntegral( curv);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cout << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";

        // Solve initial problem
        double theta= C.stk_Theta;
        time.Reset();
        int step=0;
        int iters=0;
        do
        {
            Stokes.SetupNonlinear( &Stokes.N, &Stokes.v, &cplN, lset, Stokes.t);
            cplN.Data-= (1-theta) * (Stokes.N.Data * Stokes.v.Data);
            statStokesSolver->Solve( Stokes.A.Data, Stokes.B.Data,
                               Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
            if (ProcCL::IamMaster())
                std::cout << "- Solving lin. Stokes ("<<step<<"): iter "<<statStokesSolver->GetIter()
                          <<", resid "<<statStokesSolver->GetResid()<<std::endl;
            ++step; iters+= statStokesSolver->GetIter();
        } while (statStokesSolver->GetIter() > 0);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cout << "- Solving Stokes for initialization took "<<duration<<" sec, "
                      << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<statStokesSolver->GetResid()<<'\n';

        if (C.dmc_InitialCond==2) // stationary flow without drop
        {
            C.exp_RadDrop= old_Radius;
            lset.Init( DistanceFct);
        }
        delete statStokesSolver;
      } break;

      // read from file
      case 3:
      {
        ReadEnsightP2SolCL reader( MG, false);
        reader.ReadVector( C.dmc_InitialFile+".vel", Stokes.v, Stokes.GetBndData().Vel);
        reader.ReadScalar( C.dmc_InitialFile+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        reader.ReadScalar( C.dmc_InitialFile+".scl", lset.Phi, lset.GetBndData());
      } break;

      // Use distance function
      default:
        lset.Init( DistanceFct);
    }

    const double Vol= 4./3.*M_PI*C.exp_RadDrop[0]*C.exp_RadDrop[1]*C.exp_RadDrop[2];
    double relVol = lset.GetVolume()/Vol;
    if (ProcCL::IamMaster())
        std::cout << "- Relative Volume is " << relVol << std::endl;

    // Write solution out in ensight format
    if (C.ens_EnsightOut) ensight.Write( 0.);

    // Use fractional step for solving coupled Navier-Stokes problem
    if (C.tm_Scheme)
    {

        StokesSolverParamST instatStokesParam(C);
        StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> instatStokesSolverFactory(Stokes, instatStokesParam);
        StokesSolverBaseCL* instatStokesSolver= instatStokesSolverFactory.CreateStokesSolver();


        typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
        NSSolverT nssolver( Stokes, *instatStokesSolver, C.ns_Iter, C.ns_Tol, C.ns_Reduction);

        time.Reset();
        // Coupling Navier-Stokes with Levelset
//         typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//         CouplingT cpl( Stokes, lset, nssolver, C.nonlinear);
        typedef ParPreGMResSolverCL<ParJac0CL> LsetSolverT;
        ParJac0CL jacparpc( *lidx);
        LsetSolverT gm(/*restart*/100, C.lvs_Iter, C.lvs_Tol, *lidx, jacparpc,/*rel*/true, /*acc*/ true, /*modGS*/false, LeftPreconditioning, /*parmod*/true);

        LevelsetModifyCL lsetmod( C.rpm_Freq, C.rpm_Method, 0, 1, C.lvs_VolCorrection, Vol);
        typedef RecThetaScheme2PhaseCL <StokesProblemT, LsetSolverT> CouplingT;
        CouplingT cpl( Stokes, lset, nssolver, gm, lsetmod, C.cpl_Tol, C.stk_Theta, C.lvs_Theta, C.ns_Nonlinear);

        time.Stop();
        duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cout << "- Updating discretization took "<<duration<<" sec.\n";

        // Set time step and create matrices
        cpl.SetTimeStep( C.tm_StepSize);

        for (int step= 1; step<=C.tm_NumSteps; ++step)
        {
            ParTimerCL step_time;
            step_time.Reset();
            if (ProcCL::IamMaster())
                std::cout << "=================================================================================== Schritt " << step << ":\n"
                          << "==> Solving coupled Levelset-Navier-Stokes problem ....\n";
            time.Reset();
            cpl.DoStep( C.cpl_Iter);
            time.Stop(); duration=time.GetMaxTime();
            if (ProcCL::IamMaster())
                std::cout << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec."<<std::endl;

            if ((C.ens_EnsightOut && step%C.ens_EnsightOut==0) || step==C.tm_NumSteps)
            {
                step_time.Stop();
                ensight.Write( step*C.tm_StepSize);
                step_time.Start();
            }

            step_time.Stop(); duration=step_time.GetMaxTime();
            if (ProcCL::IamMaster()){
                std::cout <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
            }
        }
        delete instatStokesSolver;
    }
    else {} // No other method yet implemented

    double min= ProcCL::GlobalMin(Stokes.p.Data.min()),
           max= ProcCL::GlobalMax(Stokes.p.Data.max());
    if (ProcCL::IamMaster())
        std::cout << "pressure min/max: "<<min<<", "<<max<<std::endl;

    if (C.qlg_Quad){
        LevelsetP2CL::const_DiscSolCL lset_sol=lset.GetSolution();
        if (ProcCL::IamMaster())
            std::cout << "Writing out velocity on quadrilateral grid"<<std::endl;
        brickout->putVector(C.qlg_FileName + std::string(".vel_norm"),
                            C.qlg_FileName + std::string(".velY"),
                            C.qlg_FileName + std::string(".velZ"),
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
        if ( DistanceFct(GetBaryCenter(*It)) <= C.ref_Width)
            It->SetRegRefMark();
    }
}


// **************************************************************
// * M A I N                                                    *
// **************************************************************
int main (int argc, char** argv)
{
  // Init parallel enviroment before using try, so error handling works correct
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
  try
  {
    DROPS::ParTimerCL alltime, time;
    double duration;
    if (DROPS::ProcCL::IamMaster())
        std::cout << "TestMzellePar: Running on "<<DROPS::ProcCL::Size()<<" processors!\n";

    // Read parameter
    if (argc!=2)
    {
        std::cout << "You have to specify one parameter:\n\t"
                  << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();

    if (DROPS::ProcCL::IamMaster())
        std::cout << C << std::endl;

    // Typedefinition of the problem
    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyNavStokesCL;

    // Create geometry on proc 0
    DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
    std::ifstream meshfile( C.dmc_MeshFile.c_str());
    if (!meshfile)
    {
        std::cout << "error while opening mesh file " << C.dmc_MeshFile << "\n";
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
    switch (C.ref_RefineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
    }

    // Create boundary-information
    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();
    if (num_bnd>10) {
        std::cout << "Increase size of BndSegs in main() for proper use!\n";
        return 1;
    }
    DROPS::BndCondT bc[10];
    DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[10];
    for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
    {
        bnd_fun[i]= (bc[i]= mgb->GetBC( i))==DROPS::DirBC ? &Inflow : &Null;
        if (DROPS::ProcCL::IamMaster())
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
    }

    // Setup problem
    MyNavStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun));

    time.Reset();
    // Distribute and refine multigrid
    lb.DoInitDistribution(DROPS::ProcCL::Master());
    for (int i=0; i<C.ref_FinestLevel; ++i)
    {
        if (DROPS::ProcCL::IamMaster())
            std::cout << "+ Refine drop "<<i<<std::endl;
        MarkDrop( mg);
        mg.Refine();
    }
    lb.DoMigration();
    time.Stop(); duration=time.GetMaxTime();
    if (DROPS::ProcCL::IamMaster())
        std::cout << " Creating and distributing of multigrid took "<<duration<<" sec."<<std::endl;

    // Check the parallel multigrid
    CheckParMultiGrid(pmg);
    if (DROPS::ProcCL::IamMaster())
        std::cout << " Number of simplices over procs:\n";
    mg.SizeInfo(std::cout);
    double tetra_ratio=lb.GetTetraBalance();
    if (DROPS::ProcCL::IamMaster())
        std::cout << " Maximal number of tetras over minimal number of tetras: "<<tetra_ratio<<std::endl;

    DROPS::Uint *numTetrasAllProc=0;
    DROPS::Uint *numFacesAllProc=0;
    DROPS::Uint *numDistFaceAllProc=0;
    if (DROPS::ProcCL::IamMaster()){
        numTetrasAllProc  = new DROPS::Uint[DROPS::ProcCL::Size()];
        numFacesAllProc   = new DROPS::Uint[DROPS::ProcCL::Size()];
        numDistFaceAllProc= new DROPS::Uint[DROPS::ProcCL::Size()];
    }
    DROPS::ProcCL::Gather(mg.GetNumTriangTetra(),      numTetrasAllProc,   DROPS::ProcCL::Master());
    DROPS::ProcCL::Gather(mg.GetNumTriangFace(),       numFacesAllProc,    DROPS::ProcCL::Master());
    DROPS::ProcCL::Gather(mg.GetNumDistributedFaces(), numDistFaceAllProc, DROPS::ProcCL::Master());

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
        std::cout << "#(master tetras in finest level): "<<allTetra<<", #(distributed Faces): "<<allDistFace<<", #(all Faces): "<<allFace<<"\n"
                  << "Proc\tTetra\tFace\t(%Distributed Faces)\n";
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            std::cout << i << '\t'<< numTetrasAllProc[i]<<'\t'<<numDistFaceAllProc[i]
                      <<'\t'<<ratioDistFace[i]<<std::endl;
        std::cout << "Ratio between max/min Tetra: "<<ratioTetra<<" max Ratio DistFace/AllFace: "<<maxRatio<<std::endl;
        delete[] numTetrasAllProc;
        delete[] numFacesAllProc;
        delete[] numDistFaceAllProc;
        delete[] ratioDistFace;
    }
    // Solve the coupled Navier-Stokes equation
    Strategy( prob);
    alltime.Stop(); duration=alltime.GetMaxTime();
    if (DROPS::ProcCL::IamMaster())
        std::cout << " The whole programm took "<<duration<<" sec."<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }   // error handling
}

