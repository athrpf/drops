//**************************************************************************
// File:    MzelleNMRParamEst.cpp                                          *
// Content: parallel solver                                                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file MzelleNMRParamEst.cpp
/// \brief Droplet in the meassurement device

// include std header for two-phase flows
#include "partests/two_phase_hdr.h"
#include "levelset/mzelle_hdr.h"

// include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include "parallel/exchange.h"
#include "parallel/logger.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "num/parsolver.h"
#include "num/parprecond.h"
#include "num/stokessolver.h"
#include "num/parstokessolver.h"
#include "num/stokessolverfactory.h"
#include "stokes/integrTime.h"
#include "levelset/adaptriang.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "out/gridOut.h"

 // include problem class
#include "navstokes/instatnavstokes2phase.h"
#include "levelset/coupling.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#endif


enum TimePart{
    T_solve_init_Stokes,
    T_solve_NS,
    T_transform_to_quad
};

DROPS::TimeStoreCL Times(3);

void SetDescriber()
{
    Times.SetDescriber(T_solve_init_Stokes, "Solving initial Stokes problem");
    Times.SetDescriber(T_solve_NS,          "Time-integration of Navier-Stokes");
    Times.SetDescriber(T_transform_to_quad, "Transforming solution to quadrilateral grid");
    Times.SetCounterDescriber("Moved MultiNodes");
}

DROPS::ParamParBrickFlowCL C;

namespace DROPS{

/// \brief Inflow
DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.exp_RadInlet*C.exp_RadInlet,
                 r2= p.norm_sq() - p[C.exp_FlowDir]*p[C.exp_FlowDir];
    ret[C.exp_FlowDir]= -(r2-s2)/s2*C.exp_InflowVel;
    return ret;
}

/// \brief Get y-distance to midpoint of the meassurements ...
double DistanceToBarycenterOfQuadGrid( const DROPS::Point3DCL& p)
{
    return p[1]-C.qlg_Barycenter[1];
}

/// \brief Get y-distance to midpoint of the droplet ...
double DistanceYToBarycenterOfDrop( const DROPS::Point3DCL& p)
{
    return p[1]-C.exp_PosDrop[1];
}


template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    ParTimerCL time;
    ParTimerCL initTimer;
    double duration;
    // droplet-information

    // Initial velocity is zero
    Stokes.InitVel( &Stokes.v, Null);
    lset.Init( EllipsoidCL::DistanceFct);

    IFInfo.Update(lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.t);

    switch (C.dmc_InitialCond)
    {
      // stationary flow with/without drop
      case 1: case 2:
      {
        StokesSolverParamST statStokesParam(C);
        statStokesParam.stk_StokesMethod= 20301;
        statStokesParam.tm_NumSteps   = 0;
        statStokesParam.tm_StepSize          = 0.;
        statStokesParam.stk_PcATol     = 0.02;
        StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> statStokesSolverFactory(Stokes, statStokesParam);
        StokesSolverBaseCL* statStokesSolver= statStokesSolverFactory.CreateStokesSolver();

        VelVecDescCL curv(&Stokes.vel_idx),
                    cplN(&Stokes.vel_idx);

        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        Stokes.SetupPrMass( &Stokes.prM, lset);
        curv.Clear();
        curv.Data=0.;
        lset.AccumulateBndIntegral( curv);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
        {
            std::cout << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";
            // Log measured duration.
            DROPS_LOGGER_SETVALUE("StokesCurvInit",duration);
        }

        //Solve initial problem
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
        {
            std::cout << "- Solving Stokes for initialization took "<<duration<<" sec, "
                      << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<statStokesSolver->GetResid()<<'\n';
            // Log measured duration
            DROPS_LOGGER_SETVALUE("SolStokesTime",duration);
            DROPS_LOGGER_SETVALUE("SolStokesIter",iters);
        }
        delete statStokesSolver;
      }break;
      case 3:
        {
            ReadEnsightP2SolCL reader( Stokes.GetMG(), false);
            reader.ReadVector( C.dmc_InitialFile+".vel", Stokes.v, Stokes.GetBndData().Vel);
            reader.ReadScalar( C.dmc_InitialFile+".pr",  Stokes.p, Stokes.GetBndData().Pr);
            reader.ReadScalar( C.dmc_InitialFile+".scl", lset.Phi, lset.GetBndData());
            if (ProcCL::IamMaster())
                std::cout << "- Initial Conditions successfull read\n";
        } break;
    }
    initTimer.Stop();
    Times.AddTime(T_solve_init_Stokes, initTimer.GetTime());
}

template<typename Coeff>
  void SolveCoupledNS(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, AdapTriangCL& adap)
{
    ParTimerCL time;
    ParTimerCL timeIntegrationTimer;
    double duration;
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    // writer in vtk-format
    typedef TwoPhaseVTKCL<StokesProblemT, LevelsetP2CL> VTKWriterT;
    VTKWriterT vtkwriter(adap.GetMG(), Stokes, lset, (C.vtk_VTKOut ? C.tm_NumSteps/C.vtk_VTKOut+1 : 0),
                         std::string(C.vtk_VTKDir + "/" + C.vtk_VTKName), C.vtk_Binary);
    // writer for ensight format
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.ens_EnsDir, C.ens_EnsCase, C.ens_GeomName, /*adaptive=*/true,
                                  (C.ens_EnsightOut? C.tm_NumSteps/C.ens_EnsightOut+1 : 0), C.ens_Binary, C.ens_MasterOut);

    if (C.vtk_VTKOut)
        vtkwriter.write();
    if (C.ens_EnsightOut)
        ensightwriter.write();

    const double Vol= EllipsoidCL::GetVolume();

    StokesSolverParamST instatStokesParam(C);
    StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> instatStokesSolverFactory(Stokes, instatStokesParam);
    StokesSolverBaseCL* instatStokesSolver= instatStokesSolverFactory.CreateStokesSolver();
    StokesSolverBaseCL& oseensolver= *instatStokesSolver;

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
    NSSolverT nssolver( Stokes, oseensolver, C.ns_Iter, C.ns_Tol, C.ns_Reduction);

    // coupling levelset NavStokes
    time.Reset();
    // Coupling Navier-Stokes with Levelset
//    typedef LinThetaScheme2PhaseCL <StokesProblemT, NSSolverT> CouplingT;
//    CouplingT cpl( Stokes, lset, nssolver, C.nonlinear);

    typedef ParPreGMResSolverCL<ParJac0CL> LsetSolverT;
    ParJac0CL jacparpc( lset.idx);
    LsetSolverT gm(/*restart*/100, C.lvs_Iter, C.lvs_Tol, lset.idx, jacparpc,/*rel*/true, /*acc*/ true, /*modGS*/false, LeftPreconditioning, /*parmod*/true);
    LevelsetModifyCL lsetmod( C.rpm_Freq, C.rpm_Method, 0, 1, C.lvs_VolCorrection, Vol);

    typedef RecThetaScheme2PhaseCL <StokesProblemT, LsetSolverT> CouplingT;
    CouplingT cpl( Stokes, lset, nssolver, gm, lsetmod, C.cpl_Tol, C.stk_Theta, C.lvs_Theta, C.ns_Nonlinear);


    time.Stop();
    duration=time.GetMaxTime();
    if (ProcCL::IamMaster())
        std::cout << "- Updating discretization took "<<duration<<" sec.\n";

    // Set time step and create matrices
    cpl.SetTimeStep( C.tm_StepSize);

    int step= 1;
    for (; step<=C.tm_NumSteps; ++step)
    {
        ParTimerCL step_time;
        step_time.Reset();
        if (ProcCL::IamMaster())
            std::cout << "=================================================================================== step " << step << ":\n"
                      << " Idx for vel  "<<Stokes.v.RowIdx->GetIdx()
                      << "\n Idx for pr   "<<Stokes.p.RowIdx->GetIdx()
                      << "\n Idx for lset "<<lset.Phi.RowIdx->GetIdx()<<std::endl;
        time.Reset();

        if (C.ref_Freq && step%C.ref_Freq==0)
        {
            if (ProcCL::IamMaster())
                std::cout << "==> Adaptive Refinement of MultiGrid"<<std::endl;

            adap.UpdateTriang( lset);

            if (C.inf_CheckMG && !ProcCL::Check( CheckParMultiGrid(adap.GetPMG())) )
                throw DROPSErrCL("MultiGrid is incorrect!");

            if (adap.WasModified() )
            {
                cpl.Update();
                // don't forget to update the pr mass/stiff matrix for the schur compl. preconditioner!!
                Stokes.prM.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);
                Stokes.SetupPrMass( &Stokes.prM, lset);
                Stokes.prA.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);
                Stokes.SetupPrStiff( &Stokes.prA, lset);
            }
        }

        if (C.inf_PrintNumUnk)
            DisplayUnks(Stokes, lset, adap.GetMG());

        if (ProcCL::IamMaster())
            std::cout << "==> Solving coupled Levelset-Navier-Stokes problem ....\n";

        VectorCL u_old= Stokes.v.Data;

        cpl.DoStep( C.cpl_Iter);

        VectorCL e (Stokes.v.Data-u_old);
        double L2NormE   = std::sqrt(Stokes.vel_idx.GetEx().ParDot(VectorCL(Stokes.M.Data*e), false, e, true, true));
        double L2NormUnew= std::sqrt(Stokes.vel_idx.GetEx().ParDot(VectorCL(Stokes.M.Data*Stokes.v.Data), false, Stokes.v.Data, true, true));
        double norm_diff = L2NormE/L2NormUnew;
        double timeDeriv = L2NormE/C.tm_StepSize;

        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cout << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec.\n"
                      << "  (Step "<<step<<") => relative L2-Norm of velocity-change: " << norm_diff
                      << ", L2-Norm of velocity: " << L2NormUnew
                      << ", time-derivative: "<<timeDeriv
                      << std::endl;
            // Store measured values
            DROPS_LOGGER_SETVALUE("NavStokesCoupledLevelset",duration);
        }

        // Write out solution
        if (C.ens_EnsightOut && step%C.ens_EnsightOut==0)
            ensightwriter.write();
        if (C.vtk_VTKOut && step%C.vtk_VTKOut==0)
            vtkwriter.write();

        // Write out droplet information
        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        step_time.Stop();
        duration=step_time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cout <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
            DROPS_LOGGER_SETVALUE("TotalStep",duration);
            // Tell the logger class that the calculation step is over now
            DROPS_LOGGER_NEXTSTEP();
        }

        if (timeDeriv<=C.stk_XFEMStab){     // use XFEMStab as criteria
            IF_MASTER
            std::cout << " ==> norm velocity has just changed by "<< norm_diff
                        << ", tolerance of "<<C.stk_XFEMStab<<" has been reached!"
                        << std::endl;
            break;
        }
    }
    if (step==C.tm_NumSteps){
        IF_MASTER
            std::cout<< " ============> ATTENTION: TOLERANCE NOT REACHED"<<std::endl;
    }

    Times.IncCounter(step);
    timeIntegrationTimer.Stop();
    Times.AddTime(T_solve_NS, timeIntegrationTimer.GetTime());
    delete instatStokesSolver;
}


template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adapt)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& mg= Stokes.GetMG();
    ParTimerCL time;
    double duration;

    // Set parameter of the surface tension
    SurfaceTensionCL::eps           = C.sft_JumpWidth;
    SurfaceTensionCL::lambda        = C.sft_RelPos;
    SurfaceTensionCL::sigma         = Stokes.GetCoeff().SurfTens;
    SurfaceTensionCL::sigma_dirt_fac= C.sft_DirtFactor;

    instat_scalar_fun_ptr sigmap  = 0;
    instat_vector_fun_ptr gsigmap = 0;
    if (C.sft_VarTension)
    {
        sigmap  = &SurfaceTensionCL::sigma_step;
        gsigmap = &SurfaceTensionCL::gsigma_step;
    }
    else
    {
        sigmap  = &SurfaceTensionCL::sigmaf;
        gsigmap = &SurfaceTensionCL::gsigma;
    }

    LevelsetP2CL lset( mg, sigmap, gsigmap, C.lvs_SD, C.lvs_CurvDiff, C.rpm_NarrowBand);
    if (C.sft_VarTension)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    mg.ElemInfo(std::cout);

    std::ofstream *infofile=0;
    if (ProcCL::IamMaster())
        infofile = new std::ofstream( string(C.ens_EnsCase + ".info").c_str());
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    if (C.inf_CheckMG && !ProcCL::Check( CheckParMultiGrid(adapt.GetPMG())) )
         throw DROPSErrCL("MultiGrid is incorrect!");

    LevelsetRepairCL lsetrepair( lset);
    adapt.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes);
    adapt.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset);
    adapt.push_back( &prrepair);

    MLIdxDescCL *vidx= &Stokes.vel_idx,
                *pidx= &Stokes.pr_idx;
    IdxDescCL   *lidx= &lset.idx;
    VecDescCL *v   = &Stokes.v,
              *p   = &Stokes.p,
              *l   = &lset.Phi;

    Stokes.CreateNumberingVel( mg.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), pidx);
    lset.CreateNumbering     ( mg.GetLastLevel(), lidx);

    // Tell matrices and vectors about the numbering
    v->SetIdx( vidx);               p->SetIdx( pidx);
    l->SetIdx( lidx);
    Stokes.b.SetIdx( vidx);         Stokes.c.SetIdx( pidx);
    Stokes.A.SetIdx(vidx, vidx);    Stokes.B.SetIdx(pidx, vidx);
    Stokes.M.SetIdx(vidx, vidx);    Stokes.N.SetIdx(vidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx); Stokes.prA.SetIdx( pidx, pidx);

    if (C.inf_PrintNumUnk)
        DisplayUnks(Stokes, lset, mg);

    //Setup initial problem
    if (ProcCL::IamMaster())
        std::cout << "=================================================================================== Init:\n"
                    << "==> Initialize Problem\n";
    InitProblemWithDrop(Stokes, lset);

    SolveCoupledNS(Stokes, lset, adapt);

    if (ProcCL::IamMaster())
        std::cout << "=================================================================================== Transform:\n";

    if (C.qlg_Quad){
        if (ProcCL::IamMaster())
            std::cout << "Write out geometry and velocity on the quadrilateral grid\n - geometry" << std::endl;
        time.Reset();
        LevelsetP2CL::const_DiscSolCL lset_sol=lset.GetSolution();
        QuadOutCL brickout( mg, &lset.idx, C.qlg_GridX, C.qlg_GridY, C.qlg_GridZ, C.qlg_Stepsize, C.qlg_Barycenter, C.qlg_Rotation);
        brickout.putGeom(C.qlg_FileName + std::string(".geo"));
        if (ProcCL::IamMaster())
            std::cout << " - (integrated) velocities"<<std::endl;
        brickout.putVector(C.qlg_FileName + std::string(".vel_norm"),
                           C.qlg_FileName + std::string(".velY"),
                           C.qlg_FileName + std::string(".velZ"),
                           Stokes.GetVelSolution());
        time.Stop();
        duration=time.GetTime();
        if (ProcCL::IamMaster())
            std::cout << " => Took "<<duration<<" sec."<<std::endl;
        Times.AddTime(T_transform_to_quad, duration);
    }

    if (ProcCL::IamMaster()){
        infofile->close();
        delete infofile;
    }
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
  try
  {
    DROPS::ParTimerCL::TestBandwidth(std::cout);
    if (argc!=2)
    {
        IF_MASTER
          std::cout << "You have to specify one parameter:\n\t"
                    << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        IF_MASTER
          std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    IF_MASTER
      std::cout << C << std::endl;

    DROPS::ParTimerCL alltime;
    SetDescriber();

    typedef DROPS::ZeroFlowCL                             CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    // Create Geometry
    DROPS::MultiGridCL      *mg=0;
    DROPS::StokesBndDataCL  *bnddata=0;
    CreateGeom(mg, bnddata, DROPS::Inflow, C.dmc_MeshFile, C.dmc_GeomType, C.dmc_BoundaryType, C.rst_Inputfile, C.exp_RadInlet);

    mg->SizeInfo(std::cout);

    // Init problem
    DROPS::EllipsoidCL::Init( C.exp_PosDrop, C.exp_RadDrop );
    DROPS::AdapTriangCL adap( *mg, C.ref_Width, 0, C.ref_FinestLevel, C.ref_RefineStrategy);
    mg->SizeInfo(std::cout);

    adap.MakeInitialTriang( DROPS::DistanceYToBarycenterOfDrop);
    MyStokesCL prob(*mg, DROPS::ZeroFlowCL(C), *bnddata);

    // Solve the problem
    Strategy( prob, adap);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cout);

    // Writeout logging-results
    if (DROPS::ProcCL::IamMaster()){
        // DROPS_LOGGER_WRITEOUT("TestMzelleAdaptPar", (DROPS::ProcCL::Size()==1));
        std::cout << C.ref_FinestLevel << std::endl;

        // build specific name for logfile.
        char refLevel[32];
        sprintf(refLevel,"%d",C.ref_FinestLevel);
        std::string refStr(&refLevel[0]);

        //std::string basename("/home/th179319/DropsTimeMeas/results/TestMzelleAdaptPar_Level");
        std::string basename("./TestMzelleAdaptPar_Level");
        basename.append(refLevel);
        basename.append("_Proc");
        //DROPS_LOGGER_WRITEOUT(basename.c_str(), (DROPS::ProcCL::Size()==1));
    }

    // free memory
    if (mg)      delete mg;
    if (bnddata) delete bnddata;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


