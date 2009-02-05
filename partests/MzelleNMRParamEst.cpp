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
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

/// \brief Get y-distance to midpoint of the meassurements ...
double DistanceToBarycenterOfQuadGrid( const DROPS::Point3DCL& p)
{
    return p[1]-C.barycenter[1];
}

/// \brief Get y-distance to midpoint of the droplet ...
double DistanceYToBarycenterOfDrop( const DROPS::Point3DCL& p)
{
    return p[1]-C.Mitte[1];
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

    ExchangeCL& ExV = Stokes.vel_idx.GetEx();
    ExchangeCL& ExP = Stokes.pr_idx.GetEx();
    ExchangeCL& ExL = lset.idx.GetEx();

    switch (C.IniCond)
    {
      // stationary flow with/without drop
      case 1: case 2:
      {
        StokesSolverParamST statStokesParam(C);
        statStokesParam.StokesMethod= 20301;
        statStokesParam.num_steps   = 0;
        statStokesParam.dt          = 0.;
        statStokesParam.pcA_tol     = 0.02;
        StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> statStokesSolverFactory(Stokes, statStokesParam);
        StokesSolverBaseCL* statStokesSolver= statStokesSolverFactory.CreateStokesSolver();
        StokesSolverBaseCL& schurSolver=*statStokesSolver;

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
            std::cerr << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";
            // Log measured duration.
            DROPS_LOGGER_SETVALUE("StokesCurvInit",duration);
        }

        //Solve initial problem
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
        {
            std::cerr << "- Solving Stokes for initialization took "<<duration<<" sec, "
                      << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<schurSolver.GetResid()<<'\n';
            // Log measured duration
            DROPS_LOGGER_SETVALUE("SolStokesTime",duration);
            DROPS_LOGGER_SETVALUE("SolStokesIter",iters);
        }
        delete statStokesSolver;
      }break;
      case 3:
        {
            ReadEnsightP2SolCL reader( Stokes.GetMG(), false);
            reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel, ExV);
            reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr,  ExP);
            reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData()     ,  ExL);
            if (ProcCL::IamMaster())
                std::cerr << "- Initial Conditions successfull read\n";
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
    VTKWriterT vtkwriter(adap.GetMG(), Stokes, lset, (C.vtk ? C.num_steps/C.vtk+1 : 0),
                         std::string(C.vtkDir + "/" + C.vtkName), C.vtkBinary);
    // writer for ensight format
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.EnsDir, C.EnsCase, C.geomName, /*adaptive=*/true,
                                  (C.ensight? C.num_steps/C.ensight+1 : 0), C.binary, C.masterOut);

    if (C.vtk)
        vtkwriter.write();
    if (C.ensight)
        ensightwriter.write();

    const double Vol= EllipsoidCL::GetVolume();
    double relVol = lset.GetVolume()/Vol;

    StokesSolverParamST instatStokesParam(C);
    StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> instatStokesSolverFactory(Stokes, instatStokesParam);
    StokesSolverBaseCL* instatStokesSolver= instatStokesSolverFactory.CreateStokesSolver();
    StokesSolverBaseCL& oseensolver= *instatStokesSolver;

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
    NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red);

    // coupling levelset NavStokes
    time.Reset();
    // Coupling Navier-Stokes with Levelset
//    typedef LinThetaScheme2PhaseCL <StokesProblemT, NSSolverT> CouplingT;
//    CouplingT cpl( Stokes, lset, nssolver, C.nonlinear);

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
                      << " Idx for vel  "<<Stokes.v.RowIdx->GetIdx()
                      << "\n Idx for pr   "<<Stokes.p.RowIdx->GetIdx()
                      << "\n Idx for lset "<<lset.Phi.RowIdx->GetIdx()<<std::endl;
        time.Reset();

        if (C.ref_freq && step%C.ref_freq==0)
        {
            if (ProcCL::IamMaster())
                std::cerr << "==> Adaptive Refinement of MultiGrid"<<std::endl;

            adap.UpdateTriang( lset);

            if (C.checkMG && !Check( CheckParMultiGrid(adap.GetPMG())) )
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

        if (C.printNumUnk)
            DisplayUnks(Stokes, lset, adap.GetMG());

        if (ProcCL::IamMaster())
            std::cerr << "==> Solving coupled Levelset-Navier-Stokes problem ....\n";

        VectorCL u_old= Stokes.v.Data;

        cpl.DoStep( C.cpl_iter);

        double rel_norm_diff=  Stokes.vel_idx.GetEx().Norm( VectorCL(Stokes.v.Data-u_old), true)
                             / Stokes.vel_idx.GetEx().Norm( Stokes.v.Data, true);
        double norm_vel = Stokes.vel_idx.GetEx().Norm( Stokes.v.Data, true);

        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cerr << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec.\n"
                      << "  => Relative Difference of velocity to previous solution: " << rel_norm_diff
                      << ", norm of velocity-vector is "<< norm_vel
                      << std::endl;
            // Store measured values
            DROPS_LOGGER_SETVALUE("NavStokesCoupledLevelset",duration);
        }

        // Write out solution
        if (C.ensight && step%C.ensight==0)
            ensightwriter.write();
        if (C.vtk && step%C.vtk==0)
            vtkwriter.write();

        // Write out droplet information
        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        // Reparametrization of levelset function
        if (C.RepFreq && step%C.RepFreq==0)
        {
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Reparametrization\n"
                          << "- rel. Volume: " << relVol << std::endl;
            time.Reset();
            lset.ReparamFastMarching( C.RepMethod, C.RepMethod==3);
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cerr << "- Reparametrization took "<<duration<<" sec."<<std::endl;
                DROPS_LOGGER_SETVALUE("Reparametrization",duration);
            }
        }

        if (ProcCL::IamMaster())
            std::cerr << "- rel. Volume: " << relVol << std::endl;

        if (C.VolCorr)
        {
            time.Reset();
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Adjust volume ...\n";
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            lset.Phi.Data+= dphi;
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cerr << "- Lifting level-set function took "<< duration << " sec\n"
                          << "- Volume correction " << dphi <<", new rel. Volume is " << relVol << std::endl;
                DROPS_LOGGER_SETVALUE("VolumeCorrection",duration);
            }
        }

        step_time.Stop();
        duration=step_time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cerr <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
            DROPS_LOGGER_SETVALUE("TotalStep",duration);
            // Tell the logger class that the calculation step is over now
            DROPS_LOGGER_NEXTSTEP();
        }

        if (rel_norm_diff<=C.XFEMStab){     // use XFEMStab as criteria
            IF_MASTER
              std::cerr << " ==> relative difference in velocity has just changed by "<< rel_norm_diff
                        << ", tolerance of "<<C.XFEMStab<<" has been reached!"
                        << std::endl;
            break;
        }
    }
    timeIntegrationTimer.Stop();
    Times.AddTime(T_solve_NS, timeIntegrationTimer.GetTime());
    delete instatStokesSolver;
}


template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& mg= pmg.GetMG();
    ParTimerCL time;
    double duration;

    // Set parameter of the surface tension
    SurfaceTensionCL::eps           = C.st_jumpWidth;
    SurfaceTensionCL::lambda        = C.st_relPos;
    SurfaceTensionCL::sigma         = Stokes.GetCoeff().SurfTens;
    SurfaceTensionCL::sigma_dirt_fac= C.st_red;

    instat_scalar_fun_ptr sigmap  = 0;
    instat_vector_fun_ptr gsigmap = 0;
    if (C.st_var)
    {
        sigmap  = &SurfaceTensionCL::sigma_step;
        gsigmap = &SurfaceTensionCL::gsigma_step;
    }
    else
    {
        sigmap  = &SurfaceTensionCL::sigmaf;
        gsigmap = &SurfaceTensionCL::gsigma;
    }

    LevelsetP2CL lset( mg, sigmap, gsigmap, C.theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff, C.NarrowBand);
    if (C.st_var)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    AdapTriangCL adapt( pmg, lb, C.ref_width, 0, C.ref_flevel);
    //typedef double (*DistFctT)(const DROPS::Point3DCL&);
    //const DistFctT &distfnct=::EllipsoidCL::DistanceFct;
    //adapt.MakeInitialTriang(distfnct);
    adapt.MakeInitialTriang(DistanceYToBarycenterOfDrop);
    mg.ElemInfo(std::cerr);

    std::ofstream *infofile=0;
    if (ProcCL::IamMaster())
        infofile = new std::ofstream( string(C.EnsCase + ".info").c_str());

    IFInfo.Init(infofile);

    if (C.checkMG && !Check( CheckParMultiGrid(pmg)) )
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

    if (C.printNumUnk)
        DisplayUnks(Stokes, lset, pmg.GetMG());

    //Setup initial problem
    if (ProcCL::IamMaster())
        std::cerr << "=================================================================================== Init:\n"
                    << "==> Initialize Problem\n";
    InitProblemWithDrop(Stokes, lset);

    SolveCoupledNS(Stokes, lset, adapt);

    if (C.quad){
        if (ProcCL::IamMaster())
            std::cerr << "Write out geometry and velocity on the quadrilateral grid\n - geometry" << std::endl;
        time.Reset();
        LevelsetP2CL::const_DiscSolCL lset_sol=lset.GetSolution();
        QuadOutCL brickout( mg, &lset.idx, C.gridX, C.gridY, C.gridZ, C.stepsize, C.barycenter, C.rotation);
        brickout.putGeom(C.quadFileName + std::string(".geo"));
        if (ProcCL::IamMaster())
            std::cerr << " - (integrated) velocities"<<std::endl;
        brickout.putVector(C.quadFileName + std::string(".vel_norm"),
                           C.quadFileName + std::string(".velY"),
                           C.quadFileName + std::string(".velZ"),
                           Stokes.GetVelSolution());
        time.Stop();
        duration=time.GetTime();
        if (ProcCL::IamMaster())
            std::cerr << " => Took "<<duration<<" sec."<<std::endl;
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
  DROPS::ProcCL Proc(&argc, &argv);
  try
  {
    if (argc!=2)
    {
        IF_MASTER
          std::cerr << "You have to specify one parameter:\n\t"
                    << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        IF_MASTER
          std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    IF_MASTER
      std::cerr << C << std::endl;

    DROPS::ParTimerCL alltime;
    SetDescriber();

    typedef DROPS::ZeroFlowCL                             CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    // Create Geometry
    DROPS::ParMultiGridCL   *pmg=0;
    DROPS::MultiGridCL      *mg=0;
    DROPS::LoadBalHandlerCL *lb=0;
    DROPS::StokesBndDataCL  *bnddata=0;
    CreateGeom(mg, pmg, lb, bnddata, DROPS::Inflow, C.meshfile, C.refineStrategy, newMZelle, C.r_inlet);

    mg->SizeInfo(std::cerr);

    // Init problem
    EllipsoidCL::Init( C.Mitte, C.Radius );
    MyStokesCL prob(*mg, DROPS::ZeroFlowCL(C), *bnddata);

    // Solve the problem
    Strategy( prob, *pmg, *lb);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cerr);

    // Writeout logging-results
    if (DROPS::ProcCL::IamMaster()){
        // DROPS_LOGGER_WRITEOUT("TestMzelleAdaptPar", (DROPS::ProcCL::Size()==1));
        std::cout << C.ref_flevel << std::endl;

        // build specific name for logfile.
        char refLevel[32];
        sprintf(refLevel,"%d",C.ref_flevel);
        std::string refStr(&refLevel[0]);

        //std::string basename("/home/th179319/DropsTimeMeas/results/TestMzelleAdaptPar_Level");
        std::string basename("./TestMzelleAdaptPar_Level");
        basename.append(refLevel);
        basename.append("_Proc");
        DROPS_LOGGER_WRITEOUT(basename.c_str(), (DROPS::ProcCL::Size()==1));
    }

    // free memory
    if (pmg)     delete pmg;
    if (mg)      delete mg;
    if (lb)      delete lb;
    if (bnddata) delete bnddata;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

