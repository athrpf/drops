//**************************************************************************
// File:    TwoPhaseFlowPar.cpp                                            *
// Content: parallel solver for a 2-phase problem                          *
// Author:  Sven Gross, Joerg Grande, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   04. September 2008                                             *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TwoPhaseFlowPar.cpp
/// \brief Parallel solver for 2-phase problems

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
#include "stokes/integrTime.h"
#include "levelset/adaptriang.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "out/gridOut.h"
#include "parallel/parmgserialization.h"

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


DROPS::ParamParBrickFlowCL C;

/// \brief Inflow for Mzelle-Experiment
DROPS::SVectorCL<3> InflowMzelle( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

/// \brief Inflow for Brickflow-Experiment
DROPS::SVectorCL<3> InflowBrick( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = p[0]*(2*C.r_inlet-p[0]) / (C.r_inlet*C.r_inlet),
                 z = p[2]*(2*C.r_inlet-p[2]) / (C.r_inlet*C.r_inlet);

    ret[1]= x * z * C.Anstroem * (1-C.inflow_ampl*std::cos(2*M_PI*C.inflow_freq*t));
    return ret;
}

namespace DROPS{

template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    if (ProcCL::IamMaster())
        std::cerr << "=================================================================================== Init:\n"
                    << "==> Initialize Problem\n";

    ParTimerCL time;
    double duration;

    ExchangeCL& ExV = Stokes.GetEx(Stokes.velocity);
    ExchangeCL& ExP = Stokes.GetEx(Stokes.pressure);
//     ExchangeCL& ExL = lset.GetEx();
    MultiGridCL& mg=pmg.GetMG();

    switch (C.IniCond)
    {
      // zero flow
      case 0:
      {
        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        DisplayUnks(Stokes, lset, mg);

        // tell parmultigrid about the unknowns
        pmg.AttachTo(0, &Stokes.v); pmg.AttachTo(&Stokes.v, &Stokes.GetBndData().Vel);
        pmg.AttachTo(1, &Stokes.p); pmg.AttachTo(&Stokes.p, &Stokes.GetBndData().Pr);
        pmg.AttachTo(2, &lset.Phi); pmg.AttachTo(&lset.Phi, &lset.GetBndData());

        // Init velocity and level-set
        Stokes.InitVel( &Stokes.v, Null);
        lset.Init( EllipsoidCL::DistanceFct);

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);
      } break;

      // stationary flow with/without drop
      case 1: case 2:
      {
        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);

        // tell parmultigrid about the unknowns
        pmg.AttachTo(0, &Stokes.v); pmg.AttachTo(&Stokes.v, &Stokes.GetBndData().Vel);
        pmg.AttachTo(1, &Stokes.p); pmg.AttachTo(&Stokes.p, &Stokes.GetBndData().Pr);
        pmg.AttachTo(2, &lset.Phi); pmg.AttachTo(&lset.Phi, &lset.GetBndData());

        // Init velocity and level-set
        Stokes.InitVel( &Stokes.v, Null);
        if (C.IniCond==2){
            lset.Init( &One);
        }
        else{
            lset.Init( EllipsoidCL::DistanceFct);
        }

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        // Setting up solvers
        typedef ParDummyPcCL<ExchangeCL> SPcT;
        SPcT ispc(ExP);
        typedef ParJac0CL<ExchangeCL>  APcPcT;
        APcPcT Apcpc(ExV);
        typedef ParPCGSolverCL<APcPcT,ExchangeCL> ASolverT;
        ASolverT Asolver( 500, 0.02, ExV, Apcpc, /*relative=*/ true, /*accur*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver/*, &std::cerr*/);
        typedef ParInexactUzawaCL<APcT, SPcT, APC_SYM, ExchangeCL, ExchangeCL> OseenSolverT;
        OseenSolverT schurSolver( Apc, ispc, ExV, ExP, C.outer_iter, C.outer_tol, 0.1);

        VelVecDescCL curv(&Stokes.vel_idx),
                     cplN(&Stokes.vel_idx);

        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        curv.Data=0.;
        lset.AccumulateBndIntegral( curv);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
        {
            std::cerr << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";
            // Log measured duration.
            DROPS_LOGGER_SETVALUE("StokesCurv",duration);
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
            DROPS_LOGGER_SETVALUE("SolStokesTime",duration);
            DROPS_LOGGER_SETVALUE("SolStokesIter",iters);
        }
      }break;

      case 3:
      {
        // Create and number (on master) indices
        IdxDescCL read_vel_idx(3,3,0,0), read_pr_idx(1,0,0,0), read_lset_idx(1,1,0,0);

        if (ProcCL::IamMaster()){
            if (C.checkMG){
                std::ofstream sanity_file("ser_sanity.txt");
                if (mg.IsSane(sanity_file))
                    std::cerr << "Multigrid is sane on master\n";
                else
                    std::cerr << "Multigrid is not sane on master\n ==> Check right dimension of the computational domain!";
            }
            Stokes.CreateNumberingVel( mg.GetLastLevel(), &read_vel_idx,  false);
            Stokes.CreateNumberingPr ( mg.GetLastLevel(), &read_pr_idx,   false);
            lset.CreateNumbering     ( mg.GetLastLevel(), &read_lset_idx, false);
        }

        // Allocate memory for unknowns
        VecDescCL read_vel_vec(&read_vel_idx), read_pr_vec(&read_pr_idx), read_lset_vec(&read_lset_idx);

        // Tell parmultigrid about unknowns
        pmg.AttachTo(0, &read_vel_vec);  pmg.AttachTo(&read_vel_vec, &Stokes.GetBndData().Vel);
        pmg.AttachTo(1, &read_pr_vec);   pmg.AttachTo(&read_pr_vec, &Stokes.GetBndData().Pr);
        pmg.AttachTo(2, &read_lset_vec); pmg.AttachTo(&read_lset_vec, &lset.GetBndData());

        // Read DOF
        if (ProcCL::IamMaster()){
            ReadDOF(mg, &read_vel_vec,  std::string(C.ser_dir+"velocity"));
            ReadDOF(mg, &read_pr_vec,   std::string(C.ser_dir+"pressure"));
            ReadDOF(mg, &read_lset_vec, std::string(C.ser_dir+"level-set"));
        }

        // Distribute MG
        lb.DoMigration();
        if (C.checkMG && !Check( CheckParMultiGrid(pmg)) )
            throw DROPSErrCL("MultiGrid is incorrect!");

        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        read_vel_idx.TriangLevel = Stokes.vel_idx.TriangLevel;
        read_pr_idx.TriangLevel  = Stokes.pr_idx.TriangLevel;
        read_lset_idx.TriangLevel= lset.idx.TriangLevel;

        pmg.HandleNewIdx(&read_vel_idx,  &Stokes.v);
        pmg.HandleNewIdx(&read_pr_idx,   &Stokes.p);
        pmg.HandleNewIdx(&read_lset_idx, &lset.Phi);
        pmg.DelAllUnkRecv();
        pmg.DeleteRecvBuffer();

        Stokes.DeleteNumberingVel(&read_vel_idx);
        Stokes.DeleteNumberingPr( &read_pr_idx);
        lset.DeleteNumbering(     &read_lset_idx);

        pmg.AttachTo(0, &Stokes.v); pmg.AttachTo(1, &Stokes.p); pmg.AttachTo(2, &lset.Phi);

        if (ProcCL::IamMaster())
                std::cerr << "- Initial Conditions successfull read\n";
      } break;
    }
}

template<typename Coeff>
  void SolveCoupledNS(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, AdapTriangCL& adap)
{
    ParTimerCL time;
    double duration;
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    // writer in vtk-format
    typedef TwoPhaseVTKCL<StokesProblemT, LevelsetP2CL> VTKWriterT;
    VTKWriterT vtkwriter(adap.GetMG(), Stokes, lset, (C.vtk ? C.num_steps/C.vtk+1 : 0),
                         std::string(C.vtkDir + "/" + C.vtkName), C.vtkBinary);
    // writer for ensight format
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.ensDir, C.ensCase, C.geomName, /*adaptive=*/true,
                                  (C.ensight? C.num_steps/C.ensight+1 : 0), C.binary, C.masterOut);

    if (C.vtk)
        vtkwriter.write();
    if (C.ensight)
        ensightwriter.write();

    const double Vol= EllipsoidCL::GetVolume();
    double relVol = lset.GetVolume()/Vol;

    ExchangeCL& ExV = Stokes.GetEx(Stokes.velocity);
    ExchangeCL& ExP = Stokes.GetEx(Stokes.pressure);
    ExchangeCL& ExL = lset.GetEx();

    // linear solvers
    typedef ParDummyPcCL<ExchangeCL>  SPcT;
    SPcT ispc(ExP);
    typedef ParJac0CL<ExchangeCL>  APcPcT;
    APcPcT Apcpc(ExV);
    typedef ParPreGMResSolverCL<APcPcT, ExchangeCL>    ASolverT;        // GMRes-based APcT
    ASolverT Asolver( /*restart*/    100, /*iter*/ 500,  /*tol*/   2e-5, ExV, Apcpc,
                       /*relative=*/ true, /*acc*/ true, /*modGS*/ false, RightPreconditioning, /*mod4par*/true);

//     ASolverT Asolver( C.navstokes_pc_restart, C.navstokes_pc_iter, C.navstokes_pc_reltol, ExV, Apcpc,
//                         /*relative=*/ true, /*acc*/ true, /*modGS*/false, RightPreconditioning, /*mod4par*/false);

    typedef SolverAsPreCL<ASolverT> APcT;
    APcT Apc( Asolver/*, &std::cerr*/);

    // stokes solver
    typedef ParInexactUzawaCL<APcT, SPcT, APC_OTHER, ExchangeCL, ExchangeCL> OseenSolverT;
    OseenSolverT oseensolver( Apc, ispc, ExV, ExP, C.outer_iter, C.outer_tol, 0.2, 500, &std::cerr);

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT, OseenSolverT> NSSolverT;
    NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red, &std::cerr);

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

        cpl.DoStep( C.cpl_iter);

        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cerr << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec."<<std::endl;
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
            lset.ReparamFastMarching( ExL, C.RepMethod, C.RepMethod==3);
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
    }
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

    LevelsetP2CL lset( mg, sigmap, gsigmap, C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff, C.NarrowBand);
    if (C.st_var)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    AdapTriangCL adapt( pmg, lb, C.ref_width, 0, C.ref_flevel);
    adapt.MakeInitialTriang(EllipsoidCL::DistanceFct);

    std::ofstream *infofile=0;
    if (ProcCL::IamMaster())
        infofile = new std::ofstream( string(C.ensCase + ".info").c_str());

    IFInfo.Init(infofile);

    if (C.checkMG && !Check( CheckParMultiGrid(pmg)) )
         throw DROPSErrCL("MultiGrid is incorrect!");

    LevelsetRepairCL lsetrepair( lset, pmg);
    adapt.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes, pmg);
    adapt.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, pmg);
    adapt.push_back( &prrepair);

    IdxDescCL *vidx= &Stokes.vel_idx,
              *pidx= &Stokes.pr_idx,
              *lidx= &lset.idx;
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

//     // Tell parmultigrid about unknowns
    pmg.AttachTo(0, v); pmg.AttachTo(v, &Stokes.GetBndData().Vel);
    pmg.AttachTo(1, p); pmg.AttachTo(p, &Stokes.GetBndData().Pr);
    pmg.AttachTo(2, l); pmg.AttachTo(l, &lset.GetBndData());

    if (C.printNumUnk)
        DisplayUnks(Stokes, lset, pmg.GetMG());

    //Setup initial problem
    InitProblemWithDrop(Stokes, lset, pmg, lb);

    SolveCoupledNS(Stokes, lset, adapt);

    if (C.quad){
        if (ProcCL::IamMaster())
            std::cerr << "Write out geometry and velocity on the quadrilateral grid\n - geometry" << std::endl;
        time.Reset();
        LevelsetP2CL::const_DiscSolCL lset_sol=lset.GetSolution();
        QuadOutCL brickout( mg, &lset.idx, C.gridX, C.gridY, C.gridZ, C.stepsize, C.barycenter, C.rotation);
        brickout.putGeom("wrap.geo");
        if (ProcCL::IamMaster())
            std::cerr << " - (integrated) velocities"<<std::endl;
        brickout.putVector("wrap.vel_norm", "wrap.velY", "wrap.velZ", Stokes.GetVelSolution(), &lset_sol);
        time.Stop();
        duration=time.GetTime();
        if (ProcCL::IamMaster())
            std::cerr << " => Took "<<duration<<" sec."<<std::endl;
    }

    if (ProcCL::IamMaster()){
        infofile->close();
        delete infofile;
    }
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  DROPS::ProcInitCL procinit(&argc, &argv);.
  DROPS::ParMultiGridInitCL pmginit;
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

    typedef DROPS::ZeroFlowCL                             CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    // Create Geometry
    DROPS::ParMultiGridCL   *pmg=0;
    DROPS::MultiGridCL      *mg=0;
    DROPS::LoadBalHandlerCL *lb=0;
    DROPS::StokesBndDataCL  *bnddata=0;
    GeometryType geomtype= (C.GeomType==0 ? newMZelle : (C.meshfile=="none" ? newBrick : fileBrick) );
    CreateGeom(mg, pmg, lb, bnddata, (C.GeomType==0 ? InflowMzelle : InflowBrick), C.meshfile, C.refineStrategy, geomtype, C.r_inlet);

    // Init problem
    EllipsoidCL::Init( C.Mitte, C.Radius );
    MyStokesCL prob(*mg, DROPS::ZeroFlowCL(C), *bnddata);

    // Solve the problem
    Strategy( prob, *pmg, *lb);    // do all the stuff

    // free memory
    delete mg;
    delete lb;
    delete bnddata;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


