//**************************************************************************
// File:    TestBrickflowPar.cpp                                           *
// Content: parallel solver for a simple 2-phase problem                   *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestBrickflowPar.cpp
/// \brief Testing parallel solver for a simple 2-phase problem

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


enum TimePart{
    T_create_pmg,
    T_disc_init,
    T_create_ex,
    T_solve_NS_init
};

DROPS::TimeStoreCL Times(4);

void SetDescriber()
{
    Times.SetDescriber(T_create_pmg,    "Initialization of parallel multigrid");
    Times.SetDescriber(T_disc_init,     "Discretizing Stokes/Curv for initial velocities");
    Times.SetDescriber(T_create_ex,     "Create ExchangeCL");
    Times.SetDescriber(T_solve_NS_init, "Solving NavStokes for initial velocities");
    Times.SetCounterDescriber("Moved MultiNodes");
}

DROPS::ParamParBrickFlowCL C;

const char line[] ="------------------------------------------------------------";

/****************************************************************************
* P R I N T  M U L T I G R I D  I N  A S C I I  F I L E                     *
****************************************************************************/


bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream check(dat);
    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);
    check.close();
    return DROPS::ProcCL::Check(pmg_sane && mg_sane);
}

/// \brief Display number of unknowns
void DisplayNumUnknownsSparse(const DROPS::MultiGridCL& MG, const DROPS::VecDescCL& x, const std::string& name)
/// accumulated unknowns: accumulation of the unknowns of all processors. Some unknowns
///   are counted multiple due to the overlapping of edges and vertices
/// global unknowns: all unknowns are just once counted
{
    const DROPS::Ulint acc_num_unk = DROPS::ProcCL::GlobalSum(x.Data.size()),
                       glo_num_unk = x.RowIdx->GetGlobalNumUnknowns(MG);
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    if (DROPS::ProcCL::IamMaster())
        std::cerr << "  + Number of "<<name<<" DOF with index "<<idx<<" (accumulated/global):  "
                  <<acc_num_unk<< "/" <<glo_num_unk<< std::endl;
}


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = (p[0]/C.r_inlet)*(p[0]/C.r_inlet)-1,
                 z = (p[2]/C.r_inlet)*(p[2]/C.r_inlet)-1;

//     Inflow with waves
//     const double freq= 50, ampl= 0.1;
//     ret[1]= x * z * C.Anstroem * (1-ampl*std::cos(2*M_PI*freq*t));  // Rohr

    // constant inflow
    ret[1]= x * z * C.Anstroem;

    return ret;
}

// distance function of a droplet droplet
double DistanceFct1( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.Mitte-p;
    const double avgRad = cbrt(C.Radius[0]*C.Radius[1]* C.Radius[2]);
    d/= C.Radius/avgRad;
    return d.norm()-avgRad;
}

// middle
double DistanceFct( const DROPS::Point3DCL& p)
{
    return p[1]-0.015;
}

namespace DROPS{

template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    ParTimerCL time;
    double duration;

    MultiGridCL& mg=pmg.GetMG();

    switch (C.IniCond)
    {
      // zero flow
      case 0:
      {
        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        DisplayUnks(Stokes, lset, mg);

        // Initial velocity is zero
        Stokes.InitVel( &Stokes.v, Null);
        // Initial levelset function is distance to the drop
        lset.Init( DistanceFct1);

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);
      } break;
      // stationary flow with/without drop
      case 1: case 2:
      {
        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        // Initial velocity is zero
        Stokes.InitVel( &Stokes.v, Null);
        // Initial levelset function is distance to the drop
        lset.Init( DistanceFct1);

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        StokesSolverParamST statStokesParam(C);
        statStokesParam.StokesMethod= 20301;
        statStokesParam.num_steps   = 0;
        statStokesParam.dt          = 0.;
        statStokesParam.pcA_tol     = 0.02;
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
            statStokesSolver->Solve( Stokes.A.Data, Stokes.B.Data,
                               Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
            if (ProcCL::IamMaster())
                std::cerr << "- Solving lin. Stokes ("<<step<<"): iter "<<statStokesSolver->GetIter()
                            <<", resid "<<statStokesSolver->GetResid()<<std::endl;
            ++step; iters+= statStokesSolver->GetIter();
        } while (statStokesSolver->GetIter() > 0);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
        {
            std::cerr << "- Solving Stokes for initialization took "<<duration<<" sec, "
                        << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<statStokesSolver->GetResid()<<'\n';
            DROPS_LOGGER_SETVALUE("SolStokesTime",duration);
            DROPS_LOGGER_SETVALUE("SolStokesIter",iters);
        }
        delete statStokesSolver;
      }break;

      case 3:
      {
        // Create and number (on master) indices
        MLIdxDescCL read_vel_idx( vecP2_FE), read_pr_idx( P1_FE);
        IdxDescCL read_lset_idx( P2_FE);

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

        // Read DOF
        if (ProcCL::IamMaster()){
            ReadDOF(mg, &read_vel_vec,  std::string(C.ser_dir+"velocity"));
            ReadDOF(mg, &read_pr_vec,   std::string(C.ser_dir+"pressure"));
            ReadDOF(mg, &read_lset_vec, std::string(C.ser_dir+"level-set"));
        }

        // Distribute MG
        lb.DoMigration();
        if (C.checkMG && !ProcCL::Check( CheckParMultiGrid(pmg)) )
            throw DROPSErrCL("MultiGrid is incorrect!");

        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        /// \todo Hier muessen noch die Level der Indices gesetzt werden. Warum?
//         read_vel_idx.TriangLevel = Stokes.vel_idx.TriangLevel;
//         read_pr_idx.TriangLevel  = Stokes.pr_idx.TriangLevel;
//         read_lset_idx.TriangLevel= lset.idx.TriangLevel;

        pmg.HandleNewIdx(&read_vel_idx,  &Stokes.v);
        pmg.HandleNewIdx(&read_pr_idx,   &Stokes.p);
        pmg.HandleNewIdx(&read_lset_idx, &lset.Phi);
        pmg.DelAllUnkRecv();
        pmg.DeleteRecvBuffer();

        Stokes.DeleteNumbering(&read_vel_idx);
        Stokes.DeleteNumbering( &read_pr_idx);
        lset.DeleteNumbering(     &read_lset_idx);

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

    // writer for vtk-format
    typedef TwoPhaseVTKCL<StokesProblemT, LevelsetP2CL> VTKWriterT;
    VTKWriterT vtkwriter(adap.GetMG(), Stokes, lset, (C.vtk ? C.num_steps/C.vtk+1 : 0),
                         std::string(C.vtkDir + "/" + C.vtkName), C.vtkBinary);
    // writer for ensight format
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.EnsDir, C.EnsCase, C.geomName, /*adaptive=*/true,
                                  (C.ensight? C.num_steps/C.ensight+1 : 0), C.binary, C.masterOut);

    // Serialization
    typedef TwoPhaseSerializationCL<StokesProblemT, LevelsetP2CL> SerializationT;
    SerializationT serializer(adap.GetMG(), C.ser_dir, Stokes, lset, C.overwrite, C.num_steps/(C.overwrite+1));

    if (C.vtk)
        vtkwriter.write();
    if (C.ensight)
        ensightwriter.write();

    const double Vol= 4./3.*M_PI*C.Radius[0]*C.Radius[1]*C.Radius[2];
    double relVol = lset.GetVolume()/Vol;

    StokesSolverParamST instatStokesParam(C);
    StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> instatStokesSolverFactory(Stokes, instatStokesParam);
    StokesSolverBaseCL* instatStokesSolver= instatStokesSolverFactory.CreateStokesSolver();

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
    NSSolverT nssolver( Stokes, *instatStokesSolver, C.ns_iter, C.ns_tol, C.ns_red);

    // coupling levelset NavStokes
    time.Reset();
    // Coupling Navier-Stokes with Levelset
//     typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
    typedef LinThetaScheme2PhaseCL <StokesProblemT, NSSolverT> CouplingT;
    CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear, true);
//     typedef RecThetaScheme2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//     CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear, false);


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

            if (C.checkMG && !ProcCL::Check( CheckParMultiGrid(adap.GetPMG())) )
                throw DROPSErrCL("MultiGrid is incorrect!");

            if (adap.WasModified() )
            {
                cpl.Update();
                // don't forget to update the pr mass/stiff matrix for the schur compl. preconditioner!!
                // / \todo (of) Ueberfluessig!
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
            DROPS_LOGGER_SETVALUE("NavStokesCoupledLevelset", duration);
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

        // Correction of the volume
        if (C.VolCorr)
        {
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Adjust volume ...\n";
            time.Reset();
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            lset.Phi.Data+= dphi;
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cerr << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< " took "<<duration<<" sec"<<std::endl;
                DROPS_LOGGER_SETVALUE("VolumeCorrection",duration);
            }
        }

        // Serialization
        if (C.serialization && step%C.serialization==0)
        {
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Serialize data ...\n";
            time.Reset();
            serializer.Serialize(step);
            time.Stop(); duration=time.GetMaxTime();
            if (ProcCL::IamMaster())
                std::cerr << "- Serialization took " <<duration<< " sec"<<std::endl;
        }

        step_time.Stop();
        duration=step_time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cerr <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
            DROPS_LOGGER_SETVALUE("TotalStep",duration);
            DROPS_LOGGER_NEXTSTEP();
        }
    }
    delete instatStokesSolver;
}


template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adapt)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    MultiGridCL& mg= Stokes.GetMG();
    ParTimerCL time;

    // Set parameter of the surface tension
    SurfaceTensionCL::eps           = C.st_jumpWidth;
    SurfaceTensionCL::lambda        = C.st_relPos;
    SurfaceTensionCL::sigma         = Stokes.GetCoeff().SurfTens;
    SurfaceTensionCL::sigma_dirt_fac= C.st_red;

//    instat_vector_fun_ptr gsigma= &(SurfaceTensionCL::grad_sm_step);
    LevelsetP2CL lset( mg, &SurfaceTensionCL::sigma_step, &SurfaceTensionCL::gsigma_step,
                       C.theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff, C.NarrowBand);

    DisplayDetailedGeom(mg);

    LevelsetRepairCL lsetrepair( lset);
    adapt.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes);
    adapt.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset);
    adapt.push_back( &prrepair);

    /// \todo (of) Testen von SF_ImprovedLBVar!!!
    lset.SetSurfaceForce( SF_ImprovedLB);

    std::ofstream *infofile=0;
    if (ProcCL::IamMaster())
        infofile = new std::ofstream( string(C.EnsCase + ".info").c_str());

    IFInfo.Init(infofile);

    //Setup initial problem
    std::cerr << "=================================================================================== Init:\n"
              << "==> Initialize Problem\n";
    InitProblemWithDrop(Stokes, lset, adapt.GetPMG(), adapt.GetLb());

    SolveCoupledNS(Stokes, lset, adapt);

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
    std::cerr << C << std::endl;

    DROPS::ParTimerCL alltime;
    SetDescriber();

    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::ParTimerCL time;

    // Create Geometry
    DROPS::MultiGridCL      *mg=0;
    DROPS::StokesBndDataCL  *bnddata=0;

    CreateGeom(mg, bnddata, Inflow, C.meshfile, C.GeomType, C.bnd_type, C.deserialization_file, C.r_inlet);

    // Init problem
    EllipsoidCL::Init( C.Mitte, C.Radius );
    DROPS::AdapTriangCL adap( *mg, C.ref_width, 0, C.ref_flevel, C.refineStrategy);

    adap.MakeInitialTriang( EllipsoidCL::DistanceFct);

    MyStokesCL prob( *mg, ZeroFlowCL(C), *bnddata, DROPS::P1_FE, 0.0);

    Strategy( prob, adap);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cout);

    // Writeout logging-results
    if (DROPS::ProcCL::IamMaster()){
        // DROPS_LOGGER_WRITEOUT("TestBrickflowPar", (DROPS::ProcCL::Size()==1));

        // build specific name for logfile.
        char refLevel[32];
        sprintf(refLevel,"%d",C.ref_flevel);
        std::string refStr(&refLevel[0]);
        std::string basename("./TestBrickflowPar_Level_");
       // std::string basename("/home/th179319/DropsTimeMeas/results/TestBrickflowPar_Rad3_Level_");
        // std::string basename("/home/th179319/DropsTimeMeas/results/TestBrickflowPar_Level_");
        basename.append(refLevel);
        basename.append("_Proc");
        DROPS_LOGGER_WRITEOUT(basename.c_str(), true);
    }

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
