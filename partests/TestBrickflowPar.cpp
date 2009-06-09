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
        std::cout << "  + Number of "<<name<<" DOF with index "<<idx<<" (accumulated/global):  "
                  <<acc_num_unk<< "/" <<glo_num_unk<< std::endl;
}


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = (p[0]/C.exp_RadInlet)*(p[0]/C.exp_RadInlet)-1,
                 z = (p[2]/C.exp_RadInlet)*(p[2]/C.exp_RadInlet)-1;

//     Inflow with waves
//     const double freq= 50, ampl= 0.1;
//     ret[1]= x * z * C.Anstroem * (1-ampl*std::cos(2*M_PI*freq*t));  // Rohr

    // constant inflow
    ret[1]= x * z * C.exp_InflowVel;

    return ret;
}

// distance function of a droplet droplet
double DistanceFct1( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.exp_PosDrop-p;
    const double avgRad = cbrt(C.exp_RadDrop[0]*C.exp_RadDrop[1]* C.exp_RadDrop[2]);
    d/= C.exp_RadDrop/avgRad;
    return d.norm()-avgRad;
}

// middle
double DistanceFct( const DROPS::Point3DCL& p)
{
    return p[1]-0.015;
}

namespace DROPS{

template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, ParMultiGridCL& pmg)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    ParTimerCL time;
    double duration;

    MultiGridCL& mg=pmg.GetMG();

    // Create indices
    CreateIdxAndAssignIdx(Stokes, lset, mg);

    switch (C.dmc_InitialCond)
    {
      // zero flow
      case 0:
      {
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
        // Initial velocity is zero
        Stokes.InitVel( &Stokes.v, Null);
        // Initial levelset function is distance to the drop
        lset.Init( DistanceFct1);

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

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
            DROPS_LOGGER_SETVALUE("StokesCurv",duration);
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
            DROPS_LOGGER_SETVALUE("SolStokesTime",duration);
            DROPS_LOGGER_SETVALUE("SolStokesIter",iters);
        }
        delete statStokesSolver;
      }break;

      case -1:
      {
        MLIdxDescCL* pidx= &Stokes.pr_idx;
        // Read velocity, pressure and levelset
        ReadFEFromFile( lset.Phi, mg, C.dmc_InitialFile+"levelset");
        ReadFEFromFile( Stokes.v, mg, C.dmc_InitialFile+"velocity");
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        if (Stokes.UsesXFEM()) {
            VecDescCL pneg( pidx), ppos( pidx);
            ReadFEFromFile( pneg, mg, C.dmc_InitialFile+"pressureNeg");
            ReadFEFromFile( ppos, mg, C.dmc_InitialFile+"pressurePos");
            P1toP1X ( pidx->GetFinest(), Stokes.p.Data, pidx->GetFinest(), ppos.Data, pneg.Data, lset.Phi, mg);
        }
        else{
            ReadFEFromFile( Stokes.p, mg, C.dmc_InitialFile+"pressure");
        }

        if (ProcCL::IamMaster())
            std::cout << "- Initial Conditions successfull read\n";
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
    VTKWriterT vtkwriter(adap.GetMG(), Stokes, lset, (C.vtk_VTKOut ? C.tm_NumSteps/C.vtk_VTKOut+1 : 0),
                         std::string(C.vtk_VTKDir + "/" + C.vtk_VTKName), C.vtk_Binary);
    // writer for ensight format
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.ens_EnsDir, C.ens_EnsCase, C.ens_GeomName, /*adaptive=*/true,
                                  (C.ens_EnsightOut? C.tm_NumSteps/C.ens_EnsightOut+1 : 0), C.ens_Binary, C.ens_MasterOut);

    // Serialization
    typedef TwoPhaseStoreCL<StokesProblemT> SerializationT;
    SerializationT serializer(adap.GetMG(), Stokes, lset, 0, C.rst_Outputfile);

    if (C.vtk_VTKOut)
        vtkwriter.write();
    if (C.ens_EnsightOut)
        ensightwriter.write();

    const double Vol= 4./3.*M_PI*C.exp_RadDrop[0]*C.exp_RadDrop[1]*C.exp_RadDrop[2];
    double relVol = lset.GetVolume()/Vol;

    StokesSolverParamST instatStokesParam(C);
    StokesSolverFactoryCL<StokesProblemT, StokesSolverParamST> instatStokesSolverFactory(Stokes, instatStokesParam);
    StokesSolverBaseCL* instatStokesSolver= instatStokesSolverFactory.CreateStokesSolver();

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
    NSSolverT nssolver( Stokes, *instatStokesSolver, C.ns_Iter, C.ns_Tol, C.ns_Reduction);

    // coupling levelset NavStokes
    time.Reset();
    // Coupling Navier-Stokes with Levelset
//     typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
    typedef LinThetaScheme2PhaseCL <StokesProblemT, NSSolverT> CouplingT;
    CouplingT cpl( Stokes, lset, nssolver, C.stk_Theta, C.ns_Nonlinear, true);
//     typedef RecThetaScheme2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//     CouplingT cpl( Stokes, lset, nssolver, C.stk_theta, C.nonlinear, false);


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
                // / \todo (of) Ueberfluessig!
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

        cpl.DoStep( C.cpl_Iter);

        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cout << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec."<<std::endl;
            // Store measured values
            DROPS_LOGGER_SETVALUE("NavStokesCoupledLevelset", duration);
        }

        // Write out solution
        if (C.ens_EnsightOut && step%C.ens_EnsightOut==0)
            ensightwriter.write();
        if (C.vtk_VTKOut && step%C.vtk_VTKOut==0)
            vtkwriter.write();

        // Write out droplet information
        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        // Reparametrization of levelset function
        if (C.rpm_Freq && step%C.rpm_Freq==0)
        {
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster())
                std::cout << "\n==> Reparametrization\n"
                          << "- rel. Volume: " << relVol << std::endl;

            time.Reset();
            lset.ReparamFastMarching( C.rpm_Method, false, false, C.rpm_Method==3);
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cout << "- Reparametrization took "<<duration<<" sec."<<std::endl;
                DROPS_LOGGER_SETVALUE("Reparametrization",duration);
            }
        }
        if (ProcCL::IamMaster())
            std::cout << "- rel. Volume: " << relVol << std::endl;

        // Correction of the volume
        if (C.lvs_VolCorrection)
        {
            if (ProcCL::IamMaster())
                std::cout << "\n==> Adjust volume ...\n";
            time.Reset();
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            lset.Phi.Data+= dphi;
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cout << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< " took "<<duration<<" sec"<<std::endl;
                DROPS_LOGGER_SETVALUE("VolumeCorrection",duration);
            }
        }

        // Serialization
        if (C.rst_Serialization && step%C.rst_Serialization==0)
        {
            if (ProcCL::IamMaster())
                std::cout << "\n==> Serialize data ...\n";
            time.Reset();
            serializer.Write();
            time.Stop(); duration=time.GetMaxTime();

            if (ProcCL::IamMaster())
                std::cout << "- Serialization took " <<duration<< " sec"<<std::endl;
        }

        step_time.Stop();
        duration=step_time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cout <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
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
    SurfaceTensionCL::eps           = C.sft_JumpWidth;
    SurfaceTensionCL::lambda        = C.sft_RelPos;
    SurfaceTensionCL::sigma         = Stokes.GetCoeff().SurfTens;
    SurfaceTensionCL::sigma_dirt_fac= C.sft_DirtFactor;

//    instat_vector_fun_ptr gsigma= &(SurfaceTensionCL::grad_sm_step);
    LevelsetP2CL lset( mg, &SurfaceTensionCL::sigma_step, &SurfaceTensionCL::gsigma_step,
                       C.lvs_Theta, C.lvs_SD, -1, C.lvs_Iter, C.lvs_Tol, C.lvs_CurvDiff, C.rpm_NarrowBand);

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
        infofile = new std::ofstream( string(C.ens_EnsCase + ".info").c_str());
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    //Setup initial problem
    std::cout << "=================================================================================== Init:\n"
              << "==> Initialize Problem\n";
    InitProblemWithDrop(Stokes, lset, adapt.GetPMG());

    SolveCoupledNS(Stokes, lset, adapt);

    if (ProcCL::IamMaster()){
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
    std::cout << C << std::endl;

    DROPS::ParTimerCL alltime;
    SetDescriber();

    typedef DROPS::ZeroFlowCL                             CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::ParTimerCL time;

    // Create Geometry
    DROPS::MultiGridCL      *mg=0;
    DROPS::StokesBndDataCL  *bnddata=0;

    CreateGeom(mg, bnddata, Inflow, C.dmc_MeshFile, C.dmc_GeomType, C.dmc_BoundaryType, C.rst_Inputfile, C.exp_RadInlet);

    DROPS::ParMultiGridCL::Instance().AttachTo(*mg);
    DROPS::CheckParMultiGrid(DROPS::ParMultiGridCL::Instance());

    // Init problem
    DROPS::EllipsoidCL::Init( C.exp_PosDrop, C.exp_RadDrop );
    DROPS::AdapTriangCL adap( *mg, C.ref_Width, 0, C.ref_FinestLevel, ((C.rst_Inputfile == "none") ? 1 : -1)*C.ref_RefineStrategy);

    if (C.rst_Inputfile == "none")
        adap.MakeInitialTriang( DROPS::EllipsoidCL::DistanceFct);

    MyStokesCL prob( *mg, DROPS::ZeroFlowCL(C), *bnddata, DROPS::P1_FE, 0.0);

    Strategy( prob, adap);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cout);

    // Writeout logging-results
    if (DROPS::ProcCL::IamMaster()){
        // DROPS_LOGGER_WRITEOUT("TestBrickflowPar", (DROPS::ProcCL::Size()==1));

        // build specific name for logfile.
        char refLevel[32];
        sprintf(refLevel,"%d",C.ref_FinestLevel);
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


