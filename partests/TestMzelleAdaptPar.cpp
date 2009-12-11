/// \file TestMzelleAdaptPar.cpp
/// \brief Testing parallel solvers for the two-phase flow in the measurement cell with adaptive refinements
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
#include "levelset/surfacetension.h"

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
    T_create_geom,
    T_write_vel,
    T_create_ex,
    T_solve_NS_init
};

DROPS::TimeStoreCL Times(4);

void SetDescriber()
{
    Times.SetDescriber(T_create_geom,   "Writing geometry in quadrilateral grid");
    Times.SetDescriber(T_write_vel,     "Integrate velocity");
    Times.SetDescriber(T_create_ex,     "Create ExchangeCL");
    Times.SetDescriber(T_solve_NS_init, "Solving NavStokes for initial velocities");
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

template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    ParTimerCL time;
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
            schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
                               Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
            if (ProcCL::IamMaster())
                std::cout << "- Solving lin. Stokes ("<<step<<"): iter "<<schurSolver.GetIter()
                            <<", resid "<<schurSolver.GetResid()<<std::endl;
            ++step; iters+= schurSolver.GetIter();
        } while (schurSolver.GetIter() > 0);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
        {
            std::cout << "- Solving Stokes for initialization took "<<duration<<" sec, "
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
            reader.ReadVector( C.dmc_InitialFile+".vel", Stokes.v, Stokes.GetBndData().Vel);
            reader.ReadScalar( C.dmc_InitialFile+".pr",  Stokes.p, Stokes.GetBndData().Pr);
            reader.ReadScalar( C.dmc_InitialFile+".scl", lset.Phi, lset.GetBndData());
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
                cpl.Update();
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
    }
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
    SurfaceTensionDataCL::eps           = C.sft_JumpWidth;
    SurfaceTensionDataCL::lambda        = C.sft_RelPos;
    SurfaceTensionDataCL::sigma         = Stokes.GetCoeff().SurfTens;
    SurfaceTensionDataCL::sigma_dirt_fac= C.sft_DirtFactor;

    instat_scalar_fun_ptr sigmap  = 0;
    instat_vector_fun_ptr gsigmap = 0;
    if (C.sft_VarTension)
    {
        sigmap  = &SurfaceTensionDataCL::sigma_step;
        gsigmap = &SurfaceTensionDataCL::gsigma_step;
    }
    else
    {
        sigmap  = &SurfaceTensionDataCL::sigmaf;
        gsigmap = &SurfaceTensionDataCL::gsigma;
    }
    SurfaceTensionCL sf( sigmap, gsigmap);
    LevelsetP2CL lset( mg, sf, C.lvs_SD, C.lvs_CurvDiff, C.rpm_NarrowBand);
    if (C.sft_VarTension)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

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

    if (C.qlg_Quad){
        if (ProcCL::IamMaster())
            std::cout << "Write out geometry and velocity on the quadrilateral grid\n - geometry" << std::endl;
        time.Reset();
        LevelsetP2CL::const_DiscSolCL lset_sol=lset.GetSolution();
        QuadOutCL brickout( mg, &lset.idx, C.qlg_GridX, C.qlg_GridY, C.qlg_GridZ, C.qlg_Stepsize, C.qlg_Barycenter, C.qlg_Rotation);
        ParTimerCL timerQuad;
        timerQuad.Start();
        brickout.putGeom(C.qlg_FileName + std::string(".geo"));
        timerQuad.Stop();
        duration= timerQuad.GetTime();
        Times.AddTime(T_create_geom, duration);
        if (ProcCL::IamMaster())
            std::cout << " - (integrated) velocities"<<std::endl;
        timerQuad.Start();
        brickout.putVector(C.qlg_FileName + std::string(".vel_norm"),
                           C.qlg_FileName + std::string(".velY"),
                           C.qlg_FileName + std::string(".velZ"),
                           Stokes.GetVelSolution());
        timerQuad.Stop();
        duration= timerQuad.GetTime();
        Times.AddTime(T_write_vel, duration);
        time.Stop();
        duration=time.GetTime();

        if (ProcCL::IamMaster())
            std::cout << " => Took "<<duration<<" sec."<<std::endl;
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

    // Create Geometry
    DROPS::MultiGridCL      *mg=0;
    DROPS::StokesBndDataCL  *bnddata=0;
    CreateGeom(mg, bnddata, DROPS::Inflow, C.dmc_MeshFile, C.dmc_GeomType, C.dmc_BoundaryType, C.rst_Inputfile, C.exp_RadInlet);

    // Init problem
    DROPS::EllipsoidCL::Init( C.exp_PosDrop, C.exp_RadDrop );
    DROPS::AdapTriangCL adap( *mg, C.ref_Width, 0, C.ref_FinestLevel, C.ref_RefineStrategy);
    mg->SizeInfo(std::cout);

    adap.MakeInitialTriang( DROPS::EllipsoidCL::DistanceFct);

    MyStokesCL prob( *mg, DROPS::ZeroFlowCL(C), *bnddata);

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
        DROPS_LOGGER_WRITEOUT(basename.c_str(), (DROPS::ProcCL::Size()==1));
    }
    if (DROPS::ProcCL::IamMaster()){
        std::cout << "Statistics: TimeQuadGeom " << Times.GetTime(T_create_geom)
                  << " TimeQuadInt " << Times.GetTime(T_write_vel) << std::endl;
    }

    // free memory
    if (mg)      delete mg;
    if (bnddata) delete bnddata;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


