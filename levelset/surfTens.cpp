/// \file surfTens.cpp
/// \brief effect of surface tension
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "num/nssolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include <fstream>


DROPS::ParamMesszelleNsCL C;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

namespace DROPS // for Strategy
{

class PSchur_PCG_CL: public PSchurSolverCL<PCG_SsorCL>
{
  private:
    SSORPcCL   _ssor;
    PCG_SsorCL _PCGsolver;
  public:
    PSchur_PCG_CL( MLMatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double omega= 1.)
        : PSchurSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol),
          _ssor( omega), _PCGsolver(_ssor, inner_iter, inner_tol)
        {}
};

class ISPSchur_PCG_CL: public PSchurSolver2CL<PCGSolverCL<SSORPcCL>, PCGSolverCL<ISPreCL> >
{
  public:
    typedef PCGSolverCL<SSORPcCL> innerSolverT;
    typedef PCGSolverCL<ISPreCL>  outerSolverT;

  private:
    SSORPcCL     ssor_;
    innerSolverT innerSolver_;
    outerSolverT outerSolver_;

  public:
    ISPSchur_PCG_CL(ISPreCL& Spc, int outer_iter, double outer_tol,
                                  int inner_iter, double inner_tol)
        : PSchurSolver2CL<innerSolverT, outerSolverT>(
              innerSolver_, outerSolver_, outer_iter, outer_tol
          ),
          innerSolver_( ssor_, inner_iter, inner_tol),
          outerSolver_( Spc, outer_iter, outer_tol)
         {}
};

void Strategy( InstatNavierStokes2PhaseP2P1CL& Stokes, const LsetBndDataCL& lsbnd)
// flow control
{
    typedef InstatNavierStokes2PhaseP2P1CL StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    // Levelset-Disc.: Crank-Nicholson
    sigma= C.sft_SurfTension;
    SurfaceTensionCL sf( sigmaf, 0);
    LevelsetP2CL lset( MG, lsbnd, sf, C.lvs_SD, C.lvs_CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);
    const double Vol= EllipsoidCL::GetVolume();
    std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    MG.SizeInfo( std::cout);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";
    Stokes.A.SetIdx(vidx, vidx);
    Stokes.B.SetIdx(pidx, vidx);
    Stokes.M.SetIdx(vidx, vidx);
    Stokes.N.SetIdx(vidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx);
    Stokes.prA.SetIdx( pidx, pidx);

    Stokes.InitVel( &Stokes.v, ZeroVel);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    PSchur_PCG_CL   schurSolver( Stokes.prM.Data, C.stk_OuterIter, C.stk_OuterTol, C.stk_InnerIter, C.stk_InnerTol);
    VelVecDescCL curv( vidx);

    if (C.dmc_InitialCond != 0)
    {
        // solve stationary problem for initial velocities
        TimerCL time;
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cout << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        time.Reset();
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
            Stokes.v.Data, Stokes.p.Data, VectorCL( Stokes.b.Data + curv.Data), Stokes.c.Data);
        time.Stop();
        std::cout << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
    }
    curv.Clear();
    lset.AccumulateBndIntegral( curv);

    // Initialize Ensight6 output
    std::string ensf( C.ens_EnsDir + "/" + C.ens_EnsCase);
    Ensight6OutCL ensight( C.ens_EnsCase + ".case", C.tm_NumSteps + 1);
    ensight.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(),        "Maesszelle", ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),           "Levelset",   ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar( Stokes.GetPrSolution(),       "Pressure",   ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector( Stokes.GetVelSolution(),      "Velocity",   ensf + ".vel", true));
    ensight.Register( make_Ensight6Vector( Stokes.GetVelSolution( curv), "Curvature",  ensf + ".crv", true));

    ensight.Write();

    ISPreCL ispc( Stokes.prA.Data, Stokes.prM.Data, 1./C.tm_StepSize, C.stk_Theta);
    ISPSchur_PCG_CL ISPschurSolver( ispc,  C.stk_OuterIter, C.stk_OuterTol, C.stk_InnerIter, C.stk_InnerTol);
    ISPschurSolver.SetTol( C.stk_OuterTol);

    typedef NSSolverBaseCL<StokesProblemT> SolverT;
    SolverT navstokessolver(Stokes, ISPschurSolver);

    typedef GMResSolverCL<SSORPcCL> LsetSolverT;
    SSORPcCL ssorpc;
    GMResSolverCL<SSORPcCL> gm( ssorpc, 100, C.lvs_Iter, C.lvs_Tol);
    LevelsetModifyCL lsetmod( 0, 0, 0, 0, 0, 0);

    LinThetaScheme2PhaseCL<StokesProblemT, LsetSolverT>
        cpl( Stokes, lset, navstokessolver, gm, lsetmod, C.stk_Theta, C.lvs_Theta, 0.);
    cpl.SetTimeStep( C.tm_StepSize);

    for (int step= 1; step<=C.tm_NumSteps; ++step)
    {
        std::cout << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( C.cpl_Iter);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);

        ensight.Write( step*C.tm_StepSize);
        std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        if (C.lvs_VolCorrection)
        {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cout << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        }

        if (C.rpm_Freq && step%C.rpm_Freq==0)
        {
            lset.Reparam( C.rpm_Method);
            curv.Clear();
            lset.AccumulateBndIntegral( curv);

            ensight.Write( (step+0.1)*C.tm_StepSize);
            std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.lvs_VolCorrection)
            {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cout << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc>2)
    {
        std::cout << "You have to specify at most one parameter:\n\t"
                  << argv[0] << " [<param_file>]" << std::endl;
        return 1;
    }
    std::ifstream param;
    if (argc>1)
        param.open( argv[1]);
    else
        param.open( "surfTens.param");
    if (!param)
    {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cout << C << std::endl;

    typedef DROPS::InstatNavierStokes2PhaseP2P1CL    MyStokesCL;

    const double L= 3e-3;
    DROPS::Point3DCL orig(-L), e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 2*L;

    const int n= std::atoi( C.dmc_MeshFile.c_str());
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, n, n, n);

    const DROPS::BndCondT bc[6]=
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel};

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    MyStokesCL prob(builder, DROPS::TwoPhaseFlowCoeffCL(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun));

    DROPS::MultiGridCL& mg = prob.GetMG();

    DROPS::EllipsoidCL::Init( C.exp_PosDrop, C.exp_RadDrop);
    for (int i=0; i<C.ref_FinestLevel; ++i)
    {
        DROPS::MarkInterface( DROPS::EllipsoidCL::DistanceFct, C.ref_Width, mg);
        mg.Refine();
    }
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;

    DROPS::Strategy( prob, lsbnd);  // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cout << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

