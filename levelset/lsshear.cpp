/// \file lsshear.cpp
/// \brief drop in shear flow
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt;  SC RWTH Aachen:

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
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/surfacetension.h"
#include <fstream>

const double      delta_t= 0.01;
const DROPS::Uint num_steps= 50;
const int         FPsteps= -1;

// du/dt - q*u - nu*laplace u + Dp = f - okn
//                          -div u = 0
//                               u = u0, t=t0

// Randdaten: x=0, x=1, y=0, y=1:  Dirichlet 0
//            z=0 und x<0.5        Neumann   0   (aus Impl.gruenden: Dir.)
//            z=0 und x>0.5        Inflow Dirichlet  parabol.
//            z=1 und x<0.5        Inflow Dirichlet  parabol.
//            z=1 und x>0.5        Neumann   0   (aus Impl.gruenden: Dir.)


// Tropfendaten:
DROPS::Point3DCL Mitte(0.5);
double           Radius= 0.2;

DROPS::SVectorCL<3> Parabol( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    if (p[0]<0.5)
        ret[2]= 4*p[0]*(p[0]-0.5);
    else
        ret[2]= 4*(1-p[0])*(p[0]-0.5);
    return ret;
}

double DistanceFct( const DROPS::Point3DCL& p)
{
    return (Mitte-p).norm()-Radius;
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }

namespace DROPS // for Strategy
{

class Uzawa_PCG_CL : public UzawaSolverCL<PCG_SsorCL>
{
  private:
    SSORPcCL   _ssor;
    PCG_SsorCL _PCGsolver;
  public:
    Uzawa_PCG_CL( MLMatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1., double omega=1.)
        : UzawaSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol, tau),
          _ssor( omega), _PCGsolver(_ssor, inner_iter, inner_tol)
        {}
};

class PSchur_GSPCG_CL: public PSchurSolverCL<PCG_SgsCL>
{
  private:
    SGSPcCL   _sgs;
    PCG_SgsCL _PCGsolver;
  public:
    PSchur_GSPCG_CL( MLMatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol)
        : PSchurSolverCL<PCG_SgsCL>( _PCGsolver, M, outer_iter, outer_tol),
          _PCGsolver( _sgs, inner_iter, inner_tol)
        {}
    PCG_SgsCL& GetPoissonSolver() { return _PCGsolver; }
};

template<class StokesProblemT>
void Strategy( StokesProblemT& Stokes, const BndDataCL<>& lsbnd, double inner_iter_tol)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();
    SurfaceTensionCL sf( sigmaf);
    LevelsetP2CL lset( MG, lsbnd, sf, 0.1);

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    VelVecDescCL* v= &Stokes.v;
    VecDescCL*    p= &Stokes.p;
    VelVecDescCL* b= &Stokes.b;
    VecDescCL* c= &Stokes.c;
    VelVecDescCL cpl_M;
    MLMatDescCL* A= &Stokes.A;
    MLMatDescCL* B= &Stokes.B;
    MLMatDescCL* M= &Stokes.M;
    MLMatDescCL prM;

    TimerCL time;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( MG.GetLastLevel(), pidx);
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    MG.SizeInfo( std::cout);
    b->SetIdx( vidx);
    c->SetIdx( pidx);
    cpl_M.SetIdx( vidx);
    p->SetIdx( pidx);
    v->SetIdx( vidx);
    std::cout << "Anzahl der Druck-Unbekannten: " << p->Data.size() << std::endl;
    std::cout << "Anzahl der Geschwindigkeitsunbekannten: " << v->Data.size() << std::endl;
    A->Reset();
    B->Reset();
    M->Reset();
    A->SetIdx(vidx, vidx);
    B->SetIdx(pidx, vidx);
    M->SetIdx(vidx, vidx);
    Stokes.N.SetIdx(vidx, vidx);
    prM.SetIdx( pidx, pidx);
    time.Reset();
    time.Start();
    Stokes.SetupSystem1( A, M, b, b, &cpl_M, lset, Stokes.t);
    Stokes.SetupSystem2( B, c, lset, Stokes.t);
    Stokes.SetupPrMass( &prM, lset);
    time.Stop();
    std::cout << time.GetTime() << " seconds for setting up all systems!" << std::endl;

    Stokes.InitVel( v, ZeroVel);
    lset.SetupSystem( Stokes.GetVelSolution(), delta_t);

    Uint meth;
    std::cout << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
    time.Reset();

    double outer_tol;
    std::cout << "tol = "; std::cin >> outer_tol;

    {
        std::cout << "Computing initial velocity..." << std::endl;
        PSchur_GSPCG_CL schurSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);
        schurSolver.Solve( A->Data, B->Data, v->Data, p->Data, b->Data, c->Data);
    }

    // Initialize Ensight6 output
    std::string ensf( "ensight/shear");
    Ensight6OutCL ensight( "shear.case", num_steps + 1);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   "shear flow field", ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",         ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",         ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",         ensf + ".vel", true));
    ensight.Write();
    typedef GMResSolverCL<SSORPcCL> LSetSolver;
    SSORPcCL ssorpc;
    LSetSolver gm( ssorpc, 100, 1000, 1e-7);
    LevelsetModifyCL lsetmod( 0, 0, 0, 0, 0, 0);
    if (meth)
    {
//        typedef PSchur_PCG_CL StokesSolverT;
        typedef PSchur_GSPCG_CL StokesSolverT;
        PSchur_GSPCG_CL StokesSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);
        typedef NSSolverBaseCL<StokesProblemT> SolverT;
        SolverT dummyFP( Stokes, StokesSolver);
        LinThetaScheme2PhaseCL<LSetSolver>
            cpl( Stokes, lset, dummyFP, gm, lsetmod, /*theta*/ 0.5, 0.5, /*nonlinear*/ 0.);
        cpl.SetTimeStep( delta_t);

        for (Uint step= 1; step<=num_steps; ++step)
        {
            std::cout << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( FPsteps);
            ensight.Write( step*delta_t);
        }
    }
    else // Uzawa
    {
        double tau;
        Uint inner_iter;
        tau=  0.5*delta_t;
        std::cout << "#PCG steps = "; std::cin >> inner_iter;
        typedef Uzawa_PCG_CL StokesSolverT;
        StokesSolverT uzawaSolver( prM.Data, 5000, outer_tol, inner_iter, inner_iter_tol, tau);
        typedef NSSolverBaseCL<StokesProblemT> SolverT;
        SolverT dummyFP( Stokes, uzawaSolver);
        LinThetaScheme2PhaseCL<LSetSolver>
            cpl( Stokes, lset, dummyFP, gm, lsetmod, /*theta*/ 0.5, 0.5, /*nonlinear*/ 0.);
        cpl.SetTimeStep( delta_t);

        for (Uint step= 1; step<=num_steps; ++step)
        {
            std::cout << "============= Schritt " << step << ":\n";
            cpl.DoStep( FPsteps);
            ensight.Write( step*delta_t);
        }
        std::cout << "Iterationen: " << uzawaSolver.GetIter()
                  << "\tNorm des Res.: " << uzawaSolver.GetResid() << std::endl;
    }

    std::cout << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=4)
    {
        std::cout << "You have to specify three parameters:\n\tlsshear <inner_iter_tol> <num_subdiv> <surf.tension>" << std::endl;
        return 1;
    }

    double inner_iter_tol= std::atof(argv[1]);
    int sub_div= std::atoi(argv[2]);
    sigma= std::atof(argv[3]);
    std::cout << "inner iter tol:  " << inner_iter_tol << std::endl;
    std::cout << "sub divisions:   " << sub_div << std::endl;
    std::cout << "surface tension: " << sigma << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    typedef DROPS::InstatNavierStokes2PhaseP2P1CL MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, sub_div, sub_div, sub_div);

    const bool IsNeumann[6]=
        {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel,  &Parabol, &Parabol };
    // parabol. Einstroembedingungen bei z=0 und z=1

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    DROPS::TwoPhaseFlowCoeffCL coeff( 1, 1, 1, 1, 0, DROPS::Point3DCL(0.));

    MyStokesCL prob(brick, coeff, DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    Strategy(prob, lsbnd, inner_iter_tol);
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cout << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
