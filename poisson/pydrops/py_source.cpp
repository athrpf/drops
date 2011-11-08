//********************************************************************************
// File:    py_source.cpp                                                         *
// Content: poisson with convection and SUPG stabilization + interface for python *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers, Liang Zhang*
//          IGPM RWTH Aachen                                                      *
//			Maka Karalashvili, Hans Pirnay                            *
//          AVT.PT                                                                *
// Version: 0.1                                                                   *
// History: begin - Jul, 01 2008 modified 08.2011                                 *
//********************************************************************************


// STL
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

// DROPS
#include "poisson/poisson.h"
#include "num/solver.h"
#include "poisson/integrTime.h"

//#include "matlabconnect.hpp"
#include "functors.hpp"
#include "drops_utils.hpp"
#include "py_coeff_dp_stat.hpp"
#include "py_journalist.hpp"

using std::string;

/**
   This code has the same functionality as poisson/ipfilm.cpp, but is extended by a matlab interface.
   It solves a convection-diffusion equation on a brick-shaped domain with certain bnd conditions
   modeling the (effective) heat/mass transfer in a flat film. Bnd conditions:
   - Dirichlet bc at x=0 (inflow)
   - Natural   bc at interface (Wall for heat and free surface for mass)
   - all other bnds: homogeneous natural bc

   The program can be called as a mex-function from matlab with the following syntax:

   [MaxIter, Csol] = source_dp( C0, B_in, F, alpha, B_interface, uN, amol, xl, yl, zl, nx, ny, nz, dt, nt, Theta, Tol, Iter,
   Flag_pr, Flag_bc, Flag_SUPG)

   scalar parameters:

   xl, yl, zl:     length of brick in x-/y-/z-direction [in mm]
   nx, ny, nz:     number of intervals in x-/y-/z-direction, i.e. Nxyz = (nx+1) x (ny+1) x (nz+1) grid points
   dt, nt:         length and number of time steps
   uN:             maximum velocity for Nusselt
   Theta:          parameter of one-step theta-scheme, controlling implicitness of time discretization
   Tol, Iter:      stopping criterion for the iterative solver (GMRES)
   a:          	  diffusion parameter a = lambda / (rho c) [in SI]
   //hflux:          bnd condition for heating: hflux = qh / (-lambda) [in SI]
   //                                       or  hflux = T_wall [K] if Dirichlet bc is used (cf. Flag)        //??
   Flag_pr:        to control if you solve the adjoint problem;
   Flag_bc:        to control the boundary type of the interface(in mass transport;
   Flag_SUPG:      to turn on stabilization, 1 with stabilization, o without.
   MaxIter (output):        number of iterations spent in the iterative solver (maximum over all time steps)

   matrix parameters:
   ------------------
   C0:     Nxyz x 1         initial temperature distribution T0(x,y,z)
   B_in:   Nyz  x Nt        temperature at inflow (x=0): T_in(y,z,t)
   F:      Nxyz x Nt        rhs term F(x,y,z,t)
   B_interface:             Interface boundary condition

   Csol (output):   Nxyz x Nt        solution of the conv-diff problem

   Here Nt = nt+1,    Nyz = (ny+1) x (nz+1),    Nxyz = (nx+1) x Nyz.

   As one can see, the columns of the matrices correspond to the different time steps.
   Note, that all the values have to be passed in entities W, mm
**/

namespace DROPS
{
  template<class Coeff>
  void InstatStrategy(Journalist& jnlst, DROPS::scalar_instat_fun_ptr& initial, InstatPoissonP1CL<Coeff>& Poisson, SolutionContainer* sol_container, int nt, double dt, double theta, double tol, int maxiter, int Flag_SUPG)
  {
    typedef InstatPoissonP1CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    IdxDescCL& idx= Poisson.idx;
    VecDescCL& x= Poisson.x;
    VecDescCL& b= Poisson.b;
    MatDescCL& A= Poisson.A;
    MatDescCL& M= Poisson.M;
    MatDescCL& U= Poisson.U;

    //SUPG
    MatDescCL& M_SD=Poisson.M_SD;
    MatDescCL& U_SD=Poisson.U_SD;
    VecDescCL& vM_SD=Poisson.vM_SD;
    VecDescCL& vU_SD=Poisson.vU_SD;

    idx.Set(1);
    // erzeuge Nummerierung zu diesem Index
    Poisson.CreateNumbering(MG.GetLastLevel(), &idx);

    // Vektoren mit Index idx
    x.SetIdx(&idx);
    b.SetIdx(&idx);
    // Matrizen mit Index idx (Zeilen und Spalten)
    A.SetIdx(&idx, &idx);
    M.SetIdx(&idx, &idx);
    U.SetIdx(&idx, &idx);

    //SUPG
    vM_SD.SetIdx(&idx);
    vU_SD.SetIdx(&idx);
    M_SD.SetIdx(&idx, &idx);
    U_SD.SetIdx(&idx, &idx);

    jnlst << "Number of Unknowns: " << (int)x.Data.size() << "\n";
    jnlst << "Theta: " <<  theta << "\n";
    jnlst << "Tolerance GMRES: " << tol <<"\t";
    jnlst << "max. Num. GMRES-Iterations: " << maxiter << "\n";

    // stationaerer Anteil
    Poisson.SetupInstatSystem(A, M, Poisson.t);

    //parameters about solver?
    SSORPcCL pc(1.0);
    typedef GMResSolverCL<SSORPcCL> SolverT;
    SolverT solver(pc, 500, maxiter, tol);

    // Zeitdiskretisierung mit one-step-theta-scheme
    // theta=1 -> impl. Euler
    // theta=0.5 -> Crank-Nicholson
    InstatPoissonThetaSchemeCL<MyPoissonCL, SolverT>
      ThetaScheme(Poisson, solver, theta, true, Flag_SUPG==1);  //first bool(convection), second bool(stabilization)
    ThetaScheme.SetTimeStep(dt);

    int MaxIter= 0;
    double MaxRes= 0., average= 0.;
    Poisson.Init(x, initial, 0);

    // main loop: return values?
    for (int step=1;step<=nt;step++)
      {
	ThetaScheme.DoStep(x);
	if (solver.GetResid()>tol) {
	  std::stringstream ss;
	  ss << "The solver tolerance could not be met.\n"
	     << "residual: " << solver.GetResid() << ", tol = " << tol
	     << "\nnumiter: " << solver.GetIter() << ", maxiter = " << maxiter;
	  throw(DROPS::DROPSErrCL(ss.str().c_str()));
	}
	average+= solver.GetIter();
	if (MaxIter<=solver.GetIter())
	  MaxIter= solver.GetIter();
	if (MaxRes<=solver.GetResid())
	  MaxRes= solver.GetResid();
	sol_container->set_solution(Poisson, Poisson.t);
      }
    average/= nt;
    jnlst << "Num. Iterations in average: " << average << "\n";
    jnlst << "max. Num. Iterations: " << MaxIter << "\n";
    jnlst << "max. res. Norm: " << MaxRes << "\n";
    jnlst << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    A.Reset(); M.Reset(); U.Reset();
    b.Reset(); U_SD.Reset(); M_SD.Reset();
    vU_SD.Reset(); vM_SD.Reset();
  }

} // end of namespace DROPS

// create Problem (Multigrid, BndData and Coeff) and call Strategy(...)
void source_dp(Journalist& jnlst, MassTransferBrick& brick, PoissonCoeffCL& pcl, SolutionContainer* sol_container, DROPS::scalar_instat_fun_ptr initial, int nt, double dt, double theta, double tol, int iter, int Flag_pr, int Flag_bc, int Flag_SUPG)
{
  DROPS::InstatPoissonP1CL<PoissonCoeffCL> prob(brick.get_brick(), pcl, brick.get_bdata(), Flag_pr & AdjFlagC);           //Adjoint problem
  DROPS::InstatStrategy(jnlst, initial, prob, sol_container, nt, dt, theta, tol, iter, Flag_SUPG);
  return;
}
