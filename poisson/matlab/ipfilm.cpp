/// \file ipfilm.cpp
/// \brief poisson with convection + interface for matlab
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt; SC RWTH Aachen:

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

#include "mex.h"

#include "../instatpoisson.h"
#include "../../num/solver.h"
#include "../integrTime.h"
#include "../../out/output.h"
#include <fstream>

const int AdjFlagC= 32,
          DirHeatBCFlagC= 16;

/**
This code has the same functionality as poisson/ipfilm.cpp, but is extended by a matlab interface.
It solves a convection-diffusion equation on a brick-shaped domain with certain bnd conditions
modeling the (effective) heat transfer in a flat film. Bnd conditions:
- Dirichlet bc at x=0 (inflow)
- Natural   bc at y=0 (heating)
- all other bnds: homogeneous natural bc

The program can be called as a mex-function from matlab with the following syntax:

  [MaxIter, Tsol] = ipfilm( T0, T_in, F, alpha, heat, xl, yl, zl, nx, ny, nz, dt, nt, Theta, Tol, Iter, Flag)

scalar parameters:
------------------
  xl, yl, zl:     length of brick in x-/y-/z-direction [in mm]
  nx, ny, nz:     number of intervals in x-/y-/z-direction, i.e. Nxyz = (nx+1) x (ny+1) x (nz+1) grid points
  dt, nt:         length and number of time steps
  Theta:          parameter of one-step theta-scheme, controlling implicitness of time discretization
  Tol, Iter:      stopping criterion for the iterative solver (GMRES)
  alpha:          diffusion parameter alpha = lambda / (rho c) [in SI]
  heat:           bnd condition for heating: heat = qh / (-lambda) [in SI]
                                         or  heat = T_wall [K] if Dirichlet bc is used (cf. Flag)
  Flag:           used for the adjoint problem and for the type of heating bc:
                  if Flag & AdjFlagC, then the adjoint operator is discretized.
                  if Flag & DirHeatBCFlagC, then heating is modeled by a Dirichlet bc (instead of natural bc).

  MaxIter (output):        number of iterations spent in the iterative solver (maximum over all time steps)

matrix parameters:
------------------
  T0:     Nxyz x 1         initial temperature distribution T0(x,y,z)
  T_in:   Nyz  x Nt        temperature at inflow (x=0): T_in(y,z,t)
  F:      Nxyz x Nt        rhs term F(x,y,z,t)

  Tsol (output):   Nxyz x nt        solution of the conv-diff problem,
                                    temperature distribution T(x,y,z,t) omitting initial time step

Here Nt = nt+1,    Nyz = (ny+1) x (nz+1),    Nxyz = (nx+1) x Nyz.

As one can see, the columns of the matrices correspond to the different time steps.

**/

extern void _main();

inline int rd( double d) { return static_cast<int>( d+0.5); } // rounding

class ParamCL
{
  public:
    double lx, ly, lz;
    int    nx, ny, nz, nt; // n=number of intervalls
    double dt, Heat, rho, mu, cp, alpha;
    std::string EnsDir, EnsCase;

    ParamCL()
      : lx(100), ly(0.3), lz(1), nx(8), ny(2), nz(2), nt(50), // in mm
        dt(0.02), Heat(5960), rho(866), mu(1.732e-3), cp(1500), alpha(0.26/1500/866), // in SI
        EnsDir("ensight"), EnsCase("FilmTemp")
      {}
} C;

class MatlabConnectCL
{ // holds the Matlab input matrices and Matlab output parameters
  private:
    int Ny, Nz, Nyz, Nxyz; // N=number of points
    double dx, dy, dz;

    double *T0, *T_in, *F,   // input  matrices: initial temp (Nxyz x 1),
                             //     inflow temp (Nyz x nt+1), rhs (Nxyz x nt+1)
           *T3D, *MaxIter;   // output matrices: temp solution (Nxyz x nt),
                             //     max. iterations of solver (1 x 1)
  public:
    int GetNum( double t)                  const { return rd(t/C.dt); }
    int GetNum( const DROPS::Point3DCL& p) const { return rd(p[0]/dx)*Nyz + rd(p[1]/dy)*Nz + rd(p[2]/dz); }
    template<int Bnd>
    int GetNumOnBnd( const DROPS::Point3DCL& p) const
    {
        switch(Bnd)
        {
            case 0: // y-z plane
                return rd(p[1]/dy)*Nz + rd(p[2]/dz);
            case 1: // x-z plane
                return rd(p[0]/dx)*Nz + rd(p[2]/dz);
            case 2: // x-y plane
                return rd(p[0]/dx)*Ny + rd(p[1]/dy);
            default:
                throw DROPS::DROPSErrCL( "MatlabConnectCL::GetNumOnBnd: unknown boundary!\n");
        }
    }

    double GetInitial( const DROPS::Point3DCL& p) const
      { return T0[GetNum(p)]; };
    double GetInflow( const DROPS::Point3DCL& p, double t) const
      { return T_in[GetNumOnBnd<0>(p) + GetNum(t)*Nyz]; };
    double GetRhs( const DROPS::Point3DCL& p, double t) const
      { return F[GetNum(p) + GetNum(t)*Nxyz]; };

  	template<class P1EvalT>
  	void SetSol3D( const P1EvalT& sol, double t)
  	{
  		const int num= (GetNum( t)-1)*Nxyz; // omit initial time step in output
  		double *out= T3D+num;

  		DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
  		{
  			out[GetNum( sit->GetCoord())]= sol.val( *sit);
  		}
  	}

  	void Init( const ParamCL& P, mxArray *plhs[], const mxArray *prhs[])
  	{
    	Ny= P.ny+1; Nz= P.nz+1; Nyz=Ny*Nz; Nxyz= Nyz*(P.nx+1);
    	dx= P.lx/P.nx; dy= P.ly/P.ny; dz= P.lz/P.nz;

		// Check to make sure the first input arguments are double matrices.
	    if(mxGetPi(prhs[0])!=NULL)
        	mexErrMsgTxt("Input T0 must be a double matrix.");
	    if(mxGetPi(prhs[1])!=NULL)
            mexErrMsgTxt("Input T_in must be a double matrix.");
	    if(mxGetPi(prhs[2])!=NULL)
	  	    mexErrMsgTxt("Input F must be a double matrix.");

	 	// Check the dimensions of the input matrices.
	 	if (mxGetM(prhs[0]) != Nxyz || mxGetN(prhs[0])!= 1)
	 	    mexErrMsgTxt("Input T0 has wrong dimensions.");
	 	if (mxGetM(prhs[1]) != Nyz  || mxGetN(prhs[1])!= P.nt+1)
	 	    mexErrMsgTxt("Input T_in has wrong dimensions.");
	 	if (mxGetM(prhs[2]) != Nxyz || mxGetN(prhs[2])!= P.nt+1)
	 	    mexErrMsgTxt("Input F has wrong dimensions.");

	 	// Get the matrix input arguments.
	 	T0=   mxGetPr(prhs[0]);
	 	T_in= mxGetPr(prhs[1]);
	 	F=    mxGetPr(prhs[2]);

	 	// Allocate memory for output arguments.
	 	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	 	plhs[1] = mxCreateDoubleMatrix(Nxyz, P.nt, mxREAL); // w/o initial time step

	 	// Set the output pointer to the output arguments.
	 	MaxIter = mxGetPr(plhs[0]);
	 	T3D =     mxGetPr(plhs[1]);
  	}
} MC;

// Solve the convection-diffusion equation
// du/dt - nu*laplace u + Vel grad u + q*u = f

class PoissonCoeffCL
{
  public:
    // static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double alpha(const DROPS::Point3DCL&, double)
      { return 1; }
    static double f(const DROPS::Point3DCL& p, double t)
    {
    	return MC.GetRhs(p,t);
    }
//    DROPS::P1EvalCL<double, const DROPS::NoBndDataCL<>, const DROPS::VecDescCL> f;
    static double Initial(const DROPS::Point3DCL& p, double t)
    {
        return MC.GetInitial(p);
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL& p, double)
    {
    	DROPS::Point3DCL ret;
        const double d= p[1]/C.ly,
            u= C.rho*9.81*C.ly*C.ly/2/C.mu*1e-3;
        ret[0]= u*(2-d)*d; // Nusselt
        return ret;
    }
};

double Zero(const DROPS::Point3DCL&, double) { return 0.0; }
double HeatFlux(const DROPS::Point3DCL&, double) { return C.Heat*1e-3; }
double HeatTemp(const DROPS::Point3DCL&, double) { return C.Heat; }

double Inflow(const DROPS::Point3DCL& p, double t) { return MC.GetInflow(p,t); }

namespace DROPS
{

template<class Coeff>
void Strategy(PoissonP1CL<Coeff>& Poisson,
  double theta, double tol, int maxiter, int EnsightFlag)
{
  typedef PoissonP1CL<Coeff> MyPoissonCL;

  MultiGridCL& MG= Poisson.GetMG();
  MLIdxDescCL& idx= Poisson.idx;
  VecDescCL& x= Poisson.x;
  VecDescCL& b= Poisson.b;
  MatDescCL& A= Poisson.A;
  MatDescCL& M= Poisson.M;
  MatDescCL& U= Poisson.U;

  idx.SetFE( P1_FE);
  // erzeuge Nummerierung zu diesem Index
  Poisson.CreateNumbering( MG.GetLastLevel(), &idx);

  // Vektoren mit Index idx
  x.SetIdx(&idx);
  b.SetIdx(&idx);
  // Matrizen mit Index idx (Zeilen und Spalten)
  A.SetIdx(&idx, &idx);
  M.SetIdx(&idx, &idx);
  U.SetIdx(&idx, &idx);

  mexPrintf("Anzahl der Unbekannten: %d\n", x.Data.size());
  mexPrintf("Theta: %g\n", theta);
  mexPrintf("Toleranz GMRES: %g\t", tol);
  mexPrintf("max. Anzahl GMRES-Iterationen: %d\n", maxiter);
  Point3DCL pIF; pIF[1]=C.ly; // Punkt auf Phasengrenze
  const double vel= norm( PoissonCoeffCL::Vel(pIF,0)),
      cfl= vel*C.dt/(C.lx/C.nx),
      Re= 9.81*std::pow(C.ly,3)*1e-9/3/C.mu/C.mu*C.rho*C.rho,
      Pr= C.mu/(C.alpha*C.rho);
  mexPrintf("Geschwindigkeit Phasengrenze: %g [mm/s]\tentspricht CFL = %g\nRe = %g, Pr = %g\n", vel, cfl, Re, Pr);
  mexPrintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");


  // stationaerer Anteil
  Poisson.SetupInstatSystem(A, M, Poisson.t);

  SSORPcCL pc(1.0);
  typedef GMResSolverCL<SSORPcCL> SolverT;
  SolverT solver(pc, 500, maxiter, tol);

  // Zeitdiskretisierung mit one-step-theta-scheme
  // theta=1 -> impl. Euler
  // theta=0.5 -> Crank-Nicholson
  InstatPoissonThetaSchemeCL<MyPoissonCL, SolverT>
    ThetaScheme(Poisson, solver, theta, true);
  const double nu= C.alpha*1e6;
  ThetaScheme.SetTimeStep(C.dt, nu);

  int MaxIter= 0;
  double MaxRes= 0., average= 0.;
  Poisson.Init(x, PoissonCoeffCL::Initial, 0);

  for (int step=1;step<=C.nt;step++)
  {
    ThetaScheme.DoStep(x);
    //mexPrintf("t= %g\n", Poisson.t);
    //mexPrintf("Iterationen: %d", solver.GetIter());
    //mexPrintf("    Norm des Residuums: %g\n", solver.GetResid());
    average+= solver.GetIter();
    if (MaxIter<=solver.GetIter())
      MaxIter= solver.GetIter();
    if (MaxRes<=solver.GetResid())
      MaxRes= solver.GetResid();
    MC.SetSol3D(Poisson.GetSolution(), Poisson.t);
  }
  average/= C.nt;
  mexPrintf("Anzahl Iterationen durchschnittlich: %g\n", average);
  mexPrintf("Anzahl Iterationen maximal: %d\n", MaxIter);
  mexPrintf("Norm des max. Residuums: %g\n", MaxRes);
  mexPrintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

  A.Reset(); M.Reset(); U.Reset();
  b.Reset();
}

} // end of namespace DROPS


static void ipdrops( double theta, double tol, int iter, int Flag)
{ // create Problem (Multigrid, BndData and Coeff) and call Strategy(...)
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= C.lx;
    e2[1]= C.ly;
    e3[2]= C.lz;

    typedef DROPS::PoissonP1CL<PoissonCoeffCL>
      InstatPoissonOnBrickCL;
    typedef InstatPoissonOnBrickCL MyPoissonCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, C.nx, C.ny, C.nz);

    mexPrintf("\nDelta t = %g", C.dt);
    mexPrintf("\nAnzahl der Zeitschritte = %d\n", C.nt);

    // bnd cond: x=0/lx, y=0/ly, z=0/lz
    const bool HeatBCType= !(Flag & DirHeatBCFlagC);
    const bool isneumann[6]=
      { false, true,   // Gamma_in, Gamma_out
        HeatBCType, true,    // Gamma_h (wall), Gamma_r (surface)
        true, true };  // Gamma_r, Gamma_r
    const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
      { &Inflow, &Zero, HeatBCType ? &HeatFlux : &HeatTemp, &Zero, &Zero, &Zero};

    DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata, Flag & AdjFlagC);
    DROPS::MultiGridCL& mg = prob.GetMG();

//    for (int count=1; count<=brick_ref; count++)
//    {
//      MarkAll(mg);
//      mg.Refine();
//    }
    // mg.SizeInfo(std::cout);
    DROPS::Strategy(prob, theta, tol, iter, Flag);

    return;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double Tol, theta;
  int Iter, Flag;

  // Check for proper number of arguments.
  if(nrhs!=17)
    mexErrMsgTxt("(T0, T_in, F, alpha, heat, xl, yl, zl, nx, ny, nz, dt, nt, Theta, Tol, Iter, Flag) as input required.");
  if(nlhs!=2)
    mexErrMsgTxt("Temperature field (3D instat.) and maximum number of GMRES-iterations as output required.");

  // Check to make sure the last input arguments are scalar.
  for (int index=3; index<nrhs; index++)
    if(!mxIsDouble(prhs[index]) ||
      mxGetN(prhs[index])*mxGetM(prhs[index])!=1)
    {
    	mexPrintf("Error in %dth input argument!\n", index+1);
        mexErrMsgTxt("Input must be a scalar.");
    }

  // Get the scalar input arguments.
  C.alpha = mxGetScalar(prhs[3]);
  C.Heat = mxGetScalar(prhs[4]);
  C.lx = mxGetScalar(prhs[5]);
  C.ly = mxGetScalar(prhs[6]);
  C.lz = mxGetScalar(prhs[7]);
  C.nx = rd( mxGetScalar(prhs[8]));
  C.ny = rd( mxGetScalar(prhs[9]));
  C.nz = rd( mxGetScalar(prhs[10]));
  C.dt = mxGetScalar(prhs[11]);
  C.nt = rd( mxGetScalar(prhs[12]));
  theta = mxGetScalar(prhs[13]);
  Tol =  mxGetScalar(prhs[14]);
  Iter = rd( mxGetScalar(prhs[15]));
  Flag = rd( mxGetScalar(prhs[16]));

  // Set the input matrices and output parameters.
  MC.Init( C, plhs, prhs);

  // Call the subroutine.
  ipdrops(theta, Tol, Iter, Flag);

  return;
}
