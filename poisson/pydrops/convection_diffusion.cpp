#include "poisson/poisson.h"
#include "num/solver.h"
#include "integrTime.h"
#include "out/output.h"
#include <fstream>
#include <string>
#include <sstream>

// Solve the convection-diffusion equation
// du/dt - nu*laplace u + Vel grad u + q*u = f

class PoissonCoeffCL
{

};

double Zero(const DROPS::Point3DCL&, double) { return 0.0; }


namespace DROPS
{

  template<class Coeff>
  void Strategy(InstatPoissonP1CL<Coeff>& Poisson, double theta, double tol, int maxiter, int Flag_SUPG)
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

    mexPrintf("Number of Unknowns: %d\n", x.Data.size());
    mexPrintf("Theta: %g\n", theta);
    mexPrintf("Tolerance GMRES: %g\t", tol);
    mexPrintf("max. Num. GMRES-Iterations: %d\n", maxiter);
    /*  Point3DCL pIF; pIF[1]=C.ly; // Punkt auf Phasengrenze
	const double vel= norm( PoissonCoeffCL::Vel(pIF,0)),
	cfl= vel*C.dt/(C.lx_/C.nx),
	Re= (9.81*1e3)*std::pow(C.ly,3)/(3*C.nu*C.nu), //
	Pr= C.nu/C.a_mol;                              //??
	mexPrintf("Geschwindigkeit Phasengrenze: %g [mm/s]\tentspricht CFL = %g\n", vel, cfl);
	mexPrintf("Geschwindigkeit Phasengrenze: %g [mm/s]\tentspricht CFL = %g\nRe = %g, Pr = %g\n", vel, cfl, Re, Pr);
	mexPrintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");*/

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
    ThetaScheme.SetTimeStep(C.dt_);

    int MaxIter= 0;
    double MaxRes= 0., average= 0.;
    Poisson.Init(x, PoissonCoeffCL::Initial, 0);

    for (int step=1;step<=C.nt_;step++)
      {
	ThetaScheme.DoStep(x);
	average+= solver.GetIter();
	if (MaxIter<=solver.GetIter())
	  MaxIter= solver.GetIter();
	if (MaxRes<=solver.GetResid())
	  MaxRes= solver.GetResid();
	MC.SetSol3D(Poisson.GetSolution(), Poisson.t);
      }
    average/= C.nt_;
    mexPrintf("Num. Iterations in average: %g\n", average);
    mexPrintf("max. Num. Iterations: %d\n", MaxIter);
    mexPrintf("max. res. Norm: %g\n", MaxRes);
    mexPrintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    A.Reset(); M.Reset(); U.Reset();
    b.Reset(); U_SD.Reset(); M_SD.Reset();
    vU_SD.Reset(); vM_SD.Reset();
  }

} // end of namespace DROPS

// Create Geometry, boundary condition, then call strategy
static void source_dp( double theta, double tol, int iter, int Flag_pr, int Flag_bc, int Flag_SUPG)
{ // create Problem (Multigrid, BndData and Coeff) and call Strategy(...)
  try
    {
      DROPS::Point3DCL null(0.0);
      DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
      e1[0]= C.lx_;
      e2[1]= C.ly_;
      e3[2]= C.lz_;

      typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL>
	InstatPoissonOnBrickCL;
      typedef InstatPoissonOnBrickCL MyPoissonCL;

      DROPS::BrickBuilderCL brick(null, e1, e2, e3, C.nx_, C.ny_, C.nz_);

      mexPrintf("Delta t = %g\n", C.dt_);
      mexPrintf("Num. of timesteps = %d\n", C.nt_);

      const bool InterfaceBCType= !(Flag_bc & DirBCFlagC);
      const bool isneumann[6]=
	{ false, true,              // Gamma_in, Gamma_out
	  true,  InterfaceBCType,   // Gamma_h (wall), Gamma_r (surface)
	  true,  true };            // Gamma_r, Gamma_r
      //Interface function creation
      const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
	{ &Inflow, &Zero, &Zero, InterfaceBCType ? &InterfaceFlux : &InterfaceValue, &Zero, &Zero};

      DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);
      MyPoissonCL prob(brick, PoissonCoeffCL(), bdata, Flag_pr & AdjFlagC);           //Adjoint problem
      DROPS::MultiGridCL& mg = prob.GetMG();

      //prepare MC
      MC.SetMG(&mg);
      MatlabConnectCL::ClearMaps();
      MatlabConnectCL::setFaceMap();
      MatlabConnectCL::setTetraMap();

      DROPS::Strategy(prob, theta, tol, iter, Flag_SUPG);

      return;
    }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}

//the main function
void convection_diffusion(IpfilmParamCL& param instat_scalar_fun* solution);
{

try
{
        DROPS::Point3DCL null(0.0);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= C.lx_;
        e2[1]= C.ly_;
        e3[2]= C.lz_;
        //create geometry
        DROPS::MultiGridCL* mg= 0;
        DROPS::PoissonBndDataCL* bdata = 0;

        //only for measuring cell, not used here
        double r = 1;
        std::string serfile = "none";
        //In poissonP1.cpp we use builddomain function
        DROPS::BrickBuilderCL brick(null, e1, e2, e3, C.nx_, C.ny_, C.nz_);

        const bool isneumann[6]=
        { false, true,              // Gamma_in, Gamma_out
          true,  false,             // Gamma_h (wall), Gamma_r (surface)
          true,  true };            // Gamma_r, Gamma_r
        const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { param.c_in_, &Zero, &Zero, param.c_surface_, &Zero, &Zero};

        DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);
        MyPoissonCL prob(brick, PoissonCoeffCL(), bdata, Flag_pr & AdjFlagC);           //Adjoint problem
        DROPS::MultiGridCL& mg = prob.GetMG();

        // Setup the problem
        DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), *bdata);


        // Refine the grid
        // Create new tetrahedra
        for ( int ref=1; ref <= param.refinement_; ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( *mg);
            mg->Refine();
        }
        mg->SizeInfo(cout);

        // Solve the problem
        DROPS::Strategy( prob);
        std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
        delete mg;
        delete bdata;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }

  return;
}
