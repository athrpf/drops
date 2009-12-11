/// \file ipfilm.cpp
/// \brief test program for the instat. poisson-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "poisson/instatpoisson.h"
#include "poisson/integrTime.h"
#include "geom/multigrid.h"
#include "num/solver.h"
#include "out/ensightOut.h"

class ParamCL
{
  public:
    double dx, dy, dz;
    int    nx, ny, nz;
    double Heat;
    double rho, mu, cp, lambda;
    std::string EnsDir, EnsCase;

    ParamCL()
      : dx(100), dy(0.3), dz(1), nx(8), ny(2), nz(2),  // in mm
        Heat(5960), rho(866), mu(1.732e-3), cp(1500), lambda(0.26), // in SI
        EnsDir("ensight"), EnsCase("FilmTemp")
      {}
} C;

// du/dt - nu*laplace u + Vel grad u + q*u = f

class PoissonCoeffCL
{
  public:
    // static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double alpha(const DROPS::Point3DCL&, double)
      { return 1; }
    static double f(const DROPS::Point3DCL& p, double t)
    {
//        return 0;
        const double u= C.rho*9.81*C.dy*C.dy/2/C.mu*1e-3;
        return std::cos((p[0] + t*u)/C.dx*2*M_PI);
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL& p, double)
    {
        DROPS::Point3DCL ret;
        const double d= p[1]/C.dy,
            u= C.rho*9.81*C.dy*C.dy/2/C.mu*1e-3;
        ret[0]= u*(2-d)*d; // Nusselt
        return ret;
    }
};


double Zero(const DROPS::Point3DCL&, double) { return 0.0; }
double Heat(const DROPS::Point3DCL&, double) { return C.Heat/C.lambda*1e-3; }


namespace DROPS
{

template<class Coeff>
void Strategy(InstatPoissonP1CL<Coeff>& Poisson, double dt, int time_steps,
  double nu, double theta, double tol, int maxiter)
{
  typedef InstatPoissonP1CL<Coeff> MyPoissonCL;

  MultiGridCL& MG= Poisson.GetMG();
  MLIdxDescCL& idx= Poisson.idx;
  VecDescCL&   x= Poisson.x;
  VecDescCL&   b= Poisson.b;
  MLMatDescCL& A= Poisson.A;
  MLMatDescCL& M= Poisson.M;
  MLMatDescCL& U= Poisson.U;

  idx.SetFE( P1_FE);
  Poisson.CreateNumbering( MG.GetLastLevel(), &idx);

  x.SetIdx(&idx);
  b.SetIdx(&idx);
  A.SetIdx(&idx, &idx);
  M.SetIdx(&idx, &idx);
  U.SetIdx(&idx, &idx);

  std::cout << "Anzahl der Unbekannten: " <<  Poisson.x.Data.size()
    << std::endl;
  Poisson.SetupInstatSystem(A, M, 0);

  SSORPcCL pc(1.0);
  typedef GMResSolverCL<SSORPcCL> SolverT;
  SolverT solver(pc, 50, maxiter, tol);
  InstatPoissonThetaSchemeCL<MyPoissonCL, SolverT>
    ThetaScheme(Poisson, solver, theta, true);

  ThetaScheme.SetTimeStep(dt, nu);

  Ensight6OutCL  ens(C.EnsCase+".case", time_steps+1);
  const std::string filename= C.EnsDir + "/" + C.EnsCase;
  ens.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(), "Film",       filename + ".geo"));
  ens.Register( make_Ensight6Scalar( Poisson.GetSolution(), "Temperatur", filename + ".tp", true));
  ens.Write();

  double average= 0.0;
  for (int step=1;step<=time_steps;step++)
  {
    ThetaScheme.DoStep(x);
    std::cout << "t= " << Poisson.t << std::endl;
    std::cout << "Iterationen: " << solver.GetIter()
      << "\tNorm des Residuums: " << solver.GetResid() << std::endl;
    average+= solver.GetIter();
    ens.Write( step*dt);
  }
  average/= time_steps;
  std::cout << "Anzahl der Iterationen im Durchschnitt: " << average
            << std::endl;
}

} // end of namespace DROPS



int main()
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= C.dx;
    e2[1]= C.dy;
    e3[2]= C.dz;

    typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL>
      InstatPoissonOnBrickCL;
    typedef InstatPoissonOnBrickCL MyPoissonCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, C.nx, C.ny, C.nz);

    double dt= 0.0;
    int time_steps= 0, brick_ref= 0;

    std::cout << "\nDelta t = "; std::cin >> dt;
    std::cout << "\nAnzahl der Zeitschritte = "; std::cin >> time_steps;
    std::cout << "\nAnzahl der Verfeinerungen = "; std::cin >> brick_ref;

    // bnd cond: x=0/dx, y=0/dy, z=0/dz
    const bool isneumann[6]=
      { false, true, true, true, true, true };
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
      { &Zero, &Zero, &Heat, &Zero, &Zero, &Zero};

    DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();

    for (int count=1; count<=brick_ref; count++)
    {
      MarkAll(mg);
      mg.Refine();
    }
    mg.SizeInfo(std::cout);

    // Diffusionskoeffizient
    //double nu= 6.303096;
    double nu= C.lambda/C.rho/C.cp*1e6;
//    double nu= 1e-3;
    // one-step-theta scheme
    double theta= 1;

    // Daten fuer den Loeser
    double tol= 1.0e-10;
    int maxiter= 5000;

    DROPS::Strategy(prob, dt, time_steps, nu, theta, tol, maxiter);

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}



