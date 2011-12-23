/// \file  prepy_product.cpp
/// \brief Setup problem for py_product function
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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany, AVT.PT RWTH Aachen
 */

// include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "geom/geomselect.h"
#include "misc/bndmap.h"                //include function container

// include numeric computing
#include "num/fe.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include "poisson/integrTime.h"
#include "num/poissonsolverfactory.h"

// include problem class
#include "misc/params.h"
#include "poisson/poissonCoeff.h"      // Coefficient-Function-Container poissonCoeffCL
#include "poisson/poisson.h"      // setting up the Poisson problem
#include "num/bndData.h"

// include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

// include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "num/parsolver.h"              // various parallel solvers
#include "num/parprecond.h"             // various parallel preconditioners
#endif

#include "pyconnect.h"
#include "pdefunction.h"
#include <math.h>

using namespace std;

typedef DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > PoissonProblem;

//const char line[] ="----------------------------------------------------------------------------------\n";

class PyScalarProductConnector {
public:
  void set_properties(const DROPS::ParamCL& P, PoissonProblem* prob_)
  {
    prob = prob_;
    nx = P.get<int>("DomainCond.nx");
    ny = P.get<int>("DomainCond.ny");
    nz = P.get<int>("DomainCond.nz");
    nt = P.get<int>("DomainCond.nt");
    double lx = P.get<double>("DomainCond.lx");
    double ly = P.get<double>("DomainCond.ly");
    double lz = P.get<double>("DomainCond.lz");
    double tmax = P.get<double>("DomainCond.tmax");
    dx = lx/(nx-1); dy = ly/(ny-1); dz = lz/(nz-1);
    dt = nt>1 ? tmax/(nt-1) : 1.0;
    //cout << "set_properties tmax, it = " << tmax << ","<< nt << ","<< dt <<","<< endl;
  }

  void SetProductFun(const PdeFunction* pdefun1_, const PdeFunction* pdefun2_) {
    pdefun1 = pdefun1_;
    pdefun2 = pdefun2_;
  }

  const PdeFunction* pdefun1;
  const PdeFunction* pdefun2;

  void getnum(const DROPS::Point3DCL& p, double t, int& ix, int& iy, int& iz, int& it)
  {
    ix=rd(p[0]/dx); iy=rd(p[1]/dy); iz=rd(p[2]/dz); it=rd(t/dt);
    //cout << "getnum t, it = " << t << "," << it << ","<< dt <<","<< endl;
  }
  PoissonProblem* prob;

  int nx, ny, nz, nt; // number of grid points
  double dx, dy, dz, dt;
private:
};

PyScalarProductConnector PySpC;

double fun1(const DROPS::Point3DCL& p, double t)
{
  int ix, iy, iz, it;
  PySpC.getnum(p, t, ix, iy, iz, it);
  (*(PySpC.pdefun1))(ix, iy, iz, it);
}
double fun2(const DROPS::Point3DCL& p, double t)
{
  int ix, iy, iz, it;
  PySpC.getnum(p, t, ix, iy, iz, it);
  PySpC.pdefun2->operator()(ix, iy, iz, it);
}


//DROPS::ParamCL P_sp; // Parameter object for scalar product - does this one really have to be global?

namespace DROPS {
  void ScalarProductSetup( PoissonProblem& Poisson, ParamCL& P)
  {
    // the triangulation
    MultiGridCL& mg= Poisson.GetMG();

    Poisson.idx.SetFE( P1_FE);                                  // set quadratic finite elements ??? P1 is quadratic ???
    //see class for explanation: template didnt work
    if ( PoissonSolverFactoryHelperCL().MGUsed(P))
      Poisson.SetNumLvl ( mg.GetNumLevel());
    Poisson.CreateNumbering( mg.GetLastLevel(), &Poisson.idx);  // number vertices and edges
    Poisson.b.SetIdx( &Poisson.idx);                            // tell b about numbering
    Poisson.x.SetIdx( &Poisson.idx);                            // tell x about numbering
    Poisson.A.SetIdx( &Poisson.idx, &Poisson.idx);              // tell A about numbering
    Poisson.M.SetIdx( &Poisson.idx, &Poisson.idx);              // tell M about numbering
    Poisson.U.SetIdx( &Poisson.idx, &Poisson.idx);

    /* std::vector<size_t> UnkOnProc( 1);
       UnkOnProc[0]  = Poisson.x.Data.size();
       IdxT numUnk   = Poisson.x.Data.size(),
       numAccUnk= Poisson.x.Data.size(); */

    Poisson.SetupInstatSystem( Poisson.A, Poisson.M, 0.0, 0 );
  }

}
/** Set up the matrices for calculating scalar products.
 *
 *  Takes number of grid points and lengths of the box
 */
int setup_sp_matrices(int nx, int ny, int nz, int nt, double lx, double ly, double lz, double tmax, bool h1)
{
  DROPS::ParamCL P;
  try
    {
      P.put<int>("DomainCond.nx", nx);
      P.put<int>("DomainCond.ny", ny);
      P.put<int>("DomainCond.nz", nz);
      P.put<int>("DomainCond.nt", nt);
      P.put<double>("DomainCond.lx", lx);
      P.put<double>("DomainCond.ly", ly);
      P.put<double>("DomainCond.lz", lz);
      P.put<double>("DomainCond.tmax", tmax);
      P.put<int>("DomainCond.RefineSteps", 0);
      stringstream MeshFile;
      MeshFile << lx << "x" << ly << "x" << lz << "@" << nx-1 << "x" << ny-1 << "x" << nz-1; // meshfile takes number of intervals, not grid points
      P.put<string>("DomainCond.MeshFile",MeshFile.str().c_str());
      if (h1) {
	P.put<string>("PoissonCoeff.Diffusion", "One");
      } else {
	P.put<string>("PoissonCoeff.Diffusion", "Zero");
      }
      P.put<string>("PoissonCoeff.Source", "Zero");
      P.put<string>("PoissonCoeff.Solution", "Zero");
      P.put<string>("PoissonCoeff.InitialVal", "Zero");
      P.put<string>("PoissonCoeff.Reaction", "Zero");
      P.put<int>("PoissonCoeff.Convection", 0);
      P.put<int>("Poisson.Method", 303);

      std::cout << P << std::endl;

      // set up data structure to represent a poisson problem
      // ---------------------------------------------------------------------
      std::cout << line << "Set up data structure to represent a Poisson problem ...\n";

      //create geometry
      DROPS::MultiGridCL* mg= 0;
      DROPS::PoissonBndDataCL* bdata = 0;

      int geomtype = 1; double r = 1; std::string serfile = "none";
      DROPS::BuildDomain( mg, MeshFile.str(), geomtype, serfile, r);

      std::string boundaryfuncs = "Zero!Zero!Zero!Zero!Zero!Zero";
      std::string boundarytype  = "0!21!21!0!21!21";
      DROPS::BuildBoundaryData( mg, bdata, boundarytype, boundaryfuncs); // here, memory is lost!

      // Setup the problem
      PoissonProblem* prob = new PoissonProblem( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), *bdata);
      ScalarProductSetup(*prob, P);

      PySpC.set_properties(P, prob);    //Initalize PythonConnector
    }  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

#include "pypdefunction.h"
using namespace boost::python;
numeric::array numpy_scalar_product(numeric::array& v, numeric::array& w) {
  typedef PdeFunction* PdeFunPtr;
  PdeFunction* vf = new PyPdeFunction(v);
  PdeFunction* wf = new PyPdeFunction(w);

  int nx=PySpC.nx, ny=PySpC.ny, nz=PySpC.nz, nt=PySpC.nt;
  assert(vf->get_dimensions(nx, ny, nz, nt));
  assert(wf->get_dimensions(nx, ny, nz, nt));

  PySpC.SetProductFun(vf, wf); //initalize two functions

  // croate output object
  npy_intp* output_dim = new npy_intp[1];
  output_dim[0] = nt;
  PyArrayObject* retval = (PyArrayObject*) PyArray_New(&PyArray_Type, 1, output_dim, PyArray_DOUBLE, NULL, NULL, 0, NPY_C_CONTIGUOUS, NULL);
  delete[] output_dim;
  object obj(handle<>((PyObject*)retval));
  double* solution_ptr = (double*)retval->data;

  PoissonProblem* prob = PySpC.prob;
  for (int timestep=0; timestep<nt; ++timestep) {
    solution_ptr[timestep] = DROPS::Py_product(prob->GetMG(), prob->idx, prob->A, prob->M, fun1, fun2, timestep*PySpC.dt, false);
    std::cout<<"The result of py_product in timestep " << timestep << " is "<<solution_ptr[timestep]<<std::endl;
  }
  array solution = extract<numeric::array>(obj);

  delete vf;
  delete wf;
  return solution;
}
