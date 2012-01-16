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
typedef boost::shared_ptr<PoissonProblem> PoissonProblemPtr;

using namespace boost::python;

class PyScalarProductConnector {
public:
  /** Constructor
   *
   *  Takes the number of grid points and lengths of the 4-dimensional box.
   */
  PyScalarProductConnector(int nx_, int ny_, int nz_, int nt_,
			   double lx_, double ly_, double lz_, double tmax_, bool h1_)
    : nx(nx_), ny(ny_), nz(nz_), nt(nt_),
      lx(lx_), ly(ly_), lz(lz_), tmax(tmax_)
  {
    //std::cout << "lx = " << lx << std::endl;
    dx = lx/(nx_-1);
    dy = ly/(ny_-1);
    dz = lz/(nz_-1);
    dt = nt_>1 ? tmax/(nt_-1) : 1.0;
    h1 = h1_;
    // TODO: timing would be nice ...
    setup_sp_matrices();
  }

  void setup_sp_matrices();

  numeric::array numpy_scalar_product(numeric::array v, numeric::array w) const;

private:
  int nx, ny, nz, nt; // number of grid points
  double lx, ly, lz, tmax;
  double dx, dy, dz, dt;
  bool h1;
  PoissonProblemPtr prob;
};


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

    Poisson.SetupInstatSystem( Poisson.A, Poisson.M, 0.0, 0 );
  }

}
/** Set up the matrices for calculating scalar products.
 *
 *  Takes number of grid points and lengths of the box
 */
void PyScalarProductConnector::setup_sp_matrices()
{
  try
    {
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
      //std::cout << P << std::endl;

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
      DROPS::BuildBoundaryData( mg, bdata, boundarytype, boundaryfuncs); // here, bdata is allocated, and we own it!

      // Setup the problem
      prob = PoissonProblemPtr(new PoissonProblem( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), *bdata));
      ScalarProductSetup(*prob, P);

      delete bdata;
    }  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

#include "pypdefunction.h"

numeric::array PyScalarProductConnector::numpy_scalar_product(numeric::array v, numeric::array w) const
{
  PdeFunction::ConstPtr vf(new PyPdeFunction(v));
  PdeFunction::ConstPtr wf(new PyPdeFunction(w));

  int tmp_nx=nx, tmp_ny=ny, tmp_nz=nz, tmp_nt=nt;
  if(!vf->get_dimensions(tmp_nx, tmp_ny, tmp_nz, tmp_nt)) {
    std::stringstream ss;
    ss << "wrong dimensions in first function. Should be " << nx << ", " << ny<< ", " << nz << ", " << nt<< ", but are in fact " << tmp_nx << ", " << tmp_ny<< ", " << tmp_nz << ", " << tmp_nt;
    throw (PyDropsErr(ss.str().c_str()));
  }
  if(!wf->get_dimensions(tmp_nx, tmp_ny, tmp_nz, tmp_nt)) {
    std::stringstream ss;
    ss << "wrong dimensions in second function. Should be " << nx << ", " << ny<< ", " << nz << ", " << nt<< ", but are in fact " << tmp_nx << ", " << tmp_ny<< ", " << tmp_nz << ", " << tmp_nt;
    throw (PyDropsErr(ss.str().c_str()));
  }

  DropsScalarProdFunction::ConstPtr fun1(new DropsScalarProdFunction(vf, dx, dy, dz, dt));
  DropsScalarProdFunction::ConstPtr fun2(new DropsScalarProdFunction(wf, dx, dy, dz, dt));

  // croate output object
  npy_intp* output_dim = new npy_intp[1];
  output_dim[0] = nt;
  PyArrayObject* retval = (PyArrayObject*) PyArray_New(&PyArray_Type, 1, output_dim, PyArray_DOUBLE, NULL, NULL, 0, NPY_C_CONTIGUOUS, NULL);
  delete[] output_dim;
  object obj(handle<>((PyObject*)retval));
  double* solution_ptr = (double*)retval->data;

  for (int timestep=0; timestep<nt; ++timestep) {
    solution_ptr[timestep] = DROPS::Py_product(prob->GetMG(), prob->idx, prob->A, prob->M, *fun1, *fun2, timestep*dt, false); //
    std::cout<<"The result of py_product in timestep " << timestep << " is "<<solution_ptr[timestep]<<std::endl;
  }
  array solution = extract<numeric::array>(obj);
  return solution;
}
