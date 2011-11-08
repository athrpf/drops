// Python (must be included as the very first)
#include "Python.h"
#include "numpy/arrayobject.h"	/* NumPy header */

// STL
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

// DROPS
#include "poisson/poisson.h"
#include "num/solver.h"

#include "py_source.hpp"
#include "py_utils.hpp"
#include "py_interpolation.hpp"
#include "functors.hpp"
#include "py_journalist.hpp"
#include "drops_utils.hpp"

#include "py_coeff_gradient.hpp"


static PyObject* drops_gradient_stat(PyObject *self, PyObject *args)
{
  PyArrayObject *u = NULL;
  PyArrayObject *psi = NULL;

  double lx, ly, lz;
  int nx, ny, nz;
  double tol;
  int maxiter;

  int h1_gradient_int = 1;
  int debug = 0;
  bool h1_gradient;

  int pyarg_retval = PyArg_ParseTuple(args,"O!O!dddiiid|ii",
				      &PyArray_Type, &u,
				      &PyArray_Type, &psi,
				      &lx, &ly, &lz, &nx, &ny, &nz,
				      &tol, &maxiter,
				      &h1_gradient_int, &debug);

  PrintLevel print_level = (debug ? JNLST_PRINTLEVEL_ALL : JNLST_PRINTLEVEL_NONE);
  Journalist jnlst(print_level);

  if (!pyarg_retval) {
    std::cout << "Inputs could not be parsed in gradient_l2_stat" << std::endl;
    return NULL;
  }
  h1_gradient = (bool) h1_gradient_int;
  if (h1_gradient) {
    jnlst << "Solving Gradient problem in H1-Norm\n";
  } else {
    jnlst << "Solving Gradient problem in L2-Norm\n";
  }
  int Nx = nx+1;
  int Ny = ny+1;
  int Nz = nz+1;
  int Nt = 1;

  // Check the dimensions of the input matrices.
  if (!check_dimensions_4(u, Nx, Ny, Nz, Nt, "u")) {return NULL;}
  if (!check_dimensions_4(psi, Nx, Ny, Nz, Nt, "psi")) {return NULL;}

  Zero zero;
  MassTransferBrick brick(nx, ny, nz, lx, ly, lz, zero, zero);

  //typedef DROPS::InstatPoissonP1CL<CoeffStatGradientCoeff> CoeffStatPoisson;
  //CoeffStatPoisson prob(brick.get_brick(), CoeffStatGradientCoeff(h1_gradient), brick.get_bdata());
  FixedValue alpha(h1_gradient ? 1.0 : 0.0);
  Zero f;
  Zero sta_coeff;
  Velocity v(1.,0.);
  PoissonCoeffCL pcl(alpha, f, v, sta_coeff);
  DROPS::PoissonP1CL<PoissonCoeffCL> prob(brick.get_brick(), pcl, brick.get_bdata());


  DROPS::MultiGridCL& MG= prob.GetMG();
  DROPS::MLIdxDescCL& idx= prob.idx;
  DROPS::VecDescCL& x= prob.x;
  DROPS::VecDescCL& b= prob.b;
  DROPS::MLMatDescCL& A= prob.A;
  DROPS::MLMatDescCL& M= prob.M;
  DROPS::VecDescCL dummy;
  DROPS::VecDescCL b1,b2;
  idx.Set(1, 0, 0, 0);

  // erzeuge Nummerierung zu diesem Index
  prob.CreateNumbering(MG.GetLastLevel(), &idx);

  // Vektoren mit Index idx
  x.SetIdx(&idx);
  b.SetIdx(&idx);
  dummy.SetIdx( &idx);
  b1.SetIdx(&idx);
  b2.SetIdx(&idx);
  // Matrizen mit Index idx (Zeilen und Spalten)
  A.SetIdx(&idx, &idx);
  M.SetIdx(&idx, &idx);
  prob.SetupInstatSystem(A, M, prob.x.t);

  PyInterpolationFunctor u_fun(u, nx, ny, nz, 0, lx, ly, lz, 0.0);
  PyInterpolationFunctor psi_fun(psi, nx, ny, nz, 0, lx, ly, lz, 0.0);

  // Set up RHS
  prob.SetupL2ProjGrad(b1, u_fun, psi_fun, NULL);
  prob.SetupInstatRhs( b2, b2, prob.t, dummy, prob.x.t); // Randwerte
  b.Data=b1.Data+b2.Data;

  //initialise solver
  DROPS::SSORPcCL pc(1.0);
  DROPS::PCG_SsorCL solver(pc, maxiter, tol);

  // Set up matrix
  DROPS::MatrixCL K;
  K.LinComb(1., A.Data, 1., M.Data);
  //solve
  solver.Solve(K, x.Data, b.Data);
  if (!py_check_convergence(solver, tol)) { return NULL;}

  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = nx+1;
  c_sol_dim[1] = ny+1;
  c_sol_dim[2] = nz+1;
  c_sol_dim[3] = 1;
  PyArrayObject* c_sol_py = (PyArrayObject*) PyArray_SimpleNew(4, c_sol_dim, PyArray_DOUBLE);
  SolutionContainer* sol_container = new
    PySolutionContainer(c_sol_py, nx, ny, nz, 0, lx, ly, lz, 0.0);
  sol_container->set_solution(prob, 0.0);

  A.Reset(); M.Reset();
  b.Reset();
  PyObject* retval = Py_BuildValue("N", PyArray_Return(c_sol_py));
  delete sol_container;
  return retval;
}

PyObject* drops_sensitivity_stat(PyObject *self, PyObject *args)
{
  PyArrayObject *u = NULL;
  PyArrayObject *a_tilde = NULL;

  double lx, ly, lz;
  int nx, ny, nz;
  double tol;
  int maxiter;

  int pyarg_retval = PyArg_ParseTuple(args,"O!O!dddiiidi",
				      &PyArray_Type, &u,
				      &PyArray_Type, &a_tilde,
				      &lx, &ly, &lz, &nx, &ny, &nz,
				      &tol, &maxiter);
  if (!pyarg_retval) {
    std::string msg = "Inputs could not be parsed in sensitivity_l2_stat";
    PyErr_SetString(PyExc_Exception, msg.c_str());
    return NULL;
  }

  int Nx = nx+1;
  int Ny = ny+1;
  int Nz = nz+1;
  int Nt = 1;

  // Check the dimensions of the input matrices.
  if (!check_dimensions_4(u, Nx, Ny, Nz, Nt, "u")) {return NULL;}
  if (!check_dimensions_4(a_tilde, Nx, Ny, Nz, Nt, "psi")) {return NULL;}

  Zero zero;
  MassTransferBrick brick(nx, ny, nz, lx, ly, lz, zero, zero);

  //typedef DROPS::InstatPoissonP1CL<CoeffStatGradientCoeff> CoeffStatPoisson;
  //CoeffStatPoisson prob(brick.get_brick(), CoeffStatGradientCoeff(true), brick.get_bdata());
  FixedValue alpha(1.0);
  Zero f;
  Zero sta_coeff;
  Velocity v(1.,0.);
  PoissonCoeffCL pcl(alpha, f, v, sta_coeff);
  DROPS::PoissonP1CL<PoissonCoeffCL> prob(brick.get_brick(), pcl, brick.get_bdata());

  DROPS::MultiGridCL& MG= prob.GetMG();
  DROPS::MLIdxDescCL& idx= prob.idx;
  DROPS::VecDescCL& x= prob.x;
  DROPS::VecDescCL& b= prob.b;
  DROPS::MLMatDescCL& A= prob.A;
  DROPS::MLMatDescCL& M= prob.M;

  DROPS::VecDescCL cplA, dummy;

  idx.Set(1, 0, 0, 0);

  // erzeuge Nummerierung zu diesem Index
  prob.CreateNumbering(MG.GetLastLevel(), &idx);

  // Vektoren mit Index idx
  x.SetIdx(&idx);
  b.SetIdx(&idx);
  cplA.SetIdx( &idx);
  dummy.SetIdx( &idx);
  // Matrizen mit Index idx (Zeilen und Spalten)
  A.SetIdx(&idx, &idx);
  M.SetIdx(&idx, &idx);
  prob.SetupInstatSystem(A, M, prob.x.t);
  prob.SetupInstatRhs( cplA, dummy, 0, dummy, 0);
  PyInterpolationFunctor u_fun(u, nx, ny, nz, 0, lx, ly, lz, 0.0);
  PyInterpolationFunctor a_tilde_fun(a_tilde, nx, ny, nz, 0, lx, ly, lz, 0.0);
  prob.SetupGradSrc(b, u_fun, a_tilde_fun, NULL);
  b.Data += cplA.Data;

  //initialise solver
  DROPS::SSORPcCL pc(1.0);
  DROPS::PCG_SsorCL solver(pc, maxiter, tol);
  //solve
  solver.Solve(A.Data, x.Data, b.Data);
  if (!py_check_convergence(solver, tol)) { return NULL;}

  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = Nx;
  c_sol_dim[1] = Ny;
  c_sol_dim[2] = Nz;
  c_sol_dim[3] = Nt;
  PyArrayObject* c_sol_py = (PyArrayObject*) PyArray_SimpleNew(4, c_sol_dim, PyArray_DOUBLE);
  SolutionContainer* sol_container = new
    PySolutionContainer(c_sol_py, nx, ny, nz, 0, lx, ly, lz, 0.0);
  sol_container->set_solution(prob, 0.0);

  A.Reset(); M.Reset();
  b.Reset();
  PyObject* retval = Py_BuildValue("N", PyArray_Return(c_sol_py));
  delete sol_container;
  return retval;
}
