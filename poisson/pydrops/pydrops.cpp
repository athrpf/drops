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
#include "poisson/integrTime.h"

//#include "matlabconnect.hpp"
#include "functors.hpp"
#include "py_functors.hpp"
#include "py_utils.hpp"
#include "py_interpolation.hpp"
#include "drops_utils.hpp"
#include "py_coeff_dp_stat.hpp"
#include "py_journalist.hpp"

using std::string;

#include "py_source.hpp"
#include "py_coeff_dp_stat.hpp"
#include "py_kdelta_psi.hpp"

static PyObject* drops_source(PyObject *self, PyObject *args)
{
  PyArrayObject *c0 = NULL;
  PyArrayObject *b_in = NULL;
  PyArrayObject *F = NULL;
  PyArrayObject *py_alpha = NULL;
  PyArrayObject *b_interface = NULL;
  double u_N, amol, lx, ly, lz;
  int nx, ny, nz;
  double dt;
  int nt;
  double theta, tol;
  int maxiter, flag_pr, flag_bc, flag_supg;
  int debug = 0;

  int pyarg_retval = PyArg_ParseTuple(args, "O!O!O!O!O!dddddiiididdiiii|i",
				      &PyArray_Type, &c0,
				      &PyArray_Type, &b_in,
				      &PyArray_Type, &F,
				      &PyArray_Type, &py_alpha,
				      &PyArray_Type, &b_interface,
				      &u_N, &amol, &lx, &ly, &lz, &nx, &ny, &nz, &dt, &nt, &theta, &tol, &maxiter, &flag_pr, &flag_bc, &flag_supg, debug);
  if (!pyarg_retval) {
    std::stringstream ss;
    ss << "Inputs could not be parsed in source problem.\n"
       << "The function signature is\n\n"
      "\tsource(c0, b_in, F, alpha, b_interface, uN, amol, \n"
      "\t\tlx, ly, lz, nx, ny, nz, dt, nt, theta, tol, maxiter, flag_pr, flag_bc, flag_supg, debug=False)";
    PyErr_SetString(PyExc_Exception, ss.str().c_str());
    return NULL;
  }
  PrintLevel print_level = (debug ? JNLST_PRINTLEVEL_ALL : JNLST_PRINTLEVEL_NONE);
  Journalist jnlst(print_level);

  // Get the scalar input arguments.
  if (flag_pr) {
    u_N = - u_N;
  }

  int Nx   = nx+1;
  int Ny   = ny+1;
  int Nz   = nz+1;
  int Nt   = nt+1;

  // Check the dimensions of the input matrices.
  if (!check_dimensions_3(c0, Nx, Ny, Nz, "C0")) {return NULL;}
  if (!check_dimensions_4(F, Nx, Ny, Nz, Nt, "F")) {return NULL;}
  if (!check_dimensions_3(b_in, Ny, Nz, Nt, "b_in")) {return NULL;}
  if (!check_dimensions_3(b_interface, Nx, Nz, Nt, "b_in")) {return NULL;}
  if (!check_dimensions_4(py_alpha, Nx, Ny, Nz, Nt, "alpha")) {return NULL;}

  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = Nx;
  c_sol_dim[1] = Ny;
  c_sol_dim[2] = Nz;
  c_sol_dim[3] = Nt;
  PyArrayObject* c_sol_py = (PyArrayObject*) PyArray_SimpleNew(4, c_sol_dim, PyArray_DOUBLE);
  //  double* c_sol    = (double*) c_sol_py->data;

  // initialize to recognizable value to be able to catch uninitialized values
  PyArrayIterObject* c_sol_it = (PyArrayIterObject*)PyArray_IterNew((PyObject*)c_sol_py);
  while (c_sol_it->index<c_sol_it->size) {
    *((double*)c_sol_it->dataptr) = -1.2345;
    PyArray_ITER_NEXT(c_sol_it);
  }

  // set first timestep to initial value - watch out for adjoint problem specifics
  int timedim = 3;
  c_sol_it = (PyArrayIterObject*)PyArray_IterAllButAxis((PyObject*)c_sol_py, &timedim);
  if (!flag_pr) {
    while (c_sol_it->index<c_sol_it->size) {
      *((double*)c_sol_it->dataptr) = *((double*)PyArray_GetPtr(c0, c_sol_it->coordinates));
      PyArray_ITER_NEXT(c_sol_it);
    }
  }
  else { // adjoint problem!
    npy_intp* coord;
    while (c_sol_it->index<c_sol_it->size) {
      coord = c_sol_it->coordinates;
      *((double*)PyArray_GETPTR4(c_sol_py, coord[0], coord[1], coord[2], Nt-1)) = *((double*)PyArray_GETPTR3(c0, coord[0],coord[1],coord[2]));
      PyArray_ITER_NEXT(c_sol_it);
    }
  }

  if (!PyArray_ISCONTIGUOUS(c_sol_py)) {
    std::string msg = "Solution array is not contiguous!";
    PyErr_SetString(PyExc_Exception, msg.c_str());
    return NULL;}

  PyBoundaryInterpolationFunctor
    inlet_concentration(b_in, nx, ny, nz, nt,
			lx, ly, lz, dt*nt, 0, flag_pr);
  PyBoundaryInterpolationFunctor
    interface_concentration(b_interface, nx, ny, nz, nt,
			    lx, ly, lz, dt*nt, 3, flag_pr);
  PyInterpolationFunctor
    initial(c0, nx, ny, nz, 0, lx, ly, lz, 0.0, false, true);

  double tmax = dt*nt;
  PyInterpolationFunctor
    alpha_eff(py_alpha, nx, ny, nz, nt, lx, ly, lz, tmax, flag_pr);
  FixedValue alpha_mol(amol);
  PyInterpolationFunctor
    rhs(F, nx, ny, nz, nt, lx, ly, lz, tmax, flag_pr);
  AddFunctor alpha(1.0, alpha_mol, 1.0, alpha_eff);
  Velocity v(ly, u_N);
  Stability_Coeff sta_coeff(alpha, v, nx, lx);

  MassTransferBrick
    brick(nx, ny, nz, lx, ly, lz, inlet_concentration, interface_concentration);

  PoissonCoeffCL pcl(alpha, rhs, v, sta_coeff);

  SolutionContainer* sol_container = new
    PySolutionContainer(c_sol_py, nx, ny, nz, nt, lx, ly, lz, tmax, flag_pr);

  // Call the subroutine.
  try {
    source_dp(jnlst, brick, pcl, sol_container, initial, nt, dt, theta, tol, maxiter, flag_pr, flag_bc, flag_supg);
  }
  catch (DROPS::DROPSErrCL err) {
    std::stringstream ss;
    ss << "In source_dp: ";
    err.what(ss);
    PyErr_SetString(PyExc_Exception, ss.str().c_str());
    return NULL;
  }
  delete sol_container;
  PyObject* retval = Py_BuildValue("N", PyArray_Return(c_sol_py));
  return retval;
}

static PyObject* drops_coeff_stat(PyObject *self, PyObject *args)
{
  PyArrayObject *b_in = NULL;
  PyArrayObject *py_F = NULL;
  PyArrayObject *py_alpha = NULL;
  PyArrayObject *b_interface = NULL;
  double lx, ly, lz;
  int nx, ny, nz;
  double tol;
  int maxiter;
  int debug = 0;

  int pyarg_retval = PyArg_ParseTuple(args, "O!O!O!O!dddiiidi|i",
				      &PyArray_Type, &b_in,
				      &PyArray_Type, &py_F,
				      &PyArray_Type, &py_alpha,
				      &PyArray_Type, &b_interface,
				      &lx, &ly, &lz, &nx, &ny, &nz, &tol, &maxiter, &debug);
  if (!pyarg_retval) {
    std::cout << "Failed to read input arguments.\n";
    return NULL;
  }

  PrintLevel print_level = (debug ? JNLST_PRINTLEVEL_ALL : JNLST_PRINTLEVEL_NONE);
  Journalist jnlst(print_level);

  int Nx   = nx+1;
  int Ny   = ny+1;
  int Nz   = nz+1;
  int Nt   = 1;

  // Check the dimensions of the input matrices.
  if (!check_dimensions_4(py_F, Nx, Ny, Nz, Nt, "F")) {return NULL;}
  if (!check_dimensions_3(b_in, Ny, Nz, Nt, "b_in")) {return NULL;}
  if (!check_dimensions_3(b_interface, Nx, Nz, Nt, "b_in")) {return NULL;}
  if (!check_dimensions_4(py_alpha, Nx, Ny, Nz, Nt, "alpha")) {return NULL;}

  // bark if the diffusion coefficient is too low
  PyArrayIterObject* py_alpha_iter = (PyArrayIterObject*)PyArray_IterNew((PyObject*)py_alpha);
  //  PyArrayIterObject* py_F_iter = (PyArrayIterObject*)PyArray_IterNew((PyObject*)py_F);
  double minD = 100., dat=0.;
  while (py_alpha_iter->index<py_alpha_iter->size) {
    dat = *((double*)py_alpha_iter->dataptr);
    if (minD>dat)
      minD = dat;
    PyArray_ITER_NEXT(py_alpha_iter);
  }
  double min_allowed_D = 1e-7; // this is completely arbitrary
  if (minD<min_allowed_D) {
    std::stringstream ss;
    ss << "The diffusion coefficient has reached a value of " << minD << " and thereby fallen below the critical value " << min_allowed_D << "\nThe problem cannot be solved.\n";
    PyErr_SetString(PyExc_Exception, ss.str().c_str());
    return NULL;
  }

  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = Nx;
  c_sol_dim[1] = Ny;
  c_sol_dim[2] = Nz;
  c_sol_dim[3] = Nt;
  PyArrayObject* c_sol_py = (PyArrayObject*) PyArray_SimpleNew(4, c_sol_dim, PyArray_DOUBLE);
  //  double* c_sol    = (double*) c_sol_py->data;

  // initialize to recognizable value to be able to catch uninitialized values
  PyArrayIterObject* c_sol_it = (PyArrayIterObject*)PyArray_IterNew((PyObject*)c_sol_py);
  while (c_sol_it->index<c_sol_it->size) {
    *((double*)c_sol_it->dataptr) = -1.2345;
    PyArray_ITER_NEXT(c_sol_it);
  }

  if (!PyArray_ISCONTIGUOUS(c_sol_py)) {
    std::string msg = "Solution array is not contiguous!";
    PyErr_SetString(PyExc_Exception, msg.c_str());
    return NULL;
  }

  PyBoundaryInterpolationFunctor
    inlet_concentration(b_in, nx, ny, nz, 0, lx, ly, lz, 0.0, 0);
  PyBoundaryInterpolationFunctor
    interface_concentration(b_interface, nx, ny, nz, 0, lx, ly, lz, 0.0, 3);
  MassTransferBrick
    brick(nx, ny, nz, lx, ly, lz, inlet_concentration, interface_concentration);

  double tmax = 0.0;
  int nt = 0;
  PyInterpolationFunctor
    alpha_eff(py_alpha, nx, ny, nz, nt, lx, ly, lz, tmax);
  PyInterpolationFunctor
    rhs(py_F, nx, ny, nz, nt, lx, ly, lz, tmax);
  Velocity v(ly, 0.0);
  Stability_Coeff sta_coeff(alpha_eff, v, nx, lx);

  PoissonCoeffCL pcl(alpha_eff, rhs, v, sta_coeff);

  SolutionContainer* sol_container = new
    PySolutionContainer(c_sol_py, nx, ny, nz, nt, lx, ly, lz, tmax);

  // Call the subroutine.
  try {
    DROPS::coeff_stat(jnlst, brick, pcl, sol_container, tol, maxiter);
  }
  catch (DROPS::DROPSErrCL err) {
    std::stringstream ss;
    ss << "In coeff_stat: ";
    err.what(ss);
    PyErr_SetString(PyExc_Exception, ss.str().c_str());
    return NULL;
  }
  delete sol_container;

  PyObject* retval = Py_BuildValue("N", PyArray_Return(c_sol_py));
  return retval;
}


static PyObject* drops_kdelta_scalarprod(PyObject *self, PyObject *args)
{
  PyArrayObject *F = NULL;
  PyArrayObject *py_alpha = NULL;
  PyArrayObject *py_dalpha = NULL;
  PyArrayObject *py_u = NULL;
  double lx, ly, lz;
  int nx, ny, nz;
  double tol;
  int maxiter;
  int debug = 0;

  int pyarg_retval = PyArg_ParseTuple(args, "O!O!O!O!dddiiidi|i",
				      &PyArray_Type, &F,
				      &PyArray_Type, &py_alpha,
				      &PyArray_Type, &py_dalpha,
				      &PyArray_Type, &py_u,
				      &lx, &ly, &lz, &nx, &ny, &nz, &tol, &maxiter, &debug);
  if (!pyarg_retval) {
    std::cout << "Failed to read input arguments.\n";
    return NULL;
  }

  PrintLevel print_level = (debug ? JNLST_PRINTLEVEL_ALL : JNLST_PRINTLEVEL_NONE);
  Journalist jnlst(print_level);

  int Nx   = nx+1;
  int Ny   = ny+1;
  int Nz   = nz+1;
  int Nt   = 1;

  // Check the dimensions of the input matrices.
  if (!check_dimensions_4(F, Nx, Ny, Nz, Nt, "F")) {return NULL;}
  if (!check_dimensions_4(py_u, Nx, Ny, Nz, Nt, "u")) {return NULL;}
  if (!check_dimensions_4(py_dalpha, Nx, Ny, Nz, Nt, "dalpha")) {return NULL;}
  if (!check_dimensions_4(py_alpha, Nx, Ny, Nz, Nt, "alpha")) {return NULL;}

  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = Nx;
  c_sol_dim[1] = Ny;
  c_sol_dim[2] = Nz;
  c_sol_dim[3] = Nt;
  PyArrayObject* c_sol_py = (PyArrayObject*) PyArray_SimpleNew(4, c_sol_dim, PyArray_DOUBLE);
  //  double* c_sol    = (double*) c_sol_py->data;

  // initialize to recognizable value to be able to catch uninitialized values
  PyArrayIterObject* c_sol_it = (PyArrayIterObject*)PyArray_IterNew((PyObject*)c_sol_py);
  while (c_sol_it->index<c_sol_it->size) {
    *((double*)c_sol_it->dataptr) = -1.2345;
    PyArray_ITER_NEXT(c_sol_it);
  }

  if (!PyArray_ISCONTIGUOUS(c_sol_py)) {
    return NULL;
  }

  Zero zero;
  MassTransferBrick
    brick(nx, ny, nz, lx, ly, lz, zero, zero);

  double tmax = 0.0;
  int nt = 0;
  PyInterpolationFunctor
    alpha_eff(py_alpha, nx, ny, nz, nt, lx, ly, lz, tmax);
  PyInterpolationFunctor
    rhs(F, nx, ny, nz, nt, lx, ly, lz, tmax);
  Velocity v(ly, 0.0);
  //Stability_Coeff sta_coeff(alpha_eff, v, nx, lx);

  PoissonCoeffCL pcl(alpha_eff, rhs, v, zero);

  PyInterpolationFunctor
    dalpha(py_dalpha, nx, ny, nz, nt, lx, ly, lz, tmax);
  PyInterpolationFunctor
    u(py_u, nx, ny, nz, nt, lx, ly, lz, tmax);

  // Call the subroutine.
  double dretval;
  try {
    dretval = DROPS::compute_kdelta(jnlst, brick, pcl, dalpha, u, tol, maxiter);
  } catch (DROPS::DROPSErrCL err) {
    std::stringstream ss;
    ss << "In compute_kdelta: ";
    err.what(ss);
    PyErr_SetString(PyExc_Exception, ss.str().c_str());
    return NULL;
  }

  PyObject* retval = Py_BuildValue("d", dretval);
  return retval;
}

#include "py_coeff_gradient.cpp"
#include "py_L2_alt_scalar_prod.cpp"

static PyMethodDef DropsMethods[] = {
  {"source",  drops_source, METH_VARARGS, "source-inverse problem."},
  {"coeff_stat", drops_coeff_stat, METH_VARARGS,
   "stationary coefficient-inverse problem"},
  {"gradient_stat", drops_gradient_stat, METH_VARARGS,
   "Gradient of stationary coefficient-inverse problem"},
  {"sensitivity_stat", drops_sensitivity_stat, METH_VARARGS,
   "Sensitivity of stationary coefficient-inverse problem"},
  {"scalarprod", drops_L2_scalar_prod, METH_VARARGS,
   "Scalar product"},
  {"kdelta_scalarprod", drops_kdelta_scalarprod, METH_VARARGS,
   "K_delta applied to the adjoint solution psi"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initdrops(void)
{
  (void) Py_InitModule("drops", DropsMethods);
  import_array();		/* Initialize the Numarray module. */
  /* A segfault will occur if I use numarray without this.. */

}
