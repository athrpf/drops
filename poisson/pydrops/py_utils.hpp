
#ifndef __PY_UTILS_HPP__
#define __PY_UTILS_HPP__

#include "Python.h"
#include "numpy/arrayobject.h"	/* NumPy header */

#include <sstream>

#include "num/solver.h"

/** Check convergence of DROPS solver */
template <class SolverT>
bool py_check_convergence(DROPS::PCGSolverCL<SolverT> solver, double tol)
{
    double res_norm = solver.GetResid();
  if (res_norm > tol) {
    std::stringstream ss;
    ss << "The residual norm returned by the solver is " << res_norm << " but desired tolerance was " << tol << ". Exiting...\n";
    PyErr_SetString(PyExc_Exception, ss.str().c_str());
    return false;
  }
  return true;
}

bool check_dimensions_3(PyArrayObject* obj, int N1, int N2, int N3, const std::string name);

bool check_dimensions_4(PyArrayObject* obj, int N1, int N2, int N3, int N4, const std::string name);

#endif
