
#include "Python.h"
#include "numpy/arrayobject.h"	/* NumPy header */

#include "poisson/poisson.h"
#include "poisson/integrTime.h"

#include <fstream>

#include "py_utils.hpp"
#include "py_param.hpp"
#include "functors.hpp"
#include "py_coeff_gradient.hpp"


static PyObject* drops_L2_scalar_prod(PyObject* self, PyObject* args)
{
  PyArrayObject* u=NULL;
  PyArrayObject* v=NULL;
  double lx, ly, lz, tmax;
  double theta;
  int h1_int;
  bool h1;

  // get the input arguments
  int pyarg_retval = PyArg_ParseTuple(args, "O!O!dddddi",
				      &PyArray_Type, &u,
				      &PyArray_Type, &v,
				      &lx, &ly, &lz, &tmax,
				      &theta,
				      &h1_int);
  if (!pyarg_retval) {
    PyErr_SetString(PyExc_TypeError, "The input arguments are not correct.\nThe correct function call is sp = scalarprod(u,v,lx,ly,lz,tmax,theta,h1)");
    return NULL;
  }
  h1 = (bool) h1_int;

  // dimensions
  int Nx_u, Nx_v;
  int Ny_u, Ny_v;
  int Nz_u, Nz_v;
  int Nt_u, Nt_v;
  Nx_u = u->dimensions[0];
  Nx_v = v->dimensions[0];
  Ny_u = u->dimensions[1];
  Ny_v = v->dimensions[1];
  Nz_u = u->dimensions[2];
  Nz_v = v->dimensions[2];
  Nt_u = u->dimensions[3];
  Nt_v = v->dimensions[3];
  if (Nx_u!=Nx_v || Ny_u!=Ny_v || Nz_u!=Nz_v || Nt_u!=Nt_v) {
    PyErr_SetString(PyExc_ValueError, "Vectors u and v do not have the same dimensions");
  }
  int nx = Nx_u-1, ny=Ny_u-1, nz=Nz_u-1, nt=Nt_u-1;

  DROPS::Point3DCL null(0.0);
  DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
  e1[0]= lx;
  e2[1]= ly;
  e3[2]= lz;

  DROPS::BrickBuilderCL brick(null, e1, e2, e3, nx, ny, nz);

  Zero zero = Zero();
  const bool isneumann[6]=
    { false, true,         // Gamma_in, Gamma_out
      true,  false,        // Gamma_h (wall), Gamma_r (surface)
      true,  true};       // Gamma_r, Gamma_r

  const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
    {zero, zero, zero, zero, zero, zero};
  DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);
  CoeffStatGradientCoeff pcl(h1);
  DROPS::PoissonP1CL<CoeffStatGradientCoeff> prob(brick, pcl, bdata, 0);

  DROPS::MultiGridCL& MG= prob.GetMG();
  DROPS::MLIdxDescCL& idx= prob.idx;
  DROPS::MLMatDescCL& A = prob.A;
  DROPS::MLMatDescCL& M = prob.M;

  idx.Set(1);
  prob.CreateNumbering(MG.GetLastLevel(), &idx);

  // Matrizen mit Index idx (Zeilen und Spalten)
  A.SetIdx(&idx, &idx);
  M.SetIdx(&idx, &idx);

  prob.SetupInstatSystem(A, M, prob.x.t);

  PyInterpolationFunctor u_fun(u, nx, ny, nz, nt, lx, ly, lz, tmax);
  PyInterpolationFunctor v_fun(v, nx, ny, nz, nt, lx, ly, lz, tmax);

  double scalar_prod = 0;
  if (nt==0) {
    scalar_prod = DROPS::Py_product(MG, idx, A, M, u_fun, v_fun, 0.0, h1);
  }
  else {
    double dt = tmax/nt, t;
    double *scalar_prod_t = new double[nt+1];
    for (int time_step=0; time_step<=nt; ++time_step) {
      t = dt*time_step;
      scalar_prod_t[time_step] = DROPS::Py_product(MG, idx, A, M, u_fun, v_fun, t, h1);
    }
    for (int time_step=0; time_step<nt; ++time_step) {
      scalar_prod += dt*((1-theta)*scalar_prod_t[time_step] + theta*scalar_prod_t[time_step+1]);
    }
    delete[] scalar_prod_t;
  }
  PyObject* retval = Py_BuildValue("d", scalar_prod);
  return retval;
}
