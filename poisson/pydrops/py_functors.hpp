//
/* This file holds the functors that need python objects */

#ifndef __PY_FUNCTORS_HPP__
#define __PY_FUNCTORS_HPP__

#include "drops_utils.hpp"
#include "functors.hpp"
#include "py_interpolation.hpp"

class PyInterpolationFunctor
{
private:
  PyArrayObject *u_;

  /** Number of intervals (not grid points!) in each dimension */
  int nx_, ny_, nz_, nt_;

  /** Length of the brick-shaped geometry in each dimension */
  double lx_, ly_, lz_, tmax_;

  /** Adjoint flag: If true, the time is inverted */
  bool invert_time_;
  bool initial_;

  /** Default constructor is not implemented */
  PyInterpolationFunctor();
public:
  /**
   * Constructor
   */
  PyInterpolationFunctor(PyArrayObject* u, int nx, int ny, int nz, int nt,
			 double lx, double ly, double lz, double tmax,
			 bool invert_time=false, bool initial=false)
  : u_(u),
    nx_(nx), ny_(ny), nz_(nz), nt_(nt),
    lx_(lx), ly_(ly), lz_(lz), tmax_(tmax),
    invert_time_(invert_time), initial_(initial)
  {}

  /** Evaluate the interpolation function */
  double operator()(const DROPS::Point3DCL& p, double t) const
  {
    if (invert_time_) {
      t = tmax_ - t;
    }
    double *t_ptr;
    t_ptr = &t;
    if (initial_) {
      t_ptr = NULL;
    }
    return trilinear_interpolation(p, t_ptr, u_, lx_, ly_, lz_, tmax_, nx_, ny_, nz_, nt_);
  }
};

class PyBoundaryInterpolationFunctor
{
private:
  PyArrayObject *u_;

  /** Number of intervals (not grid points!) in each dimension */
  int nx_, ny_, nz_, nt_;

  /** Length of the brick-shaped geometry in each dimension */
  double lx_, ly_, lz_, tmax_;

  /** Adjoint flag: If true, the time is inverted */
  bool invert_time_;

  /** indices of Point3D that will be used as indices - derived from interface_idx*/
  int idx1_, idx2_;
  int n1_, n2_;
  double l1_, l2_;

  /** Default constructor is not implemented */
  PyBoundaryInterpolationFunctor();
public:
  PyBoundaryInterpolationFunctor(PyArrayObject* u, int nx, int ny, int nz, int nt,
				 double lx, double ly, double lz, double tmax,
				 int interface_idx,
				 bool invert_time=false)
    : u_(u),
      nx_(nx), ny_(ny), nz_(nz), nt_(nt),
      lx_(lx), ly_(ly), lz_(lz), tmax_(tmax)
  {
    switch (interface_idx/2) {
    case 0:
      idx1_ = 1;
      idx2_ = 2;
      l1_   = ly_;
      l2_   = lz_;
      n1_   = ny_;
      n2_   = nz_;
      break;
    case 1:
      idx1_ = 0;
      idx2_ = 2;
      l1_   = lx_;
      l2_   = lz_;
      n1_   = nx_;
      n2_   = nz_;
      break;
    case 2:
      // this should not happen!
      throw(1);
      break;
    default:
      // this should not happen!
      throw(1);
    }
  }

  /** Evaluate the interpolation function */
  double operator()(const DROPS::Point3DCL& p, double t) const
  {
    if (invert_time_) {
      t = tmax_ - t;
    }
    double* t_ptr;
    t_ptr = &t;
    return bilinear_interpolation(p[idx1_], p[idx2_], t_ptr, u_, l1_, l2_, tmax_, n1_, n2_, nt_);

  }
};

class PySolutionContainer : public SolutionContainer
{
private:
  PyArrayObject* u_;
  int nx_, ny_, nz_, nt_;
  double lx_, ly_, lz_, tmax_;
  double dx_, dy_, dz_, dt_;
  bool invert_time_;

public:
  PySolutionContainer(PyArrayObject* u, int nx, int ny, int nz, int nt,
		      double lx, double ly, double lz, double tmax,
		      bool invert_time=false)
    :
    u_(u),
    nx_(nx), ny_(ny), nz_(nz), nt_(nt),
    lx_(lx), ly_(ly), lz_(lz), tmax_(tmax),
    invert_time_(invert_time)
  {
    dx_ = lx_/nx_;
    dy_ = ly_/ny_;
    dz_ = lz_/nz_;
    if (nt_>0) {
      dt_ = tmax_/nt_;
    } else {
      dt_ = 0.0;
    }
  }

  virtual void set_solution(const DROPS::PoissonP1CL<PoissonCoeffCL>& poisson, double t)
  {
    const DROPS::PoissonP1CL<PoissonCoeffCL>::const_DiscSolCL sol = poisson.GetSolution();
    if (!invert_time_) {
      DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
	{
	  const DROPS::Point3DCL p = sit->GetCoord();
	  if (t==0.0) { // necessary for the special case of stationary problems (there, dt=0)
	    *((double*)PyArray_GETPTR4(u_, rd(p[0]/dx_), rd(p[1]/dy_), rd(p[2]/dz_), 0)) = sol.val(*sit);
	  } else {
	    *((double*)PyArray_GETPTR4(u_, rd(p[0]/dx_), rd(p[1]/dy_), rd(p[2]/dz_), rd(t/dt_))) = sol.val(*sit);
	  }
	}
    }
    else {
      DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
	{
	  // Adjoint: time moves the other way...
	  const DROPS::Point3DCL p = sit->GetCoord();
	  *((double*)PyArray_GETPTR4(u_, rd(p[0]/dx_), rd(p[1]/dy_), rd(p[2]/dz_), rd((tmax_-t)/dt_))) = sol.val(*sit);
	}
    }
  }
};

#endif
