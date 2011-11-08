#include "py_interpolation.hpp"

void get_index(double val, double lb, double ub, int n, int& idx, double& frac)
{
  if (n==0)
    {
      idx = 0;
      frac = 0.0;
      return;
    }
  double dx = (ub-lb)/n;
  double frac1 = (val-lb)/dx;
  idx = (int)frac1;
  frac = frac1 - idx;

  if (idx==n && frac<1.e-10) {
    idx = n-1;
    frac = 1.;
  }
}

double bilinear_interpolation(double x, double y, const double* t, const PyArrayObject* ptr,
			      double xl, double yl, double tmax,
			      int nx, int ny, int nt)
{
  double x_frac, y_frac;
  int x0, y0;
  get_index(x, 0., xl, nx, x0, x_frac);
  get_index(y, 0., yl, ny, y0, y_frac);

  assert(x0>-1);
  assert(y0>-1);
  assert(x0<nx);
  assert(y0<ny);

  int x1, y1;
  x1 = x0 + 1;
  y1 = y0 + 1;

  double w1, w2;
  if (t) {
    int t0;
    double t_frac;
    get_index(*t, 0., tmax, nt, t0, t_frac);
    // assert(abs(t_frac) < 1.e-10)
    w1 = (1.-y_frac)* (*(double*)PyArray_GETPTR3(ptr, x0, y0, t0)) + y_frac*(*((double*)PyArray_GETPTR3(ptr, x0, y1, t0)));
    w2 = (1.-y_frac)*(*(double*)PyArray_GETPTR3(ptr, x1, y0, t0)) + y_frac*(*(double*)PyArray_GETPTR3(ptr, x1, y1, t0));
  }
  else {
    w1 = (1.-y_frac)*(*(double*)PyArray_GETPTR2(ptr, x0, y0)) + y_frac*(*(double*)PyArray_GETPTR2(ptr, x0, y1));
    w2 = (1.-y_frac)*(*(double*)PyArray_GETPTR2(ptr, x1, y0)) + y_frac*(*(double*)PyArray_GETPTR2(ptr, x1, y1));
  }

  return w1*(1.-x_frac) + w2*x_frac;
}

double trilinear_interpolation(const DROPS::Point3DCL& p, const double* t, const PyArrayObject* ptr,
			       double xl, double yl, double zl, double tmax,
			       int nx, int ny, int nz, int nt)
{
  double x_frac, y_frac, z_frac;
  int x0, y0, z0;
  get_index(p[0], 0., xl, nx, x0, x_frac);
  get_index(p[1], 0., yl, ny, y0, y_frac);
  get_index(p[2], 0., zl, nz, z0, z_frac);

  assert(x0>-1);
  assert(y0>-1);
  assert(z0>-1);
  assert(x0<nx);
  assert(y0<ny);
  assert(z0<nz);

  int x1, y1, z1;
  x1 = x0 + 1;
  y1 = y0 + 1;
  z1 = z0 + 1;

  double i1, i2, j1, j2, w1, w2, retval;
  if (t) {
    int t0;
    double t_frac;
    get_index(*t, 0., tmax, nt, t0, t_frac);
    assert(t0>-1);
    assert(t0<=nt);
    i1 = (1.-z_frac)*(*(double*)PyArray_GETPTR4(ptr, x0, y0, z0, t0)) + z_frac*(*(double*)PyArray_GETPTR4(ptr, x0, y0, z1, t0));
    i2 = (1.-z_frac)*(*(double*)PyArray_GETPTR4(ptr, x0, y1, z0, t0)) + z_frac*(*(double*)PyArray_GETPTR4(ptr, x0, y1, z1, t0));
    j1 = (1.-z_frac)*(*(double*)PyArray_GETPTR4(ptr, x1, y0, z0, t0)) + z_frac*(*(double*)PyArray_GETPTR4(ptr, x1, y0, z1, t0));
    j2 = (1.-z_frac)*(*(double*)PyArray_GETPTR4(ptr, x1, y1, z0, t0)) + z_frac*(*(double*)PyArray_GETPTR4(ptr, x1, y1, z1, t0));
  }
  else {
    i1 = (1.-z_frac)*(*(double*)PyArray_GETPTR3(ptr, x0, y0, z0)) + z_frac*(*(double*)PyArray_GETPTR3(ptr, x0, y0, z1));
    i2 = (1.-z_frac)*(*(double*)PyArray_GETPTR3(ptr, x0, y1, z0)) + z_frac*(*(double*)PyArray_GETPTR3(ptr, x0, y1, z1));
    j1 = (1.-z_frac)*(*(double*)PyArray_GETPTR3(ptr, x1, y0, z0)) + z_frac*(*(double*)PyArray_GETPTR3(ptr, x1, y0, z1));
    j2 = (1.-z_frac)*(*(double*)PyArray_GETPTR3(ptr, x1, y1, z0)) + z_frac*(*(double*)PyArray_GETPTR3(ptr, x1, y1, z1));
  }

  w1 = i1*(1.-y_frac) + i2*y_frac;
  w2 = j1*(1.-y_frac) + j2*y_frac;

  retval = w1*(1.-x_frac) + w2*x_frac;
  return retval;
}
