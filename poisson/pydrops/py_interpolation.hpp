
#ifndef __PY_INTERPOLATION__
#define __PY_INTERPOLATION__

#include "Python.h"
#include "numpy/arrayobject.h"	/* NumPy header */

#include "poisson/poisson.h"

/** Find the grid indices between which an interpolation point lies.
 *  n is the number of intervals, not the number of grid points!
 *  Returns the index for the grid point right before val in idx.
 *  frac is the fraction of the interval, after which val lies:
 *  val = value(idx) + frac*(value(idx+1) - value(idx))*/
void get_index(double val, double lb, double ub, int n, int& idx, double& frac);

/** Bilinear interpolation for boundary grid points */
double bilinear_interpolation(double x, double y, const double* t, const PyArrayObject* ptr,
			      double xl, double yl, double tmax,
			      int nx, int ny, int nt);

/** Trilinear interpolation for grid points (see Wikipedia...):
 * Return the value of ptr for the grid.
 * If t=NULL, the time variable is omitted.
 */
double trilinear_interpolation(const DROPS::Point3DCL& p, const double* t, const PyArrayObject* ptr,
			       double xl, double yl, double zl, double tmax,
			       int nx, int ny, int nz, int nt);

#endif
