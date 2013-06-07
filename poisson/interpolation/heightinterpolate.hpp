// Jan Becker, AVT.PT, RWTH Aachen
// 2012-08-07

#pragma once

#include "periodicdata.hpp"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

/**
 * Functions for interpolating with phase shift.
 *
 * These functions take a point in space x and time t and interpolate the value at this point using either a grid
 * of values y or a spline.
 */
double phaseinterp1d(double x, double t, double* y, const PeriodicData& pd);
double phaseinterp2d(double x1, double x2, double t, double* y, const PeriodicData& pd);
double phaseinterp1d(double x, double t, gsl_interp_accel * acc, gsl_spline * y, const PeriodicData& pd);
