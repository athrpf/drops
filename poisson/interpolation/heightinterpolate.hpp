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
double phaseinterpolate(double x, double t, double* y, const PeriodicData& pd);

double phaseinterpolate(double x, double t, gsl_interp_accel * acc, gsl_spline * y, const PeriodicData& pd);
