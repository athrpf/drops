#include <cmath>
#include "boxinterp.hpp"
#include "periodicdata.hpp"
#include "heightinterpolate.hpp"

double get_equivalent_phase(double x, double t, const PeriodicData& pd)
{
  double L = (pd.NX-1)*pd.deltaX; // lenght of reference domain
  double p1 = phase(x,t,pd.c); // Phase
  double p2 = phaseModuloWavelenghtLiangData(p1,L); // Calculate an equivalent phase, that lies in the reference domain.
  return p2;
}

double phaseinterpolate(double x, double t, double* h, const PeriodicData& pd)
{
  double p2 = get_equivalent_phase(x, t, pd);
  int nX = boxnumber(p2,pd.deltaX,0); // Find the relevant boxnumber for p2
  if(nX==pd.NX) nX=(nX-1); // we sit directly on the last grid point - step back
  double x0 = (nX-1)*pd.deltaX;
  return ((h[nX]-h[nX-1])/pd.deltaX)*(p2-x0) + h[nX-1];
}

double phaseinterpolate(double x, double t, gsl_interp_accel * acc, gsl_spline * heightspline, const PeriodicData& pd)
{
  double p2 = get_equivalent_phase(x, t, pd);
  return gsl_spline_eval(heightspline, p2, acc);
}
