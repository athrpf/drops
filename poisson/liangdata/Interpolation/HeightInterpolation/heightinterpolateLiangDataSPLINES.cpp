//#include <cmath>
//#include <math.h>
//#include "../boxinterpLiangData.hpp"
//#include "../PERIODICADDEDphysicaldata.hpp"
//#include "heightinterpolateLiangData.hpp"
//static PhysicalData pd;
//#include <gsl/gsl_spline.h>
//#include <gsl/gsl_interp.h>
//#include <gsl/gsl_errno.h>

double HeightInterpolLiangData(double x, double t, gsl_interp_accel * acc, gsl_spline * heightspline) { 
  

  double L;//lenght of reference-domain
  L=(pd.NX-1)*pd.deltaX;//lenght of reference domain

  double p1=phase(x,t,pd.c);//phase
  double p2=phaseModuloWavelenghtLiangData(p1,L);//Calculate an equivalent phase, that lies in the reference domain.
  
  

    return gsl_spline_eval(heightspline, p2, acc);
}

