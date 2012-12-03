#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include "../../PERIODICADDEDphysicaldata.hpp"
static PhysicalData pd;
#include "../../boxinterpLiangData.cpp"
#include "../createlevelLiang.cpp"
#include "../DiscreteLevelSetToDiscreteHeight.cpp"
#include "../heightinterpolateLiangDataSPLINES.cpp"

int main() {

      double* level; 
      double* level_tmp = createlevelLiang();
      level = DiscreteLevelSetToDiscreteHeight(level_tmp);
      double * xd = new double[pd.NX];
      for(int i=0; i<pd.NX; i++)
      {
            xd[i]=static_cast <double>(i)*pd.deltaX;
      }
      //Create spline:
      static gsl_interp_accel * acc;
      static gsl_spline * heightspline;
      acc = gsl_interp_accel_alloc();
      heightspline = gsl_spline_alloc(gsl_interp_cspline, pd.NX);
      gsl_spline_init(heightspline, xd, level, pd.NX);

double x=5.5*pd.NX;
double t=0.;
double height = HeightInterpolLiangData(x, t, acc, heightspline);
std::cout << height << std::endl;

return 0;
}
