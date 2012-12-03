#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "../../liangdata/Interpolation/PERIODICADDEDphysicaldata.hpp"
static PhysicalData pd;
#include "../../liangdata/Interpolation/boxinterpLiangData.cpp"

#include "../../liangdata/Interpolation/HeightInterpolation/createlevelLiang.cpp"

#include "../../liangdata/Interpolation/HeightInterpolation/DiscreteLevelSetToDiscreteHeight.cpp"

#include "../../liangdata/Interpolation/VelInterpolation/functionjaninterpolateLiangData.cpp"

namespace Jan {

    
    static bool first_call = true;
    static double* u;
    static double* v;
    static double* w;
    static double* level;

  static double Dmol;

    void setup()
    {
      //Create an array of Liang's level-set-values, which is converted to
      //an array of discrete height-values by use of the
      //function "DiscreteLevelSetToDiscreteHeight":
      double* level_tmp = createlevelLiang();
      level = DiscreteLevelSetToDiscreteHeight(level_tmp);
      //Create array of discrete reference-points on x-axis that are, in addition to the array "level" needed to 
      //create a spline for the reference-heigt-profile: 
      double * xd = new double[pd.X];
      for(int i=0; i<pd.NX; i++)
      {
            xd[i]=i*pd.deltaX;
      }
      //Create spline:
      static gsl_spline ** heightspline;
      static gsl_interp_accel * acc;
      acc = gsl interp_accel_alloc();
      heightspline = new gsl_spline;
      heightspline = gsl_splie_alloc(gsl_interp_cspline, pd.NX);
      gsl_spline_init(heightspline, xd, level, pd.NX);
      
      // read data of the velocity-field (and the level-set)
      int NumCoords = pd.NX*pd.NY;
      u = new double[NumCoords];
      v = new double[NumCoords];
      w = new double[NumCoords];
      std::string velocity_filename = "../../liangdata/Interpolation/DataForPoissonCoeff/VelocityWithLevelSetFinal.txt";
      std::ifstream ufile;
      ufile.open(velocity_filename.c_str(), std::fstream::in);
      double curr_u, curr_v, curr_w, curr_level;
      for (int k=0; k<NumCoords; k++)
      {
         ufile >> curr_u >> curr_v >> curr_w >> curr_level;
         u[k] = curr_u;
         v[k] = curr_v;
         w[k] = curr_w;
      }
       first_call = false;
    }
#include "heightinterpolateLiangDataSPLINES.cpp"
double HeightInterpolLiangData(double x, double t) { 

  double L;//lenght of reference-domain
  L=(pd.NX-1)*pd.deltaX;//lenght of reference domain

  double p1=phase(x,t,pd.c);//phase
  double p2=phaseModuloWavelenghtLiangData(p1,L);//Calculate an equivalent phase, that lies in the reference domain.
  
  

    return gsl_spline_eval(heightspline, p2, acc);
}


int main() {
    double h;
    for(int i=0; i<pd.NX; i++)
    {
            x_temp=i*pd.deltaX;
            h=HeightInterpolLiangData(x_temp, 0)
            std::cout << h << std::endl;
    }

return 0;

}
