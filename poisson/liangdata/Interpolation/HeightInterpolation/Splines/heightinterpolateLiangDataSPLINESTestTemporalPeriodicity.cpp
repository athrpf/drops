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
//#include "../createlevelLiang.cpp"
#include "../DiscreteLevelSetToDiscreteHeight.cpp"
#include "../heightinterpolateLiangDataSPLINES.cpp"

//The following function was at first copied from "../createlevelLiang.cpp" and then the 
//path to LevelSetLiangData.txt was updated.
double* createlevelLiang() {

  double* levelLiang; //[4141]

  // allocate
  int NumCoordsLiang = (pd.NX-1)*pd.NY;
  levelLiang = new double[NumCoordsLiang];

  //Please mind: The following relative path to the source-data "LevelSetLiangData.txt" referrs to the 
  //directory where fullDNS.cpp is located.
  std::string level_filename = "../../DataForPoissonCoeff/LevelSetLiangData.txt";
  std::ifstream levelfile;
  levelfile.open(level_filename.c_str(), std::fstream::in);
  double curr_level;
  for (int k=0; k<NumCoordsLiang; k++) {
    levelfile >> curr_level;
    levelLiang[k] = curr_level;
  }
  

  

return levelLiang;


}

int main() {

      double* level; 
      double* level_tmp = createlevelLiang();
      level = DiscreteLevelSetToDiscreteHeight(level_tmp);
      double * xd = new double[pd.NX];
        for(int i=0; i<pd.NX; i++){
            xd[i]=static_cast <double>(i)*pd.deltaX;
        }
        
      //Create spline:
      static gsl_interp_accel * acc;
      static gsl_spline * heightspline;
      acc = gsl_interp_accel_alloc();
      heightspline = gsl_spline_alloc(gsl_interp_cspline, pd.NX);
      gsl_spline_init(heightspline, xd, level, pd.NX);

      double L = pd.deltaX * (pd.NX-1);//Spatial period.
      double T = L/pd.c;//Temporal period.
      double t_0=2.;
      double t;
      double inct=0.01;
      double x=0.5;
      double height;
      double diff;
      double errmax=1.e-15;
      int kmax=100;//#periods
      int maxtimesteps=1000;
        for(int k=0; k<=kmax; k++){
            for(int i=0; i<=maxtimesteps; i++){
                t=t_0 + inct*static_cast <double>(i);
                height = HeightInterpolLiangData(x, t + static_cast <double>(k)*T, acc, heightspline);
                diff = height - HeightInterpolLiangData(x, t, acc, heightspline);
                if(diff<-errmax || diff > errmax){
                   std::cout << height << " period k is "<< k << 
                   " index i is " << i << "   diff = " << diff  << std::endl;
                }
            }
        }

return 0;
}
