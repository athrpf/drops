//The executable of this code creates four double-valued column-.txt-files: heightLiang.txt, xcoarsecoordsSPLINEFinrGrid.txt, heightSPLINEFineGrid.txt and xcoordsSPLINEFineGrid.txt. Plotting heightLiang.txt vs. xcoarsecoordsSPLINEFinrGrid.txt visualizes the linear interpolation of the reference-height-pofile; plotting heightSPLINEFineGrid.txt vs. xcoordsSPLINEFineGrid.txt visualizes the (periodic) interpolation of the reference-height-profile via splines. Comparison of these two plots gives a rough visual estimation of the quality of the spline-interpolation and the linear interpolation respectively.


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
      level = DiscreteLevelSetToDiscreteHeight(level_tmp);//Also creates file heightliang.txt, 
      //that contains pd.NX discrete height-values. These are used below to create a spline.
      double * xd = new double[pd.NX];
      std::ofstream xcoarsecoordsfile;
      xcoarsecoordsfile.open("xcoarsecoordsSPLINEFineGrid.txt");
      for(int i=0; i<pd.NX; i++){
            xd[i]=static_cast <double>(i)*pd.deltaX;
            xcoarsecoordsfile << xd[i] << std::endl;
      }
      xcoarsecoordsfile.close();
      //Create spline:
      static gsl_interp_accel * acc;
      static gsl_spline * heightspline;
      acc = gsl_interp_accel_alloc();
      heightspline = gsl_spline_alloc(gsl_interp_cspline, pd.NX);
      gsl_spline_init(heightspline, xd, level, pd.NX);



      std::ofstream heightfile;
      heightfile.open("heightSPLINEFineGrid.txt");

      std::ofstream xcoordsfile;
      xcoordsfile.open("xcoordsSPLINEFineGrid.txt");


      double t=0.;
      double height;
      int NumPointsFineGrid = 3000; 
      double x;
      double incx = 0.1*pd.deltaX;
        for(int i=0; i<NumPointsFineGrid; i++){
        x= incx*static_cast<double>(i);
        height = HeightInterpolLiangData(x, t, acc, heightspline);
        xcoordsfile << x << std::endl;
        heightfile << height << std::endl;
        }
      
      heightfile.close();
      xcoordsfile.close();

return 0;
}
