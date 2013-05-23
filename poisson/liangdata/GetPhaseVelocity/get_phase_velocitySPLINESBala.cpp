#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include "../Interpolation/PERIODICADDEDphysicaldataBala.hpp"
static PhysicalDataBala pd;
#include "../Interpolation/boxinterpLiangData.cpp"
#include "../Interpolation/HeightInterpolation/heightinterpolateLiangDataSPLINES.cpp"
#include "../Interpolation/VelInterpolation/functionjaninterpolateLiangData.cpp"

    static bool first_call = true;
    static double* u;
    static double* v;
    static double* h;
    static gsl_interp_accel * acc;
    static gsl_spline * heightspline;
    
    void setup(void);
    double RefInterface( double x);
    double RefVelU( double x, double y);
    double RefVelV( double x, double y);
    double diff_x(double a_x, double a_x_plus_dx, double dx);
    
    void setup()
    {
      h = new double[pd.NX];
      std::string height_filename = "../Interpolation/DataForPoissonCoeff/Bala/bala_height.dat";
      std::ifstream heightfile;
      heightfile.open(height_filename.c_str(), std::fstream::in);
      double curr_h;
      for (int k=0; k<pd.NX; k++) {
         heightfile >> curr_h;
         h[k] = curr_h;
      }
      std::vector<double> xd(pd.NX);
      for(int i=0; i<pd.NX; i++) {
         xd[i]=i*pd.deltaX;
      }
      //Create spline:
      acc = gsl_interp_accel_alloc();
      heightspline = gsl_spline_alloc(gsl_interp_cspline, pd.NX);
      gsl_spline_init(heightspline, &xd[0], h, pd.NX);
      
      // read data of the velocity-field:
      int NumCoords = pd.NX*pd.NY;
      u = new double[NumCoords];
      v = new double[NumCoords];
      std::string velocity_filename = "../Interpolation/DataForPoissonCoeff/Bala/bala_velocity.dat";
      std::ifstream ufile;
      ufile.open(velocity_filename.c_str(), std::fstream::in);
      double curr_u, curr_v;
      for (int k=0; k<NumCoords; k++) {
         ufile >> curr_u >> curr_v;
         u[k] = curr_u;
         v[k] = curr_v;
      }
      first_call = false;
    }

    double RefInterface( double x){
      if (first_call)
        setup();

      double retval = HeightInterpolLiangData(x,0, acc, heightspline);
      return retval;
    }

    double RefVelU( double x, double y){
        if (first_call)
            setup();
     
        double retval = VelInterpolLiangData(x,y,0, u);
        return retval;
    }
    
    double RefVelV( double x, double y){
        if (first_call)
            setup();
     
        double retval = VelInterpolLiangData(x,y,0, v);
        return retval;
    }
    
    double diff_x(double a_x, double a_x_plus_dx, double dx){
    
        double retval = (a_x_plus_dx-a_x)/dx;
        return retval; 
    
    }
    
int main() {

    int refpoints=pd.NX;
    //double L=(pd.NX-1)*pd.deltaX;
    double inc_x=pd.deltaX;
    //double inc_x=L/(refpoints-1);
    double x_0;
    double h;
    double h_plus;
    double dx=1.e-13;
    double u;
    double v;
    double v_w;
    double h_x;
    double epsilon=10.e-2;
    double c;
    double fehlerfaktor;
    double relfehlerfaktor;
    
    
    std::ofstream cfile;
    cfile.open("phasevelocityBalaSLPINES.txt"); 
  
    for(int i=0; i<refpoints-1; i++){
    
        x_0 = (i+0.5)*inc_x;
        h = RefInterface(x_0);
        h_plus = RefInterface(x_0+dx);
        h_x = diff_x(h, h_plus, dx);
        u = RefVelU(x_0, h); 
        v = RefVelV(x_0, h);
        v_w = RefVelV(x_0, 0.);
        if (h_x<epsilon && h_x>-epsilon){
            c = u + (v_w - v)/h_x; 
            //fehlerfaktor=sqrt(1. + 1./(h_x*h_x) + (v*v)/(h_x*h_x*h_x*h_x));
            //relfehlerfaktor=fehlerfaktor/c;
            std::cout << "Absolute-value of the derivative is quite small! c = " << c <<//" Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << 
            std::endl;
            cfile << "Absolute-value of the derivative is quite small! c = " << c <<//" Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << 
            std::endl;
        }
        else { 
            c = u + (v_w - v)/h_x;
            fehlerfaktor=sqrt(1. + 1./(h_x*h_x) + (v*v)/(h_x*h_x*h_x*h_x));
            relfehlerfaktor=fehlerfaktor/c;
            std::cout << "Calculated value for phase-velocity is " << c <<// " Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << 
            std::endl;
            cfile << "Calculated value for phase-velocity is " << c <<// " Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << 
            std::endl;
            
        }
    }
    
   cfile.close();
 
return 0;
}   





