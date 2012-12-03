#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../Interpolation/PERIODICADDEDphysicaldata.hpp"
static PhysicalData pd;
#include "../Interpolation/boxinterpLiangData.cpp"
double* createlevelLiang() {
  double* levelLiang; //[4141]
  int NumCoordsLiang = (pd.NX-1)*pd.NY;
  levelLiang = new double[NumCoordsLiang];
  std::string level_filename = "../Interpolation/DataForPoissonCoeff/LevelSetLiangData.txt";
  std::ifstream levelfile;
  levelfile.open(level_filename.c_str(), std::fstream::in);
  double curr_level;
  for (int k=0; k<NumCoordsLiang; k++) {
    levelfile >> curr_level;
    levelLiang[k] = curr_level;
  }
return levelLiang;
}
//#include "../Interpolation/HeightInterpolation/createlevelLiang.cpp"
#include "../Interpolation/HeightInterpolation/DiscreteLevelSetToDiscreteHeight.cpp"
#include "../Interpolation/HeightInterpolation/heightinterpolateLiangData.cpp"
#include "../Interpolation/VelInterpolation/functionjaninterpolateLiangData.cpp"




    
    static bool first_call = true;
    static double* u;
    static double* v;
    static double* w;
    static double* level;

    void setup(void);
    double RefInterface( double x);
    double RefVelU( double x, double y);
    double RefVelV( double x, double y);
    double diff_x(double a_x, double a_x_plus_dx, double dx);

    void setup()
    {
       //Create an array of Liang's level-set-values, which is converted to
       //an array of discrete height-values by use of the
       //function "DiscreteLevelSetToDiscreteHeight" (The latter is used in
       //the function "RefInterface" below):
      double* level_tmp = createlevelLiang();
      level = DiscreteLevelSetToDiscreteHeight(level_tmp);
       // read data of the velocity-field (and the level-set)
       int NumCoords = pd.NX*pd.NY;
       u = new double[NumCoords];
       v = new double[NumCoords];
       w = new double[NumCoords];
       std::string velocity_filename = "../Interpolation/DataForPoissonCoeff/VelocityWithLevelSetFinal.txt";
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

    double RefInterface( double x){
      if (first_call)
        setup();

        double retval = HeightInterpolLiangData(x,0, level);
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
    double h_x;
    double epsilon=10.e-2;
    double c;
    double fehlerfaktor;
    double relfehlerfaktor;
    
    
    std::ofstream cfile;
    cfile.open("phasevelocity.txt"); 
  
    for(int i=0; i<refpoints-1; i++){
    
        x_0 = (i+0.5)*inc_x;
        h = RefInterface(x_0);
        h_plus = RefInterface(x_0+dx);
        h_x = diff_x(h, h_plus, dx);
        u = RefVelU(x_0, h); 
        v = RefVelV(x_0, h);
        if (h_x<epsilon && h_x>-epsilon){
            c = u - v/h_x; 
            fehlerfaktor=sqrt(1. + 1./(h_x*h_x) + (v*v)/(h_x*h_x*h_x*h_x));
            relfehlerfaktor=fehlerfaktor/c;
            std::cout << "Absolute-value of the derivative is quite small! c = " << c <<" Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << std::endl;
            cfile << "Absolute-value of the derivative is quite small! c = " << c <<" Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << std::endl;
        }
        else { 
            c = u - v/h_x;
            fehlerfaktor=sqrt(1. + 1./(h_x*h_x) + (v*v)/(h_x*h_x*h_x*h_x));
            relfehlerfaktor=fehlerfaktor/c;
            std::cout << "Calculated value for phase-velocity is " << c << " Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << std::endl;
            cfile << "Calculated value for phase-velocity is " << c << " Fehlerfaktor = "<< fehlerfaktor << " RelFehlerfaktor = "<< relfehlerfaktor << std::endl;
            
        }
    }
    
   cfile.close();
 
return 0;
}   





