//As a test the interpolation-function VelInterpolLiangData shall be evaluated on a grid, that is much finer than that one used for bilinear interpolation of Liang's data itself.
#include <cmath>
#include <math.h>
#include "../functionjaninterpolateLiangData.hpp"
#include "../../PERIODICADDEDphysicaldata.hpp"
#include <fstream>
#include <sstream>
static PhysicalData pd;

int main() {

  double *u; //[4182]
  double *v; //[4182]
  double *w; //[4182]
  double *level; //[4182]

  // allocate
  int NumCoords = pd.NX*pd.NY;
  u = new double[NumCoords];
  v = new double[NumCoords];
  w = new double[NumCoords];
  level = new double[NumCoords];


  std::string velocity_filename = "../../DataForPoissonCoeff/VelocityWithLevelSetFinal.txt";
  

  std::ifstream ufile;

  ufile.open(velocity_filename.c_str(), std::fstream::in);
  
  //Einlesen des Geschwindigkeitsfeldes:
  double curr_u, curr_v, curr_w, x, y, curr_level;
  for (int k=0; k<NumCoords; k++) {
    ufile >> curr_u >> curr_v >> curr_w >> curr_level;
   
    u[k] = curr_u;
    v[k] = curr_v;
    w[k] = curr_w;
    level[k] = curr_level;
  }

//Evaluate the interpolation-function on the fine grid and export the computed velocity-data and the coordinates of the fine grid via the files "vel.txt" resp. "coords.txt".
int MX=1000; // #interpolation-points in x-direction
int MY=100; // #interpolation-points in y-direction
double testdeltaX=(pd.NX-1)*pd.deltaX/(MX-1); //increment in x-direction for fine-grid;
double testdeltaY=(pd.NY-1)*pd.deltaY/(MY-1);//increment in y-direction for fine-grid;
double X;
double Y;
double t=0;
double ucurrVelInterpol;
double vcurrVelInterpol;
double wcurrVelInterpol;
double levelcurrVelInterpol;


 
 std::ofstream xcoordsfile;
 xcoordsfile.open("xcoordsPeriodicADDED.txt");

 std::ofstream ycoordsfile;
 ycoordsfile.open("ycoordsPeriodicADDED.txt");
 
 std::ofstream uvelfile;
 uvelfile.open("uvelLiangPeriodicADDED.txt");

 std::ofstream vvelfile;
 vvelfile.open("vvelLiangPeriodicADDED.txt");
 
 std::ofstream wvelfile;
 wvelfile.open("wvelLiangPeriodicADDED.txt");
 
 std::ofstream levelfile;
 levelfile.open("levelLiangPeriodicADDED.txt");

      for(int i=0; i<MY; i++) {

           

           for(int k=0; k<MX; k++) {

               Y = i*testdeltaY;
               X = k*testdeltaX;
                
               ucurrVelInterpol = VelInterpolLiangData(X, Y, t, u); 
               vcurrVelInterpol = VelInterpolLiangData(X, Y, t, v);
               wcurrVelInterpol = VelInterpolLiangData(X, Y, t, w);
               levelcurrVelInterpol = VelInterpolLiangData(X, Y, t, level);
               
               uvelfile << ucurrVelInterpol << std::endl; 
               vvelfile << vcurrVelInterpol << std::endl;
               wvelfile << wcurrVelInterpol << std::endl;
               levelfile << levelcurrVelInterpol << std::endl;

               xcoordsfile << X << std::endl;
               ycoordsfile << Y << std::endl;

           }
           

       }

 uvelfile.close();

 vvelfile.close();
 
 wvelfile.close();
 
 xcoordsfile.close();

 ycoordsfile.close();
 
 levelfile.close();
    
  
 return 0;

}
