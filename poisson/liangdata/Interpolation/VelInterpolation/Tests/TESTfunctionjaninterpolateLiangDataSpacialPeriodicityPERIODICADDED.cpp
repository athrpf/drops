#include <iostream>
#include <fstream>
#include <sstream>
#include "../functionjaninterpolateLiangData.hpp"
#include "GridpointCoordinatesFromGridpointNumber.hpp"
#include "../../PERIODICADDEDphysicaldata.hpp"
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

 
  double offsetX=0.;
  double offsetY=0.;
  double L=(pd.NX-1)*pd.deltaX;
  double t=0.;
  double deltaU;
  double deltaV;
  double deltaW;
  double deltaLevel;
  double XPer;
  double X;
  double Y;
  int NumPeriods=50;
  double errormax = 1.e-13;

  
    


 for(int k=0; k<NumCoords; k++){
   
   Y = gridpointY(k+1, pd.NX, pd.deltaY, offsetY);
   X = gridpointX(k+1, pd.NX, pd.deltaX, offsetX); 
    for(int j=1; j<NumPeriods; j++){
      XPer = X + j*L;
      deltaU = VelInterpolLiangData(X, Y, t, u)-VelInterpolLiangData(XPer, Y, t, u); 
      deltaV = VelInterpolLiangData(X, Y, t, v)-VelInterpolLiangData(XPer, Y, t, v);
      deltaW = VelInterpolLiangData(X, Y, t, w)-VelInterpolLiangData(XPer, Y, t, w);  
      deltaLevel = VelInterpolLiangData(X, Y, t, level)-VelInterpolLiangData(XPer, Y, t, level);
       if(deltaU > errormax || deltaU < -errormax || deltaV > errormax || deltaV < -errormax ||  deltaW > errormax || deltaW < -errormax || deltaLevel > errormax || deltaLevel < -errormax){
          std::cout << "Error for k=   " << k <<"   and j=   " << j << "   deltaU=   " << deltaU << "   deltaV=   " << deltaV << "   deltaW=   " << deltaW << "   deltaLevel=   " << deltaLevel << "   U=  " <<VelInterpolLiangData(X, Y, t, u) << "   V=  " << VelInterpolLiangData(X, Y, t, v) << "   W=  " << VelInterpolLiangData(X, Y, t, w) << "   Level=  " << VelInterpolLiangData(X, Y, t, level) << std::endl;
          
       }
    }
  }    
    
 
  
    
    
  



  return 0;
}
