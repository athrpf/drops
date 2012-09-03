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


  std::string velocity_filename = "VelocityWithLevelSetFinal.txt";
  

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
  double T=L/pd.c; //temporal period
  double t=0;
  double deltaU;
  double deltaV;
  double deltaW;
  double deltaLevel;
  double X;
  double Y;
  int NumPeriods=24;
  double errormax=1.e-13;
  
  
    



  for(int k=0; k<NumCoords; k++){
   X = gridpointX(k+1, pd.NX, pd.deltaX, offsetX);
   Y = gridpointY(k+1, pd.NX, pd.deltaY, offsetY);
   //std::cout << k << "(X,Y) is "  << "(" << X << "," << Y << ")" << std::endl;
    for(int j=1; j<NumPeriods; j++){
      deltaU = VelInterpolLiangData(X, Y, t, u)-VelInterpolLiangData(X, Y, t+j*T, u); 
      deltaV = VelInterpolLiangData(X, Y, t, v)-VelInterpolLiangData(X, Y, t+j*T, v); 
      deltaW = VelInterpolLiangData(X, Y, t, w)-VelInterpolLiangData(X, Y, t+j*T, w); 
      deltaLevel = VelInterpolLiangData(X, Y, t, level)-VelInterpolLiangData(X, Y, t+j*T, level);
       if(deltaU > errormax || deltaU < - errormax || deltaV > errormax|| deltaV < - errormax || deltaW > errormax || deltaW < - errormax || deltaLevel > errormax|| deltaLevel < - errormax){
          std::cout << "Error for k=   " << k <<"   and j=   " << j << "   deltaU=   " << deltaU << "   deltaV=   " << deltaV << "   deltaW=   " << deltaW << "   deltaLevel=   " << deltaLevel << "   U=  " <<VelInterpolLiangData(X, Y, t, u) << "   V=  " << VelInterpolLiangData(X, Y, t, v) << "   w=  " << VelInterpolLiangData(X, Y, t, w) << "Level="<< VelInterpolLiangData(X, Y, t, level) << std::endl;
          
       }
    }
  } 
    
 
  
  return 0;
}

  
    



