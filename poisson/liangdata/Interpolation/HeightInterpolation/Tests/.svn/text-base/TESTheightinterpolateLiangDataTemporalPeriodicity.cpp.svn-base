#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include "../../PERIODICADDEDphysicaldata.hpp"
#include "../heightinterpolateLiangData.hpp"
static PhysicalData pd;

int main() {

  double *h; //[102]
  
 

  // allocate
 
  h = new double[102];
  
  std::string velocity_filename = "heightFIRST.txt";
  

  std::ifstream ufile;

  ufile.open(velocity_filename.c_str(), std::fstream::in);
  
  //Einlesen des Höhenprofils:
  double curr_h;
  for (int k=0; k<102; k++) {
    ufile >> curr_h;
   
    h[k] = curr_h;
     //std::cout << h[k] << std::endl;
    
  }
  
  double L=(pd.NX-1)*pd.deltaX;
  double T=L/pd.c; //temporal period
  double t=0;
  double deltah;
  double X;
  int NumPeriods=24;
  double errormax=1.e-13;
 

  for(int k=0; k<102; k++){
     X = k*pd.deltaX;
     for(int j=1; j<NumPeriods; j++){
      deltah = HeightInterpolLiangData(X, t, h)-HeightInterpolLiangData(X, t+j*T, h); 
     
        if(deltah > errormax || deltah < -errormax){
          std::cout << "Error for k=   " << k << " and j= "<< j << "   deltah=   " << deltah << std::endl;
        
        }
       
      }  
      
  }     
   

  return 0;
}
