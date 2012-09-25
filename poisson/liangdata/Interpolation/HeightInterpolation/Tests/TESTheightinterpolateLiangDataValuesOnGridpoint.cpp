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

  
  double offsetX = 0.;
  double t=0.;
  double deltah;
  double X;
  double errormax = 1.e-16;

  for(int k=0; k<102; k++){
     X = k*pd.deltaX;
     
    
     deltah = HeightInterpolLiangData(X, t, h)-h[k]; 
     
       if(deltah > errormax || deltah < -errormax){
       std::cout << "Error for k=   " << k << "   deltah=   " << deltah <<  h[k] << "   h[k]=  " << h[k] << std::endl;
        
       }
  }    

  return 0;
}
