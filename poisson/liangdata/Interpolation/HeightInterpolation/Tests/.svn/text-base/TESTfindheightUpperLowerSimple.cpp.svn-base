#include "../../PERIODICADDEDphysicaldata.hpp"
#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
static PhysicalData pd;

int main() {

  double *levelLiang; //[4141]

  // allocate
  int NumCoordsLiang = (pd.NX-1)*pd.NY;
  levelLiang = new double[NumCoordsLiang];


  std::string velocity_filename = "LevelSetLiangData.txt";
  std::ifstream ufile;
  ufile.open(velocity_filename.c_str(), std::fstream::in);
  double curr_level;
  for (int k=0; k<NumCoordsLiang; k++) {
    ufile >> curr_level;
    levelLiang[k] = curr_level;
    
  }
 
  
double h[pd.NX];
double lminus;
double lplus;
double yminus;
double yplus;

for (int i=0; i<pd.NX-1; i++){
  
    for (int j=0; j<pd.NY-1; j++){
      yminus=j*pd.deltaY; yplus=(j+1)*pd.deltaY;
      lminus=levelLiang[i*pd.NY+j]; lplus=levelLiang[i*pd.NY+j+1];
      if(lminus <=0 && lplus >=0) {
        h[i] = yminus - lminus*((yplus - yminus)/(lplus - lminus)); 
        
      }
      
    } 


}   


 
// for (int i=0; i<pd.NX-1; i++){

  // std::cout << "h["<<i<<"]="<<h[i]<< std::endl;

// }

double *ylower;
double *yupper;

// allocate
  
  ylower = new double[pd.NX-1];
  yupper = new double[pd.NX-1];


for (int i=0; i<pd.NX-1; i++){
  
    for (int j=0; j<pd.NY-1; j++){
      lminus=levelLiang[i*pd.NY+j]; lplus=levelLiang[i*pd.NY+j+1];
      
      if(lminus <=0 && lplus >=0) {
         ylower[i] = j*pd.deltaY; yupper[i] = (j+1)*pd.deltaY;
      }
      
    } 
 }    

 for (int i=0; i<pd.NX-1; i++){
     if(h[i] < ylower[i] || h[i] > yupper[i]) { std::cout << "  error for i = " << i << " h[i]-ylower[i] =  " << h[i] - ylower[i] << " h[i] -   yupper[i]=  " << h[i] - yupper[i] << std::endl;
     }
  }

return 0;

}




 
 
 
 
