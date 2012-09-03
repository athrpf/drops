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
 
  
double h[pd.NX-1];
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
double LevelSetMatrix[41][101];

// allocate
  
  ylower = new double[101];
  yupper = new double[101];
  
  
//Construction LevelSetMatrix. LevelSetMatrix[i][j] denotes the LevelSet at x=j*deltaX and y=i*deltaY.   
for (int i=0; i<41; i++){
  
    for (int j=0; j<101; j++){
      LevelSetMatrix[i][j]=levelLiang[j*41+i]; 
      
    } 
 }      


for (int j=0; j<101; j++){
  
    for (int i=0; i<40; i++){
      lminus=LevelSetMatrix[i][j]; lplus=LevelSetMatrix[i+1][j];
      
      if(lminus <=0 && lplus >=0) {
         ylower[j] = i*pd.deltaY; yupper[j] = (i+1)*pd.deltaY;
      }
      
    } 
 }    

 for (int j=0; j<101; j++){
     if(h[j] < ylower[j] || h[j] > yupper[j]) { std::cout << "  error for j = " << j << " h[j]-ylower[j] =  " << h[j] - ylower[j] << " h[j] -   yupper[j]=  " << h[j] - yupper[j] << std::endl;
     }
  }

return 0;

}




 
 
 
 
