#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include "../createlevelLiang.hpp"
#include "../DiscreteLevelSetToDiscreteHeight.hpp"
#include "../../PERIODICADDEDphysicaldata.hpp"
#include "../heightinterpolateLiangData.hpp"
static PhysicalData pd;

int main() {

  
  
  double L=(pd.NX-1)*pd.deltaX;
  double T=L/pd.c; //temporal period
  double t=0;
  double deltah;
  double X;
  int NumPeriods=24;
  double errormax=1.e-13;
  
  //Create an array of Liang's level-set-values, which is converted to an array of discrete height-values by use of the 
  //function "DiscreteLevelSetToDiscreteHeight":
  createlevelLiang();
  DiscreteLevelSetToDiscreteHeight(createlevelLiang());
    
 

  for(int k=0; k<102; k++){
     X = k*pd.deltaX;
     for(int j=1; j<NumPeriods; j++){
      deltah = HeightInterpolLiangData(X, t, DiscreteLevelSetToDiscreteHeight(createlevelLiang()))-HeightInterpolLiangData(X, t+j*T, DiscreteLevelSetToDiscreteHeight(createlevelLiang())); 
     
        if(deltah > errormax || deltah < -errormax){
          std::cout << "Error for k=   " << k << " and j= "<< j << "   deltah=   " << deltah << std::endl;
        
        }
       
      }  
      
  }     
   

  return 0;
}
