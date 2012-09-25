#include <iostream>
#include "boxinterp.hpp"

int main(){

  int NX = 366;
  double deltaX=0.00002;
  double offsetX=0.00001;
  double x;
  int i;
  
  std::cout << "(k+1) should be equal to boxnumber. An additional error report will occur, if this is not the case." << std::endl;
  
  for(int k=0; k<NX; k++){
      
               x = offsetX + k*deltaX;   i = boxnumber(x, deltaX, offsetX);
             
          std::cout << "boxnumber for k=" << k << " is :" << i << std::endl;
          
          if((k+1)-i < -1.e-16 || (k+1)-i > 1.e-16){
          std::cout << "Error for k=" << k << " :" << "(k+1)-i= " << (k+1)-i << ", not equal to zero" << std::endl;
          }
  }
  
  return 0;
}
