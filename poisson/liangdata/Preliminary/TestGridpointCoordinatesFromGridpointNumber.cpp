#include <iostream>
#include <math.h>
#include "GridpointCoordinatesFromGridpointNumber.hpp"



int main(){

int NX = 366;
int NY = 30;
double X;
double Y;
double deltaX=0.00002;
double deltaY=0.00002;
double offsetX=0.00001;
double offsetY=0.00001;
double errorX;
double errorY;





for(int j=0; j<NY; j++){
   
    for(int k=0; k<NX; k++){

      X = offsetX +  k*deltaX;     
      Y = offsetY + j*deltaY;

      errorX = gridpointX(j*NX + (k+1), NX, deltaX, offsetX)-X;
      errorY = gridpointY(j*NX + (k+1), NX, deltaY, offsetY)-Y;


       if(errorX > 1.e-15 || errorX < -1.e-15 || errorY > 1.e-15 || errorY < -1.e-15){
          std::cout << "Error for k=   " << k <<"   and j=   " << j << "   errorX=   " << errorX << "   errorY=   " << errorY << "   X=  " << X << "   Y=  " << Y << std::endl;
          
       }
    }
  }    

std::cout << gridpointX(366, NX, deltaX, offsetX) << std::endl;

return 0;

}
