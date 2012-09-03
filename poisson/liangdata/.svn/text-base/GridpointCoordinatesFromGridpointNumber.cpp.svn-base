#include <math.h>
#include "GridpointCoordinatesFromGridpointNumber.hpp"



//Calculates the x-coordinate of a gridpoint for a given gridpoint-# 
double gridpointX(int gridpointnumber, int NX, double deltaX, double offsetX){

   double Xcoordinate;
   
   
   while(gridpointnumber-NX>0){
          gridpointnumber=gridpointnumber-NX;
      
   }

   if(gridpointnumber-NX==0){
      Xcoordinate=offsetX + (NX-1)*deltaX;
   }
   else{
      Xcoordinate= offsetX + (gridpointnumber-1)*deltaX; 
   }

   
   return Xcoordinate;
   
   }



//Calculates the y-coordinate of a gridpoint for a given gridpoint-# 
double gridpointY(int gridpointnumber, int NX, double deltaY, double offsetY){
   
   double Ycoordinate;
   int k = 0;
   
   while(gridpointnumber-NX>0){
          gridpointnumber=gridpointnumber-NX;
          k=k+1;
   }

   Ycoordinate = offsetY + k*deltaY;

   return Ycoordinate;
   
   }




