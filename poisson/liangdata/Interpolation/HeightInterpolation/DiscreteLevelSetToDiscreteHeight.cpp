//#include "../PERIODICADDEDphysicaldata.hpp"
//#include <cmath>
//#include <math.h>
//#include <fstream>
//#include <sstream>
//#include "DiscreteLevelSetToDiscreteHeight.hpp"
//static PhysicalData pd;


 
/**Create a discrete one-dimensional height-profile (that is stored in a 
  * .txt - file) from Liang's level-set (given in Liang's format). 
  * The height-profile is cunstructed in such a way, that it has equal values 
  * at both ends (the last array value is taken as equal to the first one). 
  *Additionally the function returns the height-profile as an array h.
  */
double* DiscreteLevelSetToDiscreteHeight(const double* levelLiang) {
  
  double* h = new double[pd.NX];
  double lminus;
  double lplus;
  double yminus;
  double yplus;

  std::ofstream heightfile;
  heightfile.open("heightLiang.txt"); 
  for (int i=0; i<pd.NX-1; i++) {  
       for (int j=0; j<pd.NY-1; j++) {
    
         yminus=j*pd.deltaY; yplus=(j+1)*pd.deltaY;
         lminus=levelLiang[i*pd.NY+j]; lplus=levelLiang[i*pd.NY+j+1];

         if(lminus <=0 && lplus >=0) {
           h[i] = yminus - lminus*((yplus - yminus)/(lplus - lminus)); 
           heightfile << h[i] << std::endl;
         }
      
       } 

   }    

   h[pd.NX-1]=h[0];

   heightfile << h[pd.NX-1] << std::endl;

   heightfile.close();
   
   return h;


}




 
 
 
 
