#include "../../PERIODICADDEDphysicaldata.hpp"
#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include "../heightinterpolateLiangData.hpp"
static PhysicalData pd;

int main() {

  double *coarseheight; //[102]
  
  // allocate
  coarseheight = new double[pd.NX];


  std::string velocity_filename = "heightFIRST.txt";
  std::ifstream ufile;
  ufile.open(velocity_filename.c_str(), std::fstream::in);
  double curr_height;
  for (int k=0; k<pd.NX; k++) {
    ufile >> curr_height;
    coarseheight[k] = curr_height;
  }
  
 std::ofstream heightfile;
 heightfile.open("heightFineGrid.txt"); 
 
 std::ofstream xfile;
 xfile.open("xFineGrid.txt"); 
 
int M=10000; // #intervals of the finer grid
double L=(pd.NX-1)*pd.deltaX;
double newdeltaX=L/M; // increment in x-direction for the finer grid
double t=0.;  
double curr_x;
double curr_h;

for (int i=0; i<M+1; i++){
  
 curr_x=i*newdeltaX;
 curr_h = HeightInterpolLiangData(curr_x, t, coarseheight);   
 heightfile << curr_h << std::endl; 
 xfile << curr_x << std::endl;   
} 

   

heightfile.close();
xfile.close();

return 0;

}




 
 
 
 
