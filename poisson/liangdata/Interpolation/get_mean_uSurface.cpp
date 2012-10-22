#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>


int main() {
  //double *x; //[4182]
  //double *y; //[4182]
  //double *z; //[4182]
  double *u; //[4182]
  //double *v; //[4182]
  //double *w; //[4182]
  double *level; //[4182]

  // allocate
  int NumCoords = 4141;
  //x = new double[NumCoords];
  //y = new double[NumCoords];
  //z = new double[NumCoords];
  u = new double[NumCoords];
  //v = new double[NumCoords];
  //w = new double[NumCoords];
  level = new double[NumCoords];


  std::string velocity_filename = "NEWuniformgrid.out";
  

  std::ifstream ufile;

  ufile.open(velocity_filename.c_str(), std::fstream::in);
  
  //Einlesen des Geschwindigkeitsfeldes und des level-sets:
  double curr_x, curr_y, curr_z, curr_u, curr_v, curr_w, curr_level;
  for (int k=0; k<NumCoords; k++) {
    ufile >> curr_x >> curr_y >> curr_z >> curr_u >> curr_v >> curr_w >> curr_level;
    //x[k] = curr_u;
    //y[k] = curr_v;
    //z[k] = curr_w;
    u[k] = curr_u;
    //v[k] = curr_v;
    //w[k] = curr_w;
    level[k] = curr_level;
  }
  
  std::cout << "Level-set-data is:" << std::endl; 
  
  for (int k=1; k<NumCoords; k++) {
    
    
       
       std::cout << level[k] << std::endl; 
       
    
   
  }
  std::cout << "It´s reasonable?" << std::endl; 
  
  double SumUSurface=0.;
  int NU=0;
  double USurf_curr;
    
  for (int k=1; k<NumCoords; k++) {
    
    if(level[k]>0 && level[k-1]<0){
       USurf_curr=(u[k-1] - (u[k]-u[k-1])*level[k-1]/(level[k]-level[k-1]));
       std::cout << "USurf_curr is " << "  " << USurf_curr << std::endl; 
       SumUSurface = SumUSurface + USurf_curr;
       NU=NU+1;
    }
   
  }
  double USurfaceMean;
  USurfaceMean = SumUSurface/NU;
  std::cout << "Mean value for x-component of the velocity-field at fluid-surface is "  << USurfaceMean << std::endl;
  
  return 0;
}