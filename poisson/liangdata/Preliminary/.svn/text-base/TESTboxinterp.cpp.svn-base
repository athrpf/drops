#include <iostream>
#include <fstream>
#include <sstream>
#include "functionjaninterpolate.hpp"
#include "boxinterp.hpp"

int main(){

  double *u; //[10980]
  double *v; //[10980]
  double *w; //[10980]
  double *h;

  // allocate
  int NumCoords = 10980;
  u = new double[NumCoords];
  v = new double[NumCoords];
  w = new double[NumCoords];
  h = new double[NumCoords];


  std::string velocity_filename = "velocity.csv";
  std::string coordinate_filename = "coord.csv";
  std::string h_filename = "height.csv";

  std::ifstream ufile;
  std::ifstream cfile;
  std::ifstream hfile;
  ufile.open(velocity_filename.c_str(), std::fstream::in);
  cfile.open(coordinate_filename.c_str(), std::fstream::in);
  hfile.open(h_filename.c_str(), std::fstream::in);

  // Einlesen des Geschwindigkeitsfeldes:
  double curr_u, curr_v, curr_w, x, y, curr_h;
  for (int k=0; k<NumCoords; ++k) {
    ufile >> curr_u >> curr_v >> curr_w;
    cfile >> x >> y;
    hfile >> curr_h;
    //std::cout << "at (x,y)=(" << x << "," << y<< "): " << curr_u << ", " << curr_v << ", " << curr_w <<std::endl;
    u[k] = curr_u;
    v[k] = curr_v;
    w[k] = curr_w;
    h[k] = curr_h;
  }

  // tests
  int NX = 366;
  double deltaX=0.00002;
  double deltaY=0.00002;
  double offsetX=0.00001;
  double offsetY=0.00001;
  double L = NX*deltaX;
  double c = 0.018;
  double T = L/c; // temporal period
  double x;
  int i;
  std::cout << VelInterpol(3.2, 0.0003, 4.7, u) << std::endl;
  std::cout << VelInterpol(3.2, 0.0003, 4.7, v) << std::endl;
  std::cout << VelInterpol(3.2, 0.0003, 4.7, w) << std::endl;
  
  
  
      for(int k=0, k<NX, k++){
      
               x = deltaX + k*deltaX;
    
      
               i = boxnumber(x, deltaX, offsetX);
             
      
          if((k+1)-i < 1.e-16 || (k+1)-i > 1.e-16){
          std::cout << "Error for k=" << k << " :" << "(k+1)-i= " << (k+1)-i << std::endl;
          }
      }
  
  for(int i=0, i<NX, i++){
    std::cout << "should be zero, because the interpolatet velocity-field at a gridpoint should be the same as the given velocity at this point:" << std::endl;
    double udeltavalgridpoint= u[i]-VelInterpol(offsetX + i*deltaX, offsetY, 0, u);
    std::cout << "delta u =" << udeltagridpoint << std::endl;
    double vdeltavalgridpoint= v[i]-VelInterpol(offsetX + i*deltaX, offsetY, 0, v);
    std::cout << "delta v =" << vdeltagridpoint << std::endl;
  
    std::cout << "should be zero because of boundary-condition at the wall:" << std::endl;
    double wallu=VelInterpol(offsetX + i*deltaX, 0, 0, u);
    std::cout << "wallu =" << wallu << std::endl;
    double wallv=VelInterpol(offsetX + i*deltaX, 0, 0, v);
    std::cout << "wallu =" << wallu << std::endl;
  }
  for(int k=0, k<10, k++){
    for(int i=0, i<Numcoords, i++){
      std::cout << "should be zero because of temporal periodicity (period T=L/c)  " << std::endl;
      double udeltaperiodT = VelInterpol(offsetX + i*deltaX, offsetY, 0, u) - VelInterpol(offsetX + i*deltaX, offsetY, k*T, u);
      std::cout << "udeltaperiod" << k <<"T =" << udeltaperiodT << std::endl;
      double vdeltaperiodT = VelInterpol(offsetX + i*deltaX, offsetY, 0, v) - VelInterpol(offsetX + i*deltaX, offsetY, k*T, v);
      std::cout << "vdeltaperiod" << k <<"T=" << vdeltaperiodT << std::endl;
    }
  }
  for(int k=0, k<10, k++){
    for(int i=0, i<Numcoords, i++){
      std::cout << "should be zero because of spacial periodicity (period L)  " << std::endl;
      double udeltaperiodL = VelInterpol(offsetX + i*deltaX, offsetY, 0, u) - VelInterpol(offsetX + i*deltaX + K*L, deltaY + offsetY, 0, u);
      std::cout << "udeltaperiod" << k <<"L =" << udeltaperiodL << std::endl;
      double vdeltaperiodL = VelInterpol(offsetX + i*deltaX, offsetY, 0, v)- VelInterpol(offsetX + i*deltaX + K*L, deltaY + offsetY, 0, v);
      std::cout << "vdeltaperiod" << k <<"L =" << vdeltaperiodL << std::endl;
    }
  }
  return 0;
}
