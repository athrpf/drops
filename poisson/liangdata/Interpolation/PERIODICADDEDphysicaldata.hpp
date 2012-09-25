// Jan Becker
// 2012-07-10

#pragma once

struct PhysicalData {

  double c; //phase-velocity
  double deltaX;//increment x-direction
  double deltaY;//increment y-direction
  int NY;//#gridpoints y-direction
  int NX;//#gridpoints x-direction

  PhysicalData() 
      : c(0.018), deltaX(0.000209), deltaY(0.0001), NY(41), NX(102) 
  {}
  
};
