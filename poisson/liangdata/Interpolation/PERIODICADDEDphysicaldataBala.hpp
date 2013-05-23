// Jan Becker
// 2013-05-22

#pragma once

struct PhysicalDataBala {

  double c; //phase-velocity
  double deltaX;//increment x-direction
  double deltaY;//increment y-direction
  int NY;//#gridpoints y-direction
  int NX;//#gridpoints x-direction

  PhysicalDataBala() 
      : c(0.385), deltaX(0.000110849885544), deltaY(2.10314439544e-06), NY(100), NX(301) 
  {}
  
};
