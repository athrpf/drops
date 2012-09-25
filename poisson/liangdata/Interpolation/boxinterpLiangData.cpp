// Jan Becker
// 2012-06-19

//#include <math.h>
//#include "boxinterpLiangData.hpp"

//Enumerate the squares (needed for the bilinear interpolation) in x- resp. in y-direction
int boxnumber(double x, double delta, double offset)
{
  double q=(x-offset)/delta + 1.e-14;
  int Rq=(int)(floor(q));
  int m=Rq+1;
  return m;
}

  //Calculate a phase-value within [0,L], that is equivalent to the original one.
    double phaseModuloWavelenghtLiangData(double phi, double L)
  {
     
    double pL=phi - floor(phi/L + 1.e-14)*L;
   if(phi > - 1.e-15 && phi <  1.e-15) return 0; 
   else if(pL > - 1.e-15 && pL  < 1.e-15 ) return L;
   else return pL;
  }

//Calculate space- ant time-depentent phase (c designates the phase-velocity)
double phase(double x, double t, double c)
{

  double p=x-c*t;
  return p;
}

// definition interpolation-function
double BilinearInterpol(double UQ11, double UQ21, double UQ12, double UQ22, double x1, double y1, double x2, double y2, double X, double Y)
{

  double UR1=((x2-X)/(x2-x1))*UQ11 + ((X-x1)/(x2-x1))*UQ21;
  double UR2=((x2-X)/(x2-x1))*UQ12 + ((X-x1)/(x2-x1))*UQ22;
  double UXY=((y2-Y)/(y2-y1))*UR1 + ((Y-y1)/(y2-y1))*UR2;

  return UXY;
}


