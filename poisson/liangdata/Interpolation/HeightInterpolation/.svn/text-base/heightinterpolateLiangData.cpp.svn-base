#include <cmath>
#include <math.h>
#include "../boxinterpLiangData.hpp"
#include "../PERIODICADDEDphysicaldata.hpp"
#include "heightinterpolateLiangData.hpp"
static PhysicalData pd;

double HeightInterpolLiangData(double x, double t, double* h) { 

  double L;//lenght of reference-domain
  L=(pd.NX-1)*pd.deltaX;//lenght of reference domain

  double p1=phase(x,t,pd.c);//phase
  double p2=phaseModuloWavelenghtLiangData(p1,L);//Calculate an equivalent phase, that lies in the reference domain.
  
  //Find the relevant boxnumber for p2
  int nX=boxnumber(p2,pd.deltaX,0);
  


  if(nX==pd.NX) nX=(nX-1);
  
  double X1=(nX-1)*pd.deltaX;
  double X2=nX*pd.deltaX;


    return(((h[nX]-h[nX-1])/pd.deltaX)*(p2-X1) + h[nX-1]);
}

