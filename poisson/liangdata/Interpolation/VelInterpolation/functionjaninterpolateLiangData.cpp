#include <cmath>
#include <math.h>
#include "functionjaninterpolateLiangData.hpp"
#include "boxinterpLiangData.hpp"
#include "physicaldata.hpp"

static PhysicalData pd;

//Interpolate a given component of a velocity-field for suitable spacial- and temporal coordinates (x and t are arbitrary doubles, but y has to be chosen in [0,W] with W as defined below. The gridpoints do not have an offset.)
double VelInterpolLiangData(double x, double y, double t, double* u){

  double L;//lenght of reference-domain
  double W;//height of reference-domain
  W=(pd.NY-1)*pd.deltaY;//height of reference domain
  L=(pd.NX-1)*pd.deltaX;//lenght of reference domain

  double p1=phase(x,t,pd.c);//phase
  double p2=phaseModuloWavelenghtLiangData(p1,L);//Calculate an equivalent phase, that lies in the reference domain.
  
  //Find the relevant boxnumbers for p2 and y ("In which box of the grid is the point (x,p2) located?")
  int nX=boxnumber(p2,pd.deltaX,0);
  int nY=boxnumber(y,pd.deltaY,0);

    
  
  //Once we know the relevant boxnumbers we can calculate the gridpoint-numbers at which we are going to use the given velocity-data (via *u) for interpolation

  if(nX==pd.NX) nX=(nX-1);
  if(nY==pd.NY) nY=(nY-1);

  int N1 = (nY-1)*pd.NX + nX;
  int N2 = N1 + 1;
  int N3 = N1 + pd.NX;
  int N4 = N2 + pd.NX;
  
  //Additionally we need the coordinates of the corners of the box (p2,y) is located in.
  double X1 = (nX-1)*pd.deltaX;
  double Y1 = (nY-1)*pd.deltaY;
  double X2 = nX*pd.deltaX;
  double Y2 = nY*pd.deltaY;


  //The relevant gridpointnumbers N1, N2, N3, N4 needed for bilinear interpolation are enumerated as follows. "1": bottom left corner,  "2": bottom right corner, "3": top left corner, "4": top right corner. The associated coordinates are (in the same order): (X1,Y1), (X2,Y1), (X1,Y2), (X2,Y2). 

    return BilinearInterpol(u[N1-1],u[N2-1],u[N3-1],u[N4-1],X1,Y1,X2,Y2,p2,y);
 
    
}

