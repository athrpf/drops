#include <cmath>
#include <fstream>
#include <sstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "interpolation/boxinterp.hpp"
#include "interpolation/periodicdata.hpp"
#include "interpolation/heightinterpolate.hpp"
#include "poissonCoeff.h"

extern DROPS::ParamCL P;
PeriodicData pd;

namespace PeriodicALE {

  static bool first_call = true;
  static double* u;
  static double* v;
  static double* level;

  static double Dmol;
  static gsl_spline * heightspline;
  static gsl_interp_accel * acc;

  void setup()
  {
    pd.NX = P.get<int>("PeriodicData.nx");
    pd.NY = P.get<int>("PeriodicData.ny");
    pd.deltaX = P.get<double>("PeriodicData.dx");
    pd.deltaY = P.get<double>("PeriodicData.dy");
    pd.c = P.get<double>("PeriodicData.PhaseVelocity");
    std::string height_filename = P.get<std::string>("PeriodicData.HeighFileName"); // ./liangdata/Interpolation/DataForPoissonCoeff/Bala/bala_height.dat";
    std::string velocity_filename = P.get<std::string>("PeriodicData.VelocityFileName"); // "./liangdata/Interpolation/DataForPoissonCoeff/Bala/bala_velocity.dat";

    level = new double[pd.NX];
    std::ifstream heightfile;
    heightfile.open(height_filename.c_str(), std::fstream::in);
    for (int k=0; k<pd.NX; k++)
      heightfile >> level[k];

    // Create spline:
    double * xd = new double[pd.NX];
    for(int i=0; i<pd.NX; i++)
      xd[i] = i*pd.deltaX;
    acc = gsl_interp_accel_alloc();
    heightspline = gsl_spline_alloc(gsl_interp_cspline, pd.NX);
    gsl_spline_init(heightspline, xd, level, pd.NX);

    // read data of the velocity-field:
    int NumCoords = pd.NX*pd.NY;
    u = new double[NumCoords];
    v = new double[NumCoords];
    std::ifstream ufile;
    ufile.open(velocity_filename.c_str(), std::fstream::in);
    for (int k=0; k<NumCoords; k++)
      ufile >> u[k] >> v[k];
    first_call = false;
  }

  double Interface( const DROPS::Point3DCL& p, double t)
  {
    if (first_call)
      setup();

    //Cubic interpolation of the height-profile via the previously created splines of the reference-height-profile:
    double retval = phaseinterp1d(p[0],t, acc, heightspline, pd);
    return retval;
  }

  //Periodic flowfield translates with a constant phase-velocity in x-direction
  DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double t)
  {
    if (first_call)
      setup();
    DROPS::Point3DCL ret;
    double p0 = p[0], p1 = p[1];
    ret[0] = phaseinterp2d(p0, p1, t, u, pd);
    ret[1] = phaseinterp2d(p0, p1, t, v, pd);
    ret[2] = 0.;
    return ret;
  }

  double InitialValue(const DROPS::Point3DCL&, double)
  {
    return P.get<double>("BoundaryConcentration.cInitial");
  }

  double SurfaceConcentration(const DROPS::Point3DCL& p, double){
    double einlauf = P.get<double>("BoundaryConcentration.Einlauf");
    double inlet_concentration = P.get<double>("BoundaryConcentration.cInlet");
    double surface_concentration = P.get<double>("BoundaryConcentration.cSurface");
    if (p[0]<einlauf) {
      return p[0]*(surface_concentration-inlet_concentration)/einlauf;
    }
    return surface_concentration;
  }

  double Diffusion(const DROPS::Point3DCL&, double)
  {
    return P.get<double>("FullDNS.Dmol");
  }

  static DROPS::RegisterScalarFunction regscaa("PeriodicALE_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)    );
  static DROPS::RegisterScalarFunction regscaint("PeriodicALE_Interface",  DROPS::instat_scalar_fun_ptr(Interface)    );
  static DROPS::RegisterVectorFunction regscav("PeriodicALE_Velocity",     DROPS::instat_vector_fun_ptr(Flowfield)    );
  static DROPS::RegisterScalarFunction regscaf("PeriodicALE_Surface",       DROPS::instat_scalar_fun_ptr(SurfaceConcentration)       );
  //static DROPS::RegisterScalarFunction regscas("PeriodicALE_Inter",         DROPS::instat_scalar_fun_ptr(BInter)     );
  static DROPS::RegisterScalarFunction regscai("PeriodicALE_InitialVal",   DROPS::instat_scalar_fun_ptr(InitialValue) );
}//end of namespace
