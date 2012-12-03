#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "./liangdata/Interpolation/PERIODICADDEDphysicaldata.hpp"
static PhysicalData pd;

#include "./liangdata/Interpolation/boxinterpLiangData.cpp"

// #include "./liangdata/Interpolation/HeightInterpolation/createlevelLiang.hpp"
#include "./liangdata/Interpolation/HeightInterpolation/createlevelLiang.cpp"

// #include "./liangdata/Interpolation/HeightInterpolation/DiscreteLevelSetToDiscreteHeight.hpp"
#include "./liangdata/Interpolation/HeightInterpolation/DiscreteLevelSetToDiscreteHeight.cpp"

// #include "./liangdata/Interpolation/HeightInterpolation/heightinterpolateLiangData.hpp"
//#include "./liangdata/Interpolation/HeightInterpolation/heightinterpolateLiangData.cpp"

// #include "./liangdata/Interpolation/VelInterpolation/functionjaninterpolateLiangData.hpp"
#include "./liangdata/Interpolation/VelInterpolation/functionjaninterpolateLiangData.cpp"

#include "poissonCoeff.h"

extern DROPS::ParamCL P;

namespace Jan {

    static PhysicalData pd;
    static bool first_call = true;
    static double* u;
    static double* v;
    static double* w;
    static double* level;

  static double Dmol;

    void setup()
    {
      //Create an array of Liang's level-set-values, which is converted to
      //an array of discrete height-values by use of the
      //function "DiscreteLevelSetToDiscreteHeight":
      double* level_tmp = createlevelLiang();
      level = DiscreteLevelSetToDiscreteHeight(level_tmp);
      //Create array of discrete reference-points on x-axis that are, in addition to the array "level" needed to 
      //create a spline for the reference-heigt-profile: 
      double * xd = new double[pd.X];
      for(int i=0; i<pd.NX; i++)
      {
            xd[i]=i*pd.deltaX;
      }
      //Create spline:
      static gsl_spline ** heightspline;
      static gsl_interp_accel * acc;
      acc = gsl interp_accel_alloc();
      heightspline = new gsl_spline;
      heightspline = gsl_splie_alloc(gsl_interp_cspline, pd.NX);
      gsl_spline_init(heightspline, xd, level, pd.NX);
      
      // read data of the velocity-field (and the level-set)
      int NumCoords = pd.NX*pd.NY;
      u = new double[NumCoords];
      v = new double[NumCoords];
      w = new double[NumCoords];
      std::string velocity_filename = "./liangdata/Interpolation/DataForPoissonCoeff/VelocityWithLevelSetFinal.txt";
      std::ifstream ufile;
      ufile.open(velocity_filename.c_str(), std::fstream::in);
      double curr_u, curr_v, curr_w, curr_level;
      for (int k=0; k<NumCoords; k++)
      {
         ufile >> curr_u >> curr_v >> curr_w >> curr_level;
         u[k] = curr_u;
         v[k] = curr_v;
         w[k] = curr_w;
      }
       first_call = false;
    }
#include "./liangdata/Interpolation/HeightInterpolation/heightinterpolateLiangDataSPLINES.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double Interface( const DROPS::Point3DCL& p, double t){
      if (first_call)
        setup();

        //Cubic interpolation of the height-profile via the previously created splines of the reference-height-profile:
        double retval = HeightInterpolLiangData(p[0],t);
        return retval;
    }

    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t){

        DROPS::Point3DCL ref(0.);
        ref[0] = p[0];
        ref[1] = 0.2 * p[1]/Interface(p, t);
        ref[2] = p[2];
        return ref;
    }

    //Periodic flowfield translates with a constant phase-velocity in x-direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double t){

        if (first_call)
            setup();
        DROPS::Point3DCL ret;
        ret[0]=VelInterpolLiangData(p[0],p[1],t, u);
        ret[1]=VelInterpolLiangData(p[0],p[1],t, v);
        ret[2]=VelInterpolLiangData(p[0],p[1],t, w);
        return ret;
    }


    double Source(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    
    double InitialValue(const DROPS::Point3DCL&, double){
        return P.get<double>("JanOptions.cInitial");
    }

    double SurfaceConcentration(const DROPS::Point3DCL& p, double){
        double einlauf = P.get<double>("JanOptions.Einlauf");
        double inlet_concentration = P.get<double>("JanOptions.cInlet");
        double surface_concentration = P.get<double>("JanOptions.cSurface");
        if (p[0]<einlauf) {
            return p[0]*(surface_concentration-inlet_concentration)/einlauf;
        }
        return surface_concentration;
    }

    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }

    double Diffusion(const DROPS::Point3DCL&, double){
      return P.get<double>("FullDNS.Dmol");
    }

    /// \brief Solution
    double Solution( const DROPS::Point3DCL&, double)
    {
        return 1e-4;
    }

    //double BInter(const DROPS::Point3DCL&, double){
      //  return 0.01;
    //}

    static DROPS::RegisterScalarFunction regscaq("JanALE_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)     );
    static DROPS::RegisterScalarFunction regscaa("JanALE_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)    );
    static DROPS::RegisterScalarFunction regscaint("JanALE_Interface",  DROPS::instat_scalar_fun_ptr(Interface)    );
    static DROPS::RegisterVectorFunction regscav("JanALE_Velocity",     DROPS::instat_vector_fun_ptr(Flowfield)    );
    static DROPS::RegisterScalarFunction regscaf("JanALE_Surface",       DROPS::instat_scalar_fun_ptr(SurfaceConcentration)       );
    //static DROPS::RegisterScalarFunction regscas("JanALE_Inter",         DROPS::instat_scalar_fun_ptr(BInter)     );
    static DROPS::RegisterScalarFunction regscas("JanALE_Solution",     DROPS::instat_scalar_fun_ptr(Solution)     );
    static DROPS::RegisterScalarFunction regscai("JanALE_InitialVal",   DROPS::instat_scalar_fun_ptr(InitialValue) );
}//end of namespace

