#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include "poisson/liangdata/Interpolation/PERIODICADDEDphysicaldata.hpp"
#include "poisson/liangdata/Interpolation/HeightInterpolation/createlevelLiang.hpp"
#include "poisson/liangdata/Interpolation/HeightInterpolation/DiscreteLevelSetToDiscreteHeight.hpp"
#include "poisson/liangdata/Interpolation/HeightInterpolation/heightinterpolateLiangData.hpp"
#include "poisson/liangdata/Interpolation/VelInterpolation/functionjaninterpolateLiangData.hpp"
#include "poisson/poissonCoeff.h"

extern DROPS::ParamCL P;

namespace Jan {

    static PhysicalData pd;
    static bool first_call = true;
    static double* u;
    static double* v;
    static double* w;
    static double* level;
    
    void setup()
    {
       //Create an array of Liang's level-set-values, which is converted to an array of discrete height-values by use of the 
       //function "DiscreteLevelSetToDiscreteHeight" (The latter is used in the function "Interface" below):
       createlevelLiang();
       DiscreteLevelSetToDiscreteHeight(createlevelLiang());
       // read data of the velocity-field (and the level-set)
       int NumCoords = pd.NX*pd.NY;
       u = new double[NumCoords];
       v = new double[NumCoords];
       w = new double[NumCoords];
       level = new double[NumCoords];
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
       level[k] = curr_level;
       } 
       first_call = false;
    }
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double Interface( const DROPS::Point3DCL& p, double t){
       
        if (first_call)
            setup();
        
        //Linear interpolation of the previosly (see "setup" above) generated discrete height-profile:
        double retval = HeightInterpolLiangData(p[0],t, DiscreteLevelSetToDiscreteHeight(createlevelLiang()));
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
        return 1.e-5; 
    }

    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.e-9; 
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
    static DROPS::RegisterScalarFunction regscaf("JanALE_Source",       DROPS::instat_scalar_fun_ptr(Source)       );
    //static DROPS::RegisterScalarFunction regscas("JanALE_Inter",         DROPS::instat_scalar_fun_ptr(BInter)     );
    static DROPS::RegisterScalarFunction regscas("JanALE_Solution",     DROPS::instat_scalar_fun_ptr(Solution)     );
    static DROPS::RegisterScalarFunction regscai("JanALE_InitialVal",   DROPS::instat_scalar_fun_ptr(InitialValue) );
}//end of namespace

