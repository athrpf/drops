/// \file
/// \brief parameters for surfactant PDE on the interface.

#ifndef DROPS_SURF_PARAMS_H
#define DROPS_SURF_PARAMS_H

#include "misc/params.h"

namespace DROPS
{

/// \brief Parameter class for the surfactant model
class ParamSurfactantCL: public ParamBaseCL
{
  private:
    void RegisterParams();

  public:
    /// \name Level Set
    //@{
    double lset_tol;                            ///< tolerance for Level Set solver
    int    lset_iter;                           ///< max. number of iterations for Level Set solver
    double lset_SD,                             ///< streamline diffusion parameter
           CurvDiff;                            ///< smoothing of Level Set function before curvature term discretization
    int    VolCorr;                             ///< volume correction (0=no)
    //@}
    /// \name Time discretization
    //@{
    double dt;                                  ///< time step size
    int    num_steps;                           ///< number of timesteps
    double theta_surf, lset_theta;              ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson
    //@}

    /// \name Material data
    ///
    /// D = drop,    F =  surrounding fluid
    //@{
    double muI;                                 ///< diffusion coefficient on the interface
    //@}

    ///\name Experimental setup
    //@{
    Point3DCL Radius;                           ///< radii of the ellipsoidal droplet
    Point3DCL Mitte;                            ///< position of the droplet
    Point3DCL Velocity;                         ///< Velocity field depending on testcase
    //@}

    ///\name Adaptive refinement
    //@{
    int    ref_freq,                            ///< frequency (after how many time steps grid should be adapted)
           ref_flevel;                          ///< finest level of refinement
    double ref_width;                           ///< width of refinement zone around zero-level
    //@}
    ///\name Reparametrization
    //@{
    int    RepFreq,                             ///< frequency (after how many time steps grid should be reparametrized)
           RepMethod;                           ///< 0/1 = fast marching without/with modification of zero level set
    //@}

    int    TestCase,                            ///< 0: Laplace-Beltrami on sphere, v=0
           cdiv;                                ///< #divisions of initial cube

    string EnsCase,                             ///< name of Ensight Case, "none"= no output
           EnsDir;                              ///< local directory for Ensight files



    ParamSurfactantCL ()                       { RegisterParams(); }
    ParamSurfactantCL (const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

} // end of namespace DROPS

#endif




