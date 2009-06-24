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
    double lvs_Tol;                            ///< tolerance for Level Set solver
    int    lvs_Iter;                           ///< max. number of iterations for Level Set solver
    double lvs_Theta;                          ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson
    double lvs_SD,                             ///< streamline diffusion parameter
           lvs_CurvDiff;                       ///< smoothing of Level Set function before curvature term discretization
    int    lvs_VolCorr;                        ///< volume correction (0=no)
    //@}
    /// \name Time discretization
    //@{
    double tm_StepSize;                         ///< time step size
    int    tm_NumSteps;                         ///< number of timesteps
    //@}

    ///\name Experimental setup
    //@{
    Point3DCL exp_Radius;                       ///< radii of the ellipsoidal droplet
    Point3DCL exp_PosDrop;                      ///< position of the droplet
    Point3DCL exp_Velocity;                     ///< Velocity field depending on testcase
    //@}

    ///\name Adaptive refinement
    //@{
    int    ref_Freq,                            ///< frequency (after how many time steps grid should be adapted)
           ref_FinestLevel;                     ///< finest level of refinement
    double ref_Width;                           ///< width of refinement zone around zero-level
    //@}
    ///\name Reparametrization
    //@{
    int    rpm_Freq,                            ///< frequency (after how many time steps grid should be reparametrized)
           rpm_Method;                          ///< 0/1 = fast marching without/with modification of zero level set
    //@}

  ///\name Surfactant Transport
  //@{
    double surf_Theta;     ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson
    int    surf_Iter;      ///< iterations of solver for surfactant equation
    double surf_Tol;       ///< tolerance of solver for surfactant equation
    double surf_OmitBound; ///< omit dof with small mass on the interface
    double surf_Visc;      ///< diffusion coefficient on the interface
  //@}

    int    TestCase,                            ///< 0: Laplace-Beltrami on sphere, v=0
           cdiv;                                ///< \# divisions of initial cube

    string EnsCase,                             ///< name of Ensight Case, "none"= no output
           EnsDir;                              ///< local directory for Ensight files

    ParamSurfactantCL ()                       { RegisterParams(); }
    ParamSurfactantCL (const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

} // end of namespace DROPS

#endif




