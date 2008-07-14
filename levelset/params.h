/// \file
/// \brief parameters for two-phase flow problems.

#ifndef DROPS_LSET_PARAMS_H
#define DROPS_LSET_PARAMS_H

#include "misc/params.h"

namespace DROPS
{

/// \brief Parameter class for the problem case
/// 'Droplet in measurement device', Stokes flow
class ParamMesszelleCL: public ParamBaseCL
{
  private:
    void RegisterParams();

  public:
    /// \name Stokes
    //@{
    int    StokesMethod;                        ///< solver for the Stokes problems
    double inner_tol,                           ///< tolerance for Stokes solver
           outer_tol;
    int    inner_iter,                          ///< max. number of iterations for Stokes solver
           outer_iter;
    //@}
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
    double theta, lset_theta;                   ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson
    //@}
    /// \name Coupling
    //@{
    double cpl_tol;                             ///< tolerance for the coupling
    int    cpl_iter;                            ///< max. number of iterations for the fixed-point iteration
    double cpl_stab;                            ///< Laplace-Beltrami-stabilization
    double cpl_proj;                            ///< 1 = perform a projection step before FP
    //@}

    /// \name Material data
    ///
    /// D = drop,    F =  surrounding fluid
    //@{
    double sigma,                               ///< surface tension coefficient
           rhoD,                                ///< density
           rhoF,
           muD,                                 ///< dynamic viscosity
           muF,
           sm_eps;                              ///< width of smooth transition zone for jumping coefficients
    //@}
    ///\name Experimental setup
    //@{
    Point3DCL Radius;                           ///< radii of the ellipsoidal droplet
    Point3DCL Mitte;                            ///< position of the droplet
    double Anstroem,                            ///< max. inflow velocity (parabolic profile)
           r_inlet;                             ///< radius at inlet of measuring device
    int    flow_dir;                            ///< flow direction (x/y/z = 0/1/2)
    Point3DCL g;                                ///< gravity
    double inflow_freq,                         ///< inflow frequence
           inflow_ampl;                         ///< inflow amplitude
    int    bnd_type;                            ///< boundary type: 1= hom. Dirichlet, 2= inhom. Dirichlet bnd-data for in-/outflow, 3= tube/canal
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


    double XFEMStab;                            ///< threshold for discarding ext. dofs parameter, default 0.1

    int    IniCond,                             ///< initial condition (0=Zero, 1/2= stat. flow with/without droplet, -1= read from file)
           GeomType;                            ///< specifies the used geometry (0=ReadMeshBuilder, 1=BrickBuilder)

    string IniData,                             ///< file prefix when reading data for initial condition
           EnsCase,                             ///< name of Ensight Case
           EnsDir,                              ///< local directory for Ensight files
           meshfile,                            ///< mesh file (created by GAMBIT, FLUENT/UNS format) or dimensions of a cuboid (e.g. 2x3x4@5x6x7)
           serialization_file,                  ///< writes multigrid to serialisation files, special value "none" to ignore
           deserialization_file;                ///< reads multigrid from deserialization files, special value "none" to ignore



    ParamMesszelleCL()                        { RegisterParams(); }
    ParamMesszelleCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the problem case
/// 'Droplet in measurement device', Navier-Stokes flow
class ParamMesszelleNsCL: public ParamMesszelleCL
{
  private:
    void RegisterParams();

  public:
    /// \name Navier-Stokes
    //@{
    int    scheme;                              ///< time discretization scheme: 1=lintheta-scheme, 2=rectheta-scheme, 3=theta-scheme, 4=operator splitting, 5=crank-nicolson-scheme
    double nonlinear;                           ///< magnitude of nonlinear term
    double ns_tol,                              ///< Tolerance of the Navier-Stokes-solver
           ns_red;                              ///< The Oseen-residual is reduced by this factor (<1.0)
    int    ns_iter;                             ///< Maximal number of iterations of the solver
    //@}

    ParamMesszelleNsCL()
      : ParamMesszelleCL() { RegisterParams(); }
    ParamMesszelleNsCL( const string& filename)
      : ParamMesszelleCL() { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

class ParamFilmCL: public ParamBaseCL
{ // y = Filmnormal, x = Ablaufrichtung
  private:
    void RegisterParams();

  public:
    double inner_tol, outer_tol,                // Parameter der Loeser
           ns_tol, ns_red, nonlinear,           // fuer Flow & Levelset
           lset_tol, lset_SD, cpl_tol;
    int    inner_iter, outer_iter, 
           ns_iter, lset_iter;
    int    cpl_iter;                            // Kopplung Levelset/Flow: Anzahl Fixpunkt-Iterationen

    double dt;                                  // Zeitschrittweite
    int    num_steps;                           // Anzahl Zeitschritte
    double theta, lset_theta;                   // 0=FwdEuler, 1=BwdEuler, 0.5=CN

    double sigma,                               // Oberflaechenspannung
           CurvDiff,                            // num. Glaettung Kruemmungstermberechnung
           rhoF, rhoG, muF, muG,                // Stoffdaten: Dichte/Viskositaet
           sm_eps,                              // Glaettungszone fuer Dichte-/Viskositaetssprung
           PumpAmpl, PumpFreq;                  // Frequenz und Amplitude der Anregung

    Point3DCL g;                                // Schwerkraft
    double    Filmdicke;                        // Filmdicke
    Point3DCL mesh_res,                         // Gitteraufloesung und
              mesh_size;                        // Gittergroesse in x-/y-/z-Richtung
    int       VolCorr,                          // Volumenkorrektur (0=false)
              IniCond;                          // Anfangsbedingung (0=Null, 1= stat. flow, -1= read from file )

    int    ref_flevel, ref_freq;                // Parameter fuer
    double ref_width;                           // adaptive Verfeinerung

    int    RepFreq, RepMethod;                  // Parameter fuer Reparametrisierung

    string EnsCase,                             // Ensight Case,
           EnsDir,                              // lok.Verzeichnis, in das die geom/vec/scl-files abgelegt werden
           IniData,
           BndCond;

    ParamFilmCL()                        { RegisterParams(); }
    ParamFilmCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

} // end of namespace DROPS

#endif




