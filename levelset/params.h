/// \file
/// \brief parameters for two-phase flow and other problems.
/// \author Oliver Fortmeier, SC; Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser, IGPM

#ifndef DROPS_LSET_PARAMS_H
#define DROPS_LSET_PARAMS_H

#include "misc/params.h"

namespace DROPS
{

/// \brief Paramater class for ensight parameter
class ParamEnsightCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Ensight
  //@{
    int    ensight;                             ///< Ensight output
    string EnsCase;                             ///< name of Ensight Case
    string EnsDir;                              ///< local directory for Ensight files
    string geomName;                            ///< name for the geometry
    int    masterOut;                           ///< only master writes out ensight files
    int    binary;                              ///< write out ensight files in binary format
  //@}

  public:
    ParamEnsightCL()                        { RegisterParams(); }
    ParamEnsightCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for VTK parameter
class ParamVTKCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name VTK
  //@{
    int    vtk;                                 ///< VTK output
    string vtkDir;                              ///< local directory for vtk files
    string vtkName;                             ///< name of vtk files
    int    vtkBinary;                           ///< write out ensight files in binary format
  //@}

  public:
    ParamVTKCL()                        { RegisterParams(); }
    ParamVTKCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for general information outputs
class ParamInfoCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
    /// \name Information about runs
    //@{
    int printSize;              ///< Should the size of the mg for all MG be displayed
    int printNumUnk;            ///< Print number of unknowns
    int checkMG;                ///< check multigrid for sanety
    //@}
  public:
    ParamInfoCL()                        { RegisterParams(); }
    ParamInfoCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for output on a quadrilateral grid
class ParamQuadCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name quadrilateral output
  //@{
    int quad;                                   ///< Write out quadrilateral grid
    int gridX, gridY, gridZ;                    ///< number of gridpoints
    Point3DCL stepsize;                         ///< Stepsize in each direction
    Point3DCL barycenter;                       ///< Barycenter of the grid
    Point3DCL rotation;                         ///< Rotation
    string quadFileName;                        ///< name of the result-file
  //@}

  public:
    ParamQuadCL()                        { RegisterParams(); }
    ParamQuadCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for (de-)serialization of a (parallel) multigrid
class ParamMGSerCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name (de-)serialization of a multigrid
  //@{
    int    serialization;                 ///< Perform serialization
    int    overwrite;                     ///< Overwrite old output-files or create new for every step.

    string ser_dir,                       ///< writes multigrid to serialisation files, special value "none" to ignore
           deserialization_file;          ///< reads multigrid from deserialization files, special value "none" to ignore
  //@}

  public:
    ParamMGSerCL()                        { RegisterParams(); }
    ParamMGSerCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for creating a brick
class ParamBrickCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Initialization of a brick
  //@{
    int basicref_x;             ///< basic refinements in x-direction
    int basicref_y;             ///< basic refinements in y-direction
    int basicref_z;             ///< basic refinements in z-direction
    Point3DCL orig;             ///< origin of the brick
    Point3DCL dim;              ///< dx, dy, dz
  //@}

  public:
    ParamBrickCL()                        { RegisterParams(); }
    ParamBrickCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for reparametrization
class ParamReparamCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name reparametrize level-set function
  //@{
    int    RepFreq;             ///< number of timesteps before reparametrize the level-set function
    int    RepMethod;           ///< method of reparametrize level-set function
    double NarrowBand;          ///< Narrow-Band method for the Euclidian method (e.g. 0.5 : all dof with <50% of the maximal level-set value are considered)
    double MinGrad,             ///< minimal allowed norm of the gradient of phi
           MaxGrad;             ///< maximal allowed norm of the gradient of phi
  //@}

  public:
    ParamReparamCL()                        { RegisterParams(); }
    ParamReparamCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for Stokes equation
class ParamStokesCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name parameter for the Stokes equation
  //@{
    int    StokesMethod;                        ///< solver for the Stokes problems
    double inner_tol,                           ///< tolerance for Stokes solver
           outer_tol;
    int    inner_iter,                          ///< max. number of iterations for Stokes solver
           outer_iter;
    int    pcA_iter;                            ///< max. number of iterations for the preconditionier
    double pcA_tol,                             ///< tolerance for the preconditioner of A-block
           pcS_tol;                             ///< tolerance for the preconditioner of Schur complement
    double XFEMStab;                            ///< threshold for discarding ext. dofs parameter, default 0.1
  //@}
  public:
    ParamStokesCL()                        { RegisterParams(); }
    ParamStokesCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for adaptive refinement (with level-set function)
class ParamAdaptRefCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name reparametrize level-set function
  //@{
    int    ref_freq;            ///< number of timesteps before adaptive refinement
    int    ref_flevel;          ///< finest level in the near of the phase boundary
    double ref_width;           ///< domain of refinement
  //@}

  public:
    ParamAdaptRefCL()                        { RegisterParams(); }
    ParamAdaptRefCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for Navier-Stokes-Solver
class ParamNavStokesCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name parameter for the Navier-Stokes
  //@{
    double nonlinear;           ///< magnitude of nonlinear term
    double ns_tol,              ///< Tolerance of the Navier-Stokes-solver
           ns_red;              ///< The Oseen-residual is reduced by this factor (<1.0)
    int    ns_iter;             ///< Maximal number of iterations of the solver
    //@}
  public:
    ParamNavStokesCL()                        { RegisterParams(); }
    ParamNavStokesCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for time discretisation
class ParamTimeDiscCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name parameter for the time discretisation
  //@{
    double dt;                                  ///< time step size
    int    num_steps;                           ///< number of timesteps
    double theta;                               ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson
    int    scheme;                              ///< 1=lintheta-scheme, 2=rectheta-scheme, 3=theta-scheme, 4=operator splitting, 5=crank-nicolson-scheme
  //@}
  public:
    ParamTimeDiscCL()                        { RegisterParams(); }
    ParamTimeDiscCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for Level-Set
class ParamLevelSetCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name parameter for the Level-Set function
  //@{
    double lset_tol;                            ///< tolerance for Level Set solver
    int    lset_iter;                           ///< max. number of iterations for Level Set solver
    double lset_SD,                             ///< streamline diffusion parameter
           CurvDiff;                            ///< smoothing of Level Set function before curvature term discretization
    int    VolCorr;                             ///< volume correction (0=no)
  //@}
  public:
    ParamLevelSetCL()                        { RegisterParams(); }
    ParamLevelSetCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for coupling between level set function and Navier-Stokes equation
class ParamCouplingCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name parameter for the coupling between level set function and Navier-Stokes equation
  //@{
    double cpl_tol;                             ///< tolerance for the coupling
    int    cpl_iter;                            ///< max. number of iterations for the fixed-point iteration
    double cpl_stab;                            ///< Laplace-Beltrami-stabilization
    double cpl_proj;                            ///< 1 = perform a projection step before FP
  //@}
  public:
    ParamCouplingCL()                        { RegisterParams(); }
    ParamCouplingCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the material data
class ParamMaterialDataCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Material data
  ///
  /// D = drop,    F =  surrounding fluid
  //@{
    double rhoD,                                ///< density
           rhoF,
           muD,                                 ///< dynamic viscosity
           muF,
           sm_eps;                              ///< width of smooth transition zone for jumping coefficients
  //@}
  public:
    ParamMaterialDataCL()                        { RegisterParams(); }
    ParamMaterialDataCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the experimental data
class ParamExperimentalDataCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
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
  //@}
  public:
    ParamExperimentalDataCL()                        { RegisterParams(); }
    ParamExperimentalDataCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the surface tension
class ParamSurfaceTensionCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  ///\name Surface Force
  //@{
    int st_var;                                 ///< variable surface tension: 0= off, 1=on
    double sigma,                               ///< surface tension coefficient
           st_jumpWidth,                        ///< jump width
           st_relPos,                           ///< position of the jump between upper edge (lambda=0) and barycenter (lambda=1)
           st_red;                              ///< red. factor for surface tension (due to contamination) in lower droplet)
  //@}
  public:
    ParamSurfaceTensionCL()                        { RegisterParams(); }
    ParamSurfaceTensionCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the mass transport
class ParamTransportCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  ///\name Mass Transport
  //@{
    int    transp_do;                           ///< mass transport on (1) or off (0)
    double transp_theta;                        ///< time integration theta scheme
    int    transp_iter;
    double transp_tol,
           transp_diffPos,                      ///< diffusion coefficient (pos. part)
           transp_diffNeg,                      ///< diffusion coefficient (neg. part)
           transp_H,                            ///< Henry number cneg(\f$\infty\f$) = H*cpos(\f$\infty\f$)
           transp_cPos,                         ///< initial concentration (pos. part)
           transp_cNeg;                         ///< initial concentration (neg. part)
  //@}
  public:
    ParamTransportCL()                        { RegisterParams(); }
    ParamTransportCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for domain conditions like boundary cond. and geometry
class ParamDomainCondCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name reparametrize level-set function
  //@{
    int    IniCond,                             ///< initial condition (0=Zero, 1/2= stat. flow with/without droplet, -1= read from file)
           GeomType;                            ///< specifies the used geometry (0=ReadMeshBuilder, 1=BrickBuilder)
    string IniData,                             ///< file prefix when reading data for initial condition
           meshfile;                            ///< mesh file (created by GAMBIT, FLUENT/UNS format) or dimensions of a cuboid (e.g. 2x3x4\@5x6x7)
    int    bnd_type;                            ///< boundary type: 1= hom. Dirichlet, 2= inhom. Dirichlet bnd-data for in-/outflow, 3= tube/canal
  //@}

  public:
    ParamDomainCondCL()                        { RegisterParams(); }
    ParamDomainCondCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the problem case
/// 'Droplet in measurement device', Stokes flow
class ParamMesszelleNsCL:
        public ParamEnsightCL,
        public ParamReparamCL,
        public ParamStokesCL,
        public ParamTimeDiscCL,
        public ParamTransportCL,
        public ParamLevelSetCL,
        public ParamCouplingCL,
        public ParamMaterialDataCL,
        public ParamExperimentalDataCL,
        public ParamSurfaceTensionCL,
        public ParamAdaptRefCL,
        public ParamNavStokesCL,
        public ParamDomainCondCL,
        public ParamMGSerCL

{
  private:
    void RegisterParams();
  public:
    ParamMesszelleNsCL() { RegisterParams(); }
    ParamMesszelleNsCL( const string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamTimeDiscCL::rp_.ReadParams( file);
        ParamTransportCL::rp_.ReadParams( file);
        ParamLevelSetCL::rp_.ReadParams( file);
        ParamCouplingCL::rp_.ReadParams( file);
        ParamEnsightCL::rp_.ReadParams( file);
        ParamMaterialDataCL::rp_.ReadParams( file);
        ParamExperimentalDataCL::rp_.ReadParams( file);
        ParamReparamCL::rp_.ReadParams( file);
        ParamStokesCL::rp_.ReadParams( file);
        ParamNavStokesCL::rp_.ReadParams( file);
        ParamSurfaceTensionCL::rp_.ReadParams( file);
        ParamAdaptRefCL::rp_.ReadParams( file);
        ParamDomainCondCL::rp_.ReadParams( file);
        ParamMGSerCL::rp_.ReadParams( file);
    }
};

class ParamFilmCL: public ParamBaseCL
{ // y = Filmnormal, x = Ablaufrichtung
  private:
    void RegisterParams();

  public:
    int    StokesMethod;                        // solver for the Stokes problems
    double inner_tol, outer_tol,                // Parameter der Loeser
           ns_tol, ns_red, nonlinear,           // fuer Flow & Levelset
           lset_tol, lset_SD, cpl_tol;
    int    inner_iter, outer_iter,
           ns_iter, lset_iter;
    int    pcA_iter;                            // max. number of iterations for the preconditionier
    double pcA_tol,                             // tolerance for the preconditioner
           pcS_tol;
    int    cpl_iter;                            // Kopplung Levelset/Flow: Anzahl Fixpunkt-Iterationen
    double XFEMStab;                            ///< threshold for discarding ext. dofs parameter, default 0.1

    double dt;                                  // Zeitschrittweite
    int    num_steps;                           // Anzahl Zeitschritte
    double theta, lset_theta;                   // 0=FwdEuler, 1=BwdEuler, 0.5=CN

    double sigma,                               // Oberflaechenspannung
           CurvDiff,                            // num. Glaettung Kruemmungstermberechnung
           rhoF, rhoG, muF, muG,                // Stoffdaten: Dichte/Viskositaet
           sm_eps,                              // Glaettungszone fuer Dichte-/Viskositaetssprung
           PumpAmpl, PumpFreq,                  // Frequenz und Amplitude der Anregung
           AmplZ;                               // Amplitude in z-Richtung der initialen Phasengrenze

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
           BndCond,
           serialization_file,                  ///< writes multigrid to serialisation files, special value "none" to ignore
           deserialization_file;                ///< reads multigrid from deserialization files, special value "none" to ignore

    ParamFilmCL()                        { RegisterParams(); }
    ParamFilmCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

} // end of namespace DROPS

#endif




