/// \file params.h
/// \brief parameters for two-phase flow and other problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

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
    int    ens_EnsightOut;                          ///< Ensight output
    string ens_EnsCase;                             ///< name of Ensight Case
    string ens_EnsDir;                              ///< local directory for Ensight files
    string ens_GeomName;                            ///< name for the geometry
    int    ens_MasterOut;                           ///< only master writes out ensight files
    int    ens_Binary;                              ///< write out ensight files in binary format
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
    int    vtk_VTKOut;                           ///< VTK output
    string vtk_VTKDir;                           ///< local directory for vtk files
    string vtk_VTKName;                          ///< name of vtk files
    int    vtk_Binary;                           ///< write out ensight files in binary format
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
    int inf_PrintSize;              ///< Should the size of the mg for all MG be displayed
    int inf_PrintNumUnk;            ///< Print number of unknowns
    int inf_CheckMG;                ///< check multigrid for sanety
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
    int qlg_Quad;                                   ///< Write out quadrilateral grid
    int qlg_GridX, qlg_GridY, qlg_GridZ;            ///< number of gridpoints
    Point3DCL qlg_Stepsize;                         ///< Stepsize in each direction
    Point3DCL qlg_Barycenter;                       ///< Barycenter of the grid
    Point3DCL qlg_Rotation;                         ///< Rotation
    string qlg_FileName;                            ///< name of the result-file
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
    int    rst_Serialization;                 ///< Perform serialization
    int    rst_Overwrite;                     ///< Overwrite old output-files or create new for every step.

    string rst_Outputfile,                    ///< writes multigrid to serialisation files, special value "none" to ignore
           rst_Inputfile;                     ///< reads multigrid from deserialization files, special value "none" to ignore
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
    int brk_BasicRefX;             ///< basic refinements in x-direction
    int brk_BasicRefY;             ///< basic refinements in y-direction
    int brk_BasicRefZ;             ///< basic refinements in z-direction
    Point3DCL brk_orig;            ///< origin of the brick
    Point3DCL brk_dim;             ///< dx, dy, dz
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
    int    rpm_Freq;            ///< number of timesteps before reparametrize the level-set function
    int    rpm_Method;          ///< method of reparametrize level-set function
    double rpm_NarrowBand;      ///< Narrow-Band method for the Euclidian method (e.g. 0.5 : all dof with <50% of the maximal level-set value are considered)
    double rpm_MinGrad,         ///< minimal allowed norm of the gradient of phi
           rpm_MaxGrad;         ///< maximal allowed norm of the gradient of phi
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
    int    stk_StokesMethod;                      ///< solver for the Stokes problems
    double stk_InnerTol,                          ///< tolerance for Stokes solver
           stk_OuterTol;
    int    stk_InnerIter,                         ///< max. number of iterations for Stokes solver
           stk_OuterIter;
    int    stk_PcAIter;                           ///< max. number of iterations for the preconditionier
    double stk_PcATol,                            ///< tolerance for the preconditioner of A-block
           stk_PcSTol;                            ///< tolerance for the preconditioner of Schur complement
    double stk_XFEMStab;                          ///< threshold for discarding ext. dofs parameter, default 0.1
    double stk_Theta;                             ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson

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
    int    ref_Freq;               ///< number of timesteps before adaptive refinement
    int    ref_FinestLevel;        ///< finest level in the near of the phase boundary
    int    ref_CoarsestLevel;      ///< coarsest level in the near of the phase boundary
    double ref_Width;              ///< domain of refinement
    int    ref_RefineStrategy;     ///< algorithm to determine the load balancing graph for refinement
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
    double ns_Nonlinear;         ///< magnitude of nonlinear term
    double ns_Tol,               ///< Tolerance of the Navier-Stokes-solver
           ns_Reduction;         ///< The Oseen-residual is reduced by this factor (<1.0)
    int    ns_Iter;              ///< Maximal number of iterations of the solver
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
    double tm_StepSize;                            ///< time step size
    int    tm_NumSteps;                            ///< number of timesteps
    int    tm_Scheme;                              ///< 1=lintheta-scheme, 2=rectheta-scheme, 3=theta-scheme, 4=operator splitting, 5=crank-nicolson-scheme
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
    double lvs_Tol;                           ///< tolerance for Level Set solver
    int    lvs_Iter;                          ///< max. number of iterations for Level Set solver
    double lvs_SD,                            ///< streamline diffusion parameter
           lvs_CurvDiff;                      ///< smoothing of Level Set function before curvature term discretization
    int    lvs_VolCorrection;                 ///< volume correction (0=no)
    double lvs_Theta;                         ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson

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
    double cpl_Tol;                             ///< tolerance for the coupling
    int    cpl_Iter;                            ///< max. number of iterations for the fixed-point iteration
    double cpl_Stab;                            ///< Laplace-Beltrami-stabilization
    double cpl_Projection;                      ///< 1 = perform a projection step before FP
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
    double mat_DensDrop,                        ///< density
           mat_DensFluid,
           mat_ViscDrop,                        ///< dynamic viscosity
           mat_ViscFluid,
           mat_SmoothZone;                      ///< width of smooth transition zone for jumping coefficients
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
    Point3DCL exp_RadDrop;                         ///< radii of the ellipsoidal droplet
    Point3DCL exp_PosDrop;                         ///< position of the droplet
    double    exp_InflowVel,                       ///< max. inflow velocity (parabolic profile)
              exp_RadInlet;                        ///< radius at inlet of measuring device
    int       exp_FlowDir;                         ///< flow direction (x/y/z = 0/1/2)
    Point3DCL exp_Gravity;                         ///< gravity
    double    exp_InflowFreq,                      ///< inflow frequence
              exp_InflowAmpl;                      ///< inflow amplitude
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
    int sft_VarTension;                          ///< variable surface tension: 0= off, 1=on
    double sft_SurfTension,                      ///< surface tension coefficient
           sft_JumpWidth,                        ///< jump width
           sft_RelPos,                           ///< position of the jump between upper edge (lambda=0) and barycenter (lambda=1)
           sft_DirtFactor;                       ///< red. factor for surface tension (due to contamination) in lower droplet)
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
    int    trp_DoTransp;                     ///< mass transport on (1) or off (0)
    double trp_Theta;                        ///< time integration theta scheme
    int    trp_Iter;
    double trp_Tol,
           trp_DiffPos,                      ///< diffusion coefficient (pos. part)
           trp_DiffNeg,                      ///< diffusion coefficient (neg. part)
           trp_H,                            ///< Henry number cneg(\f$\infty\f$) = H*cpos(\f$\infty\f$)
           trp_IniCPos,                      ///< initial concentration (pos. part)
           trp_IniCNeg,                      ///< initial concentration (neg. part)
           trp_NitschePenalty,
           trp_NitscheXFEMStab;              ///< threshold for discarding ext. dofs parameter, default 0.1
  //@}
  public:
    ParamTransportCL()                        { RegisterParams(); }
    ParamTransportCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for surfactant transport
class ParamSurfactantTransportCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  ///\name Surfactant Transport
  //@{
    int    surf_DoTransp;  ///< surfactant transport on (1) or off (0) 0.5=CrankNicholson
    double surf_Theta;     ///< 0=FwdEuler, 1=BwdEuler, 0.5=CrankNicholson
    int    surf_Iter;      ///< iterations of solver for surfactant equation
    double surf_Tol;       ///< tolerance of solver for surfactant equation
    double surf_OmitBound; ///< omit dof with small mass on the interface
    double surf_Visc;      ///< diffusion coefficient on the interface
  //@}

  public:
    ParamSurfactantTransportCL () { RegisterParams(); }
    ParamSurfactantTransportCL (const string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams (file);
    }
};

/// \brief Parameter class for domain conditions like boundary cond. and geometry
class ParamDomainCondCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name reparametrize level-set function
  //@{
    int    dmc_InitialCond,                         ///< initial condition (0=Zero, 1/2= stat. flow with/without droplet, -1= read from file)
           dmc_GeomType;                            ///< specifies the used geometry (0=ReadMeshBuilder, 1=BrickBuilder)
    string dmc_InitialFile,                         ///< file prefix when reading data for initial condition
           dmc_MeshFile;                            ///< mesh file (created by GAMBIT, FLUENT/UNS format) or dimensions of a cuboid (e.g. 2x3x4\@5x6x7)
    int    dmc_BoundaryType;                        ///< boundary type: 1= hom. Dirichlet, 2= inhom. Dirichlet bnd-data for in-/outflow, 3= tube/canal
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
        public ParamSurfactantTransportCL,
        public ParamLevelSetCL,
        public ParamCouplingCL,
        public ParamMaterialDataCL,
        public ParamExperimentalDataCL,
        public ParamSurfaceTensionCL,
        public ParamAdaptRefCL,
        public ParamNavStokesCL,
        public ParamDomainCondCL,
        public ParamMGSerCL,
        public ParamVTKCL

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
        ParamSurfactantTransportCL::rp_.ReadParams( file);
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
    int    stk_StokesMethod;                        // solver for the Stokes problems
    double stk_InnerTol, stk_OuterTol,              // Parameter der Loeser
           ns_Tol, ns_Reduction, ns_Nonlinear,      // fuer Flow & Levelset
           lvs_Tol, lvs_SD, cpl_Tol;
    int    stk_InnerIter, stk_OuterIter,
           ns_Iter, lvs_Iter;
    int    stk_PcAIter;                             // max. number of iterations for the preconditionier
    double stk_PcATol,                              // tolerance for the preconditioner
           stk_PcSTol;
    int    cpl_Iter;                                // Kopplung Levelset/Flow: Anzahl Fixpunkt-Iterationen
    double stk_XFEMStab;                           ///< threshold for discarding ext. dofs parameter, default 0.1

    double tm_StepSize;                             // Zeitschrittweite
    int    tm_NumSteps;                             // Anzahl Zeitschritte
    double stk_Theta, lvs_Theta;                    // 0=FwdEuler, 1=BwdEuler, 0.5=CN

    double mat_SurfTension,                         // Oberflaechenspannung
           lvs_CurvDiff,                            // num. Glaettung Kruemmungstermberechnung
           mat_DensFluid, mat_DensGas,              // Stoffdaten: Dichte/Viskositaet
           mat_ViscFluid, mat_ViscGas,
           mat_SmoothZone,                          // Glaettungszone fuer Dichte-/Viskositaetssprung
           exp_PumpAmpl, exp_PumpFreq,              // Frequenz und Amplitude der Anregung
           exp_Ampl_zDir;                           // Amplitude in z-Richtung der initialen Phasengrenze

    Point3DCL exp_Gravity;                          // Schwerkraft
    double    exp_Thickness;                        // Filmdicke
    Point3DCL mcl_MeshResolution,                   // Gitteraufloesung und
              mcl_MeshSize;                         // Gittergroesse in x-/y-/z-Richtung
    int       lvs_VolCorrection,                    // Volumenkorrektur (0=false)
              mcl_InitialCond;                      // Anfangsbedingung (0=Null, 1= stat. flow, -1= read from file )

    int    ref_FinestLevel, ref_Freq,               // Parameter fuer
           ref_CoarsestLevel;
    double ref_Width;                               // adaptive Verfeinerung

    int    rpm_Freq, rpm_Method;                    // Parameter fuer Reparametrisierung

    string mcl_EnsightCase,                         // Ensight Case,
           mcl_EnsightDir,                          // lok.Verzeichnis, in das die geom/vec/scl-files abgelegt werden
           mcl_InitialFile,
           mcl_BndCond,
           mcl_SerializationFile,                   ///< writes multigrid to serialisation files, special value "none" to ignore
           mcl_DeserializationFile;                 ///< reads multigrid from deserialization files, special value "none" to ignore

    ParamFilmCL()                        { RegisterParams(); }
    ParamFilmCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

} // end of namespace DROPS

#endif





