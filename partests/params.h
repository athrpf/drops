/***************************************************************************
*  File:    params.h                                                       *
*  Content: Reading of Parameterfiles for parallel programms               *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           10.11.2005                                             *
*  last modified:   10.11.2005                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file partests/params.h
/// \brief Reading parameters for TestRefPar, TestStokesPar, TestPoissonPar,
///  TestInstatStokesPar and TestExchange

#ifndef DROPS_PAR_PARAMS_H
#define DROPS_PAR_PARAMS_H

#include "misc/params.h"
#include "levelset/params.h"

namespace DROPS
{

/// \brief Paramater class for loadbalancing parameter
class ParamLoadBalCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Loadbalancing
  //@{
    int refineStrategy;                         ///< Algorithm of calculate the loabal-graph for refinement
    double quality;                             ///< quality-parameter for ParMetis
  //@}

  public:
    ParamLoadBalCL()                        { RegisterParams(); }
    ParamLoadBalCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};


/// \brief Paramater class for ensight parameter
class ParamEnsightCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Ensight
  //@{
    int    ensight;                             ///< Ensight output
    string ensCase;                             ///< name of Ensight Case
    string ensDir;                              ///< local directory for Ensight files
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
    int printNumUnk;            ///< Print number of unknwons
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
    int quad;                                   ///< Write out qudrilateral grid
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

/// \brief Parameter class for serialization of a parallel multigrid
class ParamMGSerCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name serialization of a multigrid
  //@{
    int serialization;     /// Perform serialization
    int overwrite;         /// Overwrite old output-files or create new for every step.
    string ser_dir;        /// Output-Dir.
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
    int RepFreq;                ///< number of timesteps before reparametrize the level-set function
    int RepMethod;              ///< method of reparametrize level-set function
    double NarrowBand;          ///< Narrow-Band method for the Euclidian method (e.g. 0.5 : all dof with <50% of the maximal level-set value are considered)
  //@}

  public:
    ParamReparamCL()                        { RegisterParams(); }
    ParamReparamCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for adaptive refinement (with level-set function)
class ParamApdaptRefCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name reparametrize level-set function
  //@{
    int ref_freq;               ///< number of timesteps before adaptive refinement
    int ref_flevel;             ///< finest level in the near of the phase boundary
    double ref_width;           ///< domain of refinement
  //@}

  public:
    ParamApdaptRefCL()                        { RegisterParams(); }
    ParamApdaptRefCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the problem case TestRefPar
class ParamParRefCL :
        public virtual ParamBaseCL,
        public virtual ParamBrickCL
{
  private:
    void RegisterParams();

  public:
  /// \name Refining
  //@{
    int init_cond;              ///< init multigrid as brick (0) or read from serialization file
    int refined;                ///< number of refinements, that has been done to create multigrid in case of init_cond==1
    int markall;                ///< number of refinements of all tetraeder
    int markdrop;               ///< number of refinements of the drop
    int markcorner;             ///< number or refinements around the corner (0,0,0)
    int markingproc;            ///< number of the proc that should mark its tetraeder
    int Strategy;               ///< Strategy for refinement
  //@}
  /// \name Coarsening
  //@{
    int coarsedrop;             ///< number of coarsening tetraeder of the drop
    int coarseall;              ///< number of coarsening all tetraeder
    int unmarkingproc;          ///< number of the proc that should unmark its tetraeder
  //@}
  /// \name Loadbalancing
  //@{
    int refineStrategy;         ///< Algorithm of calculate the loabal-graph for refinement
    int coarseStrategy;         ///< Algorithm of calculate the loabal-graph for coarsening
    int middleMig;              ///< Do a loadbalancing step with adaptive graph partitioning after refinement and before coarsening
  //@}
  /// \name Misc
  //@{
    int printSize;              ///< Should the size of the mg for all MG be displayed
    int printPMG;               ///< Should the debuginfo of the ParMultiGridCL printed into a file
    int printGEO;               ///< Should the mesh be printed in Geomview format
    int printTime;              ///< Should the time for each part be printed
    int checkRef;               ///< Check parallel multigrid after every refine
    int checkMig;               ///< Check parallel multigrid after every migrate
    int checkDDD;               ///< do a DDD-GCC (Global Consistency Check)
    string init_pre;            ///< prefix of files, where the serializated multigrid can be found
  //@}

    ParamParRefCL()                        { RegisterParams(); }
    ParamParRefCL( const string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamBrickCL::rp_.ReadParams( file);
    }
};

/// \brief Parameter class for the problem case TestStokesPar
class ParamParStokesCL :
        public virtual ParamBaseCL,
        public virtual ParamEnsightCL,
        public virtual ParamVTKCL
{
  private:
    void RegisterParams();

  public:
  /// \name Stokes-Coeffs
  //@{
    double nu;                                  ///< Coeff of the diffusion
  //@}
  /// \name Refining
  //@{
    int basicref_x, basicref_y, basicref_z;     ///< number of basic refinements of the brick
    double dx, dy, dz;                          ///< Dimension of the brick
    int refall;                                 ///< number of refinements of all tetraeder
  //@}
  /// \name Loadbalancing
  //@{
    int refineStrategy;                         ///< Algorithm of calculate the loabal-graph for refinement
  //@}
  /// \name Solver
  //@{
    double relax;                               ///< Relax-Parameter for preconditioner
    int    pc_iter;                             ///< max iterations for solver within the approximate schur complement matrix
    double pc_rel_tol;                          ///< relative tolerance for solver within the approximate schur complement matrix
    int    inner_iter;                          ///< max iterations for the inner solver
    int    outer_iter;                          ///< max iterations for the outer solver
    double inner_tol;                           ///< toleration for the inner solver
    double outer_tol;                           ///< toleration for the outer solver
    int    restart;                             ///< if inner solver==GMRES, dimension of Krylov subspace, solver==GCR , truncation
    int    relative;                            ///< meassure resid relative
    double reduction;                           ///< residual reduction for innexact Uzawa algorithm
    int    accur;                               ///< use accurate solver
  //@}
  /// \name Misc
  //@{
    int printInfo;                              ///< Display information
  //@}

    ParamParStokesCL()                        { RegisterParams(); }
    ParamParStokesCL( const string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamVTKCL::rp_.ReadParams( file);
        ParamEnsightCL::rp_.ReadParams( file);
    }
};

/// \brief Parameter class for the problem case TestInstatStokesPar
class ParamParInstatStokesCL : public virtual ParamParStokesCL
{
 private:
    void RegisterParams();

 public:
  /// \name Experiment Data
  //@{
    double    inflowVel;                        ///< max. inflow velocity (parabolic profile)
    Point3DCL g;                                ///< gravity
    double    frequence;                        ///< frequence of inflow
    double    ampl;                             ///< amplitute of inflow
  //@}
  /// \name Time
  //@{
    int timesteps;                              ///< number of timesteps
    double theta;                               ///< for theta-scheme
    double stepsize;                            ///< timestep-size
  //@}

    ParamParInstatStokesCL() : ParamParStokesCL()
        { RegisterParams(); }
    ParamParInstatStokesCL( const string& filename) : ParamParStokesCL()
        { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the problem case TestPoissonPar
class ParamParPoissonCL : public ParamBrickCL
{
  private:
    void RegisterParams();

  public:
  /// \name Poisson-Coeffs
  //@{
    double nu;                                  ///< Coeff of the diffusion
  //@}
  /// \name Refining
  //@{
    int refall;                                 ///< number of refinements of all tetraeder or number of adaptive steps
    int markdrop;                               ///< number of refinements of tetraeder around drop
    int markcorner;                             ///< number of refinements of tetraeder around corner
    int adaptiv;                                ///< Should the adaptive-strategy be used
  //@}
  /// \name Loadbalancing
  //@{
    int refineStrategy;                         ///< Algorithm of calculate the loabal-graph for refinement
    int transferUnks;                           ///< Should the unknowns be transfered too
  //@}
  /// \name Solver
  //@{
    int solver;                                 ///< Solver
    int precond;                                ///< Preconditioner for solver
    double relax;                               ///< Relax-Parameter for preconditioner
    int pciter;                                 ///< number of maximal iterations for BiCGStab as preconditioner
    double pctol;                               ///< tolerance of BiCGStab as preconditioner
    int iter;                                   ///< max iterations for the solver
    double tol;                                 ///< toleration for the solver
    int restart;                                ///< if solver==GMRES, dimension of Krylov subspace, solver==GCR, truncation
    int useMGS;                                 ///< if solver=GMRES, use the modified Gramm-Schmidt
    int relative;                               ///< use relative resid
    int accur;                                  ///< use accur variant of solver, if exists
    int modified;                               ///< use modified variant for better scalability
    int preCondMeth;                            ///< use left or right preconditioning for GMRES
  //@}
  /// \name Misc
  //@{
    int printGEO;                               ///< print the multigrid as geomview file
    int printMG;                                ///< print the multigrid as ascii file
    int printSize;                              ///< print information about the multigrid
    int printUnknowns;                          ///< print information about number of unknowns
    int check;                                  ///< check multigrid for sanety
  //@}
  /// \name Ensight
  //@{
    int ensight;
    string EnsCase;             ///< name of Ensight Case
    string EnsDir;              ///< local directory for Ensight files
    string geomName;            ///< name for the geometry
    string varName;             ///< name of the variable
  //@}

    ParamParPoissonCL() : ParamBrickCL()       { RegisterParams(); }
    ParamParPoissonCL( const string& filename) : ParamBrickCL() {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamBrickCL::rp_.ReadParams(file);
    }
};


/// \brief Parameter class for the problem case TestExchangePar
class ParamParExchangeCL : public ParamBaseCL
{
  private:
    void RegisterParams();

  public:
  /// \name Unknowns
  //@{
    int numsV1,                                 ///< dof on vertices for first index
        numsE1,                                 ///< dof on edges for first index
        numsT1;                                 ///< dof on tetras for first index
    int numsV2,                                 ///< dof on vertices for second index
        numsE2,                                 ///< dof on edges for second index
        numsT2;                                 ///< dof on tetras for second index
  //@}
  /// \name Misc
  //@{
    int printMG;                                ///< print the MultiGrid into file
    int printEx;                                ///< print the ExchangeCL into file
    int timeMeas;                               ///< measure time for accumulation
    int tests;                                  ///< number of runs of accumulate
    int multruns;                               ///< number of runs for checking multiple runs
    int checkMG;                                ///< Check pmg for sanity
    int printMsgSize;                           ///< display the message size between all procs
  //@}
  /// \name Refining
  //@{
    int basicref_x, basicref_y, basicref_z;     ///< number of basic refinements of the brick
    double dx, dy, dz;                          ///< Dimension of the brick
    int refall;                                 ///< number of refinements of all tetraeder or number of adaptive steps
    int refineStrategy;                         ///< Algorithm of calculate the loabal-graph for refinement
    int migEveryTime;                           ///< Migration after every refinement
  //@}
    ParamParExchangeCL()                        { RegisterParams(); }
    ParamParExchangeCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }

};

/// \brief Parameter class for the problem case TestNavStokesPar
class ParamParNavStokesCL : public virtual ParamParInstatStokesCL
{
  private:
    void RegisterParams();

  public:
    /// \name NavierStokes
    //@{
    double reduction;                           ///< reduction within adaptive fixpoint iteration
    int    nav_iter;                            ///< maximal iterations of fixpoint steps
    double nav_tol;                             ///< tolerance for fixpoint iteration
    int    markTop;
    //@}

    ParamParNavStokesCL() : ParamParInstatStokesCL()
        { RegisterParams(); }

    ParamParNavStokesCL( const string& filename) : ParamParInstatStokesCL()
        { RegisterParams();
          std::ifstream file(filename.c_str());
          ParamParInstatStokesCL::rp_.ReadParams( file);
          rp_.ReadParams( file); }
};

/// \brief Parameter class for the problem case TestBrickflowPar
class ParamParBrickFlowCL :
        public ParamMesszelleNsCL,
        public ParamLoadBalCL,
        public ParamEnsightCL,
        public ParamVTKCL,
        public ParamInfoCL,
        public ParamQuadCL,
        public ParamMGSerCL,
        public ParamReparamCL
{
  private:
    void RegisterParams();
    typedef ParamMesszelleNsCL base;

  public:
    ParamParBrickFlowCL()
      : ParamMesszelleNsCL(), ParamLoadBalCL(), ParamEnsightCL(), ParamVTKCL(),
        ParamInfoCL(), ParamQuadCL()
      { RegisterParams(); }

    ParamParBrickFlowCL( const string& filename) : ParamMesszelleNsCL()
       { RegisterParams();
         std::ifstream file(filename.c_str());
         ParamMesszelleNsCL::rp_.ReadParams( file);
         ParamVTKCL::rp_.ReadParams( file);
         ParamLoadBalCL::rp_.ReadParams( file);
         ParamInfoCL::rp_.ReadParams( file);
         ParamQuadCL::rp_.ReadParams( file);
         ParamEnsightCL::rp_.ReadParams( file);
         ParamMGSerCL::rp_.ReadParams( file);
         ParamReparamCL::rp_.ReadParams( file);
       }
};

/// \brief Parameter class for the problem case TestMzellePar
class ParParamMesszelleNsCL :
        public ParamMesszelleNsCL,
        public ParamLoadBalCL,
        public ParamEnsightCL,
        public ParamVTKCL,
        public ParamInfoCL,
        public ParamQuadCL,
        public ParamReparamCL
{
  private:
    void RegisterParams();

  public:
    ParParamMesszelleNsCL()
      : ParamMesszelleNsCL(), ParamLoadBalCL(), ParamEnsightCL(), ParamVTKCL(),
        ParamInfoCL(), ParamQuadCL()
        { RegisterParams(); }

    ParParamMesszelleNsCL( const string& filename) : ParamMesszelleNsCL()
       { RegisterParams();
         std::ifstream file(filename.c_str());
         rp_.ReadParams( file);
         ParamMesszelleNsCL::rp_.ReadParams( file);
         ParamVTKCL::rp_.ReadParams( file);
         ParamLoadBalCL::rp_.ReadParams( file);
         ParamInfoCL::rp_.ReadParams( file);
         ParamQuadCL::rp_.ReadParams( file);
         ParamEnsightCL::rp_.ReadParams( file);
         ParamReparamCL::rp_.ReadParams( file);
       }
};

/// \brief Parameter class for the FilmCalculation
class ParamParFilmCL : public ParamBaseCL
{
  private:
    void RegisterParams();

  public:
  /// \name Geometry
  //@{
    Point3DCL Dim;                              ///< Dimension of the brick
    int nx;                                     ///< Number of gridpoints in x direction
    int ny;                                     ///< Number of gridpoints in y direction
    int nz;                                     ///< Number of gridpoints in z direction
    int Refine;                                 ///< Number of regular refinements
  //@}
  /// \name Material data
  //@{
    double KineticViscosity;                    ///< Convection
    double HeatConductivity;                    ///< Diffusion
    double InflowTemperature;                   ///< Inflow temperature
    double WallTemperature;                     ///< Temperature at the wall
    Point3DCL g;                                ///< Gravity
  //@}

  /// \name Solver parameter
  //@{
    double Relax;                               ///< Relaxation of Jacobi PC
    int    Restart;                             ///< Restart of GMRES
    int    Iter;                                ///< Maximal iterations
    double Tol;                                 ///< Tolerance
    int    UseMGS;                              ///< use standard or modified Gramm-Schmidt
  //@}
  /// \name Ensight
  //@{
    int Ensight;                                ///< Ensight output
    string EnsCase;                             ///< name of Ensight Case
    string EnsDir;                              ///< local directory for Ensight files
    string geomName;                            ///< name for the geometry
    string varName;                             ///< name of the variable
  //@}
    ParamParFilmCL()                        { RegisterParams(); }
    ParamParFilmCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the problem case TestMGSerPar
class ParamParSerCL :
        virtual public ParamMGSerCL,
        virtual public ParamBrickCL,
        virtual public ParamVTKCL
{
  private:
    void RegisterParams();

  public:
  /// \name Refining
  //@{
    int markall;                ///< number of refinements of all tetraeder
    int markdrop;               ///< number of refinements of the drop
    int markcorner;             ///< number or refinements around the corner (0,0,0)
    int markingproc;            ///< number of the proc that should mark its tetraeder
  //@}

    int mode;                   ///< mode=0: write serialized multigrid into a file, mode=1: read serialized multigrid
    int unknowns;               ///< test with or without unknowns

    ParamParSerCL()                        { RegisterParams(); }
    ParamParSerCL( const string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamMGSerCL::rp_.ReadParams( file);
        ParamBrickCL::rp_.ReadParams( file);
        ParamVTKCL::rp_.ReadParams( file);
    }
};

/// \brief Parameter class for the problem cas TestFilmStokesPar.cpp
class ParamParFilmStokesCL :
        virtual public ParamLoadBalCL,
        virtual public ParamEnsightCL,
        virtual public ParamVTKCL,
        virtual public ParamInfoCL,
        virtual public ParamBrickCL,
        virtual public ParamReparamCL,
        virtual public ParamApdaptRefCL
{ // y = Filmnormal, x = Ablaufrichtung
  protected:
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
    int       VolCorr,                          // Volumenkorrektur (0=false)
              IniCond;                          // Anfangsbedingung (0=Null, 1= stat. flow, -1= read from file )

    string IniData;

    ParamParFilmStokesCL() { RegisterParams(); }
    ParamParFilmStokesCL( const string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamLoadBalCL::rp_.ReadParams( file);
        ParamEnsightCL::rp_.ReadParams( file);
        ParamVTKCL::rp_.ReadParams( file);
        ParamInfoCL::rp_.ReadParams( file);
        ParamBrickCL::rp_.ReadParams( file);
        ParamReparamCL::rp_.ReadParams( file);
        ParamApdaptRefCL::rp_.ReadParams( file);
    }
};
} // end of namespace DROPS

#endif
