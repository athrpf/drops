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

/// \brief Parameter class for load balancing parameter
class ParamLoadBalCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Load balancing
  //@{
    double quality;                             ///< quality-parameter for ParMetis
  //@}

  public:
    ParamLoadBalCL()                        { RegisterParams(); }
    ParamLoadBalCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
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
    int markall;                ///< number of refinements of all tetrahedra
    int markdrop;               ///< number of refinements of the drop
    int markcorner;             ///< number or refinements around the corner (0,0,0)
    int markingproc;            ///< number of the proc that should mark its tetrahedra
    int Strategy;               ///< Strategy for refinement
  //@}
  /// \name Coarsening
  //@{
    int coarsedrop;             ///< number of coarsening tetrahedra of the drop
    int coarseall;              ///< number of coarsening all tetrahedra
    int unmarkingproc;          ///< number of the proc that should unmark its tetrahedra
  //@}
  /// \name Load balancing
  //@{
    int refineStrategy;         ///< Algorithm to calculate the load balancing graph for refinement
    int coarseStrategy;         ///< Algorithm to calculate the load balancing graph for coarsening
    int middleMig;              ///< Do a load balancing step with adaptive graph partitioning after refinement and before coarsening
  //@}
  /// \name Misc
  //@{
    int printSize;              ///< Should the size of the mg for all MG be displayed
    int printPMG;               ///< Should the debug info of the ParMultiGridCL printed into a file
    int printGEO;               ///< Should the mesh be printed in Geomview format
    int printTime;              ///< Should the time for each part be printed
    int checkRef;               ///< Check parallel multigrid after every refine
    int checkMig;               ///< Check parallel multigrid after every migrate
    int checkDDD;               ///< do a DDD-GCC (Global Consistency Check)
    string init_pre;            ///< prefix of files, where the serialized multigrid can be found
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
        public virtual ParamVTKCL,
        public virtual ParamStokesCL
{
  private:
    void RegisterParams();

  public:
  /// \name Stokes coefficients
  //@{
    double nu;                                  ///< Coeff of the diffusion
  //@}
  /// \name Refining
  //@{
    int basicref_x, basicref_y, basicref_z;     ///< number of basic refinements of the brick
    double dx, dy, dz;                          ///< Dimension of the brick
    int refall;                                 ///< number of refinements of all tetrahedra
  //@}
  /// \name Load balancing
  //@{
    int refineStrategy;                         ///< Algorithm to calculate the load balancing graph for refinement
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
        ParamStokesCL::rp_.ReadParams( file);
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
    double    frequence;                        ///< frequency of inflow
    double    ampl;                             ///< amplitude of inflow
  //@}
  /// \name Time
  //@{
    int timesteps;                              ///< number of time steps
    double theta;                               ///< for theta-scheme
    double stepsize;                            ///< time step size
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
  /// \name Poisson Coefficients
  //@{
    double nu;                                  ///< Coeff of the diffusion
  //@}
  /// \name Refining
  //@{
    int refall;                                 ///< number of refinements of all tetrahedra or number of adaptive steps
    int markdrop;                               ///< number of refinements of tetrahedra around drop
    int markcorner;                             ///< number of refinements of tetrahedra around corner
    int adaptiv;                                ///< Should the adaptive-strategy be used
  //@}
  /// \name Load balancing
  //@{
    int refineStrategy;                         ///< Algorithm to calculate the load balancing graph for refinement
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
    int useMGS;                                 ///< if solver=GMRES, use the modified Gram-Schmidt
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
    int check;                                  ///< check multigrid for sanity
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
    int refall;                                 ///< number of refinements of all tetrahedra or number of adaptive steps
    int refineStrategy;                         ///< Algorithm to calculate the load balancing graph for refinement
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
    double reduction;                           ///< reduction within adaptive fixed point iteration
    int    nav_iter;                            ///< maximal iterations of fixed point steps
    double nav_tol;                             ///< tolerance for fixed point iteration
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
        public ParamQuadCL,
        public ParamInfoCL
{
  private:
    void RegisterParams();
    typedef ParamMesszelleNsCL base;

  public:
    ParamParBrickFlowCL()
      : ParamMesszelleNsCL(), ParamLoadBalCL(), ParamQuadCL(), ParamInfoCL()
      { RegisterParams(); }

    ParamParBrickFlowCL( const string& filename) : ParamMesszelleNsCL()
       { RegisterParams();
         std::ifstream file(filename.c_str());
         ParamMesszelleNsCL::rp_.ReadParams( file);
         ParamLoadBalCL::rp_.ReadParams( file);
         ParamVTKCL::rp_.ReadParams( file);
         ParamQuadCL::rp_.ReadParams( file);
         ParamInfoCL::rp_.ReadParams( file);
       }
};

/// \brief Parameter class for the problem case TestMzellePar
class ParParamMesszelleNsCL :
        public ParamMesszelleNsCL,
        public ParamLoadBalCL,
        public ParamQuadCL
{
  private:
    void RegisterParams();

  public:
    ParParamMesszelleNsCL()
      : ParamMesszelleNsCL(), ParamLoadBalCL(), ParamQuadCL()
        { RegisterParams(); }

    ParParamMesszelleNsCL( const string& filename) : ParamMesszelleNsCL()
       { RegisterParams();
         std::ifstream file(filename.c_str());
         rp_.ReadParams( file);
         ParamMesszelleNsCL::rp_.ReadParams( file);
         ParamLoadBalCL::rp_.ReadParams( file);
         ParamQuadCL::rp_.ReadParams( file);
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
    int nx;                                     ///< Number of grid points in x direction
    int ny;                                     ///< Number of grid points in y direction
    int nz;                                     ///< Number of grid points in z direction
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
    int    UseMGS;                              ///< use standard or modified Gram-Schmidt
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
    int markall;                ///< number of refinements of all tetrahedra
    int markdrop;               ///< number of refinements of the drop
    int markcorner;             ///< number or refinements around the corner (0,0,0)
    int markingproc;            ///< number of the proc that should mark its tetrahedra
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

} // end of namespace DROPS

#endif
