/// \file
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Maxim Larin, Patrick Esser, IGPM

#ifndef DROPS_MGSOLVER_H
#define DROPS_MGSOLVER_H

#include "misc/problem.h"
#include "num/solver.h"

#include <list>
#include <cstring>

namespace DROPS
{
/*******************************************************************
*   M G L e v e l D a t a  C L                                     *
*******************************************************************/
/// \brief Represents the data for one level of the triangulation
/** Contains the IdxDescCL for velocity and pressure
    It also stores the matrices of the generalized Stokes problem
    including mass matrices, saddle-point matrices (Stokes/NavStokes)
    and a matrix pointer to the actual A-block*/
/*******************************************************************
*   M G L e v e l D a t a  C L                                     *
*******************************************************************/
struct MGLevelDataCL  // data for one triang level
{
    IdxDescCL Idx,        ///< index description for velocity
              IdxPr;      ///< index description for pressure
    MatDescCL A,          ///< Stokes A matrix
              B,          ///< Stokes B matrix
              BT;         ///< Stokes B^T matrix
    MatDescCL Mpr,        ///< pressure mass matrix
              Mvel;       ///< velocity mass matrix
    MatDescCL P,          ///< prolongation matrix for velocity
              PPr;        ///< prolongation matrix pressure
    MatDescCL AN;         ///< Navier Stokes A-Block matrix (A + alpha*N)
    MatrixCL* ABlock;     ///< pointer to the matrix which is used by MultiGrid methods
};


/*******************************************************************
*   M G D a t a  C L                                               *
*******************************************************************/
/// \brief Represents the data for all levels of the triangulation
/** Contains a list of MGLevelDataCL
    NOTE: Assumes, that the levels are stored in an ascending order
    (first=coarsest, last=finest) */
/*******************************************************************
*   M G D a t a  C L                                               *
********************************************************************/
class MGDataCL : public std::list<MGLevelDataCL>
{
  private:
    bool StokesMG_;                                                      ///< internal variable which allows the using of B/B^T etc

  public:
    MGDataCL(int n=0) : std::list<MGLevelDataCL>(n), StokesMG_(false) {} ///< creates a new MGDataCL with n entrys
    void RemoveCoarseResetFinest();                                      ///< removes all entrys exept the last
    bool StokesMG() {return StokesMG_;}                                  ///< returns an internal variable which allows the using of B/B^T etc
    void SetStokesMG(bool full) {StokesMG_=full;}                        ///< sets an internal variable which allows the using of B/B^T etc
};

typedef MGDataCL::iterator       MGDataIterCL;                           ///< iterator for the MGDataCL
typedef MGDataCL::const_iterator const_MGDataIterCL;                     ///< constant iterator for the MGDataCL

/// checks a part of the MGDataCL: A_coarse= PT * A_fine * P on each level
void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end);

/**
\brief  Multigrid method, V-cycle, beginning from level 'fine'

 numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
 after 'numLevel' visited levels or if number of unknowns <= 'numUnknDirect'
 If one of the parameters is -1, it will be neglected.
 If the coarsest level 'begin' has been reached, the direct solver is used too.
 NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)

 \param begin         coarsest level
 \param fine          actual level
 \param x             approximation of the solution
 \param b             right hand side
 \param Smoother      multigrid smoother
 \param smoothSteps   number of smoothing steps
 \param Solver        coarse grid/direct solver with relative residual measurement
 \param numLevel      number of vidited levels
 \param numUnknDirect minimal number of unknowns for the direct solver */
template<class SmootherCL, class DirectSolverCL>
void MGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b,
          const SmootherCL& Smoother, Uint smoothSteps,
          DirectSolverCL& Solver, int numLevel, int numUnknDirect);

/**
\brief Uses MGM for solving to tolerance tol or until maxiter iterations are reached.

 The error is measured as two-norm of dx for residerr=false, of Ax-b for residerr=true.
 sm controls the number of smoothing steps, lvl the number of used levels */
template<class SmootherCL, class DirectSolverCL>
void MG(const MGDataCL& MGData, const SmootherCL&, DirectSolverCL&, VectorCL& x, const VectorCL& b,
        int& maxiter, double& tol, const bool residerr= true, Uint sm=1, int lvl=-1);


/*******************************************************************
*   M G S o l v e r  C L                                           *
*******************************************************************/
/// \brief MultiGrid solver for a single matrix problem
/** Uses a Multigrid structure for a single matrix, e.g.
    a poisson problem or the A-block of a (navier-)stokes problem */
/*******************************************************************
*   M G S o l v e r  C L                                           *
********************************************************************/
template<class SmootherT, class DirectSolverT>
class MGSolverCL : public SolverBaseCL
{
  private:
    const MGDataCL&   mgdata_;           ///< multigrid hierarchy
    const SmootherT&  smoother_;         ///< multigrid smoother
    DirectSolverT&    directSolver_;     ///< coarse grid solver with relative residual measurement
    const bool        residerr_;         ///< controls the error measuring: false : two-norm of dx, true: two-norm of residual
    Uint              smoothSteps_;      ///< number of smoothing steps
    int               usedLevels_;       ///< number of used levels (-1 = all)

  public:
    /// constructor for MGSolverCL
    /** \param mgdata     multigrid hierarchy in ascending order (first=coarsest, last=finest)
        \param sm         multigrid smoother
        \param ds         coarse grid solver with relative residual measurement
        \param maxiter    maximal iteration number
        \param tol        stopping criterion
        \param residerr   controls the error measuring: false : two-norm of dx, true: two-norm of residual
        \param smsteps    number of smoothing steps
        \param lvl        number of used levels (-1 = all) */
    MGSolverCL( const MGDataCL& mgdata, const SmootherT& sm, DirectSolverT& ds, int maxiter,
                double tol, const bool residerr= true, Uint smsteps= 1, int lvl= -1 )
        : SolverBaseCL(maxiter,tol), mgdata_(mgdata), smoother_(sm), directSolver_(ds),
          residerr_(residerr), smoothSteps_(smsteps), usedLevels_(lvl) {}

    /// solve function: calls the MultiGrid-routine
    void Solve(const MatrixCL& /*A*/, VectorCL& x, const VectorCL& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MG( mgdata_, smoother_, directSolver_, x, b, _iter, _res, residerr_, smoothSteps_, usedLevels_);
    }
};

/**
\brief Multigrid method for saddle point problems, V-cycle, beginning from level 'fine'

 numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
 after 'numLevel' visited levels or if number of unknowns <= 'numUnknDirect'
 If one of the parameters is -1, it will be neglected.
 If the coarsest level 'begin' has been reached, the direct solver is used too.
 NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)
 \param begin         coarsest level
 \param fine          actual level
 \param u             velocity
 \param p             pressure
 \param b             rhs for velocity
 \param c             rhs for pressure
 \param Smoother      multigrid smoother
 \param smoothSteps   number of smoothing steps
 \param cycleSteps    number of cycle steps
 \param Solver        coarse grid/direct solver with relative residual measurement
 \param numLevel      number of vidited levels
 \param numUnknDirect minimal number of unknowns for the direct solver */
template<class StokesSmootherCL, class StokesDirectSolverCL>
void StokesMGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& u, VectorCL& p,
                const VectorCL& b, const VectorCL& c,  const StokesSmootherCL& Smoother, Uint smoothSteps,
                Uint cycleSteps, StokesDirectSolverCL& Solver, int numLevel, int numUnknDirect);

/// checks Stokes-MG-Data
void CheckStokesMGData( const_MGDataIterCL begin, const_MGDataIterCL end);

///< dummy SmootherCL for StokesMGM
class DummySmootherCL
{
  public:
    template <typename Mat, typename Vec>
    void Apply( const Mat&, const Mat&, const Mat&, const Mat&,
                Vec&, Vec&, const Vec&, const Vec& ) const
    {}
};

/*******************************************************************
*   P V a n k a S m o o t h e r  C L                               *
*******************************************************************/
/// \brief Vanka smoother for StokesMGM
/** Iterates over the pressure unknowns and construct local saddle point
    problems. These local problems are solved via an schur method with
    different approximations of A^{-1} or via a LR decomposition
    Method from "A comparative study for efficient terative solvers for
    generalized Stokes equations", Larin, Reusken. */
/*******************************************************************
*   P V a n k a S m o o t h e r  C L                               *
********************************************************************/
class PVankaSmootherCL
{
  private:
    int vanka_method_;           ///< method for solving the local problems: 0, 3, 4: schur method with diagonal/Gauss-Seidel/symmetrical Gauss-Seidel approximation of A^{-1}, 2: LR decomposition of the local problems
    double tau_;                 ///< relaxation parameter
    typedef std::vector<size_t>  NodeListVelT;

    template <typename Mat>                 ///< constructs local saddle point matrix
    DMatrixCL<double>  SetupLocalProblem (Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B) const;
    template <typename Mat, typename Vec>
    void DiagSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;
    template <typename Mat, typename Vec>
    void GSSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;
    template <typename Mat, typename Vec>
    void LRSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;

  public:
/**
\brief constructor of the PVankaSmootherCL
 \param vanka_method solving methods for local problems: 0, 3, 4 schur method with diagonal/Gauss-Seidel/symmetrical Gauss-Seidel approximation of A^{-1}, 2: LR decomposition of the local problems
 \param tau          relaxation parameter */
    PVankaSmootherCL(int vanka_method=0, double tau=1.0) : vanka_method_(vanka_method), tau_(tau) {}
    template <typename Mat, typename Vec>
    void Apply( const Mat& A, const Mat& B, const Mat& BT, const Mat&,
                Vec& u, Vec& p, const Vec& f, const Vec& g ) const;
    void SetVankaMethod( int method) {vanka_method_ = method;}  ///< change the method for solving the local problems
    int  GetVankaMethod()            {return vanka_method_;}    ///< get the number of the method for solving the local problems
};

/*******************************************************************
*   B S S m o o t h e r  C L                                       *
*******************************************************************/
/// \brief Braess-Sarazin smoother for StokesMGM
/** Method from "A comparative study for efficient terative solvers for
    generalized Stokes equations", Larin, Reusken. */
/*******************************************************************
*   B S S m o o t h e r  C L                                       *
********************************************************************/
class BSSmootherCL
{
  private:
    int    maxit_;
    double red_, omega_;

  public:
/**
\brief Constructor of the BSSmootherCL
 \param maxit  number of maximal iterations for the inner CG method
 \param red    stopping criterion for the inner CG method
 \param omega  scaling parameter for the diagonal of A */
    BSSmootherCL( int maxit = 20, double red = 2e-1, double omega = 2.0) : maxit_(maxit), red_(red), omega_(omega) {};
    template <typename Mat, typename Vec>
    void Apply( const Mat& A, const Mat& B, const Mat&, const Mat& M,
                Vec& u, Vec& p, const Vec& f, const Vec& g ) const;
};

//===================================
// definition of template functions
//===================================

} // end of namespace DROPS

#include "num/MGsolver.tpp"

#endif
