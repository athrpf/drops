/// \file MGsolver.h
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Maxim Larin, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_MGSOLVER_H
#define DROPS_MGSOLVER_H

#include "misc/problem.h"
#include "num/solver.h"

#include <list>
#include <cstring>

namespace DROPS
{
/**
\brief  Multigrid method, V-cycle, beginning from level 'fine'

 numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
 after 'numLevel' visited levels or if number of unknowns <= 'numUnknDirect'
 If one of the parameters is -1, it will be neglected.
 If the coarsest level 'begin' has been reached, the direct solver is used too.
 NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)

 \param begin         coarsest level
 \param fine          actual level
 \param P             prolongation
 \param x             approximation of the solution
 \param b             right hand side
 \param Smoother      multigrid smoother
 \param smoothSteps   number of smoothing steps
 \param Solver        coarse grid/direct solver with relative residual measurement
 \param numLevel      number of vidited levels
 \param numUnknDirect minimal number of unknowns for the direct solver */
template<class SmootherCL, class DirectSolverCL, class ProlongationIteratorT>
void MGM( const MLMatrixCL::const_iterator& begin, const MLMatrixCL::const_iterator& fine,
          const ProlongationIteratorT& P, VectorCL& x, const VectorCL& b,
          const SmootherCL& Smoother, Uint smoothSteps,
          DirectSolverCL& Solver, int numLevel, int numUnknDirect);

/**
\brief Uses MGM for solving to tolerance tol or until maxiter iterations are reached.

 The error is measured as two-norm of dx for residerr=false, of Ax-b for residerr=true.
 sm controls the number of smoothing steps, lvl the number of used levels */
template<class SmootherCL, class DirectSolverCL, class ProlongationT>
void MG(const MLMatrixCL& MGData, const ProlongationT& Prolong, const SmootherCL&,
        DirectSolverCL&, VectorCL& x, const VectorCL& b, int& maxiter, double& tol,
        const bool residerr= true, Uint sm=1, int lvl=-1);


/*******************************************************************
*   M G S o l v e r  C L                                           *
*******************************************************************/
/// \brief MultiGrid solver for a single matrix problem
/** Uses a Multigrid structure for a single matrix, e.g.
    a poisson problem or the A-block of a (navier-)stokes problem */
/*******************************************************************
*   M G S o l v e r  C L                                           *
********************************************************************/
template<class SmootherT, class DirectSolverT, class ProlongationT= MLMatrixCL>
class MGSolverCL : public SolverBaseCL
{
  private:
    ProlongationT     P;                 ///< prolongation
    const SmootherT&  smoother_;         ///< multigrid smoother
    DirectSolverT&    directSolver_;     ///< coarse grid solver with relative residual measurement
    const bool        residerr_;         ///< controls the error measuring: false : two-norm of dx, true: two-norm of residual
    Uint              smoothSteps_;      ///< number of smoothing steps
    int               usedLevels_;       ///< number of used levels (-1 = all)

  public:
    /// constructor for MGSolverCL
    /** \param sm         multigrid smoother
        \param ds         coarse grid solver with relative residual measurement
        \param maxiter    maximal iteration number
        \param tol        stopping criterion
        \param residerr   controls the error measuring: false : two-norm of dx, true: two-norm of residual
        \param smsteps    number of smoothing steps
        \param lvl        number of used levels (-1 = all) */
    MGSolverCL( const SmootherT& sm, DirectSolverT& ds, int maxiter,
                double tol, const bool residerr= true, Uint smsteps= 1, int lvl= -1 )
        : SolverBaseCL(maxiter,tol), smoother_(sm), directSolver_(ds),
          residerr_(residerr), smoothSteps_(smsteps), usedLevels_(lvl) {}

    ProlongationT* GetProlongation() { return &P; }
    /// solve function: calls the MultiGrid-routine
    void Solve(const MLMatrixCL& A, VectorCL& x, const VectorCL& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MG( A, P, smoother_, directSolver_, x, b, _iter, _res, residerr_, smoothSteps_, usedLevels_);
    }
    void Solve(const MatrixCL&, VectorCL&, const VectorCL&)
    {
        throw DROPSErrCL( "MGSolverCL::Solve: need multilevel data structure\n");
    }

};

/// checks multigrid structure
void CheckMGData( const MLMatrixCL& A, const MLMatrixCL& P);

/**
\brief Multigrid method for saddle point problems, V-cycle, beginning from level 'fine'

 numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
 after 'numLevel' visited levels or if number of unknowns <= 'numUnknDirect'
 If one of the parameters is -1, it will be neglected.
 If the coarsest level 'begin' has been reached, the direct solver is used too.
 NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)
 \param beginA        A on coarsest level
 \param fineA         A on actual level
 \param fineB         B on actual level
 \param fineBT        B^T on actual level
 \param fineprM       pressure mass matrix on actual level
 \param PVel          prolongation for P2
 \param PPr           prolongation for P1
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
template<class StokesSmootherCL, class StokesDirectSolverCL, class ProlongItT1, class ProlongItT2>
void StokesMGM( const MLMatrixCL::const_iterator& beginA,  const MLMatrixCL::const_iterator& fineA,
                const MLMatrixCL::const_iterator& fineB,   const MLMatrixCL::const_iterator& fineBT,
                const MLMatrixCL::const_iterator& fineprM, const ProlongItT1& PVel,
                const ProlongItT2& PPr, VectorCL& u, VectorCL& p, const VectorCL& b,
                const VectorCL& c, const StokesSmootherCL& Smoother, Uint smoothSteps, Uint cycleSteps,
                StokesDirectSolverCL& Solver, int numLevel, int numUnknDirect);

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
    void DiagSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v, size_t id2) const;
    template <typename Mat, typename Vec>
    void GSSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;
    template <typename Mat, typename Vec>
    void LRSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;

    const MLIdxDescCL* idx_;   ///< If not zero, the standard- and the extended- pressure-dof in a vertex are treated as a block system (only for the DiagSmoother.

  public:
/**
\brief constructor of the PVankaSmootherCL
 \param vanka_method solving methods for local problems: 0, 3, 4 schur method with diagonal/Gauss-Seidel/symmetrical Gauss-Seidel approximation of A^{-1}, 2: LR decomposition of the local problems
 \param tau          relaxation parameter
 \param idx          optional extended index for 2x2-systems with pressure & extended pressure */
    PVankaSmootherCL(int vanka_method=0, double tau=1.0, const MLIdxDescCL* idx= 0)
        : vanka_method_(vanka_method), tau_(tau), idx_( idx) {}
    template <typename Mat, typename Vec>
    void Apply( const Mat& A, const Mat& B, const Mat& BT, const Mat&,
                Vec& u, Vec& p, const Vec& f, const Vec& g ) const;
    void SetRelaxation (double tau) { tau_= tau; }
    void SetVankaMethod( int method) {vanka_method_ = method;}  ///< change the method for solving the local problems
    int  GetVankaMethod()            {return vanka_method_;}    ///< get the number of the method for solving the local problems

    void Setidx (const MLIdxDescCL* idx) { idx_= idx; }
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
