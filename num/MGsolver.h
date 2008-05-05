//**************************************************************************
// File:    MGsolver.h                                                     *
// Content: classes that constitute the poisson-problem with MG-solver     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 16 2001                                         *
//**************************************************************************

#ifndef DROPS_MGSOLVER_H
#define DROPS_MGSOLVER_H

#include "misc/problem.h"
#include "num/solver.h"

#include <list>


namespace DROPS
{



struct MGLevelDataCL  // data for one triang level
{
    IdxDescCL Idx;    // index description
    MatDescCL A, P;   // stiffness matrix / prolongation
};

class MGDataCL : public std::list<MGLevelDataCL>
{
  public:
    MGDataCL() {};
    void RemoveCoarseResetFinest() {
        if (this->empty()) return;
        //RemoveCoarse
        MGDataCL::iterator it=this->end();
        --it;
        this->erase(this->begin(), it);
        //ResetFinest
        this->begin()->A.Data.clear();
        this->begin()->P.Data.clear();
    }
};

typedef MGDataCL::iterator       MGDataIterCL;
typedef MGDataCL::const_iterator const_MGDataIterCL;


// Multigrid method, V-cycle, beginning from level 'fine'
// numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
//      after 'numLevel' visited levels or if #Unknowns <= 'numUnknDirect'
// If one of the parameters is -1, it will be neglected.
// If the coarsest level 'begin' has been reached, the direct solver is used too.
// NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)
template<class SmootherCL, class DirectSolverCL>
void MGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b,
          const SmootherCL& Smoother, Uint smoothSteps,
          DirectSolverCL& Solver, int numLevel, int numUnknDirect);


// checks  A_coarse= PT * A_fine * P on each level
void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end);


// Uses MGM for solving to tolerance tol or until maxiter iterations are reached.
// The error is measured as two-norm of dx for residerr=false, of Ax-b for residerr=true.
template<class SmootherCL, class DirectSolverCL>
void MG(const MGDataCL& MGData, const SmootherCL&, DirectSolverCL&, VectorCL& x, const VectorCL& b,
        int& maxiter, double& tol, const bool residerr= true, Uint sm=1, int lvl=-1);


template<class SmootherT, class DirectSolverT>
class MGSolverCL : public SolverBaseCL
{
  private:
    const MGDataCL&      mgdata_;
    const SmootherT&     smoother_;
    DirectSolverT&       directSolver_;
    const bool residerr_;
    Uint smoothSteps_;
    Uint usedLevels_;

  public:
    MGSolverCL( const MGDataCL& mgdata, const SmootherT& sm, DirectSolverT& ds, int maxiter, double tol,
               const bool residerr= true, Uint smsteps= 1, int lvl= -1 )
        : SolverBaseCL(maxiter,tol), mgdata_(mgdata), smoother_(sm), directSolver_(ds),
          residerr_(residerr), smoothSteps_(smsteps), usedLevels_(lvl) {}

    void Solve(const MatrixCL& /*A*/, VectorCL& x, const VectorCL& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MG( mgdata_, smoother_, directSolver_, x, b, _iter, _res, residerr_, smoothSteps_, usedLevels_);
    }
};

//===================================
// definition of template functions
//===================================

template<class SmootherCL, class DirectSolverCL>
void
MGM(const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b,
    const SmootherCL& Smoother, Uint smoothSteps,
    DirectSolverCL& Solver, int numLevel, int numUnknDirect)
// Multigrid method, V-cycle. If numLevel==0 or #Unknowns <= numUnknDirect,
// the direct solver Solver is used.
// If one of the parameters is -1, it will be neglected.
// If MGData.begin() has been reached, the direct solver is used too.
{
    const_MGDataIterCL coarse= fine;
    --coarse;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) )
       || fine==begin)
    { // use direct solver
        // We use relative residual-measurement, as otherwise the accuracy and
        // correctness depend on the scaling of the matrix and geometry.
        // This has bitten us in e.g. levelset/mzelle_instat.cpp.
        const double r0= norm( b - fine->A.Data*x);
        Solver.SetTol( 1e-5*r0);
        Solver.Solve( fine->A.Data, x, b);
//        const double r1= norm( b - fine->A.Data*x);
//        std::cerr << "MGM: direct solver: iterations: " << Solver.GetIter()
//                  << "\treduction: " << r1/r0 << '\n';
        return;
    }
    VectorCL d(coarse->Idx.NumUnknowns), e(coarse->Idx.NumUnknowns);
    // presmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( fine->A.Data, x, b);
    // restriction of defect
    d= transp_mul( fine->P.Data, VectorCL( b - fine->A.Data*x));
    // calculate coarse grid correction
    MGM( begin, coarse, e, d, Smoother, smoothSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    // add coarse grid correction
    x+= fine->P.Data * e;
    // postsmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( fine->A.Data, x, b);
}


template<class SmootherCL, class DirectSolverCL>
void MG(const MGDataCL& MGData, const SmootherCL& smoother, DirectSolverCL& solver, 
        VectorCL& x, const VectorCL& b,
        int& maxiter, double& tol, const bool residerr, Uint sm, int lvl)
{
    const_MGDataIterCL finest= --MGData.end();
    double resid= -1;
    double old_resid;
    VectorCL tmp;
    if (residerr == true) {
        resid= norm( b - finest->A.Data * x);
        //std::cerr << "initial residual: " << resid << '\n';
    }
    else
        tmp.resize( x.size());

    int it;
    for (it= 0; it<maxiter; ++it) {
        if (residerr == true) {
            if (resid <= tol) break;
        }
        else tmp= x;
        MGM( MGData.begin(), finest, x, b, smoother, sm, solver, lvl, -1);
        if (residerr == true) {
            old_resid= resid;
            resid= norm( b - finest->A.Data * x);
//            std::cerr << "iteration: " << it  << "\tresidual: " << resid;
//            std::cerr << "\treduction: " << resid/old_resid;
//            std::cerr << '\n';
        }
        else if ((resid= norm( tmp - x)) <= tol) break;
    }
    maxiter= it;
    tol= resid;
}

} // end of namespace DROPS

#endif
