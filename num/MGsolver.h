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

typedef std::list<MGLevelDataCL> MGDataCL;
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
void MG(const MGDataCL& MGData, VectorCL& x, const VectorCL& b, 
        int& maxiter, double& tol, const bool residerr= true);


class MGPreCL
{
  private:
    MGDataCL& A_;
    Uint iter_;

  public:
    MGPreCL( MGDataCL& A, Uint iter)
        :A_( A), iter_( iter)
    {}

    template <class Mat, class Vec>
    void
    Apply( const Mat&, Vec& x, const Vec& r) const;
};


// MG
class MGSolverCL : public SolverBaseCL
{
  private:
    const MGDataCL& _mgdata;

  public:
    MGSolverCL( const MGDataCL& mgdata, int maxiter, double tol )
        : SolverBaseCL(maxiter,tol), _mgdata(mgdata) {}

    void Solve(const MatrixCL& /*A*/, VectorCL& x, const VectorCL& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MG( _mgdata, x, b, _iter, _res);
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
        Solver.Solve( fine->A.Data, x, b);
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


template <class Mat, class Vec>
void
MGPreCL::Apply( const Mat&, Vec& x, const Vec& r) const
{
    Uint sm=  2; // how many smoothing steps?
    int lvl= -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
    SSORsmoothCL smoother( omega);  // Symmetric-Gauss-Seidel with over-relaxation
    SSORPcCL directpc; PCG_SsorCL solver( directpc, 200, 1e-12);
    for (Uint i= 0; i < iter_; ++i)
        MGM( A_.begin(), --A_.end(), x, r, smoother, sm, solver, lvl, -1);
}

} // end of namespace DROPS

#endif
