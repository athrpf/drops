//**************************************************************************
// File:    nssolver.h                                                     *
// Content: nonlinear solvers for the Navier-Stokes equation               *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 20 2001                                           *
//**************************************************************************

#include "num/stokessolver.h"
#include "misc/problem.h"

#ifndef DROPS_NSSOLVER_H
#define DROPS_NSSOLVER_H

namespace DROPS
{

template <class NavStokesT, class SolverT>
class AdaptFixedPtDefectCorrCL
/****************************************************************************
* adaptive fixedpoint defect correction (TUREK p. 187f)
* for the Navier-Stokes equation. The NS problem is of type NavStokesT,
* the inner problems of Stokes-type are solved via a SolverT-solver.
* After the run, the NS class contains the nonlinear part N / cplN belonging
* to the iterated solution.
****************************************************************************/
{
  private:
    NavStokesT& _NS;
    SolverT&    _solver;
    MatrixCL    _AN;
    
    int         _maxiter, _iter;
    double      _tol, _res, _red;
  
  public:
    AdaptFixedPtDefectCorrCL( NavStokesT& NS, SolverT& solver, int maxiter, double tol, double reduction= 100.)
        : _NS( NS), _solver( solver), _maxiter( maxiter), _iter(-1), _tol( tol), _res(-1.), _red( reduction) {}

    ~AdaptFixedPtDefectCorrCL() {}
    
    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    void SetReduction( double red) { _red= red; }
    
    double   GetResid()        const { return _res; }
    int      GetIter ()        const { return _iter; }
    SolverT& GetStokesSolver() const { return _solver; }

    // solves the system   [A + alpha*N] v + BT p = b + alpha*cplN
    //                                 B v        = c
    // (param. alpha is used for time integr. schemes)
    void Solve( const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
                VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha= 1.);
};

template<class NavStokesT, class SolverT>
class FixedPtDefectCorrCL
/****************************************************************************
* fixedpoint defect correction (TUREK p. 187f w/o adaption)
* for the Navier-Stokes equation. The NS problem is of type NavStokesT,
* the inner problems of Stokes-type are solved via a SolverT-solver.
* After the run, the NS class contains the nonlinear part N / cplN belonging
* to the iterated solution.
****************************************************************************/
{
  private:
    NavStokesT& _NS;
    SolverT&    _solver;
    MatrixCL    _AN;
    
    int         _maxiter, _iter;
    double      _tol, _res, _red;
  
  public:
    FixedPtDefectCorrCL( NavStokesT& NS, SolverT& solver, int maxiter, double tol, double reduction= 100.)
        : _NS( NS), _solver( solver), _maxiter( maxiter), _iter(-1), _tol( tol), _res(-1.), _red( reduction) {}

    ~FixedPtDefectCorrCL() {}
    
    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( Uint iter ) { _maxiter= iter; }
    void SetReduction( double red) { _red= red; }
    
    double GetResid()          const { return _res; }
    Uint   GetIter ()          const { return _iter; }
    SolverT& GetStokesSolver() const { return _solver; }

    // solves the system   [A + alpha*N] v + BT p = b + alpha*cplN
    //                                 B v        = c
    // (param. alpha is used for time integr. schemes)
    void Solve( const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
                VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha= 1.);
};


//==================================================
//        derived classes for easier use
//==================================================

template <class NavStokesT>
class AFPDeCo_Uzawa_PCG_CL: public AdaptFixedPtDefectCorrCL<NavStokesT, Uzawa_PCG_CL>
{
  private:
    Uzawa_PCG_CL _uzawaSolver;
  
  public:
    AFPDeCo_Uzawa_PCG_CL( NavStokesT& NS, MatrixCL& M, int fp_maxiter, double fp_tol, int stokes_maxiter,
                          int poiss_maxiter, double poiss_tol, double reduction= 100.)
        : AdaptFixedPtDefectCorrCL<NavStokesT, Uzawa_PCG_CL>( NS, _uzawaSolver, fp_maxiter, fp_tol, reduction),
          _uzawaSolver( M, stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};

template <class NavStokesT>
class AFPDeCo_Schur_PCG_CL: public AdaptFixedPtDefectCorrCL<NavStokesT, Schur_PCG_CL>
{
  private:
    Schur_PCG_CL _schurSolver;
  
  public:
    AFPDeCo_Schur_PCG_CL( NavStokesT& NS, int fp_maxiter, double fp_tol, int stokes_maxiter,
                          int poiss_maxiter, double poiss_tol, double reduction= 100.)
        : AdaptFixedPtDefectCorrCL<NavStokesT, Schur_PCG_CL>( NS, _schurSolver, fp_maxiter, fp_tol, reduction),
          _schurSolver( stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};

template <class NavStokesT>
class FPDeCo_Uzawa_CG_CL: public FixedPtDefectCorrCL<NavStokesT, Uzawa_CG_CL>
{
  private:
    Uzawa_CG_CL _uzawaSolver;
  
  public:
    FPDeCo_Uzawa_CG_CL( NavStokesT& NS, MatrixCL& M, int fp_maxiter, double fp_tol, int stokes_maxiter,
                         int poiss_maxiter, double poiss_tol, double reduction= 100.)
        : FixedPtDefectCorrCL<NavStokesT, Uzawa_CG_CL>( NS, _uzawaSolver, fp_maxiter, fp_tol, reduction),
          _uzawaSolver( M, stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};

template <class NavStokesT>
class FPDeCo_Uzawa_SGSPCG_CL: public FixedPtDefectCorrCL<NavStokesT, Uzawa_SGSPCG_CL>
{
  private:
    Uzawa_SGSPCG_CL _uzawaSolver;
  
  public:
    FPDeCo_Uzawa_SGSPCG_CL( NavStokesT& NS, MatrixCL& M, int fp_maxiter, double fp_tol, int stokes_maxiter,
                         int poiss_maxiter, double poiss_tol, double reduction= 100.)
        : FixedPtDefectCorrCL<NavStokesT, Uzawa_SGSPCG_CL>( NS, _uzawaSolver, fp_maxiter, fp_tol, reduction),
          _uzawaSolver( M, stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};

template <class NavStokesT>
class FPDeCo_Uzawa_PCG_CL: public FixedPtDefectCorrCL<NavStokesT, Uzawa_PCG_CL>
{
  private:
    Uzawa_PCG_CL _uzawaSolver;
  
  public:
    FPDeCo_Uzawa_PCG_CL( NavStokesT& NS, MatrixCL& M, int fp_maxiter, double fp_tol, int stokes_maxiter,
                         int poiss_maxiter, double poiss_tol, double reduction= 100.)
        : FixedPtDefectCorrCL<NavStokesT, Uzawa_PCG_CL>( NS, _uzawaSolver, fp_maxiter, fp_tol, reduction),
          _uzawaSolver( M, stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};

template <class NavStokesT>
class FPDeCo_Schur_PCG_CL: public FixedPtDefectCorrCL<NavStokesT, Schur_PCG_CL>
{
  private:
    Schur_PCG_CL _schurSolver;
  
  public:
    FPDeCo_Schur_PCG_CL( NavStokesT& NS, int fp_maxiter, double fp_tol, int stokes_maxiter,
                         int poiss_maxiter, double poiss_tol, double reduction= 1000.)
        : FixedPtDefectCorrCL<NavStokesT, Schur_PCG_CL>( NS, _schurSolver, fp_maxiter, fp_tol, reduction),
          _schurSolver( stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};


//=================================
//     template definitions
//=================================


template<class NavStokesT, class SolverT>
void
AdaptFixedPtDefectCorrCL<NavStokesT, SolverT>::Solve( 
    const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
    VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha)
{
        VecDescCL v_omw( v.RowIdx);
        VectorCL d( v.Data.size()), e( p.size()),
                 w( v.Data.size()), q( p.size());
        double omega= 1.; // initial value (no damping)
        _iter= 0;
        for(;;++_iter) // ever
        {
            _NS.SetupNonlinear(&_NS.N, &v, &cplN);
            _AN.LinComb( 1., A, alpha, _NS.N.Data);
            
            // calculate defect:
            d= _AN*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
            e= B*v.Data - c;
	    std::cerr << _iter << ": res = " << (_res= sqrt( dot( d, d) + dot( e, e)) ) << std::endl; 
            if (_res < _tol || _iter>=_maxiter)
                break;
            
            // solve correction:
            double outer_tol= _res/_red;
            
            _solver.SetTol( outer_tol);
            _solver.Solve( _AN, B, w, q, d, e);
            
            // calculate adaption:
            _NS.N.Data.clear();
            v_omw.Data= v.Data - omega*w;
            _NS.SetupNonlinear( &_NS.N, &v_omw, &cplN);
            
            d= A*w + alpha*(_NS.N.Data*w) + transp_mul( B, q);
            e= B*w;
            omega= dot( d, VectorCL( A*v.Data + _NS.N.Data*VectorCL( alpha*v.Data)
                + transp_mul( B, p) - b - alpha*cplN.Data))
                + dot( e, VectorCL( B*v.Data - c));
            omega/= norm_sq( d) + norm_sq( e);
            std::cerr << "omega = " << omega << std::endl;
            
            // update solution:
            v.Data-= omega*w;
            p     -= omega*q;
        }
}

template<class NavStokesT, class SolverT>
void FixedPtDefectCorrCL<NavStokesT, SolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
    VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha)
{
        VectorCL d( v.Data.size()), e( p.size()),
                 w( v.Data.size()), q( p.size());
        _iter= 0;
        for(;;++_iter) // ever
        {
            _NS.SetupNonlinear(&_NS.N, &v, &cplN);
            _AN.LinComb( 1., A, alpha, _NS.N.Data);
           
            // calculate defect:
            d= _AN*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
            e= B*v.Data - c;
            
            std::cerr << _iter << ": res = " << (_res= sqrt( norm_sq( d) + norm_sq(e))) << std::endl; 
            if (_res < _tol || _iter>=_maxiter)
                break;
            
            // solve correction:
            double outer_tol= _res/_red;
            
            _solver.SetTol( outer_tol);
            _solver.Solve( _AN, B, w, q, d, e);
            
            // update solution:
            v.Data-= w;
            p     -= q;
        }
}


}    // end of namespace DROPS

#endif
