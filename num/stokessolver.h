//**************************************************************************
// File:    stokessolver.h                                                 *
// Content: solvers for the Stokes problem                                 *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 20 2001                                           *
//**************************************************************************

#ifndef _STOKESSOLVER_H_
#define _STOKESSOLVER_H_

#include "num/solver.h"

namespace DROPS
{

template <class PoissonSolverT>
class UzawaSolverCL
{
  private:
    PoissonSolverT& _poissonSolver;
    MatrixCL&       _M; 
    
    double _tau, _tol, _res;
    int    _maxiter, _iter;
    
  public:
    UzawaSolverCL( PoissonSolverT& solver, MatrixCL& M, double tol, int maxiter, double tau= 1.)
        : _poissonSolver( solver), _M( M), _tau( tau),
          _tol( tol), _res( -1.), _maxiter( maxiter), _iter( -1)    {}

    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    
    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c);
    // solves the system   A v + BT p = b
    //                     B v        = c
};

template <class PoissonSolverT>
class SchurSolverCL
{
  private:
    PoissonSolverT&  _poissonSolver;
    
    double _tol, _res;
    int    _maxiter, _iter;
    
  public:
    SchurSolverCL( PoissonSolverT& solver, double tol, int maxiter)
        : _poissonSolver( solver),
          _tol( tol), _res( -1.), _maxiter( maxiter), _iter( -1)    {}
    
    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    
    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c);
    // solves the system   A v + BT p = b
    //                     B v        = c
};    
    
template <class PoissonSolverT>
class PSchurSolverCL
{
  private:
    typedef SsorMassPcCL<SchurComplMatrixCL> SchurPcT;
    SchurPcT _schurPc;
    PoissonSolverT&  _poissonSolver;
    
    double _tol, _res;
    int    _maxiter, _iter;
    
  public:
    PSchurSolverCL( PoissonSolverT& solver, MatrixCL& M, double tol, int maxiter)
        : _schurPc( M), _poissonSolver( solver),
          _tol( tol), _res( -1.), _maxiter( maxiter), _iter( -1)    {}
    
    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    
    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c);
    // solves the system   A v + BT p = b
    //                     B v        = c
};    


//==================================================
//        derived classes for easier use
//==================================================

typedef PCGSolverCL<VectorCL, double, SsorPcCL<VectorCL, double> > 
        PCG_T;

typedef PCGSolverCL<VectorCL, double, ImprovedSsorPcCL<VectorCL, double> > 
        IPCG_T;

typedef PCGSolverCL<VectorCL, double, SGSPcCL<VectorCL, double> > 
        GSPCG_T;

class Uzawa_PCG_CL: public UzawaSolverCL<PCG_T>
{
  private:
    PCG_T _PCGsolver;
  public:
    Uzawa_PCG_CL( MatrixCL& M, double outer_tol, int outer_iter, double inner_tol, int inner_iter, double tau= 1.)
        : UzawaSolverCL<PCG_T>( _PCGsolver,  M, outer_tol, outer_iter, tau), 
          _PCGsolver( inner_tol, inner_iter, SsorPcCL<VectorCL, double>( 1.))
        {}
};    


class Uzawa_IPCG_CL
{
  private:
    IPCG_T _M_IPCGsolver;
    IPCG_T _A_IPCGsolver;

    MatrixCL&       _M; 
    
    double _tau, _tol, _res;
    int    _maxiter, _iter;

  public:
    Uzawa_IPCG_CL(MatrixCL& M, double outer_tol, int outer_iter, double inner_tol, int inner_iter, double tau= 1.)
        :  
          _M_IPCGsolver( inner_tol, inner_iter, ImprovedSsorPcCL<VectorCL, double>( 1.) ),
          _A_IPCGsolver( inner_tol, inner_iter, ImprovedSsorPcCL<VectorCL, double>( 1.) ),
          _M(M), _tau(tau), _tol(outer_tol), _res(-1.), _maxiter(outer_iter), _iter(-1)
        { _M_IPCGsolver.GetPc().Init(_M); }

    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    
    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }

    // Always call this when A has changed, before Solve()!
    void Init_A_Pc(MatrixCL& A) { _A_IPCGsolver.GetPc().Init(A); }

    inline void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c);
};    
    
    
class Schur_PCG_CL: public SchurSolverCL<PCG_T>
{
  private:
    PCG_T _PCGsolver;
  public:
    Schur_PCG_CL(double outer_tol, int outer_iter, double inner_tol, int inner_iter)
        : SchurSolverCL<PCG_T>( _PCGsolver, outer_tol, outer_iter),
          _PCGsolver( inner_tol, inner_iter, SsorPcCL<VectorCL, double>( 1.))
        {}
};    
    
class PSchur_PCG_CL: public PSchurSolverCL<PCG_T>
{
  private:
    PCG_T _PCGsolver;
  public:
    PSchur_PCG_CL( MatrixCL& M, double outer_tol, int outer_iter, double inner_tol, int inner_iter)
        : PSchurSolverCL<PCG_T>( _PCGsolver, M, outer_tol, outer_iter),
          _PCGsolver( inner_tol, inner_iter, SsorPcCL<VectorCL, double>( 1.))
        {}
};    
    
class PSchur_IPCG_CL: public PSchurSolverCL<IPCG_T>
{
  private:
    IPCG_T _PCGsolver;
  public:
    PSchur_IPCG_CL( MatrixCL& M, double outer_tol, int outer_iter, double inner_tol, int inner_iter)
        : PSchurSolverCL<IPCG_T>( _PCGsolver, M, outer_tol, outer_iter),
          _PCGsolver( inner_tol, inner_iter, ImprovedSsorPcCL<VectorCL, double>( 1.))
        {}
    IPCG_T& GetPoissonSolver() { return _PCGsolver; }
};    
    
class PSchur_GSPCG_CL: public PSchurSolverCL<GSPCG_T>
{
  private:
    GSPCG_T _PCGsolver;
  public:
    PSchur_GSPCG_CL( MatrixCL& M, double outer_tol, int outer_iter, double inner_tol, int inner_iter)
        : PSchurSolverCL<GSPCG_T>( _PCGsolver, M, outer_tol, outer_iter),
          _PCGsolver( inner_tol, inner_iter, SGSPcCL<VectorCL, double>() )
        {}
    GSPCG_T& GetPoissonSolver() { return _PCGsolver; }
};    
    
// TODO: (P)Schur_MG_CL    



//=================================
//     template definitions
//=================================

template <class PoissonSolverT>
void UzawaSolverCL<PoissonSolverT>::Solve
    ( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c)

{
    VectorCL v_corr(v.size()),
             p_corr(p.size()),
             res1(v.size()),
             res2(p.size());

    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter)
    {
//        _poissonSolver.Solve( _M, p_corr, res2= B*v - c);
        z_xpay(res2, B*v, -1.0, c);
        _poissonSolver.Solve(_M, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(_tau, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2(res1, A*v, 1.0, transp_mul(B,p), -1.0, b);
        res1_norm= res1.norm2();
        res2_norm= res2.norm2();

        if (res1_norm + res2_norm < tol)
        {
            _res= ::sqrt( res1_norm + res2_norm );
            return;
        }   

        if( (_iter%output)==0 )
            std::cerr << "step " << _iter << ": norm of 1st eq= " << ::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << ::sqrt( res2_norm) << std::endl;

        _poissonSolver.Solve( A, v_corr, res1);
        v-= v_corr;
    }
    _res= ::sqrt( res1_norm + res2_norm );
}


template <class PoissonSolverT>
void SchurSolverCL<PoissonSolverT>::Solve
    ( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c)
{
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p

    VectorCL rhs= -c;
    {
        VectorCL tmp( v.size());
        _poissonSolver.Solve( A, tmp, b);
        std::cerr << "Iterationen: " << _poissonSolver.GetIter() << "    Norm des Residuums: " << _poissonSolver.GetResid() << std::endl;
        rhs+= B*tmp;
    }
    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;   
    double tol= _tol;     
    CG( SchurComplMatrixCL( A, B, _poissonSolver.GetTol(), 1.), p, rhs, iter, tol);
    std::cerr << "Iterationen: " << iter << "    Norm des Residuums: " << tol << std::endl;
    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;

    tol= _poissonSolver.GetTol();
    _poissonSolver.SetTol( _tol);  // same tolerance as for pressure
    _poissonSolver.Solve( A, v, b - transp_mul(B, p));
    _poissonSolver.SetTol( tol);   // reset old tolerance, so that nothing has changed
    std::cerr << "Iterationen: " << _poissonSolver.GetIter() << "    Norm des Residuums: " << _poissonSolver.GetResid() << std::endl;
    std::cerr << "-----------------------------------------------------" << std::endl;
}

inline void Uzawa_IPCG_CL::Solve
    ( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c)

{
    VectorCL v_corr(v.size()),
             p_corr(p.size()),
             res1(v.size()),
             res2(p.size());

    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter)
    {
//        _poissonSolver.Solve( _M, p_corr, res2= B*v - c);
        z_xpay(res2, B*v, -1.0, c);
        _M_IPCGsolver.Solve(_M, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(_tau, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2(res1, A*v, 1.0, transp_mul(B,p), -1.0, b);
        res1_norm= res1.norm2();
        res2_norm= res2.norm2();

        if (res1_norm + res2_norm < tol)
        {
            _res= ::sqrt( res1_norm + res2_norm );
            return;
        }   

        if( (_iter%output)==0 )
            std::cerr << "step " << _iter << ": norm of 1st eq= " << ::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << ::sqrt( res2_norm) << std::endl;

        _A_IPCGsolver.Solve( A, v_corr, res1);
//        v-= v_corr;
        axpy(-1.0, v_corr, v);
    }
    _res= ::sqrt( res1_norm + res2_norm );
}

template <class PoissonSolverT>
void PSchurSolverCL<PoissonSolverT>::Solve
    ( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, VectorCL& b, VectorCL& c)
{
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p

    VectorCL rhs= -c;
    {
        VectorCL tmp( v.size());
        _poissonSolver.Solve( A, tmp, b);
        std::cerr << "Iterationen: " << _poissonSolver.GetIter() << "    Norm des Residuums: " << _poissonSolver.GetResid() << std::endl;
        rhs+= B*tmp;
    }
    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;   
    double tol= _tol;     
    PCG( SchurComplMatrixCL( A, B, _poissonSolver.GetTol(), 1.), p, rhs, _schurPc, iter, tol);
    std::cerr << "Iterationen: " << iter << "    Norm des Residuums: " << tol << std::endl;
    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;

    tol= _poissonSolver.GetTol();
    _poissonSolver.SetTol( _tol);  // same tolerance as for pressure
    _poissonSolver.Solve( A, v, b - transp_mul(B, p));
    _poissonSolver.SetTol( tol);   // reset old tolerance, so that nothing has changed
    std::cerr << "Iterationen: " << _poissonSolver.GetIter() << "    Norm des Residuums: " << _poissonSolver.GetResid() << std::endl;
    std::cerr << "-----------------------------------------------------" << std::endl;
}


}    // end of namespace DROPS

#endif
