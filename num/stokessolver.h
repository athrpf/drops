//**************************************************************************
// File:    stokessolver.h                                                 *
// Content: solvers for the Stokes problem                                 *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 20 2001                                           *
//**************************************************************************

#ifndef DROPS_STOKESSOLVER_H
#define DROPS_STOKESSOLVER_H

#include "num/solver.h"
#include "num/MGsolver.h"

namespace DROPS
{

//=============================================================================
//  The Stokes solvers solve systems of the form
//    A v + BT p = b
//    B v        = c
//=============================================================================

template <typename PoissonSolverT>
class SchurSolverCL : public SolverBaseCL
{
  private:
    PoissonSolverT& _poissonSolver;

  public:
    SchurSolverCL (PoissonSolverT& solver, int maxiter, double tol)
        : SolverBaseCL(maxiter,tol), _poissonSolver(solver) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};


template <typename PoissonSolverT>
class PSchurSolverCL : public SolverBaseCL
{
  private:
    PoissonSolverT&        _poissonSolver;
    PreGSOwnMatCL<P_SSOR0> _schurPc;

  public:
    PSchurSolverCL (PoissonSolverT& solver, MatrixCL& M, int maxiter, double tol)
        : SolverBaseCL(maxiter,tol), _poissonSolver(solver), _schurPc(M) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};


template <typename PoissonSolverT>
class UzawaSolverCL : public SolverBaseCL
{
  private:
    PoissonSolverT& _poissonSolver;
    MatrixCL&       _M;
    double          _tau;

  public:
    UzawaSolverCL (PoissonSolverT& solver, MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : SolverBaseCL(maxiter,tol), _poissonSolver(solver), _M(M), _tau(tau) {}

    double GetTau()            const { return _tau; }
    void   SetTau( double tau)       { _tau= tau; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <typename PoissonSolverT, typename PoissonSolver2T>
class UzawaSolver2CL : public SolverBaseCL
{
  private:
    PoissonSolverT&  poissonSolver_;
    PoissonSolver2T& poissonSolver2_;
    MatrixCL& M_;
    double    tau_;

  public:
    UzawaSolver2CL (PoissonSolverT& solver, PoissonSolver2T& solver2,
                    MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : SolverBaseCL( maxiter, tol), poissonSolver_( solver), poissonSolver2_( solver2),
          M_( M), tau_( tau) {}

    double GetTau()            const { return tau_; }
    void   SetTau( double tau)       { tau_= tau; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

class Uzawa_IPCG_CL : public SolverBaseCL
{
  private:
    PCG_SsorDiagCL _M_IPCGsolver;
    PCG_SsorDiagCL _A_IPCGsolver;
    MatrixCL&      _M;
    double         _tau;

  public:
    Uzawa_IPCG_CL(MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1.)
        : SolverBaseCL(outer_iter,outer_tol),
          _M_IPCGsolver( SSORDiagPcCL(1.), inner_iter, inner_tol ),
          _A_IPCGsolver( SSORDiagPcCL(1.), inner_iter, inner_tol ),
          _M(M), _tau(tau)
        { _M_IPCGsolver.GetPc().Init(_M); }

    // Always call this when A has changed, before Solve()!
    void Init_A_Pc(MatrixCL& A) { _A_IPCGsolver.GetPc().Init(A); }

    inline void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                       const VectorCL& b, const VectorCL& c);
};


//=============================================================================
//  Derived classes for easier use
//=============================================================================

class Uzawa_CG_CL : public UzawaSolverCL<CGSolverCL>
{
  private:
    CGSolverCL _CGsolver;
  public:
    // XXX M is not used!! Perhaps modify UzawaSolverCL. Or just live with it.
    Uzawa_CG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1.)
        : UzawaSolverCL<CGSolverCL>( _CGsolver, M, outer_iter, outer_tol, tau),
          _CGsolver( inner_iter, inner_tol)
        {}
};

class Uzawa_PCG_CL : public UzawaSolverCL<PCG_SsorCL>
{
  private:
    PCG_SsorCL _PCGsolver;
  public:
    Uzawa_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1.)
        : UzawaSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol, tau),
          _PCGsolver(SSORPcCL(1.), inner_iter, inner_tol)
        {}
};

class Uzawa_SGSPCG_CL : public UzawaSolverCL<PCG_SgsCL>
{
  private:
    PCG_SgsCL _PCGsolver;
  public:
    Uzawa_SGSPCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1.)
        : UzawaSolverCL<PCG_SgsCL>( _PCGsolver, M, outer_iter, outer_tol, tau),
          _PCGsolver(SGSPcCL(1.), inner_iter, inner_tol)
        {}
};

class Schur_PCG_CL: public SchurSolverCL<PCG_SsorCL>
{
  private:
    PCG_SsorCL _PCGsolver;
  public:
    Schur_PCG_CL(int outer_iter, double outer_tol, int inner_iter, double inner_tol)
        : SchurSolverCL<PCG_SsorCL>( _PCGsolver, outer_iter, outer_tol),
          _PCGsolver(SSORPcCL(1.), inner_iter, inner_tol)
        {}
};


class PSchur_PCG_CL: public PSchurSolverCL<PCG_SsorCL>
{
  private:
    PCG_SsorCL _PCGsolver;
  public:
    PSchur_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol)
        : PSchurSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol),
          _PCGsolver(SSORPcCL(1.), inner_iter, inner_tol)
        {}
};


class PSchur_IPCG_CL: public PSchurSolverCL<PCG_SsorDiagCL>
{
  private:
    PCG_SsorDiagCL _PCGsolver;
  public:
    PSchur_IPCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol)
        : PSchurSolverCL<PCG_SsorDiagCL>( _PCGsolver, M, outer_iter, outer_tol),
          _PCGsolver(SSORDiagPcCL(1.), inner_iter, inner_tol)
        {}
    PCG_SsorDiagCL& GetPoissonSolver() { return _PCGsolver; }
};


class PSchur_GSPCG_CL: public PSchurSolverCL<PCG_SgsCL>
{
  private:
    PCG_SgsCL _PCGsolver;
  public:
    PSchur_GSPCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol)
        : PSchurSolverCL<PCG_SgsCL>( _PCGsolver, M, outer_iter, outer_tol),
          _PCGsolver(SGSPcCL(), inner_iter, inner_tol)
        {}
    PCG_SgsCL& GetPoissonSolver() { return _PCGsolver; }
};


class PSchur_MG_CL: public PSchurSolverCL<MGSolverCL>
{
  private:
    MGSolverCL _MGsolver;
  public:
    PSchur_MG_CL( MatrixCL& M,      int outer_iter, double outer_tol, 
                  MGDataCL& MGData, int inner_iter, double inner_tol )
        : PSchurSolverCL<MGSolverCL>( _MGsolver, M, outer_iter, outer_tol ),
          _MGsolver( MGData, inner_iter, inner_tol )
        {}
};

class Uzawa_MG_CL : public UzawaSolver2CL<PCG_SsorCL, MGSolverCL>
{
  private:
    PCG_SsorCL PCGsolver_;
    MGSolverCL MGsolver_;

  public:
    Uzawa_MG_CL(MatrixCL& M,      int outer_iter, double outer_tol,
                MGDataCL& MGData, int inner_iter, double inner_tol, double tau= 1.)
        : UzawaSolver2CL<PCG_SsorCL, MGSolverCL>( PCGsolver_, MGsolver_, M,
                                                  outer_iter, outer_tol, tau),
          PCGsolver_( SSORPcCL( 1.), inner_iter, inner_tol),
          MGsolver_( MGData, inner_iter, inner_tol)
        {}
};


//=============================================================================
//  SchurComplMatrixCL
//=============================================================================

template<typename>
class SchurComplMatrixCL;

template<typename T>
VectorCL operator*(const SchurComplMatrixCL<T>&, const VectorCL&);


template<class PoissonSolverT>
class SchurComplMatrixCL
{
  private:
    PoissonSolverT& solver_;
    const MatrixCL& A_;
    const MatrixCL& B_;

  public:
    SchurComplMatrixCL(PoissonSolverT& solver, const MatrixCL& A, const MatrixCL& B)
        : solver_( solver), A_( A), B_( B) {}

    friend VectorCL
    operator*<>(const SchurComplMatrixCL<PoissonSolverT>&, const VectorCL&);
};


template<class PoissonSolverT>
VectorCL operator*(const SchurComplMatrixCL<PoissonSolverT>& M, const VectorCL& v)
{
    VectorCL x( M.A_.num_cols());

    M.solver_.Solve( M.A_, x, transp_mul( M.B_, v));
//    std::cerr << "iterations: " << M.solver_.GetIter()
//              << "\tresidual: " << M.solver_.GetResid() << std::endl;
    return M.B_*x;
}


//=============================================================================
//  The "Solve" functions
//=============================================================================

template <class PoissonSolverT>
void UzawaSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    VectorCL v_corr(v.size()),
             p_corr(p.size()),
             res1(v.size()),
             res2(p.size());

    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= res2.norm2();
        _poissonSolver.SetTol( std::sqrt( res2_norm)/20.0);
        _poissonSolver.Solve(_M, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(_tau, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2(res1, A*v, 1.0, transp_mul(B,p), -1.0, b);
        res1_norm= res1.norm2();
        if (res1_norm + res2_norm < tol) {
            _res= ::sqrt( res1_norm + res2_norm );
            return;
        }
        if( (_iter%output)==0)
            std::cerr << "step " << _iter << ": norm of 1st eq= " << ::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << ::sqrt( res2_norm) << std::endl;

        _poissonSolver.SetTol( std::sqrt( res1_norm)/20.0);
        _poissonSolver.Solve( A, v_corr, res1);
        v-= v_corr;
    }
    _res= ::sqrt( res1_norm + res2_norm );
}

template <class PoissonSolverT, class PoissonSolver2T>
void UzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::Solve(
    const MatrixCL& A, const MatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    VectorCL v_corr( v.size()),
             p_corr( p.size()),
             res1( v.size()),
             res2( p.size());
    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= res2.norm2();
        poissonSolver_.SetTol( std::sqrt( res2_norm)/20.0);
        poissonSolver_.Solve( M_, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(tau_, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2( res1, A*v, 1.0, transp_mul( B, p), -1.0, b);
        res1_norm= res1.norm2();
        if (res1_norm + res2_norm < tol) {
            _res= ::sqrt( res1_norm + res2_norm);
            return;
        }
        if( (_iter%output)==0)
            std::cerr << "step " << _iter << ": norm of 1st eq= " << ::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << ::sqrt( res2_norm) << std::endl;

        poissonSolver2_.SetTol( std::sqrt( res1_norm)/20.0);
        poissonSolver2_.Solve( A, v_corr, res1);
        v-= v_corr;
    }
    _res= ::sqrt( res1_norm + res2_norm );
}

template <class PoissonSolverT>
void SchurSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    VectorCL rhs= -c;
    {
        VectorCL tmp( v.size());
        _poissonSolver.Solve( A, tmp, b);
        std::cerr << "iterations: " << _poissonSolver.GetIter()
                  << "\tresidual: " << _poissonSolver.GetResid() << std::endl;
        rhs+= B*tmp;
    }
    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;
    double tol= _tol;
    CG( SchurComplMatrixCL<PoissonSolverT>( _poissonSolver, A, B), p, rhs, iter, tol);
    std::cerr << "iterations: " << iter << "\tresidual: " << tol << std::endl;
    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;

    const double old_tol= _poissonSolver.GetTol();
    _poissonSolver.SetTol( _tol);      // same tolerance as for pressure
    _poissonSolver.Solve( A, v, b - transp_mul(B, p));
    _poissonSolver.SetTol( old_tol);   // reset old tolerance, so that nothing has changed
    std::cerr << "Iterationen: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;

    _iter= iter+_poissonSolver.GetIter();
    _res= tol + _poissonSolver.GetResid();
std::cerr << "Real residuals are: "
          << (A*v+transp_mul(B, p)-b).norm() << ", "
          << (B*v-c).norm() << std::endl;
    std::cerr << "-----------------------------------------------------" << std::endl;
}

inline void Uzawa_IPCG_CL::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)

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
void PSchurSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    VectorCL rhs= -c;
    {
        VectorCL tmp( v.size());
        _poissonSolver.Solve( A, tmp, b);
        std::cerr << "iterations: " << _poissonSolver.GetIter()
                  << "\tresidual: " << _poissonSolver.GetResid() << std::endl;
        rhs+= B*tmp;
    }
    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;
    double tol= _tol;
    PCG( SchurComplMatrixCL<PoissonSolverT>( _poissonSolver, A, B), p, rhs, _schurPc, iter, tol);
    std::cerr << "iterations: " << iter << "\tresidual: " << tol << std::endl;
    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;

    const double old_tol= _poissonSolver.GetTol();
    _poissonSolver.SetTol( _tol);      // same tolerance as for pressure
    _poissonSolver.Solve( A, v, b - transp_mul(B, p));
    _poissonSolver.SetTol( old_tol);   // reset old tolerance, so that nothing has changed
    std::cerr << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;

    _iter= iter+_poissonSolver.GetIter();
    _res= std::sqrt( tol*tol + _poissonSolver.GetResid()*_poissonSolver.GetResid());
    std::cerr << "-----------------------------------------------------" << std::endl;
}


} // end of namespace DROPS

#endif
