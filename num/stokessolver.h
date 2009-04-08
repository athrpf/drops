/// \file stokessolver.h
/// \brief Solvers for the Stokes problem
/// \author Volker Reichelt, Sven Gross, Joerg Grande, Patrick Esser, IGPM

#ifndef DROPS_STOKESSOLVER_H
#define DROPS_STOKESSOLVER_H

#include "num/solver.h"
#include "num/MGsolver.h"
#include "stokes/integrTime.h"
#ifdef _PAR
#  include "parallel/exchange.h"
#endif
namespace DROPS
{

//=============================================================================
//  The Stokes solvers solve systems of the form
//    A v + BT p = b
//    B v        = c
//=============================================================================

// What every iterative stokes-solver should have
class StokesSolverBaseCL: public SolverBaseCL
{
  public:
    StokesSolverBaseCL (int maxiter, double tol, bool rel= false, std::ostream* output= 0)
        : SolverBaseCL(maxiter, tol, rel, output){}

    virtual void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c) = 0;
    virtual void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c) = 0;
};

template <typename PoissonSolverT>
class SchurSolverCL : public StokesSolverBaseCL
{
  private:
    PoissonSolverT& _poissonSolver;
    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    SchurSolverCL (PoissonSolverT& solver, int maxiter, double tol)
        : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <typename PoissonSolverT>
class PSchurSolverCL : public StokesSolverBaseCL
{
  private:
    PoissonSolverT&        _poissonSolver;
    PreGSOwnMatCL<P_SSOR0> _schurPc;

    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    PSchurSolverCL (PoissonSolverT& solver, MatrixCL& M, int maxiter, double tol)
        : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver), _schurPc(M) {}
    PSchurSolverCL (PoissonSolverT& solver, MLMatrixCL& M, int maxiter, double tol)
    : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver), _schurPc( M.GetFinest()) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <typename InnerSolverT, typename OuterSolverT>
class PSchurSolver2CL : public StokesSolverBaseCL
{
  private:
    InnerSolverT& innerSolver_;
    OuterSolverT& outerSolver_;

    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    PSchurSolver2CL( InnerSolverT& solver1, OuterSolverT& solver2,
                    int maxiter, double tol)
        : StokesSolverBaseCL( maxiter, tol), innerSolver_( solver1),
          outerSolver_( solver2) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);

    InnerSolverT& GetInnerSolver() { return innerSolver_; }
};

template <typename PoissonSolverT>
class UzawaSolverCL : public StokesSolverBaseCL
{
  private:
    PoissonSolverT& _poissonSolver;
    MatrixCL&       _M;
    double          _tau;

    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    UzawaSolverCL (PoissonSolverT& solver, MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver), _M(M), _tau(tau) {}
    UzawaSolverCL (PoissonSolverT& solver, MLMatrixCL& M, int maxiter, double tol, double tau= 1.)
    : StokesSolverBaseCL(maxiter,tol), _poissonSolver(solver), _M(M.GetFinest()), _tau(tau) {}


    double GetTau()            const { return _tau; }
    void   SetTau( double tau)       { _tau= tau; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <typename PoissonSolverT, typename PoissonSolver2T>
class UzawaSolver2CL : public StokesSolverBaseCL
{
  private:
    PoissonSolverT&  poissonSolver_;
    PoissonSolver2T& poissonSolver2_;
    MatrixCL& M_;
    double    tau_;

    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    UzawaSolver2CL (PoissonSolverT& solver, PoissonSolver2T& solver2,
                    MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : StokesSolverBaseCL( maxiter, tol), poissonSolver_( solver), poissonSolver2_( solver2),
          M_( M), tau_( tau) {}

    double GetTau()            const { return tau_; }
    void   SetTau( double tau)       { tau_= tau; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

class Uzawa_IPCG_CL : public StokesSolverBaseCL
{
  private:
    SSORDiagPcCL   _ssordiagpc;
    PCG_SsorDiagCL _M_IPCGsolver;
    PCG_SsorDiagCL _A_IPCGsolver;
    MatrixCL&      _M;
    double         _tau;

    template <typename Mat, typename Vec>
    void doSolve( const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c);

  public:
    Uzawa_IPCG_CL(MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1.)
        : StokesSolverBaseCL(outer_iter,outer_tol),
          _M_IPCGsolver( _ssordiagpc, inner_iter, inner_tol ),
          _A_IPCGsolver( _ssordiagpc, inner_iter, inner_tol ),
          _M(M), _tau(tau)
        { _M_IPCGsolver.GetPc().Init(_M); }

    // Always call this when A has changed, before Solve()!
    void Init_A_Pc(MatrixCL& A) { _A_IPCGsolver.GetPc().Init(A); }

    inline void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                       const VectorCL& b, const VectorCL& c);
    inline void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};


//=============================================================================
// Ready-to-use solver-class with inexact Uzawa method InexactUzawa.
//=============================================================================

// Characteristics of the preconditioner for the A-block
enum InexactUzawaApcMethodT { APC_OTHER, APC_SYM, APC_SYM_LINEAR };

//=============================================================================
// Inexact Uzawa-method from "Fast Iterative Solvers for Discrete Stokes
// Equations", Peters, Reichelt, Reusken, Chapter 3.3.
// The preconditioner Apc for A must be "good" (MG-like) to guarantee
// convergence.
//=============================================================================
template <class ApcT, class SpcT, InexactUzawaApcMethodT ApcMeth= APC_OTHER>
  class InexactUzawaCL: public StokesSolverBaseCL
{
  private:
    ApcT& Apc_;
    SpcT& Spc_;
    double innerreduction_;
    int    innermaxiter_;

  public:
    InexactUzawaCL(ApcT& Apc, SpcT& Spc, int outer_iter, double outer_tol,
        double innerreduction= 0.3, int innermaxiter= 500)
        :StokesSolverBaseCL( outer_iter, outer_tol),
         Apc_( Apc), Spc_( Spc), innerreduction_( innerreduction), innermaxiter_( innermaxiter)
    {}
    inline void
    Solve(const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
          const VectorCL& b, const VectorCL& c);
    inline
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};


// Use a Krylow-method (from num/solver.h) with the standard-interface of
// the Stokessolvers in this file.
template <class SolverT>
class BlockMatrixSolverCL: public StokesSolverBaseCL
{
  private:
    SolverT& solver_;
#ifdef _PAR
    ExchangeBlockCL exBlock_;
#endif

  public:
#ifndef _PAR
    BlockMatrixSolverCL( SolverT& solver)
        : StokesSolverBaseCL(-1, -1.0), solver_( solver) {}
#else
    BlockMatrixSolverCL( SolverT& solver, const IdxDescCL& vel_idx, const IdxDescCL& pr_idx)
        : StokesSolverBaseCL(-1, -1.0), solver_( solver)
    {
        exBlock_.AttachTo( vel_idx);
        exBlock_.AttachTo( pr_idx);
    }
#endif

// We overwrite these functions.
    void   SetTol     (double tol) { solver_.SetTol( tol); }
    void   SetMaxIter (int iter)   { solver_.SetMaxIter( iter); }
    void   SetRelError(bool rel)   { solver_.SetRelError( rel); }

    double GetTol     () const { return solver_.GetTol(); }
    int    GetMaxIter () const { return solver_.GetMaxIter(); }
    double GetResid   () const { return solver_.GetResid(); }
    int    GetIter    () const { return solver_.GetIter(); }
    bool   GetRelError() const { return solver_.GetRelError(); }

    void
    Solve(const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
          const VectorCL& b, const VectorCL& c) {
        BlockMatrixCL M( &A, MUL, &B, TRANSP_MUL, &B, MUL);
        VectorCL rhs( M.num_rows());
        rhs[std::slice( 0, M.num_rows( 0), 1)]= b;
        rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]= c;
        VectorCL x( M.num_cols());
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
#ifndef _PAR
        solver_.Solve( M, x, rhs);
#else
        solver_.Solve( M, x, rhs, exBlock_);
#endif
        v= x[std::slice( 0, M.num_cols( 0), 1)];
        p= x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)];
    }
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        MLBlockMatrixCL M( &A, MUL, &B, TRANSP_MUL, &B, MUL);
        VectorCL rhs( M.num_rows());
        rhs[std::slice( 0, M.num_rows( 0), 1)]= b;
        rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]= c;
        VectorCL x( M.num_cols());
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
#ifndef _PAR
        solver_.Solve( M, x, rhs);
#else
        solver_.Solve( M, x, rhs, exBlock_);
#endif
        v= x[std::slice( 0, M.num_cols( 0), 1)];
        p= x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)];
    }
};

// Upper block-triangular preconditioning strategy in BlockPreCL
struct UpperBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (PC1T& pc1, PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc2.Apply( /*dummy*/ B, p, c);
        Vec b2( b);
        b2-= transp_mul( B, p);
        pc1.Apply( A, v, b2);
    }
};

// Block-diagonal preconditioning strategy in BlockPreCL
struct DiagBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (PC1T& pc1, PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc1.Apply( A, v, b);
        pc2.Apply( /*dummy*/ B, p, c);
   }
};

// Lower block-triangular preconditioning strategy in BlockPreCL
struct LowerBlockPreCL
{
    template <class PC1T, class PC2T, class Mat, class Vec>
    static void
    Apply (PC1T& pc1, PC2T& pc2, const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) {
        pc1.Apply( A, v, b);
        Vec c2( c);
#ifdef _PAR
        Assert(pc1.RetAcc(), DROPSErrCL("LowerBlockPreCL::Apply: Accumulation is missing"), DebugParallelNumC);
#endif
        c2-= B*v;
        pc2.Apply( /*dummy*/ B, p, c2);
   }
};

// With BlockShapeT= DiagBlockPreCL, this is the diagonal PC ( pc1^(-1) 0 \\ 0 pc2^(-1) ),
// else it is block-triangular:
// Upper... ( pc1^(-1) B^T \\ 0 pc2^(-1) ), resp. lower: ( pc1^(-1) 0 \\ B pc2^(-1) )
template <class PC1T, class PC2T, class BlockShapeT= DiagBlockPreCL>
class BlockPreCL
{
  private:
    PC1T& pc1_; // Preconditioner for A.
    PC2T& pc2_; // Preconditioner for S.

  public:
    BlockPreCL (PC1T& pc1, PC2T& pc2)
        : pc1_( pc1), pc2_( pc2) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
        BlockShapeT::Apply( pc1_, pc2_, A, B, v, p, b, c);
    }

    template <typename Mat, typename Vec>
    void
    Apply(const BlockMatrixBaseCL<Mat>& A, Vec& x, const Vec& b) const {
        VectorCL b0( b[std::slice( 0, A.num_rows( 0), 1)]);
        VectorCL b1( b[std::slice( A.num_rows( 0), A.num_rows( 1), 1)]);
        VectorCL x0( A.num_cols( 0));
        VectorCL x1( A.num_cols( 1));
        BlockShapeT::Apply( pc1_, pc2_, *A.GetBlock( 0), *A.GetBlock( 2), x0, x1, b0, b1);
        x[std::slice( 0, A.num_cols( 0), 1)]= x0;
        x[std::slice( A.num_cols( 0), A.num_cols( 1), 1)]= x1;
    }

#ifdef _PAR
    /// \brief Check if the preconditioned vector is accumulated
    bool RetAcc() const {
        Assert( pc1_.RetAcc()==pc2_.RetAcc(), DROPSErrCL("BlockPreCL::RetAcc: Preconditioners do not match"),
                DebugNumericC);
        return pc1_.RetAcc();
    }

    /// \brief Check if the diagonal of the matrix needs to be computed
    bool NeedDiag() const { return pc1_.NeedDiag() || pc2_.NeedDiag(); }

    /// \brief Set accumulated diagonal of a matrix, that is needed by most of the preconditioners
    template<typename Mat>
    void SetDiag(const Mat& A)
    {
        pc1_.SetDiag(*A.GetBlock( 0));
        pc2_.SetDiag(/*dummy*/ *(A.GetBlock( 3)!=0 ? A.GetBlock( 3) : A.GetBlock( 1)));
    }
    const PC1T& GetPC1() const { return pc1_; }
          PC1T& GetPC1()       { return pc1_; }
    const PC2T& GetPC2() const { return pc2_; }
          PC2T& GetPC2()       { return pc2_; }

#endif
};

// StokesSolverBaseCL can be used as a preconditioner for the methods in solver.h
class StokesSolverAsPreCL
{
  private:
    StokesSolverBaseCL& solver_;
    mutable std::ostream* output_;

  public:
    StokesSolverAsPreCL( StokesSolverBaseCL& solver, int max_iter, std::ostream* output= 0)
        : solver_( solver), output_( output) {solver.SetMaxIter(max_iter);}

    void
    Apply(const BlockMatrixCL& M, VectorCL& x, const VectorCL& rhs) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)),
                 b(rhs[std::slice( 0, M.num_rows( 0), 1)]),
                 c(rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]);
        solver_.Solve(*M.GetBlock(0), *M.GetBlock(1), v, p, b, c);
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
        if (output_ != 0)
            *output_<< "StokesSolverAsPreCL: iterations: " << solver_.GetIter()
                    << "\trelative residual: " << solver_.GetResid() << std::endl;
    }
    void
    Apply(const MLBlockMatrixCL& M, VectorCL& x, const VectorCL& rhs) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)),
                 b(rhs[std::slice( 0, M.num_rows( 0), 1)]),
                 c(rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]);
        solver_.Solve(*M.GetBlock(0), *M.GetBlock(1), v, p, b, c);
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
        if (output_ != 0)
            *output_<< "StokesSolverAsPreCL: iterations: " << solver_.GetIter()
                    << "\trelative residual: " << solver_.GetResid() << std::endl;
    }
};

//=============================================================================
//  SchurComplMatrixCL
//=============================================================================

template<typename, typename>
class SchurComplMatrixCL;

template<typename T, typename Mat>
VectorCL operator*(const SchurComplMatrixCL<T, Mat>&, const VectorCL&);


template<class PoissonSolverT, typename Mat>
class SchurComplMatrixCL
{
  private:
    PoissonSolverT& solver_;
    const Mat& A_;
    const Mat& B_;

  public:
    SchurComplMatrixCL(PoissonSolverT& solver, const Mat& A, const Mat& B)
        : solver_( solver), A_( A), B_( B) {}

    friend VectorCL
    operator*<>(const SchurComplMatrixCL<PoissonSolverT, Mat>&, const VectorCL&);
};


template<class PoissonSolverT, class Mat>
VectorCL operator*(const SchurComplMatrixCL<PoissonSolverT, Mat>& M, const VectorCL& v)
{
    VectorCL x( M.A_.num_cols());
    M.solver_.Solve( M.A_, x, transp_mul( M.B_, v));
//    std::cout << "> inner iterations: " << M.solver_.GetIter()
//              << "\tresidual: " << M.solver_.GetResid() << std::endl;
    return M.B_*x;
}


//=============================================================================
// ApproximateSchurComplMatrixCL
// BApc^{-1}B^T, where Apc is a preconditioner for A.
//=============================================================================
template<typename, typename>
class ApproximateSchurComplMatrixCL;

template<typename T, typename MatT>
VectorCL operator*(const ApproximateSchurComplMatrixCL<T, MatT>&, const VectorCL&);

template<typename APC, typename Mat>
class ApproximateSchurComplMatrixCL
{
  private:
    const Mat& A_;
    APC& Apc_;
    const Mat& B_;

  public:
    ApproximateSchurComplMatrixCL(const Mat& A, APC& Apc, const Mat& B)
        : A_( A), Apc_( Apc), B_( B) {}

    friend VectorCL
    operator*<>(const ApproximateSchurComplMatrixCL<APC, Mat>&, const VectorCL&);
};

template<typename APC, typename Mat>
VectorCL operator*(const ApproximateSchurComplMatrixCL<APC, Mat>& M, const VectorCL& v)
{
    VectorCL x( 0.0, M.B_.num_cols());
    VectorCL r= transp_mul( M.B_, v);
    M.Apc_.Apply( M.A_, x, r);
    return M.B_*x;
}

//=============================================================================
//  The "Solve" functions
//=============================================================================

template <class PoissonSolverT>
template <typename Mat, typename Vec>
void UzawaSolverCL<PoissonSolverT>::doSolve(
    const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
{
    Vec v_corr(v.size()),
        p_corr(p.size()),
        res1(v.size()),
        res2(p.size());

    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= norm_sq( res2);
        _poissonSolver.SetTol( std::sqrt( res2_norm)/20.0);
        _poissonSolver.Solve(_M, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(_tau, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2(res1, A*v, 1.0, transp_mul(B,p), -1.0, b);
        res1_norm= norm_sq( res1);
        if (res1_norm + res2_norm < tol) {
            _res= std::sqrt( res1_norm + res2_norm );
            return;
        }
        if( (_iter%output)==0)
            std::cout << "step " << _iter << ": norm of 1st eq= " << std::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << std::sqrt( res2_norm) << std::endl;

        _poissonSolver.SetTol( std::sqrt( res1_norm)/20.0);
        _poissonSolver.Solve( A, v_corr, res1);
        v-= v_corr;
    }
    _res= std::sqrt( res1_norm + res2_norm );
}

template <class PoissonSolverT>
void UzawaSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT>
void UzawaSolverCL<PoissonSolverT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT, class PoissonSolver2T>
template <typename Mat, typename Vec>
void UzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::doSolve(
    const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
{
    Vec v_corr( v.size()),
        p_corr( p.size()),
        res1( v.size()),
        res2( p.size());
    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= norm_sq( res2);
        poissonSolver_.SetTol( std::sqrt( res2_norm)/20.0);
        poissonSolver_.Solve( M_, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(tau_, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2( res1, A*v, 1.0, transp_mul( B, p), -1.0, b);
        res1_norm= norm_sq( res1);
        if (res1_norm + res2_norm < tol) {
            _res= std::sqrt( res1_norm + res2_norm);
            return;
        }
        if( (_iter%output)==0)
            std::cout << "step " << _iter << ": norm of 1st eq= " << std::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << std::sqrt( res2_norm) << std::endl;

        poissonSolver2_.SetTol( std::sqrt( res1_norm)/20.0);
        poissonSolver2_.Solve( A, v_corr, res1);
        v-= v_corr;
    }
    _res= std::sqrt( res1_norm + res2_norm );
}

template <class PoissonSolverT, class PoissonSolver2T>
void UzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::Solve(
    const MatrixCL& A, const MatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT, class PoissonSolver2T>
void UzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT>
template <typename Mat, typename Vec>
void SchurSolverCL<PoissonSolverT>::doSolve(
    const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    Vec _tmp;
    Vec rhs(-c);
    if (_tmp.size() != v.size())
        _tmp.resize( v.size());
    _poissonSolver.Solve( A, _tmp, b);
    std::cout << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;
    rhs+= B*_tmp;

    std::cout << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;
    double tol= _tol;
    CG( SchurComplMatrixCL<PoissonSolverT, Mat>( _poissonSolver, A, B), p, rhs, iter, tol);
    std::cout << "iterations: " << iter << "\tresidual: " << tol << std::endl;
    std::cout << "pressure has been solved! Now solving velocities..." << std::endl;

    _poissonSolver.Solve( A, v, Vec(b - transp_mul(B, p)));
    std::cout << "Iterationen: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;

    _iter= iter+_poissonSolver.GetIter();
    _res= tol + _poissonSolver.GetResid();
/* std::cout << "Real residuals are: "
          << norm( A*v+transp_mul(B, p)-b) << ", "
          << norm( B*v-c) << std::endl; */
    std::cout << "-----------------------------------------------------" << std::endl;
}

template <class PoissonSolverT>
void SchurSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT>
void SchurSolverCL<PoissonSolverT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <typename Mat, typename Vec>
inline void Uzawa_IPCG_CL::doSolve(
    const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)

{
    Vec v_corr(v.size()),
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
        res1_norm= norm_sq( res1);
        res2_norm= norm_sq( res2);

        if (res1_norm + res2_norm < tol)
        {
            _res= std::sqrt( res1_norm + res2_norm );
            return;
        }

        if( (_iter%output)==0 )
            std::cout << "step " << _iter << ": norm of 1st eq= " << std::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << std::sqrt( res2_norm) << std::endl;

        _A_IPCGsolver.Solve( A, v_corr, res1);
//        v-= v_corr;
        axpy(-1.0, v_corr, v);
    }
    _res= std::sqrt( res1_norm + res2_norm );
}

inline void Uzawa_IPCG_CL::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

inline void Uzawa_IPCG_CL::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT>
template <typename Mat, typename Vec>
void PSchurSolverCL<PoissonSolverT>::doSolve(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    Vec _tmp;
    Vec rhs( -c);
    if (_tmp.size() != v.size())
        _tmp.resize( v.size());
    _poissonSolver.Solve( A, _tmp, b);
    std::cout << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;
    rhs+= B*_tmp;

    std::cout << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;
    double tol= _tol;
    PCG( SchurComplMatrixCL<PoissonSolverT, Mat>( _poissonSolver, A, B), p, rhs, _schurPc, iter, tol);
    std::cout << "iterations: " << iter << "\tresidual: " << tol << std::endl;
    std::cout << "pressure has been solved! Now solving velocities..." << std::endl;

    _poissonSolver.Solve( A, v, Vec( b - transp_mul(B, p)));
    std::cout << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;

    _iter= iter+_poissonSolver.GetIter();
    _res= std::sqrt( tol*tol + _poissonSolver.GetResid()*_poissonSolver.GetResid());
    std::cout << "-----------------------------------------------------" << std::endl;
}

template <class PoissonSolverT>
void PSchurSolverCL<PoissonSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    doSolve( A, B, v, p, b, c);
}

template <class PoissonSolverT>
void PSchurSolverCL<PoissonSolverT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    doSolve( A, B, v, p, b, c);
}

template <typename InnerSolverT, typename OuterSolverT>
template <typename Mat, typename Vec>
void PSchurSolver2CL<InnerSolverT, OuterSolverT>::doSolve(
    const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    Vec tmp_;
    Vec rhs( -c);
    if (tmp_.size() != v.size()) tmp_.resize( v.size());
    innerSolver_.Solve( A, tmp_, b);
    std::cout << "rhs     : iterations: " << innerSolver_.GetIter()
              << "\tresidual: " << innerSolver_.GetResid() << std::endl;
    rhs+= B*tmp_;

    outerSolver_.SetTol( _tol);
    outerSolver_.SetMaxIter( _maxiter);
    outerSolver_.Solve( SchurComplMatrixCL<InnerSolverT, Mat>( innerSolver_, A, B), p, rhs);
    std::cout << "pressure: iterations: " << outerSolver_.GetIter()
              << "\tresidual: " << outerSolver_.GetResid() << std::endl;

    innerSolver_.Solve( A, v, Vec( b - transp_mul(B, p)));
    std::cout << "velocity: iterations: " << innerSolver_.GetIter()
              << "\tresidual: " << innerSolver_.GetResid() << std::endl;

    _iter= innerSolver_.GetIter() + outerSolver_.GetIter();
    _res= std::sqrt( std::pow( innerSolver_.GetResid(), 2)
                   + std::pow( outerSolver_.GetResid(), 2));
    std::cout << "-----------------------------------------------------" << std::endl;
}

template <typename InnerSolverT, typename OuterSolverT>
void PSchurSolver2CL<InnerSolverT, OuterSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}

template <typename InnerSolverT, typename OuterSolverT>
void PSchurSolver2CL<InnerSolverT, OuterSolverT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    doSolve( A, B, v, p, b, c);
}


//-----------------------------------------------------------------------------
// UzawaPCG: The basic scheme is identical to PCG, but two extra recursions for
//     InexactUzawa are performed to avoid an evaluation of a preconditioner
//     there.
//     Also, the matrix is fixed: BQ_A^{-1}B^T and we assume x=zbar=zhat=0
//     upon entry to this function.
//     See "Fast iterative solvers for Stokes equation", Peters, Reichelt,
//     Reusken, Ch. 3.3 Remark 3.
//-----------------------------------------------------------------------------
template <typename APC, typename Mat, typename Vec, typename SPC>
bool
UzawaPCG(const APC& Apc, const Mat& A, const Mat& B,
    Vec& x, Vec& zbar, Vec& zhat, const Vec& b, const SPC& M,
    int& max_iter, double& tol)
{
    const size_t n= x.size();
    Vec p(n), z(n), q(n), d(n), e(n), r= b;
    Vec q1( B.num_cols()), q2( B.num_cols());
    double rho, rho_1= 0.0, resid= norm_sq( r);

    tol*= tol;
    if (resid<=tol)
    {
        tol= std::sqrt( resid);
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        M.Apply( B/*dummy*/, z, r);
        rho= dot( r, z);
        if (i == 1)
            p= z;
        else
            z_xpay(p, z, (rho/rho_1), p); // p= z + (rho/rho_1)*p;

        // q= A*p;
        q1= transp_mul( B, p);
        q2= 0.0;
        Apc.Apply( A, q2, q1);
        q= B*q2;

        const double alpha= rho/dot( p, q);
        axpy(alpha, p, x);                // x+= alpha*p;
        axpy(-alpha, q, r);               // r-= alpha*q;
        zbar+= alpha*q1;
        zhat+= alpha*q2;
        resid= norm_sq( r);
        if (resid<=tol)
        {
            tol= std::sqrt( resid);
            max_iter= i;
            return true;
        }
        rho_1= rho;
    }
    tol= std::sqrt(resid);
    return false;
}

//=============================================================================
// Inexact Uzawa-method from "Fast Iterative Solvers for Discrete Stokes
// Equations", Peters, Reichelt, Reusken, Chapter 3.3.
// The preconditioner Apc for A must be "good" (MG-like) to guarantee
// convergence.
// Due to theory (see paper), we should use 0.2 < innerred < 0.7.
//=============================================================================
template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
InexactUzawa(const Mat& A, const Mat& B, Vec& xu, Vec& xp, const Vec& f, const Vec& g,
    PC1& Apc, PC2& Spc,
    int& max_iter, double& tol,
    InexactUzawaApcMethodT apcmeth= APC_OTHER,
    double innerred= 0.3, int innermaxiter= 500)
{
    VectorCL ru( f - A*xu - transp_mul( B, xp));
    VectorCL rp( g - B*xu);
    VectorCL w( f.size());
    VectorCL z( g.size());
    VectorCL z2( g.size());
    VectorCL zbar( f.size());
    VectorCL zhat( f.size());
    VectorCL du( f.size());
    VectorCL c( g.size());
    ApproximateSchurComplMatrixCL<PC1, Mat>* asc= apcmeth == APC_SYM_LINEAR ? 0 :
        new ApproximateSchurComplMatrixCL<PC1, Mat>( A, Apc, B);
    double innertol;
    int inneriter;
    int pr_iter_cumulative= 0;
    double resid0= std::sqrt( norm_sq( ru) + norm_sq( rp));
    double resid= resid0;
    std::cout << "residual (2-norm): " << resid
              << "\tres-impuls: " << norm( ru)
              << "\tres-mass: " << norm( rp)
              << '\n';
    if (resid <= tol) { // The fixed point iteration between levelset and Stokes
        tol= resid;     // equation uses this to determine convergence.
        max_iter= 0;
        delete asc;
        return true;
    }
    for (int k= 1; k <= max_iter; ++k) {
        w= 0.0;
        Apc.Apply( A, w, ru);
        c= B*w - rp;
        z= 0.0;
        z2= 0.0;
        inneriter= innermaxiter;
        switch (apcmeth) {
          case APC_SYM_LINEAR:
            zbar= 0.0;
            zhat= 0.0;
            innertol= innerred*norm( c);
            UzawaPCG( Apc, A, B, z, zbar, zhat, c, Spc, inneriter, innertol);
            break;
          case APC_SYM:
            innertol= innerred*norm( c);
            PCG( *asc, z, c, Spc, inneriter, innertol);
            break;
          default:
            std::cout << "WARNING: InexactUzawa: Unknown apcmeth; using GMRes.\n";
            // fall through
          case APC_OTHER:
            innertol= innerred; // GMRES can do relative tolerances.
            GMRES( *asc, z, c, Spc, /*restart*/ inneriter, inneriter, innertol,
                /*relative errors*/ true, /*don't check 2-norm*/ false);
              break;
        }
        if (apcmeth != APC_SYM_LINEAR) {
            zbar= transp_mul( B, z);
            zhat= 0.0;
            Apc.Apply( A, zhat, zbar);
        }
        pr_iter_cumulative+= inneriter;
        std::cout << "pr solver: iterations: " << pr_iter_cumulative
                  << "\tresid: " << innertol << '\n';
        du= w - zhat;
        xu+= du;
        xp+= z;
        ru-= A*du + zbar; // z_xpaypby2(ru, ru, -1.0, A*du, -1.0, zbar);
        rp= g - B*xu;
        resid= std::sqrt( norm_sq( ru) + norm_sq( rp));
        std::cout << "residual reduction (2-norm): " << resid/resid0
                  << "\nresidual (2-norm): " << resid
                  << "\tres-impuls: " << norm( ru)
                  << "\tres-mass: " << norm( rp)
                  << '\n';
        if (resid <= tol) { // absolute errors
            tol= resid;
            max_iter= k;
            delete asc;
            return true;
        }
    }
    tol= resid;
    delete asc;
    return false;
}


template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const MatrixCL& A, const MatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    _res=  _tol;
    _iter= _maxiter;
    InexactUzawa( A, B, v, p, b, c, Apc_, Spc_, _iter, _res, Apcmeth, innerreduction_, innermaxiter_);
}

template <class ApcT, class SpcT, InexactUzawaApcMethodT Apcmeth>
  inline void
  InexactUzawaCL<ApcT, SpcT, Apcmeth>::Solve( const MLMatrixCL& A, const MLMatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    _res=  _tol;
    _iter= _maxiter;
    InexactUzawa( A, B, v, p, b, c, Apc_, Spc_, _iter, _res, Apcmeth, innerreduction_, innermaxiter_);
}

//-----------------------------------------------------------------------------
// UzawaPCG: The basic scheme is identical to PCG, but two extra recursions for
//     InexactUzawa are performed to avoid an evaluation of a preconditioner
//     there.
//     Also, the matrix is fixed: BQ_A^{-1}B^T and we assume x=zbar=zhat=0
//     upon entry to this function.
//     See "Fast iterative solvers for Stokes equation", Peters, Reichelt,
//     Reusken, Ch. 3.3 Remark 3.
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
UzawaCGEff(const Mat& A, const Mat& B, Vec& xu, Vec& xp, const Vec& f, const Vec& g,
    PC1& Apc, PC2& Spc,
    int& max_iter, double& tol)
{
    double err= std::sqrt( norm_sq( f - ( A*xu + transp_mul( B, xp))) + norm_sq( g - B*xu));
    const double err0= err;
    Vec rbaru( f - (A*xu + transp_mul(  B, xp)));
    Vec rbarp( g - B*xu);
    Vec ru( f.size());
    Apc.Apply( A, ru, rbaru);
    Vec rp( B*ru - rbarp);
    Vec a( f.size()), b( f.size()), s( f.size()), pu( f.size()), qu( f.size());
    Vec z( g.size()), pp( g.size()), qp( g.size()), t( g.size());
    double alpha= 0.0, initialbeta=0.0, beta= 0.0, beta0= 0.0, beta1= 0.0;
    for (int i= 0; i < max_iter; ++i) {
        z= 0.0;
        Spc.Apply( B, z, rp);
        a= A*ru;
        beta1= dot(a, ru) - dot(rbaru,ru) + dot(z,rp);
        if (i==0) initialbeta= beta1;
        if (beta1 <= 0.0) {throw DROPSErrCL( "UzawaCGEff: Matrix is not spd.\n");}
        // This is for fair comparisons of different solvers:
        err= std::sqrt( norm_sq( f - (A*xu + transp_mul( B, xp))) + norm_sq( g - B*xu));
        std::cout << "relative residual (2-norm): " << err/err0
                  << "\t(problem norm): " << std::sqrt( beta1/initialbeta) << '\n';
//        if (beta1/initialbeta <= tol*tol) {
//            tol= std::sqrt( beta1/initialbeta);
        if (err/err0 <= tol) {
            tol= err/err0;
            max_iter= i;
            return true;
        }
        if (i != 0) {
            beta= beta1/beta0;
            z_xpay( pu, ru, beta, pu); // pu= ru + beta*pu;
            z_xpay( pp, z, beta, pp);  // pp= z  + beta*pp;
            z_xpay( b, a, beta, b);    // b= a + beta*b;
        }
        else {
            pu= ru;
            pp= z;
            b= a; // A*pu;
        }
        qu= b + transp_mul( B, pp);
        qp= B*pu;
        s= 0.0;
        Apc.Apply( A, s, qu);
        z_xpay( t, B*s, -1.0, qp); // t= B*s - qp;
        alpha= beta1/(dot( s, b) - dot(qu, pu) + dot(t, pp));
        axpy( alpha, pu, xu); // xu+= alpha*pu;
        axpy( alpha, pp, xp); // xp+= alpha*pp;
        axpy( -alpha, s, ru); // ru-= alpha*s;
        axpy( -alpha, t, rp); // rp-= alpha*t;
        axpy( -alpha, qu, rbaru);// rbaru-= alpha*qu;
        // rbarp-= alpha*qp; rbarp is not needed in this algo.
        beta0= beta1;
    }
    tol= std::sqrt( beta1);
    return false;
}

//-----------------------------------------------------------------------------
// CG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// max_iter - number of iterations performed before tolerance was reached
//      tol - residual after the final iteration
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
UzawaCG(const Mat& A, const Mat& B, Vec& u, Vec& p, const Vec& b, const Vec& c,
        PC1& Apc, PC2& Spc, int& max_iter, double& tol)
{
    double err= std::sqrt( norm_sq( b - (A*u + transp_mul( B, p))) + norm_sq( c - B*u));
    const double err0= err;
    Vec ru( b - ( A*u + transp_mul(  B, p)));
    Vec rp( c - B*u);
    Vec s1( b.size()); // This is r2u...
    Apc.Apply( A, s1, ru);
    Vec s2( B*s1 - rp);
    Vec r2p( c.size());
    Spc.Apply( B, r2p, s2);
    double rho0= dot( s1, VectorCL( A*s1 - ru)) + dot( r2p, s2);
    const double initialrho= rho0;
//    std::cout << "UzawaCG: rho: " << rho0 << '\n';
    if (rho0<=0.0) throw DROPSErrCL("UzawaCG: Matrix is not spd.\n");
//    tol*= tol*rho0*rho0; // For now, measure the relative error.
//    if (rho0<=tol) {
//        tol= std::sqrt( rho0);
//        max_iter= 0;
//        return true;
//    }
    Vec pu= s1; // s1 is r2u.
    Vec pp= r2p;
    Vec qu( A*pu + transp_mul( B, pp));
    Vec qp= B*pu;
    double rho1= 0.0;
    Vec t1( b.size());
    Vec t2( c.size());
    for (int i= 1; i<=max_iter; ++i) {
        Apc.Apply( A, t1, qu);
        z_xpay( t2, B*t1, -1.0, qp); // t2= B*t1 - qp;
        const double alpha= rho0/( dot(pu, VectorCL( A*t1 - qu)) + dot( pp, t2));
        axpy(alpha, pu, u);  // u+= alpha*pu;
        axpy(alpha, pp, p);  // p+= alpha*pp;
        axpy( -alpha, qu, ru);
        axpy( -alpha, qp, rp);
        s1= 0.0;
        Apc.Apply( A, s1, ru);
        z_xpay( s2, B*s1, -1.0, rp); // s2= B*s1 - rp;
        //axpy( -alpha, t1, s1); // kann die beiden oberen Zeilen ersetzen,
        //axpy( -alpha, t2, s2); // Algorithmus wird schneller, aber Matrix bleibt nicht spd
        r2p= 0.0;
        Spc.Apply( B, r2p, s2);
        rho1= dot( s1, VectorCL( A*s1 - ru)) + dot( r2p, s2);
//        std::cout << "UzawaCG: rho: " << rho1 << '\n';
        if (rho1<=0.0) throw DROPSErrCL("UzawaCG: Matrix is not spd.\n");
        // This is for fair comparisons of different solvers:
        err= std::sqrt( norm_sq( b - (A*u + transp_mul( B, p))) + norm_sq( c - B*u));
        std::cout << "relative residual (2-norm): " << err/err0
                << "\t(problem norm): " << std::sqrt( rho1/initialrho)<< '\n';
//        if (rho1 <= tol) {
//            tol= std::sqrt( rho1);
        if (err <= tol*err0) {
            tol= err/err0;
            max_iter= i;
            return true;
        }
        const double beta= rho1/rho0;
        z_xpay( pu, s1, beta, pu); // pu= s1 + beta*pu; // s1 is r2u.
        z_xpay( pp, r2p, beta, pp); // pp= r2p + beta*pp;
        qu= (A*pu) + transp_mul( B, pp);
        qp= (B*pu);
        rho0= rho1;
        t1= 0.0;
    }
    tol= std::sqrt( rho1);
    return false;
}

template <typename PC1, typename PC2>
class UzawaCGSolverEffCL : public StokesSolverBaseCL
{
  private:
    PC1& Apc_;
    PC2& Spc_;

  public:
    UzawaCGSolverEffCL (PC1& Apc, PC2& Spc, int maxiter, double tol)
        : StokesSolverBaseCL(maxiter, tol), Apc_( Apc), Spc_( Spc) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        UzawaCGEff( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        UzawaCGEff( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
};

template <typename PC1, typename PC2>
class UzawaCGSolverCL : public StokesSolverBaseCL
{
  private:
    PC1& Apc_;
    PC2& Spc_;

  public:
    UzawaCGSolverCL (PC1& Apc, PC2& Spc, int maxiter, double tol)
        : StokesSolverBaseCL(maxiter, tol), Apc_( Apc), Spc_( Spc) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        UzawaCG( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        UzawaCG( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
};

//important for ScaledMGPreCL
template <class ProlongationT>
double
EigenValueMaxMG(const MLMatrixCL& A, const ProlongationT& P, VectorCL& x, int iter, double  tol= 1e-3)
{
    MLMatrixCL::const_iterator finest= A.GetFinestIter();
    typename ProlongationT::const_iterator finestP= P.GetFinestIter();
    Uint   sm   =  1; // how many smoothing steps?
    int    lvl  = -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
    double l= -1.0;
    double l_old= -1.0;
    VectorCL tmp( x.size());
    VectorCL z( x.size());
//    JORsmoothCL smoother( omega); // Jacobi
//    GSsmoothCL smoother( omega); // Gauss-Seidel
//    SGSsmoothCL smoother( omega); // symmetric Gauss-Seidel
//    SORsmoothCL smoother( omega); // Gauss-Seidel with over-relaxation
    SSORsmoothCL smoother( omega); // symmetric Gauss-Seidel with over-relaxation
//    CGSolverCL  solver( 200, tol); //CG-Verfahren
    SSORPcCL directpc; PCG_SsorCL solver( directpc, 200, 1e-15);
    x/= norm( x);
    std::cout << "EigenValueMaxMG:\n";
    for (int i= 0; i<iter; ++i) {
        tmp= 0.0;
        MGM( A.begin(), finest, finestP, tmp, A*x, smoother, sm, solver, lvl, -1);
        z= x - tmp;
        l= dot( x, z);
        std::cout << "iteration: " << i  << "\tlambda: " << l << "\trelative_change= : " << (i==0 ? -1 : std::fabs( (l-l_old)/l_old)) << '\n';
        if (i > 0 && std::fabs( (l-l_old)/l_old) < tol) break;
        l_old= l;
        x= z/norm( z);
    }
    std::cout << "maximal value for lambda: " << l << '\n';
    return l;
}

//scaled MG as preconditioner
template <class ProlongationT= MLMatrixCL>
class ScaledMGPreCL
{
  private:
    const ProlongationT& P_;
    Uint iter_;
    double s_;
    mutable SSORPcCL directpc_;
    mutable PCG_SsorCL solver_;
    Uint sm_; // how many smoothing steps?
    int lvl_; // how many levels? (-1=all)
    double omega_; // relaxation parameter for smoother
    mutable SSORsmoothCL smoother_;  // Symmetric-Gauss-Seidel with over-relaxation

  public:
    ScaledMGPreCL( const ProlongationT& P, Uint iter, double s= 1.0)
        : P_(P), iter_( iter), s_( s), solver_( directpc_, 200, 1e-12), sm_( 1),
         lvl_( -1), omega_( 1.0), smoother_( omega_)
    {}

    inline void
    Apply( const MLMatrixCL& A, VectorCL& x, const VectorCL& r) const;
    inline void
    Apply( const MatrixCL&, VectorCL&, const VectorCL&) const
    {
        throw DROPSErrCL( "ScaledMGPreCL::Apply: need multilevel data structure\n");
    }
};

template <class ProlongationT>
inline void
ScaledMGPreCL<ProlongationT>::Apply( const MLMatrixCL& A, VectorCL& x, const VectorCL& r) const
{
    x= 0.0;
    MLMatrixCL::const_iterator finest = A.GetFinestIter();
    typename ProlongationT::const_iterator finestP= P_.GetFinestIter();
    const double oldres= norm(r - A*x);
    for (Uint i= 0; i < iter_; ++i) {
        VectorCL r2( r - A*x);
        VectorCL dx( x.size());
        MGM( A.begin(), finest, finestP, dx, r2, smoother_, sm_, solver_, lvl_, -1);
        x+= dx*s_;
    }
    const double res= norm(r - A*x);
    std::cout << "ScaledMGPreCL: it: " << iter_ << "\treduction: " << (oldres==0.0 ? res : res/oldres) << '\n';
}

template<class SmootherT, class PVelT= MLMatrixCL, class PPrT= MLMatrixCL>
class StokesMGSolverCL: public StokesSolverBaseCL
{
  private:
    PVelT PVel_;
    PPrT  PPr_;
    const MLMatrixCL    &prM_;
    const SmootherT&    smoother_;
    StokesSolverBaseCL& directSolver_;
    Uint  smoothSteps_;
    int   usedLevels_;
    MLMatrixCL BT_;
    size_t BVersion_;

    void UpdateBT( const MLMatrixCL& B)
    {
        BT_.resize( B.size());
        MLMatrixCL::iterator BT = BT_.begin();
        for ( MLMatrixCL::const_iterator it = B.begin(); it != B.end(); ++it)
        {
            transpose (*it, *BT);
            ++BT;
        }
        BVersion_ = B.Version();
    }

  public:
    StokesMGSolverCL( const MLMatrixCL& prM, const SmootherT& smoother, StokesSolverBaseCL& ds,
                      Uint iter_vel, double tol, bool rel= false, Uint sm = 2, int lvl = -1)
      : StokesSolverBaseCL(iter_vel, tol, rel), prM_(prM), smoother_(smoother), directSolver_(ds),
        smoothSteps_(sm), usedLevels_(lvl), BVersion_( 0) {}
    ~StokesMGSolverCL() {}

    PVelT* GetPVel() { return &PVel_;}
    PPrT*  GetPPr()  { return &PPr_; }
    void
    Solve(const MatrixCL& /*A*/, const MatrixCL& /*B*/, VectorCL&, VectorCL&, const VectorCL&, const VectorCL&)
    {
        throw DROPSErrCL( "StokesMGSolverCL::Solve: need multilevel data structure\n");
    }

    void
    Solve(const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c) {
// define MG parameters for the first diagonal blockS
        if (B.Version() != BVersion_) UpdateBT( B);
        int nit=_maxiter;
        double actualtol = 1;
        Uint   wc   = 1;   // how many W-cycle steps? (1=V-cycle)

// define initial approximation
        MLMatrixCL::const_iterator A_end   ( A.GetFinestIter());
        MLMatrixCL::const_iterator B_end   ( B.GetFinestIter());
        MLMatrixCL::const_iterator BT_end  ( BT_.GetFinestIter());
        MLMatrixCL::const_iterator prM_end ( prM_.GetFinestIter());
        typename PVelT::const_iterator PVel( PVel_.GetFinestIter());
        typename PPrT::const_iterator  PPr ( PPr_.GetFinestIter());

        const double runorm0= norm_sq( A * v + transp_mul( B, p) - b);
        const double rpnorm0= norm_sq( B * v - c);
        const double resid0= std::sqrt(runorm0+rpnorm0);
        double resid;
        for (int j=0; j<_maxiter; ++j)
        {
            StokesMGM( A.begin(), A_end, B_end, BT_end, prM_end, PVel, PPr, v, p, b, c, smoother_, smoothSteps_, wc, directSolver_, usedLevels_, -1);
            const double runorm= norm_sq( A * v + transp_mul(B, p ) - b);
            const double rpnorm= norm_sq( B * v - c);
            resid= std::sqrt(runorm+rpnorm);
            if (rel_)
                actualtol= resid/resid0;
            else
                actualtol= resid;
            std::cout << "P2P1:StokesMGSolverCL: residual = " << actualtol << std::endl;
            if (actualtol<=_tol)
            {
                nit= j+1;
                break;
            }
        }
        std::cout << "StokesMGM: actual residual = " << actualtol
                  << "  after " << nit << " iterations " << std::endl;
        _res = actualtol;
        _iter= nit;
   }
};

///\brief Preconditioner for the (Navier-) Stokes operator
///
/// Vertex-based multiplikative Schwartz-method. The implementation is based on the Vanka-smoother class for
/// Stokes-multigrid in Drops. The effect of Apply is one smoother-iteration with 0 as starting value.
class VankaPreCL
{
  private:
    PVankaSmootherCL smoother_;
    mutable MatrixCL BT_;
    mutable size_t BVersion_;

    void UpdateBT (const MatrixCL& B) const {
        transpose (B, BT_);
        BVersion_ = B.Version();
    }

  public:
    VankaPreCL (const MLIdxDescCL* idx= 0) : smoother_( 0, 1., idx), BVersion_( 0) {}

    void Setidx (const MLIdxDescCL* idx) { smoother_.Setidx( idx); }

    void
    Apply (const BlockMatrixCL& M, VectorCL& x, const VectorCL& rhs) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)),
                 b(rhs[std::slice( 0, M.num_rows( 0), 1)]),
                 c(rhs[std::slice( M.num_rows( 0), M.num_rows( 1), 1)]);
        if (M.GetBlock( 1)->Version() != BVersion_) UpdateBT( *M.GetBlock( 1));
        smoother_.Apply( *M.GetBlock( 0), *M.GetBlock( 1), BT_, /*dummy*/ *M.GetBlock( 0), v, p, b, c);
        x[std::slice( 0, M.num_cols( 0), 1)]= v;
        x[std::slice( M.num_cols( 0), M.num_cols( 1), 1)]= p;
    }
    void
    Apply (const MLBlockMatrixCL& M, VectorCL& x, const VectorCL& rhs) const {
        BlockMatrixCL MM( &M.GetBlock( 0)->GetFinest(), MUL, &M.GetBlock( 1)->GetFinest(), TRANSP_MUL,
                          &M.GetBlock( 1)->GetFinest(), MUL);
        this->Apply( MM, x, rhs);
    }
};

///\brief Preconditioner for the Schur-complement of the (Navier-) Stokes operator
///
/// Vertex-based multiplikative Schwartz-method. The implementation is based on the Vanka-smoother class for
/// Stokes-multigrid in Drops. The effect of Apply is one smoother-iteration with 0 as starting value. The
/// right-hand side for the momentum equations is set to 0. Thus, the action of the Schur-complement is approximated.
class VankaSchurPreCL
{
  private:
    PVankaSmootherCL smoother_;
    mutable MatrixCL BT_;
    mutable size_t BVersion_;
    mutable BlockMatrixCL M;

    void UpdateBT (const MatrixCL& B) const {
        transpose (B, BT_);
        BVersion_ = B.Version();
    }

  public:
    VankaSchurPreCL (const MLIdxDescCL* idx= 0)
        : smoother_( 0, 1., idx), BVersion_( 0), M( 0, MUL, 0, TRANSP_MUL, 0, MUL) {}

    void Setidx (const MLIdxDescCL* idx)   { smoother_.Setidx( idx); }
    void SetAB   (const MatrixCL* A, const MatrixCL* B) { M.SetBlock( 0, A); M.SetBlock( 1, B); M.SetBlock( 2, B); }

    void
    Apply (const MatrixCL& /*S*/, VectorCL& x, const VectorCL& rhs) const {
        VectorCL v(M.num_cols( 0)), p(M.num_cols( 1)), b( M.num_rows( 0));
        if (M.GetBlock( 1)->Version() != BVersion_) UpdateBT( *M.GetBlock( 1));
        smoother_.Apply( *M.GetBlock( 0), *M.GetBlock( 1), BT_, /*dummy*/ *M.GetBlock( 0), v, p, b, rhs);
        x= p;
    }
    void
    Apply (const MLMatrixCL& /*S*/, VectorCL& x, const VectorCL& rhs) const {
        this->Apply( MatrixCL(), x, rhs);
    }
};

} // end of namespace DROPS

#endif
