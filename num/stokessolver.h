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
#include "stokes/integrTime.h"

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
    VectorCL        _tmp;

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
    VectorCL               _tmp;

  public:
    PSchurSolverCL (PoissonSolverT& solver, MatrixCL& M, int maxiter, double tol)
        : SolverBaseCL(maxiter,tol), _poissonSolver(solver), _schurPc(M) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <typename InnerSolverT, typename OuterSolverT>
class PSchurSolver2CL : public SolverBaseCL
{
  private:
    InnerSolverT& innerSolver_;
    OuterSolverT& outerSolver_;
    VectorCL      tmp_;

  public:
    PSchurSolver2CL( InnerSolverT& solver1, OuterSolverT& solver2,
                    int maxiter, double tol)
        : SolverBaseCL( maxiter, tol), innerSolver_( solver1),
          outerSolver_( solver2) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
                
    InnerSolverT& GetInnerSolver() { return innerSolver_; }
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
// Ready-to-use solver-class with inexact Uzawa method InexactUzawa. 
//=============================================================================
// Inexact Uzawa-method from "Fast Iterative Solvers for Discrete Stokes
// Equations", Peters, Reichelt, Reusken, Chapter 3.3.
// The preconditioner Apc for A must be "good" (MG-like) to guarantee
// convergence.
//=============================================================================
class InexactUzawa_CL: public SolverBaseCL
{
  private:
    SSORPCG_PreCL Apc_;
    ISPreCL& Spc_;

  public:
    InexactUzawa_CL(ISPreCL& Spc, int outer_iter, double outer_tol)
        :SolverBaseCL( outer_iter, outer_tol),
	 Apc_( 1000, 0.2), Spc_( Spc)
    {}
    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};




// One recursive step of Lanzcos' algorithm for computing an ONB (q1, q2, q3,...)
// of the Krylovspace of A for a given starting vector r. This is a three term
// recursion, computing the next q_i from the two previous ones.
// See Arnold Reusken, "Numerical methods for elliptic partial differential equations",
// p. 148.
//
// This is the reformulation from ../num/solver.h for the Stokes-equation.
//
// Returns false for 'lucky breakdown' (see below), true in the generic case.
template <typename Mat, typename Vec>
bool
LanczosStep_SP(const Mat& A, const Mat& B,
               const Vec& qu0, const Vec& qpr0, const Vec& qu1, const Vec& qpr1,
               Vec& qu2, Vec& qpr2,
               double& a1,
               const double b0, double& b1)
{
    qu2.raw()= (A*qu1).raw() + transp_mul( B, qpr1).raw() - b0*qu0.raw();
    qpr2.raw()= (B*qu1).raw() - b0*qpr0.raw();
    a1= qu2*qu1 + qpr2*qpr1;
    qu2.raw()-= a1*qu1.raw();
    qpr2.raw()-= a1*qpr1.raw();
    b1= std::sqrt( qu2.norm2() + qpr2.norm2());
    // Lucky breakdown; the Krylov-space K up to q1 is A-invariant. Thus,
    // the correction dx needed to solve A(x0+dx)=b is in this space and
    // the Minres-algo will terminate with the exact solution in the
    // following step.
    if (b1 < 1e-15) return false;
    qu2.raw()*= 1./b1;
    qpr2.raw()*= 1./b1;
    return true;
}

template <typename Mat, typename Vec>
class LanczosONB_SPCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;
    const Mat* A;
    const Mat* B;

  public:
    SBufferCL<Vec, 3> qu;
    SBufferCL<Vec, 3> qpr;
    double a0;
    SBufferCL<double, 2> b;

    LanczosONB_SPCL()
      :A( 0), B( 0) {}

    void // Sets up initial values and computes q0.
    new_basis(const Mat& A_, const Mat& B_, const Vec& ru0, const Vec& rpr0) {
        A= &A_;
        B= &B_;
        qu[-1].resize( ru0.size(), 0.);
        qpr[-1].resize( rpr0.size(), 0.);
        norm_r0_= std::sqrt( ru0.norm2() + rpr0.norm2());
        qu[0].resize( ru0.size(), 0.); qu[0].raw()= ru0.raw()*(1./norm_r0_);
        qpr[0].resize( rpr0.size(), 0.); qpr[0].raw()= rpr0.raw()*(1./norm_r0_);
        qu[1].resize( ru0.size(), 0.);
        qpr[1].resize( rpr0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= LanczosStep_SP( *A, *B, qu[-1], qpr[-1], qu[0], qpr[0],
                                      qu[1], qpr[1], a0, b[-1], b[0]);
    }

    double norm_r0() const {
        return norm_r0_; }
    bool
    breakdown() const {
        return !nobreakdown_; }
    // Computes new q_i, a_i, b_1, q_{i+1} in q0, a0, b0, q1 and moves old
    // values to qm1, bm1.
    bool
    next() {
        qu.rotate(); qpr.rotate(); b.rotate();
        return (nobreakdown_= LanczosStep_SP( *A, *B, qu[-1], qpr[-1], qu[0], qpr[0],
                                              qu[1], qpr[1], a0, b[-1], b[0]));
    }
};

// One recursive step of the preconditioned Lanzcos algorithm for computing an ONB of
// a Krylovspace. This is a three term
// recursion, computing the next q_i from the two previous ones.
// See Arnold Reusken, "Numerical methods for elliptic partial differential equations",
// p. 153.
//
// This is the PMinres-method from ../num/solver.h adapted to the
// Stokes-equations. Keep in sync!
// Returns false for 'lucky breakdown' (see below), true in the generic case.
template <typename Mat, typename Vec, typename PreCon>
bool
PLanczosStep_SP(const Mat& A, const Mat& B,
                const PreCon& M,
                const Vec& qu1, const Vec& qpr1, Vec& qu2, Vec& qpr2,
                const Vec& tu0, const Vec& tpr0, const Vec& tu1, const Vec& tpr1,
                Vec& tu2, Vec& tpr2,
                double& a1,
                const double b0, double& b1)
{
    tu2.raw()= (A*qu1).raw() + transp_mul( B, qpr1).raw() - b0*tu0.raw();
    tpr2.raw()= (B*qu1).raw() - b0*tpr0.raw();
    a1= tu2*qu1 + tpr2*qpr1;
    tu2.raw()+= (-a1)*tu1.raw();
    tpr2.raw()+= (-a1)*tpr1.raw();
    M.Apply( A, B, qu2, qpr2, tu2, tpr2);
    const double b1sq= qu2*tu2 + qpr2*tpr2;
    Assert( b1sq >= 0.0, "PLanczosStep_SP: b1sq is negative!\n", DebugNumericC);
    b1= std::sqrt( b1sq);
    if (b1 < 1e-15) return false;
    tu2.raw()*= 1./b1;
    tpr2.raw()*= 1./b1;
    qu2.raw()*= 1./b1;
    qpr2.raw()*= 1./b1;
    return true;
}

template <typename Mat, typename Vec, typename PreCon>
class PLanczosONB_SPCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;
    const Mat* A;
    const Mat* B;
    const PreCon& M;

  public:
    SBufferCL<Vec, 2> qu;
    SBufferCL<Vec, 2> qpr;
    SBufferCL<Vec, 3> tu;
    SBufferCL<Vec, 3> tpr;
    double a0;
    SBufferCL<double, 2> b;

    PLanczosONB_SPCL(const PreCon& M_)
        :A( 0), B( 0), M( M_) {}

    void // Sets up initial values and computes q0.
    new_basis(const Mat& A_, const Mat& B_,  const Vec& r0u, const Vec& r0pr) {
        A= &A_;
        B= &B_;
        tu[-1].resize( r0u.size(), 0.);
        tpr[-1].resize( r0pr.size(), 0.);
        qu[-1].resize( r0u.size(), 0.);
        qpr[-1].resize( r0pr.size(), 0.);
        M.Apply( *A, *B, qu[-1], qpr[-1], r0u, r0pr);
        norm_r0_= std::sqrt( qu[-1]*r0u + qpr[-1]*r0pr);
        tu[0].resize( r0u.size(), 0.); tu[0].raw()= r0u.raw()*(1./norm_r0_);
        tpr[0].resize( r0pr.size(), 0.); tpr[0].raw()= r0pr.raw()*(1./norm_r0_);
        qu[0].resize( r0u.size(), 0.); qu[0].raw()= qu[-1].raw()*(1./norm_r0_);
        qpr[0].resize( r0pr.size(), 0.); qpr[0].raw()= qpr[-1].raw()*(1./norm_r0_);
        tu[1].resize( r0u.size(), 0.);
        tpr[1].resize( r0pr.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= PLanczosStep_SP( *A, *B, M, qu[0], qpr[0], qu[1], qpr[1],
                          tu[-1], tpr[-1], tu[0], tpr[0], tu[1], tpr[1],
                          a0, b[-1], b[0]);
    }

    double
    norm_r0() const {
        return norm_r0_; }
    bool
    breakdown() const {
        return !nobreakdown_; }
    // Computes new q_i, t_i, a_i, b_1, q_{i+1} in q0, t_0, a0, b0, q1 and moves old
    // values to qm1, tm1, bm1.
    bool
    next() {
        qu.rotate(); qpr.rotate(); tu.rotate(); tpr.rotate(); b.rotate();
        return (nobreakdown_= PLanczosStep_SP( *A, *B, M, qu[0], qpr[0], qu[1], qpr[1],
                                  tu[-1], tpr[-1], tu[0], tpr[0], tu[1], tpr[1],
                                  a0, b[-1], b[0]));
    }
};


//-----------------------------------------------------------------------------
// PMINRES: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// See Arnold Reusken, "Numerical methods for elliptic partial differential
// equations", pp. 149 -- 154
//
// This is the PMinres-method from ../num/solver.h adapted to the
// Stokes-equations. Keep in sync!
//
// Upon successful return, output arguments have the following values:
//
//        x - approximate solution to Ax = rhs
// max_iter - number of iterations performed before tolerance was reached
//      tol - 2-norm of the last correction dx to x after the final iteration.
//-----------------------------------------------------------------------------
template <typename Mat, typename Vec, typename Lanczos>
bool
PMINRES_SP(const Mat& /*A*/, const Mat& /*B*/,
           Vec& u, Vec& pr, const Vec& /*rhsu*/, const Vec& /*rhspr*/,
           Lanczos& q, int& max_iter, double& tol)
{
    Vec dxu( u.size()); // Vector for updating x, thus the name.
    Vec dxpr( pr.size()); // Vector for updating x, thus the name.
//    double err= 0.0; // Sink gcc-warnings.
    double res= 0.0; // Sink gcc-warnings.

    const double norm_r0= q.norm_r0();
    bool lucky= q.breakdown();
    SBufferCL<double, 3> c;
    SBufferCL<double, 3> s;
    SBufferCL<SVectorCL<3>, 3> r;
    SBufferCL<Vec, 3> pu;
    SBufferCL<Vec, 3> ppr;
    pu[0].resize( u.size()); pu[1].resize( u.size()); pu[2].resize( u.size());
    ppr[0].resize( pr.size()); ppr[1].resize( pr.size()); ppr[2].resize( pr.size());
    SBufferCL<SVectorCL<2>, 2> b;

    if (norm_r0 < tol) { // Needed as a stopping criterion in the Stokes-levelset-fp-iteration.
            tol= norm_r0;
            max_iter= 0;
            return true;
        }

    for (int k=1; k<=max_iter; ++k) {
        switch (k) {
          case 1:
            // Compute r1
            GMRES_GeneratePlaneRotation( q.a0, q.b[0], c[0], s[0]);
            r[0][0]= std::sqrt( q.a0*q.a0 + q.b[0]*q.b[0]);
            // Compute p1
            // p[0]= q.q[0]/r[0][0];
            pu[0].raw()= q.qu[0].raw()/r[0][0];
            ppr[0].raw()= q.qpr[0].raw()/r[0][0];
            // Compute b11
            b[0][0]= 1.; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation(b[0][0], b[0][1], c[0], s[0]);
            break;
          case 2:
            // Compute r2
            r[0][0]= q.b[-1]; r[0][1]= q.a0; r[0][2]= q.b[0];
            GMRES_ApplyPlaneRotation( r[0][0], r[0][1], c[-1], s[-1]);
            GMRES_GeneratePlaneRotation( r[0][1], r[0][2], c[0], s[0]);
            GMRES_ApplyPlaneRotation( r[0][1], r[0][2], c[0], s[0]);
            // Compute p2
            // p[0]= (q.q[0] - r[0][0]*p[-1])/r[0][1];
            pu[0].raw()= (q.qu[0].raw() - r[0][0]*pu[-1].raw())/r[0][1];
            ppr[0].raw()= (q.qpr[0].raw() - r[0][0]*ppr[-1].raw())/r[0][1];
            // Compute b22
            b[0][0]= b[-1][1]; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation( b[0][0], b[0][1], c[0], s[0]);
            break;
          default:
            r[0][0]= 0.; r[0][1]= q.b[-1]; r[0][2]= q.a0;
            double tmp= q.b[0];
            GMRES_ApplyPlaneRotation( r[0][0], r[0][1], c[-2], s[-2]);
            GMRES_ApplyPlaneRotation( r[0][1], r[0][2], c[-1], s[-1]);
            GMRES_GeneratePlaneRotation( r[0][2], tmp, c[0], s[0]);
            GMRES_ApplyPlaneRotation( r[0][2], tmp, c[0], s[0]);
            // p[0]= (q.q[0] - r[0][0]*p[-2] -r[0][1]*p[-1])/r[0][2];
            pu[0].raw()= (q.qu[0].raw() - r[0][0]*pu[-2].raw() -r[0][1]*pu[-1].raw())*(1./r[0][2]);
            ppr[0].raw()= (q.qpr[0].raw() - r[0][0]*ppr[-2].raw() -r[0][1]*ppr[-1].raw())*(1./r[0][2]);
            b[0][0]= b[-1][1]; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation( b[0][0], b[0][1], c[0], s[0]);
        }
        dxu.raw()= (norm_r0*b[0][0])*pu[0].raw();
        dxpr.raw()= (norm_r0*b[0][0])*ppr[0].raw();
        u.raw()+= dxu.raw();
        pr.raw()+= dxpr.raw();

        // This is for fair comparisons of different solvers:
//        err= std::sqrt( (rhsu - (A*u + transp_mul( B, pr))).norm2() + (rhspr - B*u).norm2());
        res= std::fabs( norm_r0*b[0][1]);
//        std::cerr << "PMINRES: residual: " << res << '\t' << " 2-residual: " << err << '\n';
        if (res<=tol || lucky==true) {
            tol= res;
            max_iter= k;
            return true;
        }
        q.next();
        if (q.breakdown()) {
            lucky= true;
            std::cerr << "MINRES: lucky breakdown\n";
        }
        c.rotate(); s.rotate(); r.rotate(); pu.rotate(); ppr.rotate(); b.rotate();
    }
    tol= res;
    return false;
}


// Preconditioned MINRES solver for the Stokes-equations.
template <typename Lanczos>
class PMResSPCL : public SolverBaseCL
{
  private:
    Lanczos* q_;

  public:
    PMResSPCL(Lanczos& q, int maxiter, double tol)
      :SolverBaseCL( maxiter,tol), q_( &q)
    {}

    Lanczos*&       GetONB ()       { return q_; }
    const Lanczos*& GetONB () const { return q_; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c)
    {
        q_->new_basis( A, B, b - (A*v + transp_mul( B, p)), c - B*v);
        _res=  _tol;
        _iter= _maxiter;
        PMINRES_SP( A, B, v, p, b, c, *q_, _iter, _res);
	std::cerr << "PMResSPCL::Solve: iterations: " << GetIter() << "\tresidual: " << GetResid() << '\n';
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c, int& numIter, double& resid) const
    {
        q_->new_basis( A, B, b - (A*v + transp_mul( B, p)), c - B*v);
        PMINRES_SP( A, B, v, p, b, c, *q_, numIter, resid);
    }
};


// The identity-preconditioner for PMinres.
// For testing. If you really want plain Minres, use
// LanczosONB_SPCL instead of PLanczosONB_SPCL<..., IdPreCL>;
// it is about 4% faster.
class IdPreCL
{
  private:

  public:
    template <typename Mat, typename Vec>
    void
    Apply(const Mat& /*A*/, const Mat& /*B*/, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
        v= b; p= c;
    }
};


class DiagPCGPreCL
{
  private:
    mutable PCG_SsorCL PA_; // Preconditioner for A.
    mutable PCG_SsorCL PS_; // Preconditioner for S.
    const MatrixCL&   M_; // Preconditioner for S.

  public:
    DiagPCGPreCL( const MatrixCL& M)
      :PA_( SSORPcCL(), 8, 1e-20), PS_( SSORPcCL(), 3, 1e-20), M_( M) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
//        PA_.SetMaxIter( 500); PA_.SetTol( (b - A*v).norm()*1e-4);
        PA_.Solve( A, v, b);
//        std::cerr << PA_.GetIter() << '\t' << PA_.GetResid() << '\n';
        PS_.Solve( M_, p, c);
    }
};


class DiagMGPreCL
{
  private:
    const MGDataCL& A_; // Preconditioner for A.
    const MatrixCL& M_; // Preconditioner for S.
    Uint iter_vel_;

  public:
    DiagMGPreCL(const MGDataCL& A, const MatrixCL& M, Uint iter_vel)
      :A_( A), M_( M), iter_vel_( iter_vel) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
//        PA_.SetMaxIter( 1); PA_.SetTol( (bb - K.A_*u).norm()*1e-4);
        Uint   sm   =  2; // how many smoothing steps?
        int    lvl  = -1; // how many levels? (-1=all)
        double omega= 1.; // relaxation parameter for smoother
        SORsmoothCL smoother( omega);  // Gauss-Seidel with over-relaxation
        SSORPcCL P1;
        PCG_SsorCL solver( P1, 200, 1e-12);
        for (DROPS::Uint i=0; i<iter_vel_; ++i)
            MGM( A_.begin(), --A_.end(), v, b, smoother, sm, solver, lvl, -1);
        P1.Apply( M_, p, c);
    }
};


//=============================================================================
// Preconditioner for Minres-like solvers: SSOR-steps for A, ISPreCL for S. 
//=============================================================================
class Minres_SSOR_IS_PreCL
{
  private:
    SSORPcCL Apc_; // Preconditioner for A.
    const ISPreCL& Spc_; // Preconditioner for S.
    Uint iter_vel_;

  public:
    Minres_SSOR_IS_PreCL(const ISPreCL& Spc, Uint iter_vel= 1)
      :Apc_( 1.0), Spc_( Spc), iter_vel_( iter_vel) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, const Mat& B, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
        for (Uint i=0; i<iter_vel_; ++i)
	    Apc_.Apply( A, v, b);
	Spc_.Apply( B, p, c); // B is just a dummy.
    }
};


//=============================================================================
//  Derived classes for easier use
//=============================================================================

class Uzawa_CG_CL : public UzawaSolverCL<CGSolverCL>
{
  private:
    CGSolverCL _CGsolver;
  public:
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


/*
// Reimplementation of PSchur_PCG_CL with PSchurSolver2CL;
// this is only a check for the new class.
class PSchur2_PCG_CL: public PSchurSolver2CL<PCG_SsorCL,
                                             PCGSolverCL< PreGSOwnMatCL<P_SSOR0> > >
{
  private:
    PCG_SsorCL PCGsolver_;
    PCGSolverCL< PreGSOwnMatCL<P_SSOR0> > PCGsolver2_;

  public:
    PSchur2_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol)
        : PSchurSolver2CL<PCG_SsorCL, PCGSolverCL< PreGSOwnMatCL<P_SSOR0> > >( PCGsolver_, PCGsolver2_, outer_iter, outer_tol),
          PCGsolver_( SSORPcCL( 1.), inner_iter, inner_tol),
          PCGsolver2_( PreGSOwnMatCL<P_SSOR0>( M), outer_iter, outer_tol)
        {}
};
*/

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


class MinresSPCL : public PMResSPCL<LanczosONB_SPCL<MatrixCL, VectorCL> >
{
  private:
    LanczosONB_SPCL<MatrixCL, VectorCL> q_;

  public:
    MinresSPCL(int maxiter, double tol)
        :PMResSPCL<LanczosONB_SPCL<MatrixCL, VectorCL> >( q_, maxiter, tol),
         q_()
    {}
};

class PMinresSP_DiagPCG_CL : public PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, DiagPCGPreCL> >
{
  private:
    DiagPCGPreCL pre_;
    PLanczosONB_SPCL<MatrixCL, VectorCL, DiagPCGPreCL> q_;

  public:
    PMinresSP_DiagPCG_CL(const MatrixCL& M, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, DiagPCGPreCL> >( q_, maxiter, tol),
         pre_( M), q_( pre_)
    {}
};

class PMinresSP_DiagMG_CL : public PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, DiagMGPreCL> >
{
  private:
    DiagMGPreCL pre_;
    PLanczosONB_SPCL<MatrixCL, VectorCL, DiagMGPreCL> q_;

  public:
    PMinresSP_DiagMG_CL(const MGDataCL& A, const MatrixCL& M, int iter_vel, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, DiagMGPreCL> >( q_, maxiter, tol),
         pre_( A, M, iter_vel), q_( pre_)
    {}
};

//=============================================================================
// PMinres solver for the instationary Stokes-Equations without MG.
//=============================================================================
class PMinresSP_Diag_CL: public PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, Minres_SSOR_IS_PreCL> >
{
  private:
    Minres_SSOR_IS_PreCL pre_;
    PLanczosONB_SPCL<MatrixCL, VectorCL, Minres_SSOR_IS_PreCL> q_;

  public:
    PMinresSP_Diag_CL(const ISPreCL& Spc, int iter_vel, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, Minres_SSOR_IS_PreCL> >( q_, maxiter, tol),
         pre_( Spc, iter_vel), q_( pre_)
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
    static VectorCL x( M.A_.num_cols());
    if (x.size() != M.A_.num_cols())
    {
//        std::cerr << "> vector resized: old size was " << x.size();
        x.resize( M.A_.num_cols() );
//        std::cerr << ", new size is " << x.size() << '\n';
    }
    
    M.solver_.Solve( M.A_, x, transp_mul( M.B_, v));
//    std::cerr << "> inner iterations: " << M.solver_.GetIter()
//              << "\tresidual: " << M.solver_.GetResid() << std::endl;
    return M.B_*x;
}


//=============================================================================
// ApproximateSchurComplMatrixCL
// BApc^{-1}B^T, where Apc is a preconditioner for A.
//=============================================================================
template<typename>
class ApproximateSchurComplMatrixCL;

template<typename T>
VectorCL operator*(const ApproximateSchurComplMatrixCL<T>&, const VectorCL&);

template<class APC>
class ApproximateSchurComplMatrixCL
{
  private:
    const MatrixCL& A_;
    APC& Apc_;
    const MatrixCL& B_;

  public:
    ApproximateSchurComplMatrixCL(const MatrixCL& A, APC& Apc, const MatrixCL& B)
        : A_( A), Apc_( Apc), B_( B) {}

    friend VectorCL
    operator*<>(const ApproximateSchurComplMatrixCL<APC>&, const VectorCL&);
};

template<class APC>
VectorCL operator*(const ApproximateSchurComplMatrixCL<APC>& M, const VectorCL& v)
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
    if (_tmp.size() != v.size())
        _tmp.resize( v.size());
    _poissonSolver.Solve( A, _tmp, b);
    std::cerr << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;
    rhs+= B*_tmp;

    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;
    double tol= _tol;
    CG( SchurComplMatrixCL<PoissonSolverT>( _poissonSolver, A, B), p, rhs, iter, tol);
    std::cerr << "iterations: " << iter << "\tresidual: " << tol << std::endl;
    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;

    _poissonSolver.Solve( A, v, b - transp_mul(B, p));
    std::cerr << "Iterationen: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;

    _iter= iter+_poissonSolver.GetIter();
    _res= tol + _poissonSolver.GetResid();
/* std::cerr << "Real residuals are: "
          << (A*v+transp_mul(B, p)-b).norm() << ", "
          << (B*v-c).norm() << std::endl; */
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
    if (_tmp.size() != v.size())
        _tmp.resize( v.size());
    _poissonSolver.Solve( A, _tmp, b);
    std::cerr << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;
    rhs+= B*_tmp;

    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= _maxiter;
    double tol= _tol;
    PCG( SchurComplMatrixCL<PoissonSolverT>( _poissonSolver, A, B), p, rhs, _schurPc, iter, tol);
    std::cerr << "iterations: " << iter << "\tresidual: " << tol << std::endl;
    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;

    _poissonSolver.Solve( A, v, b - transp_mul(B, p));
    std::cerr << "iterations: " << _poissonSolver.GetIter()
              << "\tresidual: " << _poissonSolver.GetResid() << std::endl;

    _iter= iter+_poissonSolver.GetIter();
    _res= std::sqrt( tol*tol + _poissonSolver.GetResid()*_poissonSolver.GetResid());
    std::cerr << "-----------------------------------------------------" << std::endl;
}

template <typename InnerSolverT, typename OuterSolverT>
void PSchurSolver2CL<InnerSolverT, OuterSolverT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
// solve:       S*p = B*(A^-1)*b - c   with SchurCompl. S = B A^(-1) BT
//              A*u = b - BT*p
{
    VectorCL rhs= -c;
    if (tmp_.size() != v.size()) tmp_.resize( v.size());
    innerSolver_.Solve( A, tmp_, b);
    std::cerr << "rhs     : iterations: " << innerSolver_.GetIter()
              << "\tresidual: " << innerSolver_.GetResid() << std::endl;
    rhs+= B*tmp_;

    outerSolver_.SetTol( _tol);
    outerSolver_.SetMaxIter( _maxiter);
    outerSolver_.Solve( SchurComplMatrixCL<InnerSolverT>( innerSolver_, A, B), p, rhs);
    std::cerr << "pressure: iterations: " << outerSolver_.GetIter()
              << "\tresidual: " << outerSolver_.GetResid() << std::endl;

    innerSolver_.Solve( A, v, b - transp_mul(B, p));
    std::cerr << "velocity: iterations: " << innerSolver_.GetIter()
              << "\tresidual: " << innerSolver_.GetResid() << std::endl;

    _iter= innerSolver_.GetIter() + outerSolver_.GetIter();
    _res= std::sqrt( std::pow( innerSolver_.GetResid(), 2)
                   + std::pow( outerSolver_.GetResid(), 2));
    std::cerr << "-----------------------------------------------------" << std::endl;
}


//=============================================================================
// Inexact Uzawa-method from "Fast Iterative Solvers for Discrete Stokes
// Equations", Peters, Reichelt, Reusken, Chapter 3.3.
// The preconditioner Apc for A must be "good" (MG-like) to guarantee
// convergence.
//=============================================================================
template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
InexactUzawa(const Mat& A, const Mat& B, Vec& xu, Vec& xp, const Vec& f, const Vec& g,
    PC1& Apc, PC2& Spc,
    int& max_iter, double& tol)
{
    VectorCL ru= f - A*xu - transp_mul( B, xp);
    VectorCL w( f.size());
    VectorCL z( g.size());
    VectorCL a( f.size());
    VectorCL xuneu( f.size());
    VectorCL b( f.size());
    ApproximateSchurComplMatrixCL<PC1> asc( A, Apc, B);
    PCGSolverCL<PC2> pcgsolver( Spc, 100, 0.5);
    const double resid0= std::sqrt( ru.norm2() + (g - B*xu).norm2());
    double resid= 0.0;
    std::cerr << "residual (2-norm): " << resid0 << '\n';
    if (resid0<=tol) { // The fixed point iteration between levelset and Stokes
        tol= resid;    // equation uses this to determine convergence.
        max_iter= 0;
        return true;
    }
    for (int k= 1; k<=max_iter; ++k) {
        w= 0.0;
        Apc.Apply( A, w, ru);
        w+= xu;
        z= 0.0;
//        std::cerr << "B*w : " << (B*w).norm() << "\tprojektion auf 1: "
//                  << (B*w)*VectorCL( 1.0/std::sqrt( (double)g.size()), g.size())
//                  << "\n";
        // Due to theory (see paper) we must use relative error of about 0.5. z==0.
	pcgsolver.SetTol( 0.5*(B*w - g).norm());
        pcgsolver.Solve( asc, z, B*w - g);
//        std::cerr << "pcgsolver: iterations: " << pcgsolver.GetIter() 
//                  << "\tresid: " << pcgsolver.GetResid() << '\n';

        b= transp_mul( B, z);
        a= 0.0;
        Apc.Apply( A, a, b);
        z_xpay( xuneu, w, -1.0, a); // xuneu= w - a;
        xp+= z;
        z_xpaypby2(ru, ru, -1.0, A*(xuneu - xu), -1.0, transp_mul( B, z)); // ru-= A*(xuneu - xu) + transp_mul( B, z);
        xu= xuneu;
        resid= std::sqrt( (f - A*xu - transp_mul( B, xp)).norm2() + (g - B*xu).norm2());
//        std::cerr << "relative residual (2-norm): " << resid/resid0 
//                  << "\tv: " << (f - A*xu - transp_mul( B, xp)).norm()
//                  << "\tp: " << (g - B*xu).norm()
//                  << '\n';
/*
        if (resid<=tol*resid0) { // relative errors
            tol= resid0==0.0 ? 0.0 : resid/resid0;
            max_iter= k;
            return true;
        }
*/
        if (resid<=tol) { // absolute errors
            std::cerr << "relative residual (2-norm): " << resid/resid0 
                      << "\tv: " << (f - A*xu - transp_mul( B, xp)).norm()
                      << "\tp: " << (g - B*xu).norm()
                      << '\n';
            tol= resid;
            max_iter= k;
            return true;
        }

    }
    tol= resid;
    return false;
}


} // end of namespace DROPS

#endif
