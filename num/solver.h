//**************************************************************************
// File:    solver.h                                                       *
// Content: iterative solvers                                              *
// Authors: Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************


#ifndef DROPS_SOLVER_H
#define DROPS_SOLVER_H

#include <vector>
#include "misc/container.h"
#include "num/spmat.h"

namespace DROPS
{

//*****************************************************************************
//
//  Gauss-Seidel type methods for preconditioning
//
//*****************************************************************************

//=============================================================================
//  Template magic for the selection of the preconditioner
//=============================================================================

// Available methods
enum PreMethGS
{
    P_JAC,     // Jacobi
    P_GS,      // Gauss-Seidel
    P_SGS,     // symmetric Gauss-Seidel
    P_SGS0,    // symmetric Gauss-Seidel with initial vector 0
    P_JOR,     // P_JAC with over-relaxation
    P_SOR,     // P_GS with over-relaxation
    P_SSOR,    // P_SGS with over-relaxation
    P_SSOR0,   // P_SGS0 with over-relaxation
    P_SGS0_D,  // P_SGS0 using SparseMatDiagCL
    P_SSOR0_D, // P_SSOR0 using SparseMatDiagCL
    P_DUMMY    // identity
};

// Base methods
enum PreBaseGS { PB_JAC, PB_GS, PB_SGS, PB_SGS0, PB_DUMMY };

// Properties of the methods
template <PreMethGS PM> struct PreTraitsCL
{
    static const PreBaseGS BaseMeth= PreBaseGS(PM<8 ? PM%4 : (PM==10 ? PB_DUMMY : PB_SGS0) );
    static const bool      HasOmega= (PM>=4 && PM<8) || PM==9;
    static const bool      HasDiag=  PM==8 || PM==9;
};

// Used to make a distinct type from each method
template <PreBaseGS> class PreDummyCL {};


//=============================================================================
//  Implementation of the methods
//=============================================================================

// One step of the Jacobi method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_JAC>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    Vec          y(x.size());

    for (size_t i=0, nz=0; i<n; ++i)
    {
        double aii, sum= b[i];
        for (const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        if (HasOmega)
            y[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            y[i]= sum/aii;
    }

    std::swap(x,y);
}


// One step of the Gauss-Seidel/SOR method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_GS>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    double aii, sum;

    for (size_t i=0, nz=0; i<n; ++i) {
        sum= b[i];
        const size_t end= A.row_beg( i+1);
        for (; A.col_ind( nz) != i; ++nz) // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
            sum-= A.val( nz)*x[A.col_ind( nz)];
        aii= A.val( nz++);
        for (; nz<end; ++nz)
            sum-= A.val( nz)*x[A.col_ind( nz)];
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_SGS>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    double aii, sum;

    for (size_t i=0, nz=0; i<n; ++i) {
        sum= b[i];
        const size_t end= A.row_beg( i+1);
        for (; A.col_ind( nz) != i; ++nz) // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
            sum-= A.val( nz)*x[A.col_ind( nz)];
        aii= A.val( nz++);
        for (; nz<end; ++nz)
            sum-= A.val( nz)*x[A.col_ind( nz)];
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
    for (size_t i= n, nz= A.row_beg( n); i>0; ) { // XXX: Rearrange this loop as the preceding one.
        --i;
        double aii, sum= b[i];
        for (const size_t beg= A.row_beg(i); nz>beg; ) {
            --nz;
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        }
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector 0
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_SGS0>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        double sum= b[i];
        size_t j= A.row_beg(i);
        for ( ; A.col_ind(j) < i; ++j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= omega*sum/A.val(j);
        else
            x[i]= sum/A.val(j);
    }

    for (size_t i=n; i>0; )
    {
        --i;
        double sum= 0;
        size_t j= A.row_beg(i+1)-1;
        for ( ; A.col_ind(j) > i; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/A.val(j);
        else
            x[i]+= sum/A.val(j);
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector 0,
// uses SparseMatDiagCL for the location of the diagonal
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_SGS0>&, const MatrixCL& A, Vec& x, const Vec& b, const SparseMatDiagCL& diag, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        double sum= b[i];
        for (size_t j= A.row_beg(i); j < diag[i]; ++j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= omega*sum/A.val(diag[i]);
        else
            x[i]= sum/A.val(diag[i]);
    }

    for (size_t i=n; i>0; )
    {
        --i;
        double sum= 0;
        for (size_t j= A.row_beg(i+1)-1; j > diag[i]; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/A.val(diag[i]);
        else
            x[i]+= sum/A.val(diag[i]);
    }
}


template <bool HasOmega, typename Vec, typename Real>
void
SolveGSstep(const PreDummyCL<PB_DUMMY>&, const SparseMatBaseCL<Real>&, Vec& x, const Vec& b, double)
{
    x=b;
}

//=============================================================================
//  Preconditioner classes
//=============================================================================

// TODO: Referenzen benutzen statt Daten kopieren. Init ueberdenken.

// Preconditioners without own matrix
template <PreMethGS PM, bool HasDiag= PreTraitsCL<PM>::HasDiag> class PreGSCL;

// Simple preconditioners
template <PreMethGS PM>
class PreGSCL<PM,false>
{
  private:
    double _omega;

  public:
    PreGSCL (double om= 1.0) : _omega(om) {}

    template <typename Vec>
    void Apply(const MatrixCL& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A, x, b, _omega);
    }
};


// Preconditioner with SparseMatDiagCL
template <PreMethGS PM>
class PreGSCL<PM,true>
{
  private:
    const SparseMatDiagCL* _diag;
    double                 _omega;

  public:
    PreGSCL (double om= 1.0) : _diag(0), _omega(om) {}
    PreGSCL (const PreGSCL& p) : _diag(p._diag ? new SparseMatDiagCL(*(p._diag)) : 0), _omega(p._omega) {}
    ~PreGSCL() { delete _diag; }

    void Init(const MatrixCL& A)
    {
        delete _diag; _diag=new SparseMatDiagCL(A);
    }

    template <typename Vec>
    void Apply(const MatrixCL& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A, x, b, *_diag, _omega);
    }
};


// Preconditioners with own matrix
template <PreMethGS PM, bool HasDiag= PreTraitsCL<PM>::HasDiag>
class PreGSOwnMatCL;

// Simple preconditioners
template <PreMethGS PM>
class PreGSOwnMatCL<PM,false>
{
  private:
    const MatrixCL& _M;
    double          _omega;

  public:
    PreGSOwnMatCL (const MatrixCL& M, double om= 1.0) : _M(M), _omega(om) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), _M, x, b, _omega);
    }
};


// Preconditioner with SparseMatDiagCL
template <PreMethGS PM>
class PreGSOwnMatCL<PM,true>
{
  private:
    const MatrixCL&        _M;
    const SparseMatDiagCL* _diag;
    double                 _omega;

  public:
    PreGSOwnMatCL (const MatrixCL& M, double om= 1.0) : _M(M), _diag(0), _omega(om) {}
    PreGSOwnMatCL (const PreGSOwnMatCL&); // not defined
    ~PreGSOwnMatCL() { delete _diag; }

    void Init(const MatrixCL& A)
    {
        delete _diag; _diag=new SparseMatDiagCL(A);
    }

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), _M, x, b, *_diag, _omega);
    }
};


//*****************************************************************************
//
//  Conjugate gradients (CG, PCG) and GMRES
//
//*****************************************************************************

//=============================================================================
//  Implementation of the methods
//=============================================================================

//-----------------------------------------------------------------------------
// CG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//
//        x - approximate solution to Ax = b
// max_iter - number of iterations performed before tolerance was reached
//      tol - residual after the final iteration
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec>
bool
CG(const Mat& A, Vec& x, const Vec& b, int& max_iter, double& tol)
{
    Vec r= A*x - b;
    Vec d= -r;
    double resid= r.norm2();

    tol*= tol;

    if (resid<=tol)
    {
        tol= sqrt(resid);
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        const Vec    Ad= A*d;
        const double delta= Ad*d;
        const double alpha= resid/delta;
        double       beta= resid;

        axpy(alpha, d, x);  // x+= alpha*d;
        axpy(alpha, Ad, r); // r+= alpha*Ad;

        resid= r.norm2();
        if (resid<=tol)
        {
            tol= sqrt(resid);
            max_iter= i;
            return true;
        }
        beta= resid / beta;
        d= beta*d-r;
    }
    tol= sqrt(resid);
    return false;
}


//-----------------------------------------------------------------------------
// PCG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// Upon successful return, output arguments have the following values:
//
//        x - approximate solution to Ax = b
// max_iter - number of iterations performed before tolerance was reached
//      tol - residual after the final iteration
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PreCon>
bool
PCG(const Mat& A, Vec& x, const Vec& b, const PreCon& M,
    int& max_iter, double& tol)
{
    const size_t n= x.size();
    Vec p(n), z(n), q(n), r= b - A*x;
    double rho, rho_1= 0.0, resid= r.norm2();

    tol*= tol;

    if (resid<=tol)
    {
        tol= sqrt(resid);
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        M.Apply(A, z, r);
        rho= r*z;
        if (i == 1)
            p= z;
        else
            z_xpay(p, z, (rho/rho_1), p); // p= z + (rho/rho_1)*p;

        q= A*p;
        const double alpha= rho/(p*q);
        axpy(alpha, p, x);                // x+= alpha*p;
        axpy(-alpha, q, r);               // r-= alpha*q;

        resid= r.norm2();
        if (resid<=tol)
        {
            tol= sqrt(resid);
            max_iter= i;
            return true;
        }
        rho_1= rho;
    }
    tol= sqrt(resid);
    return false;
}


//-----------------------------------------------------------------------------
// GMRES:
//-----------------------------------------------------------------------------

inline void GMRES_GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
    if (dy == 0.0)
    {
        cs = 1.0;
        sn = 0.0;
    }
    else
    {
        const double r=std::sqrt(dx*dx+dy*dy);
        cs = dx/r;
        sn = dy/r;
    }
}


inline void GMRES_ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
    const double tmp = cs*dx + sn*dy;
    dy = -sn*dx + cs*dy;
    dx = tmp;
}


template <typename Mat, typename Vec>
void GMRES_Update(Vec &x, int k, const Mat &H, const Vec &s, const std::vector<Vec> &v)
{
    Vec y(s);

    // Backsolve:
    for ( int i=k; i>=0; --i )
    {
        y[i] /= H(i,i);
        for ( int j=i-1; j>=0; --j )
            y[j] -= H(j,i) * y[i];
    }

    for ( int i=0; i<=k; ++i )
        x += y[i] * v[i];
}


template <typename Mat, typename Vec, typename PreCon>
bool
GMRES(const Mat& A, Vec& x, const Vec& b, const PreCon& M,
      int m, int& max_iter, double& tol)
{
    DMatrixCL<double> H(m,m);
    Vec               s(m), cs(m), sn(m), w(b.size()), r(b.size());
    std::vector<Vec>  v(m);
    double            beta, normb, resid;

    for ( int i=0; i<m; ++i )
        v[i].resize(b.size());

    M.Apply(A, r, b-A*x);
    beta = r.norm();

    M.Apply(A, w, b);
    normb=w.norm();
    if (normb == 0.0) normb=1;

    resid = beta/normb;
    if (resid<=tol)
    {
        tol = resid;
        max_iter = 0;
        return true;
    }

    int j=1;
    while (j <= max_iter)
    {
        v[0] = r * (1.0 / beta);
        s = 0.0;
        s[0] = beta;

        for ( int i=0; i<m-1 && j<=max_iter; ++i, ++j )
        {
            M.Apply(A, w, A*v[i]);
            for ( int k=0; k<=i; ++k )
            {
                H(k, i) = w * v[k];
                w -= H(k, i) * v[k];
            }
            H(i+1, i) = w.norm();
            v[i+1] = w * (1.0 / H(i+1,i));

            for ( int k=0; k<i; ++k )
                GMRES_ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);

            GMRES_GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
            GMRES_ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
            GMRES_ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

            resid = std::abs(s[i+1])/normb;
            if (resid<=tol)
            {
                GMRES_Update(x, i, H, s, v);
                tol = resid;
                max_iter = j;
                return true;
            }
        }

        GMRES_Update(x, m - 2, H, s, v);
        M.Apply(A, r, b-A*x);
        beta = r.norm();
        resid = beta/normb;
        if (resid<=tol)
        {
            tol = resid;
            max_iter = j;
            return true;
        }
    }

    tol = resid;
    return false;
}


// One recursive step of Lanzcos' algorithm for computing an ONB (q1, q2, q3,...)
// of the Krylovspace of A for a given starting vector r. This is a three term
// recursion, computing the next q_i from the two previous ones.
// See Arnold Reusken, "Numerical methods for elliptic partial differential equations",
// p. 148.
// Returns false for 'lucky breakdown' (see below), true in the generic case.
template <typename Mat, typename Vec>
bool
LanczosStep(const Mat& A,
            const Vec& q0, const Vec& q1, Vec& q2,
            double& a1,
            const double b0, double& b1)
{
    q2.raw()= (A*q1).raw() - b0*q0.raw();
    a1= q2*q1;
    q2.raw()-= a1*q1.raw();
    b1= q2.norm();
    // Lucky breakdown; the Krylov-space K up to q1 is A-invariant. Thus,
    // the correction dx needed to solve A(x0+dx)=b is in this space and
    // the Minres-algo will terminate with the exact solution in the
    // following step.
    if (std::fabs(b1) < 1e-15) return false;
    q2/= b1;
    return true;
}

template <typename Mat, typename Vec>
class LanczosONBCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;

  public:
    const Mat& A;
    SBufferCL<Vec, 3> q;
    double a0;
    SBufferCL<double, 2> b;

    // Sets up initial values and computes q0.
    LanczosONBCL(const Mat& A_, const Vec& r0)
      :A( A_) {
        q[-1].resize( r0.size(), 0.);
        norm_r0_= r0.norm();
        q[0].resize( r0.size(), 0.); q[0]= r0/norm_r0_;
        q[1].resize( r0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= LanczosStep( A, q[-1], q[0], q[1], a0, b[-1], b[0]);
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
        q.rotate(); b.rotate();
        return (nobreakdown_= LanczosStep( A, q[-1], q[0], q[1], a0, b[-1], b[0]));
    }
};


// One recursive step of the preconditioned Lanzcos algorithm for computing an ONB of
// a Krylovspace. This is a three term
// recursion, computing the next q_i from the two previous ones.
// See Arnold Reusken, "Numerical methods for elliptic partial differential equations",
// p. 153.
// Returns false for 'lucky breakdown' (see below), true in the generic case.
template <typename Mat, typename Vec, typename PreCon>
bool
PLanczosStep(const Mat& A,
             const PreCon& M,
             const Vec& q1, Vec& q2,
             const Vec& t0, const Vec& t1, Vec& t2,
             double& a1,
             const double b0, double& b1)
{
    t2.raw()= (A*q1).raw() - b0*t0.raw();
    a1= t2*q1;
    t2.raw()-= a1*t1.raw();
    M.Apply( A, q2, t2);
    b1= std::sqrt( q2*t2);
    if (b1 < 1e-15) return false;
    t2.raw()*= 1./b1;
    q2.raw()*= 1./b1;
    return true;
}

template <typename Mat, typename Vec, typename PreCon>
class PLanczosONBCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;

  public:
    const Mat& A;
    const PreCon& M;
    SBufferCL<Vec, 2> q;
    SBufferCL<Vec, 3> t;
    double a0;
    SBufferCL<double, 2> b;

    // Sets up initial values and computes q0.
    PLanczosONBCL(const Mat& A_, const PreCon& M_, const Vec& r0)
      :A( A_), M( M_) {
        t[-1].resize( r0.size(), 0.);
        q[-1].resize( r0.size(), 0.); M.Apply( A, q[-1], r0);
        norm_r0_= std::sqrt( q[-1]*r0);
        t[0].resize( r0.size(), 0.); t[0]= r0/norm_r0_;
        q[0].resize( r0.size(), 0.); q[0]= q[-1]/norm_r0_;
        t[1].resize( r0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= PLanczosStep( A, M, q[0], q[1], t[-1], t[0], t[1], a0, b[-1], b[0]);
    }

    double norm_r0() const {
        return norm_r0_; }
    bool
    breakdown() const {
        return !nobreakdown_; }
    // Computes new q_i, t_i, a_i, b_1, q_{i+1} in q0, t_0, a0, b0, q1 and moves old
    // values to qm1, tm1, bm1.
    bool
    next() {
        q.rotate(); t.rotate(); b.rotate();
        return (nobreakdown_= PLanczosStep( A, M, q[0], q[1], t[-1], t[0], t[1], a0, b[-1], b[0]));
    }
};

//-----------------------------------------------------------------------------
// PMINRES: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// See Arnold Reusken, "Numerical methods for elliptic partial differential
// equations", pp. 149 -- 154
//
// Upon successful return, output arguments have the following values:
//
//        x - approximate solution to Ax = rhs
// max_iter - number of iterations performed before tolerance was reached
//      tol - 2-norm of the last correction dx to x after the final iteration.
//-----------------------------------------------------------------------------
template <typename Mat, typename Vec, typename Lanczos>
bool
PMINRES(const Mat& A, Vec& x, const Vec& rhs, Lanczos& q, int& max_iter, double& tol)
{
    Vec dx= rhs - A*x; // First, the residual, later used as vector for
                       // updating x, thus the name.
    const double resid0= dx.norm();
    double err= resid0*resid0;

    tol*= tol;
//    if (err<=tol) {
//        tol= sqrt( err);
//        max_iter= 0;
//        return true;
//    }
    const double norm_r0= q.norm_r0();
    bool lucky= q.breakdown();
    SBufferCL<double, 3> c;
    SBufferCL<double, 3> s;
    SBufferCL<SVectorCL<3>, 3> r;
    SBufferCL<Vec, 3> p;
    p[0].resize( x.size()); p[1].resize( x.size()); p[2].resize( x.size());
    SBufferCL<SVectorCL<2>, 2> b;
    
    for (int k=1; k<=max_iter; ++k) {
        switch (k) {
          case 1:
            // Compute r1
            GMRES_GeneratePlaneRotation( q.a0, q.b[0], c[0], s[0]);
            r[0][0]= std::sqrt( q.a0*q.a0 + q.b[0]*q.b[0]);
            // Compute p1
            // p[0]= q.q[0]/r[0][0];
            p[0].raw()= q.q[0].raw()/r[0][0];
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
            p[0].raw()= (q.q[0].raw() - r[0][0]*p[-1].raw())/r[0][1];
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
            p[0].raw()= (q.q[0].raw() - r[0][0]*p[-2].raw() -r[0][1]*p[-1].raw())*(1/r[0][2]);
            b[0][0]= b[-1][1]; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation( b[0][0], b[0][1], c[0], s[0]);
        }
        dx.raw()= (norm_r0*b[0][0])*p[0].raw();
//std::cout << "q.q[0] " << q.q[0] << " q.a0 " << q.a0 << " q.b[0] " << q.b[0] << " c " << c[0] << " s " << s[0] << " r " << r[0] << " p " << p[0] << " b " << b[0]
//          << std::endl;
//std::cout << "dx\n" << dx;
        x.raw()+= dx.raw();
        err= dx.norm2();
//        err= (rhs - A*x).norm2();
        if (err<=tol || lucky==true) {
            tol= sqrt( err);
            max_iter= k;
            return true;
        }
        q.next();
        if (q.breakdown()) {
            lucky= true;
            Comment( "MINRES: lucky breakdown\n", ~0);
        }
        c.rotate(); s.rotate(); r.rotate(); p.rotate(); b.rotate();
    }
    tol= sqrt( err);
    return false;
}


template <typename Mat, typename Vec>
bool
MINRES(const Mat& A, Vec& x, const Vec& rhs, int& max_iter, double& tol)
{
    Vec dx= rhs - A*x;
    LanczosONBCL<Mat, Vec> q( A, dx);
    return PMINRES( A,  x, rhs, q, max_iter, tol);
}


//=============================================================================
//  Drivers
//=============================================================================

// What every iterative solver should have
class SolverBaseCL
{
  protected:
    int    _maxiter, _iter;
    double _tol, _res;

    SolverBaseCL (int maxiter, double tol)
        : _maxiter(maxiter), _iter(-1), _tol(tol), _res(-1.) {}

  public:
    void   SetTol    (double tol) { _tol= tol; }
    void   SetMaxIter(int iter)   { _maxiter= iter; }

    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }
};


// Bare CG solver
class CGSolverCL : public SolverBaseCL
{
  public:
    CGSolverCL(int maxiter, double tol) : SolverBaseCL(maxiter,tol) {}

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        CG(A, x, b, _iter, _res);
    }
    template <typename Vec>
    void Solve(const MatrixCL& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        CG(A, x, b, numIter, resid);
    }
};


// With preconditioner
template <typename PC>
class PCGSolverCL : public SolverBaseCL
{
  private:
    PC _pc;

  public:
    PCGSolverCL(const PC& pc, int maxiter, double tol)
        : SolverBaseCL(maxiter,tol), _pc(pc) {}

    PC&       GetPc ()       { return _pc; }
    const PC& GetPc () const { return _pc; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        PCG(A, x, b, _pc, _iter, _res);
    }
    template <typename Vec>
    void Solve(const MatrixCL& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        PCG(A, x, b, _pc, numIter, resid);
    }
};

// Bare MINRES solver
class MResSolverCL : public SolverBaseCL
{
  public:
    MResSolverCL(int maxiter, double tol) : SolverBaseCL( maxiter,tol) {}

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MINRES( A, x, b, _iter, _res);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        MINRES( A, x, b, numIter, resid);
    }
};

// Preconditioned MINRES solver
template <typename Lanczos>
class PMResSolverCL : public SolverBaseCL
{
  private:
    Lanczos& q_;

  public:
    PMResSolverCL(Lanczos& q, int maxiter, double tol)
      :SolverBaseCL( maxiter,tol), q_( q) {}

    Lanczos&       GetONB ()       { return q_; }
    const Lanczos& GetONB () const { return q_; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        PMINRES( A, x, b, q_, _iter, _res);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        PMINRES( A, x, b, q_, numIter, resid);
    }
};

// GMRES
template <typename PC>
class GMResSolverCL : public SolverBaseCL
{
  private:
    PC  pc_;
    int restart_;

  public:
    GMResSolverCL(const PC& pc, int restart, int maxiter, double tol)
        : SolverBaseCL(maxiter,tol), pc_(pc), restart_(restart) {}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }
    int       GetRestart () const { return restart_; }

    template <typename Vec>
    void Solve(const MatrixCL& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        GMRES(A, x, b, pc_, restart_, _iter, _res);
    }
    template <typename Vec>
    void Solve(const MatrixCL& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        GMRES(A, x, b, pc_, restart_, numIter, resid);
    }
};


//=============================================================================
//  Typedefs
//=============================================================================

typedef PreGSCL<P_SSOR>    SSORsmoothCL;
typedef PreGSCL<P_SOR>     SORsmoothCL;
typedef PreGSCL<P_SGS>     SGSsmoothCL;
typedef PreGSCL<P_JOR>     JORsmoothCL;
typedef PreGSCL<P_GS>      GSsmoothCL;
typedef PreGSCL<P_SGS0>    SGSPcCL;
typedef PreGSCL<P_SSOR0>   SSORPcCL;
typedef PreGSCL<P_SSOR0_D> SSORDiagPcCL;
typedef PreGSCL<P_DUMMY>   DummyPcCL;

typedef PCGSolverCL<SGSPcCL>      PCG_SgsCL;
typedef PCGSolverCL<SSORPcCL>     PCG_SsorCL;
typedef PCGSolverCL<SSORDiagPcCL> PCG_SsorDiagCL;

//*****************************************************************************
//
//  Non-iterative methods: Gauss solver with pivoting
//
//*****************************************************************************

template<size_t _Dim>
void
gauss_pivot(SMatrixCL<_Dim, _Dim>& A, SVectorCL<_Dim>& b)
{
    double max;
    size_t ind_max;
    size_t p[_Dim];
    SVectorCL<_Dim> b2= b;
    for (size_t i=0; i<_Dim; ++i) p[i]= i;

    for (size_t i=0; i<_Dim-1; ++i)
    {
        max= fabs(A(p[i], i));
        ind_max= i;
        for (size_t l=i+1; l<_Dim; ++l)
            if (fabs(A(p[l], i))>max)
            {
                max= fabs(A(p[l], i));
                ind_max= l;
            }
        if (max == 0.0) throw DROPSErrCL("gauss_pivot: Matrix is singular.");
        if (i!=ind_max) std::swap(p[i], p[ind_max]);
        const double piv= A(p[i], i);
        for (size_t j=i+1; j<_Dim; ++j)
        {
            const double fac= A(p[j], i)/piv;
            b2[p[j]]-= fac*b2[p[i]];
            for (size_t k=i+1; k<_Dim; ++k)
                A(p[j], k)-= fac* A(p[i], k);
        }
    }

    for (int i=_Dim-1; i>=0; --i)
    {
        for (int j=_Dim-1; j>i; --j)
            b2[p[i]]-= A(p[i], j)*b2[p[j]];
        b2[p[i]]/= A(p[i], i);
        b[i]= b2[p[i]];
    }
}

} // end of namespace DROPS

#endif
