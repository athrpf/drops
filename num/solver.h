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
    P_JAC,     //  0 Jacobi
    P_GS,      //  1 Gauss-Seidel
    P_SGS,     //  2 symmetric Gauss-Seidel
    P_SGS0,    //  3 symmetric Gauss-Seidel with initial vector 0
    P_JOR,     //  4 P_JAC with over-relaxation
    P_SOR,     //  5 P_GS with over-relaxation
    P_SSOR,    //  6 P_SGS with over-relaxation
    P_SSOR0,   //  7 P_SGS0 with over-relaxation
    P_SGS0_D,  //  8 P_SGS0 using SparseMatDiagCL
    P_SSOR0_D, //  9 P_SSOR0 using SparseMatDiagCL
    P_DUMMY,   // 10 identity
    P_GS0,     // 11 Gauss-Seidel with initial vector 0
    P_JAC0     // 12 Jacobi with initial vector 0
};

// Base methods
enum PreBaseGS { PB_JAC, PB_GS, PB_SGS, PB_SGS0, PB_DUMMY, PB_GS0, PB_JAC0 };

// Properties of the methods
template <PreMethGS PM> struct PreTraitsCL
{
    static const PreBaseGS BaseMeth= PreBaseGS(PM<8 ? PM%4
                                     : (PM==10 ? PB_DUMMY : (PM==11 ? PB_GS0 : (PM==12 ? PB_JAC0 : PB_SGS0))) );
    static const bool      HasOmega= (PM>=4 && PM<8) || PM==9 || PM==11 || PM==12;
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

// One step of the Jacobi method with start vector 0
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_JAC0>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    size_t nz;

    for (size_t i= 0; i < n; ++i) {
        nz= A.row_beg( i);
        for (const size_t end= A.row_beg( i+1); A.col_ind( nz) != i && nz < end; ++nz); // empty loop
        if (HasOmega)
            x[i]= /*(1.-omega)*x[i]=0  + */ omega*b[i]/A.val( nz);
        else
            x[i]= b[i]/A.val( nz);
    }
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

// One step of the Gauss-Seidel/SOR method with start vector x
template <bool HasOmega, typename Vec>
void
SolveGSstep(const PreDummyCL<PB_GS0>&, const MatrixCL& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    double aii, sum;

    for (size_t i=0, nz=0; i<n; ++i) {
        sum= b[i];
        const size_t end= A.row_beg( i+1);
        for (; A.col_ind( nz) != i; ++nz) // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
            sum-= A.val( nz)*x[A.col_ind( nz)];
        aii= A.val( nz);
        nz= end;
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
    for (size_t i= n, nz= A.row_beg( n); i>0; ) { // This is safe: Without diagonal entry, Gauss-Seidel would explode anyway.
        --i;
        double aii, sum= b[i];
        const size_t beg= A.row_beg( i);
        for (; A.col_ind( --nz) != i; ) {
            sum-= A.val( nz)*x[A.col_ind( nz)];
        }
        aii= A.val( nz);
        for (; nz>beg; ) {
            --nz;
            sum-= A.val( nz)*x[A.col_ind( nz)];
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


template <bool HasOmega, typename Vec, typename Mat>
void
SolveGSstep(const PreDummyCL<PB_DUMMY>&, const Mat&, Vec& x, const Vec& b, double)
{
    x= b;
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

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec& x, const Vec& b) const
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

class MultiSSORPcCL
// do multiple SSOR-steps
{
  private:
    double _omega;
    int   _num;

  public:
    MultiSSORPcCL(double om= 1.0, int num= 1) : _omega(om), _num(num) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec& x, const Vec& b) const
    {
        // one SSOR0-step
        SolveGSstep<PreTraitsCL<P_SSOR0>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<P_SSOR0>::BaseMeth>(), A, x, b, _omega);
        // _num-1 SSOR-steps
        for (int i=1; i<_num; ++i)
            SolveGSstep<PreTraitsCL<P_SSOR>::HasOmega,Vec>(PreDummyCL<PreTraitsCL<P_SSOR>::BaseMeth>(), A, x, b, _omega);
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
//      tol - (relative, see next parameter) 2-norm of the residual after the
//    final iteration
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol.
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec>
bool
CG(const Mat& A, Vec& x, const Vec& b, int& max_iter, double& tol,
    bool measure_relative_tol= false)
{
    Vec r( A*x - b);
    Vec d( -r);
    double normb= norm( b), res, resid= norm_sq( r);

    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    if ((res= std::sqrt( resid)/normb) <= tol)
    {
        tol= res;
        max_iter= 0;
        return true;
    }

    for (int i= 1; i <= max_iter; ++i)
    {
        const Vec    Ad= A*d;
        const double delta= dot( Ad, d);
        const double alpha= resid/delta;
        double       beta= resid;

        axpy(alpha, d, x);  // x+= alpha*d;
        axpy(alpha, Ad, r); // r+= alpha*Ad;

        resid= norm_sq( r);
        if ((res= std::sqrt( resid)/normb) <= tol)
        {
            tol= res;
            max_iter= i;
            return true;
        }
        beta= resid / beta;
        d= beta*d-r;
    }
    tol= res;
    return false;
}


//-----------------------------------------------------------------------------
// PCG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// Upon successful return, output arguments have the following values:
//
//        x - approximate solution to Ax = b
// max_iter - number of iterations performed before tolerance was reached
//      tol - 2-norm of the (relative, see below) residual after the final iteration
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol.
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PreCon>
bool
PCG(const Mat& A, Vec& x, const Vec& b, const PreCon& M,
    int& max_iter, double& tol, bool measure_relative_tol= false)
{
    const size_t n= x.size();
    Vec p(n), z(n), q(n), r( b - A*x);
    double rho, rho_1= 0.0, normb= norm( b), resid;

    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    if ((resid= norm( r)/normb) <= tol)
    {
        tol= resid;
        max_iter= 0;
        return true;
    }

    for (int i= 1; i <= max_iter; ++i)
    {
        M.Apply(A, z, r);
        rho= dot( r, z);
        if (i == 1)
            p= z;
        else
            z_xpay(p, z, (rho/rho_1), p); // p= z + (rho/rho_1)*p;

        q= A*p;
        const double alpha= rho/dot( p, q);
        axpy( alpha, p, x);                // x+= alpha*p;
        axpy( -alpha, q, r);               // r-= alpha*q;

        if ((resid= norm( r)/normb)<= tol)
        {
            tol= resid;
            max_iter= i;
            return true;
        }
        rho_1= rho;
    }
    tol= resid;
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

//-----------------------------------------------------------------------------
// GMRES: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//
// x - approximate solution to Ax = b
// max_iter - number of iterations performed before tolerance was reached
// tol - 2-norm of the (relative, see below) preconditioned residual after the
//     final iteration
//
// measure_relative_tol - If true, stop if |M^(-1)( b - Ax)|/|M^(-1)b| <= tol,
//     if false, stop if |M^(-1)( b - Ax)| <= tol.
// calculate2norm - If true, the unpreconditioned (absolute) residual is
//     calculated in every step for debugging. This is rather expensive.
//-----------------------------------------------------------------------------
template <typename Mat, typename Vec, typename PreCon>
bool
GMRES(const Mat& A, Vec& x, const Vec& b, const PreCon& M,
      int /*restart parameter*/ m, int& max_iter, double& tol,
      bool measure_relative_tol= true, bool calculate2norm= false)
{
    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    DMatrixCL<double> H( m, m);
    Vec               s( m), cs( m), sn( m), w( b.size()), r( b.size());
    std::vector<Vec>  v( m);
    double            beta, normb, resid;

    for (int i= 0; i < m; ++i)
        v[i].resize( b.size());

    M.Apply( A, r, Vec( b-A*x));
    beta= norm( r);

    M.Apply( A, w, b);
    normb= norm( w);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    resid = beta/normb;
    if (resid <= tol) {
        tol= resid;
        max_iter= 0;
        return true;
    }

    int j= 1;
    while (j <= max_iter) {
        v[0]= r*(1.0/beta);
        s= 0.0;
        s[0]= beta;

        int i;
        for (i= 0; j <= max_iter; ++i, ++j) {
            M.Apply( A, w, A*v[i]);
            for (int k= 0; k <= i; ++k ) {
                H( k, i)= dot( w, v[k]);
                w-= H( k, i)*v[k];
            }

            if (i == m - 1) break;

            H( i + 1, i)= norm( w);
            v[i + 1]= w*(1.0/H( i + 1, i));

            for (int k= 0; k < i; ++k)
                GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

            GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
            GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
            GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

            resid= std::abs( s[i+1])/normb;
            if (calculate2norm == true) { // debugging aid
                Vec y( x);
                GMRES_Update( y, i, H, s, v);
                double resid2= norm( Vec( b - A*y));
                std::cerr << "GMRES: absolute residual 2-norm: " << resid2
                          << "\tabsolute preconditioned residual 2-norm: "
                          << std::fabs( s[i+1]) << '\n';
            }
            if (resid <= tol) {
                GMRES_Update( x, i, H, s, v);
                tol= resid;
                max_iter= j;
                return true;
            }
        }
        GMRES_Update( x, i, H, s, v);
        M.Apply( A, r, Vec( b - A*x));
        beta= norm( r);
        resid= beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= j;
            return true;
        }
    }
    tol= resid;
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
    q2= A*q1 - b0*q0;
    a1= dot( q2, q1);
    q2-= a1*q1;
    b1= norm( q2);
    // Lucky breakdown; the Krylov-space K up to q1 is A-invariant. Thus,
    // the correction dx needed to solve A(x0+dx)=b is in this space and
    // the Minres-algo will terminate with the exact solution in the
    // following step.
    if (b1 < 1e-15) return false;
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
    const Mat* A;
    SBufferCL<Vec, 3> q;
    double a0;
    SBufferCL<double, 2> b;

    LanczosONBCL()
      :A( 0) {}

    void // Sets up initial values and computes q0.
    new_basis(const Mat& A_, const Vec& r0) {
        A= &A_;
        q[-1].resize( r0.size(), 0.);
        norm_r0_= norm( r0);
        q[0].resize( r0.size(), 0.); q[0]= r0/norm_r0_;
        q[1].resize( r0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= LanczosStep( *A, q[-1], q[0], q[1], a0, b[-1], b[0]);
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
        return (nobreakdown_= LanczosStep( *A, q[-1], q[0], q[1], a0, b[-1], b[0]));
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
    t2= A*q1 - b0*t0;
    a1= dot( t2, q1);
    t2-= a1*t1;
    M.Apply( A, q2, t2);
    const double b1sq= dot( q2, t2);
    Assert( b1sq >= 0.0, "PLanczosStep: b1sq is negative!\n", DebugNumericC);
    b1= std::sqrt( b1sq);
    if (b1 < 1e-15) return false;
    t2*= 1./b1;
    q2*= 1./b1;
    return true;
}

template <typename Mat, typename Vec, typename PreCon>
class PLanczosONBCL
{
  private:
    bool nobreakdown_;
    double norm_r0_;

  public:
    const Mat* A;
    const PreCon& M;
    SBufferCL<Vec, 2> q;
    SBufferCL<Vec, 3> t;
    double a0;
    SBufferCL<double, 2> b;

    // Sets up initial values and computes q0.
    PLanczosONBCL(const PreCon& M_)
      :A( 0), M( M_) {}

    void // Sets up initial values and computes q0.
    new_basis(const Mat& A_, const Vec& r0) {
        A= &A_;
        t[-1].resize( r0.size(), 0.);
        q[-1].resize( r0.size(), 0.); M.Apply( *A, q[-1], r0);
        norm_r0_= std::sqrt( dot( q[-1], r0));
        t[0].resize( r0.size(), 0.); t[0]= r0/norm_r0_;
        q[0].resize( r0.size(), 0.); q[0]= q[-1]/norm_r0_;
        t[1].resize( r0.size(), 0.);
        b[-1]= 0.;
        nobreakdown_= PLanczosStep( *A, M, q[0], q[1], t[-1], t[0], t[1], a0, b[-1], b[0]);
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
        return (nobreakdown_= PLanczosStep( *A, M, q[0], q[1], t[-1], t[0], t[1], a0, b[-1], b[0]));
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
//      tol - (relative, see below) residual b - Ax measured in the (M^-1 ., .)-
//     inner-product-norm.
// measure_relative_tol - If true, stop if (M^(-1)( b - Ax), b - Ax)/(M^(-1)b, b) <= tol,
//     if false, stop if (M^(-1)( b - Ax), b - Ax) <= tol.
//-----------------------------------------------------------------------------
template <typename Mat, typename Vec, typename Lanczos>
bool
PMINRES(const Mat&, Vec& x, const Vec&, Lanczos& q, int& max_iter, double& tol,
    bool measure_relative_tol= false)
{
    Vec dx( x.size());
    const double norm_r0= q.norm_r0();
    double normb= std::fabs( norm_r0);
    double res= norm_r0;
    bool lucky= q.breakdown();
    SBufferCL<double, 3> c;
    SBufferCL<double, 3> s;
    SBufferCL<SVectorCL<3>, 3> r;
    SBufferCL<Vec, 3> p;
    p[0].resize( x.size()); p[1].resize( x.size()); p[2].resize( x.size());
    SBufferCL<SVectorCL<2>, 2> b;

    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

    if ((res= norm_r0/normb) <= tol) {
        tol= res;
        max_iter= 0;
        return true;
    }

    for (int k= 1; k <= max_iter; ++k) {
        switch (k) {
          case 1:
            // Compute r1
            GMRES_GeneratePlaneRotation( q.a0, q.b[0], c[0], s[0]);
            r[0][0]= std::sqrt( q.a0*q.a0 + q.b[0]*q.b[0]);
            // Compute p1
            // p[0]= q.q[0]/r[0][0];
            p[0]= q.q[0]/r[0][0];
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
            p[0]= (q.q[0] - r[0][0]*p[-1])/r[0][1];
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
            p[0]= (q.q[0] - r[0][0]*p[-2] -r[0][1]*p[-1])*(1/r[0][2]);
            b[0][0]= b[-1][1]; b[0][1]= 0.;
            GMRES_ApplyPlaneRotation( b[0][0], b[0][1], c[0], s[0]);
        }
        dx= norm_r0*b[0][0]*p[0];
        x+= dx;

        res= std::fabs( norm_r0*b[0][1])/normb;
//        std::cerr << "PMINRES: residual: " << res << '\n';
        if (res<= tol || lucky==true) {
            tol= res;
            max_iter= k;
            return true;
        }
        q.next();
        if (q.breakdown()) {
            lucky= true;
            std::cerr << "PMINRES: lucky breakdown\n";
        }
        c.rotate(); s.rotate(); r.rotate(); p.rotate(); b.rotate();
    }
    tol= res;
    return false;
}


template <typename Mat, typename Vec>
bool
MINRES(const Mat& A, Vec& x, const Vec& rhs, int& max_iter, double& tol,
    bool measure_relative_tol= false)
{
    LanczosONBCL<Mat, Vec> q;
    q.new_basis( A, Vec( rhs - A*x));
    return PMINRES( A,  x, rhs, q, max_iter, tol, measure_relative_tol);
}


//*****************************************************************
// BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = b 
// using the Preconditioned BiConjugate Gradient Stabilized method
//
// BiCGSTAB follows the algorithm described on p. 27 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence or breakdown within
// max_iter iterations (false). In cases of breakdown a message is printed
// std::cerr and the iteration returns the approximate solution found. 
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol. ( |.| is the euclidean norm.)
//  
//*****************************************************************
template <class Mat, class Vec, class Preconditioner>
bool
BICGSTAB( const Mat& A, Vec& x, const Vec& b,
    const Preconditioner& M, int& max_iter, double& tol,
    bool measure_relative_tol= true)
{
    double rho_1= 0.0, rho_2= 0.0, alpha= 0.0, beta= 0.0, omega= 0.0;
    Vec p( x.size()), phat( x.size()), s( x.size()), shat( x.size()),
        t( x.size()), v( x.size());

    double normb= norm( b);
    Vec r( b - A*x);
    Vec rtilde= r;

    if (normb == 0.0 || measure_relative_tol == false) normb = 1.0;
  
    double resid= norm( r)/normb;
    if (resid <= tol) {
        tol= resid;
        max_iter= 0;
        return true;
    }

    for (int i= 1; i <= max_iter; ++i) {
        rho_1= dot( rtilde, r);
        if (rho_1 == 0.0) {
            tol = norm( r)/normb;
            max_iter= i;
            std::cerr << "BiCGSTAB: Breakdown with rho_1 = 0.\n";
            return false;
        }
        if (i == 1) p= r;
        else {
            beta= (rho_1/rho_2)*(alpha/omega);
            p= r + beta*(p - omega*v);
        }
        M.Apply( A, phat, p);
        v= A*phat;
        alpha= rho_1/dot( rtilde, v);
        s= r - alpha*v;
        if ((resid= norm( s)/normb) < tol) {
            x+= alpha*phat;
            tol= resid;
            max_iter= i;
            return true;
        }
        M.Apply( A, shat, s);
        t= A*shat;
        omega= dot( t, s)/dot( t, t);
        x+= alpha*phat + omega*shat;
        r= s - omega*t;

        rho_2= rho_1;
        if ((resid= norm( r)/normb) < tol) {
            tol= resid;
            max_iter= i;
            return true;
        }
        if (omega == 0.0) {
            tol= norm( r)/normb;
            max_iter= i;
            std::cerr << "BiCGSTAB: Breakdown with omega_ = 0.\n";
            return false;
        }
    }
    tol= resid;
    return false;
}

//*****************************************************************
// GCR
//
// GCR solves the unsymmetric linear system Ax = b 
// using the Preconditioned Generalized Conjugate Residuals method;
//
// The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within
// max_iter iterations (false).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
// measure_relative_tol - If true, stop if |b - Ax|/|b| <= tol,
//     if false, stop if |b - Ax| <= tol. ( |.| is the euclidean norm.)
//  
//*****************************************************************
template <class Mat, class Vec, class Preconditioner>
bool
GCR( const Mat& A, Vec& x, const Vec& b, const Preconditioner& M,
    int /*truncation parameter m*/, int& max_iter, double& tol,
    bool measure_relative_tol= true)
{
    Vec r( b - A*x);
    std::vector<Vec> s( 1), v( 1); // Positions s[0], v[0] are unused below.
    double normb= norm( b);
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;
    double resid= -1.0;

    for (int k= 0; k < max_iter; ++k) {
        if ((resid= norm( r)/normb) < tol) {
            tol= resid;
            max_iter= k;
            return true;
        }
        s.push_back( Vec( b.size()));
        M.Apply( A, s[k+1], r);
        v.push_back( A*s[k+1]);
        for (int i= 1; i <= k; ++i) {
            const double alpha= dot( v[k+1], v[i]);
            v[k+1]-= alpha*v[i];
            s[k+1]-= alpha*s[i];
        }
        const double beta= norm( v[k+1]);
        v[k+1]/= beta;
        s[k+1]/= beta;
        const double gamma= dot( r, v[k+1]);
        x+= gamma*s[k+1];
        r-= gamma*v[k+1];
    }
    tol= resid;
    return false;
}


//=============================================================================
//  Drivers
//=============================================================================

// What every iterative solver should have
class SolverBaseCL
{
  protected:
    int         _maxiter;
    mutable int _iter;
    double         _tol;
    mutable double _res;
    bool rel_;

    SolverBaseCL (int maxiter, double tol, bool rel= false)
        : _maxiter( maxiter), _iter( -1), _tol( tol), _res( -1.), rel_( rel)  {}

  public:
    void   SetTol     (double tol) { _tol= tol; }
    void   SetMaxIter (int iter)   { _maxiter= iter; }
    void   SetRelError(bool rel)   { rel_= rel; }

    double GetTol     () const { return _tol; }
    int    GetMaxIter () const { return _maxiter; }
    double GetResid   () const { return _res; }
    int    GetIter    () const { return _iter; }
    bool   GetRelError() const { return rel_; }

};


// Bare CG solver
class CGSolverCL : public SolverBaseCL
{
  public:
    CGSolverCL(int maxiter, double tol, bool rel= false)
        : SolverBaseCL(maxiter, tol, rel) {}

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        CG(A, x, b, _iter, _res, rel_);
    }
    template <typename Mat, typename Vec>
    void Solve(const MatrixCL& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        CG(A, x, b, numIter, resid, rel_);
    }
};


// With preconditioner
template <typename PC>
class PCGSolverCL : public SolverBaseCL
{
  private:
    PC _pc;

  public:
    PCGSolverCL(const PC& pc, int maxiter, double tol, bool rel= false)
        : SolverBaseCL(maxiter, tol, rel), _pc(pc) {}

    PC&       GetPc ()       { return _pc; }
    const PC& GetPc () const { return _pc; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        PCG(A, x, b, _pc, _iter, _res, rel_);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        PCG(A, x, b, _pc, numIter, resid, rel_);
    }
};

// Bare MINRES solver
class MResSolverCL : public SolverBaseCL
{
  public:
    MResSolverCL(int maxiter, double tol, bool rel= false)
        : SolverBaseCL( maxiter, tol, rel) {}

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MINRES( A, x, b, _iter, _res, rel_);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        MINRES( A, x, b, numIter, resid, rel_);
    }
};

// Preconditioned MINRES solver
template <typename Lanczos>
class PMResSolverCL : public SolverBaseCL
{
  private:
    Lanczos& q_;

  public:
    PMResSolverCL(Lanczos& q, int maxiter, double tol, bool rel= false)
      :SolverBaseCL( maxiter, tol, rel), q_( q) {}

    Lanczos&       GetONB ()       { return q_; }
    const Lanczos& GetONB () const { return q_; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        q_.new_basis( A, Vec( b - A*x));
        PMINRES( A, x, b, q_, _iter, _res, rel_);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        q_.new_basis( A, Vec( b - A*x));
        PMINRES( A, x, b, q_, numIter, resid, rel_);
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
    GMResSolverCL(const PC& pc, int restart, int maxiter, double tol, bool relative= true)
        : SolverBaseCL( maxiter, tol, relative), pc_(pc), restart_(restart) {}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }
    int       GetRestart () const { return restart_; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        GMRES(A, x, b, pc_, restart_, _iter, _res, rel_);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        GMRES(A, x, b, pc_, restart_, numIter, resid, rel_);
    }
};

// BiCGStab
template <typename PC>
class BiCGStabSolverCL : public SolverBaseCL
{
  private:
    PC pc_;

  public:
    BiCGStabSolverCL(const PC& pc, int maxiter, double tol, bool relative= true)
        : SolverBaseCL( maxiter, tol, relative), pc_( pc){}

    PC&       GetPc ()       { return pc_; }
    const PC& GetPc () const { return pc_; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        BICGSTAB( A, x, b, pc_, _iter, _res, rel_);
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        BICGSTAB(A, x, b, pc_, numIter, resid, rel_);
    }
};

// GCR
template <typename PC>
class GCRSolverCL : public SolverBaseCL
{
  private:
    PC  pc_;
    int truncate_; // no effect atm.

  public:
    GCRSolverCL(const PC& pc, int truncate, int maxiter, double tol, bool relative= true)
        : SolverBaseCL( maxiter, tol, relative), pc_( pc), truncate_( truncate) {}

    PC&       GetPc      ()       { return pc_; }
    const PC& GetPc      () const { return pc_; }
    int       GetTruncate() const { return truncate_; }

    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        GCR( A, x, b, pc_, truncate_, _iter, _res, rel_);
        std::cerr << "GCRSolverCL iterations: " << GetIter()
                  << "\tresidual: " << GetResid() << std::endl;
    }
    template <typename Mat, typename Vec>
    void Solve(const Mat& A, Vec& x, const Vec& b, int& numIter, double& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        GCR( A, x, b, pc_, truncate_, numIter, resid, rel_);
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
typedef PreGSCL<P_JAC0>    JACPcCL;
typedef PreGSCL<P_SGS0>    SGSPcCL;
typedef PreGSCL<P_SSOR0>   SSORPcCL;
typedef PreGSCL<P_SSOR0_D> SSORDiagPcCL;
typedef PreGSCL<P_DUMMY>   DummyPcCL;
typedef PreGSCL<P_GS0>     GSPcCL;

typedef PCGSolverCL<SGSPcCL>      PCG_SgsCL;
typedef PCGSolverCL<SSORPcCL>     PCG_SsorCL;
typedef PCGSolverCL<SSORDiagPcCL> PCG_SsorDiagCL;


//=============================================================================
// Krylov-methods as preconditioner. 
// Needed for InexactUzawa. These shall be replaced by an MG-preconditioner.
//=============================================================================
class SSORPCG_PreCL
{
  private:
    mutable PCGSolverCL<SSORPcCL> solver_;

  public:
    SSORPCG_PreCL(int maxiter= 500, double reltol= 0.02)
        : solver_( SSORPcCL( 1.0), maxiter, reltol, /*relative residuals*/ true) {}                                                                // real tol for

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, Vec& x, const Vec& b) const {
        x= 0.0;
        solver_.Solve( A, x, b);
    }
};


template <class PC>
class GMRes_PreCL : public SolverBaseCL
{
  private:
    PC& pc_;
    int restart_;

  public:
    GMRes_PreCL(PC& pc, int maxiter= 500, int restart= 250, double reltol= 0.02)
        : SolverBaseCL( maxiter, reltol, /*relative residual*/true), pc_( pc), restart_( restart)
    {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, Vec& x, const Vec& b) const {
        _res= _tol;
        _iter= _maxiter;
        GMRES( A, x, b, pc_, restart_, _iter, _res,
            /*relative errors!*/ true, /*don't check 2-norm*/ false);
//        std::cerr << "GMRes_PreCL iterations: " << GetIter()
//                  << "\trelative residual: " << GetResid() << std::endl;
    }
};

typedef GMRes_PreCL<DummyPcCL>  DummyGMRes_PreCL;
typedef GMRes_PreCL<JACPcCL>    JACGMRes_PreCL;
typedef GMRes_PreCL<GSPcCL>     GSGMRes_PreCL;
typedef GMRes_PreCL<SSORPcCL>   SSORGMRes_PreCL;


template <class PC>
class BiCGStab_PreCL : public SolverBaseCL
{
  private:
    PC& pc_;

  public:
    BiCGStab_PreCL(PC& pc, int maxiter= 500, double reltol= 0.02)
        : SolverBaseCL( maxiter, reltol, /*relative residual*/true), pc_( pc)
    {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& A, Vec& x, const Vec& b) const {
        _res= _tol;
        _iter= _maxiter;
        BICGSTAB( A, x, b, pc_, _iter, _res,
            /*relative errors!*/ true);
        std::cerr << "BiCGStab_PreCL iterations: " << GetIter()
                  << "\trelative residual: " << GetResid() << std::endl;
    }
};

typedef BiCGStab_PreCL<DummyPcCL> DummyBiCGStab_PreCL;
typedef BiCGStab_PreCL<JACPcCL>   JACBiCGStab_PreCL;
typedef BiCGStab_PreCL<GSPcCL>    GSBiCGStab_PreCL;
typedef BiCGStab_PreCL<SSORPcCL>  SSORBiCGStab_PreCL;


//*****************************************************************************
//
//  Non-iterative methods: Gauss solver with pivoting
//
//*****************************************************************************

template<Uint _Dim>
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
        max= std::fabs(A(p[i], i));
        ind_max= i;
        for (size_t l=i+1; l<_Dim; ++l)
            if (std::fabs(A(p[l], i))>max)
            {
                max= std::fabs(A(p[l], i));
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
