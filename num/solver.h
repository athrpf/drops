//**************************************************************************
// File:    solver.h                                                       *
// Content: iterative solvers                                              *
// Authors: Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************


#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "spmat.h"

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
    P_JAC,    // Jacobi
    P_GS,     // Gauss-Seidel
    P_SGS,    // symmetric Gauss-Seidel
    P_SGS0,   // symmetric Gauss-Seidel with initial vector 0
    P_JOR,    // P_JAC with over-relaxation
    P_SOR,    // P_GS with over-relaxation
    P_SSOR,   // P_SGS with over-relaxation
    P_SSOR0,  // P_SGS0 with over-relaxation
    P_SGS0_D, // P_SGS0 using SparseMatDiagCL
    P_SSOR0_D // P_SSOR0 using SparseMatDiagCL
};

// Base methods
enum PreBaseGS { PB_JAC, PB_GS, PB_SGS, PB_SGS0 };

// Properties of the methods
template <PreMethGS PM> struct PreTraitsCL
{
    static const PreBaseGS BaseMeth= PreBaseGS(PM<8 ? PM%4 : PB_SGS0);
    static const bool      HasOmega= (PM>=4 && PM<8) || PM==10;
    static const bool      HasDiag=  PM>=8;
};

// Used to make a distinct type from each method
template <PreBaseGS> class PreDummyCL {};


//=============================================================================
//  Implementation of the methods
//=============================================================================

// One step of the Jacobi method with start vector x
template <bool HasOmega, typename Vec, typename Real>
void
SolveGSstep(const PreDummyCL<PB_JAC>&, const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();
    Vec          y(x.size());

    for (size_t i=0, nz=0; i<n; ++i)
    {
        Real aii, sum= b[i];
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

    swap(x,y);
}


// One step of the Gauss-Seidel/SOR method with start vector x
template <bool HasOmega, typename Vec, typename Real>
void
SolveGSstep(const PreDummyCL<PB_GS>&, const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0, nz=0; i<n; ++i)
    {
        Real aii, sum= b[i];
        for (const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }
}


// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector x
template <bool HasOmega, typename Vec, typename Real>
void
SolveGSstep(const PreDummyCL<PB_SGS>&, const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0, nz=0; i<n; ++i)
    {
        Real aii, sum= b[i];
        for (const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        if (HasOmega)
            x[i]= (1.-omega)*x[i]+omega*sum/aii;
        else
            x[i]= sum/aii;
    }

    for (size_t i=n, nz=A.row_beg(n); i>0; )
    {
        --i;
        Real aii, sum= b[i];
        for (const size_t beg= A.row_beg(i); nz>beg; )
        {
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
template <bool HasOmega, typename Vec, typename Real>
void
SolveGSstep(const PreDummyCL<PB_SGS0>&, const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        Real sum= b[i];
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
        Real sum= 0;
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
template <bool HasOmega, typename Vec, typename Real>
void
SolveGSstep(const PreDummyCL<PB_SGS0>&, const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, const SparseMatDiagCL& diag, double omega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        Real sum= b[i];
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
        Real sum= 0;
        for (size_t j= A.row_beg(i+1)-1; j > diag[i]; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/A.val(diag[i]);
        else
            x[i]+= sum/A.val(diag[i]);
    }
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

    template <typename Vec, typename Real>
    void Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec,Real>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A, x, b, _omega);
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
    PreGSCL (const PreGSCL& p) : _diag(new SparseMatDiagCL(*(p._diag))), _omega(p._omega) {}
    ~PreGSCL() { delete _diag; }

    template <typename Real>
    void Init(const SparseMatBaseCL<Real>& A)
    {
        delete _diag; _diag=new SparseMatDiagCL(A);
    }

    template <typename Vec, typename Real>
    void Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec,Real>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), A, x, b, *_diag, _omega);
    }
};


// Preconditioners with own matrix
template <PreMethGS PM, typename Real, bool HasDiag= PreTraitsCL<PM>::HasDiag>
class PreGSOwnMatCL;

// Simple preconditioners
template <PreMethGS PM, typename Real>
class PreGSOwnMatCL<PM,Real,false>
{
  private:
    const SparseMatBaseCL<Real>& _M;
    double                       _omega;

  public:
    PreGSOwnMatCL (const SparseMatBaseCL<Real>& M, double om= 1.0) : _M(M), _omega(om) {}
    PreGSOwnMatCL (const PreGSOwnMatCL&); // not defined

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec,Real>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), _M, x, b, _omega);
    }
};


// Preconditioner with SparseMatDiagCL
template <PreMethGS PM, typename Real>
class PreGSOwnMatCL<PM,Real,true>
{
  private:
    const SparseMatBaseCL<Real>& _M;
    const SparseMatDiagCL*       _diag;
    double                       _omega;

  public:
    PreGSOwnMatCL (const SparseMatBaseCL<Real>& M, double om= 1.0)
      : _M(M), _diag(0), _omega(om) {}
    PreGSOwnMatCL (const PreGSOwnMatCL&); // not defined
    ~PreGSOwnMatCL() { delete _diag; }

    void Init(const SparseMatBaseCL<Real>& A)
    {
        delete _diag; _diag=new SparseMatDiagCL(A);
    }

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& x, const Vec& b) const
    {
        SolveGSstep<PreTraitsCL<PM>::HasOmega,Vec,Real>(PreDummyCL<PreTraitsCL<PM>::BaseMeth>(), _M, x, b, *_diag, _omega);
    }
};


//*****************************************************************************
//
//  Conjugate gradients - CG, PCG
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

template <typename Matrix, typename Vector, typename Real>
bool
CG(const Matrix& A, Vector& x, const Vector& b, int& max_iter, Real& tol)
{
    Vector r= A*x - b;
    Vector d= -r;
    Real resid= r.norm2();

    tol*= tol;

    if (resid<=tol)
    {
        tol= sqrt(resid);
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        const Vector Ad= A*d;
        const Real   delta= Ad*d;
        const Real   alpha= resid/delta;
        Real         beta= resid;

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

template <typename Matrix, typename Vector, typename PreCon, typename Real>
bool
PCG(const Matrix& A, Vector& x, const Vector& b, const PreCon& M,
    int& max_iter, Real& tol)
{
    const size_t n= x.size();
    Vector p(n), z(n), q(n), r= b - A*x;
    Real rho, rho_1= 0.0, resid= r.norm2();

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
        const Real alpha= rho/(p*q);
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


//=============================================================================
//  Drivers for the CG methods
//=============================================================================

// Bare CG solver
class CGSolverCL
{
  protected:
    double _tol, _res;
    int    _maxiter, _iter;

  public:
    CGSolverCL(double tol, int maxiter)
        : _tol(tol), _res( -1.), _maxiter(maxiter), _iter( -1) {}

    void   SetTol    (double tol) { _tol= tol; }
    void   SetMaxIter(int iter)   { _maxiter= iter; }

    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }

    template <typename Vec, typename Real>
    void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        CG(A, x, b, _iter, _res);
    }
    template <typename Vec, typename Real>
    void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, int& numIter, Real& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        CG(A, x, b, numIter, resid);
    }
};


// Add the preconditioner
template <typename PC>
class PCGSolverCL : public CGSolverCL
{
  private:
    PC _pc;

  public:
    PCGSolverCL(double tol, int maxiter, const PC& pc)
        : CGSolverCL(tol,maxiter), _pc(pc) {}

    PC&       GetPc ()       { return _pc; }
    const PC& GetPc () const { return _pc; }

    template <typename Vec, typename Real>
    void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        PCG(A, x, b, _pc, _iter, _res);
    }
    template <typename Vec, typename Real>
    void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, int& numIter, Real& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        PCG(A, x, b, _pc, numIter, resid);
    }
};


//=============================================================================
//  Typedefs
//=============================================================================

typedef PreGSCL<P_SOR>     SORsmoothCL;
typedef PreGSCL<P_SGS0>    SGSPcCL;
typedef PreGSCL<P_SSOR0>   SSORPcCL;
typedef PreGSCL<P_SSOR0_D> SSORDiagPcCL;

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
