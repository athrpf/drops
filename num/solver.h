//**************************************************************************
// File:    solver.h                                                       *
// Content: iterative solvers                                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 3 2001                                          *
//                                                                         *
// Remarks: We should use the const-qualifier to make it difficult to      *
//          accidentally change the multigrid structure from anywhere      *
//          outside of the multigrid algorithms.                           *
//          Thus the pointer to user data structures should probably be    *
//          a pointer to mutable.                                          *
//**************************************************************************


#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "num/spmat.h"
#include <algorithm>

namespace DROPS
{


//**************************************************************************
// Iterative template routine - CG                                         *
// The return value indicates convergence within max_iter (input)          *
// iterations (true), or no convergence within max_iter iterations (false).*
//                                                                         *
// Upon successful return, output arguments have the following values:     *
//                                                                         *
//        x - approximate solution to Ax = b                               *
// max_iter - the number of iterations performed before the                *
//            tolerance was reached                                        *
//      tol - the residual after the final iteration                       *
//**************************************************************************
template <class Matrix, class Vector, class Real>
bool
CG(const Matrix& A, Vector& x, const Vector& b,
   int& max_iter, Real& tol)
{
    Vector r= A*x - b;
    Vector d= -r;
    Real gamma= r.norm2();

    tol*= tol;

    if (gamma<tol)
    {
        max_iter= 0;
        tol= sqrt(gamma);
        return true;
    }

    for (int it=0; it<=max_iter; ++it)
    {
        const Vector Ad= A*d;
        const Real delta= Ad*d;
/*        if (delta<1.e-20)
        {
            max_iter= 0;
            tol= sqrt(gamma);
            std::cerr << "res= " << r.norm() 
                      << ", gamma= " << gamma
                      << ", delta= " << delta << std::endl;
            return true;
        }
*/      const Real alpha= gamma/delta;
        x+= alpha*d;
        r+= alpha*Ad;
        Real beta= gamma;
//        std::cerr << "aeussere Iteration "<<it<<" mit Res.= " << r.norm() <<std::endl;
        if ( (gamma= r.norm2()) < tol ) 
        {
            max_iter= it;
            tol= sqrt(gamma);
            return true;
        }
        beta= gamma / beta;
        d*= beta;
        d-= r;
    }
    tol= sqrt(gamma);
    return false;
}


//**************************************************************************
// Iterative template routine -- PCG                                       *
// The return value indicates convergence within max_iter (input)          *
// iterations (true), or no convergence within max_iter iterations (false).*
// Upon successful return, output arguments have the following values:     *
//                                                                         *
//        x  --  approximate solution to Ax = b                            *
// max_iter  --  the number of iterations performed before the             *
//               tolerance was reached                                     *
//      tol  --  the residual after the final iteration                    *
//**************************************************************************

template < class Matrix, class Vector, class Preconditioner, class Real >
bool
PCG(const Matrix& A, Vector& x, const Vector& b,
    const Preconditioner& M, int& max_iter, Real& tol)
{
    const size_t n= x.size();
    Real resid;
    Vector p(n), z(n), q(n);
    Real alpha, beta, rho, rho_1;
    Vector r= b - A*x;
    tol*= tol;

    if ( (resid= r.norm2()) <= tol) 
    {
        tol= sqrt(resid);
        max_iter= 0;
        return true;
    }
    for (int i= 1; i <= max_iter; ++i) 
    {
        M.Apply(A, z, r);
        rho= r*z;
        if (i == 1) p= z;
        else 
        {
            beta= rho/rho_1;
            p= z + beta*p;
        }
        q= A*p;
        alpha= rho/(p*q);
        x+= alpha*p;
        r-= alpha*q;
        if ( (resid= r.norm2()) <= tol) 
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


/********************************************************************************************
*
*                           Gauss-Seidel methods
*
********************************************************************************************/


template <class Vec, class Real>
void GSStep(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b)
// performs one step of the Symmetric Gauss-Seidel method with start vector x
{
    const size_t n= A.num_rows();
    Real aii;

    for(size_t i=0, nz=0; i<n; ++i)
    {
        x[i]= b[i];
        for(const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                x[i]-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        x[i]/= aii;
    }
}


template <class Vec, class Real>
void
SGSStep(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b)
// performs one step of the Symmetric Gauss-Seidel method with start vector x
{
    const size_t n= A.num_rows();
    Real aii;

    for(size_t i=0, nz=0; i<n; ++i)
    {
        x[i]= b[i];
        for(const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                x[i]-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        x[i]/= aii;
    }

    for(long int i=n-1, nz= A.row_beg(n)-1; i>=0; --i)
    {
        x[i]= b[i];
        for(const size_t beg= A.row_beg(i); nz>=beg; --nz)
            if (col_ind(nz) != i)
                x[i]-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        x[i]/= aii;
    }
}


template <class Vec, class Real>
void
SORStep(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, Real omega)
// performs one step of the Symmetric SOR method with start vector x
{
    const size_t n= A.num_rows();
    Real aii, sum;

    for(size_t i=0, nz=0; i<n; ++i)
    {
        sum= b[i];
        for(const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        sum*= omega/aii;
	x[i]*= 1. - omega;
        x[i]+= sum;
    }

}
template <class Vec, class Real>
void
SSORStep(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, Real omega)
// performs one step of the Symmetric SOR method with start vector x
{
    const size_t n= A.num_rows();
    Real aii, sum;

    for(size_t i=0, nz=0; i<n; ++i)
    {
        sum= -b[i];
        for(const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum+= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        sum*= omega/aii;
	x[i]*= 1. - omega;
        x[i]-= sum;
    }

    for(long int i=n-1, nz=A.row_beg(n)-1; i>=0; --i)
    {
        sum= -b[i];
        for(const size_t beg= A.row_beg(i); nz>=beg; --nz)
            if (A.col_ind(nz) != i)
                sum+= A.val(nz)*x[A.col_ind(nz)];
            else
                aii= A.val(nz);
        sum*= omega/aii;
	x[i]*= 1. - omega;
        x[i]-= sum;
    }
}


template <class Vec, class Real>
void
SSORPC(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, Real omega)
// performs one step of the Symmetric SOR method with start vector x
{
    typedef typename SparseMatBaseCL<Real>::subs_type*  coliterT;
    typedef typename SparseMatBaseCL<Real>::value_type* valiterT;

    const size_t n= A.num_rows();
    Real sum, *valit;
    coliterT colit, diag,
             *diagonal= new coliterT[n];
//std::cerr << "SSOR-Preconditioner\n";
//    Assert(n == x.size(), DROPSErrCL("SSOR: incompatible dimensions\n"));
    for(size_t i=0; i<n; ++i)
    {
        x[i]= b[i];
        colit= A.GetFirstCol(i);
        valit= A.GetFirstVal(i);
        diag= diagonal[i]= std::lower_bound( colit, A.GetFirstCol(i+1), i);
//std::cerr << i <<'='<< (*diag) <<"---------"<<std::endl; 
        for( ; colit!=diag; ++colit, ++valit)
            x[i]-= (*valit)*x[*colit];

        // now: a_ii = *valit
	x[i]*= omega/(*valit);
    }

    size_t i=n;
    do
    {
        --i;
        sum= 0;
        colit= A.GetFirstCol(i+1);
        valit= A.GetFirstVal(i+1);
        --colit; --valit;
//        diag= std::lower_bound( A.GetFirstCol(i), A.GetFirstCol(i+1), i); 
        for(diag= diagonal[i] ; colit != diag; --colit, --valit)
            sum+= (*valit)*x[*colit];

        // now: a_ii = *valit
        sum*= omega/(*valit);
	x[i]*= 2. - omega;
        x[i]-= sum;
    }
    while (i>0);
    delete[] diagonal;
}


/**********************************************************************************************************
*
*   preconditioner, smoother and solver classes, embedding the foregoing methods
*
**********************************************************************************************************/


template<class Vec, class Real>
class SsorPcCL
{ 
  private:
    Real _omega;

  public:
    SsorPcCL( Real om= 1.0)
        : _omega(om)  {}
    inline void
    Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
    {
        SSORPC( A, x, b, _omega);
    }
};

template<class MatrixT>
class SsorMassPcCL
{
  private:
    typedef SparseMatBaseCL<double> MatrixCL;
    typedef VectorBaseCL<double>    VectorCL;  
  
    const MatrixCL& _M;
    
  public:
    SsorMassPcCL( const MatrixCL& M)
        : _M( M) {}
    inline void Apply( const MatrixT&, VectorCL& x, const VectorCL& b) const
    {
        SSORPC( _M, x, b, 1.);
    }
};

template<class Vec, class Real>
class ImprovedSsorPcCL
{ 
  private:
    typedef typename SparseMatBaseCL<Real>::subs_type*  coliterT;
    typedef typename SparseMatBaseCL<Real>::value_type* valiterT;

    Real _omega;
    coliterT* _diagonal;

  public:
    ImprovedSsorPcCL( Real om= 1.0)
        : _omega(om), _diagonal(0)  {}
    ~ImprovedSsorPcCL() { delete _diagonal; }
    
    void Init(const SparseMatBaseCL<Real>& A);
      // call Init before calls to Apply, do so when matrix has changed
    void Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const;
};

template<class Vec, class Real>
class SGSPcCL
{ 
  private:
    typedef typename SparseMatBaseCL<Real>::subs_type*  coliterT;
    typedef typename SparseMatBaseCL<Real>::value_type* valiterT;

  public:
    void Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const;
};

template<class Vec, class Real>
class WSGSSmootherCL
{
  private:
    Real _omega;
    
  public:
    WSGSSmootherCL( Real om=1.0) : _omega(om) {}
    inline void Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
    {
        SSORStep( A, x, b, _omega);
    }
};


template<class Vec, class Real>
class WGSSmootherCL
{
  private:
    Real _omega;
    
  public:
    WGSSmootherCL( Real om=1.0) : _omega(om) {}
    inline void Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
    {
        SORStep( A, x, b, _omega);
    }
};


template<class Vec, class Real>
class CGSolverCL
{
  private:
    Real _tol, _res;
    int  _maxiter, _iter;
  
  public:
    CGSolverCL(Real tol, int maxiter)
        : _tol(tol), _res( -1.), _maxiter(maxiter), _iter( -1) {}
    
    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    
    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }

    inline void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        CG(A, x, b, _iter, _res);
    }
    inline void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, int& numIter, Real& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        CG(A, x, b, numIter, resid);
    }
};
  
template<class Vec, class Real, class PC>
class PCGSolverCL
{
  private:
    Real _tol, _res;
    int  _maxiter, _iter;
    PC   _pc;
  
  public:
    PCGSolverCL(Real tol, int maxiter, const PC pc)
        : _tol(tol), _res( -1.), _maxiter(maxiter), _iter( -1), _pc(pc) {}
    
    void SetTol      ( double tol) { _tol= tol; }
    void SetMaxIter  ( int iter  ) { _maxiter= iter; }
    
    double GetTol    () const { return _tol; }
    int    GetMaxIter() const { return _maxiter; }
    double GetResid  () const { return _res; }
    int    GetIter   () const { return _iter; }
    PC&    GetPc     ()       { return _pc; }

    inline void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        PCG(A, x, b, _pc, _iter, _res);
    }
    inline void Solve(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b, int& numIter, Real& resid) const
    {
        resid=   _tol;
        numIter= _maxiter;
        PCG(A, x, b, _pc, numIter, resid);
    }
};
  

/***********************************************************************************************************
*
*                                     other methods
*
***********************************************************************************************************/

template<Uint _Dim>
void
gauss_pivot(SMatrixCL<_Dim, _Dim>& A, SVectorCL<_Dim>& b)
{
    double max;
    Uint ind_max;
    Uint p[_Dim];
    SVectorCL<_Dim> b2= b;
    for (Uint i=0; i<_Dim; ++i) p[i]= i;

    for (Uint i=0; i<_Dim-1; ++i)
    {
        max= fabs(A(p[i], i));
        ind_max= i;
        for (Uint l=i+1; l<_Dim; ++l)
            if (fabs(A(p[l], i))>max)
            {
                max= fabs(A(p[l], i));
                ind_max= l;
            }
        if (max == 0.0) throw DROPSErrCL("gauss_pivot: Matrix is singular.");
        if (i!=ind_max) std::swap(p[i], p[ind_max]);
        const double piv= A(p[i], i);
        for (Uint j=i+1; j<_Dim; ++j)
        {
            const double fac= A(p[j], i)/piv;
            b2[p[j]]-= fac*b2[p[i]];
            for (Uint k=i+1; k<_Dim; ++k)
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


template<class Matrix, class Vector, class Real>
void
cgVerf(const Matrix& A, Vector& x, const Vector& b,
       int& iter, Real& eps)
{
    Vector r= b-A*x, p= r, c(b.size());
    Real alpha, beta, d, dhelp;
    eps*= eps;

    r= b-A*x;
    p= r;
    d= r*r;
    do
    {
        c= A*p;		// c ist Hilfsvektor, um Rechenzeit zu sparen
        alpha= d/(p*c);	// berechne alpha, x, r, beta und p neu
        x+= alpha*p;
        r-= alpha * c; 	// hier mutiert r von riMinus1 zu ri
        dhelp= r*r;
        beta= dhelp/d;
        d= dhelp;
        p= r + (beta*p);
        ++iter;
    }
    while (r.norm2() > eps);
    eps= sqrt(r.norm2());
}

//========================================
//        template definitions
//========================================

template <class Vec, class Real>
void ImprovedSsorPcCL<Vec,Real>::Init(const SparseMatBaseCL<Real>& A)
{
    const size_t n= A.num_rows();
    if (_diagonal) delete _diagonal;
    _diagonal= new coliterT[n];
    
    for(size_t i=0; i<n; ++i)
    {
        _diagonal[i]= std::lower_bound( A.GetFirstCol(i), A.GetFirstCol(i+1), i);
    }
}

template <class Vec, class Real>
void ImprovedSsorPcCL<Vec,Real>::Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
{
    const size_t n= A.num_rows();
    Real sum, *valit;
    coliterT colit, diag;
//std::cerr << "SSOR-Preconditioner\n";
//    Assert(n == x.size(), DROPSErrCL("SSOR: incompatible dimensions\n"));
    for (size_t i=0; i<n; ++i)
    {
        x[i]= b[i];
        colit= A.GetFirstCol(i);
        valit= A.GetFirstVal(i);
//std::cerr << i <<'='<< (*diag) <<"---------"<<std::endl; 
        for (diag= _diagonal[i]; colit!=diag; ++colit, ++valit)
            x[i]-= (*valit)*x[*colit];

        // now: a_ii = *valit
	x[i]*= _omega/(*valit);
    }

    size_t i=n;
    do
    {
        --i;
        sum= 0;
        colit= A.GetFirstCol(i+1);
        valit= A.GetFirstVal(i+1);
        --colit; --valit;
        diag= _diagonal[i];
        for (; colit != diag; --colit, --valit)
            sum+= (*valit)*x[*colit];

        // now: a_ii = *valit
        sum*= _omega/(*valit);
	x[i]*= 2. - _omega;
        x[i]-= sum;
    }
    while (i>0);
}

template <class Vec, class Real>
void SGSPcCL<Vec,Real>::Apply(const SparseMatBaseCL<Real>& A, Vec& x, const Vec& b) const
{
    const size_t n= A.num_rows();
    Real sum, *valit;
    coliterT colit, diag;
//std::cerr << "SSOR-Preconditioner\n";
//    Assert(n == x.size(), DROPSErrCL("SSOR: incompatible dimensions\n"));
    for (size_t i=0; i<n; ++i)
    {
        x[i]= b[i];
        colit= A.GetFirstCol(i);
        valit= A.GetFirstVal(i);
        diag= std::lower_bound( colit, A.GetFirstCol(i+1), i);
//std::cerr << i <<'='<< (*diag) <<"---------"<<std::endl; 
        for (; colit!=diag; ++colit, ++valit)
            x[i]-= (*valit)*x[*colit];

        // now: a_ii = *valit
	x[i]/= *valit;
    }

    size_t i=n;
    do
    {
        --i;
        sum= 0;
        colit= A.GetFirstCol(i+1);
        valit= A.GetFirstVal(i+1);
        diag= std::lower_bound( A.GetFirstCol(i), colit, i);
        --colit; --valit;
        for (; colit != diag; --colit, --valit)
            sum+= (*valit)*x[*colit];

        // now: a_ii = *valit
        sum/= *valit;
        x[i]-= sum;
    }
    while (i>0);
}

} // end of namespace DROPS

#endif
