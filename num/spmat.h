//**************************************************************************
// File:     spmat.h                                                       *
// Content:  sparse matrix in compressed row format                        *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.2                                                           *
//**************************************************************************


#ifndef DROPS_SPMAT_H
#define DROPS_SPMAT_H

#include <iostream>
#include <valarray>
#include <numeric>
#include <map>
#include <misc/utils.h>

namespace DROPS
{

//*****************************************************************************
//
//  V e c t o r B a s e C L :   base class for numerical vectors
//                              based on std::valarray
//
//*****************************************************************************
template <typename T>
class VectorBaseCL: public std::valarray<T>
{
private:
#if (DROPSDebugC & DebugNumericC)
    void AssertDim (__UNUSED__ const std::valarray<T>& va,
                    __UNUSED__ const char msg[]) const
      { Assert( size()==va.size(), msg, DebugNumericC); }
#endif

public:
    typedef T value_type;
    typedef std::valarray<T> base_type;

    // ctors
    VectorBaseCL()                      : base_type()       {}
#ifdef VALARRAY_BUG
    VectorBaseCL (size_t s)             : base_type(T(),s)  {}
#else
    VectorBaseCL (size_t s)             : base_type(s)      {}
#endif
    VectorBaseCL (T c, size_t s)        : base_type(c, s)   {}
    VectorBaseCL (const T* tp, size_t s): base_type(tp, s)  {}
    template <class V>
      VectorBaseCL (V exptpl)           : base_type(exptpl) {}

    // assignment
    template <class VT>
      VectorBaseCL& operator=(const VT& v) {
        *static_cast<base_type*>( this)= v; return *this; }

    // computed assignment
    template <class VT>
    VectorBaseCL& operator+=(const VT& v) {
        *static_cast<base_type*>( this)+= v; return *this; }
    template <class VT>
    VectorBaseCL& operator-=(const VT& v) {
        *static_cast<base_type*>( this)-= v; return *this; }
    template <class VT>
    VectorBaseCL& operator*=(const VT& v) {
        *static_cast<base_type*>( this)*= v; return *this; }
    template <class VT>
    VectorBaseCL& operator/=(const VT& v) {
        *static_cast<base_type*>( this)/= v; return *this; }

#if (DROPSDebugC & DebugNumericC)
    // For checking new code; the following functions inhibit the expression
    // template mechanism of std::valarray, but do some checks before the
    // operations.

    // element access
    T  operator [] (size_t s) const;
    T& operator [] (size_t s);

    // assignment
    VectorBaseCL& operator=(const std::valarray<T>& va);
    VectorBaseCL& operator= (const VectorBaseCL& v);

    // computed assignment
    VectorBaseCL& operator+= (const VectorBaseCL& v);
    VectorBaseCL& operator-= (const VectorBaseCL& v);
    VectorBaseCL& operator*= (const VectorBaseCL& v);
    VectorBaseCL& operator/= (const VectorBaseCL& v);

    friend VectorBaseCL operator+ (const VectorBaseCL& v, const VectorBaseCL& w);
    friend VectorBaseCL operator- (const VectorBaseCL& v, const VectorBaseCL& w);
#endif
};


template <class VT1, class VT2>
  inline typename VT1::value_type
  dot(const VT1& v, const VT2& w)
{
#if (DROPSDebugC & DebugNumericC)
        v.AssertDim(w, "dot: incompatible dimensions");
#endif
    typedef typename VT1::value_type valueT;
    valueT ret= valueT();
    const size_t s= v.size();
    for (size_t i= 0; i<s; ++i) ret+= v[i]*w[i];
    return ret;
}

template <class VT>
  inline typename VT::value_type
  norm_sq(const VT& v)
{
    typedef typename VT::value_type valueT;
    valueT ret= valueT();
    const size_t s= v.size();
    for (size_t i= 0; i<s; ++i) ret+= v[i]*v[i];
    return ret;
}

template <class VT>
  inline typename VT::value_type
  norm(const VT& v)
{
    return std::sqrt( norm_sq( v));
}

template <class VT>
  inline typename VT::value_type
  supnorm(const VT& v)
{
    typedef typename VT::value_type valueT;
    const size_t s= v.size();
    valueT ret= s==0 ? valueT() : std::abs( v[0]);
    for (size_t i= 1; i < s; ++i) {
        const valueT t= std::abs( v[i]);
        if (t > ret) ret= t;
    }
    return ret;
}

template <typename T>
  inline void
  axpy(T a, const VectorBaseCL<T>& x, VectorBaseCL<T>& y)
{
    Assert(x.size()==y.size(), "axpy: incompatible dimensions", DebugNumericC);
    y+= a*x;
}

template <typename T>
  inline void
  z_xpay(VectorBaseCL<T>& z, const VectorBaseCL<T>& x, T a, const VectorBaseCL<T>& y)
{
    Assert(z.size()==x.size() && z.size()==y.size(),
        "z_xpay: incompatible dimensions", DebugNumericC);
    z= x + a*y;
}

template <typename T>
  inline void
  z_xpaypby2 (VectorBaseCL<T>& z, const VectorBaseCL<T>& x,
    T a, const VectorBaseCL<T>& y, T b, const VectorBaseCL<T>& y2)
{
    Assert(z.size()==x.size() && z.size()==y.size() && z.size()==y2.size(),
        "z_xpaypby2: incompatible dimensions", DebugNumericC);
    z= x + a*y + b*y2;
}


#if (DROPSDebugC & DebugNumericC)
template <typename T>
T  VectorBaseCL<T>::operator[](size_t s) const
{
    Assert(s<size(), "VectorBaseCL []: index out of bounds", DebugNumericC);
    return (*static_cast<base_type*>( this))[s];
}

template <typename T>
T& VectorBaseCL<T>::operator[](size_t s)
{
    Assert(s<size(), "VectorBaseCL []: index out of bounds", DebugNumericC);
    return (*static_cast<base_type*>( this))[s];
}

// assignment
template <typename T>
VectorBaseCL<T>& VectorBaseCL<T>::operator=(const std::valarray<T>& va)
{
    AssertDim(va, "VectorBaseCL =: incompatible dimensions");
    *static_cast<base_type*>( this)= va;
    return *this;
}

template <typename T>
VectorBaseCL<T>& VectorBaseCL<T>::operator=(const VectorBaseCL<T>& v)
{
    AssertDim(v, "VectorBaseCL =: incompatible dimensions");
    *static_cast<base_type*>( this)= v;
    return *this;
}

// computed assignment
template <typename T>
VectorBaseCL<T>& VectorBaseCL<T>::operator+=(const VectorBaseCL<T>& v)
{
    AssertDim(v, "VectorBaseCL +=: incompatible dimensions");
    *static_cast<base_type*>( this)+= v;
    return *this;
}

template <typename T>
VectorBaseCL<T>& VectorBaseCL<T>::operator-=(const VectorBaseCL<T>& v)
{
    AssertDim(v, "VectorBaseCL -=: incompatible dimensions");
    *static_cast<base_type*>( this)-= v;
    return *this;
}

template <typename T>
VectorBaseCL<T>& VectorBaseCL<T>::operator*=(const VectorBaseCL<T>& v)
{
    AssertDim(v, "VectorBaseCL *=: incompatible dimensions");
    *static_cast<base_type*>( this)*= v;
    return *this;
}

template <typename T>
VectorBaseCL<T>& VectorBaseCL<T>::operator/=(const VectorBaseCL<T>& v)
{
    AssertDim(v, "VectorBaseCL /=: incompatible dimensions");
    *static_cast<base_type*>( this)/= v;
    return *this;
}

template <typename T>
VectorBaseCL<T> operator+(const VectorBaseCL<T>& v, const VectorBaseCL<T>& w)
{
    v.AssertDim(w, "VectorBaseCL + VectorBaseCL: incompatible dimensions");
    return VectorBaseCL( static_cast<typename VectorBaseCL<T>::base_type>( v)
        + static_cast<typename VectorBaseCL<T>::base_type>( w));
}

template <typename T>
VectorBaseCL<T> operator-(const VectorBaseCL<T>& v, const VectorBaseCL<T>& w)
{
    v.AssertDim(w, "VectorBaseCL - VectorBaseCL: incompatible dimensions");
    return VectorBaseCL( static_cast<typename VectorBaseCL<T>::base_type>( v)
        - static_cast<typename VectorBaseCL<T>::base_type>( w));
}
#endif


// Get the address of the first element in a valarray
// ("&x[0]" doesn't work, because "operator[] const" only returns a value)
template <typename T>
  inline const T*
  Addr(const std::valarray<T>& x)
    { return &(const_cast<std::valarray<T>&>(x)[0]); }

template <typename T>
  inline T*
  Addr(std::valarray<T>& x)
    { return &(x[0]); }


//*****************************************************************************
//
//  S p a r s e M a t B u i l d e r C L :  used for setting up sparse matrices
//
//*****************************************************************************
template <typename T>
class SparseMatBaseCL;


template <typename T>
class SparseMatBuilderCL
{
public:
    typedef T                        valueT;
    typedef SparseMatBaseCL<T>       spmatT;
    typedef std::pair<size_t,valueT> entryT;
    typedef std::map<size_t,valueT>  couplT;

private:
    size_t  _rows;
    size_t  _cols;
    spmatT* _mat;
    couplT* _coupl;

public:
    SparseMatBuilderCL(spmatT* mat, size_t rows, size_t cols)
        : _rows(rows), _cols(cols), _mat(mat), _coupl(new couplT[_rows])
        { Assert( _coupl!=0, "SparseMatBuilderCL: out of memory", ~0); }
    ~SparseMatBuilderCL() { delete[] _coupl; }

    T& operator() (size_t i, size_t j)
    {
        Assert(i<_rows && j<_cols, "SparseMatBuilderCL (): index out of bounds", DebugNumericC);
        return _coupl[i][j];
    }

    void Build();
};


template <typename T>
void SparseMatBuilderCL<T>::Build()
{
    size_t nz= 0;
    for (size_t i= 0; i<_rows; ++i)
        for (typename couplT::const_iterator it= _coupl[i].begin(), end= _coupl[i].end(); it != end; ++it)
    // TODO: was machen wir mit Eintraegen kleiner 1e-18 ?
            if (it->second != 0)
                ++nz;

    _mat->resize(_rows, _cols, nz);

    nz= 0;
    for (size_t i=0; i<_rows; ++i)
    {
        _mat->_rowbeg[i]= nz;
        for (typename couplT::const_iterator it= _coupl[i].begin(), end= _coupl[i].end(); it != end; ++it)
            if (it->second != 0)
            {
                _mat->_colind[nz]= it->first;
                _mat->_val[nz]=    it->second;
                ++nz;
            }
        // the col_ind-entries in each row are sorted, as they were stored sorted in the map
    }
    _mat->_rowbeg[_rows]= nz;

    Assert( nz == _mat->num_nonzeros(), "SparseMatBuilderCL::Build: wrong count of nonzeros", ~0);

    delete[] _coupl;
    _coupl= 0;
}


//*****************************************************************************
//
//  S p a r s e M a t B a s e C L :  row based sparse matrix format
//                                   use SparseMatBuilderCL for setting up!
//
//*****************************************************************************

template <typename T>
class SparseMatBaseCL
{
private:
    size_t _rows;                   // global dimensions
    size_t _cols;

    std::valarray<size_t> _rowbeg;  // (_rows+1 entries, last entry must be <=_nz) index of first non-zero-entry in _val belonging to the row given as subscript
    std::valarray<size_t> _colind;  // (_nz elements) column-number of corresponding entry in _val
    std::valarray<T>      _val;     // nonzero-entries

    void resize (size_t rows, size_t cols, size_t nz)
        { _rows=rows; _cols=cols; _rowbeg.resize(rows+1); _colind.resize(nz); _val.resize(nz); }
    void resize_rows (size_t rows)            { _rows=rows; _rowbeg.resize(rows+1); }
    void resize_cols (size_t cols, size_t nz) { _cols=cols; _colind.resize(nz); }
    void resize_val  (size_t nz)              { _val.resize(nz); }

public:
    typedef T value_type;

    SparseMatBaseCL () {}
    // default copy-ctor, assignment-op, dtor
    SparseMatBaseCL (size_t rows, size_t cols, size_t nz)
        : _rows(rows), _cols(cols), _rowbeg(rows+1), _colind(nz), _val(nz) {}
    SparseMatBaseCL (size_t rows, size_t cols, size_t nz,
                     const T* valbeg , const size_t* rowbeg, const size_t* colindbeg)
        : _rows(rows), _cols(cols), _rowbeg(rowbeg, rows+1), _colind(colindbeg, nz), _val(valbeg, nz) {}

    const T*      raw_val() const { return Addr( _val); }
    T*            raw_val()       { return &_val[0]; }
    const size_t* raw_row() const { return Addr( _rowbeg); }
    size_t*       raw_row()       { return &_rowbeg[0]; }
    const size_t* raw_col() const { return Addr( _colind); }
    size_t*       raw_col()       { return &_colind[0]; }

    size_t num_rows     () const { return _rows; }
    size_t num_cols     () const { return _cols; }
    size_t num_nonzeros () const { return _val.size(); }

    size_t row_beg (size_t i) const { return _rowbeg[i]; }
    size_t col_ind (size_t i) const { return _colind[i]; }
    T      val     (size_t i) const { return _val[i]; }

    const size_t* GetFirstCol(size_t i) const { return Addr(_colind)+_rowbeg[i]; }
    const T*      GetFirstVal(size_t i) const { return Addr(_val)+_rowbeg[i]; }

    inline T  operator() (size_t i, size_t j) const;

    SparseMatBaseCL& operator*= (T c) { _val*= c; return *this; }
    SparseMatBaseCL& operator/= (T c) { _val/= c; return *this; }

    SparseMatBaseCL& LinComb (double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&);

    void clear() { resize(0,0,0); }

    friend class SparseMatBuilderCL<T>;
};

template <typename T>
T SparseMatBaseCL<T>::operator() (size_t i, size_t j) const
{
    Assert(i<num_rows() && j<num_cols(), "SparseMatBaseCL (): index out of bounds", DebugNumericC);
    const size_t *pos= std::lower_bound( GetFirstCol(i), GetFirstCol(i+1), j);
    // lower_bound returns the iterator to the next column entry, if col j is not found
    return (pos != GetFirstCol(i+1) && *pos==j) ? _val[pos-GetFirstCol(0)] : T();
}


//**********************************************************************************
//
//  S p a r s e M a t D i a g C L :  stores location of diagonal of a sparse matrix
//
//**********************************************************************************

class SparseMatDiagCL
{
  private:
    std::valarray<size_t> _diagpos;

  public:
    template <typename T> SparseMatDiagCL (const SparseMatBaseCL<T>& A)
      : _diagpos(A.num_rows())
    {
        const size_t n=A.num_rows();

        for (size_t i=0; i<n; ++i)
            _diagpos[i]=std::lower_bound( A.GetFirstCol(i), A.GetFirstCol(i+1), i) - A.GetFirstCol(0);
    }

    size_t operator[] (size_t i) const { return _diagpos[i]; }
};


//*****************************************************************************
//
//  I/O for vectors and matrices
//
//*****************************************************************************

// Human readable output
template <typename T>
std::ostream& operator << (std::ostream& os, const VectorBaseCL<T>& v)
{
    os << v.size() << "      ";
    for (size_t i=0; i<v.size(); ++i) os << v[i] << ' ';
    return os << std::endl;
}


// Readable by in
template <typename T>
void out (std::ostream& os, const VectorBaseCL<T>& v)
{
    os << v.size() << '\n';
    for (size_t i=0; i<v.size(); ++i) os << v[i] << ' ';
    os << std::endl;
}


// Read vector from a stream
template <typename T>
void in (std::istream& is, VectorBaseCL<T>& v)
{
    size_t s;
    is >> s;
    v.resize(s);
    for (size_t i=0; i<s; ++i) is >> v[i];
}


// Human/Matlab readable output
template <typename T>
std::ostream& operator << (std::ostream& os, const SparseMatBaseCL<T>& A)
{
    const size_t M = A.num_rows();

    os << "% " << M << 'x' << A.num_cols() << ' ' << A.num_nonzeros() << " nonzeros\n";

    for (size_t row=0; row<M; ++row)
        for (size_t col=A.row_beg(row), rowend=A.row_beg(row+1); col<rowend; ++col)
            os << row+1 << ' ' << A.col_ind(col)+1 << ' ' << A.val(col) << '\n';

    return os << std::flush;
}


// Can be read by "in", see below
template <typename T>
void out (std::ostream& os, const SparseMatBaseCL<T>& A)
{
    os << A.num_rows() << ' ' << A.num_cols() << ' ' << A.num_nonzeros() << '\n';
    for (size_t row=0; row<=A.num_rows(); ++row) os << A.row_beg(row) << ' ';
    os << '\n';
    for (size_t nz=0; nz<A.num_nonzeros(); ++nz) os << A.col_ind(nz) << ' ';
    os << '\n';
    for (size_t nz=0; nz<A.num_nonzeros(); ++nz) os << A.val(nz) << ' ';
    os << std::endl;
}


// Read a sparse matrix from a stream
template <typename T>
void in (std::istream& is, SparseMatBaseCL<T>& A)
{
    size_t numrows, numcols, numnz;
    is >> numrows >> numcols >> numnz;
    A.resize(numrows, numcols, numnz);
    for (size_t row=0; row<=numrows; ++row) is >> A.row_beg(row);
    for (size_t nz=0; nz<numnz; ++nz) is >> A.col_ind(nz);
    for (size_t nz=0; nz<numnz; ++nz) is >> A.val(nz);
}


//*****************************************************************************
//
//  Matrix- and vector-operations
//
//*****************************************************************************

template <typename T>
SparseMatBaseCL<T>& SparseMatBaseCL<T>::LinComb (double coeffA, const SparseMatBaseCL<T>& A,
                                                 double coeffB, const SparseMatBaseCL<T>& B)
{
    Assert( A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(),
            "LinComb: incompatible dimensions", DebugNumericC);

    _rows=A.num_rows();
    _cols=A.num_cols();

    // The new sparsity pattern is computed by merging the lists of _colind (removing the duplicates) row by row.

    // Calculate the entries of _rowbeg (we need the number of nonzeros and get the rest for free)
    _rowbeg.resize(A.num_rows()+1);
    size_t i=0, iA=0, iB=0;

    for (size_t row=1; row<=A.num_rows(); ++row)                           // for each row
    {
        while ( iA < A.row_beg(row) )                                      // process every entry in A
        {
            while ( iB < B.row_beg(row) && B.col_ind(iB) < A.col_ind(iA) ) // process entries in B with smaller col_ind
                { ++i; ++iB; }
            if ( iB < B.row_beg(row) && B.col_ind(iB) == A.col_ind(iA) )   // process entries in B with equal col_ind
                ++iB;
            ++i; ++iA;
        }
        while ( iB < B.row_beg(row) )                                      // process, what is left in B
            { ++iB; ++i; }
        _rowbeg[row]=i;                                                    // store result
    }

    // Calculate the entries of _colind, _val (really perform the merging of the matrices)
    _colind.resize(row_beg(num_rows()));
    _val.resize(row_beg(num_rows()));
    i=0, iA=0, iB=0;

    for (size_t row=1; row<=A.num_rows(); ++row) // same algorithm as above
    {
        while ( iA < A.row_beg(row) )
        {
            while ( iB < B.row_beg(row) && B.col_ind(iB) < A.col_ind(iA) )
            {
                _val[i]=coeffB*B._val[iB];
                _colind[i++]=B._colind[iB++];
            }
            if ( iB < B.row_beg(row) && B.col_ind(iB) == A.col_ind(iA) )
                _val[i]=coeffA*A._val[iA]+coeffB*B._val[iB++];
            else
                _val[i]=coeffA*A._val[iA];
            _colind[i++]=A._colind[iA++];
        }
        while ( iB < B.row_beg(row) )
        {
            _val[i]=coeffB*B._val[iB];
            _colind[i++]=B._colind[iB++];
        }
    }

    return *this;
}


// y= A*x
// fails, if num_rows==0.
// Assumes, that none of the arrays involved do alias.
template <typename T>
inline void
y_Ax(T* __restrict y,
     size_t num_rows, 
     const T* __restrict Aval,
     const size_t* __restrict Arow,
     const size_t* __restrict Acol,
     const T* __restrict x)
{
    T sum;
    size_t rowend;
    size_t nz= 0;
    do {
        rowend= *++Arow;
        sum= T();
        for (; nz<rowend; ++nz)
            sum+= (*Aval++)*x[*Acol++];
        (*y++)= sum;
    } while (--num_rows > 0);
}


template <typename _MatEntry, typename _VecEntry>
VectorBaseCL<_VecEntry> operator * (const SparseMatBaseCL<_MatEntry>& A, const VectorBaseCL<_VecEntry>& x)
{
    VectorBaseCL<_VecEntry> ret( A.num_rows());
    Assert( A.num_cols()==x.size(), "SparseMatBaseCL * VectorBaseCL: incompatible dimensions", DebugNumericC);
    y_Ax( &ret[0],
          A.num_rows(),
          A.raw_val(),
          A.raw_row(),
          A.raw_col(),
          Addr( x));
    return ret;
}


// y+= A^T*x
// fails, if num_rows==0.
// Assumes, that none of the arrays involved do alias.
template <typename T>
inline void
y_ATx(T* __restrict y,
     size_t num_rows, 
     const T* __restrict Aval,
     const size_t* __restrict Arow,
     const size_t* __restrict Acol,
     const T* __restrict x)
{
    size_t rowend;
    size_t nz= 0;
    do {
        rowend= *++Arow;
        const T xrow= *x++;
        for (; nz<rowend; ++nz)
            y[*Acol++]+= (*Aval++)*xrow;
    } while (--num_rows > 0);
}

template <typename _MatEntry, typename _VecEntry>
VectorBaseCL<_VecEntry> transp_mul (const SparseMatBaseCL<_MatEntry>& A, const VectorBaseCL<_VecEntry>& x)
{
    VectorBaseCL<_VecEntry> ret( A.num_cols());
    Assert( A.num_rows()==x.size(), "transp_mul: incompatible dimensions", DebugNumericC);
    y_ATx( &ret[0],
           A.num_rows(),
           A.raw_val(),
           A.raw_row(),
           A.raw_col(),
           Addr( x));
    return ret;
}


//=============================================================================
//  Typedefs
//=============================================================================

typedef VectorBaseCL<double>       VectorCL;
typedef SparseMatBaseCL<double>    MatrixCL;
typedef SparseMatBuilderCL<double> MatrixBuilderCL;

} // end of namespace DROPS

#endif
