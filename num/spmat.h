//**************************************************************************
// File:     spmat.h                                                       *
// Content:  sparse matrix in compressed row format                        *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.2                                                           *
// Remarks:  developped from TNT fcscmat.h                                 *
//**************************************************************************


#ifndef DROPS_SPMAT_H
#define DROPS_SPMAT_H

#include <iostream>
#include <valarray>
#include <map>

namespace DROPS
{

template <class T>
class SparseMatBuilderCL;

/**************************************************************************
*
*  S p a r s e M a t B a s e C L :  row based sparse matrix format
*                                   hint: use SparseMatBuilderCL for setting up
*
**************************************************************************/

template <class T>
class SparseMatBaseCL
{
public:
    typedef T      value_type;
    typedef size_t subs_type;            // type of subscripts
private:
    subs_type _rows;                     // global dimensions
    subs_type _cols;

    std::valarray<subs_type> _rowbeg;    // (_rows+1 entries, last entry must be <=_nz) index of first non-zero-entry in _val belonging to the row given as subscript
    std::valarray<subs_type> _colind;    // (_nz elements) column-number of correspnding entry in _val
    std::valarray<T>         _val;       // nonzero-entries

private:
    void resize (subs_type rows, subs_type cols, subs_type nz)
        { _rows=rows; _cols=cols; _rowbeg.resize(rows+1); _colind.resize(nz); _val.resize(nz); }
    void resize_rows (subs_type rows)
        { _rows=rows; _rowbeg.resize(rows+1); }
    void resize_cols (subs_type cols, subs_type nz)
        { _cols=cols; _colind.resize(nz); }
    void resize_val (subs_type nz)
        {  _val.resize(nz); }
public:
    SparseMatBaseCL ()   {}
    // default copy-ctor, assignment-op, dtor
    SparseMatBaseCL (subs_type rows, subs_type cols, subs_type nz)
        : _rows(rows), _cols(cols), _rowbeg(rows+1), _colind(nz), _val(nz) {}
    SparseMatBaseCL (subs_type rows, subs_type cols, subs_type nz,
                 const T* valbeg , const subs_type* rowbegbeg, const subs_type* colindbeg)
        : _rows(rows), _cols(cols), _rowbeg(rowbeg, rows+1), _colind(colindbeg, nz), _val(valbeg, nz) {}

    subs_type num_rows     () const { return _rows; }
    subs_type num_cols     () const { return _cols; }
    subs_type num_nonzeros () const { return _val.size(); }

//    subs_type& row_beg (subs_type i)       { return _rowbeg[i]; }
    subs_type  row_beg (subs_type i) const { return _rowbeg[i]; }
//    subs_type& col_ind (subs_type i)       { return _colind[i]; }
    subs_type  col_ind (subs_type i) const { return _colind[i]; }
    subs_type* GetFirstCol(subs_type i) const
              { return &(*const_cast<std::valarray<subs_type>*>(&_colind))[_rowbeg[i]]; }
    value_type* GetFirstVal(subs_type i) const
              { return &(*const_cast<std::valarray<value_type>*>(&_val))[_rowbeg[i]]; }
//    T&         val     (subs_type i)       { return _val[i]; }
    T          val     (subs_type i) const { return _val[i]; }
    inline T&   operator() (subs_type i, subs_type j);
    inline T    operator() (subs_type i, subs_type j) const;
// deprecated
//    inline T    get        (subs_type i, subs_type j) const;
    SparseMatBaseCL& operator*= (T coeff) { _val*= coeff; return *this; }
    void clear() { resize(0,0,0); }
    friend class SparseMatBuilderCL<T>;
};

//TODO: fix following access operators
template <class T>
T& SparseMatBaseCL<T>::operator() (subs_type i, subs_type j)
{
    Assert(i<num_rows() && j<num_cols(), "SparseMatBaseCL operator(): index out of bounds", DebugNumericC);
    subs_type *pos= std::lower_bound( &_colind[_rowbeg[i]], &_colind[_rowbeg[i+1]], j);
    // perhaps column entry j doesn't exist in row i
    // N o t e :  lower_bound returns the iterator to the next column entry, if col j is not found !!
    Assert(pos!=&_colind[_rowbeg[i+1]] && *pos==j, DROPSErrCL("SparseMatBaseCL operator(): invalid index"), DebugNumericC);   
    return _val[pos - &_colind[0]];
}


// std::find does not work, because we cannot obtain the address of the first element in an
// valarray, if we are not allowed to invoke non-const members ( operator[] const only returns a value...)
template <class T>
T SparseMatBaseCL<T>::operator() (subs_type i, subs_type j) const
{
    Assert(i<num_rows() && j<num_cols(), "SparseMatBaseCL operator(): index out of bounds", DebugNumericC);
    subs_type *pos= std::lower_bound( GetFirstCol(i), GetFirstCol(i+1), j);
    // N o t e :  lower_bound returns the iterator to the next column entry, if col j is not found !!
    if (pos != GetFirstCol(i+1) && *pos==j) return _val[pos-GetFirstCol(0)];
    return 0;
}

// std::find does not work, because we cannot obtain the address of the first element in an
// valarray, if we are not allowed to invoke non-const members ( operator[] const only returns a value...)
/*
template <class T>
T SparseMatBaseCL<T>::get (subs_type i, subs_type j) const
{
    Assert(i<num_rows() && j<num_cols(), "SparseMatBaseCL get(): index out of bounds", DebugNumericC);
    for(subs_type k= row_beg(i); k<row_beg(i+1); ++k)
        if (col_ind(k)==j)
            return val(k);
    return 0;
}
*/
template <class T>
void ConvexComb( SparseMatBaseCL<T>& M, const double coeffA, const SparseMatBaseCL<T>& A, 
                                        const double coeffB, const SparseMatBaseCL<T>& B)
{
    Assert( A.num_rows()==B.num_rows && A.num_cols()==B.num_cols(), 
            DROPSErrCL("SparseMatBaseCL(double,Mat,double,Mat): Matrices have different size!"), DebugNumericC);
    SparseMatBuilderCL<double> MB( &M, A.num_rows(), A.num_cols());
    size_t col, indA= 0, indB= 0;

    for (size_t row= 0, row_end= A.num_rows(); row<row_end; ++row)
    {
        while ( indA < A.row_beg(row+1) )
        {
            col= A.col_ind( indA++);
            MB(row,col)+= coeffA * A(row,col);
        }
        while ( indB < B.row_beg(row+1) )
        {
            col= B.col_ind( indB++);
            MB(row,col)+= coeffB * B(row,col);
        }
    }
    MB.Build();
}

// Human readable output, MatLab readable as well
template <class T>
std::ostream& operator << (std::ostream &os, const SparseMatBaseCL<T> &A)
{
    typedef typename SparseMatBaseCL<T>::subs_type subs_t;

    const subs_t M = A.num_rows();
    const subs_t N = A.num_cols();

    os << "% " << M << 'x' << N << ' ' << A.num_nonzeros() << " nonzeros" << std::endl;

    for (subs_t row=0; row<M; ++row)
    {
        const subs_t rowbegin = A.row_beg(row);
        const subs_t rowend   = A.row_beg(row+1);
        for (subs_t col=rowbegin; col<rowend; ++col)
            os << row+1 << ' ' << A.col_ind(col)+1 << ' ' << A.val(col) << '\n';
        os.flush();
    }

    return os;
}


// Can be read by "in", see below
template <class T>
void out( std::ostream &os, const SparseMatBaseCL<T> &A)
{
   typedef typename SparseMatBaseCL<T>::subs_type subs_t;
   
   os << A.num_rows() << " " << A.num_cols() << " " << A.num_nonzeros() << std::endl;
   for(subs_t row=0; row<=A.num_rows(); ++row) os << A.row_beg(row) << " ";
   os << "\n";
   for(subs_t nz=0; nz<A.num_nonzeros(); ++nz) os << A.col_ind(nz) << " ";
   os << "\n";
   for(subs_t nz=0; nz<A.num_nonzeros(); ++nz) os << A.val(nz) << " ";
   os << std::endl;
}


// Read a sparse matrix from a stream
template <class T>
void in( std::istream &is, SparseMatBaseCL<T> &A)
{
   typedef typename SparseMatBaseCL<T>::subs_type subs_t;
   subs_t numrows, numcols, numnz;
   is >> numrows >> numcols >> numnz;
   A.resize(numrows, numcols, numnz);
   for(subs_t row=0; row<=numrows; ++row) is >> A.row_beg(row);
   for(subs_t nz=0; nz<numnz; ++nz) is >> A.col_ind(nz);
   for(subs_t nz=0; nz<numnz; ++nz) is >> A.val(nz);
}


/**************************************************************************
*
*  S p a r s e M a t B u i l d e r C L :  used for setting up sparse matrices
*
**************************************************************************/

template <class T>
class SparseMatBuilderCL
{
public:
    typedef T                       valueT;
    typedef size_t                  subsT;            // type of subscripts
    typedef SparseMatBaseCL<T>      spmatT;
    typedef std::pair<subsT,valueT> entryT;
    typedef std::map<subsT,valueT>  couplT;

private:
    Uint    _rows;
    Uint    _cols;
    spmatT* _mat;
    couplT* _coupl;
    
public:
    SparseMatBuilderCL(spmatT* mat, subsT rows, subsT cols)
        : _rows(rows), _cols(cols), _mat(mat), _coupl(new couplT[_rows])
        { 
            Assert( _coupl!=0, DROPSErrCL("SparseMatBuilderCL: not enough memory!"), -1); 
        }
    ~SparseMatBuilderCL() { delete[] _coupl; }

    inline T&   operator() (subsT i, subsT j);
    void Build();
};

template<class T>
T& SparseMatBuilderCL<T>::operator() (subsT i, subsT j)
{
    Assert(i<_rows && j<_cols, DROPSErrCL("SparseMatBuilderCL operator(): index out of bounds"), DebugNumericC);
    return _coupl[i][j];
}

template<class T>
void SparseMatBuilderCL<T>::Build()
{
    Uint nz= 0;
    for (Uint i= 0; i<_rows; ++i)
        for (typename couplT::const_iterator it= _coupl[i].begin(), end= _coupl[i].end();
             it != end; ++it)
            if(it->second!=0) //nz+= _coupl[i].size();
                ++nz;
                    
    _mat->resize(_rows, _cols, nz);

    nz= 0;
    for (Uint i=0; i<_rows; ++i)
    {
        _mat->_rowbeg[i]= nz;
        for (typename couplT::const_iterator it= _coupl[i].begin(), end= _coupl[i].end();
             it != end; ++it)
            if (it->second != 0)
            {
                _mat->_colind[nz]= it->first;
                _mat->_val[nz]=    it->second;
                ++nz;
            }
        // the col_ind-entries in each row are sorted, as they were stored sorted in the map
   }
   _mat->_rowbeg[_rows]= nz;

   Assert( nz == _mat->num_nonzeros(), 
           DROPSErrCL("SparseMatBuilderCL: wrong counting of nonzeros in Build()"), -1);
              
   delete[] _coupl;
   _coupl= 0;
}


/**************************************************************************
*
*  V e c t o r B a s e C L :   base class for numerical vectors
*                              based on std::valarray
*
**************************************************************************/

template <class T>
class VectorBaseCL
{
private:
    std::valarray<T> _va;

public:
    typedef T      value_type;
    typedef size_t subs_type;            // type of subscripts

    // ctors
    VectorBaseCL ()                                 : _va()       {}
    VectorBaseCL (subs_type s)                      : _va(s)      {}
    VectorBaseCL (const T& val, subs_type s)        : _va(val, s) {}
    VectorBaseCL (const T* tp, subs_type s)         : _va(tp, s)  {}
    VectorBaseCL (const std::valarray<T>& va)       : _va(va)     {}
    VectorBaseCL (const VectorBaseCL& v)            : _va(v._va)  {}

    VectorBaseCL (const std::slice_array<T>& sla)   : _va(sla)    {}
    VectorBaseCL (const std::gslice_array<T>& gsla) : _va(gsla)   {}
    VectorBaseCL (const std::mask_array<T>& ma)     : _va(ma)     {}
    VectorBaseCL (const std::indirect_array<T>& ia) : _va(ia)     {}

    void resize (subs_type s, T val = T()) { _va.resize(s, val); }

    // element access
    T  operator [] (subs_type s) const { Assert( s<size(), "VectorCL: out of bounds", DebugNumericC);
                                         return _va[s]; }
    T& operator [] (subs_type s)       { Assert( s<size(), "VectorCL: out of bounds", DebugNumericC);
                                         return _va[s]; }

    // assignment
    inline VectorBaseCL<T>& operator= (const std::valarray<T>& va);
    inline VectorBaseCL<T>& operator= (const VectorBaseCL<T>& v);
    VectorBaseCL<T>& operator= (const T& val)                     { _va = val;  return *this; }
    VectorBaseCL<T>& operator= (const std::slice_array<T>& sla)   { _va = sla;  return *this; }
    VectorBaseCL<T>& operator= (const std::gslice_array<T>& gsla) { _va = gsla; return *this; }
    VectorBaseCL<T>& operator= (const std::mask_array<T>& ma)     { _va = ma;   return *this; }
    VectorBaseCL<T>& operator= (const std::indirect_array<T>& ia) { _va = ia;   return *this; }

    // computed assignment
    VectorBaseCL<T>& operator+= (const T& val)         { _va +=  val; return *this; }
    VectorBaseCL<T>& operator-= (const T& val)         { _va -=  val; return *this; }
    VectorBaseCL<T>& operator*= (const T& val)         { _va *=  val; return *this; }
    VectorBaseCL<T>& operator/= (const T& val)         { _va /=  val; return *this; }
    inline VectorBaseCL<T>& operator+= (const VectorBaseCL<T>& v);
    inline VectorBaseCL<T>& operator-= (const VectorBaseCL<T>& v);
    inline VectorBaseCL<T>& operator*= (const VectorBaseCL<T>& v);
    inline VectorBaseCL<T>& operator/= (const VectorBaseCL<T>& v);

    // unary minus
    VectorBaseCL<T> operator - () const
        { VectorBaseCL<T> ret(0.0, this->size()); ret-=*this; return ret; }

    // member functions
    subs_type size () const { return _va.size(); }
    T sum     () const { return _va.sum(); }
    T min     () const { return _va.min(); }
    T max     () const { return _va.max(); }
    T norm2   () const { return (*this)*(*this); }
    T norm    () const { return sqrt(norm2()); }
    T supnorm () const;
};


template<typename T>
inline VectorBaseCL<T>&
VectorBaseCL<T>::operator=(const std::valarray<T>& va)
{
    Assert(_va.size()==va.size(), DROPSErrCL("Vec= valarray: incompatible dimensions!"), DebugNumericC);
    _va = va;   return *this;
}

template<typename T>
inline VectorBaseCL<T>& 
VectorBaseCL<T>::operator=(const VectorBaseCL<T>& v)
{
    Assert(_va.size()==v.size(), DROPSErrCL("Vec= Vec: incompatible dimensions!"), DebugNumericC);
    _va = v._va;   return *this;
}

template<typename T>
inline VectorBaseCL<T>& 
VectorBaseCL<T>::operator+=(const VectorBaseCL<T>& v)
{
    Assert(_va.size()==v.size(), DROPSErrCL("Vec+= Vec: incompatible dimensions!"), DebugNumericC);
    _va += v._va; return *this;
}

template<typename T>
inline VectorBaseCL<T>&
VectorBaseCL<T>::operator-=(const VectorBaseCL<T>& v)
{
    Assert(_va.size()==v.size(), DROPSErrCL("Vec-= Vec: incompatible dimensions!"), DebugNumericC);
    _va -= v._va; return *this;
}

template<typename T>
inline VectorBaseCL<T>&
VectorBaseCL<T>::operator*=(const VectorBaseCL<T>& v)
{
    Assert(_va.size()==v.size(), DROPSErrCL("Vec*= Vec: incompatible dimensions!"), DebugNumericC);
    _va *= v._va; return *this;
}

template<typename T>
inline VectorBaseCL<T>&
VectorBaseCL<T>::operator/=(const VectorBaseCL<T>& v)
{
    Assert(_va.size()==v.size(), DROPSErrCL("Vec/= Vec: incompatible dimensions!"), DebugNumericC);
    _va /= v._va; return *this;
}

template<typename T>
inline T VectorBaseCL<T>::supnorm () const
{
    T ret=0, tmp;
    for (subs_type i=0; i<size(); ++i) if ((tmp=fabs(_va[i]))>ret) ret= tmp;
    return ret;
}


template<typename T>
inline VectorBaseCL<T> operator + (const VectorBaseCL<T> &v, const VectorBaseCL<T> &w)
{
    Assert( v.size()==w.size(), DROPSErrCL("vec + vec: incompatible dimensions!"), DebugNumericC);
    VectorBaseCL<T> ret(v.size());
    for (typename VectorBaseCL<T>::subs_type i=0; i<v.size(); ++i) ret[i] = v[i]+w[i];
    return ret;
}


template<typename T>
inline VectorBaseCL<T> operator - (const VectorBaseCL<T> &v, const VectorBaseCL<T> &w)
{
    Assert( v.size()==w.size(), DROPSErrCL("vec - vec: incompatible dimensions!"), DebugNumericC);
    VectorBaseCL<T> ret(v.size());
    for (typename VectorBaseCL<T>::subs_type i=0; i<v.size(); ++i) ret[i] = v[i]-w[i];
    return ret;
}


template<typename T>
inline T operator * (const VectorBaseCL<T> &v, const VectorBaseCL<T> &w)
{
    Assert( v.size()==w.size(), DROPSErrCL("vec * vec: incompatible dimensions!"), DebugNumericC);
    T sum= T();
    for(Uint i=0, len= v.size(); i<len; ++i) sum+= v[i]*w[i];
    return sum;
}


template<typename T>
inline VectorBaseCL<T> operator * (T c, const VectorBaseCL<T> &v)
{
    VectorBaseCL<T> ret(v.size());
    for (typename VectorBaseCL<T>::subs_type i=0; i<v.size(); ++i) ret[i] = c*v[i];
    return ret;
}

// z_xpaypby2(z, x, a,y,b,y2) := z= x + a*y + b*y2
template<typename T>
inline void z_xpaypby2(VectorBaseCL<T>& z,
                       const VectorBaseCL<T>& x,
                       T a, const VectorBaseCL<T>& y,
                       T b, const VectorBaseCL<T>& y2)
{
    Assert(z.size()==x.size() && z.size() == y.size() && z.size()==y2.size(),
           DROPSErrCL("z_xpaypy2: incompatible dimensions!"), DebugNumericC);

    const typename VectorBaseCL<T>::subs_type s= z.size();
    for (typename VectorBaseCL<T>::subs_type i=0; i<s; ++i)
        z[i]= x[i] + a*y[i] + b*y2[i];
}

// z_xpay(z,x,a,y) := z= x + a*y
template<typename T>
inline void z_xpay(VectorBaseCL<T>& z, const VectorBaseCL<T>& x, T a, const VectorBaseCL<T>& y)
{
    Assert( z.size()==x.size() && z.size() == y.size(), DROPSErrCL("z_xpay: incompatible dimensions!"), DebugNumericC);

    const typename VectorBaseCL<T>::subs_type s= z.size();
    for (typename VectorBaseCL<T>::subs_type i=0; i<s; ++i)
        z[i]= x[i] + a*y[i];
}

// axpy(a,x,y) :=  y+= a*x
template<typename T>
inline void axpy(const T a, const VectorBaseCL<T>& x, VectorBaseCL<T>& y)
{
    Assert( x.size()==y.size(), DROPSErrCL("axpy: incompatible dimensions!"), DebugNumericC);

    const typename VectorBaseCL<T>::subs_type s= x.size();
    for (typename VectorBaseCL<T>::subs_type i=0; i<s; ++i) y[i]+= a*x[i];
}

// Human readable output
template <class T>
std::ostream& operator << (std::ostream &os, const VectorBaseCL<T> &v)
{
    const typename VectorBaseCL<T>::subs_type s = v.size();

    os << s << "      ";
    for (typename VectorBaseCL<T>::subs_type i=0; i<s; ++i)
        os << v[i] << ' ';
    return os << std::endl;
}


// Readable by in
template <class T>
void out( std::ostream &os, const VectorBaseCL<T> &v)
{
   typedef typename VectorBaseCL<T>::subs_type subs_t;
   
   os << v.size() << std::endl;
   for(subs_t i=0; i<v.size(); ++i) os << v[i] << " ";
   os << std::endl;
}


// Read vector from a stream
template <class T>
void in( std::istream &is, VectorBaseCL<T> &v)
{
   typedef typename VectorBaseCL<T>::subs_type subs_t;
   subs_t s;
   is >> s;
   v.resize(s);
   for(subs_t i=0; i<s; ++i) is >> v[i];
}


template <class _MatEntry, class _VecEntry>
VectorBaseCL<_VecEntry> operator * (const SparseMatBaseCL<_MatEntry>& A, const VectorBaseCL<_VecEntry>& x)
{
    typedef typename SparseMatBaseCL<_MatEntry>::subs_type subs_t;
    const subs_t M(A.num_rows());
    VectorBaseCL<_VecEntry> ret(M);

    Assert( A.num_cols()==x.size(), DROPSErrCL("SparseMat * Vec: incompatible dimensions!"), DebugNumericC);
    for (subs_t row=0, nz=0; row<M; ++row)
    {
        const subs_t rowend(A.row_beg(row+1));
        _VecEntry sum= _VecEntry();
        for (; nz<rowend; ++nz)
            sum += A.val(nz)*x[A.col_ind(nz)];
        ret[row]= sum;
    }
    return ret;
}


template <class _MatEntry, class _VecEntry>
VectorBaseCL<_VecEntry> transp_mul (const SparseMatBaseCL<_MatEntry>& A, const VectorBaseCL<_VecEntry>& x)
{
    typedef typename SparseMatBaseCL<_MatEntry>::subs_type subs_t;
    const subs_t num_rows= A.num_rows(),
                 num_cols= A.num_cols();
    VectorBaseCL<_VecEntry> ret(num_cols);

    Assert( num_rows==x.size(), DROPSErrCL("transp_mul: incompatible dimensions!"), DebugNumericC);
    for (subs_t row=0, nz=0; row<num_rows; ++row)
    {
        const subs_t rowend(A.row_beg(row+1));
        const _VecEntry xrow= x[row];
        for (; nz<rowend; ++nz)
        {
            Assert( A.col_ind(nz)<num_cols, DROPSErrCL("transp_mul: out of bounds!"), DebugNumericC);
            ret[A.col_ind(nz)]+= A.val(nz)*xrow;
        }
    }
    return ret;
}


} // end of namespace DROPS

#endif
