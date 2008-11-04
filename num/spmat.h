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
#include <vector>
#include <deque>
#include <numeric>
#include <limits>
#if __GNUC__ >= 4 && !defined(__INTEL_COMPILER) 
#    include <tr1/unordered_map>
#else
#    include <map>
#endif
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
  public:
    typedef T value_type;
    typedef std::valarray<T> base_type;

    // ctors
    VectorBaseCL()                      : base_type()       {}
#ifdef VALARRAY_BUG
    VectorBaseCL (size_t s)             : base_type( T(),s) {}
#else
    VectorBaseCL (size_t s)             : base_type( s)     {}
#endif
    VectorBaseCL (T c, size_t s)        : base_type( c, s)  {}
    VectorBaseCL (const T* tp, size_t s): base_type( tp, s) {}

DROPS_DEFINE_VALARRAY_DERIVATIVE( VectorBaseCL, T, base_type)

#if (DROPSDebugC & DebugNumericC)
    using base_type::operator[];

    // For checking new code; the following functions interfere with the expression
    // template mechanism of std::valarray, but do some checks before the
    // operations.
    // element access
    T  operator[](size_t s) const;
    T& operator[](size_t s);
#endif
};


#if (DROPSDebugC & DebugNumericC)
template <typename T>
T  VectorBaseCL<T>::operator[](size_t s) const
{
    Assert(s<base_type::size(), "VectorBaseCL []: index out of bounds", DebugNumericC);
    return (*static_cast<const base_type*>( this))[s];
}

template <typename T>
T& VectorBaseCL<T>::operator[](size_t s)
{
    Assert(s<base_type::size(), "VectorBaseCL []: index out of bounds", DebugNumericC);
    return (*static_cast<base_type*>( this))[s];
}
#endif


// Get the address of the first element in a valarray
// ("&x[0]" doesn't work, because "operator[] const" only returns a value)
template <typename T>
  inline const T*
  Addr(const std::valarray<T>& x)
{
    return &(const_cast<std::valarray<T>&>(x)[0]);
}

template <typename T>
  inline T*
  Addr(std::valarray<T>& x)
{
    return &(x[0]);
}

template <class T>
  inline T
  dot(const VectorBaseCL<T>& v, const VectorBaseCL<T>& w)
{
    Assert( v.size()==w.size(), "dot: incompatible dimensions", DebugNumericC);
    return std::inner_product( Addr( v), Addr( v) + v.size(), Addr( w), T());
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


template <typename T>
  void
  permute_Vector (VectorBaseCL<T>& v, const PermutationT& p)
{
    Assert( v.size() == p.size(),
        DROPSErrCL( "permute_Vector: v and p have different dimension.\n"), DebugNumericC);

    VectorBaseCL<T> w( v);
    for (size_t i= 0; i < v.size(); ++i)
        v[p[i]]= w[i];
}


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
#if __GNUC__ >= 4 && !defined(__INTEL_COMPILER) 
    typedef std::tr1::unordered_map<size_t,valueT> couplT;
#else
    typedef std::map<size_t,valueT> couplT;
#endif

private:
    size_t  _rows;
    size_t  _cols;
    spmatT* _mat;
    bool    _reuse;
    couplT* _coupl;

public:
    SparseMatBuilderCL(spmatT* mat, size_t rows, size_t cols)
        : _rows(rows), _cols(cols), _mat(mat),
          _reuse(!(DROPSDebugC & DebugNoReuseSparseC) && mat->num_nonzeros()!=0
                 && mat->num_rows()==rows && mat->num_cols()==cols)
    {
        mat->IncrementVersion();
        if (_reuse)
        {
            Comment("SparseMatBuilderCL: Reusing OLD matrix" << std::endl, DebugNumericC);
            _coupl=0;
            _mat->_val=T();
        }
        else
        {
            Comment("SparseMatBuilderCL: Creating NEW matrix" << std::endl, DebugNumericC);
            _coupl=new couplT[_rows];
            Assert( _coupl!=0, "SparseMatBuilderCL: out of memory", ~0);
        }
    }

    ~SparseMatBuilderCL() { delete[] _coupl; }

    T& operator() (size_t i, size_t j)
    {
        Assert(i<_rows && j<_cols, "SparseMatBuilderCL (): index out of bounds", DebugNumericC);

        if (_reuse)
        {
            // search row for the correct entry
            const size_t rowend=_mat->_rowbeg[i+1]-1;

            for (size_t k=_mat->_rowbeg[i]; k<rowend; ++k)
                if (_mat->_colind[k]==j)
                    return _mat->_val[k];

            Assert(_mat->_colind[rowend]==j, "SparseMatBuilderCL (): no such index", ~0);
            return _mat->_val[rowend];
        }
        else
            return _coupl[i][j];
    }

    void Build();
};

template <typename T>
void SparseMatBuilderCL<T>::Build()
{
    if (_reuse) return;

    size_t nz= 0;
    for (size_t i= 0; i < _rows; ++i)
        nz+= _coupl[i].size();

    _mat->resize( _rows, _cols, nz);

    nz= 0;
#if __GNUC__ >= 4 && !defined(__INTEL_COMPILER) 
    typedef std::pair<size_t, T> PT;
    std::vector<PT> pv;
    for (size_t i= 0; i < _rows; ++i) {
        pv.resize( _coupl[i].size());
        _mat->_rowbeg[i]= nz;
        nz+= pv.size();
        std::copy( _coupl[i].begin(), _coupl[i].end(), pv.begin());
        // The col_ind-entries in each row are sorted.
        std::sort( pv.begin(), pv.end(), less1st<PT>());
        for (size_t k= _mat->_rowbeg[i], j= 0; k < nz; ++k, ++j) {
            _mat->_colind[k]= pv[j].first;
            _mat->_val[k]= pv[j].second;
        // std::cout << _coupl[i].load_factor() << '\t' << std::setfill('0') << std::setw(3) <<_coupl[i].size() << '\n';
        }
    }
#else
    for (size_t i=0; i<_rows; ++i)
    {
        _mat->_rowbeg[i]= nz;
        for (typename couplT::const_iterator it= _coupl[i].begin(), end= _coupl[i].end(); it != end; ++it)
        {
            _mat->_colind[nz]= it->first;
            _mat->_val[nz]=    it->second;
            ++nz;
        }
        // the col_ind-entries in each row are sorted, as they were stored sorted in the map
    }
#endif
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

    size_t version_;                // All modifications increment this. Starts with 1.

    std::valarray<size_t> _rowbeg;  // (_rows+1 entries, last entry must be <=_nz) index of first non-zero-entry in _val belonging to the row given as subscript
    std::valarray<size_t> _colind;  // (_nz elements) column-number of corresponding entry in _val
    std::valarray<T>      _val;     // nonzero-entries

    void resize_rows (size_t rows)            { _rows=rows; _rowbeg.resize(rows+1); }
    void resize_cols (size_t cols, size_t nz) { _cols=cols; _colind.resize(nz); }
    void resize_val  (size_t nz)              { _val.resize(nz); }

public:
    typedef T value_type;

    SparseMatBaseCL () : version_( 1) {}
    SparseMatBaseCL& operator= (const SparseMatBaseCL& m);    
    // default copy-ctor, dtor
    SparseMatBaseCL (size_t rows, size_t cols, size_t nz)
        : _rows(rows), _cols(cols), version_( 1), _rowbeg(rows+1), _colind(nz), _val(nz) {}
    SparseMatBaseCL (size_t rows, size_t cols, size_t nz,
                     const T* valbeg , const size_t* rowbeg, const size_t* colindbeg)
        : _rows(rows), _cols(cols), version_( 1), _rowbeg(rowbeg, rows+1), _colind(colindbeg, nz), _val(valbeg, nz) {}
    SparseMatBaseCL (const std::valarray<T>&); // Creates a square diagonal matrix.

    const T*      raw_val() const { return Addr( _val); }
    T*            raw_val()       { return &_val[0]; }
    const size_t* raw_row() const { return Addr( _rowbeg); }
    size_t*       raw_row()       { return &_rowbeg[0]; }
    const size_t* raw_col() const { return Addr( _colind); }
    size_t*       raw_col()       { return &_colind[0]; }

    std::valarray<T>&       val()       { return _val; }
    const std::valarray<T>& val() const { return _val; }

    size_t num_rows     () const { return _rows; }
    size_t num_cols     () const { return _cols; }
    size_t num_nonzeros () const { return _val.size(); }

    size_t row_beg (size_t i) const { return _rowbeg[i]; }
    size_t col_ind (size_t i) const { return _colind[i]; }
    T      val     (size_t i) const { return _val[i]; }

    void IncrementVersion() { ++version_; }
    size_t Version() const {  return version_; }

    const size_t* GetFirstCol(size_t i) const { return Addr(_colind)+_rowbeg[i]; }
    const T*      GetFirstVal(size_t i) const { return Addr(_val)+_rowbeg[i]; }

    inline T  operator() (size_t i, size_t j) const;

    SparseMatBaseCL& operator*= (T c) { IncrementVersion(); _val*= c; return *this; }
    SparseMatBaseCL& operator/= (T c) { IncrementVersion(); _val/= c; return *this; }

    SparseMatBaseCL& LinComb (double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&);
    SparseMatBaseCL& LinComb (double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&);

    void resize (size_t rows, size_t cols, size_t nz)
        { IncrementVersion(); _rows=rows; _cols=cols; _rowbeg.resize(rows+1); _colind.resize(nz); _val.resize(nz); }
    void clear() { resize(0,0,0); }

    VectorBaseCL<T> GetDiag() const;

    void permute_rows (const PermutationT&);
    void permute_columns (const PermutationT&);

    friend class SparseMatBuilderCL<T>;
};

template <typename T>
  SparseMatBaseCL<T>& SparseMatBaseCL<T>::operator= (const SparseMatBaseCL<T>& m)
{
    if (&m == this) return *this;

    IncrementVersion();
    _rows= m._rows;
    _cols= m._cols;
    _rowbeg.resize(m._rowbeg.size());
    _rowbeg= m._rowbeg;
    _colind.resize(m._colind.size());
    _colind= m._colind;
    _val.resize(m._val.size());
    _val= m._val;
    return *this;
}

template <typename T>
  SparseMatBaseCL<T>::SparseMatBaseCL(const std::valarray<T>& v)
      : _rows( v.size()), _cols( v.size()), version_( 1),
        _rowbeg( v.size() + 1), _colind( v.size()), _val( v)
{
    for (size_t i= 0; i < _rows; ++i)
        _rowbeg[i]= _colind[i]= i;
    _rowbeg[_rows]= _rows;
}

template <typename T>
T SparseMatBaseCL<T>::operator() (size_t i, size_t j) const
{
    Assert(i<num_rows() && j<num_cols(), "SparseMatBaseCL (): index out of bounds", DebugNumericC);
    const size_t *pos= std::lower_bound( GetFirstCol(i), GetFirstCol(i+1), j);
    // lower_bound returns the iterator to the next column entry, if col j is not found
    return (pos != GetFirstCol(i+1) && *pos==j) ? _val[pos-GetFirstCol(0)] : T();
}

template <typename T>
VectorBaseCL<T> SparseMatBaseCL<T>::GetDiag() const
{
    const size_t n=num_rows();
    Assert(n==num_cols(), "SparseMatBaseCL::GetDiag: no square Matrix", DebugParallelC);
    VectorBaseCL<T> diag(n);
    for (size_t i=0; i<n; ++i)
        diag[i] = (*this)(i,i);
    return diag;
}

template <typename T>
  void
  SparseMatBaseCL<T>::permute_rows (const PermutationT& p)
{
    Assert( num_rows() == p.size(),
        DROPSErrCL( "permute_rows: Matrix and Permutation have different dimension.\n"), DebugNumericC);

    IncrementVersion();
    PermutationT pi( invert_permutation( p));
    std::valarray<size_t> r( _rowbeg);
    std::valarray<size_t> c( _colind);
    std::valarray<T> v( _val);
    _rowbeg[0]= 0;
    for (size_t i= 0; i < num_rows(); ++i) {
        size_t nonzeros= r[pi[i] + 1] - r[pi[i]];
        _rowbeg[i + 1]= _rowbeg[i] + nonzeros;
        std::slice oldpos( std::slice( r[pi[i]], nonzeros, 1));
        std::slice newpos( std::slice( _rowbeg[i], nonzeros, 1));
        _val[newpos]= v[oldpos];
        _colind[newpos]= c[oldpos];
    }
}

template <typename T>
  void
  SparseMatBaseCL<T>::permute_columns (const PermutationT& p)
{
    Assert( num_cols() == p.size(),
        DROPSErrCL( "permute_columns: Matrix and Permutation have different dimension.\n"), DebugNumericC);

    IncrementVersion();
    for (size_t i= 0; i < _colind.size(); ++i)
        _colind[i]= p[_colind[i]];

    // Sort the indices in each row. This affects the internal matrix-layout only.
    typedef std::pair<size_t, T> PT;
    for (size_t r= 0; r < num_rows(); ++r) {
        std::vector<PT> pv( row_beg( r + 1) - row_beg( r));
        for (size_t i= row_beg( r), j= 0; i < row_beg( r + 1); ++i, ++j)
            pv[j]= std::make_pair( _colind[i], _val[i]);
        std::sort( pv.begin(), pv.end(), less1st<PT>());
        for (size_t i= row_beg( r), j= 0; i < row_beg( r + 1); ++i, ++j) {
            _colind[i]= pv[j].first;
            _val[i]= pv[j].second;
        }        
    }
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
//  2x2-Blockmatrix: ( A & B \\ C & D)
//
//*****************************************************************************
enum BlockMatrixOperationT { MUL, TRANSP_MUL };

template <class MatT>
class BlockMatrixBaseCL
{
  public:
    typedef BlockMatrixOperationT OperationT;

  private:
    const MatT* block_[4];
    OperationT operation_[4];

    bool block_num_rows(size_t b, size_t& nr) const;
    bool block_num_cols(size_t b, size_t& nc) const;

  public:
    BlockMatrixBaseCL( const MatT* A, OperationT Aop, const MatT* B, OperationT Bop,
        const MatT* C, OperationT Cop, const MatT* D= 0, OperationT Dop= MUL);

    size_t num_rows(size_t) const;
    size_t num_cols(size_t) const;
    size_t num_rows() const { return this->num_rows( 0) + this->num_rows( 1); }
    size_t num_cols() const { return this->num_cols( 0) + this->num_cols( 1); }

    const MatT* GetBlock( size_t b) const { return block_[b]; }

    OperationT GetOperation( size_t b) const { return operation_[b]; }
    OperationT GetTransposeOperation( size_t b) const {
        return operation_[b] == MUL ? TRANSP_MUL : MUL;
    }

    BlockMatrixBaseCL GetTranspose() const;

    VectorBaseCL<typename MatT::value_type> GetDiag() const;
};

template <class MatT>
BlockMatrixBaseCL<MatT>::BlockMatrixBaseCL( const MatT* A, OperationT Aop,
    const MatT* B, OperationT Bop, const MatT* C, OperationT Cop,
    const MatT* D, OperationT Dop)
{
    block_[0]= A; operation_[0]= Aop;
    block_[1]= B; operation_[1]= Bop;
    block_[2]= C; operation_[2]= Cop;
    block_[3]= D; operation_[3]= Dop;
}

template <class MatT>
bool
BlockMatrixBaseCL<MatT>::block_num_rows(size_t b, size_t& nr) const
{
    if (block_[b] == 0) return false;
    switch (operation_[b]) {
      case MUL:        nr= block_[b]->num_rows(); return true;
      case TRANSP_MUL: nr= block_[b]->num_cols(); return true;
      default:
        Comment("BlockMatrixBaseCL::block_num_rows: No such operation.\n", DebugNumericC);
        return false;
    }
}

template <class MatT>
bool
BlockMatrixBaseCL<MatT>::block_num_cols(size_t b, size_t& nc) const
{
    if (block_[b] == 0) return false;
    switch (operation_[b]) {
      case MUL:        nc= block_[b]->num_cols(); return true;
      case TRANSP_MUL: nc= block_[b]->num_rows(); return true;
      default:
        Comment("BlockMatrixBaseCL::block_num_cols: No such operation.\n", DebugNumericC);
        return false;
    }
}

template <class MatT>
size_t
BlockMatrixBaseCL<MatT>::num_rows(size_t block_row) const
{
    size_t ret= 0;
    bool block_found;
    switch (block_row) {
      case 0:
        block_found= block_num_rows( 0, ret) || block_num_rows( 1, ret);
        break;
      case 1:
        block_found= block_num_rows( 2, ret) || block_num_rows( 3, ret);
        break;
      default:
        Comment("BlockMatrixBaseCL::num_rows: No such block_row.\n", DebugNumericC);
        return 0;
    }
    Assert( block_found, "BlockMatrixBaseCL::num_rows: All pointers are 0.\n", DebugNumericC);
    return ret;
}

template <class MatT>
size_t
BlockMatrixBaseCL<MatT>::num_cols(size_t block_col) const
{
    size_t ret= 0;
    bool block_found;
    switch (block_col) {
      case 0:
        block_found= block_num_cols( 0, ret) || block_num_cols( 2, ret);
        break;
      case 1:
        block_found= block_num_cols( 1, ret) || block_num_cols( 3, ret);
        break;
      default:
        Comment("BlockMatrixBaseCL::num_cols: No such block_col.\n", DebugNumericC);
        return 0;
    }
    Assert( block_found, "BlockMatrixBaseCL::num_cols: All pointers are 0.\n", DebugNumericC);
    return ret;
}

template <class MatT>
BlockMatrixBaseCL<MatT>
BlockMatrixBaseCL<MatT>::GetTranspose() const
{
    return BlockMatrixBaseCL<MatT>( block_[0], GetTransposeOperation( 0),
        block_[2], GetTransposeOperation( 2),
        block_[1], GetTransposeOperation( 1),
        block_[3], GetTransposeOperation( 3));
}

template <typename MatT>
VectorBaseCL<typename MatT::value_type> BlockMatrixBaseCL<MatT>::GetDiag() const
{
    const size_t n=num_rows();
    Assert(n==num_cols(), "SparseMatBaseCL::GetDiag: no square Matrix", DebugParallelC);
    VectorBaseCL<typename MatT::value_type> diag(n);

    if (block_[0])
        diag[std::slice(0,num_rows(0),1)] = block_[0]->GetDiag();
    else
        diag[std::slice(0,num_rows(0),1)] = 1.;

    if (block_[3])
        diag[std::slice(num_rows(0),num_rows(1),1)] = block_[3]->GetDiag();
    else
        diag[std::slice(num_rows(0),num_rows(1),1)] = 1.;

    return diag;
}


//*****************************************************************************
//
//  Composition of 2 matrices
//
//*****************************************************************************
template <class MatT0, class MatT1>
class CompositeMatrixBaseCL
{
  public:
    typedef BlockMatrixOperationT OperationT;

  private:
    // The matrices are applied as block1_*block0_*v
    const MatT0* block0_;
    const MatT1* block1_;
    OperationT operation_[2];

  public:
    CompositeMatrixBaseCL( const MatT0* A, OperationT Aop, const MatT1* B, OperationT Bop);

    size_t num_rows() const;
    size_t num_cols() const;
    size_t intermediate_dim() const;

    const MatT0* GetBlock0 () const { return block0_; }
    const MatT1* GetBlock1 () const { return block1_; }
    void SetBlock0 (const MatT0* p) { block0_= p; }
    void SetBlock1 (const MatT1* p) { block1_= p; }

    OperationT GetOperation( size_t b) const { return operation_[b]; }
    OperationT GetTransposeOperation( size_t b) const {
        return operation_[b] == MUL ? TRANSP_MUL : MUL;
    }

    CompositeMatrixBaseCL<MatT1, MatT0> GetTranspose() const;
};

template <class MatT0, class MatT1>
CompositeMatrixBaseCL<MatT0, MatT1>::CompositeMatrixBaseCL( const MatT0* A, OperationT Aop,
    const MatT1* B, OperationT Bop)
{
    block0_= A; operation_[0]= Aop;
    block1_= B; operation_[1]= Bop;
}

template <class MatT0, class MatT1>
size_t
CompositeMatrixBaseCL<MatT0, MatT1>::num_rows() const
{
    switch (operation_[1]) {
      case MUL:        return block1_->num_rows();
      case TRANSP_MUL: return block1_->num_cols();
      default:
        Comment("CompositeMatrixBaseCL::num_rows: No such operation.\n", DebugNumericC);
        return (size_t)-1;
    }
}

template <class MatT0, class MatT1>
size_t
CompositeMatrixBaseCL<MatT0, MatT1>::num_cols() const
{
    switch (operation_[0]) {
      case MUL:        return block0_->num_cols();
      case TRANSP_MUL: return block0_->num_rows();
      default:
        Comment("CompositeMatrixBaseCL::num_cols: No such operation.\n", DebugNumericC);
        return (size_t)-1;
    }
}

template <class MatT0, class MatT1>
size_t
CompositeMatrixBaseCL<MatT0, MatT1>::intermediate_dim() const
{
    switch (operation_[0]) {
      case MUL:        return block0_->num_rows();
      case TRANSP_MUL: return block0_->num_cols();
      default:
        Comment("CompositeMatrixBaseCL::intermediate_dim: No such operation.\n", DebugNumericC);
        return (size_t)-1;
    }
}

template <class MatT0, class MatT1>
CompositeMatrixBaseCL<MatT1, MatT0>
CompositeMatrixBaseCL<MatT0, MatT1>::GetTranspose() const
{
    return CompositeMatrixBaseCL<MatT1,MatT0>( block1_, GetTransposeOperation( 1),
        block0_, GetTransposeOperation( 0));
}

//*****************************************************************************
//
//  Vector as diagonal matrix
//
//*****************************************************************************
template <class T>
class VectorAsDiagMatrixBaseCL
{
  private:
    const VectorBaseCL<T>* v_;

  public:
    VectorAsDiagMatrixBaseCL( const VectorBaseCL<T>* v)
        : v_( v) {}

    const VectorBaseCL<T>& GetDiag() const { return *v_; }

    size_t num_rows() const { return v_->size(); }
    size_t num_cols() const { return v_->size(); }
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
    A.resize( numrows, numcols, numnz);
    T* val= A.raw_val();
    size_t* row_beg= A.raw_row();
    size_t* col_ind= A.raw_col();

    for (size_t row=0; row<=numrows; ++row) is >> row_beg[row];
    for (size_t nz=0; nz<numnz; ++nz) is >> col_ind[nz];
    for (size_t nz=0; nz<numnz; ++nz) is >> val[nz];
}

// Read output of "operator<<".
template <typename T>
std::istream& operator>> (std::istream& in, SparseMatBaseCL<T>& A)
{
    size_t numrows, numcols, numnz;

    // Read the leading comment: % rows 'x' columns nonzeros "nonzeros\n"
    while (in && in.get() != '%') ;
    in >> numrows;
    while ( in && in.get() != 'x') ;
    in >> numcols >> numnz;
    while ( in && in.get() != '\n') ;
    if (!in)
        throw DROPSErrCL( "SparseMatBaseCL operator>>: Missing \"% rows cols nz\" comment.\n");

    SparseMatBuilderCL<T> B( &A, numrows, numcols);
    size_t r, c;
    T v;
    size_t nz;
    for (nz= 0; nz < numnz; ++nz) {
        in >> r >> c >> v;
        if (!in)
            throw DROPSErrCL( "SparseMatBaseCL operator>>: Stream corrupt.\n");
        if (r > numrows || c > numcols) {
            std::cerr << "r: " << r << "\tc: " << c << "\tv: " << v << "\tnz: " << nz << '\n';
            throw DROPSErrCL( "SparseMatBaseCL operator>>: Inconsistent data.\n");
        }
        B( r - 1, c - 1)= v; // Indices in the stream are 1-based.
    }
    B.Build();

    return in;
}

//*****************************************************************************
//
//  Matrix- and vector-operations
//
//*****************************************************************************

template <typename T>
  inline typename SparseMatBaseCL<T>::value_type
  supnorm(const SparseMatBaseCL<T>& M)
{
    typedef typename SparseMatBaseCL<T>::value_type valueT;
    const size_t nr= M.num_rows();
    valueT ret= valueT(), tmp;
    for (size_t i= 0; i < nr; ++i) {
        tmp= valueT();
        const size_t rowend= M.row_beg( i + 1);
        for (size_t j= M.row_beg( i); j < rowend; ++j)
            tmp+= std::fabs( M.val( j));
        if (tmp > ret) ret= tmp;
    }
    return ret;
}

template <typename T>
  typename SparseMatBaseCL<T>::value_type
  frobeniusnorm (const SparseMatBaseCL<T>& M)
{
    typedef typename SparseMatBaseCL<T>::value_type valueT;
    const size_t nz= M.num_nonzeros();
    valueT ret= valueT();
    for (size_t i= 0; i< nz; ++i)
            ret+= std::pow( M.val( i), 2);
    return std::sqrt( ret);
}

template <typename T>
  std::valarray<typename SparseMatBaseCL<T>::value_type>
  LumpInRows(const SparseMatBaseCL<T>& M)
{
    std::valarray<typename SparseMatBaseCL<T>::value_type>
        v( M.num_rows());
    for (size_t r= 0, nz= 0; nz < M.num_nonzeros(); ++r)
        for (; nz < M.row_beg( r + 1); )
            v[r]+= M.val( nz++);
    return v;
}

template <typename T>
  void
  ScaleRows(SparseMatBaseCL<T>& M,
    const std::valarray<typename SparseMatBaseCL<T>::value_type>& v)
{
    M.IncrementVersion();
    typename SparseMatBaseCL<T>::value_type* val= M.raw_val();
    for (size_t r= 0, nz= 0; nz < M.num_nonzeros(); ++r)
        for (; nz < M.row_beg( r + 1); )
            val[nz++]*= v[r];
}

template <typename T>
  void
  ScaleCols(SparseMatBaseCL<T>& M,
    const std::valarray<typename SparseMatBaseCL<T>::value_type>& v)
{
    M.IncrementVersion();
    typename SparseMatBaseCL<T>::value_type* val= M.raw_val();
    size_t* col= M.raw_col();
    for (size_t nz= 0; nz < M.num_nonzeros(); ++nz)
        val[nz]*= v[col[nz]];
}

/// \brief Compute the linear combination of two sparse matrices efficiently.
template <typename T>
SparseMatBaseCL<T>& SparseMatBaseCL<T>::LinComb (double coeffA, const SparseMatBaseCL<T>& A,
                                                 double coeffB, const SparseMatBaseCL<T>& B)
{
    Assert( A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(),
            "LinComb: incompatible dimensions", DebugNumericC);

    IncrementVersion();
    // Todo: Das alte Pattern wiederzuverwenden, macht mal wieder Aerger:
    // Zur Zeit (2.2008) mit der Matrix im NS-Loeser nach Gitteraenderungen, die
    // die Anzahl der Unbekannten nicht aendert. Daher schalten wir die
    // Wiederverwendung vorerst global aus.
    if (false && (!(DROPSDebugC & DebugNoReuseSparseC) && _val.size()!=0
        && _rows==A.num_rows() && _cols==A.num_cols()))
    {
        Comment("LinComb: Reusing OLD matrix" << std::endl, DebugNumericC);
        std::cerr << "LinComb: Reusing OLD matrix" << std::endl;

        size_t i=0, iA=0, iB=0;

        // same algorithm as below without writing _colind
        for (size_t row=1; row<=A.num_rows(); ++row)
        {
            while ( iA < A.row_beg(row) )
            {
                while ( iB < B.row_beg(row) && B.col_ind(iB) < A.col_ind(iA) )
                    _val[i++]=coeffB*B._val[iB++];
                if ( iB < B.row_beg(row) && B.col_ind(iB) == A.col_ind(iA) )
                    _val[i++]=coeffA*A._val[iA++]+coeffB*B._val[iB++];
                else
                    _val[i++]=coeffA*A._val[iA++];
            }
            while ( iB < B.row_beg(row) )
                _val[i++]=coeffB*B._val[iB++];
        }
        Assert( i==_val.size() && iA==A._val.size() && iB==B._val.size(), "LinComb: reuse of matrix pattern failed", DebugNumericC);
    }
    else
    {
        Comment("LinComb: Creating NEW matrix" << std::endl, DebugNumericC);
        std::cerr << "LinComb: Creating NEW matrix" << std::endl;

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
    }

    return *this;
}

/// \brief Compute the linear combination of three sparse matrices.
template <typename T>
SparseMatBaseCL<T>& SparseMatBaseCL<T>::LinComb (double coeffA, const SparseMatBaseCL<T>& A,
                                                 double coeffB, const SparseMatBaseCL<T>& B,
                                                 double coeffC, const SparseMatBaseCL<T>& C)
{
    SparseMatBaseCL<T> tmp;
    tmp.LinComb( coeffA, A, coeffB, B);
    return this->LinComb( 1.0, tmp, coeffC, C);
}

/// \brief Compute the transpose matrix of M explicitly.
/// This could be handled more memory-efficient by exploiting
/// that the column-indices in the rows of M are ascending.
template <typename T>
void
transpose (const SparseMatBaseCL<T>& M, SparseMatBaseCL<T>& Mt)
{
    SparseMatBuilderCL<T> Tr( &Mt, M.num_cols(), M.num_rows());
    for (size_t i= 0; i< M.num_rows(); ++i)
        for (size_t nz= M.row_beg( i); nz < M.row_beg( i + 1); ++nz)
            Tr( M.col_ind( nz), i)= M.val( nz);
    Tr.Build();
}   


/// \brief Compute the diagonal of B*B^T.
///
/// The commented out version computes B*M^(-1)*B^T
template <typename T>
VectorBaseCL<T>
BBTDiag (const SparseMatBaseCL<T>& B /*, const VectorBaseCL<T>& Mdiaginv*/)
{
    VectorBaseCL<T> ret( B.num_rows());
    
    T Bik;
    for (size_t i= 0; i < B.num_rows(); ++i) {
        for (size_t l= B.row_beg( i); l < B.row_beg( i + 1); ++l) {
            Bik= B.val( l);
            ret[i]+= /*Mdiaginv[B.col_ind( l)]**/ Bik*Bik;
        }
    }
    return ret;
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


template <typename _MatEntry, typename _VecEntry>
VectorBaseCL<_VecEntry>
operator*(const BlockMatrixBaseCL<SparseMatBaseCL<_MatEntry> >& A,
    const VectorBaseCL<_VecEntry>& x)
{
    VectorBaseCL<_VecEntry> x0( x[std::slice( 0, A.num_cols( 0), 1)]),
                            x1( x[std::slice( A.num_cols( 0), A.num_cols( 1), 1)]);
    VectorBaseCL<_VecEntry> r0( A.num_rows( 0)),
                            r1( A.num_rows( 1));
    const SparseMatBaseCL<_MatEntry>* mat;

    if ( (mat= A.GetBlock( 0)) != 0) {
        switch( A.GetOperation( 0)) {
          case MUL:
            r0= (*mat)*x0; break;
          case TRANSP_MUL:
            r0= transp_mul( *mat, x0); break;
        }
    }
    if ( (mat= A.GetBlock( 1)) != 0) {
        switch( A.GetOperation( 1)) {
          case MUL:
            r0+= (*mat)*x1; break;
          case TRANSP_MUL:
            r0+= transp_mul( *mat, x1); break;
        }
    }
    if ( (mat= A.GetBlock( 2)) != 0) {
        switch( A.GetOperation( 2)) {
          case MUL:
            r1= (*mat)*x0; break;
          case TRANSP_MUL:
            r1= transp_mul( *mat, x0); break;
        }
    }
    if ( (mat= A.GetBlock( 3)) != 0) {
        switch( A.GetOperation( 3)) {
          case MUL:
            r1+= (*mat)*x1; break;
          case TRANSP_MUL:
            r1+= transp_mul( *mat, x1); break;
        }
    }
    VectorBaseCL<_VecEntry> ret( A.num_rows());
    ret[std::slice( 0, A.num_rows( 0), 1)]= r0;
    ret[std::slice( A.num_rows( 0), A.num_rows( 1), 1)]= r1;
    return ret;
}

template <typename _MatEntry, typename _VecEntry>
VectorBaseCL<_VecEntry>
transp_mul(const BlockMatrixBaseCL<SparseMatBaseCL<_MatEntry> >& A,
    const VectorBaseCL<_VecEntry>& x)
{
    return A.GetTranspose()*x;
}


template <typename _MatT0, typename _MatT1, typename _VecEntry>
VectorBaseCL<_VecEntry>
operator*(const CompositeMatrixBaseCL<_MatT0, _MatT1>& A,
    const VectorBaseCL<_VecEntry>& x)
{
    VectorBaseCL<_VecEntry> tmp( A.intermediate_dim());
    switch( A.GetOperation( 0)) {
      case MUL:
        tmp= (*A.GetBlock0())*x; break;
      case TRANSP_MUL:
        tmp= transp_mul( *A.GetBlock0(), x); break;
    }
    VectorBaseCL<_VecEntry> ret( A.num_rows());
    switch( A.GetOperation( 1)) {
      case MUL:
        ret= (*A.GetBlock1())*tmp; break;
      case TRANSP_MUL:
        ret= transp_mul( *A.GetBlock1(), tmp); break;
    }
    return ret;

}

template <typename _MatT0, typename _MatT1, typename _VecEntry>
VectorBaseCL<_VecEntry>
transp_mul(const CompositeMatrixBaseCL<_MatT0, _MatT1>& A,
    const VectorBaseCL<_VecEntry>& x)
{
    return A.GetTranspose()*x;
}

template <typename _VecEntry>
VectorBaseCL<_VecEntry>
operator* (const VectorAsDiagMatrixBaseCL<_VecEntry>& A,
    const VectorBaseCL<_VecEntry>& x)
{
    return VectorBaseCL<_VecEntry>( A.GetDiag()*x);
}


template <typename _VecEntry>
VectorBaseCL<_VecEntry>
transp_mul (const VectorAsDiagMatrixBaseCL<_VecEntry>& A,
    const VectorBaseCL<_VecEntry>& x)
{
    return VectorBaseCL<_VecEntry>( A.GetDiag()*x);
}

//=============================================================================
//  Reverse Cuthill-McKee ordering 
//=============================================================================

/// \brief The first component is the degree of the vertex, the second is its number.
///
/// This is used by reverse_cuthill_mckee and its helpers.
typedef std::pair<size_t, size_t> DegVertT;
const size_t NoVert= std::numeric_limits<size_t>::max();

/// \brief Traverse the graph starting at v in breadth first preorder.
///
/// This is a helper of reverse_cuthill_mckee.  Vertices in each level
/// are numbered by increasing degree.
/// \param v The first vertex.
/// \param M Interpreted as adjacency matrix of the graph; see reverse_cuthill_mckee. Only verticed with p[v] != NoIdx are considered.
/// \param p Numbering of the vertices.
/// \param idx Next number to be used.
/// \param degree The degree of each vertex
template <typename T>
  void
  breadth_first (size_t v, const SparseMatBaseCL<T>& M, PermutationT& p, size_t& idx,
    std::vector<size_t>& degree)
{
    typedef std::deque<size_t> QueueT;
    QueueT Q;    

    // Insert v into queue and number v.
    Q.push_back( v);
    p[v]= idx++;

    std::vector<DegVertT> succ;
    size_t w;
    while (!Q.empty()) {
        v= Q.front();
        Q.pop_front();
        succ.resize( 0);
        succ.reserve( M.row_beg( v + 1) - M.row_beg( v));
        // Find all unnumbered successors of v.
        for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
            if (M.val( e) == 0.0 || (w= M.col_ind( e)) == v)
                continue;
            if (p[w] == NoVert)
                succ.push_back( std::make_pair( degree[w], w));
        }
        // Number successors by increasing degree and enter into queue.
        std::sort( succ.begin(), succ.end(), less1st<DegVertT>());
        for (size_t i= 0; i < succ.size(); ++i) {
            Q.push_back( succ[i].second);
            p[succ[i].second]= idx++;
        }
    }
}

/// \brief Find a vertex with high eccentricity.
///
/// This is a helper of reverse_cuthill_mckee. It computes (heuristically)
/// a start-vertex, the distance of which to at least one other vertex
/// is close to the diameter of of the graph. In a sense, such nodes
/// are on the boundary of the graph.
///
/// After a breadth first traversal a min-degree vertex from the last
/// level is chosen and the procedure is repeated if its eccentricity
/// is larger than that of the previous candidate.
///
/// \param V The vertices of the (sub-) graph together with their degree.
/// \param M Interpreted as adjacency matrix of the graph; see reverse_cuthill_mckee.
/// \param max_iter Maximal number of iteration. Though this heuristic will terminate
///     after at least |V| steps, this should be a small number like 5.
template <typename T>
  size_t
  pseudo_peripheral_node (std::vector<DegVertT>& V, const SparseMatBaseCL<T>& M,
    size_t& max_iter)
{
    PermutationT p_inv( M.num_rows(), NoVert); // Maps vertex v to its index in V.
    for (size_t i= 0; i < V.size(); ++i)
        p_inv[V[i].second]= i;

    // Select a minimal degree vertex.
    size_t v_old= arg_min( V.begin(), V.end(), less1st<DegVertT>())->second;

    int l_old= 0;
    size_t v, w;

    for (size_t iter= 0; iter < max_iter; ++iter) {
        // Breadth first traversal
        std::vector<int> level( V.size(), -1); // level[p_inv[v]]==-1 means: v is undiscovered.
        typedef std::deque<size_t> QueueT;
        QueueT Q, Qnext;

        // Insert v_old into queue; v_old is in level 0.
        Qnext.push_back( v_old);
        level[p_inv[v_old]]= 0;

        do {
            Q.swap( Qnext);
            Qnext.clear();
            for (QueueT::iterator it= Q.begin(); it != Q.end(); ++it) {
                v= *it;
                // Find all unnumbered successors of v.
                for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
                    if (M.val( e) == 0.0 || (w= M.col_ind( e)) == v)
                        continue;
                    if (level[p_inv[w]] == -1) {
                        Qnext.push_back( w);
                        level[p_inv[w]]= level[p_inv[v]] + 1;
                    }
                }
            }
        }
        while (!Qnext.empty());
        // Now Q contains the last level. It is not empty, if V is not empty.

        // Select a minimal degree vertex from the last level.
        std::vector<DegVertT> lastlevel;
        for (size_t i= 0; i < Q.size(); ++i) {
            lastlevel.push_back( V[p_inv[Q[i]]]);
        }
        v= arg_min( lastlevel.begin(), lastlevel.end(), less1st<DegVertT>())->second;

        if (level[p_inv[v]] <= l_old) {
            max_iter= iter;
            return v_old;
        }
        l_old= level[p_inv[v]];
        v_old= v;
    }
    return v_old;
}

/// \brief Compute a permutation such that M has (in most cases) near
/// minimal bandwith.
///
/// M_in is interpreted as adjacency matrix of a graph with vertices i
/// from [0...M.num_rows()).  Edge (i, j) exists, iff M_ij != 0.
///
/// From the linear algebra standpoint that means that x_i depends on
/// x_j, so these edges are actually the incoming edges of x_i. We use
/// this definition because we can traverse the incoming edges (rows)
/// much faster than the outgoing edges (columns). Thus, this routine
/// computes the rcm-permutation of M^T. Note, that the bandwith of A and
/// A^T is identical, so this probably does not matter. If desired the
/// function can compute the transpose of M and find the corresponding
/// rcm-permutation. This needs as much storage as an assembly of M.
///
/// \param M_in Interpreted as adjacency matrix of the graph.
/// \param p Contains for each unknown i its new number p[i].
/// \param use_indegree true: Compute the rcm-permutation of M^T
///     false: Compute the rcm-permutation of M; this is memory intensive.
template <typename T>
void
reverse_cuthill_mckee (const SparseMatBaseCL<T>& M_in, PermutationT& p,
    bool use_indegree= true)
{
    if (M_in.num_rows() == NoVert)
        throw DROPSErrCL( "reverse_cuthill_mckee: Graph is too large.\n");
    if (M_in.num_rows() != M_in.num_cols())
        throw DROPSErrCL( "reverse_cuthill_mckee: Matrix is not square.\n");

    size_t N= M_in.num_rows();
    SparseMatBaseCL<T> const* M;
    if (use_indegree == true)  M= &M_in;
    else {
        SparseMatBaseCL<T>* M2= new SparseMatBaseCL<T>();
        transpose( M_in, *M2);
        M= M2;
    }

    size_t idx= 0;
    p.assign( N, NoVert);

    std::vector<size_t> degree( N);    
    for (size_t r= 0; r < N; ++r)
        degree[r]= M->row_beg( r + 1) - M->row_beg( r);
    std::vector<DegVertT> V;
    V.reserve( N);
    for (size_t i= 0; i < N; ++i)
        V.push_back( std::make_pair( degree[i], i));

    do {
        // Find a start vertex v.
        size_t max_iter= 5;
        size_t v= pseudo_peripheral_node( V, *M, max_iter);
        std::cerr << "reverse_cuthill_mckee: p_p_n iterations: " << max_iter << '\n';

        // Traverse the graph starting with v and number the vertices.
        breadth_first( v, *M, p, idx, degree);

        V.clear();
        for (size_t i= 0; i < N; ++i) // Number all components of the graph.
            if (p[i] == NoVert)
                V.push_back( std::make_pair( degree[i], i));
    }
    while (!V.empty());

    if (use_indegree == false) delete M;

    // Revert the order of the vertices.
    for (size_t i= 0; i < N; ++i)
        p[i]= N - 1 - p[i];

    if (idx != N)
        throw DROPSErrCL( "reverse_cuthill_mckee: Could not number all unkowns.\n");
}

//=============================================================================
//  Typedefs
//=============================================================================

typedef VectorBaseCL<double>            VectorCL;
typedef SparseMatBaseCL<double>         MatrixCL;
typedef SparseMatBuilderCL<double>      MatrixBuilderCL;
typedef BlockMatrixBaseCL<MatrixCL>     BlockMatrixCL;
typedef CompositeMatrixBaseCL<MatrixCL, MatrixCL> CompositeMatrixCL;
typedef VectorAsDiagMatrixBaseCL<double>VectorAsDiagMatrixCL;

} // end of namespace DROPS

#endif

