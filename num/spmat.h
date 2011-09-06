/// \file spmat.h
/// \brief sparse matrix in compressed row format
/// \author LNM RWTH Aachen: Joerg Peters, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_SPMAT_H
#define DROPS_SPMAT_H

#define DROPS_SPARSE_MAT_BUILDER_USES_HASH_MAP (__GNUC__ >= 4 && !defined(__INTEL_COMPILER))

#include <iostream>
#include <valarray>
#include <vector>
#include <deque>
#include <numeric>
#include <limits>
#if DROPS_SPARSE_MAT_BUILDER_USES_HASH_MAP
#    include <tr1/unordered_map>
#else
#    include <map>
#endif
#include "misc/utils.h"
#include "misc/container.h"
#ifdef _PAR
# include "parallel/parallel.h"
#endif


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

    const T* raw() const { return Addr( *this); }
    T*       raw()       { return &(*this)[0]; }

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
    Assert(s<base_type::size(), "VectorBaseCL [] const: index out of bounds", DebugNumericC);
    return (*static_cast<const base_type*>( this))[s];
}

template <typename T>
T& VectorBaseCL<T>::operator[](size_t s)
{
    Assert(s<base_type::size(), "VectorBaseCL []: index out of bounds", DebugNumericC);
    return (*static_cast<base_type*>( this))[s];
}
#endif

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

/// \brief v[begin..(begin+2)]+= p[0..2]. A service function for the assembly of right-hand sides for vector-valued PDEs.
/// \return The updated vector
template <class T>
inline VectorBaseCL<T>& add_to_global_vector (VectorBaseCL<T>& v, const Point3DCL& p, size_t begin)
{
    for (int i= 0; i < 3; ++i)
        v[begin+i]+= p[i];
    return v;
}

/// \brief Use Kahan's algorithm to sum up elements 
/** This algorithm accumulates the error made by the floting point 
    arithmetics and addes this error to the sum.
    \param first iterator to the first element
    \param end   iterator behind the last element
    \param init  initiali value of the sum
    \return init + sum_{i=first}^end *i
*/
template <typename T, typename Iterator>
inline T KahanSumm( Iterator first, const Iterator& end, const T init=(T)0)
{
    T sum= init, c=T(0), t, y;
    while(first!=end){
        y  = *first++ - c;
        t  = sum + y;
        c  = (t-sum)-y;
        sum= t;
    }
    return sum;
}

/// \brief Use Kahan's algorithm to perform an inner product
template <typename T, typename Iterator>
inline T KahanInnerProd( Iterator first1, const Iterator& end1, Iterator first2, const T init=(T)0)
{
    T sum= init, c=T(0), t, y;
    while(first1!=end1){
        y  = (*first1++)*(*first2++) - c;
        t  = sum + y;
        c  = (t-sum)-y;
        sum= t;
    }
    return sum;
}

/// \brief Use Kahan's algorithm to perform an inner product on given indices
template <typename T, typename Cont, typename Iterator>
inline T KahanInnerProd( const Cont& a, const Cont&b, Iterator firstIdx, const Iterator& endIdx, const T init=(T)0)
{
    T sum= init, c=T(0), t, y;
    while(firstIdx!=endIdx){
        y  = a[*firstIdx]*b[*firstIdx] - c;
        ++firstIdx;
        t  = sum + y;
        c  = (t-sum)-y;
        sum= t;
    }
    return sum;
}

template<typename iterator>
  iterator*
  equal_split (iterator begin, iterator end, int num_threads, iterator* split_begin)
  {
    Assert( num_threads > 0, DROPSErrCL("equal_split: num_threads must be positive.\n"), DebugNumericC);
    const size_t chunk_size= (end - begin)/num_threads;
    const size_t remainder = (end - begin)%num_threads;

    iterator s= begin;
    *split_begin= begin;
    for (size_t i= 0; i < static_cast<size_t>( num_threads); ++i) {
        s+= i < remainder ? chunk_size + 1 : chunk_size;
        *++split_begin= s;
    }
    return split_begin;
  }

/// \brief Performs the std::partial_sum in parallel 
/// Can be used with any number of threads in a surrounding parallel region.
/// \param t_sum Array of size omp_num_threads(); used internally.
template<typename iterator>
void inplace_parallel_partial_sum (iterator begin, iterator end, iterator t_sum)
{
    if (begin == end)
        return;

    const Uint num_threads= omp_get_num_threads();
    const Uint tid= omp_get_thread_num();

    // Compute the part of [begin, end) to be considered in the current thread.
    const size_t chunk_size= (end - begin)/num_threads;
    const size_t remainder = (end - begin)%num_threads; // The threads [0,remainder) will have chunks of size chunk_size+1.
    const iterator t_begin= begin + tid*chunk_size + (tid < remainder ? tid : remainder);
    const iterator t_end= t_begin + chunk_size + (tid < remainder ? 1 : 0);

    // Compute the sum for each chunk in parallel
    typedef typename std::iterator_traits<iterator>::value_type value_type;
    value_type tmp= value_type();
    for (iterator p= t_begin; p < t_end; ++p)
        tmp+= *p;
    t_sum[tid]= tmp;
#   pragma omp barrier
    // Compute the prefix sums of the thread-local sums; these are the offsets
    // for the prefix sums on each thread.
    // This costs O(nthreads) ops. serially. There is a nice parallel algo to
    // compute the sum in O(lg(nthreads)) steps in parallel, but for nthreads=O(10) (and
    // end - begin = O(10000)), this does not matter.
#   pragma omp single
        std::partial_sum( t_sum, t_sum + num_threads, t_sum);

    // Add offset to the first summand for each thread and compute the final partial sums.
    if (tid > 0)
        *t_begin+= t_sum[tid - 1];
    std::partial_sum( t_begin, t_end, t_begin);
}


//*****************************************************************************
//
//  S p a r s e M a t B u i l d e r C L :  used for setting up sparse matrices
//
//*****************************************************************************
template <typename T>
class SparseMatBaseCL;

///\brief Traits for the SparseMatBuilderCL depending on the Blocks used in assembling the sparse matrix
///
/// SparseMatBaseCL<T> uses the compressed row storage format with entries of type T (e.g. T==double).
/// For building matrices for the Stokes-equations, it is useful to set-up the matrices with block-elements. For example, the deformation-tensor is representable as block-matrix with 3x3-blocks. By storing these blocks in the builder instead of the individual doubles, a factor 9 of calls to SparseMatBuilderCL::operator() (i,j) is saved. Memory savings result from the fact that the ratio of hash-value to numerical data improves by the same factor 9. For Setupsystem1_P2, the improved version is nearly twice as fast as the old version.
///
/// The use of std::map is about 10% slower than std::unordered_map + sorting entries.
///@{

/// \brief Generic traits for SparseMatbuilderCL. It must be specialized for individual block-types.
template <class BlockT>
struct BlockTraitsCL
{
    typedef BlockT block_type; ///< The data-blocks used in the builder (e.g. double, SMatrixCL<3,3>)
    typedef std::pair<size_t, block_type*> sort_pair_type; ///< Type used for sorting the rows when using hash-maps

    static const Uint num_rows= BlockT::num_rows; ///< Number of rows of one block
    static const Uint num_cols= BlockT::num_cols; ///< Number of columns of one block
    static const bool no_reuse= true; ///< False, if the sparsity-pattern can be reused with this block-type

    ///\brief Computes the beginning of double-valued rows in the matrix for one block-row
    static inline void row_begin ( size_t*, size_t); // not defined
    ///\brief Inserts one block-row in the form of double-valued rows into the matrix
    template <class Iter>
    static inline void insert_block_row (Iter, Iter, const size_t*, size_t*, double*); // not defined
    ///\brief Creates a pair for sorting the rows from a key-value pari stored in the hash-map.
    static inline sort_pair_type pair_copy (const std::pair<size_t, double>& p); // not defined
    ///\brief Return an entry, if the sparsity pattern is reused (only for block_type == double).
    static inline block_type& get_entry_reuse (const size_t*, const size_t*, size_t, double*); // not defined
};

template <>
struct BlockTraitsCL<double>
{
    typedef double block_type;
    typedef std::pair<size_t, double> sort_pair_type;

    static const Uint num_rows= 1;
    static const Uint num_cols= 1;
    static const bool no_reuse= false;

    static inline void row_begin (size_t* prev_row_end, size_t num_blocks)
        { prev_row_end[1]= /*prev_row_end[0] + */num_blocks; }
    template <class Iter>
    static inline void insert_block_row (Iter begin, Iter end, const size_t* rb, size_t* colind, double* val) {
        for (size_t j= rb[0]; begin != end; ++begin, ++j) {
            colind[j]= begin->first;
            val[j]= begin->second;
        }
    }
    static inline sort_pair_type pair_copy (const std::pair<size_t, double>& p)
        { return p; }

    static inline double& get_entry_reuse (const size_t* col_idx_begin, const size_t* col_idx_end, size_t j, double* val) {
        const size_t* pos= std::lower_bound( col_idx_begin, col_idx_end, j);
        Assert( pos != col_idx_end, "SparseMatBuilderCL (): no such index", DebugNumericC);
        return val[pos - col_idx_begin];
    }
};

template <Uint Rows, Uint Cols>
struct BlockTraitsCL< SMatrixCL<Rows, Cols> >
{
    typedef SMatrixCL<Rows, Cols> block_type;
    typedef std::pair<size_t, block_type*> sort_pair_type;

    static const Uint num_rows= Rows;
    static const Uint num_cols= Cols;
    static const bool no_reuse= true;

    static inline  void row_begin (size_t* prev_row_end, size_t num_blocks) {
        for (Uint k= 0; k < num_rows; ++k)
            prev_row_end[k + 1]= /*prev_row_end[k] +*/ num_cols*num_blocks;
    }
    template <class Iter>
    static inline void insert_block_row (Iter begin, Iter end, const size_t* rb, size_t* colind, double* val) {
        for (size_t l= 0 ; begin != end; ++begin, ++l)
            for (size_t i= 0; i < num_rows; ++i)
                for (size_t j= 0; j < num_cols; ++j) {
                    colind[j + l*num_cols + rb[i]]= begin->first*num_cols + j;
#if DROPS_SPARSE_MAT_BUILDER_USES_HASH_MAP
                    val   [j + l*num_cols + rb[i]]= (*begin->second)( i, j);
#else
                    val   [j + l*num_cols + rb[i]]= ( begin->second)( i, j);
#endif
                }
    }
    static inline sort_pair_type pair_copy (std::pair<const size_t, block_type>& p)
        { return std::make_pair( p.first, &p.second); }
    static inline block_type& get_entry_reuse (const size_t*, const size_t*, size_t, double*)
        { throw DROPSErrCL( "SparseMatBuilderCL: Cannot reuse the sparsity-pattern for SMatrixCL-blocks"); }
};

template <Uint Rows>
struct BlockTraitsCL< SDiagMatrixCL<Rows> >
{
    typedef SDiagMatrixCL<Rows> block_type;
    typedef std::pair<size_t, block_type*> sort_pair_type;

    static const Uint num_rows= Rows;
    static const Uint num_cols= Rows;
    static const bool no_reuse= true;

    static inline  void row_begin (size_t* prev_row_end, size_t num_blocks) {
        for (Uint k= 0; k < num_rows; ++k)
            prev_row_end[k + 1]= /*prev_row_end[k] +*/ num_blocks;
    }
    template <class Iter>
    static inline void insert_block_row (Iter begin, Iter end, const size_t* rb, size_t* colind, double* val) {
        for (size_t l= 0 ; begin != end; ++begin, ++l)
            for (size_t i= 0; i < num_rows; ++i) {
                colind[l + rb[i]]= begin->first*num_cols + i;
#if DROPS_SPARSE_MAT_BUILDER_USES_HASH_MAP
                val   [l + rb[i]]= (*begin->second)( i);
#else
                val   [l + rb[i]]= ( begin->second)( i);
#endif
            }
    }
    static inline sort_pair_type pair_copy (std::pair<const size_t, block_type>& p)
        { return std::make_pair( p.first, &p.second); }
    static inline block_type& get_entry_reuse (const size_t*, const size_t*, size_t, double*)
        { throw DROPSErrCL( "SparseMatBuilderCL: Cannot reuse the sparsity-pattern for SDiagMatrixCL-blocks"); }
};
///@}

/// \brief Building sparse matrices
///
/// \param T is the type of the matrix-entries
/// \param BlockT is a T-valued container-type used in the builder
template <typename T= double, typename BlockT= T>
class SparseMatBuilderCL
{
private:
    typedef BlockTraitsCL<BlockT> BlockTraitT;
    typedef typename BlockTraitT::block_type block_type;

public:
    typedef T                        valueT;
    typedef SparseMatBaseCL<T>       spmatT;
    typedef std::pair<size_t, block_type> entryT;
#if DROPS_SPARSE_MAT_BUILDER_USES_HASH_MAP
    typedef std::tr1::unordered_map<size_t, block_type> couplT;
#else
    typedef std::map<size_t, block_type> couplT;
#endif

private:
    size_t  _rows;
    size_t  _cols;
    spmatT* _mat;
    bool    _reuse;
    couplT* _coupl;

public:
    SparseMatBuilderCL(spmatT* mat, size_t rows, size_t cols, bool reuse= false)
        : _rows(rows), _cols(cols), _mat(mat), _reuse( reuse)
    {
        Assert( _rows%BlockTraitT::num_rows == 0, "SparseMatBuilderCL (): number of rows incompatible with block_type", DebugNumericC);
        Assert( _cols%BlockTraitT::num_cols == 0, "SparseMatBuilderCL (): number of columns incompatible with block_type", DebugNumericC);
        mat->IncrementVersion();
        if (_reuse)
        {
            if (BlockTraitT::no_reuse)
                throw DROPSErrCL( "SparseMatBuilderCL: Cannot reuse the pattern for block_type != double.");
            Comment("SparseMatBuilderCL: Reusing OLD matrix" << std::endl, DebugNumericC);
            _coupl=0;
            std::memset( _mat->_val, 0, mat->num_nonzeros()*sizeof( T));
        }
        else
        {
            Comment("SparseMatBuilderCL: Creating NEW matrix" << std::endl, DebugNumericC);
            _coupl= new couplT[_rows/BlockTraitT::num_rows];
            for (size_t i=0; i< _rows/BlockTraitT::num_rows; ++i) 
                _coupl[i].rehash(100);
        }
    }

    ~SparseMatBuilderCL() { if (_coupl) delete[] _coupl; }

    block_type& operator() (size_t i, size_t j)
    {
        Assert( i < _rows/BlockTraitT::num_rows && j <_cols/BlockTraitT::num_cols,
            "SparseMatBuilderCL (): index out of bounds", DebugNumericC);

        if (!_reuse)
            return _coupl[i/BlockTraitT::num_rows][j/BlockTraitT::num_cols];
        else
            return BlockTraitT::get_entry_reuse( _mat->GetFirstCol( i),  _mat->GetFirstCol( i + 1), j, _mat->GetFirstVal( i));
    }

    void Build();
};

template <typename T, typename BlockT>
void SparseMatBuilderCL<T, BlockT>::Build()
{
    if (_reuse) return;

    const size_t block_rows= _rows/BlockTraitT::num_rows;
    _mat->num_rows( _rows);
    _mat->num_cols( _cols);

    size_t* rb= _mat->raw_row();
    rb[0]= 0;

    size_t* t_sum= new size_t[omp_get_max_threads()];
#pragma omp parallel
{
#pragma omp for
    for (size_t i= 0; i < block_rows; ++i)
        BlockTraitT::row_begin( rb + i*BlockTraitT::num_rows, _coupl[i].size());

    inplace_parallel_partial_sum( rb, rb + _mat->num_rows() + 1, t_sum); 
#pragma omp barrier
#pragma omp master
    _mat->num_nonzeros( rb[_rows]);
#pragma omp barrier

#if DROPS_SPARSE_MAT_BUILDER_USES_HASH_MAP
    std::vector<typename BlockTraitT::sort_pair_type> pv;
#pragma omp for
    for (size_t i= 0; i < block_rows; ++i) {
        pv.resize( _coupl[i].size());
        std::transform( _coupl[i].begin(), _coupl[i].end(), pv.begin(), &BlockTraitT::pair_copy);
        std::sort( pv.begin(), pv.end(), less1st<typename BlockTraitT::sort_pair_type>());
        BlockTraitT::insert_block_row( pv.begin(), pv.end(), rb + i*BlockTraitT::num_rows, _mat->raw_col(), _mat->raw_val());
        // std::cout << _coupl[i].load_factor() << '\t' << std::setfill('0') << std::setw(3) <<_coupl[i].size() << '\n';
    }
}
#else
    for (size_t i= 0; i < block_rows; ++i)
        BlockTraitT::insert_block_row( _coupl[i].begin(), _coupl[i].end(), rb + i*BlockTraitT::num_rows, _mat->raw_col(), _mat->raw_val());
#endif
    delete[] t_sum;
    delete[] _coupl;
    _coupl= 0;
}

///\brief  SparseMatBaseCL: compressed row storage sparse matrix
/// Use SparseMatBuilderCL for setting up.
template <typename T>
class SparseMatBaseCL
{
private:
    size_t _rows; ///< number of rows
    size_t _cols; ///< number of columns
    size_t nnz_;  ///< number of non-zeros

    size_t version_; ///< All modifications increment this. Starts with 1.

    size_t* _rowbeg; ///< (_rows+1 entries, last entry must be <=_nz) index of first non-zero-entry in _val belonging to the row given as subscript
    size_t* _colind; ///< (nnz_ entries) column-number of corresponding entry in _val
    T*      _val;    ///< (nnz_ entries) the components of the matrix

    void num_rows (size_t rows);
    void num_cols (size_t cols);
    void num_nonzeros (size_t nnz);

public:
    typedef T value_type;

    SparseMatBaseCL (); ///< empty zero-matrix
    SparseMatBaseCL (const SparseMatBaseCL&);
    ~SparseMatBaseCL ();

    SparseMatBaseCL (size_t rows, size_t cols, size_t nnz); ///< the fields are allocated, but not initialized
    SparseMatBaseCL (size_t rows, size_t cols, size_t nnz,  ///< construct matrix from the CRS-fields
                     const T* valbeg , const size_t* rowbeg, const size_t* colindbeg);
    SparseMatBaseCL (const std::valarray<T>&); ///< Creates a square diagonal matrix.

    SparseMatBaseCL& operator= (const SparseMatBaseCL& m);

    const T*      raw_val() const { return _val; }
    T*            raw_val()       { return _val; }
    const size_t* raw_row() const { return _rowbeg; }
    size_t*       raw_row()       { return _rowbeg; }
    const size_t* raw_col() const { return _colind; }
    size_t*       raw_col()       { return _colind; }

    size_t num_rows     () const { return _rows; }
    size_t num_cols     () const { return _cols; }
    size_t num_nonzeros () const { return nnz_; }
#ifdef _PAR
    size_t num_acc_nonzeros() const { return ProcCL::GlobalSum(num_nonzeros()); }
#endif

    size_t row_beg (size_t i) const { return _rowbeg[i]; }
    size_t col_ind (size_t i) const { return _colind[i]; }
    T      val     (size_t i) const { return _val[i]; }

    void IncrementVersion() { ++version_; }      ///< Increment modification version number
    size_t Version() const  { return version_; } ///< Get modification version number

    const size_t* GetFirstCol(size_t i) const { return _colind + _rowbeg[i]; }
          size_t* GetFirstCol(size_t i)       { return _colind + _rowbeg[i]; }
    const T*      GetFirstVal(size_t i) const { return _val    + _rowbeg[i]; }
          T*      GetFirstVal(size_t i)       { return _val    + _rowbeg[i]; }

    inline T  operator() (size_t i, size_t j) const;

    SparseMatBaseCL& operator*= (T c);
    SparseMatBaseCL& operator/= (T c);

    SparseMatBaseCL& LinComb (double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&);
    SparseMatBaseCL& LinComb (double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&);
    SparseMatBaseCL& LinComb (double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&,
                              double, const SparseMatBaseCL<T>&);

    void insert_col (size_t c, const VectorBaseCL<T>& v);

    ///\brief Resize to new dimensions. The old content is lost.
    void resize (size_t rows, size_t cols, size_t nnz);
    void clear  () { resize(0,0,0); }

    VectorBaseCL<T> GetDiag()                              const;
    VectorBaseCL<T> GetLumpedDiag()                        const;
    VectorBaseCL<T> GetSchurDiag(const VectorBaseCL<T>& W) const; ///< returns diagonal of B W B^T

    void permute_rows (const PermutationT&);
    void permute_columns (const PermutationT&);

    template <class, class>
      friend class SparseMatBuilderCL;
};

template <typename T>
  void
  SparseMatBaseCL<T>::num_rows (size_t rows)
{
    _rows= rows;
    delete[] _rowbeg;
    _rowbeg= new size_t[rows + 1];
}

template <typename T>
  void
  SparseMatBaseCL<T>::num_cols (size_t cols)
{
    _cols= cols;
}

template <typename T>
  void
  SparseMatBaseCL<T>::num_nonzeros (size_t nnz)
{
    nnz_= nnz;
    delete[] _val;
    delete[] _colind;
    _val= new T[nnz];
    _colind= new size_t[nnz];
}

template <typename T>
  SparseMatBaseCL<T>::SparseMatBaseCL ()
    : _rows(0), _cols(0), nnz_( 0), version_(1), _rowbeg( new size_t[1]), _colind(0), _val(0)
{
    _rowbeg[0]= 0;
}

template <typename T>
  SparseMatBaseCL<T>::SparseMatBaseCL (const SparseMatBaseCL& m)
    : _rows( m._rows), _cols( m._cols), nnz_( m.nnz_), version_( m.version_),
      _rowbeg( new size_t[m._rows+1]), _colind( new size_t[m.num_nonzeros()]), _val(new T[m.num_nonzeros()])
{
    std::copy( m.raw_row(), m.raw_row() + m.num_rows() + 1, raw_row());
    std::copy( m.raw_col(), m.raw_col() + m.num_nonzeros(), raw_col());
    std::copy( m.raw_val(), m.raw_val() + m.num_nonzeros(), raw_val());
}

template <typename T>
  SparseMatBaseCL<T>::~SparseMatBaseCL ()
{
    delete[] _rowbeg;
    delete[] _colind;
    delete[] _val;
}

template <typename T>
  SparseMatBaseCL<T>::SparseMatBaseCL (size_t rows, size_t cols, size_t nnz)
    : _rows( rows), _cols( cols), nnz_( nnz), version_( 1),
      _rowbeg( new size_t[rows+1]), _colind( new size_t[nnz]), _val( new T[nnz])
{
    // std::memset( _rowbeg, 0, (_rows + 1)*sizeof( size_t));
    // std::memset( _colind, 0, nnz_*sizeof( size_t));
    // std::memset( _val,    0, nnz_*sizeof( T));
}

template <typename T>
  SparseMatBaseCL<T>::SparseMatBaseCL (size_t rows, size_t cols, size_t nnz,
    const T* valbeg , const size_t* rowbeg, const size_t* colindbeg)
    : _rows(rows), _cols(cols), nnz_(nnz), version_( 1),
      _rowbeg( new size_t[rows+1]), _colind( new size_t[nnz]), _val( new T[nnz])
{
    std::copy( rowbeg, rowbeg + num_rows() + 1, raw_row());
    std::copy( colindbeg, colindbeg + num_nonzeros(), raw_col());
    std::copy( valbeg,    valbeg    + num_nonzeros(), raw_val());
}

template <typename T>
  SparseMatBaseCL<T>::SparseMatBaseCL(const std::valarray<T>& v)
      : _rows( v.size()), _cols( v.size()), nnz_( v.size()), version_( 1),
        _rowbeg( new size_t[v.size() + 1]), _colind( new size_t[v.size()]), _val( new T[v.size()])
{
    for (size_t i= 0; i < _rows; ++i)
        _rowbeg[i]= _colind[i]= i;
    _rowbeg[_rows]= _rows;
    std::copy( Addr( v), Addr( v) + num_nonzeros(), raw_val());
}

template <typename T>
  SparseMatBaseCL<T>& SparseMatBaseCL<T>::operator= (const SparseMatBaseCL<T>& m)
{
    if (&m == this) return *this;

    resize( m.num_rows(), m.num_cols(), m.num_nonzeros());
    std::copy( m.raw_row(), m.raw_row() + m.num_rows() + 1, raw_row());
    std::copy( m.raw_col(), m.raw_col() + m.num_nonzeros(), raw_col());
    std::copy( m.raw_val(), m.raw_val() + m.num_nonzeros(), raw_val());
    return *this;
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
  SparseMatBaseCL<T>&
  SparseMatBaseCL<T>::operator*= (T c)
{
    IncrementVersion();
    for (size_t i= 0; i < nnz_; ++i)
        _val[i]*= c;
    return *this;
}

template <typename T>
  SparseMatBaseCL<T>&
  SparseMatBaseCL<T>::operator/= (T c)
{
    IncrementVersion();
    for (size_t i= 0; i < nnz_; ++i)
        _val[i]/= c;
    return *this;
}

template <typename T>
  void
  SparseMatBaseCL<T>::resize (size_t rows, size_t cols, size_t nnz)
{
    IncrementVersion();
    num_rows( rows);
    num_cols( cols);
    num_nonzeros( nnz);
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
VectorBaseCL<T> SparseMatBaseCL<T>::GetLumpedDiag() const
{
    const size_t n=num_rows(),
        nnz=num_nonzeros();
    Assert(n==num_cols(), "SparseMatBaseCL::GetLumpedDiag: no square Matrix", DebugParallelC);
    VectorBaseCL<T> diag(n);
    for (size_t r= 0, nz= 0; nz < nnz; ++r)
        for (; nz < _rowbeg[ r + 1]; ++nz) {
            const double v= _val[nz];
            diag[r]+= std::abs(v);
        }
    return diag;
}

template <typename T>
VectorBaseCL<T> SparseMatBaseCL<T>::GetSchurDiag( const VectorBaseCL<T>& W) const
/// In parallel, this function may not work as expected
{
    const size_t n=num_rows(),
        nnz=num_nonzeros();
    VectorBaseCL<T> diag(n);
    for (size_t r= 0, nz= 0; nz < nnz; ++r)
        for (; nz < _rowbeg[ r + 1]; ++nz) {
            const double v= _val[nz];
            diag[r]+= v*v*W[_colind[nz]];
        }
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
    SparseMatBaseCL<T> tmp( *this);

    _rowbeg[0]= 0;
    for (size_t i= 0; i < num_rows(); ++i) {
        _rowbeg[i + 1]= _rowbeg[i] + tmp.row_beg( pi[i] + 1) - tmp.row_beg( pi[i]);
        std::copy( tmp.GetFirstVal( pi[i]),  tmp.GetFirstVal( pi[i] + 1), GetFirstVal( i));
        std::copy( tmp.GetFirstCol( pi[i]),  tmp.GetFirstCol( pi[i] + 1), GetFirstCol( i));
    }
}

template <typename T>
  void
  SparseMatBaseCL<T>::permute_columns (const PermutationT& p)
{
    Assert( num_cols() == p.size(),
        DROPSErrCL( "permute_columns: Matrix and Permutation have different dimension.\n"), DebugNumericC);

    IncrementVersion();
    for (size_t i= 0; i < nnz_; ++i)
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
            std::cout << "r: " << r << "\tc: " << c << "\tv: " << v << "\tnz: " << nz << '\n';
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
/// In parallel, this function may not work as expected
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
/// In parallel, this function may not work as expected
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
/// In parallel, this function returns the distributed form of the lumped diagonal (can be accumulated afterwards)
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

template <typename T>
void ortho( VectorBaseCL<T>& v, const VectorBaseCL<T>& k) /// orthogonalize v w.r.t. k
{
    const double alpha= dot(v,k)/dot(k,k);
    v-= alpha*k;
}


template <typename T>
void ortho( VectorBaseCL<T>& v, const VectorBaseCL<T>& k, const SparseMatBaseCL<T>& Y) /// orthogonalize v w.r.t. k and inner product induced by Y
{
    const VectorBaseCL<T> Yk(Y*k);
    const double alpha= dot(v,Yk)/dot(k,Yk);
    v-= alpha*k;
}

/// \brief Compute the linear combination of two sparse matrices efficiently.
/// The new sparsity pattern is computed by merging the lists of _colind (removing the duplicates) row by row.
/// \todo Das alte Pattern wiederzuverwenden, macht mal wieder Aerger:
///   Zur Zeit (2.2008) mit der Matrix im NS-Loeser nach Gitteraenderungen, die
///   die Anzahl der Unbekannten nicht aendert. Daher schalten wir die
///   Wiederverwendung vorerst global aus.
template <typename T>
SparseMatBaseCL<T>& SparseMatBaseCL<T>::LinComb (double coeffA, const SparseMatBaseCL<T>& A,
                                                 double coeffB, const SparseMatBaseCL<T>& B)
{
    Assert( A.num_rows()==B.num_rows() && A.num_cols()==B.num_cols(),
            "LinComb: incompatible dimensions", DebugNumericC);

    IncrementVersion();
    Comment( "LinComb: Creating NEW matrix" << std::endl, DebugNumericC);
    num_rows( A.num_rows());
    num_cols( A.num_cols());
    _rowbeg[0]= 0;
    size_t* t_sum= new size_t[omp_get_max_threads()];

#pragma omp parallel
    {
    // Compute the entries of _rowbeg (that is the number of nonzeros in each row of the result)
    size_t i;
    const size_t* rA;
    const size_t* rB;
#pragma omp for
    for (size_t row= 0; row < A.num_rows(); ++row) {
        i= 0;//row_beg( row);
        rA= A.GetFirstCol( row);
        rB= B.GetFirstCol( row);
        const size_t* const rAend= A.GetFirstCol( row + 1);
        const size_t* const rBend= B.GetFirstCol( row + 1);
        for (; rA != rAend && rB != rBend; ++i)
            if (*rB < *rA)
                ++rB;
            else if (*rA < *rB)
                ++rA;
            else {
                ++rA;
                ++rB;
            }
        _rowbeg[row + 1]= i + (rAend - rA) + (rBend - rB);
    }

    inplace_parallel_partial_sum( _rowbeg, _rowbeg + num_rows() + 1, t_sum); 
#   pragma omp barrier
#    pragma omp master
        num_nonzeros( row_beg( num_rows()));
#    pragma omp barrier
    // Compute the entries of _colind, _val (actual merge).
    size_t iA, iB;

#pragma omp for
    for (size_t row= 0; row < A.num_rows(); ++row) { // same structure as above
        i=    row_beg( row);
        iA= A.row_beg( row);
        iB= B.row_beg( row);
        const size_t rAend= A.row_beg( row + 1),
                     rBend= B.row_beg( row + 1);
        for (; iA != rAend && iB != rBend; ++i)
            if (B.col_ind( iB) < A.col_ind( iA)) {
                _val[i]= coeffB*B._val[iB];
                _colind[i]= B._colind[iB++];
            }
            else if (A.col_ind( iA) < B.col_ind( iB)) {
                _val[i]= coeffA*A._val[iA];
                _colind[i]= A._colind[iA++];
            }
            else {
                _val[i]= coeffA*A._val[iA++] + coeffB*B._val[iB];
                _colind[i]= B._colind[iB++];
            }
        // At most one of A or B might have entries left.
        std::copy( B._colind + iB, B._colind + rBend,
                   std::copy( A._colind + iA, A._colind + rAend, _colind + i));
        for (; iA < rAend; ++iA, ++i)
            _val[i]= coeffA*A._val[iA];
        for (; iB < rBend; ++iB, ++i)
            _val[i]= coeffB*B._val[iB];
    }
    } //end of omp parallel
    delete [] t_sum;
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

/// \brief Compute the linear combination of four sparse matrices.
template <typename T>
SparseMatBaseCL<T>& SparseMatBaseCL<T>::LinComb (double coeffA, const SparseMatBaseCL<T>& A,
                                                 double coeffB, const SparseMatBaseCL<T>& B,
                                                 double coeffC, const SparseMatBaseCL<T>& C,
                                                 double coeffD, const SparseMatBaseCL<T>& D)
{
    SparseMatBaseCL<T> tmp;
    tmp.LinComb( coeffA, A, coeffB, B, coeffC, C);
    return this->LinComb( 1.0, tmp, coeffD, D);
}

/// \brief Inserts v as column c. The old columns [c, num_cols()) are shifted to the right.
///
/// If c > num_cols(), implicit zero-columns will show up in the matrix.
template <typename T>
inline void
SparseMatBaseCL<T>::insert_col (size_t c, const VectorBaseCL<T>& v)
{
    if (c > num_cols())
        throw DROPSErrCL( "SparseMatBaseCL<T>::insert_col: Only one column can be inserted.\n");

    const size_t numzero= std::count( Addr( v), Addr( v) + v.size(), 0.),
                 nnz= v.size() - numzero;
    size_t       zerocount= 0,
                 shift= 1;
    IncrementVersion();
    // These will become the new data-arrays of the matrix.
    std::valarray<size_t> rowbeg( num_rows() + 1);
    std::valarray<size_t> colind( num_nonzeros() + nnz);
    std::valarray<T>      val(    num_nonzeros() + nnz);
    size_t* newcbeg= Addr( colind);
    T*      newvbeg= Addr( val);

    for (size_t row= 0; row < num_rows(); ++row) {
        // Column indices: Copy the entries up to column c.
        rowbeg[row]= newcbeg - Addr( colind);
        const size_t* cbeg= GetFirstCol( row);
        const size_t* cend= GetFirstCol( row + 1);
        const size_t *cpos= std::lower_bound( cbeg, cend, c);
        newcbeg= std::copy( cbeg, cpos, newcbeg);
        // Same for the val-array.
        const T* vbeg= GetFirstVal( row);
        const T* vend= GetFirstVal( row + 1);
        const T* vpos= vbeg + (cpos - cbeg);
        newvbeg= std::copy( vbeg, vpos, newvbeg);

        // Insert c and the vector-entry.
        if ( v[row] != 0.) {
            *newcbeg++= c;
            *newvbeg++= v[row];
            shift= 1;
        }
        else {
            ++zerocount;
            shift= 0;
        }

        // Copy the rest of the row and shift column indices by shift.
        for( const size_t* p= cpos; p != cend; ++p, ++newcbeg)
            *newcbeg= *p + shift;
        // The same for the val-array; no index shifting here.
        newvbeg= std::copy( vpos, vend, newvbeg);
    }
    if (zerocount != numzero)
        throw DROPSErrCL( "SparseMatBaseCL<T>::insert_col: Inconsistent zero-counts in v.\n");

    num_cols( num_cols() + c <= _cols ? 1 : c - _cols + 1);
    // Adjust the last rowbeg-entry.
    rowbeg[num_rows()]= colind.size();
    // Copy adapted arrays into the matrix.
    std::copy( &rowbeg[0], &rowbeg[0] + num_rows() + 1, _rowbeg);
    num_nonzeros( colind.size());
    std::copy( &colind[0], &colind[0] + colind.size(), _colind);
    std::copy( &val[0],    &val[0]    + val.size(),    _val);
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
///
/// In parallel, this function may not work as expected
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
    size_t rowend, rowbeg;
    size_t nz= 0;

#pragma omp parallel for private(sum, rowend, rowbeg, nz)
    for (size_t i = 0; i < num_rows; i++){
        sum = 0.0;
        rowend = Arow[i+1];
        rowbeg = Arow[i];
        for (nz=rowbeg; nz<rowend; ++nz)
            sum += Aval[nz] * x[Acol[nz]];
        y[i] = sum;
    }
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


/// \brief returns the entry in row 'row' of A*x.
template <typename Mat, typename Vec>
inline typename Vec::value_type
mul_row (const Mat& A, const Vec& x, size_t row)
{
    typename Vec::value_type sum= typename Vec::value_type();
    for (size_t k= A.row_beg( row); k < A.row_beg( row + 1 ); ++k)
        sum+= A.val( k)*x[A.col_ind( k)];
    return sum;
}


/// \brief x+= w*A(row, .); the updated vector is returned.
template <typename Mat, typename Vec>
inline Vec&
add_row_to_vec (const Mat& A, double w, Vec& x, size_t row)
{
    for (size_t k= A.row_beg( row); k < A.row_beg( row + 1 ); ++k)
        x[A.col_ind( k)]+= w*A.val( k);
    return x;
}

/// \brief x+= w*A(., col); the updated vector is returned.
template <typename Mat, typename Vec>
inline Vec&
add_col_to_vec (const Mat& A, double w, Vec& x, size_t col)
{
    for (size_t row= 0; row < A.num_rows(); ++row) {
        x[row]+= w*A( row, col);
    }
    return x;
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
        std::cout << "reverse_cuthill_mckee: p_p_n iterations: " << max_iter << '\n';

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

//*****************************************************************************
//
//  MLSparseMatBaseCL
//
//*****************************************************************************

template <typename T>
class MLSparseMatBaseCL : public MLDataCL<SparseMatBaseCL<T> >
{
  private:
    typedef typename MLSparseMatBaseCL<T>::const_iterator ML_const_iterator;
    typedef typename MLSparseMatBaseCL<T>::iterator       ML_iterator;
  public:
    MLSparseMatBaseCL (size_t lvl= 1)
    {
#ifdef _PAR
        if (lvl>1)
            throw DROPSErrCL("MLSparseMatBaseCL::MLSparseMatBaseCL: No multilevel matrices in the parallel version, yet, sorry");
#endif
        this->resize(lvl);
    }
    size_t Version      () const { return this->GetFinest().Version();}
    size_t num_nonzeros () const { return this->GetFinest().num_nonzeros(); }
    size_t num_rows     () const { return this->GetFinest().num_rows(); }
    size_t num_cols     () const { return this->GetFinest().num_cols(); }
    MLSparseMatBaseCL<T>& LinComb (double ma, const MLSparseMatBaseCL<T>& A,
                                   double mb, const MLSparseMatBaseCL<T>& B);
    MLSparseMatBaseCL<T>& LinComb (double ma, const MLSparseMatBaseCL<T>& A,
                                   double mb, const MLSparseMatBaseCL<T>& B,
                                   double mc, const MLSparseMatBaseCL<T>& C);

    VectorBaseCL<T> GetDiag() const { return this->GetFinest().GetDiag(); }
    void clear() { for (ML_iterator it = this->begin(); it != this->end(); ++it) it->clear();}
};

template <typename T>
MLSparseMatBaseCL<T>& MLSparseMatBaseCL<T>::LinComb (double ma, const MLSparseMatBaseCL<T>& A, double mb, const MLSparseMatBaseCL<T>& B)
{
    Assert( A.size()==B.size(), "MLMatrixCL::LinComb: different number of levels", DebugNumericC);
    ML_const_iterator itA = A.begin();
    ML_const_iterator itB = B.begin();
    SparseMatBaseCL<T> mat;
    this->resize( A.size());
    for (ML_iterator it = this->begin(); it != this->end(); ++it)
    {
        mat.LinComb( ma, *itA, mb, *itB);
        *it = mat;
        ++itA;
        ++itB;
    }
    return *this;
}

template <typename T>
MLSparseMatBaseCL<T>& MLSparseMatBaseCL<T>::LinComb (double ma, const MLSparseMatBaseCL<T>& A, double mb, const MLSparseMatBaseCL<T>& B, double mc, const MLSparseMatBaseCL<T>& C)
{
    Assert( A.size()==B.size(), "MLMatrixCL::LinComb: different number of levels", DebugNumericC);
    Assert( A.size()==C.size(), "MLMatrixCL::LinComb: different number of levels", DebugNumericC);
    ML_const_iterator itA = A.begin();
    ML_const_iterator itB = B.begin();
    ML_const_iterator itC = C.begin();
    SparseMatBaseCL<T> mat;
    this->resize( A.size());
    for (ML_iterator it = this->begin(); it != this->end(); ++it)
    {
        mat.LinComb( ma, *itA, mb, *itB, mc, *itC);
        *it = mat;
        ++itA;
        ++itB;
        ++itC;
    }
    return *this;
}

template <typename _MatEntry, typename _VecEntry>
VectorBaseCL<_VecEntry> transp_mul (const MLSparseMatBaseCL<_MatEntry>& A, const VectorBaseCL<_VecEntry>& x)
{
    return transp_mul(A.GetFinest(), x);
}

template <typename _MatEntry, typename _VecEntry>
VectorBaseCL<_VecEntry> operator* (const MLSparseMatBaseCL<_MatEntry>& A, const VectorBaseCL<_VecEntry>& x)
{
    return A.GetFinest()*x;
}

//Human Readable
template <typename T>
std::ostream& operator << (std::ostream& os, const MLSparseMatBaseCL<T>& A)
{
	for (typename MLSparseMatBaseCL<T>::const_iterator it=A.begin(); it!=A.end(); ++it)
		os << *it << std::endl;
    return os;
}

//=============================================================================
//  Typedefs
//=============================================================================

typedef VectorBaseCL<double>             VectorCL;
typedef SparseMatBaseCL<double>          MatrixCL;
typedef SparseMatBuilderCL<>             MatrixBuilderCL;
typedef VectorAsDiagMatrixBaseCL<double> VectorAsDiagMatrixCL;
typedef MLSparseMatBaseCL<double>        MLMatrixCL;
} // end of namespace DROPS

#endif
