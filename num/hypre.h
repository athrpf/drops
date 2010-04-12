/// \file hypre.h
/// \brief Interface to HYPRE (high performance preconditioners), Lawrence Livermore Nat. Lab. 
/// \author LNM RWTH Aachen: Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_HYPRE_H
#define DROPS_HYPRE_H

#ifdef _HYPRE

#include "misc/utils.h"
#include "num/spmat.h"
#include "num/solver.h"
#include "num/parsolver.h"
#include "num/parprecond.h"

namespace DROPS {

/// \brief Distribution of a index, needed by HYPRE
class HypreIndexCL
{
  private:
    const IdxDescCL&   idx_;                    ///< description of finite elements
    std::vector<int>   idx_beg_;                ///< start index of each process
    VectorBaseCL<int>  global_idx_;             ///< map local index to a global index
    std::vector<int>   local_idx_;              ///< map global index to a local index (elements idx_beg_[MyRank()]..idx_beg_[MyRank()+1]-1)
    std::vector<int>   exclusive_idx_;          ///< map global index to a local exclusive index (elements idx_beg_[MyRank()]..idx_beg_[MyRank()+1]-1)
    size_t             num_exclusive_;          ///< number of exclusive indices
    VectorBaseCL<int>  exclusive_global_idx_;   ///< global indices of exclusive dof

  public:
    /// \brief Construct a HYPRE index, i.e., store information about each index
    HypreIndexCL( const IdxDescCL& idx)
        : idx_(idx), idx_beg_( ProcCL::Size()+1), global_idx_(0), local_idx_(0),
          exclusive_idx_(0), num_exclusive_(0), exclusive_global_idx_(0) { Generate(); }

    /// \brief Rebuild all arrays to represent the indices
    void Generate();
    /// \brief Clear arrays
    void Clear();

    /// \brief Get a constant list of global exclusive indices
    const int* GetGlobalExList() const { return Addr(exclusive_global_idx_); }
    /// \brief Get a list of global exclusive indices
    int* GetGlobalExList() { return Addr(exclusive_global_idx_); }
    //// \brief Get a list of global exclusive indices as a vector
    const VectorBaseCL<int>& GetGlobalExListVector() const { return exclusive_global_idx_; }

    /// \brief Get first index of this process
    inline int  GetLower() const { return idx_beg_[ProcCL::MyRank()]; }
    /// \brief Get last index of this process
    inline int  GetUpper() const { return idx_beg_[ProcCL::MyRank()+1]-1; }
    /// \brief Get number of indices this process is resposnsible for
    int  GetNumExclusive() const { return (int) num_exclusive_; }
    /// \brief Check if a local index is exclusive
    inline bool IsLocExclusive( int loc_i) const { return idx_.GetEx().IsExclusive( (size_t)loc_i); }
    /// \brief Check if a global index is exclusibe
    bool IsGlobExclusive( int glob_i) const { return glob_i>=GetLower() && glob_i<=GetUpper(); }
    /// \brief Get exclusive process of a local index
    inline int  GetLocExclusiveProc( int loc_i) const { return idx_.GetEx().GetExclusiveProc( (size_t)(loc_i)); }
    /// \brief Get global index of a local index
    int  GetGlobalIdx( int loc_i) const { return global_idx_[loc_i]; }
    /// \brief Get local index of a global index
    int  GetLocalIdx( int glob_i) const { return IsGlobExclusive(glob_i) ? local_idx_[glob_i-GetLower()] : -1;}
    /// \brief Get exclusive index of a global index
    int  GetExclusiveNum( int glob_i) const { return IsGlobExclusive(glob_i) ? exclusive_idx_[glob_i-GetLower()]: -1; }
    /// \brief Get number of neighbors
    int  GetNumNeighs() const { return idx_.GetEx().GetNumNeighs(); }
    /// \brief Get access to ExchangeCL
    const ExchangeCL& GetEx() const { return idx_.GetEx(); }
};

/// \brief Interface for a DROPS matrix to generate a HYPRE matrix
class HypreMatrixCL
{
  private:
    const MatrixCL      M_;             ///< DROPS matrix
    HypreIndexCL        colidx_;        ///< HYPRE index to respresent the column
    HypreIndexCL        rowidx_;        ///< HYPRE index to respresent the row
    HYPRE_IJMatrix      ijMat_;         ///< HYPRE matrix
    HYPRE_ParCSRMatrix  parMat_;        ///< HYPRE parallel crs matrix
    std::vector<int>    cols_,          ///< column indices of all non-zeros (the row indices are given by rowidx_)
                        ncols_,         ///< number of column entries per row
                        nOffCol_,       ///< number of column entries per row whose corresponding vector entries live on another process
                        nMyCol_;        ///< number of column entries per row whose corresponding vector entries live on this process
    std::vector<double> vals_;          ///< values of the non-zeros
    ProcCL::DatatypeT   sendMPItype_;   ///< type for sending a NonZeroST
    int                 tag_;           ///< tag used by this class

    /// \brief Datatype for sending a non-zero
    struct NonZeroST
    {
        int    row;     ///< global row
        int    col;     ///< global col
        double value;   ///< local value of the non-zero

        NonZeroST() : row(-1), col(-1), value(0.) {}                        ///< Generate a dummy element
        NonZeroST( int r, int c, double v) : row(r), col(c), value(v) {}    ///< Generate a non-zero element
        bool IsDummy() const { return row<0; }                              ///< Check if structure holds a valid entry
    };

    /// \brief Generate the MPI datatype for sending a NonZeroST
    void CreateSendMPIType();
    /// \brief Generate the HYPRE matrix 
    void Init();

  public:
    HypreMatrixCL( const MatrixCL& M, const IdxDescCL& rowIdx, const IdxDescCL& colIdx)
        : M_(M), colidx_(colIdx), rowidx_(rowIdx), sendMPItype_(ProcCL::NullDataType), tag_(1001) 
    { Init();  Generate(); }
    ~HypreMatrixCL() { ProcCL::Free(sendMPItype_); }

    /// \brief Update the matrix
    void Generate();
    /// \brief Clear arrays
    void Clear();

    /// \name Get reference to index classes
    //@{
    const HypreIndexCL& GetRowIdx() const { return rowidx_; }
    const HypreIndexCL& GetColIdx() const { return colidx_; }
    HypreIndexCL&       GetRowIdx()       { return rowidx_; }
    HypreIndexCL&       GetColIdx()       { return colidx_; }
    //@}

    /// \name Get reference to Hypre matrix
    //@{
    const HYPRE_IJMatrix& GetHypreIJMatrix() const    { return ijMat_; }
    HYPRE_IJMatrix& GetHypreIJMatrix()                { return ijMat_; }
    const HYPRE_ParCSRMatrix& operator() (void) const { return parMat_; }
    HYPRE_ParCSRMatrix& operator() (void)             { return parMat_; }
    //@}
};

/// \brief Respresent a vector in HYPRE format
class HypreVectorCL
{
  private:
    VectorCL&           v_;         ///< DROPS vector
    const HypreIndexCL& idx_;       ///< corresponding indexing
    HYPRE_IJVector      ijVec_;     ///< type of the HYPRE vector
    HYPRE_ParVector     parVec_;    ///< HYPRE vector
    
    void Init();
    
  public:
    HypreVectorCL( VectorCL& v, const HypreIndexCL& idx)
      : v_(v), idx_(idx) { Init(); Generate(); }
    ~HypreVectorCL();
    
    /// \brief generate the HYPRE vector
    void Generate();
    /// \brief generate an accumulated DROPS vector from the HYPRE vector
    void Retrieve();
    /// \brief Access to the HYPRE vector
    HYPRE_ParVector operator() (void) const { return parVec_; }
};

/// \brief Interface to the HYPRE Boomer AMG solver
class HypreAMGSolverCL : public SolverBaseCL
{
  private:
    typedef SolverBaseCL baseT;

    const IdxDescCL&     idx_;
    HYPRE_Solver         solver_;

    /// \brief Initialize the solver
    void Init();
    ParDummyPcCL dummy_;

  public:
    HypreAMGSolverCL( const IdxDescCL& idx, int maxiter= 100, double tol=1e-7)
        : baseT(maxiter, tol), idx_(idx), dummy_(idx_) { Init(); }
    ~HypreAMGSolverCL();
    
    void SetTol( double tol)   { _tol= tol; HYPRE_BoomerAMGSetTol( solver_, tol); }
    void SetMaxIter( int iter) { _maxiter= iter; HYPRE_BoomerAMGSetMaxIter( solver_, iter); }
      
    void Setup( const HypreMatrixCL& A, const HypreVectorCL& x, const HypreVectorCL& b) const;
    /// \name Solve a system of linear equations
    //@{
    void SetupAndSolve( const MatrixCL& A, VectorCL& x, const VectorCL& b, const IdxDescCL& idx) const;
    void Solve( const HypreMatrixCL& A, HypreVectorCL& x, const HypreVectorCL& b) const;
    void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b)   { SetupAndSolve(A,x,b,idx_); }
    void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b) { SetupAndSolve(A.GetFinest(),x,b,idx_); }
    ParDummyPcCL GetPC() { return dummy_; }
    //@}
};   

}    // end of namespace DROPS

#endif      // of _HYPRE
#endif // DROPS_HYPRE_H

