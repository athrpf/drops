/// \file parprecond.h
/// \brief parallel preconditioners
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

/// One important difference between parallel and seriell preconditioners is
/// that some parallel preconditioners need the accumulated diag of the matrix.
/// So the use 'SetDiag(const Vec& diag_acc)' or 'SetDiag(const Mat& A, const ExchangeCL& ex)'
/// to set the accumulated diagonal of the matrix

#ifndef _DROPS_PARPRECOND_H_
#define _DROPS_PARPRECOND_H_

#include "parallel/exchange.h"
#include "num/spmat.h"
#include "num/spblockmat.h"
#include "num/solver.h"
#include "misc/problem.h"

namespace DROPS
{
//***************************************************************************
// implemented parallel methods for preconditioning
//***************************************************************************
template <typename Mat, typename Vec>
void Jacobi(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega);

template <typename Mat, typename Vec>
void Jacobi0(const Mat&, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega);

template <typename Mat, typename Vec>
void SSOR0(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega);


//***************************************************************************
// preconditioning-classes
//***************************************************************************

// ***************************************************************************
/// \brief Base class for parallel preconditioning classes
// ***************************************************************************
class ParPreconditioningBaseCL
/** This class can be used as a base class for parallel preconditioners.
    Therefore this class stores a reference to an IdxDescCL to get a reference
    to the ExchangeCL. Since, the index within the problem classes are lying all
    the time at the same position in memory, we store a reference to the
    IdxDescCL instead of storing a reference to the ExchanegCL.
    This class can create an accumulated diagonal out of a matrix by its own.
    This diagonal is only created, if the corresponding matrix has been
    changed.*/
{
  protected:
    const IdxDescCL& idx_;                ///< index for accessing the ExchangeCL
    VectorCL   diag_;                     ///< accumulated diagonal of corresponding matrix
    size_t     mat_version_;              ///< version of the corresponding matrix
    bool       check_mat_version_;        ///< check matrix version before apply preconditioning (only if DebugNumericC is set)

  public:
    ParPreconditioningBaseCL(const IdxDescCL& idx)
      : idx_(idx), diag_(0), mat_version_(0), check_mat_version_(true) {}
    // default c-tor and de-tor

    /// \brief Set accumulated diagonal of a matrix, that is needed by most of the preconditioners
    template<typename Mat>
    void SetDiag(const Mat& A)
    {
        Comment( (A.Version() == mat_version_ ? "SetDiag, Reusing OLD diagonal\n" : "SetDiag: Creating NEW diagonal\n"), DebugNumericC);
        // Just set diagonal, if the diagonal has changed
        if (A.Version() == mat_version_)
            return;
        if (diag_.size() != A.num_rows())
            diag_.resize(A.num_rows());
        // accumulate diagonal of the matrix
        diag_= A.GetDiag();
        GetEx().Accumulate(diag_);
        // remeber version of the matrix
        mat_version_= A.Version();
    }

    /// \brief Set accumulated diagonal of a matry by the given accumulated diagonal
    void SetDiag(const VectorCL& diag)
    /** \pre diag must be accumulated*/
        { if (diag.size()!=diag_.size()) diag_.resize(diag.size()); diag_=diag; }
    /// \brief Get constant reference on accumulated diagonal of corresponding matrix
    const VectorCL& GetDiag() const { return diag_; }
    /// \brief Get reference on accumulated diagonal of corresponding matrix
    VectorCL&       GetDiag()       { return diag_; }
    /// \brief Get constant reference on exchange class
    const ExchangeCL& GetEx() const { return idx_.GetEx(); }
    /// \brief Check matrix version before apply preconditioner in Debug-Mode (DebugNumericC) default
    void CheckMatVersion() { check_mat_version_=true; }
    /// \brief Don't check matrix version before apply preconditioner
    void DoNotCheckMatVersion() { check_mat_version_=false; }
};

// ***************************************************************************
/// \brief Class for performing a preconditioning step with the identity matrix
// ***************************************************************************
class ParDummyPcCL : public ParPreconditioningBaseCL
{
  private:
    typedef ParPreconditioningBaseCL base_;

  public:
    ParDummyPcCL(const IdxDescCL& idx) : base_(idx) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc()   const { return false; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }
    /// \name Set diagonal of the matrix for consistency
    //@{
    void SetDiag(const VectorCL&) {}        // just for consistency
    template<typename Mat>
    void SetDiag(const Mat&) {}             // just for consistency
    //@}

    /// \brief Apply preconditioner: x <- b
    template <typename Mat, typename Vec>
    void Apply(const Mat& /*A*/, Vec &x, const Vec& b) const
    {
        x=b;
    }
    /// \brief Apply preconditioner of A^T: x <- b
    template <typename Mat, typename Vec>
    void transp_Apply(const Mat&, Vec &x, const Vec& b) const
    {
        x=b;
    }

};

// ***************************************************************************
/// \brief Class for performing one Step of the Jacobi-Iteration
// ***************************************************************************
class ParJacCL : public ParPreconditioningBaseCL
{
  protected:
    typedef ParPreconditioningBaseCL base_;
    using base_::diag_;
    using base_::mat_version_;
    using base_::check_mat_version_;

  protected:
    double  omega_;                                 // overrelaxion-parameter

  public:
    ParJacCL (const IdxDescCL& idx, double omega=1) : base_(idx), omega_(omega) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const   { return false; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return true; }
    /// \brief Get overrelaxation parameter
    double GetOmega() const {return omega_;}

    /// \brief Apply preconditioner: one step of the Jacobi-iteration
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJacCL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }
};

// ********************************************************************************
/// \brief Class for performing one Step of the Jacobi-Iteration with startvector 0
// ********************************************************************************
class ParJac0CL : public ParJacCL
{
  protected:
    typedef ParJacCL base_;
    using base_::mat_version_;
    using base_::diag_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:

    ParJac0CL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \brief Apply preconditioner: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJac0CL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }

    /// \brief Apply preconditioner of A^T: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void transp_Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version()  || !check_mat_version_,
               DROPSErrCL("ParJac0CL::transp_Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }
};

// ********************************************************************************
/// \brief Class for performing one step of the Jacobi-Iteration on a matrix A*A^T
///   with startvector 0
// ********************************************************************************
class ParJacNEG0CL : public ParJac0CL
{
  protected:
    typedef ParJac0CL base_;        ///< base class
    ExchangeMatrixCL  exMat_;       ///< handling of accumulating matrix-entries

    /// \brief Determine the diagonal of A*A^T
    inline void MySetDiag(const MatrixCL& A, const ExchangeCL& RowEx, const ExchangeCL& ColEx);

  public:
    /// \brief Constructor
    ParJacNEG0CL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \name Compute diagonal of matrix AA^T
    //@{
    /// \brief Set known diagonal of A*A^T as diagonal
    void SetDiag(const VectorCL& d) { base_::diag_.resize(d.size()); base_::diag_=d; }

    /// \brief Determine diagonal of a matrix A is not implemented for all MatTs
    template<typename MatT>
    void SetDiag(const MatT&){
        throw DROPSErrCL("ParJacNEG0CL::SetDiag: Not defined for that matrix type.");
    }
    //@}
};

inline void ParJacNEG0CL::MySetDiag(const MatrixCL& A, const ExchangeCL& RowEx, const ExchangeCL& ColEx)
/// This function determines the diagonal of the matrix A*A^T and stores this diagonal
/// in the base class. Therefore, the matrix A has to be accumulated. This is done by
/// the class MatrixExchangeCL.
///
/// Conservatively, we assume that the pattern of the matrix has changed, so the MatrixExchangeCL
/// is created each time, this function is called.
///
/// \todo(par) Determine when to create the communication-pattern more precisely
/// \param A     Compute diagonal of the matrix A*A^T
/// \param RowEx ExchangeCL according to the RowIdx of matrix A
/// \param ColEx ExchangeCL according to the ColIdx of matrix A
{
    // Accumulate matrix
    exMat_.BuildCommPattern(A, RowEx, ColEx);
    MatrixCL Aacc(exMat_.Accumulate(A));
    // Determine diagonal of AA^T
    base_::diag_.resize(Aacc.num_rows());
    for (size_t i = 0; i < Aacc.num_rows(); ++i)
        for (size_t nz = Aacc.row_beg(i); nz < Aacc.row_beg(i + 1); ++nz)
            base_::diag_[i] += Aacc.val(nz) * A.val(nz);
    RowEx.Accumulate(base_::diag_);
    mat_version_= A.Version();
}

/// \brief (Specialization) SetDiag for CompositeMatrixBaseCL<MatrixCL, MatrixCL>
template <>
inline void ParJacNEG0CL::SetDiag(const CompositeMatrixBaseCL<MatrixCL, MatrixCL>& BBT)
{
    MySetDiag(*BBT.GetBlock1(), BBT.GetEx1(), BBT.GetEx0());
}


// ********************************************************************************
/// \brief Class for performing one Step of the Jacobi-Iteration with startvector 0 and own matrix
// ********************************************************************************
class ParJac0OwnMatCL : public ParJac0CL
{
  private:
    typedef ParJac0CL base_;
    using base_::mat_version_;
    using base_::diag_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:
    ParJac0OwnMatCL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \brief Apply preconditioner: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }

    /// \brief Apply preconditioner of A^T: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void transp_Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Jacobi0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }

    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }

    /// \brief Set Diag makes nothing, because own matrix is used
    template<typename Mat>
    void SetDiag(const Mat&) {}

    /// \brief Set new diagonal
    template<typename Mat>
    void SetNewDiag(const Mat& M) { base_::SetDiag(M); }

    /// \brief Set new diagonal (given as an accumulated vector)
    void SetNewDiag(const VectorCL& v) { base_::SetDiag(v); }
};

// ********************************************************************************
/// \brief Class for performing an accumulation as preconditioning
// ********************************************************************************
class ParAccPcCL : public ParDummyPcCL
{
  private:
    typedef ParDummyPcCL base_;

  public:
    ParAccPcCL(const IdxDescCL& idx) : base_(idx) {}

    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }
    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const { return true; }

    /// \brief Apply preconditioner: x <- Accumulate(b)
    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec &x, const Vec& b) const
    {
        x= base_::GetEx().GetAccumulate(b);
    }
};

// ********************************************************************************
/// \brief Class for performing one Step of the blocked inexact SSOR0 step
// ********************************************************************************
class ParSSOR0CL : public ParJacCL
{
  private:
    typedef ParJacCL base_;
    using base_::diag_;
    using base_::mat_version_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:
    ParSSOR0CL (const IdxDescCL& idx, double omega=1) : base_(idx, omega) {}

    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const { return false; }
    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return true; }

    /// \brief Apply preconditioner: one step of the incomplete SSOR-iteration
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJac0CL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        SSOR0(A, diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
    }
};

// ********************************************************************************
/// \brief serial CG as preconditioner
// ********************************************************************************
class ParCGPreCL
{
  private:
    int iters_;
    double tol_;

  public:
    ParCGPreCL(int Steps, double tol) : iters_(Steps), tol_(tol) {}

    inline bool RetAcc() const {return true;}
    inline bool NeedDiag() const {return false;}
    void SetDiag(VectorCL */*diag*/) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        int iter= iters_;
        double tol=tol_;
        x=0.;
        PCG(A, x, b, SSORPcCL(), iter, tol);
        if (ProcCL::MyRank()==0)
            std::cout << "  -- CG as Precond used " << iter <<" steps\n";
    }
};

// ********************************************************************************
/// \brief serial GMRES as preconditioner
// ********************************************************************************
class ParGMResPreCL
{
  private:
    int iters_;
    int restart_;
    double tol_;

 public:
    ParGMResPreCL(int Restart, int Steps, double tol) : iters_(Steps), restart_(Restart), tol_(tol) {}

    inline bool RetAcc() const {return true;}
    inline bool NeedDiag() const {return false;}
    void SetDiag(VectorCL */*diag*/) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        int iter= iters_;
        double tol=tol_;
        GMRES(A, x, b, MultiSSORPcCL(1,1), restart_, iter, tol);
        if (ProcCL::MyRank()==0)
            std::cout << "  -- GMRES as Precond used " << iter <<" steps\n";
    }
};


//***************************************************************************
// Implementations of the methods
//***************************************************************************

/// \brief One step of a Jacobi-iteration
template <typename Mat, typename Vec>
void Jacobi(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega)
        /// \param[in]     A        local distributed coefficients-matrix of the linear equation system
        /// \param[in]     Diag     accumulated form of the diagonalelements of A
        /// \param[in,out] x        startvector (accumulated) and result (distributed)
        /// \param[in]     b        rhs (distributed form)
        /// \param[in]     omega    overrelaxion-parameter
        /// \param[in]     HasOmega flag, if this routine takes care of an overrelaxion-parameter
{
    const size_t n= A.num_rows();
    Vec          y(x.size());

    for (size_t i=0, nz=0; i<n; ++i)
    {
        double sum= b[i];
        for (const size_t end= A.row_beg(i+1); nz<end; ++nz)
            if (A.col_ind(nz) != i)
                sum-= A.val(nz)*x[A.col_ind(nz)];
        if (HasOmega)
            y[i]= (1.-omega)*x[i]+omega*sum/Diag[i];
        else
            y[i]= sum/Diag[i];
    }

    std::swap(x,y);
}

/// \brief One Step of a Jacobi-iteration with startvector 0
template <typename Mat, typename Vec>
void Jacobi0(const Mat&, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega)
        /// \param[in]     Diag     accumulated form of the diagonalelements of A
        /// \param[out]    x        result (distributed form)
        /// \param[in]     b        rhs (distributed form)
        /// \param[in]     omega    overrelaxion-parameter
        /// \param[in]     HasOmega flag, if this routine takes care of an overrelaxion-parameter

{
    const size_t n = x.size();
    for (size_t i=0; i<n; ++i)
    {
        if (HasOmega)
            x[i] = omega * b[i] / Diag[i];
        else
            x[i] = b[i] / Diag[i];
    }
}

// One step of the Symmetric-Gauss-Seidel/SSOR method with start vector 0
template <typename Mat, typename Vec>
void SSOR0(const Mat& A, const Vec& Diag, Vec& x, const Vec& b, const double omega, const bool HasOmega)
{
    const size_t n= A.num_rows();

    for (size_t i=0; i<n; ++i)
    {
        double sum= b[i];
        size_t j= A.row_beg(i);
        for ( ; A.col_ind(j) < i; ++j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= omega*sum/Diag[i];
        else
            x[i]= sum/Diag[i];
    }

    for (size_t i=n; i>0; )
    {
        --i;
        double sum= 0;
        size_t j= A.row_beg(i+1)-1;
        for ( ; A.col_ind(j) > i; --j)
            sum-= A.val(j)*x[A.col_ind(j)];
        if (HasOmega)
            x[i]= (2.-omega)*x[i]+omega*sum/Diag[i];
        else
            x[i]+= sum/Diag[i];
    }
}
} // end of namespace DROPS

#endif
