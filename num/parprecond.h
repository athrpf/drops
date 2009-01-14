//**************************************************************************
// File:    parprecond.h                                                    *
// Content: parallel preconditioners                                       *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   March, 27th 2006                                               *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file parprecond.h
/// \brief Parallel preconditioners
///
/// One important difference between parallel and seriell preconditioners is
/// that some parallel preconditioners need the accumulated diag of the matrix.
/// So the use 'SetDiag(const Vec& diag_acc)' or 'SetDiag(const Mat& A, const ExchangeCL& ex)'
/// to set the accumulated diagonal of the matrix

#ifndef _DROPS_PARPRECOND_H_
#define _DROPS_PARPRECOND_H_

#include "parallel/exchange.h"
#include "num/spmat.h"
#include "num/solver.h"

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
template <typename ExCL>
class ParPreconditioningBaseCL
/** This class can be used as a base class for parallel preconditioners.
    Therefore this class stores a pointer to an ExchangeCL. This is usefull,
    so this class can create an accumulated diagonal out of a matrix by its
    own. This diagonal is only created, if the corresponding matrix has been
    changed.*/
{
  protected:
    ExCL*    ex_;                       ///< exchanging numerical values
    VectorCL diag_;                     ///< accumulated diagonal of corresponding matrix
    size_t   mat_version_;              ///< version of the corresponding matrix
    bool     check_mat_version_;         ///< check matrix version before apply preconditioning (only if DebugNumericC is set)

  public:
    ParPreconditioningBaseCL(ExCL& ex)
      : ex_(&ex), diag_(0), mat_version_(0), check_mat_version_(true) {}
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
        ex_->Accumulate(diag_);
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
    const ExCL& GetEx() const { return *ex_; }
    /// \brief Check matrix version before apply preconditioner in Debug-Mode (DebugNumericC) default
    void CheckMatVersion() { check_mat_version_=true; }
    /// \brief Don't check matrix version before apply preconditioner
    void DoNotCheckMatVersion() { check_mat_version_=false; }
};

// ***************************************************************************
/// \brief Class for performing a preconditioning step with the identity matrix
// ***************************************************************************
template <typename ExCL>
class ParDummyPcCL : public ParPreconditioningBaseCL<ExCL>
{
  private:
    typedef ParPreconditioningBaseCL<ExCL> base_;

  public:
    ParDummyPcCL(ExCL& ex) : base_(ex) {}

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
template <typename ExCL>
class ParJacCL : public ParPreconditioningBaseCL<ExCL>
{
  private:
    typedef ParPreconditioningBaseCL<ExCL> base_;
    using base_::diag_;
    using base_::mat_version_;
    using base_::check_mat_version_;

  protected:
    double  omega_;                                 // overrelaxion-parameter

  public:
    ParJacCL (ExCL& ex, double omega=1) : base_(ex), omega_(omega) {}

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
template <typename ExCL>
class ParJac0CL : public ParJacCL<ExCL>
{
  private:
    typedef ParJacCL<ExCL> base_;
    using base_::mat_version_;
    using base_::diag_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:

    ParJac0CL (ExCL &ex, double omega=1) : base_(ex, omega) {}

    /// \brief Apply preconditioner: one step of the Jacobi-iteration with start vector 0
    template <typename Mat, typename Vec>
    void Apply(const Mat& A, Vec &x, const Vec& b) const
    {
        Assert(mat_version_==A.Version() || !check_mat_version_,
               DROPSErrCL("ParJac0CL::Apply: Diagonal of actual matrix has not been set"),
               DebugNumericC);
        Jacobi0(A,diag_, x, b, omega_, std::fabs(omega_-1.)>DoubleEpsC);
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
/// \brief Class for performing one Step of the Jacobi-Iteration with startvector 0 and own matrix
// ********************************************************************************
template <typename ExCL>
class ParJac0OwnMatCL : public ParJac0CL<ExCL>
{
  private:
    typedef ParJac0CL<ExCL> base_;
    using base_::mat_version_;
    using base_::diag_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:
    ParJac0OwnMatCL (ExCL &ex, double omega=1) : base_(ex, omega) {}

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
template <typename ExCL>
class ParAccPcCL : public ParDummyPcCL<ExCL>
{
  private:
    typedef ParDummyPcCL<ExCL> base_;
    using base_::ex_;

  public:
    ParAccPcCL(ExCL& ex) : base_(ex) {}

    /// \brief Check if the diagonal of the matrix is needed
    bool NeedDiag() const { return false; }
    /// \brief Check if return preconditioned vectors are accumulated after calling Apply
    bool RetAcc() const { return true; }

    /// \brief Apply preconditioner: x <- Accumulate(b)
    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec &x, const Vec& b) const
    {
        x=ex_->GetAccumulate(b);
    }
};

// ********************************************************************************
/// \brief Class for performing one Step of the blocked inexact SSOR0 step
// ********************************************************************************
template <typename ExCL>
class ParSSOR0CL : public ParJacCL<ExCL>
{
  private:
    typedef ParJacCL<ExCL> base_;
    using base_::diag_;
    using base_::mat_version_;
    using base_::omega_;
    using base_::check_mat_version_;

  public:
    ParSSOR0CL (ExCL& ex, double omega=1) : base_(ex, omega) {}

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
