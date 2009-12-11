/// \file parsolver.h
/// \brief Parallel basic iterative solvers
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

#ifndef _DROPS_PARSOLVER_H_
#define _DROPS_PARSOLVER_H_

#include "num/spmat.h"
#include "num/parlanczos.h"
#include "num/solver.h"
#include "parallel/parallel.h"
#include "misc/problem.h"
#include "parallel/exchange.h"
#include <algorithm>

// for file i/o
//#define FILE_OUT
#ifdef FILE_OUT
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#endif

// *****************************************************************************
// *                         T O C                                             *
// *****************************************************************************
// * Definition of parallel solver functions                                    *
// * Definition of parallel solver base classes                                 *
// *   - ParSolverBaseCL                                                        *
// *   - ParPreSolverBaseCL                                                     *
// * Definition of parallel solver classes                                      *
// *   - Conjugate Gradients (ParCGSolverCL)                                    *
// *   - preconditioned Conjugate Gradients (ParPCGSolverCL)                    *
// *   - preconditioned Generalized Minimal Residual  (ParPreGMResSolverCL)     *
// *   - preconditioned Bi-Conjugate Gradient Stabilized (ParBiCGSTABSolverCL)  *
// *   - preconditioned Generalized Conjugate Residuals (ParPreGCRSolverCL)     *
// *   - (preconditioned) Quasi Minimal Residual (ParQMRSolverCL)               *
// * Declaration of parallel solver functions                                   *
// ******************************************************************************

namespace DROPS
{

// const double BreakDownC=1e-35;

//***************************************************************************
// implemented parallel methods for solving a linear system
//***************************************************************************
template <typename Mat, typename Vec, typename ExCL>
bool ParCG(const Mat& A,Vec& x_acc,const Vec& b, const ExCL& ExX, int& max_iter, double& tol);

template <typename Mat, typename Vec, typename PreCon, typename ExCL>
        bool ParPCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M, int& max_iter, double& tol, bool measure_relative_tol=false);

template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParAccurPCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M, int& max_iter, double& tol, bool measure_relative_tol=false, std::ostream* output=0);

// Preconditioned GMRES with Gramm-Schmidt
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPreGS_GMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                    int m, int& max_iter, double& tol, bool measure_relative_tol =true);

// Preconditioned GMRES with modified Gramm-Schmidt
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPreMGS_GMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                     int m, int& max_iter, double& tol, bool measure_relative_tol =true);

// Preconditioned GMRES with unmodified Gramm-Schmidt and accurate inner product computation
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParAccurPreGS_GMRES(const Mat& A, Vec& x, const Vec& b, const ExCL& ExX, PreCon& M,
                         int m, int& max_iter, double& tol, bool measure_relative_tol_=true);

// Preconditioned GMRES
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
              int m, int& max_iter, double& tol, bool measure_relative_tol=true, bool useAcc=true,
              PreMethGMRES method=LeftPreconditioning);

// Preconditioned GMRES with modifications for better scalability.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParModGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                 int m, int& max_iter, double& tol, bool measure_relative_tol=true, bool useAcc=true,
                 bool useMGS=false, PreMethGMRES method=LeftPreconditioning);

// BiCGSTAB
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParBiCGSTAB(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M, int& max_iter, double& tol, bool measure_relative_tol=true);

// Preconditioned GCR with truncation with high accuracy
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPGCR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
             int m, int& max_iter, double& tol, bool measure_relative_tol =true, bool useAcc=true);

// Preconditioned GCR with truncation and modifications to reduce sync-points
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParModPGCR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
             int m, int& max_iter, double& tol, bool measure_relative_tol =true);

// Preconditioned GCR with truncation with higher accuracy and modifications to reduce sync-points
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParModAccurPGCR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                  int m, int& max_iter, double& tol,bool measure_relative_tol /*=true*/);

// QMR-Solver
template <typename Mat, typename Vec, typename Lanczos, typename ExCL>
bool ParQMR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, Lanczos lan,
            int& max_iter, double& tol, bool measure_relative_tol);



//***************************************************************************
//                      S O L V E R   B A S E  C L A S S E S
//***************************************************************************

// ***************************************************************************
/// \brief Parallel base solver class
// ***************************************************************************
class ParSolverBaseCL : public SolverBaseCL
/** All parallel solvers needs additionally to SolverBaseCL a class for
    exchanging numerical values. Since these information is stored as a part
    of the IdxDescCL, the parallel solvers stores a reference to the IdxDescCL
    as well.
*/
{
  private:
    const IdxDescCL* idx_;      // for getting the ExchangeCL
    bool             acc_;      // use accurate version for computing norms and inner products

  public:
    /// \brief Constructor
    ParSolverBaseCL(int maxiter, double tol, const IdxDescCL& idx, bool rel= false, bool acc= true, std::ostream* output=0)
      : SolverBaseCL(maxiter, tol, rel, output), idx_(&idx), acc_(acc) {}
    /// \brief Constructor, that does not initialize the index description
    ParSolverBaseCL(int maxiter, double tol, bool rel= false, bool acc= true, std::ostream* output=0)
      : SolverBaseCL(maxiter, tol, rel, output), idx_(0), acc_(acc) {}

    /// \brief Aks for ExchangeCL
    const ExchangeCL& GetEx() const {
        Assert(idx_, DROPSErrCL("ParSolverBaseCL::GetEx: Index not set, do you want to use an ExchangeBlockCL?"), DebugParallelNumC);
        return idx_->GetEx();
    }
    bool              Accurate()          const { return acc_; }              ///< Check if accurate version of inner products and norms are used
    void              SetAccurate(bool a)       { acc_= a; }                  ///< Set to accurate version
};

// ***************************************************************************
/// \brief Parallel preconditioned base solver class
// ***************************************************************************
template <typename PC>
class ParPreSolverBaseCL : public ParSolverBaseCL
/** See ParSolverBaseCL with an preconditioner and a flag if the residual should be measured relative.*/
{
  private:
    typedef ParSolverBaseCL base;

  protected:
    PC   &pc_;                                      // preconditioner

  public:
    typedef PC PrecondT;                            ///< Preconditioner

    /// \brief Constructor
    ParPreSolverBaseCL(int maxiter, double tol, const IdxDescCL& idx, PC& pc, bool rel=true, bool acc= true, std::ostream* output=0)
      : base(maxiter, tol, idx, rel, acc, output), pc_(pc) {}
    /// \brief Constructor, that does not initialize the index description
    ParPreSolverBaseCL(int maxiter, double tol, PC& pc, bool rel=true, bool acc= true, std::ostream* output=0)
      : base(maxiter, tol, rel, acc, output), pc_(pc) {}

    PC& GetPC()             {return pc_;}          ///< return reference on preconditioner
    const PC& GetPC() const {return pc_;}          ///< return constant reference on preconditioner

    void SetPC(PC& pc) { pc_= pc; }                ///< Set preconditioner
};

//***************************************************************************
//                      S O L V E R  C L A S S E S
//***************************************************************************

// ***************************************************************************
/// \brief Parallel CG-Solver class
// ***************************************************************************
class ParCGSolverCL : public ParSolverBaseCL
{
  private:
    typedef ParSolverBaseCL base;

  public:
    /// \brief Constructor for parallel CG-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products.*/
    ParCGSolverCL(int maxiter, double tol, const IdxDescCL& idx, bool rel= false, bool acc=true, std::ostream* output=0)
      : base(maxiter, tol, idx, rel, acc, output) {}

    /// \brief Solve a linear equation system with Conjugate-Gradients-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// CG algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        ParCG(A, x, b, base::GetEx(), base::_iter, base::_res);
    }
};


// ***************************************************************************
/// \brief Parallel preconditioned CG-Solver class
// ***************************************************************************
template <typename PC>
class ParPCGSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;

  public:
    /// \brief Constructor for the parallel preconditioned CG Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative and
        with \a acc the inner products are determined with accure variant (see ExchangeCL). */
    ParPCGSolverCL(int maxiter, double tol, const IdxDescCL &idx, PC& pc, bool rel=false, bool acc=true, std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, acc, output) {}

    /// \brief Solve a linear equation system with Conjugate Gradients-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned CG algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        if (base::Accurate())
            ParAccurPCG(A, x, b, base::GetEx(),  base::GetPC(), base::_iter, base::_res, base::rel_, base::output_);
        else
            ParPCG(A, x, b, base::GetEx(),  base::GetPC(), base::_iter, base::_res, base::rel_);
    }
};

// ***************************************************************************
/// \brief Parallel preconditioned GMRES-Solver class
// ***************************************************************************
template <typename PC>
class ParPreGMResSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;
    int          restart_;                  // number of iterations before restart
    bool         useModGS_;                 // which Gramm-Schmidt method should be used to compute Krylov basis
    PreMethGMRES method_;                   // left or right preconditioning
    bool         mod_;                      // use modified variant for better scalability


  public:
    /// \brief Constructor of the parallel preconditioned GMRES-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. After \a restart steps, a restart is performed. The ExCL
        \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative and
        with \a acc the inner products are determined with accure variant (see ExchangeCL).
        (this configuration needs less memory!). By setting \a ModGS the modified Gramm-Schmidt
        algorithm is used for the Arnoldi method.*/
    ParPreGMResSolverCL(int restart, int maxiter, double tol, const IdxDescCL& idx, PC &pc,
                        bool rel=true, bool acc=true, bool ModGS=false,
                        PreMethGMRES method=LeftPreconditioning, bool mod=true,
                        std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, acc, output),
        restart_(restart), useModGS_(ModGS), method_(method), mod_(mod) {}

    int  GetRestart()           const { return restart_; }  ///< number of iterations before restart
    void SetRestart(int restart)      { restart_=restart; } ///< set number of iterations before restart

    /// \brief Solve a linear equation system with a preconditioned Generalized Minimal Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GMRES algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;

        if (mod_)
            ParModGMRES(A, x, b, base::GetEx(), base::GetPC(), restart_,
                        base::_iter, base::_res, base::GetRelError(), base::Accurate(),
                        useModGS_, method_);
        else
            ParGMRES(A, x, b, base::GetEx(),  base::GetPC(), restart_,
                     base::_iter, base::_res, base::GetRelError(), base::Accurate(),
                     method_);
    }
};

// ***************************************************************************
/// \brief Parallel BiCGSTAB-Solver class
// ***************************************************************************
template<typename PC>
class ParBiCGSTABSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;

  public:
    /// \brief Constructor of the parallel preconditioned BiCGStab-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative. */
    ParBiCGSTABSolverCL(int maxiter, double tol, const IdxDescCL &idx, PC& pc, bool rel=true,
                        bool acc=true, std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, acc, output) {}

    /// \brief Solve a linear equation system with preconditioned Bi-Conjugate Gradient Stabilized-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec& b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned BiCGStab algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        ParBiCGSTAB(A, x, b, base::GetEx(), base::GetPC(), base::_iter, base::_res, base::GetRelError());
    }
};

// ***************************************************************************
/// \brief Parallel GCR-Solver class with preconditioning
// ***************************************************************************
template <typename PC>
class ParPreGCRSolverCL : public ParPreSolverBaseCL<PC>
{
  private:
    typedef ParPreSolverBaseCL<PC> base;
    int  trunc_;
    bool mod_;

  public:
    /// \brief Constructor of the parallel preconditioned GCR-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. \a trunc vectors are used to span the Krylov subspace. The ExCL
        \a ex is used to do parallel inner products. \a pc is
        the given preconditioner. If \a rel is given, the residual is computed relative and
        with \a acc the inner products are determined with accure variant (see ExchnageCL).
        If \a mod is set than a modified variant for computing the Krylov subspace
        is used, to reduce sync-points.
        \todo (of) <b>truncation strategy with modified GCR do not work!</b>
    */
    ParPreGCRSolverCL(int trunc, int maxiter, double tol, const IdxDescCL& idx, PC &pc, bool mod=false,
                      bool rel=true, bool acc=true, std::ostream* output=0)
      : base(maxiter, tol, idx, pc, rel, acc, output), trunc_(trunc), mod_(mod) {}
    /// \brief Constructor, that does not initialize the index description
    ParPreGCRSolverCL(int trunc, int maxiter, double tol, PC &pc, bool mod=false,
                      bool rel=true, bool acc=true, std::ostream* output=0)
      : base(maxiter, tol, pc, rel, acc, output), trunc_(trunc), mod_(mod) {}

    /// \brief Solve a linear equation system with preconditioned Generalized Conjugate Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GCR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        if (mod_){
            if (base::Accurate())
                ParModAccurPGCR(A, x, b, base::GetEx(), base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError(), base::output_);
            else
                ParModPGCR(A, x, b, base::GetEx(), base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError());
        }
        else
            ParPGCR(A, x, b, base::GetEx(), base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError(), base::Accurate());
    }

    /// \brief Solve a linear equation system with preconditioned Generalized Conjugate Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b, const ExchangeBlockCL& ex)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// preconditioned GCR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        if (mod_){
            if (base::Accurate())
                ParModAccurPGCR(A, x, b, ex, base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError(), base::output_);
            else
                ParModPGCR(A, x, b, base::GetEx(), base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError());
        }
        else
            ParPGCR(A, x, b, ex, base::GetPC(), trunc_, base::_iter, base::_res, base::GetRelError(), base::Accurate());
    }

};

// ***************************************************************************
/// \brief Parallel QMR-Solver class
// ***************************************************************************
template<typename Lanczos>
class ParQMRSolverCL : public ParSolverBaseCL
{
  private:
    Lanczos *lan_;
    typedef ParSolverBaseCL base;

  public:
    /// \brief Constructor of the parallel preconditioned QMR-Solver
    /** Tries to solve a linear equation system within \a maxiter steps with
        accuracy \a tol. The ExCL \a ex is used to do parallel inner products. The
        Lanczos-Algorithm to compute the bi-orthogonal-basis is given by a Lanczos class \a lan.
        If \a measure_relative_tol is given, the residual is computed relative.*/
    ParQMRSolverCL(int maxiter, double tol, const IdxDescCL& idx, Lanczos &lan, bool rel=true, bool acc= true, std::ostream* output=0) :
        base(maxiter, tol, idx, rel, acc, output), lan_(&lan) {}

    /// \brief Solve a linear equation system with Quasi Minimal Residuals-Method
    template <typename Mat, typename Vec>
      void Solve(const Mat& A, Vec& x, const Vec &b)
    /// Solve the linear equation system with coefficient matrix \a A and rhs \a b iterative with
    /// QMR algorithm, uses \a x as start-vector and result vector.
    /// \post x has accumulated form
    {
        base::_res  = base::_tol;
        base::_iter = base::_maxiter;
        ParQMR(A, x, b, base::GetEx(), *lan_, base::_iter, base::_res, base::GetRelError());
    }
};

//***************************************************************************
// Implementations of the methods
//***************************************************************************

/// \brief Parallel CG-Algorithm
template <typename Mat, typename Vec, typename ExCL>
bool ParCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, int& max_iter, double& tol)
    /// \param[in]     A        local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc    start vector and the solution in accumulated form
    /// \param[in]     b        rhs of the linear equation system (distributed form)
    /// \param[in]     ExX      ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] max_iter IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol      IN: tolerance for the residual, OUT: residual
    /// \return                 convergence within max_iter iterations
{
    Vec r ( b - A*x_acc ), r_acc(r);

    double alpha = ExX.ParDotAcc(r_acc,r);
    Vec p_acc (r_acc);

    tol*= tol;

    if (alpha<=tol)
    {
        if (alpha<0)
            std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of resid in CG because of accumulation!" << std::endl;
        tol= std::sqrt(alpha);
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        const Vec    v      = A*p_acc;
        const double gamma  = ProcCL::GlobalSum(dot( v, p_acc));
        const double lambda = alpha/gamma;
        double       beta   = alpha;

        x_acc += lambda * p_acc;
        r     -= lambda*v;

        alpha = ExX.Norm_sq_Acc(r_acc,r);

        if (alpha<=tol)
        {
            if (alpha<0)
                std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of resid in CG because of accumulation!" << std::endl;

            tol= std::sqrt(std::fabs(alpha));
            max_iter= i;
            return true;
        }

        p_acc=r_acc + (alpha/beta) * p_acc;
    }
    tol= std::sqrt(std::fabs(alpha));

    return false;
}


/// \brief Parallel preconditioned CG-Algorithm
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
            PreCon& M, int& max_iter, double& tol, bool measure_relative_tol)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol measure resid relative
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n= b.size();
    Vec p_acc(n), z_acc(n), q(n), r( b - A*x_acc), r_acc(r);

    double rho,
           rho_1= 0.0,
           resid= std::sqrt(ExX.ParDotAcc(r_acc,r)),
           normb= ExX.Norm(b);

    if (normb == 0.0 || measure_relative_tol == false)
        normb= 1.0;
    resid = resid/normb;


    if (resid<=tol){
        tol= resid;
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        M.Apply(A, z_acc, r);
        if (M.RetAcc())
            rho = ProcCL::GlobalSum(dot(z_acc,r));
        else
            rho= ExX.ParDotAcc( z_acc, r);

        if (i == 1)
            p_acc= z_acc;
        else
            p_acc = z_acc + (rho/rho_1)*p_acc;

        q= A*p_acc;
        const double lambda = ProcCL::GlobalSum(dot(p_acc, q));
        const double alpha  = rho/lambda;

        x_acc += alpha * p_acc;
        r     -= alpha * q;

        const double res= ExX.Norm_sq_Acc(r_acc, r);
        resid= std::sqrt(res<0 ? 0 : res) / normb;

        if (resid<=tol){
            if (res<0){
                std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of resid in PCG because of accumulation!"
                          << "\n   Please use accurate version."<<std::endl;
                resid=0;
            }

            tol= resid;
            max_iter= i;
            return true;
        }
        rho_1= rho;
    }
    tol= resid;
    return false;
}

/// \brief Parallel preconditioned CG-Algorithm with higher accuracy as ParPCG
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParAccurPCG(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
                   PreCon& M, int& max_iter, double& tol, bool measure_relative_tol, std::ostream* output)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol measure resid relative
    /// \param[in]     output               write information onto output stream
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n= b.size();
    Vec p_acc(n), z(n), z_acc(n), q(n), q_acc(n), r( b - A*x_acc), r_acc(r);


    double rho,
           rho_1= 0.0,
           resid= ExX.Norm(r, false, true, &r_acc), //std::sqrt(ExX.AccNorm_sq(r,r_acc)),      // accumulation of r
           normb= ExX.Norm(b);

    if (normb == 0.0 || measure_relative_tol == false)
        normb= 1.0;
    resid = resid/normb;

    if (resid<=tol){
        tol= resid;
        max_iter= 0;
        return true;
    }

    for (int i=1; i<=max_iter; ++i)
    {
        if (M.RetAcc()){
            M.Apply(A, z_acc, r);
            rho = ExX.ParDot(r_acc, true, z_acc, true);     //GlobalSum(ExX.LocAccDot(z_acc,r_acc));
        }
        else{
            M.Apply(A, z, r);
            rho= ExX.ParDot(z, false, r_acc, true, true, &z_acc);  //ExX.AccParDot( z, r_acc, z_acc);           // accumulation of z_acc
        }

        if (i == 1)
            p_acc= z_acc;
        else{
            p_acc = z_acc + (rho/rho_1)*p_acc;              // z_xpay(p_acc, z_acc, (rho/rho_1), p_acc);
        }

        q= A*p_acc;
        const double lambda = ExX.AccParDot(q,p_acc,q_acc); // accumulation of q
        const double alpha  = rho/lambda;

        x_acc += alpha * p_acc;                             //axpy(alpha, p_acc, x_acc)
        r     -= alpha*q;                                   //axpy(-alpha, q, r)

        resid= ExX.Norm(r, false, true, &r_acc);//std::sqrt(ExX.AccNorm_sq(r, r_acc));         // accumulation of r_acc
        resid= resid/normb;

        if (output){
            if (ProcCL::IamMaster())
                (*output) << "ParAccurPCG: "<<i<<": residual "<<resid<<std::endl;
        }
        if (resid<=tol)
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


/// \brief Computes an orthogonal vector on i vectors by the standard Gramm-Schmidt method
template <typename Vec, typename ExCL>
void StandardGrammSchmidt(DMatrixCL<double>& H,
                          Vec& w, bool acc_w, const std::vector<Vec>& v, bool acc_v,
                          int i, const ExCL& ex, Vec& tmpHCol)
/// This routine uses only one synchronization point.
/// \param H       Hessenberg matrix
/// \param w       orthogonal vector
/// \param acc_w   is w accumulated
/// \param v       vectors used to orthogonalize w
/// \param acc_v   is v accumulated
/// \param i       number of vectors v
/// \param ex      class to accumulate a vector
/// \param tmpHCol temporary vector to store a column of H (as parameter so no mem is allocated in this routine)
{
    if (acc_v&&acc_w){                          // both vectors are accumulated
        for (int k=0; k<=i; ++k)
            tmpHCol[k] = ex.LocDot(w,true,v[k],true,/*accur*/true);
    }
    else if (!acc_v&&!acc_w){                   // one of the both vectors have to be accumulated: really bad!
        for (int k=0; k<=i; ++k)
            tmpHCol[k] = ex.LocDot(w,false,v[k],false,/*accur*/false);
    }
    else        // update of w do only works on same types
        throw DROPSErrCL("StandardGrammSchmidt: Cannot do Gramm Schmidt on that kind of vectors!");

    // Syncpoint!
    ProcCL::GlobalSum(Addr(tmpHCol), H.GetCol(i), i+1);
    for (int k=0; k<=i; ++k)
        w -= H(k,i)*v[k];
}


/// \brief Computes an orthogonal vector on i vectors by the modified Gramm-Schmidt method
template <typename Vec, typename ExCL>
void ModifiedGrammSchmidt(DMatrixCL<double>& H, Vec& w, bool acc_w, const std::vector<Vec>& v, bool acc_v, int i, const ExCL& ex)
/// \param H     Hessenberg matrix
/// \param w     orthogonal vector
/// \param acc_w is w accumulated
/// \param v     vectors used to orthogonalize w
/// \param acc_v is v accumulated
/// \param i     number of vectors v
/// \param ex    class to accumulate a vector
{
    if (acc_v&&acc_w){
        for (int k=0; k<=i; ++k){
            H( k, i)= ex.ParDot( w, true, v[k], true, /*useAccur*/true);
            w-= H( k, i)*v[k];
        }
    }
    else if (!acc_v&&!acc_w){
        for (int k=0; k<=i; ++k){
            H( k, i)= ex.ParDot( w, false, v[k], false, /*useAccur*/false);
            w-= H( k, i)*v[k];
        }
    }
    else
        throw DROPSErrCL("StandardGrammSchmidt: Cannot do Gramm Schmidt on that kind of vectors!");
}

/// \brief Parallel GMRES-method.
///
/// This method is the same algorithm as the serial algorithm. For performance issues take the ParModGMRES
/// procedure.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                int m, int& max_iter, double& tol,
                bool measure_relative_tol, bool useAcc, PreMethGMRES method)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in]     m                    number of steps after a restart is performed
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     useAcc               use accur variant for performing inner products and norms or do not use accure variant
    /// \param[in]     method               left or right preconditioning (see solver.h for definition and declaration)
    /// \return  convergence within max_iter iterations
    /// \pre     the preconditioner should be able to handle a accumulated b
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    m= (m <= max_iter) ? m : max_iter; // m > max_iter only wastes memory.

    // The preconditioner can have one of the following properties by calling M.Apply(A,x,b) with b distributed:
    // 1. x has distributed form    (like Dummy, Jacobi, BlockSSOR)
    // 2. x has accumulated form    (like SolverAsPreCL)
    // Unfortunately there are a lot of differences, so we have to implement two strategies

    if (!M.RetAcc()){                                   // case 1
        if (method==LeftPreconditioning)
            throw DROPSErrCL("ParGMRES: Left preconditioning does not work!");
        /// \todo (of) Left preconditioning funktioniert nicht!
        DMatrixCL<double> H( m, m);
        Vec               s( m), cs( m), sn( m),
                          w( b.size()), r( b.size());
        std::vector<Vec>  v( m);
        double            beta, normb, resid;
        Vec z(x_acc.size()), t(x_acc.size());
        for (int i= 0; i < m; ++i)
            v[i].resize( b.size());

        if (method == RightPreconditioning)
        {
            r= b - A*x_acc;
            beta = ExX.Norm(r, false, useAcc);
            normb= ExX.Norm(b, false, useAcc);
        }
        else
        {
            M.Apply( A, r, Vec( b - A*x_acc));
            beta= ExX.Norm(r, false, useAcc);
            M.Apply( A, w, b);
            normb= ExX.Norm(w, false, useAcc);
        }
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
            for (i= 0; i<m-1 && j<=max_iter; ++i, ++j) {
                if (method == RightPreconditioning)
                {
                    M.Apply( A, w, v[i]);
                    w=A*ExX.GetAccumulate(w);
                }
                else
                    M.Apply( A, w, A*ExX.GetAccumulate(v[i]));
                for (int k= 0; k <= i; ++k ) {
                    H( k, i)= ExX.ParDot( w, false, v[k], false, useAcc);
                    w-= H( k, i)*v[k];
                }

                if (i == m - 1) break;

                H( i + 1, i)= ExX.Norm(w, false, useAcc);
                v[i + 1]= w*(1.0/H( i + 1, i));

                for (int k= 0; k < i; ++k)
                    GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

                GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

                resid= std::abs( s[i+1])/normb;
                if (resid <= tol) {
                    if (method == RightPreconditioning)
                    {
                        z=0.;
                        GMRES_Update( z, i, H, s, v);
                        M.Apply( A, t, z);
                        x_acc+=ExX.GetAccumulate(t);
                    }
                    else{
                        z=0.;
                        GMRES_Update( z, i, H, s, ExX.GetAccumulate(v));
                        x_acc += z;
                    }
                    tol= resid;
                    max_iter= j;
                    return true;
                }
            }

            if (method == RightPreconditioning)
            {
                z=0.;
                GMRES_Update( z, i-1, H, s, v);
                M.Apply( A, t, z);
                x_acc+=ExX.GetAccumulate(t);
                r= b - A*x_acc;
            }
            else
            {
                GMRES_Update( x_acc, i-1, H, s, ExX.GetAccumulate(v));
                M.Apply( A, r, Vec( b - A*x_acc));
            }
            beta=ExX.Norm(r, false, useAcc);
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
    else                                                   // case 2: Apply returns accumulated vector
    {
        DMatrixCL<double> H( m, m);
        Vec               s( m), cs( m), sn( m),
                          w( b.size()), r( b.size());
        std::vector<Vec>  v( m);
        double            beta, normb, resid;
        Vec z(x_acc.size()), t(x_acc.size());
        for (int i= 0; i < m; ++i)
            v[i].resize( b.size());

        if (method == RightPreconditioning)
        {
            Vec r_tmp(b - A*x_acc);
            beta = ExX.Norm(r_tmp, false, useAcc, &r);
            normb= ExX.Norm(b, false, useAcc);
        }
        else
        {
            M.Apply( A, r, Vec( b - A*x_acc));
            beta = ExX.Norm(r, true, true);                                 // This works only with useAcc==true
            M.Apply( A, w, b);
            normb= ExX.Norm(w, true, true);                                 // This works only with useAcc==true
        }
        if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;

        resid = beta/normb;
        if (resid <= tol) {
            tol= resid;
            max_iter= 0;
            return true;
        }

        int j= 1;
        while (j <= max_iter) {
            v[0]= r*(1.0/beta);                             // v has accumulated form because r is accumulated
            s= 0.0;
            s[0]= beta;

            int i;
            for (i= 0; i<m-1 && j<=max_iter; ++i, ++j) {
                if (method == RightPreconditioning)
                {
                    M.Apply( A, w, v[i]);                   // hopefully, preconditioner do right things with accumulated v[i]
                    w=A*w;
                    ExX.Accumulate(w);
                }
                else
                    M.Apply( A, w, A*v[i]);
                for (int k= 0; k <= i; ++k ) {
                    H( k, i)= ExX.ParDot(w, true, v[k], true, true);        // This works only with useAcc==true
                    w-= H( k, i)*v[k];
                }

                if (i == m - 1) break;

                H( i + 1, i)= ExX.Norm(w, true, true);                      // This works only with useAcc==true
                v[i + 1]= w*(1.0/H( i + 1, i));

                for (int k= 0; k < i; ++k)
                    GMRES_ApplyPlaneRotation( H(k,i), H(k + 1, i), cs[k], sn[k]);

                GMRES_GeneratePlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( H(i,i), H(i+1,i), cs[i], sn[i]);
                GMRES_ApplyPlaneRotation( s[i], s[i+1], cs[i], sn[i]);

                resid= std::abs( s[i+1])/normb;

                if (resid <= tol) {
                    if (method == RightPreconditioning)
                    {
                        z=0.;
                        GMRES_Update( z, i, H, s, v);
                        M.Apply( A, t, z);
                        x_acc+=t;
                    }
                    else
                        GMRES_Update( x_acc, i, H, s, v);
                    tol= resid;
                    max_iter= j;
                    return true;
                }
            }

            if (method == RightPreconditioning)
            {
                z=0.;
                GMRES_Update( z, i-1, H, s, v);
                M.Apply( A, t, z);
                x_acc+=t;
                r= ExX.GetAccumulate(Vec(b - A*x_acc));
            }
            else
            {
                GMRES_Update( x_acc, i-1, H, s, v);
                M.Apply( A, r, Vec( b - A*x_acc));
            }
            beta=ExX.Norm(r, true, true);                       // this works only with useAcc==true
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
}


/// \brief Parallel GMRES-method with modifications for better scalability.
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParModGMRES(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                   int m, int& max_iter, double& tol,
                   bool measure_relative_tol, bool useAcc, bool useMGS, PreMethGMRES method)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in,out] M                    Preconditioner
    /// \param[in]     m                    number of steps after a restart is performed
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     useAcc               use accur variant for performing inner products and norms or do not use accure variant
    /// \param[in]     useMGS               use modified Gramm-Schmidt ortogonalization (many more sync-points exists!)
    /// \param[in]     method               left or right preconditioning (see solver.h for definition and declaration)
    /// \return  convergence within max_iter iterations
    /// \pre     the preconditioner should be able to handle a accumulated b
{
    Assert(x_acc.size()==b.size() && x_acc.size()==ExX.GetNum(), DROPSErrCL("ParModGMRES: Incompatible dimension"), DebugParallelNumC);

    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n = x_acc.size();      // dimension

    DMatrixCL<double> H(m,m);           // upper Hessenberg-matrix
    VectorCL tmpHCol(m);
    double beta, normb, resid;

    Vec r(n), w(n), w_acc(n), r_acc(n), z_acc(n), t_acc(n);
    Vec c(m), s(m), gamma(m);

    std::vector<Vec> v_acc(m);    // basis of the krylov-subspaces
    for (int i=0; i<m; ++i)
        v_acc[i].resize(n);

    if (method == RightPreconditioning){
        r    = b - A*x_acc;
        beta = ExX.Norm(r, false, useAcc, &r_acc);
        normb= ExX.Norm(b, false, useAcc);
    }
    else{
        M.Apply(A, r, VectorCL( b-A*x_acc));
        beta = ExX.Norm(r, M.RetAcc(), useAcc, &r_acc);
        M.Apply(A, w, b);
        normb = ExX.Norm(w, M.RetAcc(), useAcc, &w_acc);
    }

    if (normb == 0. || measure_relative_tol==false) normb=1.0;

    resid = beta/normb;
    if (resid<=tol){                        // finished
        tol = resid; max_iter = 0; return true;
    }

    int j=1;                                // number of steps
    while (j<=max_iter)
    {
        v_acc[0] = r_acc * (1./beta);
        gamma    = 0.;
        gamma[0] = beta;

        int i;
        for (i=0; i<m-1 && j<=max_iter; ++i, ++j)
        {
            if (method == RightPreconditioning){
                M.Apply(A, w_acc, v_acc[i]);                // hopefully M does the right thing
                w = A*w_acc;
                w_acc = ExX.GetAccumulate(w);
            }
            else{
                M.Apply(A, w, A*v_acc[i]);
                if (!M.RetAcc())
                    w_acc = ExX.GetAccumulate(w);
                else
                    w_acc = w;
            }

            // Gramm-Schmidt ortogonalization without  update of w!
            if (!useMGS)
                StandardGrammSchmidt(H, w_acc, true, v_acc, true, i, ExX, tmpHCol);
            else
                ModifiedGrammSchmidt(H, w_acc, true, v_acc, true, i, ExX);

            H(i+1,i) = ExX.Norm(w_acc, true, useAcc);
            v_acc[i+1] = w_acc * (1.0 / H(i+1,i));

            for (int k=0; k<i; ++k)
                GMRES_ApplyPlaneRotation(H(k,i),H(k+1,i), c[k], s[k]);

            GMRES_GeneratePlaneRotation(H(i,i), H(i+1,i), c[i], s[i]);
            GMRES_ApplyPlaneRotation(H(i,i), H(i+1,i), c[i], s[i]);
            GMRES_ApplyPlaneRotation(gamma[i], gamma[i+1], c[i], s[i]);

            resid = std::abs(gamma[i+1])/normb;

            if (resid<=tol){            // finished
                if (method == RightPreconditioning){
                    z_acc=0.;
                    GMRES_Update( z_acc, i, H, gamma, v_acc);
                    M.Apply( A, t_acc, z_acc);              // hopefully M does the right thing
                    x_acc+=t_acc;
                }
                else
                    GMRES_Update(x_acc, i, H, gamma, v_acc);
                tol = resid; max_iter = j; return true;
            }
        }
        if (method == RightPreconditioning){
            z_acc=0.;
            GMRES_Update( z_acc, i-1, H, gamma, v_acc);
            M.Apply( A, t_acc, z_acc);                      // hopefully M does the right thing
            x_acc += t_acc;
            r      = b-A*x_acc;
            beta = ExX.Norm(r, false, useAcc, &r_acc);
        }
        else{
            GMRES_Update(x_acc, i-1, H, gamma, v_acc);
            M.Apply(A, r, static_cast<Vec>( b-A*x_acc));
            beta = ExX.Norm(r, M.RetAcc(), ( M.RetAcc()?true:useAcc), &r_acc);
        }


        resid = beta/normb;
        if (resid<=tol){                // finished
            tol = resid; max_iter=j; return true;
        }
    }
    tol = std::fabs(gamma[0]);
    return false;
}


/// \brief Preconditioned BiCGStab-Method with accure inner products
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
  bool ParBiCGSTAB(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX,
                   PreCon& M, int& max_iter, double& tol, bool measure_relative_tol)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in]     M                    Preconditioner
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \pre                                the preconditioner must fulfill the following condition: M.Apply(A,x,b): b is accumulated => x is accumulated
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    const size_t n = x_acc.size();
    Vec r(b-A*x_acc), r_acc(n), r0hat_acc(n),
        p_acc(n), phat_acc(n),
        v(n), v_acc(n),
        s_acc(n), shat_acc(n),
        t(n), t_acc(n), that(n), that_acc(n);

    double resid, rho=1, alpha=1, beta, omega=1, sigma;

    VectorCL dots(4), glob_dots(4);

    double normb = std::sqrt(ExX.AccNorm_sq(b));
    if ( normb==0. || !measure_relative_tol) normb= 1.;

    sigma = ExX.AccNorm_sq(r, r_acc);
    resid = std::sqrt(sigma) / normb;

    M.Apply(A, r0hat_acc, r);
    if (!M.RetAcc())
        r0hat_acc= ExX.GetAccumulate(r0hat_acc);

    for (int i=0; i<max_iter; ++i)
    {
        beta = (sigma*alpha)/(rho*omega);
        rho  = sigma;
        if (i>0)
            p_acc    = r_acc + beta*p_acc - (beta*omega)*v_acc;
        else
            p_acc    = r_acc;

        M.Apply(A, phat_acc, p_acc);

        v= A*phat_acc;

        dots[0]= ExX.LocDot(    v, false, r0hat_acc, true, true, &v_acc);
        dots[1]= ExX.LocNorm_sq(r_acc, true, true);

        ProcCL::GlobalSum(Addr(dots), Addr(glob_dots), 2);

        resid = std::sqrt(glob_dots[1])/normb;
        if (resid<tol){
            tol= resid;
            max_iter=i;
            return true;
        }

        sigma = glob_dots[0];

        if (glob_dots[0]==0.)
        {
            if (ProcCL::MyRank()==0)
                std::cout << ">>>>> BREAKDOWN of BiCGStab!" <<std::endl;
            tol = resid;
            max_iter=i;
            return false;
        }

        alpha = rho/sigma;
        s_acc = r_acc -alpha*v_acc;

        M.Apply(A, shat_acc, s_acc);
        t = A*shat_acc;

        if (!M.RetAcc()){
            M.Apply(A,that,t);
            dots[0]= ExX.LocDot(that, false, shat_acc, true, true, &that_acc);
        }
        else{
            M.Apply(A,that_acc,t);
            dots[0]= ExX.LocAccDot(that_acc, shat_acc);
        }

        dots[1]= ExX.LocAccDot(that_acc, that_acc);
        dots[2]= ExX.LocAccDot(r0hat_acc, s_acc);
        dots[3]= ExX.LocDot(t, false, r0hat_acc, true, true, &t_acc);

        ProcCL::GlobalSum(Addr(dots), Addr(glob_dots), 4);

        if (glob_dots[1]==0.)
        {
            if (ProcCL::MyRank()==0)
                std::cout << ">>>>> BREAKDOWN of BiCGStab!" <<std::endl;
            tol = resid;
            max_iter=i;
            return false;
        }

        omega = glob_dots[0]/glob_dots[1];
        sigma = glob_dots[2] - omega * glob_dots[3];

        r_acc =  s_acc - omega * t_acc;
        x_acc += alpha * phat_acc + omega * shat_acc;
    }

    tol = resid;
    return false;
}


/// \brief Parallel preconditioned general conjugate residual algorithm with truncation
/** Only the last m (truncation parameter) residual vectors are kept. If m==max_iter,
    then all residual vectors are stored */
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParModPGCR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                int m, int& max_iter, double& tol,bool measure_relative_tol /*=true*/)
    /// \param[in]     A                    coefficients of the linear equation system
    /// \param[in,out] x_acc                IN: start vector, OUT: solution of the linear equation system
    /// \param[in]     b                    rhs of the linear equation system
    /// \param[in]     ExX                  class for exchanging values for accumulation
    /// \param[in]     M                    Preconditioner
    /// \param[in]     m                    truncation parameter
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    if (m>max_iter) m=max_iter;

    // preprocessing and malloc
    Vec r(b-A*x_acc), r_acc(r), y(b.size()), z_acc(b.size());       // a=A*r

    std::vector<Vec> p_acc, Ap, Ap_acc;
    VectorCL a_loc(m+1), a(m+1), c(m), gamma_loc(2), gamma(2);  // derived from std::valarray, so all doubles are stored consecutive in memory!

    double normb = ExX.Norm(b);
    if (normb==0.0 || measure_relative_tol==false) normb=1.0;

    double resid = std::sqrt(ExX.ParDotAcc(r_acc,r))/normb;
    if (resid<tol){
        max_iter=0; tol=resid; return true;
    }
    M.Apply(A,z_acc,r_acc);

    p_acc.push_back(z_acc);
    Ap.push_back(A*z_acc);
    Ap_acc.push_back(Ap[0]);
    int last_idx=0, next_idx=1;

    for (int j=0; j<max_iter; ++j)
    {
        // Calc of (r,Ap) / (Ap,Ap)
        gamma_loc[0] = ExX.DotAcc(Ap_acc[last_idx], Ap[last_idx]);  // (Ap,Ap) and accumulation of Ap
        gamma_loc[1] = dot(r_acc, Ap[last_idx]);                    // (r,Ap)
        DROPS::ProcCL::GlobalSum(Addr(gamma_loc), Addr(gamma), 2);
        const double alpha = gamma[1] / gamma[0];

        // update of x and r
        x_acc += alpha * p_acc[last_idx];      // axpy(alpha,  p_acc[last_idx] , x_acc);
        r     -= alpha * Ap[last_idx];         // axpy(-alpha, Ap[last_idx]    , r);
        r_acc -= alpha * Ap_acc[last_idx];     // axpy(-alpha, Ap_acc[last_idx], r_acc);

        // orthogonalization of p
        M.Apply(A,z_acc,r_acc);
        y  = A*z_acc;
        c[last_idx] = gamma[0];
        int k;
        for (k=0; k<=j && k<m; ++k)
            a_loc[k] = dot(y,Ap_acc[k]);

        a_loc[k] = dot(z_acc,r_acc);                // calc of the residual

        ProcCL::GlobalSum(Addr(a_loc), Addr(a), k+1);

        if (a[k]<0)
            std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of resid in PGCR because of accumulation!" << std::endl;

        resid = std::sqrt(a[k]) / normb;
        if (true && j%20==0){
            double realresid = ExX.Norm(static_cast<Vec>(b-A*x_acc));
            if (ProcCL::IamMaster())
                std::cout << "j: "<<j<<"\tresid="<< resid<<", realresid " <<realresid << std::endl;
        }

        if (/*a[j+1]*/ resid<tol){
            max_iter=j; tol=resid; return true;
        }

        for (int i=0; i<=j && i<m; ++i)
            a[i] = a[i]/c[i];

        if (j+1<m){
            p_acc.push_back(z_acc);
            Ap.push_back(y);
        }
        else{
            p_acc[next_idx] = z_acc;
            Ap[next_idx]    = y;
        }

        for (int i=0; i<=j && i<m; ++i){
            if (i!=next_idx){
                p_acc[next_idx] -= a[i] * p_acc[i];   // axpy(-a[i],p_acc[i],p_acc[next_idx]);
                Ap[next_idx]    -= a[i] * Ap[i];      // axpy(-a[i],Ap[i],Ap[next_idx]);
            }
        }
        if (j+1<m)
            Ap_acc.push_back(Ap[j+1]);
        else
            Ap_acc[next_idx] = Ap[next_idx];

        last_idx = next_idx;
        next_idx = (next_idx+1) % m;
    }
    tol = resid;
    return false;
}

/// \brief Parallel preconditioned general conjugate residual algorithm with truncation and high accuracy
/** Only the last m (truncation parameter) residual vectors are kept. If m==max_iter,
    then all residual vectors are stored */
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParModAccurPGCR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
                     int m, int& max_iter, double& tol, bool measure_relative_tol /*=true*/,
                     std::ostream* os)
    /// \param[in]     A                    coefficients of the linear equation system
    /// \param[in,out] x_acc                IN: start vector, OUT: solution of the linear equation system
    /// \param[in]     b                    rhs of the linear equation system
    /// \param[in]     ExX                  class for exchanging values for accumulation
    /// \param[in]     M                    Preconditioner
    /// \param[in]     m                    truncation parameter
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    if (m>max_iter) m=max_iter;

    // preprocessing and malloc
    Vec r(b-A*x_acc), r_acc(r), y(b.size()), y_acc(b.size()), z_acc(b.size());

    std::vector<Vec> p_acc, Ap, Ap_acc;
    VectorCL a_loc(m+1), a(m+1), c(m), gamma_loc(2), gamma(2);  // derived from std::valarray, so all doubles are stored consecutive in memory!
    VectorCL b_acc(b.size());

    double normb = ExX.Norm(b, false, true, &b_acc);
    if (normb==0.0 || measure_relative_tol==false)
        normb=1.0;

    double resid= ExX.Norm(r, false, true, &r_acc)/normb;
    if (resid<tol){
        max_iter=0; tol=resid; return true;
    }
    M.Apply(A, z_acc, r_acc);

    p_acc.push_back(z_acc);
    Ap.push_back(A*z_acc);
    Ap_acc.push_back(Ap[0]);
    int last_idx=0, next_idx=1;
    double a_min=0;

    for (int j=0; j<max_iter; ++j)
    {
        // Calc of (r,Ap) / (Ap,Ap)
        gamma_loc[0]= ExX.LocNorm_sq( Ap[last_idx], false, true, &Ap_acc[last_idx]);
        gamma_loc[1]= ExX.LocDot( r_acc, true, Ap_acc[last_idx], true, true);
        ProcCL::GlobalSum(Addr(gamma_loc), Addr(gamma), 2);
        const double alpha = gamma[1] / gamma[0];

        // update of x and r
        x_acc += alpha * p_acc[last_idx];      // axpy(alpha,  p_acc[last_idx] , x_acc);
        r     -= alpha * Ap[last_idx];         // axpy(-alpha, Ap[last_idx]    , r);
        r_acc -= alpha * Ap_acc[last_idx];     // axpy(-alpha, Ap_acc[last_idx], r_acc);

        // orthogonalization of p
        M.Apply(A,z_acc,r_acc);
        y  = A*z_acc;
        y_acc = ExX.GetAccumulate(y);
        c[last_idx] = gamma[0];
        int k;
        for (k=0; k<=j && k<m; ++k)
            a_loc[k]= ExX.LocDot( y_acc, true, Ap_acc[k], true, true);
        a_loc[k]= ExX.LocDot( z_acc, true, r_acc, true, true);  // calc of the residual

        ProcCL::GlobalSum(Addr(a_loc), Addr(a), k+1);

        resid = std::sqrt(a[k]) / normb;
        if (os)
            IF_MASTER
                (*os) << "GCR " << j << ": res " << resid << std::endl;

        if (resid<tol){
            max_iter=j; tol=resid; return true;
        }

        for (int i=0; i<=j && i<m; ++i)
            a[i] = a[i]/c[i];

        if (j<m){
            p_acc.push_back(z_acc);
            Ap.push_back(y);
        }
        else{
            p_acc[next_idx] = z_acc;
            Ap[next_idx]    = y;
        }

        for (int i=0; i<=j && i<m; ++i){
            if (i!=next_idx){
                p_acc[next_idx] -= a[i] * p_acc[i];   // axpy(-a[i],p_acc[i],p_acc[next_idx]);
                Ap[next_idx]    -= a[i] * Ap[i];      // axpy(-a[i],Ap[i],Ap[next_idx]);
            }
        }
        if (j<m)
            Ap_acc.push_back(Ap[j+1]);
        else
            Ap_acc[next_idx] = Ap[next_idx];

        last_idx = next_idx;
        if (j<m)
            next_idx = (next_idx+1) % m;
        else{
            next_idx = (next_idx+1) % m;
            int min_idx = 0;
            a_min= std::fabs( a[0]); // m >= 1, thus this access is valid.
            for (int i= 1; i < k && i < m; ++i){
                if ( std::fabs( a[i]) < a_min) {
                    min_idx= i;
                    a_min= std::fabs( a[i]);
                }
            }
            next_idx=min_idx;
        }
    }
    tol = resid;
    return false;
}

/// \brief Parallel preconditioned general conjugate residual algorithm with truncation and high accuracy
///        with same strategy as serial version
template <typename Mat, typename Vec, typename PreCon, typename ExCL>
bool ParPGCR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, PreCon& M,
             int m, int& max_iter, double& tol, bool measure_relative_tol, bool useAcc)
    /// \param[in]     A                    coefficients of the linear equation system
    /// \param[in,out] x_acc                IN: start vector, OUT: solution of the linear equation system
    /// \param[in]     b                    rhs of the linear equation system
    /// \param[in]     ExX                  class for exchanging values for accumulation
    /// \param[in]     M                    Preconditioner
    /// \param[in]     m                    truncation parameter
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \param[in]     useAcc               use accurate or fast variant for inner products and norms
    /// \return                             convergence within max_iter iterations
{
    // Check if preconditioner needs diagonal of matrix. The preconditioner
    // only computes the diagonal new, if the matrix has changed
    if (M.NeedDiag())
        M.SetDiag(A);

    m= (m <= max_iter) ? m : max_iter;

    Vec r( b - A*x_acc);
    Vec sn( b.size()), vn( b.size()), sn_acc(b.size()), vn_acc(b.size());
    std::vector<Vec> s, v;
    std::vector<double> a( m);

    // Norm Berechnung der rechten Seite detailliert:
//         VectorCL b0( b[std::slice( 0, A.num_rows( 0), 1)]);
//         VectorCL b1( b[std::slice( A.num_rows( 0), A.num_rows( 1), 1)]);
//         double normb0= ExX.Get(0).Norm( b0, false, useAcc);
//         double normb1= ExX.Get(1).Norm( b1, false, useAcc);
//         double max_rhs=GlobalMax(supnorm(b));
    // Norm Berechnung der rechten Seite detailliert:

    double normb= ExX.Norm( b, false, useAcc);
//         const double mynorm=normb;
    if (normb == 0.0 || measure_relative_tol == false) normb= 1.0;
    double resid= ExX.Norm( r, false, useAcc)/normb;

//     IF_MASTER
//       std::cout << "Starting GCR with tol: "<<tol<<",\tnorm_rhs: "<<mynorm<<'\n'
//                 << "          norm_rhs_up: "<<normb0<<"\tnorm_rhs_low:   "<<normb1<<'\n'
//                 << "          max_rhs:     "<<max_rhs<<'\n';
    for (int k= 0; k < max_iter; ++k) {
        // compute residual for output
//             VectorCL r0( r[std::slice( 0, A.num_rows( 0), 1)]);
//             VectorCL r1( r[std::slice( A.num_rows( 0), A.num_rows( 1), 1)]);
//             double res_0=ExX.Get(0).Norm(r0,false, true);
//             double res_1=ExX.Get(1).Norm(r1,false, true);
//             double max_x  = GlobalMax(supnorm(x_acc));
//             double max_up = GlobalMax(supnorm(Vec(x_acc[std::slice( 0, A.num_rows( 0), 1)])));
//             double max_low= GlobalMax(supnorm(Vec(x_acc[std::slice( A.num_rows( 0), A.num_rows( 1), 1)])));
//             IF_MASTER
//                 std::cout << "GCR: k: " << k << "\tresidual: " << resid<<'\n'
//                         << "         \tresid_up:   "<<res_0<<'\n'
//                         << "         \tresid_low:  "<<res_1<<'\n'
//                         << "         \tsup(x):     "<<max_x<<'\n'
//                         << "         \tsup(x_up):  "<<max_up<<'\n'
//                         << "         \tsup(x_low): "<<max_low<<'\n';
        // end of residual computation for output
        if (resid < tol) {
            tol= resid;
            max_iter= k;
            return true;
        }
        M.Apply( A, sn, r);

        if (!M.RetAcc())
            ExX.Accumulate(sn);
        // sn is accumulated
        vn= A*sn;

        for (int i= 0; i < k && i < m; ++i) {
            const double alpha= ExX.ParDot( vn, false, v[i], false, useAcc);
            a[i]= alpha;
            vn-= alpha*v[i];
            sn-= alpha*s[i];
        }
        const double beta= ExX.Norm( vn, false, useAcc);
        vn/= beta;
        sn/= beta;
        const double gamma= ExX.ParDot( r, false, vn, false, useAcc);
        x_acc+= gamma*sn;

        r-= gamma*vn;
        resid= ExX.Norm( r, false, useAcc)/normb;
        if (k < m) {
            s.push_back( sn);
            v.push_back( vn);
        }
        else {
            int min_idx= 0;
            double a_min= std::fabs( a[0]); // m >= 1, thus this access is valid.
            for (int i= 1; i < k && i < m; ++i)
                if ( std::fabs( a[i]) < a_min) {
                    min_idx= i;
                    a_min= std::fabs( a[i]);
                }
                s[min_idx]= sn;
                v[min_idx]= vn;
        }
    }
    tol= resid;
    return false;
}

/// \brief Preconditioned QMR-Method
template <typename Mat, typename Vec, typename Lanczos, typename ExCL>
bool ParQMR(const Mat& A, Vec& x_acc, const Vec& b, const ExCL& ExX, Lanczos lan,
            int& max_iter, double& tol, bool measure_relative_tol)
    /// \param[in]     A                    local distributed coefficients-matrix of the linear equation system
    /// \param[in,out] x_acc                start vector and the solution in accumulated form
    /// \param[in]     b                    rhs of the linear equation system (distributed form)
    /// \param[in]     ExX                  ExchangeCL corresponding to the RowIdx of x and the ColIdx of A
    /// \param[in]     lan                  Lanczos-Class for computing a bi-orthogonal basis of Krylov subspaces
    /// \param[in,out] max_iter             IN: maximal iterations, OUT: used iterations
    /// \param[in,out] tol                  IN: tolerance for the residual, OUT: residual
    /// \param[in]     measure_relative_tol if true stop if |M^(-1)(b-Ax)|/|M^(-1)b| <= tol, else stop if |M^(-1)(b-Ax)|<=tol
    /// \return                             convergence within max_iter iterations
{
    tol*=tol;

    Vec d_acc(x_acc.size()), s(d_acc), r(b-A*x_acc), r_acc(r);

    double normb = ExX.Norm_sq(b);
    if (normb==0.0 || measure_relative_tol==false) normb=1.0;
    double norm_r = std::sqrt(std::fabs(ExX.ParDotAcc(r_acc,r)));

    if (norm_r/normb<std::sqrt(tol))
    {
        tol = std::fabs(norm_r)/std::sqrt(normb);
        max_iter=1;
        return true;
    }

    lan.Init(A,static_cast<Vec>(r*(1./norm_r)),static_cast<Vec>(r*(1./norm_r)),ExX,1);
    double tetha,
    kappa=-1,
    lambda=1,
    beta,
    rho = norm_r,
    rho_last = norm_r;
//  rho_last=lan.GetRho();
    double quot;

    for (int j=1; j<=max_iter; ++j)
    {
        if (!lan.Step()){
            tol = std::sqrt(norm_r/normb);
            max_iter=j-1;
            return false;
        }

        beta    = lan.GetBeta();
        rho_last= rho;
        rho     = lan.GetRho();
        quot    = (lambda*beta*beta+rho*rho);
        tetha   = beta*beta * (1-lambda)/quot;
        kappa   = -rho_last*beta*kappa/quot;
        lambda  = lambda*beta*beta/quot;

        d_acc   = tetha * d_acc + kappa*lan.GetAccP();
        s       = tetha * s     + kappa*lan.GetAp();
        x_acc  += d_acc;
        r      -= s;

        r_acc   = r;
        norm_r  = ExX.ParDotAcc(r_acc,r);

        if (norm_r<0)
            std::cout << "["<<ProcCL::MyRank()<<"]==> negative squared norm of resid in QMR because of accumulation!" << std::endl;
        if (norm_r/normb<tol)
        {
            tol = std::sqrt(std::fabs(norm_r)/normb);
            max_iter=j;
            return true;
        }
    }
    tol = std::sqrt(norm_r);
    return false;
}

} // end of namespace DROPS

#endif
