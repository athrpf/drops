//**************************************************************************
// File:    parstokessolver.h                                              *
// Content: parallel solvers for different problems                        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   Februar, 8th 2006                                              *
//**************************************************************************
/// \author Oliver Fortmeier, SC RWTH Aachen
/// \file parstokessolver.h
/// \brief Parallel itertative solvers for stokes problems

#ifndef DROPS_PARSTOKESSOLVER_H_
#define DROPS_PARSTOKESSOLVER_H_

#include "num/parsolver.h"
#include "num/stokessolver.h"
#include "parallel/logger.h"
#include <algorithm>

namespace DROPS
{

//=============================================================================
//  The Stokes solvers solve systems of the form
//    A v + BT p = b
//    B v        = c
//=============================================================================

// Inexact Uzawa algorithm
template <typename Mat, typename Vec, typename PC1, typename PC2, typename ExVCL, typename ExPCL>
bool ParInexactUzawa(const Mat& A, const Mat& B, Vec& xu_acc, Vec& xp_acc, const Vec& f, const Vec& g,
                  const ExVCL& exV, const ExPCL& exP, PC1& Apc, PC2& Spc,
                  int& max_iter, double& tol, int &usedinnerit, InexactUzawaApcMethodT apcmeth= APC_OTHER,
                  double innerred= 0.3, int innermaxiter= 500, bool useAcc=true, std::ostream* output=0);


// ***************************************************************************
/// \brief Parallel base solver class for Stokes problems
// ***************************************************************************
template <typename ApcT, typename SpcT, typename ExVCL, typename ExPCL>
class ParStokesSolverBaseCL : public StokesSolverBaseCL
{
  protected:
    typedef StokesSolverBaseCL base;// base class

    ExVCL* exV_;                    // exchange velocity
    ExPCL* exP_;                    // exchange pressure
    ApcT*  Apc_;                    // preconditioner for A
    SpcT*  Spc_;                    // preconditioner for Schurcomplement

  public:
    ParStokesSolverBaseCL(int maxiter, double tol, ExVCL& exv, ExPCL& exp, ApcT& apc, SpcT& spc)
      : base(maxiter, tol), exV_(&exv), exP_(&exp), Apc_(&apc), Spc_(&spc)
        /// \param maxiter maximal iterations
        /// \param tol     tolerance for residual
        /// \param exv     exchange class for velocity
        /// \param exp     exchange class for pressure
        /// \param apc     preconditioner for A-block
        /// \param spc     preconditioner for Schur complement matrix
    {}

    ExVCL& GetExV()                { return *exV_; }    ///< return reference on exchange class for velocities
    const ExVCL& GetExV() const    { return *exV_; }    ///< return constant reference on exchange class for velocities

    ExPCL& GetExP()                { return *exP_; }    ///< return reference on exchange class for pressures
    const ExPCL& GetExP() const    { return *exP_; }    ///< return constant reference on exchange class for pressures

    ApcT& GetInnerPC()             { return *Apc_; }    ///< return reference on inner preconditioner
    const ApcT& GetInnerPC() const { return *Apc_; }    ///< return constant reference on inner preconditioner

    SpcT& GetOuterPC()             { return *Spc_; }    ///< return reference on inner preconditioner
    const ApcT& GetOuterPC() const { return *Spc_; }    ///< return constant reference on inner preconditioner

};

// ***************************************************************************
/// \brief Parallel preconditioned inexact Uzawa class
// ***************************************************************************
template <typename ApcT, typename SpcT, InexactUzawaApcMethodT ApcMeth= APC_OTHER,
          typename ExVCL=ExchangeCL, typename ExPCL=ExchangeCL>
class ParInexactUzawaCL : public ParStokesSolverBaseCL<ApcT, SpcT, ExVCL, ExPCL>
{
  private:
    typedef ParStokesSolverBaseCL<ApcT, SpcT, ExVCL, ExPCL> base;

    double innerreduction_;     // reduction
    int    innermaxiter_;       // maximal inner iterations
    int    inneriter_;          // number of inneritter (accumulated)
    std::ostream* output_;      // for output

  public:
      ParInexactUzawaCL(ApcT& Apc, SpcT& Spc, ExVCL& exv, ExPCL& exp, int outer_iter, double outer_tol,
                        double innerreduction= 0.3, int innermaxiter= 500, std::ostream* output= 0)
      : base( outer_iter, outer_tol, exv, exp, Apc, Spc),
        innerreduction_( innerreduction), innermaxiter_( innermaxiter), inneriter_(0), output_(output)
    {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c)
    /// Solve the Stokes problem given by the matrices \a A and \a B and rhs \a b and \a c with the inexact Uzawa
    /// algorithm and stores the solution for the velocity and pressure in \a v and \a p
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        inneriter_ = 0;
        ParInexactUzawa( A, B, v, p, b, c, *base::exV_, *base::exP_, *base::Apc_,
                      *base::Spc_, base::_iter, base::_res, inneriter_, ApcMeth, innerreduction_, innermaxiter_, true, output_);
    }

    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                        const VectorCL& b, const VectorCL& c)
    {
        base::_res=  base::_tol;
        base::_iter= base::_maxiter;
        inneriter_ = 0;
        ParInexactUzawa( A, B, v, p, b, c, *base::exV_, *base::exP_, *base::Apc_,
                      *base::Spc_, base::_iter, base::_res, inneriter_, ApcMeth, innerreduction_, innermaxiter_, true, output_);
    }

    int GetInnerIter()    const  { return inneriter_; }
    int GetInnerMaxIter() const  { return innermaxiter_; }
};


// ***************************************************************************
/// \brief Parallel Schur complement matrix class
// ***************************************************************************
template<typename, typename, typename>
  class ParSchurComplMatrixCL;

template<typename T1, typename T2, typename T3>
  VectorCL operator*(const ParSchurComplMatrixCL<T1, T2, T3>&, const VectorCL&);

template<typename PoissonSolverT, typename MatT, typename ExVCL>
class ParSchurComplMatrixCL
/// This class can multiply B*A^(-1)*B * x.
{
  private:
    PoissonSolverT& solver_;
    const MatT& A_;
    const MatT& B_;
    const ExVCL& ex_;

  public:
    ParSchurComplMatrixCL(PoissonSolverT& solver, const MatT& A, const MatT& B, const ExVCL& ex)
      : solver_(solver), A_( A), B_( B), ex_(ex) {}

    friend VectorCL
    operator*<>(const ParSchurComplMatrixCL<PoissonSolverT, MatT, ExVCL>&, const VectorCL&);
};

template<class PoissonSolverT, typename MatT, typename ExVCL>
  VectorCL operator*(const ParSchurComplMatrixCL<PoissonSolverT, MatT, ExVCL>& M, const VectorCL& v)
{
    VectorCL x(  M.A_.num_cols() );
    M.solver_.Solve( M.A_, x, transp_mul( M.B_, v));
//     IF_MASTER
//       std::cerr << "> inner iterations: " << M.solver_.GetIter()
//                 << "\tresidual: " << M.solver_.GetResid()
//                 << std::endl;
    return M.B_*x;
}


// ***************************************************************************
/// \brief Parallel approximate Schur complement matrix class
// ***************************************************************************
template<typename, typename, typename>
  class ParApproximateSchurComplMatrixCL;

template<typename T1, typename T2, typename T3>
  VectorCL operator*(const ParApproximateSchurComplMatrixCL<T1, T2, T3>&, const VectorCL&);

template<typename APC, typename MatT, typename ExVCL>
class ParApproximateSchurComplMatrixCL
/// Assume M is a approximate of A. Then instead of performing B*A^(-1)*B^T * v this class can
/// perform the multiplication B*M^(-1)*B^T * v
{
  private:
    const MatT& A_;
    APC& Apc_;
    const MatT& B_;
    const ExVCL& ex_;

  public:
    ParApproximateSchurComplMatrixCL(const MatT& A, APC& Apc, const MatT& B, const ExVCL& ex)
      : A_( A), Apc_( Apc), B_( B), ex_(ex) {}

    friend VectorCL
    operator*<>(const ParApproximateSchurComplMatrixCL<APC, MatT, ExVCL>&, const VectorCL&);
};

template<class APC, typename MatT, typename ExVCL>
  VectorCL operator*(const ParApproximateSchurComplMatrixCL<APC, MatT, ExVCL>& M, const VectorCL& v)
{
    VectorCL x( 0.0, M.B_.num_cols());
    VectorCL r= transp_mul( M.B_, v);
    M.Apc_.Apply( M.A_, x, r);
    if (!M.Apc_.RetAcc())
        M.ex_.Accumulate(x);
    return M.B_*x;
}


//***************************************************************************
// Implementations of the methods
//***************************************************************************

template <typename Mat, typename Vec, typename PC1, typename PC2, typename ExVCL, typename ExPCL>
  bool ParInexactUzawa(const Mat& A, const Mat& B, Vec& xu_acc, Vec& xp_acc, const Vec& f, const Vec& g,
                    const ExVCL& exV, const ExPCL& exP, PC1& Apc, PC2& Spc,
                    int& max_iter, double& tol, int &usedinnerit, InexactUzawaApcMethodT apcmeth,
                    double innerred, int innermaxiter, bool useAcc, std::ostream* output)
/// \param[in]     A            Coefficient matrix for velocities
/// \param[in]     B            Coefficient matrix for coupling velocity and pressure
/// \param[in,out] xu_acc IN:   Startvector for velocity, OUT: velocity result
/// \param[in,out] xp_acc IN:   Startvector for pressure, OUT: pressure result
/// \param[in]     f            upper rhs
/// \param[in]     g            lower rhs
/// \param[in]     exV          Exchange class for velocities
/// \param[in]     exP          Exchange class for pressures
/// \param[in]     Apc          Preconditioner for A
/// \param[in]     Spc          Preconditioner for Schur-Complement
/// \param[in,out] max_iter     IN: maximal iterations, OUT: used iterations
/// \param[in,out] tol          IN: tolerance for the residual, OUT: residual
/// \param[in,out] usedinnerit  number of accumulated steps for inner solver
/// \param[in]     apcmeth      Characteristics of the preconditioner for the A-block
/// \param[in]     innerred     inner reduction
/// \param[in]     innermaxiter maximal inner solver steps
/// \param[in]     useAcc       use accur variant for performing inner products and norms or do not use accure variant
/// \param[in,out] output       Debug information
{
    if (Apc.NeedDiag())
        Apc.SetDiag(A);
    if (Spc.NeedDiag())
        Spc.SetDiag(B);
    if (apcmeth==APC_SYM_LINEAR)
        throw DROPSErrCL("InexactUzawa: Because the exact UzawaAlgorithm has not been implemented the choice \"apcmeth==APC_SYM_LINEAR\" is not valid");
    VectorCL ru( f - A*xu_acc - transp_mul( B, xp_acc));
    VectorCL rp( g - B*xu_acc);
    VectorCL w( f.size());
    VectorCL z( g.size());
    VectorCL z2( g.size());
    VectorCL zbar( f.size());
    VectorCL zhat( f.size());
    VectorCL du( f.size());
    VectorCL c( g.size());
    ParApproximateSchurComplMatrixCL<PC1, Mat, ExVCL>* asc = apcmeth == APC_SYM_LINEAR ? 0 :
            new ParApproximateSchurComplMatrixCL<PC1, Mat, ExVCL>( A, Apc, B, exV);
    double innertol;
    int inneriter;
    int pr_iter_cumulative= 0;
    double res_u = exV.Norm_sq( ru, false, useAcc);
    double res_p = exP.Norm_sq( rp, false, useAcc);
    double resid0= std::sqrt( res_u + res_p);
    double resid= resid0;
    if (output)
      IF_MASTER
        (*output) << "   o InexactUzawa 0: res " << resid
                  << ", res-impuls " << std::sqrt(res_u)
                  << ", res-mass " << std::sqrt(res_p)
                  << '\n';
    if (resid <= tol) { // The fixed point iteration between levelset and Stokes
        tol= resid;     // equation uses this to determine convergence.
        max_iter= 0;
        if (asc!=0) delete asc;
        return true;
    }
    for (int k= 1; k <= max_iter; ++k) {
        w= 0.0;
        Apc.Apply( A, w, ru);
        if (!Apc.RetAcc())
            exV.Accumulate(w);
        c= B*w - rp;            // w is accumulated
        z= 0.0;
        z2= 0.0;
        inneriter= innermaxiter;
        switch (apcmeth) {
//             case APC_SYM_LINEAR:         // this case has not been implemented yet!
//                 zbar= 0.0;
//                 zhat= 0.0;
//                 innertol= innerred*norm( c);
//                 UzawaPCG( Apc, A, B, z, zbar, zhat, c, Spc, inneriter, innertol);
//                 break;
            case APC_SYM:
                innertol= innerred*exP.Norm(c, false, useAcc);
                if (useAcc)
                    ParAccurPCG( *asc, z, c, exP, Spc, inneriter, innertol);
                else
                    ParPCG( *asc, z, c, exP, Spc, inneriter, innertol);
                break;
            default:
                std::cerr << "WARNING: InexactUzawa: Unknown apcmeth; using GMRes.\n";
            // fall through
            case APC_OTHER:
                innertol= innerred; // GMRES can do relative tolerances.
                ParModGMRES( *asc, z, c, exP, Spc, /*restart*/ inneriter, inneriter, innertol, /*relative errors*/ true, useAcc);
                break;
        }
        if (apcmeth != APC_SYM_LINEAR) {
            zbar= transp_mul( B, z);        // z is accumulated
            zhat= 0.0;
            Apc.Apply( A, zhat, zbar);
            if (!Apc.RetAcc())
                exV.Accumulate(zhat);
        }
        pr_iter_cumulative+= inneriter;
        du= w - zhat;           // w and zhat are accumulated, so du is accumulated as well
        xu_acc+= du;
        xp_acc+= z;             // z is in every case accumulated, because it is a result of a solver
        ru-= A*du + zbar;
        rp= g - B*xu_acc;
        res_u= exV.Norm_sq( ru, false, useAcc);
        res_p = exP.Norm_sq( rp, false, useAcc);
        resid= std::sqrt( res_u + res_p);
        if (output)
          IF_MASTER
            (*output) << "   o InexactUzawa "<<k<<": res " << resid
                      << ", res-impuls " << std::sqrt(res_u)
                      << ", res-mass " << std::sqrt(res_p)
                      << ", inner_iter "<<inneriter
                      << '\n';

        IF_MASTER{
            DROPS_LOGGER_ADDVALUE("NavStokesInnerIter",inneriter);
        }

        if (resid <= tol) { // absolute errors
            tol= resid;
            max_iter= k;
            usedinnerit = pr_iter_cumulative;
            if (asc!=0) delete asc;
            return true;
        }
    }
    tol= resid;
    if (asc!=0) delete asc;
    usedinnerit = pr_iter_cumulative;
    IF_MASTER
        std::cerr << "===> Warning, InexactUzawa stoped without reaching tolerance!" << std::endl;
    return false;
}
}
#endif
