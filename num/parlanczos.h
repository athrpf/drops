//**************************************************************************
// File:    parlanczos.h                                                   *
// Content: lanczos algorithms                                             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   May, 8th 2006                                                  *
//**************************************************************************
/// \author Oliver Fortmeier, SC RWTH Aachen
/// \file   parlanczos.h
/// \brief  Parallel lanczos algorithms for QMR itartive solver


#ifndef _DROPS_PARLANCZOS_H_
#define _DROPS_PARLANCZOS_H_

#include "num/spmat.h"
#include "parallel/parallel.h"
//#include "parallel/partime.h"
#include "parallel/exchange.h"
#include <algorithm>

namespace DROPS
{

/****************************************************************************
* L A N C Z O S  3  C L A S S                                                   *
****************************************************************************/
/// \brief Encapsulate three-term Lanczos Algorithm
/** This class construct  biorthogonal basises to K(A,v0) and K(A^T,w0),
    by a modified Lanczos-Algorithm, that uses only one synchronisation-
    point within one step. This class uses a three-term recurrence.
*/
/****************************************************************************
* L A N C Z O S  3  C L A S S                                               *
****************************************************************************/
template <typename Mat, typename Vec>
class ParLanczos3CL
{
  private:
    const Mat           *A_;            // Find biorthogonal basis to this Matrix
    const ExchangeCL    *ex_;           // how to exchange numerical data

    SBufferCL<Vec,2> v_acc_;            // last two basis vectors of Krylov-Subspace K(A,v0)
    SBufferCL<Vec,2> w_;                // last two basis vectors of Krylov-Subspace K(A^T,v0)
    Vec u_acc_, w_acc_;                 // u = A*\hat{v}_j  (help-vector for put together sync-point), w_acc, accumulated form of w_j

    double alpha_, alpha_next_,         // alpha=(Av,w), alpha_next_ is used for put together sync-points
           beta_, delta_;               // beta=scaling of w, delta=scaling of v
    double breakdownEps_;               // tolerance for test to zero sqrt|(v,w)|
    double val_loc_[2], val_glob_[2];   // tmp filds for global-reduce operations
    bool init_;                         // flag, if the fist step is been made
    bool breakdown_;                    // flag, if Lanczos breaks in last step down

  public:
    ParLanczos3CL(const Mat& A, const Vec &v0_acc, const Vec &w0, const ExchangeCL& ex, const double BreakDownEps=DoubleEpsC);
    ParLanczos3CL(const double BreakDownEps=DoubleEpsC);

    bool Step();                                ///< Do a Lanczos-Step
    void Init(const Mat& A, const Vec &v0_acc, const Vec &w0, const ExchangeCL& ex);

    void   SetBreakDownEps(const double eps);   ///< Set new tolerance for breakdown
    double GetBreakDownEps();                   ///< Get tolerance for breakdown

    Vec& GetAccV()  {return v_acc_[1];}         ///< Get accumulated version of basisvector of Krylov-subspace K(A,v0)
    Vec& GetW()     {return w_[1];}             ///< Get distributed version of basisvector of Krylov-subspace K(A^T,w0)
    double GetAlpha() const {return alpha_;}    ///< Get last coefficient (Av,w)
    double GetBeta()  const {return beta_;}     ///< Get last scaling-parameter for w
    double GetDelta() const {return delta_;}    ///< Get last scaling-parameter for v
    bool Breakdown() const {return breakdown_;} ///< Check if Lanczos breaks in last step down
};


/****************************************************************************
* L A N C Z O S  2  C L A S S                                               *
****************************************************************************/
/// \brief Encapsulate coupled two-term Lanczos Algorithm
/** This class construct biorthogonal basises to K(A,v0) and K(A^T,w0),
    by a modified Lanczos-Algorithm, that uses only one synchronisation-
    point within one step. This class uses a coupled two-term recurrence.
 */
/****************************************************************************
* L A N C Z O S  2  C L A S S                                               *
****************************************************************************/
template <typename Mat, typename Vec, typename ExCL>
class ParLanczos2CL
{
  private:
    const Mat           *A_;            // Find biorthogonal basis to this Matrix
    const ExCL          *ex_;           // how to exchange numerical data

    Vec v_, v_acc_, Av_;
    Vec w_, w_acc_;                     // last two basis vectors of Krylov-Subspace K(A^T,v0)
    Vec p_acc_, Ap_;
    Vec q_;

    double eps_, delta_, rho_, xi_, beta_, sigma_;
    double breakdownEps_;

    double val_loc_[4], val_glob_[4];   // tmp filds for global-reduce operations
    bool init_;                         // flag, if the fist step is been made (for computation of p,q,beta)
    bool breakdown_;                    // flag, if Lanczos breaks in last step down
    bool luckybreakdown_;               // this indicates, if the breakdown is lucky (i.e. v==0 od w==0)

  public:
    ParLanczos2CL(const Mat&, const Vec &v0, const Vec &w0, const ExCL&, const double eps, const double BreakDownEps=DoubleEpsC);
    ParLanczos2CL(const double BreakDownEps=DoubleEpsC);

    bool Step();                                            ///< Do a Lanczos-Step
    void Init(const Mat& A, const Vec &v0, const Vec &w0, const ExCL& ex, const double eps);

    void   SetBreakDownEps(const double eps);               ///< Set new tolerance for breakdown
    double GetBreakDownEps();                               ///< Get tolerance for breakdown

    Vec& GetV()    /*const*/ {return v_;}                   ///< Get distributed version of basisvector of Krylov-subspace K(A,v0)
    Vec& GetAccV() /*const*/ {return v_acc_;}               ///< Get accumulated version of basisvector of Krylov-subspace K(A,v0)
    Vec& GetW()    /*const*/ {return w_;}                   ///< Get distributed version of basisvector of Krylov-subspace K(A^T,w0)
    Vec& GetAccW() /*const*/ {return w_acc_;}               ///< Get accumulated version of basisvector of Krylov-subspace K(A^T,w0)
    Vec& GetAp()   /*const*/ {return Ap_;}                  ///< Get distributed version of A times p
    Vec& GetAccP()           {return p_acc_;}

    double GetDelta() const {return delta_;}                ///< Get last (hat{v},hat{w})
    double GetEps()   const {return eps_;}                  ///< Get last scaling-parameter of p and q
    double GetBeta()  const {return beta_;}                 ///< Get last recurrence paramter for v and w
    double GetSigma() const {return sigma_;}                ///< Get last (A*hat{v},hat{w})
    double GetRho()   const {return rho_;}                  ///< Get scaling of v   (v^T*v)
    double GetXi()    const {return xi_;}                   ///< Get scaling of w   (w^T*w)
    bool Breakdown() const {return breakdown_;}             ///< Check if Lanczos breaks in last step down
    bool LuckyBreakdown() const {return luckybreakdown_;}   ///< Check if in last step v or w is zero
};

/****************************************************************************
* P R E C O N D I T I O N E D   L A N C Z O S  2 C L A S S                  *
****************************************************************************/
/// \brief Encapsulate coupled two-term Lanczos Algorithm
/** This class construct biorthogonal basises to K(A,v0) and K(A^T,w0),
    by a modified Lanczos-Algorithm, that uses only one synchronisation-
    point within one step. This class uses a coupled two-term recurrence.
 */
/****************************************************************************
* P R E C O N D I T I O N E D   L A N C Z O S  2 C L A S S                  *
****************************************************************************/
template <typename Mat, typename Vec, typename PC, typename ExCL>
class ParPreLanczos2CL
{
  private:
    const Mat           *A_;            // Find biorthogonal basis to this Matrix
    const ExCL          *ex_;           // how to exchange numerical data
    PC                  pc_;            // Preconditioner

    Vec v_, v_acc_, Av_;
    Vec w_, w_acc_;                     // last two basis vectors of Krylov-Subspace K(A^T,v0)
    Vec p_acc_, Ap_;
    Vec q_;

    double eps_, delta_, rho_, xi_, beta_, sigma_;
    double breakdownEps_;

    double val_loc_[4], val_glob_[4];   // tmp filds for global-reduce operations
    bool init_;                         // flag, if the fist step is been made (for computation of p,q,beta)
    bool breakdown_;                    // flag, if Lanczos breaks in last step down
    bool luckybreakdown_;               // this indicates, if the breakdown is lucky (i.e. v==0 od w==0)

  public:
    ParPreLanczos2CL(const Mat&, const Vec &v0, const Vec &w0, const PC &pc, const ExCL&, const double eps, const double BreakDownEps=DoubleEpsC);
    ParPreLanczos2CL(const PC &pc, const double BreakDownEps=DoubleEpsC);

    bool Step();                                            ///< Do a Lanczos-Step
    void Init(const Mat& A, const Vec &v0, const Vec &w0, const ExCL& ex, const double eps);

    void   SetBreakDownEps(const double eps);               ///< Set new tolerance for breakdown
    double GetBreakDownEps();                               ///< Get tolerance for breakdown

    Vec& GetV()    /*const*/ {return v_;}                   ///< Get distributed version of basisvector of Krylov-subspace K(A,v0)
    Vec& GetAccV() /*const*/ {return v_acc_;}               ///< Get accumulated version of basisvector of Krylov-subspace K(A,v0)
    Vec& GetW()    /*const*/ {return w_;}                   ///< Get distributed version of basisvector of Krylov-subspace K(A^T,w0)
    Vec& GetAccW() /*const*/ {return w_acc_;}               ///< Get accumulated version of basisvector of Krylov-subspace K(A^T,w0)
    Vec& GetAp()   /*const*/ {return Ap_;}                  ///< Get distributed version of A times p
    Vec& GetAccP()           {return p_acc_;}

    double GetDelta() const {return delta_;}                ///< Get last (hat{v},hat{w})
    double GetEps()   const {return eps_;}                  ///< Get last scaling-parameter of p and q
    double GetBeta()  const {return beta_;}                 ///< Get last recurrence paramter for v and w
    double GetSigma() const {return sigma_;}                ///< Get last (A*hat{v},hat{w})
    double GetRho()   const {return rho_;}                  ///< Get scaling of v   (v^T*v)
    double GetXi()    const {return xi_;}                   ///< Get scaling of w   (w^T*w)
    bool Breakdown() const {return breakdown_;}             ///< Check if Lanczos breaks in last step down
    bool LuckyBreakdown() const {return luckybreakdown_;}   ///< Check if in last step v or w is zero
};

} // end of namespace DROPS

#include "num/parlanczos.tpp"

#endif
