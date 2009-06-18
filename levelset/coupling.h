/// \file
/// \brief coupling of levelset and (Navier-)Stokes equations
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM, Oliver Fortmeier, SC

#ifndef DROPS_COUPLING_H
#define DROPS_COUPLING_H

#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "num/MGsolver.h"
#include "num/nssolver.h"
#include <vector>
#ifdef _PAR
#include "num/parstokessolver.h"
#endif

namespace DROPS
{

template <class StokesT>
class TimeDisc2PhaseCL
{
  protected:
    StokesT&        Stokes_;
    LevelsetP2CL&   LvlSet_;

    VelVecDescCL *b_, *old_b_;        // rhs + couplings with poisson matrix A
    VelVecDescCL *cplM_, *old_cplM_;  // couplings with mass matrix M
    VelVecDescCL *cplN_, *old_cplN_;  // couplings with convection matrix N
    VecDescCL    *curv_, *old_curv_;  // curvature term
    VectorCL      rhs_, ls_rhs_;
    MLMatrixCL*   mat_;               // navier-stokes upper left block
    MatrixCL*     L_;                 // level set system

    double stk_theta_, ls_theta_, dt_;
    const double nonlinear_;

    VecDescCL    cplLB_;
    MLMatDescCL  LB_;

  public:
    TimeDisc2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, double stk_theta= 0.5, double sl_theta = 0.5, double nonlinear=1.);
    virtual ~TimeDisc2PhaseCL();

    double GetStokesTheta()    const { return stk_theta_; }
    double GetLsetTheta()      const { return ls_theta_; }
    double GetTime()           const { return Stokes_.t; }
    double GetTimeStep()       const { return dt_; }
    MLMatrixCL*       GetUpperLeftBlock ()       { return mat_; }
    const MLMatrixCL* GetUpperLeftBlock () const { return mat_; }

    /// \name Get reference on Stokes and level set classes
    //@{
    const StokesT& GetStokes() const { return Stokes_;}
    StokesT& GetStokes() { return Stokes_;}

    const LevelsetP2CL& GetLset() const { return LvlSet_; }
    LevelsetP2CL& GetLset() { return LvlSet_; }
    //@}

    virtual void SetTimeStep (double dt) {dt_= dt;}

    virtual void DoStep( int maxFPiter= -1) = 0;

    // update after grid has changed
    virtual void Update() = 0;
};

template <class StokesT, class LsetSolverT>
class LinThetaScheme2PhaseCL: public TimeDisc2PhaseCL<StokesT>
{
  private:
    typedef TimeDisc2PhaseCL<StokesT> base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_; // 1./dt*M + stk_theta*A + stab_*_theta*_dt*LB
    using base_::L_;   // 1./dt*E + ls_theta*H
    using base_::stk_theta_;
    using base_::ls_theta_;
    using base_::nonlinear_;
    using base_::dt_;
    using base_::cplLB_;
    using base_::LB_;

    StokesSolverT& solver_;
    LsetSolverT&   lsetsolver_;
    VecDescCL    *cplA_;
    bool         implCurv_;

  public:
    LinThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
    		                StokesSolverT& solver, LsetSolverT& lsetsolver, double stk_theta= 0.5, double ls_theta = 0.5,
                            double nonlinear= 1., bool implicitCurv= false);
    ~LinThetaScheme2PhaseCL();

    void SetTimeStep (double dt) { // overwrites base-class-method
        base_::SetTimeStep( dt);
    }
    void SetTimeStep( double dt, double theta) {
        base_::SetTimeStep( dt);
        if (theta >=0 ) stk_theta_= theta;
        ls_theta_ = theta;
    }

    void SolveLsNs();
    void CommitStep();

    void DoStep( int = -1);

    void Update();
};

template <class StokesT, class LsetSolverT>
class OperatorSplitting2PhaseCL : public TimeDisc2PhaseCL<StokesT>
{
  private:
    typedef TimeDisc2PhaseCL<StokesT> base_;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_; // 1./dt*M + theta*A
    using base_::L_;
    using base_::stk_theta_;
    using base_::ls_theta_;
    using base_::dt_;
    using base_::nonlinear_;

    StokesSolverBaseCL&     solver_;
    LsetSolverT&            lsetsolver_;
    SSORPcCL                pc_;
    GMResSolverCL<SSORPcCL> gm_;

    MLMatrixCL    AN_;                // A + N
    VelVecDescCL *cplA_, *old_cplA_;  // couplings with stiff matrix A

    double       alpha_;
    Uint         iter_nonlinear_;


  public:
    OperatorSplitting2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
    		StokesSolverBaseCL& solver, LsetSolverT& lsetsolver, int gm_iter, double gm_tol, double nonlinear= 1);
    ~OperatorSplitting2PhaseCL();

    void SetTimeStep( double dt) { base_::SetTimeStep(dt); }

    void InitStep( bool StokesStep= true);
    // perform fixed point iteration
    void DoStokesFPIter();
    void DoNonlinearFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

/// \brief Compute the relaxation factor in RecThetaScheme2PhaseCL by Aitken's delta-squared method.
///
/// This vector version of classical delta-squared convergence-acceleration computes the
/// relaxation factor in span{ (v, phi)^T}.
class cplDeltaSquaredPolicyCL
{
  private:
    VectorCL v_old_,  phi_old_,
             v_diff_, phi_diff_;
    double   omega_;
    bool     firststep_;
    std::ostream* output_;

  public:
    cplDeltaSquaredPolicyCL( std::ostream* output = 0) :
        omega_( 1.0), firststep_( true), output_( output) {}

    inline void Update( VecDescCL& v, VecDescCL& phi);
};

/// \brief Always uses 1 as relaxation factor in RecThetaScheme2PhaseCL
class cplFixedPolicyCL
{
  private:
    std::ostream* output_;
  public:
    cplFixedPolicyCL  ( std::ostream* output = 0) : output_( output)  {}
    inline void Update( const VecDescCL&, const VecDescCL&) {}
};

/// \brief Broyden method for nonlinear system (velocity - levelset)
///
/// Deuflhard: Newton Methods for Nonlinear Problems,
/// Affine Invariance and Adaptive Algorithms, pp 81-90
class cplBroydenPolicyCL
{
  private:
    double thetamax_;
    double kappamax_;
    double sigma0_, sigma_;
    double kappa_;
    bool   firststep_;
    double tol_;
    std::ostream* output_;

    typedef std::vector<VectorCL> VecValueT;

    VecValueT F1_, F2_, deltaF1_, deltaF2_;
    std::vector<double> gamma_;

  public:
    cplBroydenPolicyCL ( std::ostream* output = 0) : thetamax_( 0.45), kappamax_(10000), kappa_( 1.0), firststep_( true), tol_( 1e-99), output_( output) {}
    inline void Update( VecDescCL&, VecDescCL&);
};


template <class StokesT, class LsetSolverT, class RelaxationPolicyT= cplDeltaSquaredPolicyCL>
class RecThetaScheme2PhaseCL: public TimeDisc2PhaseCL<StokesT>
{
  protected:
    typedef TimeDisc2PhaseCL<StokesT> base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_;      // 1./dt*M + theta*A + stab_*_theta*_dt*LB
    using base_::stk_theta_;
    using base_::ls_theta_;
    using base_::dt_;
    using base_::nonlinear_;
    using base_::cplLB_;
    using base_::LB_;
    using base_::L_;

    VectorCL vdot_, // time derivative of v
             oldv_, // old velocity
             phidot_, // time derivate of phi
             oldphi_; // old level set

    StokesSolverT& solver_;
    LsetSolverT&   lsetsolver_;

    bool         withProj_;
    const double stab_;
    MLMatrixCL*  Mold_;
    MatrixCL*    Eold_;
    bool         trapezoid_;

#ifndef _PAR
    SSORPcCL ssorpc_;
    PCG_SsorCL Msolver_;
    ISBBTPreCL ispc_;
    GCRSolverCL<ISBBTPreCL> Ssolver_;
#else
    typedef ParJac0CL MsolverPCT;
    typedef ParPCGSolverCL<MsolverPCT> MsolverT;
    MsolverPCT MsolverPC_;
    MsolverT   Msolver_;

    typedef ISBBTPreCL SsolverPCT;
    typedef ParPreGMResSolverCL<SsolverPCT> SsolverT;
    SsolverPCT SsolverPC_;
    SsolverT   Ssolver_;
#endif

    void MaybeStabilize (VectorCL&);
    void ComputePressure ();

    void ComputeDots ();
    void EvalLsetNavStokesEquations();

  public:
    RecThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
    		             StokesSolverT& solver, LsetSolverT& lsetsolver, double stk_theta= 0.5, double ls_theta = 0.5, double nonlinear= 1.,
                         bool withProjection= false, double stab= 0.0, bool trapezoid= false);
    ~RecThetaScheme2PhaseCL();

    void SetTimeStep (double dt) { // overwrites baseclass-version
        base_::SetTimeStep( dt);
    }
    void SetTimeStep (double dt, double theta) { // for the fractional-step-method
        base_::SetTimeStep( dt);
        stk_theta_= theta;
        ls_theta_ = theta;
    }

    void InitStep();
    void DoProjectionStep(const VectorCL&);
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

template <class StokesT, class LsetSolverT, class RelaxationPolicyT= cplDeltaSquaredPolicyCL>
class CrankNicolsonScheme2PhaseCL: public RecThetaScheme2PhaseCL<StokesT, LsetSolverT, RelaxationPolicyT>
{
  private:
    typedef RecThetaScheme2PhaseCL<StokesT, LsetSolverT, RelaxationPolicyT> base_;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::mat_;
    using base_::dt_;
    double tmpdt_;
    int step_;

  public:
    CrankNicolsonScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
    		             NSSolverBaseCL<StokesT>& solver, LsetSolverT& lsetsolver, double nonlinear= 1.,
                         bool withProjection= false, double stab= 0.0);
    ~CrankNicolsonScheme2PhaseCL();

    void DoStep( int maxFPiter= -1);
    void Update();

};

template <class StokesT, class LsetSolverT, class RelaxationPolicyT= cplDeltaSquaredPolicyCL>
class FracStepScheme2PhaseCL : public RecThetaScheme2PhaseCL<StokesT, LsetSolverT, RelaxationPolicyT>
{
  private:
    static const double facdt_[3];
    static const double theta_[3];

    typedef RecThetaScheme2PhaseCL<StokesT, LsetSolverT, RelaxationPolicyT> base_;

    double dt3_;
    int step_;

  public:
    FracStepScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
    	                       NSSolverBaseCL<StokesT>& solver, LsetSolverT& lsetsolver, double nonlinear= 1, bool withProjection= false,
                               double stab= 0.0, int step = -1)
        : base_( Stokes, ls, solver, lsetsolver, 0.5, 0.5, nonlinear, withProjection, stab), step_((step >= 0) ? step%3 : 0) {}

    double GetSubTimeStep() const { return facdt_[step_]*dt3_; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep( double dt, int step = -1) {
        dt3_= dt;
        if (step>=0) step_= step%3;
    }

    void DoSubStep( int maxFPiter= -1) {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fractional Step Method: Substep " << step_ << '\n';
        base_::SetTimeStep( GetSubTimeStep(), GetSubTheta());
        base_::DoStep( maxFPiter);
        step_= (step_ + 1)%3;
    }

    void DoStep( int maxFPiter= -1) {
        DoSubStep( maxFPiter);
        DoSubStep( maxFPiter);
        DoSubStep( maxFPiter);
    }

    void Update() { base_::Update(); }
};

template <class NavStokesT, class LsetSolverT, class RelaxationPolicyT>
const double FracStepScheme2PhaseCL<NavStokesT, LsetSolverT, RelaxationPolicyT>::facdt_[3]
//  = { 1./3, 1./3, 1./3 };
//  = { 1./3, 1./3, 1./3 };
  = { 1.0 - std::sqrt( 0.5), std::sqrt( 2.0) - 1.0, 1.0 - std::sqrt( 0.5) };

template <class NavStokesT, class LsetSolverT, class RelaxationPolicyT>
const double FracStepScheme2PhaseCL<NavStokesT,LsetSolverT,RelaxationPolicyT>::theta_[3]
//  = { 1.0, 1.0, 1.0 };
//  = { 1./3, 5./6, 1./3 };
  = { 2.0 - std::sqrt( 2.0), std::sqrt( 2.0) - 1.0, 2.0 - std::sqrt( 2.0) };

} // end of namespace DROPS

#include "levelset/coupling.tpp"

#endif
