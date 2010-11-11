/// \file coupling.h
/// \brief coupling of levelset and (Navier-)Stokes equations
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_COUPLING_H
#define DROPS_COUPLING_H

#include "navstokes/instatnavstokes2phase.h"
#include "levelset/levelset.h"
#include "num/MGsolver.h"
#include "num/nssolver.h"
#include <vector>
#ifdef _PAR
#include "num/parstokessolver.h"
#endif

namespace DROPS
{

typedef InstatNavierStokes2PhaseP2P1CL StokesT;

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

    double       dt_;
    const double nonlinear_;

    LevelsetModifyCL& lsetmod_;

    VecDescCL    cplLB_;
    MLMatDescCL  LB_;

    SchurPreBaseCL* ispc_;             // pointer to preconditioner for the schur complement

  public:
    TimeDisc2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, LevelsetModifyCL& lsetmod, double nonlinear=1.);
    virtual ~TimeDisc2PhaseCL();

    double GetTime()           const { return Stokes_.v.t; }
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
    
    void SetSchurPrePtr( SchurPreBaseCL* ptr) { ispc_ = ptr; }
};

template <class LsetSolverT>
class LinThetaScheme2PhaseCL: public TimeDisc2PhaseCL
{
  private:
    typedef TimeDisc2PhaseCL base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;

    StokesSolverT& solver_;
    LsetSolverT&   lsetsolver_;

    double stk_theta_, ls_theta_;

    VecDescCL    *cplA_;
    bool         implCurv_;

  public:
    LinThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                                StokesSolverT& solver, LsetSolverT& lsetsolver,
                                LevelsetModifyCL& lsetmod, double stk_theta= 0.5, double ls_theta = 0.5,
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

    double GetStokesTheta()    const { return stk_theta_; }
    double GetLsetTheta()      const { return ls_theta_; }

};

template <class LsetSolverT>
class OperatorSplitting2PhaseCL : public TimeDisc2PhaseCL
{
  private:
    typedef TimeDisc2PhaseCL base_;

    StokesSolverBaseCL&     solver_;
    LsetSolverT&            lsetsolver_;
    SSORPcCL                pc_;
    GMResSolverCL<SSORPcCL> gm_;

    double stk_theta_, ls_theta_;

    MLMatrixCL    AN_;                // A + N
    VelVecDescCL *cplA_, *old_cplA_;  // couplings with stiff matrix A

    double       alpha_;
    Uint         iter_nonlinear_;


  public:
    OperatorSplitting2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                                StokesSolverBaseCL& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod, int gm_iter, double gm_tol, double nonlinear= 1);
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


//forward declarations
class cplDeltaSquaredPolicyCL;
class cplFixedPolicyCL;
class cplBroydenPolicyCL;

template <class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class CoupledTimeDisc2PhaseBaseCL: public TimeDisc2PhaseCL
{
  protected:
    typedef TimeDisc2PhaseCL base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;

    double dphi_;

    StokesSolverT& solver_;
    LsetSolverT&   lsetsolver_;
    double         tol_;

    bool           withProj_;
    double         stab_;
    double         alpha_;

    virtual void InitStep();
    virtual void CommitStep();
    virtual void SetupNavStokesSystem() = 0;
    virtual void SetupLevelsetSystem()  = 0;

    void DoProjectionStep( const VectorCL&);
    void MaybeStabilize  ( VectorCL&);
    void EvalLsetNavStokesEquations();
    void SetupStokesMatVec();

  public:
    CoupledTimeDisc2PhaseBaseCL( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver,
                                 LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod, double tol,
                                 double nonlinear= 1., bool withProjection= false, double stab= 0.0);
    ~CoupledTimeDisc2PhaseBaseCL();

    void DoStep( int maxFPiter= -1);

    virtual void Update();
};

template <class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class MidPointTimeDisc2PhaseCL: public CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT>
{
  protected:
    typedef CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT> base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_;      // NavStokes Operator
    using base_::L_;        // Levelset Operator
    using base_::dt_;
    using base_::nonlinear_;
    using base_::cplLB_;
    using base_::LB_;
    using base_::lsetmod_;
    using base_::alpha_;

  private:
    VectorCL oldv_,     // old velocity
             oldp_,     // old pressure
             oldphi_;   // old levelset

    bool implicitpressure_;

    void InitStep();
    void CommitStep();
    void SetupNavStokesSystem();
    void SetupLevelsetSystem();

  public:
    MidPointTimeDisc2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver,
                              LevelsetModifyCL& lsetmod, double tol, double nonlinear = 1.0,
                              bool withProjection =  false, double stab = 0.0, bool implicitpressure = false);
    ~MidPointTimeDisc2PhaseCL() {}

    void Update();

};

template <class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class SpaceTimeDiscTheta2PhaseCL: public CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT>
{
  protected:
    typedef CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT> base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_;      // NavStokes Operator
    using base_::L_;        // Levelset Operator
    using base_::dt_;
    using base_::nonlinear_;
    using base_::cplLB_;
    using base_::LB_;
    using base_::lsetmod_;
    using base_::alpha_;
    using base_::stab_;

    void ComputeDots ();
    VectorCL oldv_,     // old velocity
             oldphi_;   // old levelset
  private:
    VectorCL phidot_,   // "time derivate" of phi
             vdot_;     // "time derivate" of v

    VectorCL fixed_rhs_, fixed_ls_rhs_;

    double stk_theta_, ls_theta_;
    bool implicitpressure_;

    MLMatrixCL*  Mold_;
    MatrixCL*    Eold_;

    void InitStep();
    void CommitStep();
    void SetupNavStokesSystem();
    void SetupLevelsetSystem();

  public:
    SpaceTimeDiscTheta2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver,
                              LevelsetModifyCL& lsetmod, double tol, double stk_theta= 0.5, double ls_theta = 0.5, double nonlinear = 1.0,
                              bool withProjection =  false, double stab = 0.0, bool implicitpressure = false);
    ~SpaceTimeDiscTheta2PhaseCL();

    void Update();

    void SetTimeStep (double dt) { // overwrites baseclass-version
        base_::SetTimeStep( dt);
    }
    void SetTimeStep (double dt, double theta) { // for the fractional-step-method
        base_::SetTimeStep( dt);
        if (base_::ispc_)
            base_::ispc_->SetWeights(1.0/dt, theta);
        stk_theta_= theta;
        ls_theta_ = theta;
    }
};

template <class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class EulerBackwardScheme2PhaseCL: public CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT>
{
  protected:
    typedef CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT> base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;
    using base_::cplM_;
    using base_::cplN_;
    using base_::curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_;      // NavStokes Operator
    using base_::L_;        // Levelset Operator
    using base_::dt_;
    using base_::nonlinear_;
    using base_::cplLB_;
    using base_::LB_;
    using base_::lsetmod_;
    using base_::alpha_;


  private:
    VectorCL fixed_rhs_, fixed_ls_rhs_;

    void InitStep();
    void CommitStep();
    void SetupNavStokesSystem();
    void SetupLevelsetSystem();

  public:
    EulerBackwardScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                         StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
                         double tol, double nonlinear= 1., bool withProjection= false, double stab= 0.0);
    ~EulerBackwardScheme2PhaseCL();

    void Update();
};

template <class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class RecThetaScheme2PhaseCL: public CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT>
{
  protected:
    typedef CoupledTimeDisc2PhaseBaseCL<LsetSolverT, RelaxationPolicyT> base_;
    typedef NSSolverBaseCL<StokesT> StokesSolverT;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_;      // NavStokes Operator
    using base_::L_;        // Levelset Operator
    using base_::dt_;
    using base_::nonlinear_;
    using base_::cplLB_;
    using base_::LB_;
    using base_::lsetmod_;
    using base_::alpha_;
    using base_::stab_;

    void ComputeDots ();

    VectorCL oldv_,     // old velocity
             oldphi_;   // old levelset

  private:
    VectorCL fixed_rhs_, fixed_ls_rhs_;
    double stk_theta_, ls_theta_;

    void InitStep();
    void CommitStep();
    void SetupNavStokesSystem();
    void SetupLevelsetSystem();

    VectorCL vdot_,     // time derivative of v
             phidot_;   // time derivate of phi

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

    void ComputePressure ();

  public:
    RecThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                             StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
                             double tol, double stk_theta= 0.5, double ls_theta = 0.5, double nonlinear= 1.,
                             bool withProjection= false, double stab= 0.0);
    ~RecThetaScheme2PhaseCL();

    void SetTimeStep (double dt) { // overwrites baseclass-version
        base_::SetTimeStep( dt);
    }
    void SetTimeStep (double dt, double theta) { // for the fractional-step-method
        base_::SetTimeStep( dt);
        if (base_::ispc_)
            base_::ispc_->SetWeights(1.0/dt, theta);
        stab_ *= theta/stk_theta_;
        stk_theta_= theta;
        ls_theta_ = theta;
    }

    void Update();
};

template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class CrankNicolsonScheme2PhaseCL: public BaseMethod<LsetSolverT, RelaxationPolicyT>
{
  private:
    static const double theta_[3];
    double facdt_[3];
    typedef BaseMethod<LsetSolverT, RelaxationPolicyT> base_;
    double dt2_;
    int step_;

  public:
    CrankNicolsonScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                             NSSolverBaseCL<StokesT>& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
                             double tol, double nonlinear= 1., bool withProjection= false, double stab= 0.0, int step = -1)
        : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, 1.0, 1.0, nonlinear, withProjection, stab), step_((step >= 0) ? step%2 : 0) {}

    ~CrankNicolsonScheme2PhaseCL() {};

    double GetSubTimeStep() const { return facdt_[step_]; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep (double dt) { // overwrites baseclass-version
        dt2_= dt;
        facdt_[0] = 0.2*dt2_*dt2_; facdt_[1] = dt2_-0.2*dt2_*dt2_; facdt_[2] = dt2_;
    }

    void DoSubStep( int maxFPiter= -1) {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Crank Nicolson Method: Substep variant" << step_ << '\n';
        base_::SetTimeStep( GetSubTimeStep(), GetSubTheta());
        base_::DoStep( maxFPiter);
    }

    void DoStep( int maxFPiter= -1) {
        if (step_ == 0) {
            DoSubStep( maxFPiter);
            ++step_;
            DoSubStep( maxFPiter);
            ++step_;
        }
        else
            DoSubStep( maxFPiter);
    }

    void Update() { step_ = 0; base_::Update(); }

};

template < template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
const double CrankNicolsonScheme2PhaseCL<BaseMethod, LsetSolverT, RelaxationPolicyT>::theta_[3]
  = { 1.0, 0.5, 0.5 };

// Handbook of Numerical Analysis (Ciarlet/Lions), Volume IX: Numerical Methods for Fluids (Part 3) by R. Glowinski
template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class FracStepScheme2PhaseCL : public BaseMethod<LsetSolverT, RelaxationPolicyT>
{
  private:
    static const double facdt_[3];
    static const double theta_[3];

    typedef BaseMethod<LsetSolverT, RelaxationPolicyT> base_;

    double dt3_;
    int step_;

  public:
    FracStepScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                               NSSolverBaseCL<StokesT>& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
                               double tol, double nonlinear= 1, bool withProjection= false, double stab= 0.0, int step = -1)
        : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, 0.5, 0.5, nonlinear, withProjection, stab), step_((step >= 0) ? step%3 : 0) {}

    double GetSubTimeStep() const { return facdt_[step_]*dt3_; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep (double dt) { // overwrites baseclass-version
        dt3_= dt;
    }

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

template < template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
const double FracStepScheme2PhaseCL<BaseMethod, LsetSolverT, RelaxationPolicyT>::facdt_[3]
//  = { 1./3, 1./3, 1./3 };
//  = { 1./3, 1./3, 1./3 };
  = { 1.0 - std::sqrt( 0.5), std::sqrt( 2.0) - 1.0, 1.0 - std::sqrt( 0.5) };

template < template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
const double FracStepScheme2PhaseCL<BaseMethod, LsetSolverT, RelaxationPolicyT>::theta_[3]
//  = { 1.0, 1.0, 1.0 };
//  = { 1./3, 5./6, 1./3 };
  = { 2.0 - std::sqrt( 2.0), std::sqrt( 2.0) - 1.0, 2.0 - std::sqrt( 2.0) };


// Handbook of Numerical Analysis (Ciarlet/Lions), Volume IX: Numerical Methods for Fluids (Part 3) by R. Glowinski: Remark 10.3
template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT= cplBroydenPolicyCL>
class Frac2StepScheme2PhaseCL : public BaseMethod<LsetSolverT, RelaxationPolicyT>
{
  private:
    static const double facdt_[2];
    static const double theta_[2];

    typedef BaseMethod<LsetSolverT, RelaxationPolicyT> base_;
    using base_::oldv_;
    using base_::oldphi_;
    using base_::Stokes_;
    using base_::LvlSet_;
    double dt2_;
    int step_;

  public:
    Frac2StepScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                               NSSolverBaseCL<StokesT>& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
                               double tol, double nonlinear= 1, bool withProjection= false, double stab= 0.0, int step = -1)
        : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, 0.5, 0.5, nonlinear, withProjection, stab), step_((step >= 0) ? step%2 : 0) {}

    double GetSubTimeStep() const { return facdt_[step_]*dt2_; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep (double dt) { // overwrites baseclass-version
        dt2_= dt;
    }

    void SetTimeStep( double dt, int step = -1) {
        dt2_= dt;
        if (step>=0) step_= step%2;
    }

    void DoSubStep( int maxFPiter= -1) {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fractional Step Method: Substep " << step_ << '\n';
        base_::SetTimeStep( GetSubTimeStep(), GetSubTheta());
        base_::DoStep( maxFPiter);
        step_= (step_ + 1)%2;
    }

    void DoStep( int maxFPiter= -1) {
        VectorCL v(Stokes_.v.Data);
        VectorCL phi(LvlSet_.Phi.Data);
        double t(Stokes_.v.t);
        DoSubStep( maxFPiter);
        oldv_ = v;
        oldphi_ = phi;
        Stokes_.v.t = t;
        Stokes_.p.t = t;
        LvlSet_.Phi.t = t;
        DoSubStep( maxFPiter);
    }

    void Update() { base_::Update(); }
};

template < template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
const double Frac2StepScheme2PhaseCL<BaseMethod, LsetSolverT, RelaxationPolicyT>::facdt_[2]
  = { 1.0 - std::sqrt( 0.5), 1.0 };

template < template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
const double Frac2StepScheme2PhaseCL<BaseMethod, LsetSolverT, RelaxationPolicyT>::theta_[2]
  = { 1.0, 1.0 - std::sqrt( 0.5) };

/// \brief Compute the relaxation factor in RecThetaScheme2PhaseCL by Aitken's delta-squared method.
///
/// This vector version of classical delta-squared convergence-acceleration computes the
/// relaxation factor in span{ (v, phi)^T}.
class cplDeltaSquaredPolicyCL
{
  private:
    VectorCL v_old_, v_diff_;
    double   omega_;
    bool     firststep_;
    std::ostream* output_;

  public:
    cplDeltaSquaredPolicyCL( std::ostream* output = 0) :
        omega_( 1.0), firststep_( true), output_( output) {}

    void Update( VecDescCL& v);
};

/// \brief Always uses 1 as relaxation factor in RecThetaScheme2PhaseCL
class cplFixedPolicyCL
{
  private:
    std::ostream* output_;

  public:
    cplFixedPolicyCL (std::ostream* output = 0) : output_( output) {}
    void Update (const VecDescCL&) {}
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
    cplBroydenPolicyCL (std::ostream* output = 0) : thetamax_( 0.45), kappamax_(10000), sigma_( -1.), kappa_( 1.0), firststep_( true), tol_( 1e-99), output_( output) {}

    void Update( VecDescCL&);
};


} // end of namespace DROPS

#include "levelset/coupling.tpp"

#endif
