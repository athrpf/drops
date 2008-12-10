//**************************************************************************
// File:    coupling.h                                                     *
// Content: coupling of levelset and (Navier-)Stokes equations             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_COUPLING_H
#define DROPS_COUPLING_H

#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "num/MGsolver.h"

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
    MLMatrixCL*   mat_;               // 1./dt*M + theta*A

    double theta_, dt_;
    const double nonlinear_;

    VecDescCL    cplLB_;
    MLMatDescCL  LB_;

  public:
    TimeDisc2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, double theta= 0.5, double nonlinear=1.);
    virtual ~TimeDisc2PhaseCL();

    double GetTheta()    const { return theta_; }
    double GetTime()     const { return Stokes_.t; }
    double GetTimeStep() const { return dt_; }
    const MatrixCL* GetUpperLeftBlock () const { return mat_->GetFinest(); }

    virtual void SetTimeStep (double dt) {dt_= dt;}

    virtual void DoStep( int maxFPiter= -1) = 0;

    // update after grid has changed
    virtual void Update() = 0;
};

template <class StokesT, class SolverT>
class LinThetaScheme2PhaseCL: public TimeDisc2PhaseCL<StokesT>
{
  private:
    typedef TimeDisc2PhaseCL<StokesT> base_;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_; // 1./dt*M + theta*A + stab_*_theta*_dt*LB
    using base_::theta_;
    using base_::nonlinear_;
    using base_::dt_;
    using base_::cplLB_;
    using base_::LB_;

    SolverT&     solver_;
    VecDescCL    *cplA_;
    bool         implCurv_;

  public:
    LinThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                            SolverT& solver, double theta= 0.5,
                            double nonlinear= 1., bool implicitCurv= false);
    ~LinThetaScheme2PhaseCL();

    void SetTimeStep (double dt) { // overwrites base-class-method
        base_::SetTimeStep( dt);
        LvlSet_.SetTimeStep( dt);
    }
    void SetTimeStep( double dt, double theta) {
        base_::SetTimeStep( dt);
        if (theta >=0 ) theta_= theta;
        LvlSet_.SetTimeStep( dt, theta);
    }

    void SolveLsNs();
    void CommitStep();

    void DoStep( int = -1);

    void Update();
};

template <class StokesT, class SolverT>
class ThetaScheme2PhaseCL: public TimeDisc2PhaseCL<StokesT>
{
  private:
    typedef TimeDisc2PhaseCL<StokesT> base_;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_; // 1./dt*M + theta*A + stab_*_theta*_dt*LB
    using base_::theta_;
    using base_::nonlinear_;
    using base_::dt_;
    using base_::cplLB_;
    using base_::LB_;

    SolverT&     solver_;
    bool         withProj_;
    const double stab_;

    void MaybeStabilize (VectorCL&);
    void ComputePressure ();

  public:
    ThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                         SolverT& solver, double theta= 0.5, double nonlinear= 1.,
                         bool withProjection= false, double stab= 0.0);
    ~ThetaScheme2PhaseCL();

    void SetTimeStep (double dt) { // overwrites base-class-method
        base_::SetTimeStep( dt);
        LvlSet_.SetTimeStep( dt);
    }
    void SetTimeStep (double dt, double theta) { // for fractional-step
        base_::SetTimeStep( dt);
        theta_= theta;
        LvlSet_.SetTimeStep( dt, theta);
    }

    void InitStep();
    void DoProjectionStep(const VectorCL&);
    void DoFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

template <class StokesT, class SolverT>
class FracStepScheme2PhaseCL : public ThetaScheme2PhaseCL<StokesT,SolverT>
{
  private:
    static const double facdt_[3];
    static const double theta_[3];

    typedef ThetaScheme2PhaseCL<StokesT,SolverT> base_;

    double dt3_;
    int step_;

  public:
    FracStepScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                               SolverT& solver, double nonlinear= 1, bool withProjection= false,
                               double stab= 0.0, int step = -1)
        : base_( Stokes, ls, solver, 0.5, nonlinear, withProjection, stab), step_((step >= 0) ? step%3 : 0) {}

    double GetSubTimeStep() const { return facdt_[step_]*dt3_; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep( double dt, int step = -1) {
        dt3_= dt;
        if (step>=0) step_= step%3;
    }

    void DoSubStep( int maxFPiter= -1) {
        std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fractional Step Method: Substep " << step_ << '\n';
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

template <class NavStokesT, class SolverT>
const double FracStepScheme2PhaseCL<NavStokesT,SolverT>::facdt_[3]
//  = { 1./3, 1./3, 1./3 };
//  = { 1./3, 1./3, 1./3 };
  = { 1.0 - std::sqrt( 0.5), std::sqrt( 2.0) - 1.0, 1.0 - std::sqrt( 0.5) };

template <class NavStokesT, class SolverT>
const double FracStepScheme2PhaseCL<NavStokesT,SolverT>::theta_[3]
//  = { 1.0, 1.0, 1.0 };
//  = { 1./3, 5./6, 1./3 };
  = { 2.0 - std::sqrt( 2.0), std::sqrt( 2.0) - 1.0, 2.0 - std::sqrt( 2.0) };


template <class StokesT, class SolverT>
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
    using base_::theta_;
    using base_::dt_;
    using base_::nonlinear_;

    SolverT&                solver_;
    SSORPcCL                pc_;
    GMResSolverCL<SSORPcCL> gm_;

    MLMatrixCL    AN_;                // A + N
    VelVecDescCL *cplA_, *old_cplA_;  // couplings with stiff matrix A

    double       alpha_;
    Uint         iter_nonlinear_;


  public:
    OperatorSplitting2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
        SolverT& solver, int gm_iter, double gm_tol, double nonlinear= 1);
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


template <class StokesT, class SolverT>
class RecThetaScheme2PhaseCL: public TimeDisc2PhaseCL<StokesT>
{
  protected:
    typedef TimeDisc2PhaseCL<StokesT> base_;
    using base_::Stokes_;
    using base_::LvlSet_;
    using base_::b_;       using base_::old_b_;
    using base_::cplM_;    using base_::old_cplM_;
    using base_::cplN_;    using base_::old_cplN_;
    using base_::curv_;    using base_::old_curv_;
    using base_::rhs_;
    using base_::ls_rhs_;
    using base_::mat_; // 1./dt*M + theta*A + stab_*_theta*_dt*LB
    using base_::theta_;
    using base_::dt_;
    using base_::nonlinear_;
    using base_::cplLB_;
    using base_::LB_;

    VectorCL vdot_, // time derivative of v
             oldv_; // old velocity

    SolverT&     solver_;
    bool         withProj_;
    const double stab_;

    void MaybeStabilize (VectorCL&);
    void ComputePressure ();

    void ComputeVelocityDot ();

  public:
    RecThetaScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                         SolverT& solver, double theta= 0.5, double nonlinear= 1.,
                         bool withProjection= false, double stab= 0.0);
    ~RecThetaScheme2PhaseCL();

    void SetTimeStep (double dt) { // overwrites baseclass-version
        base_::SetTimeStep( dt);
        LvlSet_.SetTimeStep( dt);
    }
    void SetTimeStep (double dt, double theta) { // for the fractional-step-method
        base_::SetTimeStep( dt);
        theta_= theta;
        LvlSet_.SetTimeStep( dt, theta_);
    }

    void InitStep();
    void DoProjectionStep(const VectorCL&);
    void DoFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

template <class StokesT, class SolverT>
class CrankNicolsonScheme2PhaseCL: public RecThetaScheme2PhaseCL<StokesT, SolverT>
{
  private:
    typedef RecThetaScheme2PhaseCL<StokesT, SolverT> base_;
    using base_::Stokes_;
    using base_::mat_;
    using base_::dt_;
    double tmpdt_;
    int step_;

  public:
    CrankNicolsonScheme2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                         SolverT& solver, double nonlinear= 1.,
                         bool withProjection= false, double stab= 0.0);
    ~CrankNicolsonScheme2PhaseCL();

    void DoStep( int maxFPiter= -1);
    void Update();

};

} // end of namespace DROPS

#include "levelset/coupling.tpp"

#endif
