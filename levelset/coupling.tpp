/// \file
/// \brief coupling of levelset and (Navier-)Stokes equations
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM, Oliver Fortmeier, SC

#include "num/nssolver.h"

namespace DROPS
{

// ==============================================
//              TimeDisc2PhaseCL
// ==============================================

template <class StokesT>
TimeDisc2PhaseCL<StokesT>::TimeDisc2PhaseCL (StokesT& Stokes, LevelsetP2CL& ls, double theta, double nonlinear)
  : Stokes_( Stokes), LvlSet_( ls), b_( &Stokes.b), old_b_( new VelVecDescCL),
    cplM_( new VelVecDescCL), old_cplM_( new VelVecDescCL),
    cplN_( new VelVecDescCL), old_cplN_( new VelVecDescCL),
    curv_( new VelVecDescCL), old_curv_(new VelVecDescCL),
    rhs_( Stokes.b.RowIdx->NumUnknowns()), ls_rhs_( ls.Phi.RowIdx->NumUnknowns()),
    mat_( 0), theta_( theta), nonlinear_( nonlinear)
{
    Stokes_.SetLevelSet( ls);
    mat_= new MLMatrixCL( Stokes.vel_idx.size());
    LB_.Data.resize( Stokes.vel_idx.size());
}

template <class StokesT>
TimeDisc2PhaseCL<StokesT>::~TimeDisc2PhaseCL()
{
    delete mat_;
    if (old_b_ == &Stokes_.b)
        delete b_;
    else
        delete old_b_;
    delete cplM_; delete old_cplM_;
    delete cplN_; delete old_cplN_;
    delete curv_; delete old_curv_;
}

// ==============================================
//           LinThetaScheme2PhaseCL
// ==============================================

template <class StokesT, class SolverT>
LinThetaScheme2PhaseCL<StokesT,SolverT>::LinThetaScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double theta, double nonlinear, bool implicitCurv)
  : base_( Stokes, ls, theta, nonlinear), solver_( solver),
    cplA_(new VelVecDescCL), implCurv_( implicitCurv)
{
    Update();
}

template <class StokesT, class SolverT>
LinThetaScheme2PhaseCL<StokesT,SolverT>::~LinThetaScheme2PhaseCL()
{
    delete cplA_;
}

template <class StokesT, class SolverT>
void LinThetaScheme2PhaseCL<StokesT,SolverT>::SolveLsNs()
// solve decoupled level set / Navier-Stokes system
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();

    // operators are computed for old level set
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, old_b_, cplA_, cplM_, LvlSet_, Stokes_.t);
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass( &Stokes_.prM, LvlSet_);

    if (!implCurv_)
        mat_->LinComb( 1./dt_, Stokes_.M.Data, theta_, Stokes_.A.Data);
    else // semi-implicit treatment of curvature term, cf. Baensch
    {
        MLMatrixCL mat0;
        mat0.LinComb( 1./dt_, Stokes_.M.Data, theta_, Stokes_.A.Data);
        // The MatrixBuilderCL's method of determining when to reuse the pattern
        // is not save for matrix LB_
        LB_.Data.clear();
        LB_.Data.resize( Stokes_.A.Data.size());
        LB_.SetIdx( &Stokes_.vel_idx, &Stokes_.vel_idx);
        Stokes_.SetupLB( &LB_, &cplLB_, LvlSet_, Stokes_.t);
        mat_->LinComb( 1., mat0, dt_, LB_.Data);
    }
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing NavierStokes for old level set took "<< duration <<" sec.\n";
    time.Reset();

    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution() );
    LvlSet_.SetTimeStep( dt_); // does LinComb of system matrix L_
    ls_rhs_= (1./dt_)*(LvlSet_.E*LvlSet_.Phi.Data) - (1.-theta_)*(LvlSet_.H*LvlSet_.Phi.Data);

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing Levelset took " << duration << " sec.\n";
    time.Reset();

    LvlSet_.DoLinStep( ls_rhs_);

    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving Levelset took " << duration << " sec.\n";

    time.Reset();

    Stokes_.t+= dt_;

    // rhs for new level set
    curv_->Clear();
    LvlSet_.AccumulateBndIntegral( *curv_);
    Stokes_.SetupRhs1( b_, LvlSet_, Stokes_.t);

    rhs_=  (1./dt_)*(Stokes_.M.Data*Stokes_.v.Data) + theta_*b_->Data + cplA_->Data
        + (1.-theta_)*(old_b_->Data - Stokes_.A.Data*Stokes_.v.Data - nonlinear_*(Stokes_.N.Data*Stokes_.v.Data - old_cplN_->Data));

    if (!implCurv_)
        rhs_+= theta_*curv_->Data + (1.-theta_)*old_curv_->Data;
    else // semi-implicit treatment of curvature term, cf. Baensch
        rhs_+= old_curv_->Data + dt_*cplLB_.Data;

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing Rhs/Curv took "<< duration <<" sec.\n";

    time.Reset();
    solver_.Solve( *mat_, Stokes_.B.Data,
        Stokes_.v, Stokes_.p.Data,
        rhs_, *cplN_, Stokes_.c.Data, /*alpha*/ theta_*nonlinear_);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving NavierStokes: residual: " << solver_.GetResid()
              << "\titerations: " << solver_.GetIter()
              << "\ttime: " << duration << "sec\n";
}

template <class StokesT, class SolverT>
void LinThetaScheme2PhaseCL<StokesT,SolverT>::CommitStep()
{
    if (Stokes_.UsesXFEM()) { // update XFEM
        Stokes_.UpdateXNumbering( &Stokes_.pr_idx, LvlSet_);
        Stokes_.UpdatePressure( &Stokes_.p);
        Stokes_.c.SetIdx( &Stokes_.pr_idx);
        Stokes_.B.SetIdx( &Stokes_.pr_idx, &Stokes_.vel_idx);
        Stokes_.prA.SetIdx( &Stokes_.pr_idx, &Stokes_.pr_idx);
        Stokes_.prM.SetIdx( &Stokes_.pr_idx, &Stokes_.pr_idx);
        // The MatrixBuilderCL's method of determining when to reuse the pattern
        // is not save for P1X-elements.
        Stokes_.B.Data.clear();
        Stokes_.prA.Data.clear();
        Stokes_.prM.Data.clear();
        Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    }
    else
        Stokes_.SetupRhs2( &Stokes_.c, LvlSet_, Stokes_.t);

    std::swap( curv_, old_curv_);
    std::swap( cplN_, old_cplN_);
}

template <class StokesT, class SolverT>
void LinThetaScheme2PhaseCL<StokesT,SolverT>::DoStep( int)
{
    SolveLsNs();
    CommitStep();
}

template <class StokesT, class SolverT>
void LinThetaScheme2PhaseCL<StokesT,SolverT>::Update()
{
    MLIdxDescCL* const vidx= &Stokes_.vel_idx;
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();

    std::cout << "Updating discretization...\n";

    if (implCurv_)
    {
        LB_.Data.clear();
        LB_.SetIdx( vidx, vidx);
        cplLB_.SetIdx( vidx);
    }
    Stokes_.ClearMat();
    LvlSet_.ClearMat();
    if (Stokes_.UsesXFEM()) { // update XFEM
        Stokes_.UpdateXNumbering( &Stokes_.pr_idx, LvlSet_);
        Stokes_.UpdatePressure( &Stokes_.p);
        // The MatrixBuilderCL's method of determining when to reuse the pattern
        // is not save for P1X-elements.
        Stokes_.B.Data.clear();
        Stokes_.prA.Data.clear();
        Stokes_.prM.Data.clear();
    }
    // IndexDesc setzen
    cplA_->SetIdx( vidx);    cplM_->SetIdx( vidx);
    old_b_->SetIdx( vidx);
    cplN_->SetIdx( vidx);    old_cplN_->SetIdx( vidx);
    curv_->SetIdx( vidx);    old_curv_->SetIdx( vidx);
    rhs_.resize( vidx->NumUnknowns());
    ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
    Stokes_.SetIdx();

    // Diskretisierung
    LvlSet_.AccumulateBndIntegral( *old_curv_);
    LvlSet_.SetupSystem( Stokes_.GetVelSolution() );
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}

// ==============================================
//              OperatorSplitting2PhaseCL
// ==============================================

template <class StokesT, class SolverT>
OperatorSplitting2PhaseCL<StokesT,SolverT>::OperatorSplitting2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, int gm_iter, double gm_tol, double nonlinear)

  : base_( Stokes, ls, /*theta*/ 1.0 - std::sqrt( 2.)/2., nonlinear), solver_(solver),
    gm_( pc_, 100, gm_iter, gm_tol, false /*test absolute resid*/),
    cplA_( new VelVecDescCL), old_cplA_( new VelVecDescCL),
    alpha_( (1.0 - 2.0*theta_)/(1.0 - theta_))
{
#ifdef _PAR
    throw DROPSErrCL("OperatorSplitting2PhaseCL: Not parallelized, yet, sorry");
#endif
    std::cout << "theta = " << theta_ << "\talpha = " << alpha_ << std::endl;
    Update();
}

template <class StokesT, class SolverT>
OperatorSplitting2PhaseCL<StokesT,SolverT>::~OperatorSplitting2PhaseCL()
{
    delete cplA_; delete old_cplA_;
}

template <class StokesT, class SolverT>
void OperatorSplitting2PhaseCL<StokesT,SolverT>::InitStep( bool StokesStep)
{
// compute all terms that don't change during the following FP iterations

    const double fracdt_= StokesStep ? theta_*dt_ : (1-2*theta_)*dt_;

    LvlSet_.SetTimeStep( fracdt_, StokesStep ? alpha_ : 1-alpha_);
    LvlSet_.ComputeRhs( ls_rhs_);

    Stokes_.t+= fracdt_;

    if (StokesStep)
    {
        rhs_= -(1-alpha_)*(Stokes_.A.Data * Stokes_.v.Data - old_cplA_->Data)
              -nonlinear_*(Stokes_.N.Data * Stokes_.v.Data - old_cplN_->Data);
    }
    else
    {
        rhs_= -alpha_*(Stokes_.A.Data * Stokes_.v.Data - old_cplA_->Data)
              - transp_mul( Stokes_.B.Data, Stokes_.p.Data);
    }
    rhs_+= (1./fracdt_)*(Stokes_.M.Data*Stokes_.v.Data - old_cplM_->Data);

}

template <class StokesT, class SolverT>
void OperatorSplitting2PhaseCL<StokesT,SolverT>::DoStokesFPIter()
// perform fixed point iteration: Levelset / Stokes
// for fractional steps A and C
{
    TimerCL time;
    time.Reset();
    time.Start();
    const double fracdt_= theta_*dt_;
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution() );
    LvlSet_.SetTimeStep( fracdt_, alpha_);
    time.Stop();
    std::cout << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    LvlSet_.DoStep( ls_rhs_);
    time.Stop();
    std::cout << "Solving Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    time.Start();
    curv_->Clear();
    LvlSet_.AccumulateBndIntegral( *curv_);
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, b_, cplA_, cplM_, LvlSet_, Stokes_.t);
    mat_->LinComb( 1./fracdt_, Stokes_.M.Data, alpha_, Stokes_.A.Data);
    if (Stokes_.UsesXFEM()) {
        Stokes_.UpdateXNumbering( &Stokes_.pr_idx, LvlSet_);
        Stokes_.UpdatePressure( &Stokes_.p);
        Stokes_.c.SetIdx( &Stokes_.pr_idx);
        Stokes_.B.SetIdx( &Stokes_.pr_idx, &Stokes_.vel_idx);
        Stokes_.prA.SetIdx( &Stokes_.pr_idx, &Stokes_.pr_idx);
        Stokes_.prM.SetIdx( &Stokes_.pr_idx, &Stokes_.pr_idx);
        // The MatrixBuilderCL's method of determining when to reuse the pattern
        // is not save for P1X-elements.
        Stokes_.B.Data.clear();
        Stokes_.prA.Data.clear();
        Stokes_.prM.Data.clear();
        Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    }
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass ( &Stokes_.prM, LvlSet_);
    time.Stop();
    std::cout << "Discretizing Stokes/Curv took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    solver_.Solve( *mat_, Stokes_.B.Data, Stokes_.v.Data, Stokes_.p.Data,
                   VectorCL( rhs_ + (1./fracdt_)*cplM_->Data + alpha_*cplA_->Data + curv_->Data + b_->Data), Stokes_.c.Data);
    time.Stop();
    std::cout << "Solving Stokes took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void OperatorSplitting2PhaseCL<StokesT,SolverT>::DoNonlinearFPIter()
// perform fixed point iteration: Levelset / nonlinear system
// for fractional step B
{
    const double fracdt_= dt_*(1.-2.*theta_);
    TimerCL time;
    time.Reset();
    time.Start();
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution() );
    LvlSet_.SetTimeStep( fracdt_, 1-alpha_);
    time.Stop();
    std::cout << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    LvlSet_.DoStep( ls_rhs_);
    time.Stop();
    std::cout << "Solving Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    time.Start();
    curv_->Clear();
    LvlSet_.AccumulateBndIntegral( *curv_);
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, b_, cplA_, cplM_, LvlSet_, Stokes_.t);
    time.Stop();
    std::cout << "Discretizing Stokes/Curv took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    std::cout << "Starting fixed point iterations for solving nonlinear system...\n";
    iter_nonlinear_= 0;
    do
    {
        Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, cplN_, LvlSet_, Stokes_.t);
        AN_.LinComb( 1-alpha_, Stokes_.A.Data, nonlinear_, Stokes_.N.Data);
        mat_->LinComb( 1./fracdt_, Stokes_.M.Data, 1., AN_);
        gm_.Solve( *mat_, Stokes_.v.Data,
            VectorCL( rhs_ + (1./fracdt_)*cplM_->Data + (1-alpha_)*cplA_->Data
            + nonlinear_*cplN_->Data + curv_->Data + b_->Data));
        std::cout << "fp cycle " << ++iter_nonlinear_ << ":\titerations: "
                  << gm_.GetIter() << "\tresidual: " << gm_.GetResid() << std::endl;
    } while (gm_.GetIter() > 0 && iter_nonlinear_<20);
    time.Stop();
    std::cout << "Solving nonlinear system took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void OperatorSplitting2PhaseCL<StokesT,SolverT>::CommitStep()
{
    std::swap( cplM_, old_cplM_);
    std::swap( cplA_, old_cplA_);
    std::swap( cplN_, old_cplN_);
}

template <class StokesT, class SolverT>
void OperatorSplitting2PhaseCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 99;

    // ------ frac. step A ------
    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoStokesFPIter();
        if (solver_.GetIter()==0 && LvlSet_.GetSolver().GetResid()<LvlSet_.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cout << "===> frac.step A: Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();

    // ------ frac. step B ------
    InitStep( false);
    for (int i=0; i<maxFPiter; ++i)
    {
        DoNonlinearFPIter();
        if (iter_nonlinear_==1 && LvlSet_.GetSolver().GetResid()<LvlSet_.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cout << "===> frac.step B: Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();

    // ------ frac. step C ------
    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoStokesFPIter();
        if (solver_.GetIter()==0 && LvlSet_.GetSolver().GetResid()<LvlSet_.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cout << "===> frac.step C: Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();
}

template <class StokesT, class SolverT>
void OperatorSplitting2PhaseCL<StokesT,SolverT>::Update()
{
    MLIdxDescCL* const vidx= &Stokes_.vel_idx;
    TimerCL time;
    time.Reset();
    time.Start();

    std::cout << "Updating discretization...\n";
    mat_->clear();
    AN_.clear();
    Stokes_.ClearMat();
    LvlSet_.ClearMat();
    mat_->clear();
    // IndexDesc setzen
    cplM_->SetIdx( vidx);    old_cplM_->SetIdx( vidx);
    cplA_->SetIdx( vidx);    old_cplA_->SetIdx( vidx);
    cplN_->SetIdx( vidx);    old_cplN_->SetIdx( vidx);
    curv_->SetIdx( vidx);
    rhs_.resize( vidx->NumUnknowns());
    ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
    Stokes_.SetIdx();

    // Diskretisierung
    LvlSet_.SetupSystem( Stokes_.GetVelSolution() );
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, b_, old_cplA_, old_cplM_, LvlSet_, Stokes_.t);
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass( &Stokes_.prM, LvlSet_);

    time.Stop();
    std::cout << "Discretizing took " << time.GetTime() << " sec.\n";
}


// ==============================================
//              RecThetaScheme2PhaseCL
// ==============================================

template <class StokesT, class SolverT, class RelaxationPolicyT>
RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::RecThetaScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double theta, double nonlinear, bool withProjection, double stab, bool trapezoid)
  : base_( Stokes, ls, theta, nonlinear),
    solver_( solver), withProj_( withProjection), stab_( stab), Mold_( 0), trapezoid_( trapezoid && theta == 0.5),
#ifndef _PAR
    ssorpc_(), Msolver_( ssorpc_, 200, 1e-10, true),
    ispc_( &Stokes_.B.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), Stokes_.pr_idx.GetFinest(), 1.0, 0.0, 1e-4, 1e-4),
    Ssolver_( ispc_, 200, 200, 1e-10, true)
#else
    MsolverPC_(Stokes.vel_idx.GetFinest()), Msolver_(200, 1e-10, Stokes.vel_idx.GetFinest(), MsolverPC_, false, true),
    SsolverPC_(Stokes.B.Data.GetFinestPtr(), Stokes.prM.Data.GetFinestPtr(), Stokes.M.Data.GetFinestPtr(),
               Stokes.pr_idx.GetFinest(), Stokes.vel_idx.GetFinest(), 1.0, 0.0, 1e-4, 1e-4),
               Ssolver_(100, 200, 1e-10, Stokes.pr_idx.GetFinest(), SsolverPC_, true)
#endif
{
    Update();
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::~RecThetaScheme2PhaseCL()
{
    delete Mold_;
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::MaybeStabilize (VectorCL& b)
{
    if (stab_ == 0.0) return;

    cplLB_.SetIdx( &Stokes_.vel_idx);
    LB_.SetIdx( &Stokes_.vel_idx, &Stokes_.vel_idx);
    // The MatrixBuilderCL's method of determining when to reuse the pattern
    // is not save for matrix LB_
    LB_.Data.clear();
    Stokes_.SetupLB( &LB_, &cplLB_, LvlSet_, Stokes_.t);

    MLMatrixCL mat0( *mat_);
    mat_->clear();
    const double s= stab_*theta_*dt_;
    std::cout << "Stabilizing with: " << s << '\n';
    mat_->LinComb( 1., mat0, s, LB_.Data);
    b+= s*(LB_.Data*Stokes_.v.Data);
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::InitStep()
// compute all terms that don't change during the following FP iterations
{
    LvlSet_.ComputeRhs( ls_rhs_);

    std::cout << "InitStep-dt_: " << dt_ << std::endl;
    if (trapezoid_) {
        rhs_ = (1./dt_) * (Stokes_.M.Data * Stokes_.v.Data) + vdot_;
    }
    else {
        rhs_= (1./dt_)*Stokes_.v.Data;
        if (theta_ != 1.)
            rhs_+= (1. - theta_)*vdot_;
    }
    Stokes_.t+= dt_;

    if (theta_ != 0. && theta_ != 1. && (!trapezoid_))
        Stokes_.p.Data*= theta_; // Just to have a better starting-value for p.
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::DoProjectionStep( const VectorCL& /*rhscurv*/)
// perform preceding projection step
{
    std::cout << "~~~~~~~~~~~~~~~~ NO Projection step\n";
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::EvalLsetNavStokesEquations()
// perform fixed point iteration
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    time.Start();

    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution());
    LvlSet_.SetTimeStep( dt_);

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing Levelset took " << duration << " sec.\n";

    time.Reset();

    LvlSet_.DoStep( ls_rhs_);

    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving Levelset took " << duration << " sec.\n";

    time.Reset();
    time.Start();

    curv_->Clear();
    LvlSet_.AccumulateBndIntegral( *curv_);

    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, b_, b_, cplM_, LvlSet_, Stokes_.t);
    if (Stokes_.UsesXFEM()) {
        Stokes_.UpdateXNumbering( &Stokes_.pr_idx, LvlSet_);
        Stokes_.UpdatePressure( &Stokes_.p);
        Stokes_.c.SetIdx( &Stokes_.pr_idx);
        Stokes_.B.SetIdx( &Stokes_.pr_idx, &Stokes_.vel_idx);
        Stokes_.prA.SetIdx( &Stokes_.pr_idx, &Stokes_.pr_idx);
        Stokes_.prM.SetIdx( &Stokes_.pr_idx, &Stokes_.pr_idx);
        // The MatrixBuilderCL's method of determining when to reuse the pattern
        // is not save for P1X-elements.
        Stokes_.B.Data.clear();
        Stokes_.prA.Data.clear();
        Stokes_.prM.Data.clear();
        Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    }
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass ( &Stokes_.prM, LvlSet_);
    VectorCL b2( curv_->Data + b_->Data);
    double alpha = nonlinear_;
    if (trapezoid_) {
        mat_->LinComb( 1./dt_, Stokes_.M.Data, 1./dt_, *Mold_, 1.0, Stokes_.A.Data);
        b2 += rhs_ + (1./dt_) * (Stokes_.M.Data * oldv_); /* only if time-dep DirBC:+ (1./dt_)*cplM_->Data + coupling of M with vdot_new*/
    }
    else {
        mat_->LinComb( 1./dt_, Stokes_.M.Data, theta_, Stokes_.A.Data);
        b2 *= theta_;
        b2 += Stokes_.M.Data*rhs_; /* only if time-dep DirBC:+ (1./dt_)*cplM_->Data + coupling of M with vdot_new*/
        alpha *= theta_;
    }
    MaybeStabilize( b2);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing NavierStokes/Curv took "<<duration<<" sec.\n";

    time.Reset();
    solver_.Solve( *mat_, Stokes_.B.Data,
        Stokes_.v, Stokes_.p.Data,
        b2, *cplN_, Stokes_.c.Data, alpha);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving NavierStokes: residual: " << solver_.GetResid()
              << "\titerations: " << solver_.GetIter()
              << "\ttime: " << duration << "s\n";
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::CommitStep()
{
    std::swap( b_, old_b_);
    std::swap( cplM_, old_cplM_);
    std::swap( cplN_, old_cplN_);
    std::swap( curv_, old_curv_);

    if (theta_ != 0.) {
        VectorCL vdot1( (1./dt_)*(Stokes_.v.Data - oldv_));
        if (trapezoid_) {
            vdot_ = Stokes_.M.Data * vdot1 + (*Mold_) * vdot1 - vdot_;
            delete Mold_;
            Mold_ = new MLMatrixCL( Stokes_.M.Data);
        }
        else {
            if (theta_ != 1.) {
                vdot_= vdot1 - (1. - theta_)*vdot_;

                Stokes_.p.Data*= 1./theta_;
                vdot_*= 1./theta_;
            }
        }
    }
    else {
        ComputePressure();
        ComputeVelocityDot();
    }
    oldv_= Stokes_.v.Data;

//     static int mycount( 1);
//     if ( mycount++ % 20 == 1) {
//         VectorCL pp( Stokes_.p.Data);
//         ComputePressure();
//         std::cout << "pressure difference: " << norm( Stokes_.p.Data - pp)/norm( pp) << '\n';
//         Stokes_.p.Data= pp;
//
//         if (theta_ != 1.) { // Implicit Euler does not need and calculate vdot_.
//             VectorCL vd( vdot_);
//             ComputeVelocityDot();
//             std::cout << "vdot difference: " << norm( vdot_ - vd)/norm( vd) << '\n';
//             vdot_= vd;
//         }
//     }
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 99;

    InitStep();
    RelaxationPolicyT relax;
    for (int i=0; i<maxFPiter; ++i)
    {
        std::cout << "~~~~~~~~~~~~~~~~ FP-Iter " << i+1 << '\n';
        const VectorCL v( Stokes_.v.Data);
        const VectorCL phi( LvlSet_.Phi.Data);
        EvalLsetNavStokesEquations();
        if (solver_.GetIter()==0 && LvlSet_.GetSolver().GetResid()<LvlSet_.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cout << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
        Stokes_.v.Data   = v   - Stokes_.v.Data;
        LvlSet_.Phi.Data = phi - LvlSet_.Phi.Data;

        // quasi newton method: relax computes the update vector
        relax.Update( Stokes_.v, LvlSet_.Phi);

        Stokes_.v.Data   = v   - Stokes_.v.Data;
        LvlSet_.Phi.Data = phi - LvlSet_.Phi.Data;
    }
    CommitStep();
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::Update()
{
    MLIdxDescCL* const vidx= &Stokes_.vel_idx;
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    time.Start();

    std::cout << "Updating discretization...\n";
    if (withProj_ || stab_!=0.0)
    {
        LB_.Data.clear();
        LB_.SetIdx( vidx, vidx);
        cplLB_.SetIdx( vidx);
    }
    Stokes_.ClearMat();
    LvlSet_.ClearMat();
    mat_->clear();
    // IndexDesc setzen
    b_->SetIdx( vidx);       old_b_->SetIdx( vidx);
    cplM_->SetIdx( vidx);    old_cplM_->SetIdx( vidx);
    cplN_->SetIdx( vidx);    old_cplN_->SetIdx( vidx);
    curv_->SetIdx( vidx);    old_curv_->SetIdx( vidx);
    rhs_.resize( vidx->NumUnknowns());
    ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
    vdot_.resize( vidx->NumUnknowns());
    oldv_.resize( vidx->NumUnknowns());
    oldv_= Stokes_.v.Data;
    Stokes_.SetIdx();

    cplLB_.SetIdx( vidx);
    LB_.Data.clear();
    LB_.SetIdx( vidx, vidx);
    Stokes_.SetupLB( &LB_, &cplLB_, LvlSet_, Stokes_.t);

    // Diskretisierung
    LvlSet_.AccumulateBndIntegral( *old_curv_);
    LvlSet_.SetupSystem( Stokes_.GetVelSolution() );
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, old_b_, old_b_, old_cplM_, LvlSet_, Stokes_.t);
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);

    // Vorkonditionierer
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass( &Stokes_.prM, LvlSet_);

    if (trapezoid_) {
        ComputeVelocityDot();
    }
    else {
        // initialer Druck
        if (theta_ != 1.) {
            ComputePressure();
            ComputeVelocityDot();
        }
    }

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}


template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::ComputePressure ()
{
    VectorCL b2( old_b_->Data + old_curv_->Data
        - Stokes_.A.Data*Stokes_.v.Data
        + nonlinear_*(old_cplN_->Data - Stokes_.N.Data*Stokes_.v.Data));
    VectorCL b3( b2.size());

#ifndef _PAR
    SchurComplMatrixCL<PCG_SsorCL, MLMatrixCL> S( Msolver_, Stokes_.M.Data, Stokes_.B.Data);
#else
    ParSchurComplMatrixCL<MsolverT, MLMatrixCL, ExchangeCL> S(Msolver_, Stokes_.M.Data, Stokes_.B.Data, Stokes_.vel_idx.GetEx());
#endif
    Msolver_.Solve( Stokes_.M.Data, b3, b2);
    std::cout << "ComputePressure: rhs: iter= " << Msolver_.GetIter() << "\tres= " << Msolver_.GetResid() << '\n';

    Msolver_.SetTol( 1e-13);

    VectorCL b4( Stokes_.B.Data*b3);
    if (Stokes_.UsesXFEM()) {
        VecDescCL Bdotv( &Stokes_.pr_idx);
        Stokes_.SetupBdotv( &Bdotv, &Stokes_.v, LvlSet_, Stokes_.t);
        b4+= Bdotv.Data;
    }
    Ssolver_.Solve( S, Stokes_.p.Data, b4);
    std::cout << "ComputePressure: pressure: iter= " << Ssolver_.GetIter() << "\tres= " << Ssolver_.GetResid() << '\n';
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::ComputeVelocityDot ()
{
    if (trapezoid_) {
        delete Mold_;
        Mold_ = new MLMatrixCL( Stokes_.M.Data);
        vdot_ = (-1.0)*( Stokes_.A.Data * Stokes_.v.Data ) + old_curv_->Data + old_b_->Data - transp_mul( Stokes_.B.Data, Stokes_.p.Data );
        return;
    }
    VectorCL b2( old_b_->Data + old_curv_->Data
        - Stokes_.A.Data*Stokes_.v.Data
        + nonlinear_*(old_cplN_->Data - Stokes_.N.Data*Stokes_.v.Data)
        - transp_mul( Stokes_.B.Data, Stokes_.p.Data));
    Msolver_.SetTol( 1e-10);
    Msolver_.Solve( Stokes_.M.Data, vdot_, b2);
    std::cout << "ComputeVelocityDot: vdot: iter= " << Msolver_.GetIter() << "\tres= " << Msolver_.GetResid() << '\n';
}


// ==============================================
//            CrankNicolsonScheme2PhaseCL
// ==============================================

template <class StokesT, class SolverT, class RelaxationPolicyT>
CrankNicolsonScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::CrankNicolsonScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double nonlinear, bool withProjection, double stab)
  : base_( Stokes, ls, solver, 1.0, nonlinear, withProjection, stab), step_(1)
{}

template <class StokesT, class SolverT, class RelaxationPolicyT>
CrankNicolsonScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::~CrankNicolsonScheme2PhaseCL()
{}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void CrankNicolsonScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::DoStep(int maxFPIter)
{
    switch (step_)
    {
        case 1 :
        {
            std::cout << "~~~~~~~~~~~~~~~~ initial backward euler step\n";
            tmpdt_=dt_;
            base_::SetTimeStep(0.2*tmpdt_*tmpdt_, 1.0);
            base_::DoStep(maxFPIter);
            base_::ComputeVelocityDot();
            base_::SetTimeStep((1.0-0.2*tmpdt_)*tmpdt_, 0.5);
            step_++;
        } break;
        case 2 :
        {
            base_::SetTimeStep(tmpdt_, 0.5);
            step_++;
        } break;
    }
    base_::DoStep(maxFPIter);
}

template <class StokesT, class SolverT, class RelaxationPolicyT>
void CrankNicolsonScheme2PhaseCL<StokesT,SolverT,RelaxationPolicyT>::Update()
{
    base_::SetTimeStep(tmpdt_, 1.0);
    base_::Update();
    step_= 1;
}

inline void cplDeltaSquaredPolicyCL::Update( VecDescCL& v, VecDescCL& phi)
{
    if (firststep_) {
        const size_t vsize = v.Data.size();
        const size_t phisize = phi.Data.size();
        v_old_.resize   ( vsize);
        phi_old_.resize ( phisize);
        v_diff_.resize  ( vsize);
        phi_diff_.resize( phisize);
        v_old_= v.Data; phi_old_= phi.Data;
        firststep_ = false;
        if (output_)
            (*output_) << "omega: " << omega_ << std::endl;
        return;
    }
    v_diff_=  v.Data - v_old_; phi_diff_= phi.Data - phi_old_;
#ifndef _PAR
    omega_*= -(dot( v_diff_, v_old_) + dot( phi_diff_, phi_old_))
            / (norm_sq( v_diff_) + norm_sq( phi_diff_));
#else
    const bool useAccur=true;
    ExchangeCL& ExVel  = v.RowIdx->GetEx();
    ExchangeCL& ExLset = phi.RowIdx->GetEx();
    omega_*=-(ExVel.ParDot( v_diff_, true, v_old_, true, useAccur)
            + ExLset.ParDot( phi_diff_, true, phi_old_, true, useAccur))
            / (ExVel.Norm_sq( v_diff_, true, useAccur) + ExLset.Norm_sq( phi_diff_, true, useAccur));
#endif
    if (output_)
        (*output_) << "omega: " << omega_ << std::endl;
    v_old_= v.Data; phi_old_= phi.Data;
    v.Data   *= omega_;
    phi.Data *= omega_;
}

inline void cplBroydenPolicyCL::Update( VecDescCL& v, VecDescCL& phi)
{
    F1_.push_back( v.Data);
    F2_.push_back( phi.Data);
#ifndef _PAR
    sigma_ = norm_sq( v.Data) + norm_sq( phi.Data);
#else
    const bool useAccur=true;
    ExchangeCL& ExVel  = v.RowIdx->GetEx();
    ExchangeCL& ExLset = phi.RowIdx->GetEx();
    sigma_ = ExVel.Norm_sq( v.Data, true, useAccur) + ExLset.Norm_sq( phi.Data, true, useAccur);
#endif

    if (sigma_ < tol_*tol_) {
        if (output_)
            (*output_) << "Solution found" << std::endl;
        //return;
    }

    if (firststep_) {
        firststep_ = false;
        sigma0_ = sigma_;
        return;
    }

    // at this point: F1_.size() >= 2
    const size_t pos = F1_.size() - 2;
    deltaF1_.push_back( VectorCL( F1_.back() - F1_[pos]));
    deltaF2_.push_back( VectorCL( F2_.back() - F2_[pos]));

    const double theta = std::sqrt(sigma_/sigma0_);
    if (theta >= thetamax_) {
        if (output_)
            (*output_) << "No convergence: theta = " << theta << std::endl;
        //return;
    }
    sigma0_ = sigma_;

    VectorCL w1( deltaF1_.back());
    VectorCL w2( deltaF2_.back());
#ifndef _PAR
    gamma_.push_back( norm_sq( w1) + norm_sq( w2));
#else
    gamma_.push_back( ExVel.Norm_sq( w1, true, useAccur) + ExLset.Norm_sq( w2, true, useAccur));
#endif

    kappa_ /= (1.0-2.0*theta);

    if (kappa_ >= kappamax_){
        if (output_)
            (*output_) << "ill-conditioned update: kappa = "<< kappa_ << std::endl;
        //return;
    }

    VectorCL v1( F1_.back()), v2( F2_.back());
#ifndef _PAR
    const double factor = 1.0 - ( dot( w1, v1) + dot( w2, v2)) / gamma_.back();
#else
    const double factor = 1.0 - ( ExVel.ParDot( w1, true, v1, true, useAccur) + ExLset.ParDot( w2, true, v2, true, useAccur)) / gamma_.back();
#endif
    v1*= factor;
    v2*= factor;

    for (int j=deltaF1_.size()-1; j>=0; --j) {
#ifndef _PAR
        const double beta = ( dot (deltaF1_[j], v1) + dot( deltaF2_[j], v2)) / gamma_[j];
#else
        const double beta = ( ExVel.ParDot (deltaF1_[j], true, v1, true, useAccur) + ExLset.ParDot( deltaF2_[j], true, v2, true, useAccur)) / gamma_[j];
#endif
        v1 -= beta*F1_[j+1];
        v2 -= beta*F2_[j+1];
    }
    v.Data   = v1;
    phi.Data = v2;
}


} // end of namespace DROPS
