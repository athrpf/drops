/// \file coupling.tpp
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

#include "num/nssolver.h"

namespace DROPS
{


// ==============================================
//           LinThetaScheme2PhaseCL
// ==============================================

template <class LsetSolverT>
LinThetaScheme2PhaseCL<LsetSolverT>::LinThetaScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
                double stk_theta, double ls_theta, double nonlinear, bool implicitCurv)
  : base_( Stokes, ls, lsetmod, nonlinear), solver_( solver), lsetsolver_( lsetsolver),
  stk_theta_( stk_theta), ls_theta_( ls_theta), cplA_(new VelVecDescCL), implCurv_( implicitCurv)
{
    Update();
}

template <class LsetSolverT>
LinThetaScheme2PhaseCL<LsetSolverT>::~LinThetaScheme2PhaseCL()
{
    delete cplA_;
}

template <class LsetSolverT>
void LinThetaScheme2PhaseCL<LsetSolverT>::SolveLsNs()
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
        mat_->LinComb( 1./dt_, Stokes_.M.Data, stk_theta_, Stokes_.A.Data);
    else // semi-implicit treatment of curvature term, cf. Baensch
    {
        MLMatrixCL mat0;
        mat0.LinComb( 1./dt_, Stokes_.M.Data, stk_theta_, Stokes_.A.Data);
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
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    L_->LinComb( 1./dt_, LvlSet_.E, ls_theta_, LvlSet_.H);
    ls_rhs_= (1./dt_)*(LvlSet_.E*LvlSet_.Phi.Data) - (1.-ls_theta_)*(LvlSet_.H*LvlSet_.Phi.Data);

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing Levelset took " << duration << " sec.\n";
    time.Reset();

    lsetsolver_.Solve( *L_, LvlSet_.Phi.Data, ls_rhs_);
    std::cout << "res = " << lsetsolver_.GetResid() << ", iter = " << lsetsolver_.GetIter() <<std::endl;

    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving Levelset took " << duration << " sec.\n";

    lsetmod_.maybeDoVolCorr( LvlSet_);
    lsetmod_.maybeDoReparam( LvlSet_);

    time.Reset();

    Stokes_.t+= dt_;

    // rhs for new level set
    curv_->Clear();
    LvlSet_.AccumulateBndIntegral( *curv_);
    Stokes_.SetupRhs1( b_, LvlSet_, Stokes_.t);

    rhs_=  (1./dt_)*(Stokes_.M.Data*Stokes_.v.Data) + stk_theta_*b_->Data + cplA_->Data
        + (1.-stk_theta_)*(old_b_->Data - Stokes_.A.Data*Stokes_.v.Data - nonlinear_*(Stokes_.N.Data*Stokes_.v.Data - old_cplN_->Data));

    if (!implCurv_)
        rhs_+= stk_theta_*curv_->Data + (1.-stk_theta_)*old_curv_->Data;
    else // semi-implicit treatment of curvature term, cf. Baensch
        rhs_+= old_curv_->Data + dt_*cplLB_.Data;

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing Rhs/Curv took "<< duration <<" sec.\n";

    time.Reset();
    solver_.Solve( *mat_, Stokes_.B.Data,
        Stokes_.v, Stokes_.p.Data,
        rhs_, *cplN_, Stokes_.c.Data, /*alpha*/ stk_theta_*nonlinear_);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving NavierStokes: residual: " << solver_.GetResid()
              << "\titerations: " << solver_.GetIter()
              << "\ttime: " << duration << "sec\n";
}

template <class LsetSolverT>
void LinThetaScheme2PhaseCL<LsetSolverT>::CommitStep()
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

template <class LsetSolverT>
void LinThetaScheme2PhaseCL<LsetSolverT>::DoStep( int)
{
    lsetmod_.init();
    SolveLsNs();
    CommitStep();
}

template <class LsetSolverT>
void LinThetaScheme2PhaseCL<LsetSolverT>::Update()
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
    mat_->clear();
    L_->clear();
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
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}

// ==============================================
//              OperatorSplitting2PhaseCL
// ==============================================


template <class LsetSolverT>
OperatorSplitting2PhaseCL<LsetSolverT>::OperatorSplitting2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverBaseCL& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod, int gm_iter, double gm_tol, double nonlinear)

  : base_( Stokes, ls, lsetmod, nonlinear), solver_(solver), lsetsolver_( lsetsolver),
    gm_( pc_, 100, gm_iter, gm_tol, false /*test absolute resid*/),
    stk_theta_( 1.0 - std::sqrt( 2.)/2.), ls_theta_( 1.0 - std::sqrt( 2.)/2.),
    cplA_( new VelVecDescCL), old_cplA_( new VelVecDescCL),
    alpha_( (1.0 - 2.0*stk_theta_)/(1.0 - stk_theta_))
{
    throw DROPSErrCL("OperatorSplitting2PhaseCL: Not correctly implemented, yet");
#ifdef _PAR
    throw DROPSErrCL("OperatorSplitting2PhaseCL: Not parallelized, yet, sorry");
#endif
    std::cout << "theta = " << stk_theta_ << "\talpha = " << alpha_ << std::endl;
    Update();
}

template <class LsetSolverT>
OperatorSplitting2PhaseCL<LsetSolverT>::~OperatorSplitting2PhaseCL()
{
    delete cplA_; delete old_cplA_;
}

template <class LsetSolverT>
void OperatorSplitting2PhaseCL<LsetSolverT>::InitStep( bool StokesStep)
{
// compute all terms that don't change during the following FP iterations

    const double fracdt_= StokesStep ? stk_theta_*dt_ : (1-2*stk_theta_)*dt_;

    ls_theta_ = StokesStep ? alpha_ : 1-alpha_;
    L_->LinComb( 1./fracdt_, LvlSet_.E, ls_theta_, LvlSet_.H);
    ls_rhs_= (1./fracdt_)*LvlSet_.Phi.Data;
    if (ls_theta_ != 1.) {
        VectorCL tmp( ls_rhs_.size());
        LsetSolverT gm( lsetsolver_);
        gm.Solve( LvlSet_.E, tmp, (const VectorCL)( LvlSet_.H*LvlSet_.Phi.Data));
        std::cout << "ComputeRhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
        ls_rhs_-= (1. - ls_theta_)*tmp;
    }

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

template <class LsetSolverT>
void OperatorSplitting2PhaseCL<LsetSolverT>::DoStokesFPIter()
// perform fixed point iteration: Levelset / Stokes
// for fractional steps A and C
{
    TimerCL time;
    time.Reset();
    time.Start();
    const double fracdt_= stk_theta_*dt_;
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    L_->LinComb( 1./fracdt_, LvlSet_.E, alpha_, LvlSet_.H);
    time.Stop();
    std::cout << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    lsetsolver_.Solve( *L_, LvlSet_.Phi.Data, VectorCL( LvlSet_.E*ls_rhs_));
    std::cout << "res = " << lsetsolver_.GetResid() << ", iter = " << lsetsolver_.GetIter() <<std::endl;
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
    else
        Stokes_.SetupRhs2( &Stokes_.c, LvlSet_, Stokes_.t);
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

template <class LsetSolverT>
void OperatorSplitting2PhaseCL<LsetSolverT>::DoNonlinearFPIter()
// perform fixed point iteration: Levelset / nonlinear system
// for fractional step B
{
    const double fracdt_= dt_*(1.-2.*stk_theta_);
    TimerCL time;
    time.Reset();
    time.Start();
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    L_->LinComb( 1./fracdt_, LvlSet_.E, 1.-alpha_, LvlSet_.H);
    time.Stop();
    std::cout << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();
    lsetsolver_.Solve( *L_, LvlSet_.Phi.Data, VectorCL( LvlSet_.E*ls_rhs_));
    std::cout << "res = " << lsetsolver_.GetResid() << ", iter = " << lsetsolver_.GetIter() <<std::endl;
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

template <class LsetSolverT>
void OperatorSplitting2PhaseCL<LsetSolverT>::CommitStep()
{
    std::swap( cplM_, old_cplM_);
    std::swap( cplA_, old_cplA_);
    std::swap( cplN_, old_cplN_);
}

template <class LsetSolverT>
void OperatorSplitting2PhaseCL<LsetSolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 99;

    // ------ frac. step A ------
    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoStokesFPIter();
        if (solver_.GetIter()==0 && lsetsolver_.GetResid()<lsetsolver_.GetTol()) // no change of vel -> no change of Phi
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
        if (iter_nonlinear_==1 && lsetsolver_.GetResid()<lsetsolver_.GetTol()) // no change of vel -> no change of Phi
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
        if (solver_.GetIter()==0 && lsetsolver_.GetResid()<lsetsolver_.GetTol()) // no change of vel -> no change of Phi
        {
            std::cout << "===> frac.step C: Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();
}

template <class LsetSolverT>
void OperatorSplitting2PhaseCL<LsetSolverT>::Update()
{
    MLIdxDescCL* const vidx= &Stokes_.vel_idx;
    TimerCL time;
    time.Reset();
    time.Start();

    std::cout << "Updating discretization...\n";
    AN_.clear();
    Stokes_.ClearMat();
    LvlSet_.ClearMat();
    mat_->clear();
    L_->clear();
    // IndexDesc setzen
    cplM_->SetIdx( vidx);    old_cplM_->SetIdx( vidx);
    cplA_->SetIdx( vidx);    old_cplA_->SetIdx( vidx);
    cplN_->SetIdx( vidx);    old_cplN_->SetIdx( vidx);
    curv_->SetIdx( vidx);
    rhs_.resize( vidx->NumUnknowns());
    ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
    Stokes_.SetIdx();

    // Diskretisierung
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, b_, old_cplA_, old_cplM_, LvlSet_, Stokes_.t);
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass( &Stokes_.prM, LvlSet_);

    time.Stop();
    std::cout << "Discretizing took " << time.GetTime() << " sec.\n";
}


// ==============================================
//            CoupledTimeDisc2PhaseBaseCL
// ==============================================

template <class LsetSolverT, class RelaxationPolicyT>
CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::CoupledTimeDisc2PhaseBaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod, double tol,
            double nonlinear, bool withProjection, double stab)
  : base_( Stokes, ls, lsetmod, nonlinear),
    solver_( solver), lsetsolver_( lsetsolver), tol_(tol), withProj_( withProjection), stab_( stab), alpha_( nonlinear_)
{
    Update();
}

template <class LsetSolverT, class RelaxationPolicyT>
CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::~CoupledTimeDisc2PhaseBaseCL()
{}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::MaybeStabilize (VectorCL& b)
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
    const double s= stab_*dt_;
    std::cout << "Stabilizing with: " << s << '\n';
    mat_->LinComb( 1., mat0, s, LB_.Data);
    b+= s*(LB_.Data*Stokes_.v.Data);
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::InitStep()
// compute all terms that don't change during the following FP iterations
{
    lsetmod_.init();
    dphi_ = 0;
    std::cout << "InitStep-dt_: " << dt_ << std::endl;
    Stokes_.t+= dt_;
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::DoProjectionStep( const VectorCL& /*rhscurv*/)
// perform preceding projection step
{
    std::cout << "~~~~~~~~~~~~~~~~ NO Projection step\n";
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::EvalLsetNavStokesEquations()
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
    SetupLevelsetSystem();
    LvlSet_.Phi.Data -= dphi_;
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing Levelset took " << duration << " sec.\n";

    time.Reset();

    lsetsolver_.Solve( *L_, LvlSet_.Phi.Data, ls_rhs_);
    std::cout << "res = " << lsetsolver_.GetResid() << ", iter = " << lsetsolver_.GetIter() << std::endl;

    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving Levelset took " << duration << " sec.\n";

    dphi_ = lsetmod_.maybeDoVolCorr( LvlSet_);

    time.Reset();
    time.Start();

    SetupNavStokesSystem();

    MaybeStabilize( rhs_);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing NavierStokes/Curv took "<<duration<<" sec.\n";

    time.Reset();

    solver_.Solve( *mat_, Stokes_.B.Data,
        Stokes_.v, Stokes_.p.Data,
        rhs_, *cplN_, Stokes_.c.Data, alpha_);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Solving NavierStokes: residual: " << solver_.GetResid()
              << "\titerations: " << solver_.GetIter()
              << "\ttime: " << duration << "s\n";
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::SetupStokesMatVec()
/// setup matrices A, M, B, prA, prM and vectors b+cplA, cplM, curv, c
{
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
    else
        Stokes_.SetupRhs2( &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass ( &Stokes_.prM, LvlSet_);
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::CommitStep()
{
    std::swap( b_, old_b_);
    std::swap( cplM_, old_cplM_);
    std::swap( cplN_, old_cplN_);
    std::swap( curv_, old_curv_);
    lsetmod_.maybeDoReparam( LvlSet_);
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 99;

    InitStep();
    RelaxationPolicyT relax;
#ifdef _PAR
    const bool useAccur=true;
    ExchangeCL& ExVel  = Stokes_.v.RowIdx->GetEx();
#endif
    double res_u = 0.0;
    for (int i=0; i<maxFPiter; ++i)
    {
        std::cout << "~~~~~~~~~~~~~~~~ FP-Iter " << i+1 << '\n';
        const VectorCL v( Stokes_.v.Data);
        EvalLsetNavStokesEquations();
        if (solver_.GetIter()==0 && lsetsolver_.GetResid()<lsetsolver_.GetTol()) // no change of vel -> no change of Phi
        {
            std::cout << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
        Stokes_.v.Data = v - Stokes_.v.Data;

        // quasi newton method: relax computes the update vector
        relax.Update( Stokes_.v);

#ifndef _PAR
        res_u   = std::sqrt( dot( Stokes_.v.Data, Stokes_.M.Data*Stokes_.v.Data));
#else
        res_u   = std::sqrt( ExVel.ParDot( Stokes_.v.Data, true, Stokes_.M.Data*Stokes_.v.Data, true, useAccur));
#endif

        std::cout << "residual of u: " << res_u << std::endl;

        Stokes_.v.Data = v - Stokes_.v.Data;

        if (res_u < tol_) {
            std::cout << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();
}

template <class LsetSolverT, class RelaxationPolicyT>
void CoupledTimeDisc2PhaseBaseCL<LsetSolverT,RelaxationPolicyT>::Update()
{
    MLIdxDescCL* const vidx= &Stokes_.vel_idx;

    Stokes_.ClearMat();
    LvlSet_.ClearMat();
    mat_->clear();
    L_->clear();

    // IndexDesc setzen
    b_->SetIdx( vidx);       old_b_->SetIdx( vidx);
    cplM_->SetIdx( vidx);    old_cplM_->SetIdx( vidx);
    cplN_->SetIdx( vidx);    old_cplN_->SetIdx( vidx);
    curv_->SetIdx( vidx);    old_curv_->SetIdx( vidx);

    rhs_.resize   ( vidx->NumUnknowns());
    ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
}

// ==============================================
//            MidPointTimeDisc2PhaseCL
// ==============================================

template <class LsetSolverT, class RelaxationPolicyT>
MidPointTimeDisc2PhaseCL<LsetSolverT,RelaxationPolicyT>::MidPointTimeDisc2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
            double tol, double nonlinear, bool withProjection, double stab, bool implicitpressure)
  : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, nonlinear, withProjection, stab),
    implicitpressure_( implicitpressure)
{
    if (nonlinear != 0.0)
        throw DROPSErrCL("MidPointTimeDisc2PhaseCL: Not yet implemented for Navier-Stokes equations\n");
    if (Stokes_.UsesXFEM() && !implicitpressure_) {
        std::cerr << "MidPointTimeDisc2PhaseCL: XFEM for pressure detected. This is not implemented, yet, using fully implicit pressure" << std::endl;
        implicitpressure_ = true;
    }
    Update();
}

template <class LsetSolverT, class RelaxationPolicyT>
void MidPointTimeDisc2PhaseCL<LsetSolverT,RelaxationPolicyT>::InitStep()
{
    base_::InitStep();
    if (implicitpressure_) {       // Just to have a better starting-value for p.
        Stokes_.p.Data *= 2.0;
    }
}

template <class LsetSolverT, class RelaxationPolicyT>
void MidPointTimeDisc2PhaseCL<LsetSolverT,RelaxationPolicyT>::CommitStep()
{
    base_::CommitStep();
    oldphi_ = LvlSet_.Phi.Data;
    if (implicitpressure_) {
        Stokes_.p.Data *= 0.5;
    }
    else {
        oldp_= Stokes_.p.Data;
    }
    oldv_= Stokes_.v.Data;
}

template <class LsetSolverT, class RelaxationPolicyT>
void MidPointTimeDisc2PhaseCL<LsetSolverT,RelaxationPolicyT>::Update()
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    time.Start();

    base_::Update();
    oldphi_.resize( LvlSet_.idx.NumUnknowns());
    oldphi_ = LvlSet_.Phi.Data;

    oldv_.resize( Stokes_.v.Data.size());
    oldv_= Stokes_.v.Data;

    oldp_.resize( Stokes_.p.Data.size());
    oldp_= Stokes_.p.Data;

    Stokes_.SetIdx();

    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}

template <class LsetSolverT, class RelaxationPolicyT>
void MidPointTimeDisc2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupNavStokesSystem()
{
    const VectorCL tmpphi( LvlSet_.Phi.Data);
    LvlSet_.Phi.Data += oldphi_;
    LvlSet_.Phi.Data *= 0.5;

    base_::SetupStokesMatVec(); // setup all matrices (except N) and rhs

    LvlSet_.Phi.Data = tmpphi;

    alpha_ = nonlinear_;

    mat_->LinComb( 2./dt_, Stokes_.M.Data, 1.0, Stokes_.A.Data);
    rhs_ = 2.0/dt_ * VectorCL(Stokes_.M.Data * oldv_) - Stokes_.A.Data * oldv_ + 2.0* curv_->Data + 2.0* Stokes_.b.Data;
    if (!implicitpressure_)
        rhs_ -= transp_mul(Stokes_.B.Data, oldp_);

}

template <class LsetSolverT, class RelaxationPolicyT>
void MidPointTimeDisc2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupLevelsetSystem()
{
    const VectorCL tmpv( Stokes_.v.Data);
    Stokes_.v.Data = 0.5 * VectorCL( Stokes_.v.Data + oldv_);
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    Stokes_.v.Data = tmpv;

    L_->LinComb( 2./dt_, LvlSet_.E, 1.0, LvlSet_.H);
    ls_rhs_ = 2.0/dt_ * VectorCL(LvlSet_.E * oldphi_) - LvlSet_.H * oldphi_;
}

// ==============================================
//            SpaceTimeDiscTheta2PhaseCL
// ==============================================

template <class LsetSolverT, class RelaxationPolicyT>
SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::SpaceTimeDiscTheta2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
            double tol, double stk_theta, double ls_theta, double nonlinear, bool withProjection, double stab, bool implicitpressure)
  : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, nonlinear, withProjection, stab),  stk_theta_( stk_theta), ls_theta_( ls_theta),
    implicitpressure_( implicitpressure), Mold_( 0), Eold_( 0)
{
    stab_ *= stk_theta_;
    Update();
}

template <class LsetSolverT, class RelaxationPolicyT>
SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::~SpaceTimeDiscTheta2PhaseCL()
{
    delete Eold_; delete Mold_;
}

template <class LsetSolverT, class RelaxationPolicyT>
void SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::InitStep()
{
    base_::InitStep();
    fixed_ls_rhs_ = (1./dt_) * (LvlSet_.E * oldphi_)    + phidot_;
    fixed_rhs_    = (1./dt_) * (Stokes_.M.Data * oldv_) + vdot_;
    if (!implicitpressure_)              // Just to have a better starting-value for p.
        Stokes_.p.Data *= stk_theta_;
}

template <class LsetSolverT, class RelaxationPolicyT>
void SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::CommitStep()
{
    base_::CommitStep();
    VectorCL vdot1( (1./dt_)*(Stokes_.v.Data - oldv_));
    vdot_ = stk_theta_ * (Stokes_.M.Data * vdot1) + (1.0 - stk_theta_) * ((*Mold_) * vdot1) - (1.0 - stk_theta_) * vdot_;
    delete Mold_;
    Mold_ = new MLMatrixCL( Stokes_.M.Data);

    phidot_ = (-1.0) * (LvlSet_.H * LvlSet_.Phi.Data);
    delete Eold_;
    Eold_ = new MatrixCL (LvlSet_.E);

    oldphi_ = LvlSet_.Phi.Data;
    oldv_= Stokes_.v.Data;
    if (implicitpressure_)
        vdot_ += transp_mul( Stokes_.B.Data, Stokes_.p.Data);
    else
        Stokes_.p.Data /= stk_theta_;
    vdot_ /= stk_theta_;
}

template <class LsetSolverT, class RelaxationPolicyT>
void SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::Update()
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    time.Start();

    base_::Update();

    oldphi_.resize( LvlSet_.idx.NumUnknowns());
    oldphi_ = LvlSet_.Phi.Data;
    fixed_ls_rhs_.resize( LvlSet_.idx.NumUnknowns());

    phidot_.resize( LvlSet_.idx.NumUnknowns());

    oldv_.resize( Stokes_.v.Data.size());
    oldv_= Stokes_.v.Data;
    fixed_rhs_.resize( Stokes_.v.Data.size());

    vdot_.resize( Stokes_.v.Data.size());

    Stokes_.SetIdx();

    // Diskretisierung
    LvlSet_.AccumulateBndIntegral( *old_curv_);
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, old_b_, old_b_, old_cplM_, LvlSet_, Stokes_.t);
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);

    // Vorkonditionierer
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass( &Stokes_.prM, LvlSet_);

    ComputeDots();
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}

template <class LsetSolverT, class RelaxationPolicyT>
void SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::ComputeDots ()
{
    delete Mold_;
    Mold_ = new MLMatrixCL( Stokes_.M.Data);
    vdot_ = (-1.0)*( Stokes_.A.Data * Stokes_.v.Data ) + old_curv_->Data + old_b_->Data
            + nonlinear_*(old_cplN_->Data - Stokes_.N.Data*Stokes_.v.Data);
    if (!implicitpressure_)
        vdot_ -= transp_mul( Stokes_.B.Data, Stokes_.p.Data);

    delete Eold_;
    Eold_ = new MatrixCL( LvlSet_.E);
    phidot_ = (-1.0) * (LvlSet_.H * LvlSet_.Phi.Data);
}

template <class LsetSolverT, class RelaxationPolicyT>
void SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupNavStokesSystem()
{
    base_::SetupStokesMatVec(); // setup all matrices (except N) and rhs

    alpha_ = nonlinear_ * stk_theta_;
    mat_->LinComb( stk_theta_/dt_, Stokes_.M.Data, (1.0-stk_theta_)/dt_, *Mold_, stk_theta_, Stokes_.A.Data);
    rhs_ = (1.0 - stk_theta_)*fixed_rhs_ + stk_theta_*( (1./dt_) * (Stokes_.M.Data * oldv_) + curv_->Data + b_->Data); /* TODO time-dep DirBC*/
}

template <class LsetSolverT, class RelaxationPolicyT>
void SpaceTimeDiscTheta2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupLevelsetSystem()
{
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);

    L_->LinComb( ls_theta_/dt_, LvlSet_.E, (1.0-ls_theta_)/dt_, *Eold_, ls_theta_, LvlSet_.H);
    ls_rhs_ = (1.0 - ls_theta_) * fixed_ls_rhs_ + (ls_theta_/dt_) * (LvlSet_.E * oldphi_);
}

// ==============================================
//            EulerBackwardScheme2PhaseCL
// ==============================================

template <class LsetSolverT, class RelaxationPolicyT>
EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::EulerBackwardScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
            double tol, double nonlinear, bool withProjection, double stab)
  : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, nonlinear, withProjection, stab)
{
    Update();
}

template <class LsetSolverT, class RelaxationPolicyT>
EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::~EulerBackwardScheme2PhaseCL()
{}

template <class LsetSolverT, class RelaxationPolicyT>
void EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::InitStep()
// compute all terms that don't change during the following FP iterations
{
    base_::InitStep();
    fixed_ls_rhs_ = (1./dt_)*LvlSet_.Phi.Data;
    fixed_rhs_    = (1./dt_)*Stokes_.v.Data;
}


template <class LsetSolverT, class RelaxationPolicyT>
void EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::CommitStep()
{}

template <class LsetSolverT, class RelaxationPolicyT>
void EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::Update()
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    time.Start();

    std::cout << "Updating discretization...\n";
    base_::Update();
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    fixed_ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
    fixed_rhs_.resize   ( Stokes_.v.Data.size());

    Stokes_.SetIdx();
    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}

template <class LsetSolverT, class RelaxationPolicyT>
void EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupNavStokesSystem()
{
    base_::SetupStokesMatVec(); // setup all matrices (except N) and rhs

    alpha_ = nonlinear_;
    mat_->LinComb( 1./dt_, Stokes_.M.Data, 1.0, Stokes_.A.Data);
    rhs_ = Stokes_.M.Data * fixed_rhs_ + curv_->Data + b_->Data;
}

template <class LsetSolverT, class RelaxationPolicyT>
void EulerBackwardScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupLevelsetSystem()
{
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);

    L_->LinComb( 1./dt_, LvlSet_.E, 1.0, LvlSet_.H);
    ls_rhs_ = LvlSet_.E * fixed_ls_rhs_;
}

// ==============================================
//              RecThetaScheme2PhaseCL
// ==============================================

template <class LsetSolverT, class RelaxationPolicyT>
RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::RecThetaScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, StokesSolverT& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
            double tol, double stk_theta, double ls_theta, double nonlinear, bool withProjection, double stab)
  : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, nonlinear, withProjection, stab),
    stk_theta_( stk_theta), ls_theta_( ls_theta),
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
    stab_ *= stk_theta_;
    Update();
}

template <class LsetSolverT, class RelaxationPolicyT>
RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::~RecThetaScheme2PhaseCL()
{}

template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::InitStep()
// compute all terms that don't change during the following FP iterations
{
    base_::InitStep();
    fixed_ls_rhs_ = (1./dt_)*oldphi_ + ( 1. - ls_theta_ )*phidot_;

    fixed_rhs_= (1./dt_)*oldv_ + ( 1. - stk_theta_)*vdot_;

    Stokes_.p.Data*= stk_theta_; // Just to have a better starting-value for p.
}

template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::CommitStep()
{
    base_::CommitStep();

    vdot_  = ((1./dt_)*(Stokes_.v.Data - oldv_)     - (1. - stk_theta_)*vdot_)  / stk_theta_;

    phidot_= ((1./dt_)*(LvlSet_.Phi.Data - oldphi_) - (1. - ls_theta_)*phidot_) / ls_theta_;

    oldphi_ = LvlSet_.Phi.Data;
    oldv_= Stokes_.v.Data;
    Stokes_.p.Data*= 1./stk_theta_;


//     static int mycount( 1);
//     if ( mycount++ % 20 == 1) {
//         VectorCL pp( Stokes_.p.Data);
//         ComputePressure();
//         std::cout << "pressure difference: " << norm( Stokes_.p.Data - pp)/norm( pp) << '\n';
//         Stokes_.p.Data= pp;
//
//         if (stk_theta_ != 1.) { // Implicit Euler does not need and calculate vdot_.
//             VectorCL vd( vdot_);
//             ComputeDots();
//             std::cout << "vdot difference: " << norm( vdot_ - vd)/norm( vd) << '\n';
//             vdot_= vd;
//         }
//     }
}

template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::Update()
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
    base_::Update();
    phidot_.resize( LvlSet_.idx.NumUnknowns());
    oldphi_.resize( LvlSet_.idx.NumUnknowns());
    oldphi_ = LvlSet_.Phi.Data;

    vdot_.resize( vidx->NumUnknowns());
    oldv_.resize( vidx->NumUnknowns());
    oldv_= Stokes_.v.Data;
    fixed_ls_rhs_.resize( LvlSet_.idx.NumUnknowns());
    fixed_rhs_.resize   ( Stokes_.v.Data.size());
    Stokes_.SetIdx();

    // Diskretisierung
    LvlSet_.AccumulateBndIntegral( *old_curv_);
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);
    Stokes_.SetupSystem1( &Stokes_.A, &Stokes_.M, old_b_, old_b_, old_cplM_, LvlSet_, Stokes_.t);
    Stokes_.SetupSystem2( &Stokes_.B, &Stokes_.c, LvlSet_, Stokes_.t);
    Stokes_.SetupNonlinear( &Stokes_.N, &Stokes_.v, old_cplN_, LvlSet_, Stokes_.t);

    // Vorkonditionierer
    Stokes_.SetupPrStiff( &Stokes_.prA, LvlSet_);
    Stokes_.SetupPrMass( &Stokes_.prM, LvlSet_);

    // initialer Druck
    ComputePressure();
    ComputeDots();

    time.Stop();
    duration=time.GetTime();
    std::cout << "Discretizing took " << duration << " sec.\n";
}


template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::ComputePressure ()
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

template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::ComputeDots ()
{
    VectorCL b2( old_b_->Data + old_curv_->Data
        - Stokes_.A.Data*Stokes_.v.Data
        + nonlinear_*(old_cplN_->Data - Stokes_.N.Data*Stokes_.v.Data)
        - transp_mul( Stokes_.B.Data, Stokes_.p.Data));
    Msolver_.SetTol( 1e-10);
    Msolver_.Solve( Stokes_.M.Data, vdot_, b2);
    std::cout << "ComputeDots: vdot:   iter= " << Msolver_.GetIter() << "\tres= " << Msolver_.GetResid() << std::endl;

    Msolver_.SetTol( 1e-15);
    VectorCL b3 ((-1.0) * (LvlSet_.H * LvlSet_.Phi.Data));
    Msolver_.Solve( LvlSet_.E, phidot_, b3);
    std::cout << "ComputeDots: phidot: iter= " << Msolver_.GetIter() << "\tres= " << Msolver_.GetResid() << std::endl;
}

template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupNavStokesSystem()
{
    base_::SetupStokesMatVec(); // setup all matrices (except N) and rhs

    alpha_ = nonlinear_ * stk_theta_;
    mat_->LinComb( 1./dt_, Stokes_.M.Data, stk_theta_, Stokes_.A.Data);
    rhs_ = Stokes_.M.Data * fixed_rhs_ + stk_theta_* curv_->Data + stk_theta_* b_->Data ; /* TODO time-dep DirBC*/
}

template <class LsetSolverT, class RelaxationPolicyT>
void RecThetaScheme2PhaseCL<LsetSolverT,RelaxationPolicyT>::SetupLevelsetSystem()
{
    // setup system for levelset eq.
    LvlSet_.SetupSystem( Stokes_.GetVelSolution(), dt_);

    L_->LinComb( 1./dt_, LvlSet_.E, ls_theta_, LvlSet_.H);
    ls_rhs_ = LvlSet_.E * fixed_ls_rhs_;
}


// ==============================================
//            CrankNicolsonScheme2PhaseCL
// ==============================================

template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
CrankNicolsonScheme2PhaseCL<BaseMethod, LsetSolverT,RelaxationPolicyT>::CrankNicolsonScheme2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, NSSolverBaseCL<StokesT>& solver, LsetSolverT& lsetsolver, LevelsetModifyCL& lsetmod,
            double tol, double nonlinear, bool withProjection, double stab)
  : base_( Stokes, ls, solver, lsetsolver, lsetmod, tol, 1.0, 1.0, nonlinear, withProjection, stab), step_(1)
{}

template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
CrankNicolsonScheme2PhaseCL<BaseMethod, LsetSolverT,RelaxationPolicyT>::~CrankNicolsonScheme2PhaseCL()
{}

template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
void CrankNicolsonScheme2PhaseCL<BaseMethod, LsetSolverT,RelaxationPolicyT>::DoStep(int maxFPIter)
{
    switch (step_)
    {
        case 1 :
        {
            std::cout << "~~~~~~~~~~~~~~~~ initial backward euler step\n";
            tmpdt_=dt_;
            base_::SetTimeStep(0.2*tmpdt_*tmpdt_, 1.0);
            base_::DoStep(maxFPIter);
            base_::ComputeDots();
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

template< template<class, class> class BaseMethod, class LsetSolverT, class RelaxationPolicyT>
void CrankNicolsonScheme2PhaseCL<BaseMethod,LsetSolverT,RelaxationPolicyT>::Update()
{
    base_::SetTimeStep(tmpdt_, 1.0);
    base_::Update();
    step_= 1;
}

} // end of namespace DROPS
