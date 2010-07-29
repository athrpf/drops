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
 * Copyright 2010 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_COUPLING_CPP
#define DROPS_COUPLING_CPP

#include "levelset/coupling.h"

namespace DROPS
{
// ==============================================
//              TimeDisc2PhaseCL
// ==============================================

TimeDisc2PhaseCL::TimeDisc2PhaseCL (StokesT& Stokes, LevelsetP2CL& ls, LevelsetModifyCL& lsetmod, double nonlinear)
  : Stokes_( Stokes), LvlSet_( ls), b_( &Stokes.b), old_b_( new VelVecDescCL),
    cplM_( new VelVecDescCL), old_cplM_( new VelVecDescCL),
    cplN_( new VelVecDescCL), old_cplN_( new VelVecDescCL),
    curv_( new VelVecDescCL), old_curv_(new VelVecDescCL),
    rhs_( Stokes.b.RowIdx->NumUnknowns()), ls_rhs_( ls.Phi.RowIdx->NumUnknowns()),
    mat_( new MLMatrixCL( Stokes.vel_idx.size())), L_( new MatrixCL()), nonlinear_( nonlinear), lsetmod_( lsetmod)
{
    Stokes_.SetLevelSet( ls);
    LB_.Data.resize( Stokes.vel_idx.size());
}

TimeDisc2PhaseCL::~TimeDisc2PhaseCL()
{
    delete mat_; delete L_;
    if (old_b_ == &Stokes_.b)
        delete b_;
    else
        delete old_b_;
    delete cplM_; delete old_cplM_;
    delete cplN_; delete old_cplN_;
    delete curv_; delete old_curv_;
}

void cplDeltaSquaredPolicyCL::Update( VecDescCL& v)
{
    if (firststep_) {
        const size_t vsize = v.Data.size();
        v_old_.resize   ( vsize);
        v_diff_.resize  ( vsize);
        v_old_= v.Data;
        firststep_ = false;
        if (output_)
            (*output_) << "omega: " << omega_ << std::endl;
        return;
    }
    v_diff_=  v.Data - v_old_;
#ifndef _PAR
    omega_*= -dot( v_diff_, v_old_) / norm_sq( v_diff_);
#else
    const bool useAccur=true;
    ExchangeCL& ExVel  = v.RowIdx->GetEx();
    omega_*=-ExVel.ParDot( v_diff_, true, v_old_, true, useAccur) / ExVel.Norm_sq( v_diff_, true, useAccur);
#endif
    if (output_)
        (*output_) << "omega: " << omega_ << std::endl;
    v_old_= v.Data;
    v.Data   *= omega_;
}


void cplBroydenPolicyCL::Update( VecDescCL& v)
{
    F1_.push_back( v.Data);
#ifndef _PAR
    sigma_ = norm_sq( v.Data);
#else
    const bool useAccur=true;
    ExchangeCL& ExVel  = v.RowIdx->GetEx();
    sigma_ = ExVel.Norm_sq( v.Data, true, useAccur);
#endif

    if (output_)
        (*output_) << "sigma = " << sigma_ << std::endl;
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

    const double theta = std::sqrt(sigma_/sigma0_);
    if (theta >= thetamax_) {
        if (output_)
            (*output_) << "No convergence: theta = " << theta << std::endl;
        //return;
    }
    sigma0_ = sigma_;

    VectorCL w1( deltaF1_.back());
#ifndef _PAR
    gamma_.push_back( norm_sq( w1));
#else
    gamma_.push_back( ExVel.Norm_sq( w1, true, useAccur));
#endif

    kappa_ /= (1.0-2.0*theta);

    if (output_)
        (*output_) << "kappa = " << kappa_ << std::endl;
    if (std::fabs(kappa_) >= kappamax_){
        if (output_)
            (*output_) << "ill-conditioned update: kappa = "<< kappa_ << std::endl;
        //return;
    }

    VectorCL v1( F1_.back());
#ifndef _PAR
    const double factor = 1.0 - ( dot( w1, v1)) / gamma_.back();
#else
    const double factor = 1.0 - ( ExVel.ParDot( w1, true, v1, true, useAccur)) / gamma_.back();
#endif
    v1*= factor;

    for (int j=deltaF1_.size()-2; j>=0; --j) {
#ifndef _PAR
        const double beta = ( dot (deltaF1_[j], v1))/ gamma_[j];
#else
        const double beta = ( ExVel.ParDot (deltaF1_[j], true, v1, true, useAccur)) / gamma_[j];
#endif
        v1 -= beta*F1_[j+1];
    }
    v.Data   = v1;
}

} // end of namespace drops
#endif
