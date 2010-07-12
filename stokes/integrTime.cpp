/// \file integrTime.cpp
/// \brief Stokes-preconditioner an time integration
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Hieu Nguyen; SC RWTH Aachen: Oliver Fortmeier

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

#include "stokes/integrTime.h"
#include "stokes/instatstokes2phase.h"
#include "num/stokessolver.h"

namespace DROPS
{

#ifndef _PAR
// Append the kernel of Bs as last column to Bs.
static void Regularize (MatrixCL& Bs, const IdxDescCL& rowidx, VectorCL ker0, const NEGSPcCL& spc, double regularize)
{
    if (rowidx.IsExtended())
        ker0[std::slice( rowidx.GetXidx().GetNumUnknownsStdFE(), rowidx.NumUnknowns() - rowidx.GetXidx().GetNumUnknownsStdFE(), 1)]= 0.;
    ker0*= 1./norm( ker0);
    VectorCL ker( spc.mul( Bs, ker0));
    ker*= regularize/std::sqrt( dot( ker, ker0));
    Bs.insert_col( Bs.num_cols(), ker);
}
#endif

void ISMGPreCL::MaybeInitOnes() const
{
    if (Mpr_.size() == ones_.size()) return;
    // Compute projection on constant pressure function only once.
    Uint i= 0;
    ones_.resize(0); // clear all
    ones_.resize(Mpr_.size());
    for (MLMatrixCL::const_iterator it= Mpr_.begin(); it != Mpr_.end(); ++it, ++i) {
        ones_[i].resize( it->num_cols(), 1.0/it->num_cols());
    }
}

void ISBBTPreCL::Update() const
{
    IF_MASTER
      std::cout << "ISBBTPreCL::Update: old version: " << Bversion_
                << "\tnew version: " << B_->Version() << '\n';
    delete Bs_;
    Bs_= new MatrixCL( *B_);
    Bversion_= B_->Version();

#ifndef _PAR
    VectorCL Dvelinv( 1.0/ Mvel_->GetDiag());
#else
    BBT_.SetBlock0( Bs_);
    BBT_.SetBlock1( Bs_);
    if (solver_.GetPC().NeedDiag())
        solver_.GetPC().SetDiag(BBT_);
    VectorCL Dvelinv( 1.0/ vel_idx_->GetEx().GetAccumulate(Mvel_->GetDiag()));
#endif
    ScaleCols( *Bs_, VectorCL( std::sqrt( Dvelinv)));

#ifndef _PAR
    VectorCL Dprsqrt( std::sqrt( M_->GetDiag()));
#else
    VectorCL Dprsqrt( std::sqrt( pr_idx_->GetEx().GetAccumulate( M_->GetDiag())));
#endif
    Dprsqrtinv_.resize( M_->num_rows());
    Dprsqrtinv_= 1.0/Dprsqrt;
    ScaleRows( *Bs_, Dprsqrtinv_);

#ifndef _PAR
    if (regularize_ != 0.)
        Regularize( *Bs_, *pr_idx_, Dprsqrt, spc_, regularize_);
#endif
}

#ifndef _PAR
void MinCommPreCL::Update() const
{
    std::cout << "MinCommPreCL::Update: old/new versions: " << Aversion_  << '/' << A_->Version()
        << '\t' << Bversion_ << '/' << B_->Version() << '\t' << Mversion_ << '/' << M_->Version()
        << '\t' << Mvelversion_ << '/' << Mvel_->Version() << '\n';
    delete Bs_;
    Bs_= new MatrixCL( *B_);
    Aversion_= A_->Version();
    Bversion_= B_->Version();
    Mversion_= M_->Version();
    Mvelversion_= Mvel_->Version();

    Assert( Mvel_->GetDiag().min() > 0., "MinCommPreCL::Update: Mvel_->GetDiag().min() <= 0\n", DebugNumericC);
    VectorCL Dvelsqrt( std::sqrt( Mvel_->GetDiag()));
    Dvelsqrtinv_.resize( Mvel_->num_rows());
    Dvelsqrtinv_= 1.0/Dvelsqrt;
    ScaleCols( *Bs_, Dvelsqrtinv_);

    Assert( M_->GetDiag().min() > 0., "MinCommPreCL::Update: M_->GetDiag().min() <= 0\n", DebugNumericC);
    VectorCL Dprsqrt( std::sqrt( M_->GetDiag()));
    Dprsqrtinv_.resize( M_->num_rows());
    Dprsqrtinv_= 1.0/Dprsqrt;
    ScaleRows( *Bs_, Dprsqrtinv_);

    if (regularize_ != 0.)
        Regularize( *Bs_, *pr_idx_, Dprsqrt, spc_, regularize_);
}


void BDinvBTPreCL::Update() const
{
    std::cout << "BDinvBTPreCL::Update: old/new versions: " << Lversion_  << '/' << L_->Version()
        << '\t' << Bversion_ << '/' << B_->Version() << '\t' << Mversion_ << '/' << M_->Version()
        << '\t' << Mvelversion_ << '/' << Mvel_->Version() << '\n';
    delete Bs_;
    Bs_= new MatrixCL( *B_);
    Lversion_= L_->Version();
    Bversion_= B_->Version();
    Mversion_= M_->Version();
    Mvelversion_= Mvel_->Version();

    Dvelinv_.resize( Mvel_->num_rows());
    if (lumped_)
        Dvelinv_= 1.0/LumpInRows(*Mvel_);
    else
        Dvelinv_= 1.0/L_->GetDiag();
    // note: Dvelinv_ has negative entries for P2-FE (hence, CGNE cannot be used for solving), however, B Dvelinv_ B^T has only positive diagonal entries
    VectorCL Dprsqrt( std::sqrt( M_->GetDiag()));
    Dprsqrtinv_.resize( M_->num_rows());
    Dprsqrtinv_= 1.0/Dprsqrt;

    ScaleRows( *Bs_, Dprsqrtinv_);
    DSchurinv_.resize( Dprsqrt.size());    
    DSchurinv_= 1.0/Bs_->GetSchurDiag(Dvelinv_);
    delete BDinvBT_;
    BDinvBT_= new AppSchurComplMatrixT( *L_, diagVelPc_, *Bs_);
    
    if (regularize_ != 0.)
        Regularize( *Bs_, *pr_idx_, Dprsqrt, spc_, regularize_);
}

BDinvBTPreCL::~BDinvBTPreCL() 
{
    delete Bs_;
    delete BDinvBT_;
}
#endif

} // end of namespace DROPS
