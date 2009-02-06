//**************************************************************************
// File:    integrTime.cpp                                                 *
// Content: Stokes-preconditioner an time integration                      *
// Author:  Sven Gross, Joerg Grande, Hieu Nguyen, IGPM RWTH Aachen        *
//          Oliver Fortmeiere, SC RWTH Aachen                              *
// Version: 0.1                                                            *
// History: begin - Jan, 24 2008                                           *
//**************************************************************************

#include <stokes/integrTime.h>

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

void ISBBTPreCL::Update() const
{
    IF_MASTER
      std::cerr << "ISBBTPreCL::Update: old version: " << Bversion_
                << "\tnew version: " << B_->Version() << '\n';
    delete Bs_;
    Bs_= new MatrixCL( *B_);
    Bversion_= B_->Version();

#ifndef _PAR
    VectorCL Dvelinv( 1.0/ Mvel_->GetDiag());
#else
    BBT_.SetBlock0( Bs_);
    BBT_.SetBlock1( Bs_);
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
    std::cerr << "MinCommPreCL::Update: old/new versions: " << Aversion_  << '/' << A_->Version()
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
#endif

} // end of namespace DROPS
