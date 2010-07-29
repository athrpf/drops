/// \file instatstokes2phase.tpp
/// \brief classes that constitute the 2-phase Stokes problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

namespace DROPS
{

/// \name Routines for SetupSystem2
//@{
/// \brief P2 / P0 FEs for vel/pr
void SetupSystem2_P2P0( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2 / P1 FEs (Taylor-Hood) for vel/pr
void SetupSystem2_P2P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2 / P1X FEs (X=extended) for vel/pr
void SetupSystem2_P2P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2X / P1X FEs (X=extended) for vel/pr
void SetupSystem2_P2RP1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2X / P1 FEs (X=extended) for vel/pr
void SetupSystem2_P2RP1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);

/// \brief P2 / P1D FEs for vel/pr
void SetupSystem2_P2P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t);
//@}


/// \name Routines for SetupRhs2
//@{
/// \brief P2 / P0 FEs for vel/pr
void SetupRhs2_P2P0( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t);

/// \brief P2 / P1 FEs (Taylor-Hood) for vel/pr
void SetupRhs2_P2P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t);

/// \brief P2 / P1X FEs (X=extended) for vel/pr
void SetupRhs2_P2P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL& lset, double t);

/// \brief P2 / P1D FEs for vel/pr
void SetupRhs2_P2P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t);
//@}


/// \name Routines for SetupPrMass
//@{
/// \brief P0 FEs for pr
void SetupPrMass_P0(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief P1 FEs for pr
void SetupPrMass_P1(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief P1X FEs for pr
void SetupPrMass_P1X(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);

/// \brief PD FEs for pr
void SetupPrMass_P1D(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset);
//@}


/// \name Routines for SetupPrSiff
//@{
/// \brief P1 FEs for pr
void SetupPrStiff_P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset);

/// \brief P1X FEs for pr
/// \todo: As in SetupPrMass_P1X, replace the smoothed density-function with integration
///        over the inner and outer part.
void SetupPrStiff_P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset);

/// \brief P1D FEs for pr
void SetupPrStiff_P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset);
//@}


//*****************************************************************************
//                               VelocityRepairCL
//*****************************************************************************
#ifndef _PAR
inline void VelocityRepairCL::pre_refine()
  /// do nothing
{ }
#else
inline void VelocityRepairCL::pre_refine()
  /// tell parallel multigrid about velocities
{
    GetPMG().AttachTo( &stokes_.v, &stokes_.GetBndData().Vel);
}
#endif

inline void
  VelocityRepairCL::post_refine ()
{
    VelVecDescCL loc_v;
    VelVecDescCL& v= stokes_.v;
    Uint LastLevel= stokes_.GetMG().GetLastLevel();
    match_fun match= stokes_.GetMG().GetBnd().GetMatchFun();
    IdxDescCL loc_vidx( vecP2_FE);

    loc_vidx.CreateNumbering( LastLevel, stokes_.GetMG(), stokes_.GetBndData().Vel, match);
    if (LastLevel != v.RowIdx->TriangLevel()) {
        std::cout << "LastLevel: " << LastLevel
                  << " old v->TriangLevel(): " << v.RowIdx->TriangLevel() << std::endl;
        throw DROPSErrCL( "VelocityRepairCL::post_refine: Sorry, not yet implemented.");
    }
    loc_v.SetIdx( &loc_vidx);
#ifdef _PAR
    GetPMG().HandleNewIdx(&stokes_.vel_idx, &loc_v);
#endif
    RepairAfterRefineP2( stokes_.GetVelSolution( v), loc_v);
#ifdef _PAR
    GetPMG().CompleteRepair( &loc_v);
#endif
    v.Clear();
    v.RowIdx->DeleteNumbering( stokes_.GetMG());

    stokes_.vel_idx.GetFinest().swap( loc_vidx);
    v.SetIdx( &stokes_.vel_idx);
    v.Data= loc_v.Data;
}

inline void
  VelocityRepairCL::post_refine_sequence ()
  /// Create numbering for all idx level
{
    match_fun match= stokes_.GetMG().GetBnd().GetMatchFun();
    stokes_.CreateNumberingVel( stokes_.GetMG().GetLastLevel(), &stokes_.vel_idx, *match);
}


//*****************************************************************************
//                               PressureRepairCL
//*****************************************************************************

#ifndef _PAR
inline void PressureRepairCL::pre_refine()
  /// do nothing
{ }
#else
inline void PressureRepairCL::pre_refine()
  /// tell parallel multigrid about (extended) finite element function
{
    GetPMG().AttachTo( &stokes_.p, &stokes_.GetBndData().Pr);
    if ( stokes_.UsesXFEM()){
        GetPMG().AttachTo( p1xrepair_->GetExt(), &stokes_.GetBndData().Pr);
    }
}
#endif

inline void
  PressureRepairCL::post_refine ()
{
    VecDescCL loc_p;
    IdxDescCL loc_pidx( stokes_.GetPrFE());
    VecDescCL& p= stokes_.p;
    match_fun match= stokes_.GetMG().GetBnd().GetMatchFun();

    loc_pidx.CreateNumbering( stokes_.GetMG().GetLastLevel(), stokes_.GetMG(), stokes_.GetBndData().Pr, match, &ls_.Phi);
    loc_p.SetIdx( &loc_pidx);
#ifdef _PAR
    GetPMG().HandleNewIdx(&stokes_.pr_idx, &loc_p);
#endif
    RepairAfterRefineP1( stokes_.GetPrSolution( p), loc_p);
#ifdef _PAR
    GetPMG().CompleteRepair( &loc_p);
    if ( stokes_.UsesXFEM()){
        IdxDescCL loc_xpr_idx( P1_FE);
        loc_xpr_idx.CreateNumbering( stokes_.GetMG().GetLastLevel(), stokes_.GetMG());
        VecDescCL loc_xpr( &loc_xpr_idx);
        GetPMG().HandleNewIdx( p1xrepair_->GetExt()->RowIdx, &loc_xpr);
        // swap index and data
        p1xrepair_->GetExt()->RowIdx->swap( loc_xpr_idx);
        p1xrepair_->GetExt()->Data.resize( loc_xpr.Data.size());
        p1xrepair_->GetExt()->Data= loc_xpr.Data;
        // delete this numbering
        loc_xpr_idx.DeleteNumbering( stokes_.GetMG());
    }
#endif
    p.Clear();
    p.RowIdx->DeleteNumbering( stokes_.GetMG());
    stokes_.pr_idx.GetFinest().swap( loc_pidx);
    p.SetIdx( &stokes_.pr_idx);
    p.Data= loc_p.Data;
}

inline void
  PressureRepairCL::pre_refine_sequence ()
{
    p1xrepair_= std::auto_ptr<P1XRepairCL>( new P1XRepairCL( stokes_.GetMG(), stokes_.p));
}

inline void
  PressureRepairCL::post_refine_sequence ()
  /// Create numbering for all idx level
{
    match_fun match= stokes_.GetMG().GetBnd().GetMatchFun();
    stokes_.CreateNumberingPr( stokes_.GetMG().GetLastLevel(), &stokes_.pr_idx, match, &ls_);
    (*p1xrepair_)();
    p1xrepair_.reset();
}
} // end of namespace DROPS
