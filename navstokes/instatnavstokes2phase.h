/// \file instatnavstokes2phase.h
/// \brief classes that constitute the 2-phase Navier-Stokes problem
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_INSTATNAVSTOKES2PHASE_H
#define DROPS_INSTATNAVSTOKES2PHASE_H

#include "stokes/instatstokes2phase.h"
#include "levelset/levelset.h"
#include "num/renumber.h"


namespace DROPS
{

/// problem class for instationary two-phase Navier-Stokes flow

class InstatNavierStokes2PhaseP2P1CL : public InstatStokes2PhaseP2P1CL
{
  private:
    typedef InstatStokes2PhaseP2P1CL       base_;

    void SetupNonlinear_P2 (MatrixCL& N, const VelVecDescCL* vel, VelVecDescCL* cplN, const LevelsetP2CL& lset,
        IdxDescCL& RowIdx, double t) const;

  public:
    MLMatDescCL    N;
    const LevelsetP2CL* ls_;

    InstatNavierStokes2PhaseP2P1CL(const MGBuilderCL& mgb, const TwoPhaseFlowCoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab= 0.1, FiniteElementT velFE= vecP2_FE)
        : InstatStokes2PhaseP2P1CL( mgb, coeff, bdata, prFE, XFEMstab, velFE), ls_( 0) {}
    InstatNavierStokes2PhaseP2P1CL(MultiGridCL& mg, const TwoPhaseFlowCoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab= 0.1, FiniteElementT velFE= vecP2_FE)
        : InstatStokes2PhaseP2P1CL( mg, coeff, bdata, prFE, XFEMstab, velFE), ls_( 0) {}

    /// \name Discretization
    //@{
    /// \brief Set up matrix for nonlinearity
    void SetupNonlinear(MLMatDescCL* matN, const VelVecDescCL* vel, VelVecDescCL* cplN, const LevelsetP2CL& lset, double t) const;
    MLTetraAccumulatorTupleCL& nonlinear_accu (MLTetraAccumulatorTupleCL& accus, MLMatDescCL* matN, const VelVecDescCL* vel, VelVecDescCL* cplN, const LevelsetP2CL& lset, double t) const;
    /// \brief Set up matrix for nonlinearity at the time in the base-class using the registered Levelset-object.
    void SetupNonlinear(MLMatDescCL* matN, const VelVecDescCL* vel, VelVecDescCL* cplN) const {
        this->SetupNonlinear( matN, vel, cplN, *ls_, vel->t);
    }
    void SetupNonlinear(MatrixCL& N, const VelVecDescCL* vel, VelVecDescCL* cplN, IdxDescCL& RowIdx) const {
        this->SetupNonlinear_P2( N, vel, cplN, *ls_, RowIdx, vel->t);
    }
    //@}

    /// \brief Register a Levelset-object for use in SetupNonlinear; this is needed for Navier-Stokes-solvers.
    void SetLevelSet(const LevelsetP2CL& ls) { ls_= &ls; }
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void ClearMat() { base_::ClearMat(); N.Data.clear(); }
    void SetIdx()   { base_::SetIdx(); N.SetIdx(&vel_idx, &vel_idx); }
    void SetNumVelLvl( size_t n) { base_::SetNumVelLvl( n); N.Data.resize (vel_idx.size()); }

    /// \brief Perform downwind numbering for the velocity FE-space. The permutation is returned.
    PermutationT downwind_numbering (const LevelsetP2CL& lset, IteratedDownwindCL dw);
};

} // end of namespace DROPS

#endif
