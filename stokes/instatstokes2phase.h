/// \file instatstokes2phase.h
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

#ifndef DROPS_INSTATSTOKES2PHASE_H
#define DROPS_INSTATSTOKES2PHASE_H

#include <memory>

#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "num/MGsolver.h"
#include "misc/xfem.h"

namespace DROPS
{

/// \brief Repair a P1X-vector if grid changes occur
///
/// Create such an object with the variable to be repaired before any grid-modifications.
/// Repair the linear part however you like and call the operator() to repair the extended part.
class P1XRepairCL
{
  private:
    bool UsesXFEM_;
    MultiGridCL& mg_;
    IdxDescCL idx_;
    VecDescCL ext_;
    VecDescCL& p_;

  public:
    P1XRepairCL( MultiGridCL& mg, VecDescCL& p);

    VecDescCL* GetExt() { return &ext_; }

    void operator() ();
};

/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag_P1 (const MultiGridCL& MG, VectorCL& M, IdxDescCL& RowIdx,
                       const BndDataCL<>& bnd= BndDataCL<>( 0));

/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag_P1X (const MultiGridCL& MG, VectorCL& M, IdxDescCL& RowIdx, const VecDescCL& lset,
                        const BndDataCL<>& bnd= BndDataCL<>( 0));

/// \brief Compute the main diagonal of the unscaled \f$L_2(\Omega)\f$-mass-matrix.
void SetupMassDiag_vecP2 (const MultiGridCL& MG, VectorCL& M, IdxDescCL& RowIdx,
                          const BndDataCL<Point3DCL>& bnd= BndDataCL<Point3DCL>( 0));


/// problem class for instationary two-pase Stokes flow

template <class Coeff>
class InstatStokes2PhaseP2P1CL : public ProblemCL<Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, StokesBndDataCL>       _base;
    typedef InstatStokes2PhaseP2P1CL<Coeff>         _self;
    typedef typename _base::CoeffCL                 CoeffCL;
    typedef typename _base::BndDataCL               BndDataCL;
    using                                           _base::_MG;
    using                                           _base::_Coeff;
    using                                           _base::_BndData;
    using                                           _base::GetBndData;
    using                                           _base::GetMG;

    typedef P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>   DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>   const_DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

  public:
    MLIdxDescCL  vel_idx;  ///< for velocity unknowns
    MLIdxDescCL  pr_idx;   ///< for pressure unknowns
    double       t;        ///< time
    VelVecDescCL v;        ///< velocity
    VecDescCL    p;        ///< pressure
    VelVecDescCL b;
    VecDescCL    c;
    MLMatDescCL  A,
                 B,
                 M,
                 prA,
                 prM;

  public:
    InstatStokes2PhaseP2P1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1)
        : _base(mgb, coeff, bdata), vel_idx(vecP2_FE, 1, bdata.Vel), pr_idx(prFE, 1, bdata.Pr, 0, XFEMstab), t( 0.) {}
    InstatStokes2PhaseP2P1CL( MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1)
        : _base(mg, coeff, bdata),  vel_idx(vecP2_FE, 1, bdata.Vel), pr_idx(prFE, 1, bdata.Pr, 0, XFEMstab), t( 0.) {}

    /// \name Numbering
    //@{
    /// Create/delete numbering of unknowns
    void CreateNumberingVel( Uint level, MLIdxDescCL* idx, match_fun match= 0);
    void CreateNumberingPr ( Uint level, MLIdxDescCL* idx, match_fun match= 0, const LevelsetP2CL* lsetp= 0);
    /// \brief Only used for XFEM
    void UpdateXNumbering( MLIdxDescCL* idx, const LevelsetP2CL& lset)
        {
            if (UsesXFEM()) idx->UpdateXNumbering( _MG, lset.Phi);
        }
    /// \brief Only used for XFEM
    void UpdatePressure( VecDescCL* p)
        {
            if (UsesXFEM()) p->RowIdx->GetXidx().Old2New( p);
        }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( _MG); }
    //@}

    /// \name Discretization
    //@{
    /// Returns whether extended FEM are used
    bool UsesXFEM() const { return pr_idx.GetFinest().IsExtended(); }
    /// Set up matrices A, M and rhs b (depending on phase bnd)
    void SetupSystem1( MLMatDescCL* A, MLMatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const;
    /// Set up rhs b (depending on phase bnd)
    void SetupRhs1( VecDescCL* b, const LevelsetP2CL& lset, double t) const;
    /// Set up the Laplace-Beltrami-Operator
    void SetupLB( MLMatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const;
    /// Set up matrix B and rhs c
    void SetupSystem2( MLMatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up rhs c
    void SetupRhs2( VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up the time-derivative of B times velocity
    void SetupBdotv (VecDescCL* Bdotv, const VelVecDescCL* vel, const LevelsetP2CL& lset, double t) const;
    /// Set up the mass matrix for the pressure, scaled by \f$\mu^{-1}\f$.
    void SetupPrMass( MLMatDescCL* prM, const LevelsetP2CL& lset) const;
    /// Set up the stiffness matrix for the pressure, scaled by \f$\rho^{-1}\f$.
    void SetupPrStiff(MLMatDescCL* prA, const LevelsetP2CL& lset) const;
    //@}

    /// Initialize velocity field
    void InitVel( VelVecDescCL*, instat_vector_fun_ptr, double t0= 0.) const;
    /// Smooth velocity field
    void SmoothVel( VelVecDescCL*, int num= 1, double tau=0.5);
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void ClearMat() { A.Data.clear(); B.Data.clear(); M.Data.clear(); prA.Data.clear(); prM.Data.clear(); }
    /// Set all indices
    void SetIdx();
    /// Set number of used levels
    void SetNumVelLvl( size_t n);
    void SetNumPrLvl ( size_t n);
    /// Get FE type for pressure space
    FiniteElementT GetPrFE() const { return pr_idx.GetFinest().GetFE(); }
    /// \name Get extended index (only makes sense for P1X_FE)
    //@{
    const ExtIdxDescCL& GetXidx() const { return pr_idx.GetFinest().GetXidx(); }
    ExtIdxDescCL&       GetXidx()       { return pr_idx.GetFinest().GetXidx(); }
    //@}
    /// Get pressure solution on inner/outer part (especially for P1X_FE)
    void GetPrOnPart( VecDescCL& p_part, const LevelsetP2CL& lset, bool posPart= true); // false = inner = Phi<0, true = outer = Phi>0

    /// \name Evaluate Solution
    //@{
    /// Get solution as FE-function for evaluation
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL( &p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL( &v, &GetBndData().Vel, &GetMG(), t); }

    const_DiscPrSolCL GetPrSolution( const VecDescCL& pr) const
        { return const_DiscPrSolCL( &pr, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution( const VelVecDescCL& vel) const
        { return const_DiscVelSolCL( &vel, &GetBndData().Vel, &GetMG(), t); }
    //@}
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function stokes_.v.
///
/// The actual work is done in post_refine().
template<class StokesT>
class VelocityRepairCL : public MGObserverCL
{
  private:
    StokesT& stokes_;

  public:
    VelocityRepairCL (StokesT& stokes)
        : stokes_( stokes) {}
    void pre_refine  ();
    void post_refine ();
    void pre_refine_sequence  () {}
    void post_refine_sequence ();
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function stokes_.pr.
///
/// For the P1-part, the actual work is done in post_refine().
/// For the P1X-part, a P1XRepairCL is created in pre_refine_sequence() and used in
/// post_refine_sequence(). Holding the P1XRepairCL* in an auto_ptr simplifies the use
/// of heap-memory: No memory is lost, even if successive calls of pre_refine_sequence()
/// occur without interleaved post_refine_sequence()-calls.
template<class StokesT>
class PressureRepairCL : public MGObserverCL
{
  private:
    StokesT& stokes_;
    std::auto_ptr<P1XRepairCL> p1xrepair_;
    const LevelsetP2CL& ls_;

  public:
    PressureRepairCL (StokesT& stokes, const LevelsetP2CL& ls)
        : stokes_( stokes), ls_( ls) {}
    void pre_refine  ();
    void post_refine ();
    void pre_refine_sequence  ();
    void post_refine_sequence ();
};

} // end of namespace DROPS

#include "stokes/instatstokes2phase.tpp"

#endif

