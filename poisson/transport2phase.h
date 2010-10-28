/// \file transport2phase.h
/// \brief Classes that constitute a 2-phase-transport-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Hieu Nguyen, Marcus Soemers, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_TRANSPORT2PHASE_H
#define DROPS_TRANSPORT2PHASE_H

#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/solver.h"
#include "num/bndData.h"
#include <iostream>
#include <numeric>
#include <cstring>

namespace DROPS
{

/// \brief P1-discretization and solution of the transport equation for two phase flow problems.
class TransportP1CL
{
  public:
    typedef BndDataCL<>                                       BndDataT;
    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    MLIdxDescCL idx;
    VecDescCL   c,  ///< concentration
                ct; ///< transformed concentration
    MLMatDescCL A,  ///< diffusion matrix
                M,  ///< mass matrix
                C;  ///< convection matrix
    VecDescCL   cplA,    cplM,    cplC,
                oldcplA, oldcplM, oldcplC;

  private:
    MLMatrixCL   L_;              ///< sum of matrices
    MultiGridCL& MG_;
    double       D_[2],           ///< diffusion constants
                 H_,              ///< jump of concentration at the interface
                 theta_, dt_;     ///< time scheme parameter and time step
    BndDataT&    Bnd_;            ///< Boundary condition for the concentration

    const VelBndDataT&  Bnd_v_;  ///< Boundary condition for the velocity
    VecDescCL*          v_;      ///< velocity at current time step
    LevelsetP2CL&       lset_;   ///< levelset at current time step

    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;

    void SetupInstatSystem (MatrixCL& matA, VecDescCL* cplA,
                        MatrixCL& matM, VecDescCL* cplM, MatrixCL& matC, VecDescCL* cplC,
                        IdxDescCL& RowIdx, const double time) const;

  public:
    TransportP1CL( MultiGridCL& mg, BndDataT& Bnd, const VelBndDataT& Bnd_v,
        double theta, double D[2], double H, VecDescCL* v, LevelsetP2CL& lset,
        double dt, int iter= 1000, double tol= 1e-7)
    : idx( P1_FE), MG_( mg), H_( H), theta_( theta), dt_( dt), Bnd_( Bnd),
        Bnd_v_( Bnd_v), v_( v), lset_( lset), gm_( pc_, 100, iter, tol, true)
    { std::memcpy( D_, D, 2*sizeof( double)); }

    const MultiGridCL& GetMG() const { return MG_; }
    GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }
    const BndDataT& GetBndData() const { return Bnd_; }

    /// \name Numbering
    ///@{
    void CreateNumbering( Uint level, MLIdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, MG_, Bnd_, match); }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    ///@}
    /// initialize transformed concentration function
    void Init( instat_scalar_fun_ptr, instat_scalar_fun_ptr cpos);

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetTimeStep( double dt, double theta=-1);
    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetupLocalSystem (const TetraCL&, double[4][4], double[4][4], double[4][4], const double,
        const LocalP2CL<>[4], const LocalP2CL<>[4][4], const Quad5CL<>[4]) const;
    void SetupInstatSystem ( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, const double) const;

    /// perform one time step
    void DoStep( double new_t);

    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &c, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myc) const
        { return const_DiscSolCL( &Myc, &Bnd_, &MG_); }
    ///@}

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the coupling navstokes-levelset. They should not be called by a common user.
    /// Use LevelsetP2CL::DoStep() instead.
    ///@{
    void InitStep (VectorCL&);
    void DoStep (const VectorCL&);
    void CommitStep ();
    /// The transformed concentration ct must be repaired externally; c is restored from ct.
    void Update ();
    ///@}

    /// The following member transform between the continuous an discontinuous version of c.
    ///@{
    void c2ct ();
    void ct2c ();
    ///@}
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the function c.ct.
///
/// The actual work is done in post_refine().
class TransportRepairCL : public MGObserverCL
{
  private:
    TransportP1CL& c_;
    MultiGridCL& mg_;

  public:
    TransportRepairCL (TransportP1CL& c, MultiGridCL& mg)
        : c_( c), mg_( mg) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    const IdxDescCL* GetIdxDesc() const { return c_.c.RowIdx; }
};

} // end of namespace DROPS

#endif
