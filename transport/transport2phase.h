/// \file
/// \brief Classes that constitute a 2-phase-transport-problem.

//**************************************************************************
// File:    transport2phase.h                                              *
// Content: classes that constitute the transport-problem                  *
// Author:  Joerg Grande, Sven Gross, Maxim Larin, Hieu Nguyen, IGPM RWTH Aachen*
// Version: 0.1                                                            *
// History: begin - August 2007                                            *
//**************************************************************************

#ifndef DROPS_TRANSPORT2PHASE_H
#define DROPS_TRANSPORT2PHASE_H

#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/solver.h"
#include "num/bndData.h"
#include "stokes/instatstokes2phase.h"
#include "misc/xfem.h"
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

    double t_;
    MLIdxDescCL idx;
    VecDescCL   c,  ///< concentration
                ct; ///< transformed concentration
    MLMatDescCL A,  ///< diffusion matrix
                M,  ///< mass matrix
                C;  ///< convection matrix
    VecDescCL   cplA,    cplM,    cplC, b,
                oldcplA, oldcplM, oldcplC, oldb;
    
  private:
    MLMatrixCL   L_;              ///< sum of matrices
    MultiGridCL& MG_;
    BndDataT&    Bnd_;            ///< Boundary condition for the concentration
    const VelBndDataT&  Bnd_v_;  ///< Boundary condition for the velocity
    double       D_[2],           ///< diffusion constants
                 H_,              ///< jump of concentration at the interface
                 theta_, dt_; ///< time scheme parameter, time step and time
    VecDescCL*          v_;      ///< velocity at current time step
    LevelsetP2CL&       lset_;   ///< levelset at current time step

    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;
    instat_scalar_fun_ptr   f_;  ///< right hand side
    void SetupInstatSystem (MatrixCL& matA, VecDescCL* cplA,
                        MatrixCL& matM, VecDescCL* cplM, MatrixCL& matC, VecDescCL* cplC, VecDescCL* cplb, VecDescCL* oldcplM,
                        IdxDescCL& RowIdx, const double time, const double time_old) const;

  public:
    TransportP1CL( MultiGridCL& mg, BndDataT& Bnd, const VelBndDataT& Bnd_v, 
        double theta, double D[2], double H, VecDescCL* v, LevelsetP2CL& lset,
        double t, double dt, instat_scalar_fun_ptr rhs=0, int iter= 1000, double tol= 1e-7)
    :  MG_( mg), Bnd_( Bnd), Bnd_v_( Bnd_v),theta_( theta), H_( H), v_( v), lset_( lset),
    t_( t),dt_( dt), f_(rhs), gm_( pc_, 100, iter, tol, true), idx( P1_FE)
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
    void Init( instat_scalar_fun_ptr, instat_scalar_fun_ptr cpos, double t = 0.0);

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetTimeStep( double dt, double theta=-1);
    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetupLocalSystem (const TetraCL&, double[4][4], double[4][4], double[4][4], const double,
        const LocalP2CL<>[4], const LocalP2CL<>[4][4], const Quad5CL<>[4]) const;
    void SetupInstatSystem ( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, VecDescCL&, VecDescCL&, const double, const double) const;
    
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
};

class VelTranspRepairCL_P1 : public MGObserverCL
{
  public:
    typedef BndDataCL<Point3DCL>                                          VelBndDataT;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL>    const_DiscVelSolCL;
  private:
    TransportP1CL& c_;
    MultiGridCL& mg_;
    VecDescCL& v_;
    VelBndDataT Bnd_v_;
    IdxDescCL& vidx_;

  public:
    VelTranspRepairCL_P1 (TransportP1CL& c, VecDescCL& v, MultiGridCL& mg, VelBndDataT Bnd_v, IdxDescCL& vidx)
        : c_(c), v_(v), mg_(mg), Bnd_v_(Bnd_v), vidx_(vidx) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
};

} // end of namespace DROPS

#endif
