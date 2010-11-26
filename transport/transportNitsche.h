/// \file transportNitsche.h
/// \brief Classes that constitute a 2-phase-transport-problem with Nitsche-XFEM discretization.
/// \author Trung Hieu Nguyen (small fixes: Martin Horsky, Christoph Lehrenfeld), IGPM

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

#ifndef DROPS_TRANSPORTNITSCHE_H
#define DROPS_TRANSPORTNITSCHE_H

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
  /**
   * TODO: check how and if t is used, do we need it as a class member?
   * 
   */
  
  
class TransportP1XCL
{
  public:
    typedef BndDataCL<>                                       BndDataT;
    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;

    double oldt_, t_;
    MLIdxDescCL idx, oldidx;
    VecDescCL c,           ///< concentration
              ct, oldct;   ///< transformed concentration
    MLMatDescCL A,         ///< diffusion matrix
                M,         ///< mass matrix
                C,        ///< convection matrix
                NA;       ///< Nitsche term 
    VecDescCL cplA, cplM, cplC, b; 
    
  private:
    MLMatrixCL     L_;              ///< sum of matrices
    MultiGridCL&   MG_;
    BndDataT&      Bnd_, Bndt_;     ///< Boundary condition for the concentration
    VelBndDataT    Bnd_v_;          ///< Boundary condition for the velocity
    LsetBndDataCL& Bnd_ls_;          ///< Boundary condition for the level set function
    double         theta_,          ///< time scheme parameter
                   dt_,             ///< time step 
                   lambda_,         ///< Nitsche's parameter
                   D_[2],           ///< diffusion constants
                   H_; 
    VecDescCL *v_, *oldv_;         ///< velocity at current time step and previous step
    VecDescCL &lset_,  &oldlset_;  ///< levelset at current time step and previous step
    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;
    instat_scalar_fun_ptr f_;

    void SetupInstatSystem(MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, VecDescCL*,
        IdxDescCL&, const double) const;
    void SetupNitscheSystem( MatrixCL&, IdxDescCL&)const;
    void SetupNitscheSystem( MatrixCL&, IdxDescCL&, bool)const;
    void SetupInstatMixedMassMatrix( MatrixCL&, VecDescCL*, IdxDescCL&, IdxDescCL&, const double) const;
    void SetupInstatMixedSystem(MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, VecDescCL*,
        IdxDescCL&, IdxDescCL&, const double) const;
    void SetupMixedNitscheSystem( MatrixCL&, IdxDescCL&, IdxDescCL&, const double) const;

  public:
    TransportP1XCL( MultiGridCL& mg, BndDataT& Bnd, BndDataT& Bndt, const VelBndDataT& Bnd_v, LsetBndDataCL& Bnd_ls,
        double theta, double D[2], double H, VecDescCL* v, VecDescCL* oldv, VecDescCL& lset, VecDescCL& oldlset,
        double dt, instat_scalar_fun_ptr rhs=0, int iter= 1000, double tol= 1e-7,
        double XFEMstab=0.1, double lambda=0.)
        : oldt_(v->t), t_( v->t), 
        idx( P1X_FE, 1, Bndt, 0, XFEMstab), oldidx( P1X_FE, 1, Bndt, 0, XFEMstab), 
        MG_( mg), Bnd_( Bnd), Bndt_( Bndt), Bnd_v_( Bnd_v), Bnd_ls_(Bnd_ls), 
        theta_( theta), dt_( dt), 
        lambda_(lambda), H_( H), v_( v), oldv_(oldv), 
        lset_( lset), oldlset_(oldlset), 
        gm_( pc_, 20, iter, tol, false, false, RightPreconditioning),
        f_(rhs) 
    {
        std::memcpy( D_, D, 2*sizeof( double));
        if (theta_!=1.0) std::cerr<< "TransportP1XCL::TransportP1XCL: Sorry. Cannot deal with theta != 1.0. Overwriting: Using implicit Euler now! \n";
    }

    TransportP1XCL( MultiGridCL& mg, BndDataT& Bnd, BndDataT& Bndt, const VelBndDataT& Bnd_v, LsetBndDataCL& Bnd_ls,
        VecDescCL* v, VecDescCL* oldv, VecDescCL& lset, VecDescCL& oldlset,
        DROPS::ParamMesszelleNsCL& C, instat_scalar_fun_ptr rhs=0)
        : oldt_(v->t), t_( v->t), 
        idx( P1X_FE, 1, Bndt, 0, C.trp_NitscheXFEMStab), oldidx( P1X_FE, 1, Bndt, 0, C.trp_NitscheXFEMStab), 
        MG_( mg), Bnd_( Bnd), Bndt_( Bndt), Bnd_v_( Bnd_v), Bnd_ls_(Bnd_ls), 
        theta_( C.trp_Theta), dt_( C.tm_StepSize), 
        lambda_(C.trp_NitschePenalty), H_( C.trp_HNeg/C.trp_HPos), v_( v), oldv_(oldv), 
        lset_( lset), oldlset_(oldlset), 
        gm_( pc_, 20, C.trp_Iter, C.trp_Tol, false, false, RightPreconditioning),
        f_(rhs) 
    {
        double D[2] = {C.trp_DiffPos, C.trp_DiffNeg};
        std::memcpy( D_, D, 2*sizeof( double));
        if (theta_!=1.0) std::cerr<< "TransportP1XCL::TransportP1XCL: Sorry. Cannot deal with theta != 1.0. Overwriting: Using implicit Euler now! \n";
    }




    MultiGridCL& GetMG()               { return MG_; }
    const MultiGridCL& GetMG() const   { return MG_; }
    GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }
    const BndDataT& GetBndData() const { return Bndt_; }
    
    /// \brief Only used for XFEM
    void UpdateXNumbering( const VecDescCL& lset, bool= false, bool NewSolution=true) {
         if (NewSolution)
            idx.UpdateXNumbering( MG_, lset, Bnd_ls_);
         else
            oldidx.UpdateXNumbering( MG_, lset, Bnd_ls_);
    }

    void CreateNumbering( Uint level, MLIdxDescCL* idx1, MLIdxDescCL *oldidx1, VecDescCL& lset1, VecDescCL& oldlset1) {
         idx1->CreateNumbering( level, MG_, &lset1, &Bnd_ls_);
         oldidx1->CreateNumbering( level, MG_, &oldlset1, &Bnd_ls_);
    }
    void DeleteNumbering( MLIdxDescCL* idx1) { idx1->DeleteNumbering( MG_); }
    /// initialize transformed concentration function
    void Init( instat_scalar_fun_ptr, instat_scalar_fun_ptr, double t=0.0);
    /// Set one P1X Function as the other with according scaling parameters in the positive and negative domain parts
    void TransformWithScaling(const VecDescCL& concin, VecDescCL& concout, double scalingp, double scalingn);

    void SetTimeStep( double dt, double theta=-1);
    void SetTwoStepIdx();
    void SetNewIdx();
    void SetupInstatSystem( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, VecDescCL&, const double) const;
    void SetupNitscheSystem( MLMatDescCL&) const;
    void SetupNitscheSystem( MLMatDescCL&, bool) const;
    void SetupInstatMixedMassMatrix( MLMatDescCL&, VecDescCL&, const double) const;
    void SetupInstatMixedSystem( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, VecDescCL&, const double) const;
    void SetupMixedNitscheSystem( MLMatDescCL&, const double) const;

    /// perform one time step
    void DoStep( double new_t);
    void SetRHS( instat_scalar_fun_ptr rhs) {f_= rhs;}
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &c, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myc, bool Is_ct) const
        { return (Is_ct) ? const_DiscSolCL( &Myc, &Bndt_, &MG_)
                         : const_DiscSolCL( &Myc, &Bnd_, &MG_); }
    const_DiscVelSolCL GetVelocity() const
        { return const_DiscVelSolCL( v_, &Bnd_v_, &MG_); }
    const_DiscVelSolCL GetVelocity(const VecDescCL& vel) const
        { return const_DiscVelSolCL( &vel, &Bnd_v_, &MG_); }

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the locling navstokes-levelset. They should not be called by a common user.
    /// Use LevelsetP2CL::DoStep() instead.
    ///@{
    void InitStep (VectorCL&);
    void DoStep (const VectorCL&);
    void CommitStep ();
    ///@}
    /// Get FE type for pressure space
    VecDescCL& GetLevelset ()  {return lset_;}
    const VecDescCL& GetLevelset ()  const {return lset_;}
    LsetBndDataCL& GetLevelsetBnd ()  {return Bnd_ls_;}
    const LsetBndDataCL& GetLevelsetBnd ()  const {return Bnd_ls_;}
    VecDescCL& GetOldLevelset ()  {return oldlset_;}
    const VecDescCL& GetOldLevelset ()  const {return oldlset_;}
    //@}
    void GetSolutionOnPart( VecDescCL& ct_part, bool posPart, bool Is_ct );
    double Interface_L2error() const;
    void SetupInstatRhs( VecDescCL & b,  const double time) const;
    double MeanDropConcentration();
    double GetHenry(bool positive = false){
      if (positive)
        return 1.0;
      else
        return H_;
    }
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the function c.ct.
///
/// The actual work is done in post_refine().

class TransportXRepairCL : public MGObserverCL
{
  private:
    TransportP1XCL& c_;
    std::auto_ptr<P1XRepairCL> oldp1xrepair_;

  public:
    TransportXRepairCL (TransportP1XCL& c)
        : c_( c) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  (); 
    void post_refine_sequence (); 
    const IdxDescCL* GetIdxDesc() const;
    
};

class VelTranspRepairCL : public MGObserverCL
{
  public:
    typedef BndDataCL<Point3DCL>                                          VelBndDataT;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL>    const_DiscVelSolCL;
  private:
    MultiGridCL& mg_;
    VecDescCL& v_;
    VecDescCL& lset_;
    const VelBndDataT& Bnd_v_;
    LsetBndDataCL& Bnd_ls_;
    IdxDescCL& vidx_;
    double& time_;

  public:
    VelTranspRepairCL (VecDescCL& v, MultiGridCL& mg, const VelBndDataT& Bnd_v, IdxDescCL& vidx, VecDescCL& lset, LsetBndDataCL& Bnd_ls, double t )
        :  mg_(mg), v_(v), lset_(lset), Bnd_v_(Bnd_v), Bnd_ls_(Bnd_ls), vidx_(vidx) , time_(t){}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    void SetTime(double t) {time_=t;}
    const IdxDescCL* GetIdxDesc() const;
};

} // end of namespace DROPS

#endif
