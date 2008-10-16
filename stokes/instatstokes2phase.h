/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#ifndef DROPS_INSTATSTOKES2PHASE_H
#define DROPS_INSTATSTOKES2PHASE_H

#include <memory>

#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "num/MGsolver.h"

namespace DROPS
{

/// \brief Extended index for P1X_FE-elements.
///
/// Depending on the position of the zero-level of the levelset-function
/// additional dof are needed in the vertices. Internally, a std::vector
/// with a component for each vertex is kept. If the vertex has an
/// extended dof, the index is stored in this component; otherwise,
/// NoIdx is stored.
class ExtIdxDescCL
{
  private:
    double omit_bound_; ///< constant for stabilization of XFEM, controls omission of extended DoFs 
  public:
    typedef std::vector<IdxT> ExtendedIdxT;

    IdxDescCL* Idx; ///< Pointer to the index-description.

    ExtendedIdxT   Xidx;
    ExtendedIdxT   Xidx_old;

    ExtIdxDescCL(IdxDescCL* idx, double omit_bound= 1./32.) : omit_bound_(omit_bound), Idx( idx) {}

    IdxT operator[](const IdxT i) const { return Xidx[i]; }
    IdxT GetNumUnknownsP1() const { return Xidx.size(); }

    void UpdateXNumbering(IdxDescCL*, const LevelsetP2CL&, bool NumberingChanged= false);
    void Old2New(VecDescCL*);
};

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
    ExtIdxDescCL& extidx_;
    VectorCL extData_;
    VecDescCL& p_;

  public:
    P1XRepairCL( bool UsesXFEM, MultiGridCL& mg, VecDescCL& p,
        ExtIdxDescCL& extidx);

    void operator() (const LevelsetP2CL& lset);
};

/// problem class for instationary two-pase Stokes flow

template <class Coeff>
class InstatStokes2PhaseP2P1CL : public ProblemCL<Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, StokesBndDataCL>      _base;
    typedef InstatStokes2PhaseP2P1CL<Coeff>        _self;
    typedef typename _base::CoeffCL                CoeffCL;
    typedef typename _base::BndDataCL              BndDataCL;
    using                                          _base::_MG;
    using                                          _base::_Coeff;
    using                                          _base::_BndData;
    using                                          _base::GetBndData;
    using                                          _base::GetMG;

    typedef P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>   DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>   const_DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

  private:
    FiniteElementT prFE_;        ///< controls FE type for pressure space
    ExtIdxDescCL   Xidx_;        ///< extended index for P1X_FE

  protected:
    MGDataCL       matMG_;       ///< grid level data

  public:
    IdxDescCL    vel_idx;  ///< for velocity unknowns
    IdxDescCL    pr_idx;   ///< for pressure unknowns
    double       t;        ///< time
    VelVecDescCL v;        ///< velocity
    VecDescCL    p;        ///< pressure
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    &B,
                 &prM,
                 &M;
    MatDescCL    A,
                 prA;

    InstatStokes2PhaseP2P1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1)
        : _base(mgb, coeff, bdata), prFE_(prFE), Xidx_( &pr_idx, XFEMstab), matMG_(1), vel_idx(vecP2_FE), pr_idx(prFE), t( 0.), B(matMG_.begin()->B), prM(matMG_.begin()->Mpr), M(matMG_.begin()->Mvel) {}
    InstatStokes2PhaseP2P1CL( MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata, FiniteElementT prFE= P1_FE, double XFEMstab=0.1)
        : _base(mg, coeff, bdata),  prFE_(prFE), Xidx_( &pr_idx, XFEMstab), matMG_(1), vel_idx(vecP2_FE), pr_idx(prFE), t( 0.), B(matMG_.begin()->B), prM(matMG_.begin()->Mpr), M(matMG_.begin()->Mvel) {}

    /// \name Numbering
    //@{
    /// Create/delete numbering of unknowns
    void CreateNumberingVel( Uint level, IdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, _MG, _BndData.Vel, match); }
    void CreateNumberingPr ( Uint level, IdxDescCL* idx, match_fun match= 0, const LevelsetP2CL* lsetp= 0)
        { idx->CreateNumbering( level, _MG, _BndData.Pr, match); if (lsetp && UsesXFEM()) { Xidx_.UpdateXNumbering( idx, *lsetp, true); } }
    /// \brief Only used for XFEM
    void UpdateXNumbering( IdxDescCL* idx, const LevelsetP2CL& lset, bool NumberingChanged= false)
        { if (UsesXFEM()) Xidx_.UpdateXNumbering( idx, lset, NumberingChanged); }
    void UpdatePressure( VecDescCL* p)
        { if (UsesXFEM()) Xidx_.Old2New( p); }
    void DeleteNumbering( IdxDescCL* idx)
        { idx->DeleteNumbering( _MG); }
    //@}
    /// \name Discretization
    //@{
    /// Returns whether extended FEM are used
    bool UsesXFEM() const { return prFE_==P1X_FE; }
    /// Set up matrices A, M and rhs b (depending on phase bnd)
    void SetupSystem1( MatDescCL* A, MatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const;
    /// Set up rhs b (depending on phase bnd)
    void SetupRhs1( VecDescCL* b, const LevelsetP2CL& lset, double t) const;
    /// Set up matrices A, M on an arbitrary level; needed for MG-preconditioner
    void SetupMatrices1( MatDescCL* A, MatDescCL* M, const LevelsetP2CL& lset, double t) const;
    /// Set up MG-hierarchy
    void SetupMatricesMG( MGDataCL* matMG, const LevelsetP2CL& lset, double dt, double theta) const;
    /// Set up the Laplace-Beltrami-Operator
    void SetupLB( MatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const;
    /// Set up matrix B and rhs c
    void SetupSystem2( MatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up matrix B
    void SetupMatrix2( MatDescCL* B, MatDescCL* BT) const;
    /// Set up rhs c
    void SetupRhs2( VecDescCL* c, const LevelsetP2CL& lset, double t) const;
    /// Set up the time-derivative of B times velocity
    void SetupBdotv (VecDescCL* Bdotv, const VelVecDescCL* vel,
        const ExtIdxDescCL& v_idx, const LevelsetP2CL& lset, double t) const;
    /// Set up the mass matrix for the pressure, scaled by \f$\mu^{-1}\f$.
    void SetupPrMass( MatDescCL* prM, const LevelsetP2CL& lset) const;
    /// Set up the stiffness matrix for the pressure, scaled by \f$\rho^{-1}\f$.
    void SetupPrStiff(MatDescCL* prA, const LevelsetP2CL& lset) const;
    //@}

    /// Initialize velocity field
    void InitVel( VelVecDescCL*, instat_vector_fun_ptr, double t0= 0.) const;
    /// Smooth velocity field
    void SmoothVel( VelVecDescCL*, int num= 1, double tau=0.5);
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void ClearMat() { A.Data.clear(); B.Data.clear(); M.Data.clear(); prA.Data.clear(); prM.Data.clear(); }
    /// Set all indices
    void SetIdx();

    /// Get FE type for pressure space
    FiniteElementT GetPrFE() const { return prFE_; }
    /// \name Get extended index (only makes sense for P1X_FE)
    //@{
    const ExtIdxDescCL& GetXidx() const { return Xidx_; }
    ExtIdxDescCL&       GetXidx()       { return Xidx_; }
    //@}
    /// Get MG-data structure
    MGDataCL& GetMGData() {return matMG_;}
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

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
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

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  ();
    void post_refine_sequence ();
};

} // end of namespace DROPS

#include "stokes/instatstokes2phase.tpp"

#endif

