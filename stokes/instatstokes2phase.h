/// \file
/// \brief classes that constitute the 2-phase stokes-problem

#ifndef DROPS_INSTATSTOKES2PHASE_H
#define DROPS_INSTATSTOKES2PHASE_H

#include "stokes/stokes.h"
#include "levelset/levelset.h"

namespace DROPS
{

/// problem class for instationary two pase stokes flow

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

    IdxDescCL    vel_idx;  ///< for velocity unknowns
    IdxDescCL    pr_idx;   ///< for pressure unknowns
    double       t;        ///< time
    VelVecDescCL v;        ///< velocity
    VecDescCL    p;        ///< pressure
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A, 
                 B,
                 M;
    
    InstatStokes2PhaseP2P1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mgb, coeff, bdata), vel_idx(3,3), pr_idx(1), t( 0.) {}  
    InstatStokes2PhaseP2P1CL( MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mg, coeff, bdata), vel_idx(3,3), pr_idx(1), t( 0.) {}  

    //@{
    /// Create/delete numbering of unknowns
    void CreateNumberingVel( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Vel, match); }
    void CreateNumberingPr ( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Pr, match); }
    void DeleteNumberingVel( IdxDescCL*);
    void DeleteNumberingPr ( IdxDescCL*);
    //@}
    
    /// Set up matrices A, M and rhs b (depending on phase bnd)
    void SetupSystem1( MatDescCL* A, MatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const;
    /// Set up matrices A, M on an arbitrary level; needed for MG-preconditioner
    void SetupMatrices1( MatDescCL* A, MatDescCL* M, const LevelsetP2CL& lset, double t) const;
    /// Set up matrix B and rhs c (independent of phase bnd, but c is time dependent)
    void SetupSystem2( MatDescCL* B, VecDescCL* c, double t) const;
    /// Set up rhs c (time dependent)
    void SetupRhs2( VecDescCL* c, double t) const;
    
    // Set up mass/stiffness matrix for pressure-unknowns (P1, time-independent) 
    // needed for preconditioning of the Schur complement
    void SetupPrMass( MatDescCL* prM, const LevelsetP2CL& lset) const;
    void SetupPrStiff(MatDescCL* prA, const LevelsetP2CL& lset) const;

    /// Initialize velocity field
    void InitVel( VelVecDescCL*, instat_vector_fun_ptr, double t0= 0.) const;

    //@{
    /// Get solution as FE-function
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

} // end of namespace DROPS

#include "stokes/instatstokes2phase.tpp"

#endif

