//**************************************************************************
// File:    instatstokes2phase.h                                           *
// Content: classes that constitute the 2-phase stokes-problem             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Sep, 1 2003                                            *
//**************************************************************************

#ifndef DROPS_INSTATSTOKES2PHASE_H
#define DROPS_INSTATSTOKES2PHASE_H

#include "stokes/instatstokes.h"
#include "levelset/levelset.h"

namespace DROPS
{

// InstatStokes2PhaseP2P1CL: problem class for instationary two pase stokes flow

template <class Coeff>
class InstatStokes2PhaseP2P1CL : public ProblemCL<Coeff, InstatStokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, InstatStokesBndDataCL>      _base;
    typedef InstatStokes2PhaseP2P1CL<Coeff>              _self;
    typedef typename _base::CoeffCL                      CoeffCL;
    typedef typename _base::BndDataCL                    BndDataCL;
    using                                                _base::_MG;
    using                                                _base::_Coeff;
    using                                                _base::_BndData;
    using                                                _base::GetBndData;
    using                                                _base::GetMG;

    typedef P1EvalCL<double, const InstatStokesPrBndDataCL, const VecDescCL>   DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const InstatStokesVelBndDataCL, const VelVecDescCL> DiscVelSolCL;

    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    double       t;        // time
    VelVecDescCL v;        // velocity
    VecDescCL    p;        // pressure
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A, 
                 B,
                 M;
    
    InstatStokes2PhaseP2P1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mgb, coeff, bdata), vel_idx(3,3), pr_idx(1), t( 0.) {}  
    InstatStokes2PhaseP2P1CL( MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mg, coeff, bdata), vel_idx(3,3), pr_idx(1), t( 0.) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Vel, match); }
    void CreateNumberingPr ( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Pr, match); }
    void DeleteNumberingVel( IdxDescCL*);
    void DeleteNumberingPr ( IdxDescCL*);
    
    // Set up matrices A, M and rhs b (depending on phase bnd)
    void SetupSystem1( MatDescCL* A, MatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const;
    // Set up matrices A, M on an arbitrary level; needed for MG-preconditioner
    void SetupMatrices1( MatDescCL* A, MatDescCL* M, const LevelsetP2CL& lset, double t) const;
    // Set up matrix B and rhs c (independent of phase bnd, but c is time dependent)
    void SetupSystem2( MatDescCL* B, VecDescCL* c, double t) const;

    // Set up rhs c (time dependent)
    void SetupRhs2( VecDescCL* c, double t) const;
    // Set up mass and stiffness matrix for pressure-unknowns (P1, time-independent) 
    // needed for preconditioning of the Schur complement
    void SetupPrMass( MatDescCL* prM, const LevelsetP2CL& lset, double nu1=1, double nu2=1) const;
    void SetupPrStiff(MatDescCL* prA) const;
    void InitVel( VelVecDescCL*, vector_instat_fun_ptr, double t0= 0.) const;

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution() const
        { return DiscPrSolCL( &p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution() const
        { return DiscVelSolCL( &v, &GetBndData().Vel, &GetMG(), t); }

    DiscPrSolCL GetPrSolution( const VecDescCL& pr) const
        { return DiscPrSolCL( &pr, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution( const VelVecDescCL& vel) const
        { return DiscVelSolCL( &vel, &GetBndData().Vel, &GetMG(), t); }
};

} // end of namespace DROPS

#include "stokes/instatstokes2phase.tpp"

#endif

