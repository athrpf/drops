//**************************************************************************
// File:    instatnavstokes2phase.h                                        *
// Content: classes that constitute the 2-phase navier-stokes-problem      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Dec, 16 2001                                           *
//**************************************************************************

#ifndef DROPS_INSTATNAVSTOKES2PHASE_H
#define DROPS_INSTATNAVSTOKES2PHASE_H

#include "stokes/instatstokes2phase.h"


namespace DROPS
{

template <class Coeff>
class InstatNavierStokes2PhaseP2P1CL : public InstatStokes2PhaseP2P1CL<Coeff>
{
  private:
    typedef InstatStokes2PhaseP2P1CL<Coeff>       _base;
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> _self;

  public:
    using _base::GetBndData;
    using _base::GetMG;
    using _base::_Coeff;
    using _base::_MG;
    using _base::_BndData;
    using _base::b;
    using _base::c;
    using _base::A;
    using _base::B;
    
    typedef Coeff                        CoeffCL;
    typedef typename _base::BndDataCL    BndDataCL;
    typedef typename _base::DiscVelSolCL DiscVelSolCL;
  
    MatDescCL    N;
    
    InstatNavierStokes2PhaseP2P1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : InstatStokes2PhaseP2P1CL<Coeff>( mgb, coeff, bdata) {}  

    // Set up matrix for nonlinearity; 
    // couplings with dir.bnd are accumulated in cplN,
    // so call cplN->Clear() before if only couplings are needed.
    void SetupNonlinear(MatDescCL* matN, const VelVecDescCL* vel, VelVecDescCL* cplN, const LevelsetP2CL& lset, double t) const;
};

} // end of namespace DROPS

#include "navstokes/instatnavstokes2phase.tpp"

#endif
