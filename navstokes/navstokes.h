//**************************************************************************
// File:    navstokes.h                                                    *
// Content: classes that constitute the navier-stokes-problem              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 30 2001                                         *
//**************************************************************************

#ifndef DROPS_NAVSTOKES_H
#define DROPS_NAVSTOKES_H

#include "stokes/stokes.h"


namespace DROPS
{

template <class Coeff>
class NavierStokesP2P1CL : public StokesP2P1CL<Coeff>
{
  private:
    typedef StokesP2P1CL<Coeff> _base;

  public:
    using                            _base::_MG;
    using                            _base::_BndData;
    using                            _base::b;
    using                            _base::c;
    using                            _base::A;
    using                            _base::B;
    using                            _base::t; // Initialized with 0.0 by base-class.

    typedef Coeff                     CoeffCL;
    typedef typename _base::BndDataCL BndDataCL;
    typedef typename _base::DiscVelSolCL DiscVelSolCL;
    typedef typename _base::DiscPrSolCL DiscPrSolCL;

    MatDescCL    N;
    VelVecDescCL cplN;
    VelVecDescCL cplM;

    NavierStokesP2P1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : StokesP2P1CL<Coeff>( mgb, coeff, bdata) {}  
    NavierStokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : StokesP2P1CL<Coeff>( mg, coeff, bdata) {}  

    // Set up matrix and rhs for nonlinearity: use time t1 for the velocity in N,
    // t2 for the boundary-data in the velocity unknowns
    void SetupNonlinear(MatDescCL*, const VelVecDescCL*, VelVecDescCL*, double, double) const;
    // Set up matrix for nonlinearity, use time _t
    void SetupNonlinear(MatDescCL* matN, const VelVecDescCL* velvec, VelVecDescCL* vecb) const
    { this->SetupNonlinear(matN, velvec, vecb, t, t); }

    // Set time for use with stationary NavStokes-Solvers. This shall be the new time t_old+dt!!!!!!!!!!!!!!!!!!
    void SetTime (double tt) { t= tt; }

    // Check computed solution
    void CheckSolution(const VelVecDescCL*, const VecDescCL*,
        instat_vector_fun_ptr, instat_scalar_fun_ptr, double= 0.0);
};


} // end of namespace DROPS

#include "navstokes/navstokes.tpp"

#endif
