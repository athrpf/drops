/// \file
/// \brief classes that constitute the Navier-Stokes problem
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

/// History: begin - April, 30 2001

#ifndef DROPS_NAVSTOKES_H
#define DROPS_NAVSTOKES_H

#include "stokes/stokes.h"


namespace DROPS
{

template <class Coeff>
class NavierStokesP2P1CL : public StokesP2P1CL<Coeff>
{
  private:
    typedef StokesP2P1CL<Coeff> base_;
    void SetupNonlinear_P2( MatrixCL&, const VelVecDescCL*, VelVecDescCL*, IdxDescCL&, double) const;

  public:
    using                            base_::MG_;
    using                            base_::BndData_;
    using                            base_::v;
    using                            base_::b;
    using                            base_::c;
    using                            base_::A;
    using                            base_::B;

    typedef Coeff                     CoeffCL;
    typedef typename base_::BndDataCL BndDataCL;
    typedef typename base_::DiscVelSolCL DiscVelSolCL;
    typedef typename base_::DiscPrSolCL DiscPrSolCL;
    typedef typename base_::const_DiscVelSolCL const_DiscVelSolCL;
    typedef typename base_::const_DiscPrSolCL const_DiscPrSolCL;

    MLMatDescCL  N;
    VelVecDescCL cplN;
    VelVecDescCL cplM;

    NavierStokesP2P1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : StokesP2P1CL<Coeff>( mgb, coeff, bdata) {}
    NavierStokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : StokesP2P1CL<Coeff>( mg, coeff, bdata) {}

    // Set up matrix and rhs for nonlinearity: use time t1 for the velocity in N,
    // t2 for the boundary-data in the velocity unknowns
    void SetupNonlinear(MLMatDescCL*, const VelVecDescCL*, VelVecDescCL*, double) const;
    // Set up matrix for nonlinearity, use time v.t
    void SetupNonlinear(MLMatDescCL* matN, const VelVecDescCL* velvec, VelVecDescCL* vecb) const
    { this->SetupNonlinear(matN, velvec, vecb, v.t); }
    void SetupNonlinear(MatrixCL& N, const VelVecDescCL* vel, VelVecDescCL* cplN, IdxDescCL& RowIdx) const
    { this->SetupNonlinear_P2( N, vel, cplN, RowIdx, v.t); }

    // Set time for use with stationary NavStokes-Solvers. This shall be the new time t_old+dt!!!!!!!!!!!!!!!!!!
    void SetTime (double tt) { v.t= tt; }
    // in parallel version the error is arisen in the super  class
    void SetNumVelLvl( size_t n) {base_::SetNumVelLvl(n); N.Data.resize( n);}

    // Check computed solution
    void CheckSolution(const VelVecDescCL*, MLIdxDescCL* idx, const VecDescCL*,
        instat_vector_fun_ptr, instat_scalar_fun_ptr);
};


} // end of namespace DROPS

#include "navstokes/navstokes.tpp"

#endif
