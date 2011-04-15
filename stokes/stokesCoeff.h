/// \file stokesCoeff.h
/// \brief  coefficient class for stokes problems
/// \author LNM RWTH Aachen: Eva Loch, Yuanjun Zhang

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

#ifndef STOKESCOEFF_H_
#define STOKESCOEFF_H_

#include "misc/bndmap.h"
#include "misc/container.h"
#include <sstream>

namespace DROPS{

typedef DROPS::Point3DCL (*vector_fun_ptr)       (const DROPS::Point3DCL&);
typedef double    (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);
typedef DROPS::Point3DCL (*instat_vector_fun_ptr)(const DROPS::Point3DCL&, double);
typedef DROPS::SMatrixCL<3, 3> (*instat_matrix_fun_ptr) (const DROPS::Point3DCL&, double);

class StokesFlowCoeffCL
{
  public:
  //reaction
    static instat_scalar_fun_ptr q;
  //source term
    static instat_vector_fun_ptr f;
  //solution of velocity
    static instat_vector_fun_ptr LsgVel;
  //solution of Jacobi-matrix of exact solution for velocity
    static instat_matrix_fun_ptr DLsgVel;
    //solution of pressure
    static instat_scalar_fun_ptr LsgPr;

    const double rho, nu;
    const DROPS::Point3DCL g;

    StokesFlowCoeffCL( const DROPS::ParamStokesProblemCL& C)
      : rho( C.mat_Dens),
        nu( C.mat_Visc),
        g( C.exp_Gravity){
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        q = scamap[C.stc_Reaction];
        DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
        f = vecmap[C.stc_Source];
        LsgVel = vecmap[C.stc_Solution_Vel];
        DROPS::InMatMap & matmap = DROPS::InMatMap::getInstance();
        DLsgVel = matmap[C.stc_Solution_DVel];
        DROPS::InScaMap & inscamap = DROPS::InScaMap::getInstance();
        LsgPr = inscamap[C.stc_Solution_Pr];
    }



};

instat_scalar_fun_ptr StokesFlowCoeffCL::q;
instat_vector_fun_ptr StokesFlowCoeffCL::f;
instat_vector_fun_ptr StokesFlowCoeffCL::LsgVel;
instat_matrix_fun_ptr StokesFlowCoeffCL::DLsgVel;
instat_scalar_fun_ptr StokesFlowCoeffCL::LsgPr;
}//end of namespace

#endif
