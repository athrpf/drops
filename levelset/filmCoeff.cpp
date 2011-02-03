/// \file filmCoeff.cpp
/// \brief boundary and source functions for film-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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

#include "misc/bndmap.h"
#include "levelset/params.h"

//========================================================================
//                        Functions for the film problem
//========================================================================

DROPS::Point3DCL FilmInflow( const DROPS::Point3DCL& p, double t)
{
    extern DROPS::ParamFilmCL C;
    DROPS::Point3DCL ret(0.);
    const double d= p[1]/C.exp_Thickness;
    static const double u= C.mat_DensFluid*C.exp_Gravity[0]*C.exp_Thickness*C.exp_Thickness/C.mat_ViscFluid/2;
    ret[0]= d<=1 ? (2*d-d*d)*u * (1 + C.exp_PumpAmpl*std::sin(2*M_PI*t*C.exp_PumpFreq))
                 : (C.mcl_MeshSize[1]-p[1])/(C.mcl_MeshSize[1]-C.exp_Thickness)*u;
    return ret;
}

//========================================================================
//        Registration of the function(s) in the func-container
//========================================================================

static DROPS::RegisterVectorFunction regvelfilminflow("FilmInflow", FilmInflow);


