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

extern DROPS::ParamFilmCL C;

//========================================================================
//                        Functions for the film problem
//========================================================================
namespace filminflow{

    DROPS::Point3DCL FilmInflow( const DROPS::Point3DCL& p, double t)
    {
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
}

//========================================================================
//                        Functions for the levelset function
//========================================================================
namespace filmdistance{
    double WavyDistanceFct( const DROPS::Point3DCL& p)
    {
        // wave length = 100 x film width
        const double wave= std::sin(2*M_PI*p[0]/C.mcl_MeshSize[0]),
            z= p[2]/C.mcl_MeshSize[2]*2; // z \in [-1,1]
    //    return p[1] - C.exp_Thickness * (1 + C.exp_PumpAmpl*wave);
    //    return p[1] - C.exp_Thickness * (1 + C.exp_PumpAmpl*(wave + C.exp_Ampl_zDir*std::cos(z*M_PI)));
        const double z_fac=  (1 + C.exp_Ampl_zDir/2*std::cos(z*M_PI));  // (z=+-1) 1-C.exp_Ampl_zDir <= z_fac <= 1+C.exp_Ampl_zDir (z=0)
        return p[1] - C.exp_Thickness * (1 + C.exp_PumpAmpl*wave) * z_fac;
    }

    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterScalarFunction regscafilmlset("WavyFilm", WavyDistanceFct);
}

