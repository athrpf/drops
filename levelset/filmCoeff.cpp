/// \file filmCoeff.cpp
/// \brief boundary and source functions for film-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Thorolf Schulte

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
#include "misc/params.h"

extern DROPS::ParamCL P;

//========================================================================
//                        Functions for the film problem
//========================================================================
namespace filminflow{

    DROPS::Point3DCL FilmInflow( const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ret(0.);
        const double d= p[1]/P.get<double>("Exp.Thickness");
        static const double u= P.get<double>("Mat.DensFluid")*P.get<DROPS::Point3DCL>("Exp.Gravity")[0]*P.get<double>("Exp.Thickness")*P.get<double>("Exp.Thickness")/P.get<double>("Mat.ViscFluid")/2;
        ret[0]= d<=1 ? (2*d-d*d)*u * (1 + P.get<double>("Exp.PumpAmpl")*std::sin(2*M_PI*t*P.get<int>("Exp.PumpFreq")))
                     : (P.get<DROPS::Point3DCL>("MeshSize")[1]-p[1])/(P.get<DROPS::Point3DCL>("MeshSize")[1]-P.get<double>("Exp.Thickness"))*u;
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
        const double wave= std::sin(2*M_PI*p[0]/P.get<DROPS::Point3DCL>("MeshSize")[0]),
            z= p[2]/P.get<DROPS::Point3DCL>("MeshSize")[2]*2; // z \in [-1,1]
    //    return p[1] - P.get<double>("Exp.Thickness") * (1 + P.get<double>("Exp.PumpAmpl")*wave);
    //    return p[1] - P.get<double>("Exp.Thickness") * (1 + P.get<double>("Exp.PumpAmpl")*(wave + P.get<double>("Exp.Ampl_zDir")*std::cos(z*M_PI)));
        const double z_fac=  (1 + P.get<double>("Exp.Ampl_zDir")/2*std::cos(z*M_PI));  // (z=+-1) 1-P.get<double>("Exp.Ampl_zDir") <= z_fac <= 1+P.get<double>("Exp.Ampl_zDir") (z=0)
        return p[1] - P.get<double>("Exp.Thickness") * (1 + P.get<double>("Exp.PumpAmpl")*wave) * z_fac;
    }

    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterScalarFunction regscafilmlset("WavyFilm", WavyDistanceFct);
}



//========================================================================
//                        Functions for matching function
//========================================================================
namespace filmperiodic{
    template<int A, int B>
    bool periodic_2sides( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    { 
        const DROPS::Point3DCL d= fabs(p-q),
                               L= fabs(P.get<DROPS::Point3DCL>("MeshSize"));
        
        const int D = 3 - A - B;
        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12)  // dB=dD=0 and dA=LA
          ||   (d[A] + d[D] < 1e-12 && std::abs( d[B] - L[B]) < 1e-12)  // dA=dD=0 and dB=LB
          ||   (d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12 && std::abs( d[B] - L[B]) < 1e-12);  // dD=0 and dA=LA and dB=LB
    }

    template<int A>
    bool periodic_1side( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    { 
        const int B = (A+1)%2;
        const int D = (B+1)%2;
        const DROPS::Point3DCL d= fabs(p-q), L= fabs(P.get<DROPS::Point3DCL>("MeshSize"));
        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12);
    }
    
    //========================================================================
    //        Registration of the function(s) in the func-container
    //========================================================================
    static DROPS::RegisterMatchingFunction regmatch2_xy("periodicxy", periodic_2sides<0,1>);
    static DROPS::RegisterMatchingFunction regmatch2_xz("periodicxz", periodic_2sides<0,2>);
    static DROPS::RegisterMatchingFunction regmatch2_yz("periodicyz", periodic_2sides<1,2>);
    static DROPS::RegisterMatchingFunction regmatch1_x("periodicx", periodic_1side<0>);
    static DROPS::RegisterMatchingFunction regmatch1_y("periodicy", periodic_1side<1>);
    static DROPS::RegisterMatchingFunction regmatch1_z("periodicz", periodic_1side<2>);
}
