/// \file twophaseCoeff.cpp
/// \brief boundary and source functions for the twophasedrops-type problems
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

extern DROPS::ParamMesszelleNsCL C;

//========================================================================
//          Functions for twophasedrops-executable
//========================================================================
namespace tpd_inflow{
    /// \name inflow condition
    DROPS::SVectorCL<3> InflowBrick( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double x = p[0]*(2* C.exp_RadInlet -p[0]) / (C.exp_RadInlet*C.exp_RadInlet),
                     z = p[2]*(2*C.exp_RadInlet-p[2]) / (C.exp_RadInlet*C.exp_RadInlet);

        ret[1]= x * z * C.exp_InflowVel * (1-C.exp_InflowAmpl*std::cos(2*M_PI*C.exp_InflowFreq*t));

        return ret;
    }


    ///microchannel (eindhoven)
    DROPS::SVectorCL<3> InflowChannel( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double y = p[1]*(2*25e-6-p[1]) / (25e-6*25e-6),
                     z = p[2]*(2*50e-6-p[2]) / (50e-6*50e-6);

        ret[0]= y * z * C.exp_InflowVel * (1-C.exp_InflowAmpl*std::cos(2*M_PI*C.exp_InflowFreq*t));

        return ret;
    }

    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    DROPS::SVectorCL<3> InflowCell( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double s2= C.exp_RadInlet*C.exp_RadInlet,
                     r2= p.norm_sq() - p[C.exp_FlowDir]*p[C.exp_FlowDir];
        ret[C.exp_FlowDir]= -(r2-s2)/s2*C.exp_InflowVel;

        return ret;
    }


    //========================================================================
    //                       Functions for brick_transp.cpp
    //========================================================================

    DROPS::SVectorCL<3> InflowBrickTransp (const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double x = p[0]*(2*C.exp_RadInlet-p[0]) / (C.exp_RadInlet*C.exp_RadInlet),
                     z = p[2]*(2*C.exp_RadInlet-p[2]) / (C.exp_RadInlet*C.exp_RadInlet);
        ret[1]= x * z * C.exp_InflowVel;
        return ret;
    }

    //========================================================================
    //            Registration of functions in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelbrick("InflowBrick", InflowBrick);
    static DROPS::RegisterVectorFunction regvelcell("InflowCell", InflowCell);
    static DROPS::RegisterVectorFunction regvelchannel("InflowChannel", InflowChannel);
    static DROPS::RegisterVectorFunction regvelbricktransp("InflowBrickTransp", InflowBrickTransp);
}

//========================================================================
//                       Functions for LevelSet Distance
//========================================================================
namespace levelsetdistance{

    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    double CubeDistance( const DROPS::Point3DCL& p)
    {
        double maxd = - C.exp_RadDrop[0] - C.exp_RadDrop[1]- C.exp_RadDrop[2];
        for (int i=0;i<3;i++){
          double x = std::abs(p[i] - C.exp_PosDrop[i]) - C.exp_RadDrop[i];
          if (x>maxd) maxd=x;
        }
        return maxd;
    }

    template<int i>
    double planedistance( const DROPS::Point3DCL& p)
    {
        double x=p[i]-C.exp_PosDrop[i]; 
        return x;
    }

    static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
    static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
    static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
    static DROPS::RegisterScalarFunction regscacube("CubeDistance", CubeDistance);
}

//========================================================================
//                       Functions for transport
//========================================================================
namespace transpfunctions{
    double tInitialcneg (const DROPS::Point3DCL& , double)
    {
        return C.trp_IniCNeg;
    }

    double tInitialcpos (const DROPS::Point3DCL& , double)
    {
        return C.trp_IniCPos;
    }
    static DROPS::RegisterScalarFunction regscainineg("Initialcpos", tInitialcneg);
    static DROPS::RegisterScalarFunction regscainipos("Initialcpos", tInitialcpos);
}

//========================================================================
//                       Functions for surfactants
//========================================================================
namespace surffunctions{
    /// \name Initial data and rhs for surfactant transport
    //@{
    const double a( -13./8.*std::sqrt( 35./M_PI));
    double surf_rhs (const DROPS::Point3DCL& p, double)
    {
        return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
    }
    double surf_sol (const DROPS::Point3DCL& p, double)
    {
        return 1. + std::sin( atan2( p[0] - C.exp_PosDrop[0], p[2] - C.exp_PosDrop[2]));  
    }
    //@}
    static DROPS::RegisterScalarFunction regscasurfrhs("surf_rhs", surf_rhs);
    static DROPS::RegisterScalarFunction regscasurfsol("surf_sol", surf_sol);  
}

//TODO: unification with filmCoeff stoff
//========================================================================
//                        Functions for matching function
//========================================================================
namespace filmperiodic{
    template<int A, int B>
    bool periodic_2sides( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    { 

        static bool first = true;
        static DROPS::Point3DCL dx;
        //dirty hack
        if (first){
            std::string mesh( C.dmc_MeshFile), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> dx[0] >> dx[1] >> dx[2] ;
            first = false;
        }      
        
        const DROPS::Point3DCL d= fabs(p-q),
                               L= fabs(dx);
        
        const int D = 3 - A - B;
        return (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12)  // dB=dD=0 and dA=LA
          ||   (d[A] + d[D] < 1e-12 && std::abs( d[B] - L[B]) < 1e-12)  // dA=dD=0 and dB=LB
          ||   (d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12 && std::abs( d[B] - L[B]) < 1e-12);  // dD=0 and dA=LA and dB=LB
    }

    template<int A>
    bool periodic_1side( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
    { 

        static bool first = true;
        static DROPS::Point3DCL dx;
        //dirty hack
        if (first){
            std::string mesh( C.dmc_MeshFile), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> dx[0] >> dx[1] >> dx[2] ;
            first = false;
        }      
        const int B = (A+1)%3;
        const int D = (B+1)%3;
        const DROPS::Point3DCL d= fabs(p-q), L= fabs(dx);
        bool res = (d[B] + d[D] < 1e-12 && std::abs( d[A] - L[A]) < 1e-12);
        
        DROPS::Point3DCL diff = p - q;
        if (res) std::cout << " diff = " << diff << std::endl;
        return res;
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