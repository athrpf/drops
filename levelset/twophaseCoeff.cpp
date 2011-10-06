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
#include "misc/params.h"

extern DROPS::ParamCL P;

//========================================================================
//          Functions for twophasedrops-executable (inflow)
//========================================================================
namespace tpd_inflow{
    /// \name inflow condition
    DROPS::SVectorCL<3> InflowBrick( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double x = p[0]*(2* P.get<double>("Exp.RadInlet") -p[0]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet")),
                     z = p[2]*(2*P.get<double>("Exp.RadInlet")-p[2]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet"));

        ret[1]= x * z * P.get<double>("Exp.InflowVel") * (1-P.get<double>("Exp.InflowAmpl")*std::cos(2*M_PI*P.get<double>("Exp.InflowFreq")*t));

        return ret;
    }


    ///microchannel (eindhoven)
    DROPS::SVectorCL<3> InflowChannel( const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double y = p[1]*(2*25e-6-p[1]) / (25e-6*25e-6),
                     z = p[2]*(2*50e-6-p[2]) / (50e-6*50e-6);

        ret[0]= y * z * P.get<double>("Exp.InflowVel") * (1-P.get<double>("Exp.InflowAmpl")*std::cos(2*M_PI*P.get<double>("Exp.InflowFreq")*t));

        return ret;
    }

    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    DROPS::SVectorCL<3> InflowCell( const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double s2= P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet"),
                     r2= p.norm_sq() - p[P.get<int>("Exp.FlowDir")]*p[P.get<int>("Exp.FlowDir")];
        ret[P.get<int>("Exp.FlowDir")]= -(r2-s2)/s2*P.get<double>("Exp.InflowVel");

        return ret;
    }


    //========================================================================
    //                       Functions for brick_transp.cpp
    //========================================================================

    DROPS::SVectorCL<3> InflowBrickTransp (const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret(0.);
        const double x = p[0]*(2*P.get<double>("Exp.RadInlet")-p[0]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet")),
                     z = p[2]*(2*P.get<double>("Exp.RadInlet")-p[2]) / (P.get<double>("Exp.RadInlet")*P.get<double>("Exp.RadInlet"));
        ret[1]= x * z * P.get<double>("Exp.InflowVel");
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
//          Functions for twophasedrops-executable (Volume Force)
//========================================================================
namespace tpd_volforce{
    /// \name inflow condition
    template<int D>
    DROPS::SVectorCL<3> PeriodicDropletPressure( const DROPS::Point3DCL& , double )
    {
        DROPS::SVectorCL<3> ret(0.);
        ret[D] = -P.get<DROPS::Point3DCL>("Exp.Gravity")[D];
        
        static bool first = true;
        static DROPS::Point3DCL dx;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> dx[0] >> dx[1] >> dx[2] ;
            first = false;
        }      
        
        static double voldrop = 4./3.*M_PI* P.get<DROPS::Point3DCL>("Exp.RadDrop")[0]*P.get<DROPS::Point3DCL>("Exp.RadDrop")[1]*P.get<DROPS::Point3DCL>("Exp.RadDrop")[2] ;
        static double brickvol = dx[0]* dx[1]* dx[2];
        static double volforce = P.get<double>("Mat.DensFluid") * brickvol - (P.get<double>("Mat.DensFluid") - P.get<double>("Mat.DensDrop")) * voldrop;
        ret[D] *= volforce/brickvol;
        return ret;
    }

    //========================================================================
    //            Registration of functions in the func-container
    //========================================================================
    static DROPS::RegisterVectorFunction regvelppx("PeriodicDropletPressurex", PeriodicDropletPressure<0>);
    static DROPS::RegisterVectorFunction regvelppy("PeriodicDropletPressurey", PeriodicDropletPressure<1>);
    static DROPS::RegisterVectorFunction regvelppz("PeriodicDropletPressurez", PeriodicDropletPressure<2>);
}


//========================================================================
//                       Functions for LevelSet Distance
//========================================================================
namespace levelsetdistance{

    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    double CubeDistance( const DROPS::Point3DCL& p)
    {
        double maxd = - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0] - P.get<DROPS::Point3DCL>("Exp.RadDrop")[1]- P.get<DROPS::Point3DCL>("Exp.RadDrop")[2];
        for (int i=0;i<3;i++){
          double x = std::abs(p[i] - P.get<DROPS::Point3DCL>("Exp.PosDrop")[i]) - P.get<DROPS::Point3DCL>("Exp.RadDrop")[i];
          if (x>maxd) maxd=x;
        }
        return maxd;
    }


    ///mzelle_ns_adap.cpp + mzelle_instat.cpp
    template<int i>    
    double PeriodicEllipsoidDistance( const DROPS::Point3DCL& p)
    {
      
        static bool first = true;
        static DROPS::Point3DCL dx;
        //dirty hack
        if (first){
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
            size_t idx_;
            while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx_]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> dx[0] >> dx[1] >> dx[2] ;
            first = false;
        }      
              
        DROPS::Point3DCL dp;
        dp[i] = dx[i];
        DROPS::Point3DCL ExpRadDrop = P.get<DROPS::Point3DCL>("Exp.RadDrop");
        DROPS::Point3DCL ExpPosDrop = P.get<DROPS::Point3DCL>("Exp.PosDrop");
        DROPS::Point3DCL d= p - ExpPosDrop;
        DROPS::Point3DCL d1= p + dp - ExpPosDrop;
        DROPS::Point3DCL d2= p - dp - ExpPosDrop;
        const double avgRad= cbrt(ExpRadDrop[0]*ExpRadDrop[1]*ExpRadDrop[2]);
        d/= ExpRadDrop;
        d1/= ExpRadDrop;
        d2/= ExpRadDrop;
        double dd = std::min(std::min(d.norm(),d1.norm()),d2.norm());
        return std::abs( avgRad)*dd - avgRad;        
    }

    template<int i>
    double planedistance( const DROPS::Point3DCL& p)
    {
        double x=p[i]-P.get<DROPS::Point3DCL>("Exp.PosDrop")[i];
        return x;
    }

    static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
    static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
    static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
    static DROPS::RegisterScalarFunction regscacube("CubeDistance", CubeDistance);
    static DROPS::RegisterScalarFunction regscaperellx("perEllipsoidx", PeriodicEllipsoidDistance<0>);
    static DROPS::RegisterScalarFunction regscaperelly("perEllipsoidy", PeriodicEllipsoidDistance<1>);
    static DROPS::RegisterScalarFunction regscaperellz("perEllipsoidz", PeriodicEllipsoidDistance<2>);
}

//========================================================================
//                       Functions for transport
//========================================================================
namespace transpfunctions{
    double tInitialcneg (const DROPS::Point3DCL& , double)
    {
        return P.get<double>("Transp.IniCNeg");
    }

    double tInitialcpos (const DROPS::Point3DCL& , double)
    {
        return P.get<double>("Transp.IniCPos");
    }
    static DROPS::RegisterScalarFunction regscainineg("Initialcneg", tInitialcneg);
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
        return 1. + std::sin( atan2( p[0] - P.get<DROPS::Point3DCL>("Exp.PosDrop")[0], p[2] - P.get<DROPS::Point3DCL>("Exp.PosDrop")[2]));
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
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
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
            std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
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
