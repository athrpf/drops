/// \file poissonCoeff.cpp
/// \brief boundary and source functions for the poisson-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Joerg Grande, Thorolf Schulte, Liang Zhang

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

//========================================================================
//          Functions for poissonP1.cpp and drops_statP2.cpp
//========================================================================

double Heat(const DROPS::Point3DCL&, double)
{
    extern DROPS::ParamCL P;
    return P.get<double>("Exp.Heat")/P.get<double>("Exp.Lambda")*1e-3;
}

/// \brief Interface condition for one phase film
double Interface(const DROPS::Point3DCL&, double)
{
    return 0.02;    
}

double InitFilm(const DROPS::Point3DCL&, double)
{
    return 1.0e-3;    
}


/// boundary description of a neumann problem
// uses constant function f = (-1)^seg *4.0
template<int sel>
double NeuConst( const DROPS::Point3DCL& , double ){return std::pow(-1.,sel)*4.0; }

/// boundary description of a neumann problem
// uses exp-function f =  e^(t)* e^(px + py + pz)
template<int sel>
double NeuExp( const DROPS::Point3DCL& p, double t){return std::pow(-1.,sel)*std::exp(t)*std::exp(p[0]+p[1]+p[2]); }

/// boundary description of a neumann problem
// uses polynomial function
double NeuPoly( const DROPS::Point3DCL& p, double ){return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]);}

/// \brief Nusselt velocity profile for flat film
DROPS::Point3DCL Nusselt(const DROPS::Point3DCL& p, double)
{
    extern DROPS::ParamCL P;

    static bool first = true;
    static double dx, dy;
    static double Rho, Mu;   //density, viscosity
    //dirty hack
    if (first){
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy;
        Rho = P.get<double>("Exp.Rho");
        Mu  = P.get<double>("Exp.Mu");
        first = false;
    }

    DROPS::Point3DCL ret;
    const double d= p[1]/dy,
        U= Rho*9.81*dy*dy/2/Mu;  //U=gh^2/(2*nu)      
    ret[0]= U*(2-d)*d;                       
    ret[1]=0.;
    ret[2]=0.;

    return ret;
}

static DROPS::RegisterScalarFunction regscaheat("Heat", Heat);
static DROPS::RegisterScalarFunction regscainterf("Interface", Interface);
static DROPS::RegisterScalarFunction regscainitf("InitFilm", InitFilm);
static DROPS::RegisterScalarFunction regscaconstpos("NeuConstPos", NeuConst<0>);
static DROPS::RegisterScalarFunction regscaconstneg("NeuConstNeg", NeuConst<1>);
static DROPS::RegisterScalarFunction regscaexppos("NeuExpPos", NeuExp<0>);
static DROPS::RegisterScalarFunction regscaexpneg("NeuExpNeg", NeuExp<1>);
static DROPS::RegisterScalarFunction regscapoly("NeuPoly", NeuPoly);
static DROPS::RegisterVectorFunction regvecnus("Nusselt", Nusselt);

/****************************************
 *  Example 1:                        \n*
 *   - stationary setup               \n*
 *   - constant diffusion             \n*
 *   - no convection                  \n*
 *   - no reaction                    \n*
 *   - given solution (adapted r.h.s) \n*
 ****************************************/
 namespace PoissonExample1{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection: no convection
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL&, double){ 
        return DROPS::Point3DCL(0.); 
    } 
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double){
        return 128*(p[1]*p[2]*(1.-p[1])*(1.-p[2])
                + p[0]*p[2]*(1.-p[0])*(1.-p[2])
                + p[0]*p[1]*(1.-p[0])*(1.-p[1]));
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double){
        return 1 + 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
    }
    /// \brief Initial value
    double InitialValue( const DROPS::Point3DCL& , double){
        return 1.;
    }

    static DROPS::RegisterScalarFunction regscaq("Example1_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example1_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example1_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example1_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example1_Flowfield",    Flowfield   );
    static DROPS::RegisterScalarFunction regscai("Example1_InitialValue", InitialValue);

}//end of namespace


/****************************************
 *  Example 2:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - no convection                  \n*
 *   - no reaction                    \n*
 *   - given solution (adapted r.h.s) \n*
 ****************************************/
 namespace PoissonExample2{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection: no convection
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL&, double){ 
        return DROPS::Point3DCL(0.); 
    } 
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double t){
        return (-2.0*std::exp(t)*std::exp(p[0]+p[1]+p[2])); 
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t){
        return (std::exp(t)*std::exp(p[0]+p[1]+p[2])); 
    }

    static DROPS::RegisterScalarFunction regscaq("Example2_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example2_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example2_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example2_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example2_Flowfield",    Flowfield   );

}//end of namespace


/****************************************
 *  Example 3:                        \n*
 *   - stationary setup               \n*
 *   - quasi-1D setup                 \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 *   - source = 1                     \n*
 ****************************************/
 namespace PoissonExample3{

    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL&, double){ 
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        return v; 
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL&, double){
        return 1.0; 
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double){
        return 1.0 + p[0] + (1-exp(p[0]))/exp(1.0); 
    }
    /// \brief Initial value
    double InitialValue( const DROPS::Point3DCL& , double){
        return 1.0;
    }

    static DROPS::RegisterScalarFunction regscaq("Example3_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example3_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example3_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example3_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example3_Flowfield",    Flowfield   );
    static DROPS::RegisterScalarFunction regscai("Example3_InitialValue", InitialValue);

}//end of namespace


double Solution( const DROPS::Point3DCL& p, double=0.0){
    return 1 + p[0] + (1-exp(p[0]))/exp(1.0);
}

double diffSolution( const DROPS::Point3DCL& p, double=0.0){
    return 3.0/2.0 - (p[0]-1)*(p[0]-1)/2.0;
}

static DROPS::RegisterScalarFunction regscasol("OurSolution", Solution);
static DROPS::RegisterScalarFunction regscadsol("OurDiffSolution", diffSolution);




