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
#include "poisson/poissonCoeff.h"
#include "misc/bndmap.h"
#include "misc/params.h"
#include <string>
#define PI 3.14159265

//======================================================================================================================
//                                  static members declarations of poissonCoeff class
//======================================================================================================================

namespace DROPS{
ParamCL PoissonCoeffCL::C_;
double PoissonCoeffCL::dx_;
double PoissonCoeffCL::dy_;
int PoissonCoeffCL::nx_;
int PoissonCoeffCL::ny_;
double PoissonCoeffCL::dt_;
int PoissonCoeffCL::Ref_;
instat_scalar_fun_ptr PoissonCoeffCL::q;
instat_scalar_fun_ptr PoissonCoeffCL::alpha;
instat_scalar_fun_ptr PoissonCoeffCL::f;
instat_scalar_fun_ptr PoissonCoeffCL::Solution;
instat_scalar_fun_ptr PoissonCoeffCL::InitialCondition;
instat_vector_fun_ptr PoissonCoeffCL::Vel;
instat_scalar_fun_ptr PoissonCoeffCL::interface;
}

//======================================================================================================================
//                                  special Functions for poissonP1.cpp and poissonP2.cpp
//======================================================================================================================

extern DROPS::ParamCL P;

double Heat(const DROPS::Point3DCL&, double)
{
    extern DROPS::ParamCL P;
    return P.get<double>("Exp.Heat")/P.get<double>("Exp.Lambda")*1e-3;
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
    static bool first = true;
    static double dx, dy;
    static double Nu, g; // kinematic viscosity, gravitational constant
    static bool adjoint;
    //dirty hack
    if (first){
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy;
        Nu  = P.get<double>("PoissonCoeff.KinematicViscosity");
	    g   = P.get<double>("PoissonCoeff.Gravity"); // this must not be fixed to 9.81 - maybe we use different scalings!!
        std::string adstr ("IA1Adjoint");
        std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
        adjoint = (adstr.compare(IAProbstr) == 0);
        first = false;
    }

    DROPS::Point3DCL ret;
    const double d= p[1]/dy,
        U= g*dy*dy/2/Nu;  //U=gh^2/(2*nu)
    double sgn =1.;
    if(adjoint)
        sgn=-1.;
    ret[0]= sgn*U*(2-d)*d;
    ret[1]=0.;
    ret[2]=0.;

    return ret;
}

/*
static DROPS::RegisterScalarFunction regscaheat("Heat", Heat);
static DROPS::RegisterScalarFunction regscaconstpos("NeuConstPos", NeuConst<0>);
static DROPS::RegisterScalarFunction regscaconstneg("NeuConstNeg", NeuConst<1>);
static DROPS::RegisterScalarFunction regscaexppos("NeuExpPos", NeuExp<0>);
static DROPS::RegisterScalarFunction regscaexpneg("NeuExpNeg", NeuExp<1>);
static DROPS::RegisterScalarFunction regscapoly("NeuPoly", NeuPoly);
static DROPS::RegisterVectorFunction regvecnus("Nusselt", Nusselt);*/
//cast to boost function
static DROPS::RegisterScalarFunction regscaheat("Heat", DROPS::instat_scalar_fun_ptr(Heat));
static DROPS::RegisterScalarFunction regscaconstpos("NeuConstPos", DROPS::instat_scalar_fun_ptr(NeuConst<0>));
static DROPS::RegisterScalarFunction regscaconstneg("NeuConstNeg", DROPS::instat_scalar_fun_ptr(NeuConst<1>));
static DROPS::RegisterScalarFunction regscaexppos("NeuExpPos", DROPS::instat_scalar_fun_ptr(NeuExp<0>));
static DROPS::RegisterScalarFunction regscaexpneg("NeuExpNeg", DROPS::instat_scalar_fun_ptr(NeuExp<1>));
static DROPS::RegisterScalarFunction regscapoly("NeuPoly", DROPS::instat_scalar_fun_ptr(NeuPoly));
static DROPS::RegisterVectorFunction regvecnus("Nusselt", DROPS::instat_vector_fun_ptr(Nusselt));


//======================================================================================================================
//
//               Registered functions in function container which are used in corresponding *Ex.json file
//
//======================================================================================================================

/****************************************
 *  Example 1:                        \n*
 *   - stationary setup               \n*
 *   - constant diffusion             \n*
 *   - no convection                  \n*
 *   - no reaction                    \n*
 *   - given solution (adapted r.h.s) \n*
 ****************************************/
 namespace statpoissonExample{
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
    ///Neumann boundary condition at x=1
    double Neumann(const DROPS::Point3DCL& p, double){
        return -64.*p[0]*p[1]*p[2]*(1-p[1])*(1-p[2]);
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double){
        return 1 + 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
    }
    /// \brief Initial value
    double InitialValue( const DROPS::Point3DCL& , double){
        return 1.;
    }

/*    static DROPS::RegisterScalarFunction regscaq("Example1_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example1_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example1_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example1_Diffusion",    Diffusion   );
    static DROPS::RegisterScalarFunction regscan("Example1_Neumann",      Neumann     );
    static DROPS::RegisterVectorFunction regscav("Example1_Flowfield",    Flowfield   );
    static DROPS::RegisterScalarFunction regscai("Example1_InitialValue", InitialValue);*/
    static DROPS::RegisterScalarFunction regscaq("Example1_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Example1_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Example1_Solution",     DROPS::instat_scalar_fun_ptr(Solution)    );
    static DROPS::RegisterScalarFunction regscaa("Example1_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterScalarFunction regscan("Example1_Neumann",      DROPS::instat_scalar_fun_ptr(Neumann)     );
    static DROPS::RegisterVectorFunction regscav("Example1_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
    static DROPS::RegisterScalarFunction regscai("Example1_InitialValue", DROPS::instat_scalar_fun_ptr(InitialValue));

}//end of namespace


/****************************************
 *  Example 2:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - no convection                  \n*
 *   - no reaction                    \n*
 *   - given solution (adapted r.h.s) \n*
 ****************************************/
 namespace instatpoissonExample{
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
/*    static DROPS::RegisterScalarFunction regscaq("Example2_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example2_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example2_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example2_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example2_Flowfield",    Flowfield   );*/
   static DROPS::RegisterScalarFunction regscaq("Example2_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
   static DROPS::RegisterScalarFunction regscaf("Example2_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
   static DROPS::RegisterScalarFunction regscas("Example2_Solution",     DROPS::instat_scalar_fun_ptr(Solution)    );
   static DROPS::RegisterScalarFunction regscaa("Example2_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
   static DROPS::RegisterVectorFunction regscav("Example2_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );

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
 namespace convdiffExample{

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

/*    static DROPS::RegisterScalarFunction regscaq("Example3_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("Example3_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("Example3_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("Example3_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("Example3_Flowfield",    Flowfield   );
    static DROPS::RegisterScalarFunction regscai("Example3_InitialValue", InitialValue);*/
    static DROPS::RegisterScalarFunction regscaq("Example3_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Example3_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Example3_Solution",     DROPS::instat_scalar_fun_ptr(Solution)    );
    static DROPS::RegisterScalarFunction regscaa("Example3_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterVectorFunction regscav("Example3_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
    static DROPS::RegisterScalarFunction regscai("Example3_InitialValue", DROPS::instat_scalar_fun_ptr(InitialValue));

}//end of namespace
/****************************************
 *  Example 4:                        \n*
 *   - stationary setup               \n*
 *   - quasi-1D setup                 \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 *   - source = 1                     \n*
 ****************************************/
 namespace statSUPGExample{
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
        return 0.01;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double)
    {
        double D = Diffusion(p, 0.);
        return p[0] - (1 - exp(p[0]/D))/(1 - exp(1./D));
    }
/*    static DROPS::RegisterScalarFunction regscaq("SUPG_Reaction",     Reaction    );
    static DROPS::RegisterScalarFunction regscaf("SUPG_Source",       Source      );
    static DROPS::RegisterScalarFunction regscas("SUPG_Solution",     Solution    );
    static DROPS::RegisterScalarFunction regscaa("SUPG_Diffusion",    Diffusion   );
    static DROPS::RegisterVectorFunction regscav("SUPG_Flowfield",    Flowfield   );*/
    static DROPS::RegisterScalarFunction regscaq("SUPG_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("SUPG_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("SUPG_Solution",     DROPS::instat_scalar_fun_ptr(Solution));
    static DROPS::RegisterScalarFunction regscaa("SUPG_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterVectorFunction regscav("SUPG_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
}//end of namespace

/****************************************
 *  Example 5:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace instatSUPGExample{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }

    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 0.01;
    }

    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL&, double){
    //(1, 0, 0)^T
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double t){
    // 1 - (1 + D)e^(-t-y)
        double alpha;
        alpha = Diffusion(p, t);
        double ret=0.;
        ret = 1. - (1 + alpha) * exp(-t-p[1]);
        return ret;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        double alpha;
        alpha = Diffusion(p, t);
        return exp(-t-p[1]) + p[0] - (1 - exp(p[0]/alpha))/(1 - exp(1./alpha));
    }
    static DROPS::RegisterScalarFunction regscaq("instatSUPG_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaa("instatSUPG_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterScalarFunction regscaf("instatSUPG_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("instatSUPG_Solution",     DROPS::instat_scalar_fun_ptr(Solution)    );
    static DROPS::RegisterVectorFunction regscav("instatSUPG_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
}//end of namespace
/****************************************
 *  Example 6:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace AdjointExample{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL&, double){
        DROPS::Point3DCL v(0.);
        v[0] = 1.;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double t){
        return -exp(-t)*exp(-p[0]);
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        return exp(-t)*exp(-p[0]);
    }
    static DROPS::RegisterScalarFunction regscaq("Adjoint_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Adjoint_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Adjoint_Solution",     DROPS::instat_scalar_fun_ptr(Solution)    );
    static DROPS::RegisterScalarFunction regscaa("Adjoint_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterVectorFunction regscav("Adjoint_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
}//end of namespace
/****************************************
 *  Example 7:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace ALEExample1{
    //refH, Mag, paraX and paraT are used to change the free surface functions
    double refH   = 0.2;
    double Mag    = 0.25;
    double paraX  = 6. * PI;
    double paraT  = 10. * PI;
    //Free surface
    double Interface( const DROPS::Point3DCL& p, double t)
    {
        double h= refH + refH * Mag * sin ( paraX * p[0]  + paraT * t );
        return h;
    }
    //Transform the physical to coordinates to reference coordinates
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t)
    {
        DROPS::Point3DCL ref(0.);
        ref[0] = p[0];
        ref[1] = refH * p[1]/Interface(p, t);
        ref[2] = p[2];
        return ref;
    }
    //Gradx, Grady, Grad1 and Grad2 are used for source term
    double Gradx( const DROPS::Point3DCL& ref, double t)   //a=h_x
    {
        return ref[1]* paraX * Mag * cos(paraX * ref[0]  + paraT * t);
    }
    double Grady( const DROPS::Point3DCL& ref, double t)   //b=h_y
    {
        return 1. + Mag * sin(paraX * ref[0]  + paraT * t);
    }
    double Grad1( const DROPS::Point3DCL& ref, double t)   //\nabla_y(a/b)
    {
        double ay=paraX * Mag * cos(paraX * ref[0]  + paraT * t);
        return ay/Grady(ref, t);
    }
    double Grad2( const DROPS::Point3DCL& ref, double t)    //\nabla_x(a/b)
    {
        double ax     = -ref[1]* paraX * paraX * Mag * sin(paraX * ref[0]  + paraT * t);
        double bx     = paraX * Mag * cos(paraX * ref[0]  + paraT * t);
        return (ax*Grady(ref, t) - Gradx(ref, t)*bx)/(Grady(ref,t)*Grady(ref, t));
    }
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double t){
        DROPS::Point3DCL ref = TransBack(p, t);
        DROPS::Point3DCL v(0.);
        v[0] = 1.0;
        v[1] = ref[1] * Mag * paraT * cos(paraX * ref[0]  + paraT * t);
        return v;
    }
    DROPS::Point3DCL ALEFlowfield(const DROPS::Point3DCL&, double){
        DROPS::Point3DCL v(0.);
        v[0] = 1.;
        v[1] = 0.;
        return v;
    }

    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        DROPS::Point3DCL ref = TransBack(p, t);
        return exp(t)*exp(ref[0] + ref[1] + ref[2]);
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double t){
        double alpha;
        alpha = Diffusion(p, t);
        DROPS::Point3DCL ref = TransBack(p, t);
        double a= Gradx(ref,t);
        double b= Grady(ref,t);
        double c= Grad1(ref,t);
        double d= Grad2(ref,t);
        double sol= Solution(p, t);
        double timederiv = sol;
        double conv      = (1. - a/b)*sol;
        double diff1     = -alpha*(1.-d - 2.*a/b + a/b*c + a*a/b/b)*sol;
        double diff2     = -alpha/(b*b)*sol;
        double diff3     = -alpha*sol;
        return timederiv + conv + diff1 +diff2 +diff3;
    }
    static DROPS::RegisterScalarFunction regscaq("TestALE_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)     );
    static DROPS::RegisterScalarFunction regscaf("TestALE_Source",       DROPS::instat_scalar_fun_ptr(Source)       );
    static DROPS::RegisterScalarFunction regscaa("TestALE_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)    );
    static DROPS::RegisterScalarFunction regscaint("TestALE_Interface",  DROPS::instat_scalar_fun_ptr(Interface)    );
    static DROPS::RegisterScalarFunction regscas("TestALE_Solution",     DROPS::instat_scalar_fun_ptr(Solution)     );
    static DROPS::RegisterVectorFunction regscav("TestOrigin_Velocity",  DROPS::instat_vector_fun_ptr(Flowfield)    );
    static DROPS::RegisterVectorFunction regscaalev("TestALE_Velocity",  DROPS::instat_vector_fun_ptr(ALEFlowfield) );
}//end of namespace

 namespace ALE{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    double Interface( const DROPS::Point3DCL& p, double t)
    {
        double h;
        h = 0.2 + 0.05 * sin ( p[0] * 0.1 * PI - 20. * PI * t );
        return h;
    }
    DROPS::Point3DCL TransBack(const DROPS::Point3DCL &p, double t)
    {
        DROPS::Point3DCL ref(0.);
        ref[0] = p[0];
        ref[1] = 0.2 * p[1]/Interface(p, t);
        ref[2] = p[2];
        return ref;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double t){
        DROPS::Point3DCL ref = TransBack(p, t);
        DROPS::Point3DCL ret;
        const double d= ref[1]/0.2,
            U= 200.;  //U=gh^2/(2*nu)
        ret[0]= U*(2-d)*d;
        ret[1]=0.;
        ret[2]=0.;
        return ret;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    double BInter(const DROPS::Point3DCL&, double){
        return 0.01;
    }
    /// \brief Solution ====================================careful

    static DROPS::RegisterScalarFunction regscaq("ALE_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("ALE_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscaint("ALE_Interface",  DROPS::instat_scalar_fun_ptr(Interface)   );
    static DROPS::RegisterVectorFunction regscav("ALE_Velocity",     DROPS::instat_vector_fun_ptr(Flowfield)   );
    static DROPS::RegisterScalarFunction regscas("ALE_Inter",         DROPS::instat_scalar_fun_ptr(BInter)     );
}//end of namespace


//======================================================================================================================
//
//               Examples for incremental approach for IA1adjoint, IA1adjoint2, IA2sensitivity, IA2Grad
//
//======================================================================================================================
 namespace PoissonAdjoint{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL&, double){
        DROPS::Point3DCL v(0.);
        v[0] = 1.;
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double t){
        return -exp(-t)*exp(-p[0]);
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        return exp(-t)*exp(-p[0]);
    }
    static DROPS::RegisterScalarFunction regscaq("Adjoint_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Adjoint_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Adjoint_Solution",     DROPS::instat_scalar_fun_ptr(Solution));
    static DROPS::RegisterScalarFunction regscaa("Adjoint_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterVectorFunction regscav("Adjoint_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
}//end of namespace
 namespace PoissonAdjoint2{
    double Pi =3.1415926;
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection: constant flow in x direction
    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double){
        DROPS::Point3DCL v(0.);
        v[0] = (2.-p[1])*p[1];
        return v;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double t){
        double time=1.-t;
        double t_deriv = (1.-cos(2*Pi*p[0]))*sin(Pi/2.*(1.-p[1]))*Pi/2.*sin(Pi/2.*time);
        double conv    = - Flowfield(p, time)[0]*2*Pi*sin(2*Pi*p[0])*sin(Pi/2.*(1.-p[1]))*cos(Pi/2.*time);
        double diff    = - 4*Pi*Pi*cos(2*Pi*p[0])*sin(Pi/2.*(1.-p[1]))*cos(Pi/2.*time) + (1.-cos(2*Pi*p[0]))*Pi*Pi/4.*sin(Pi/2.*(1.-p[1]))*cos(Pi/2.*time);
        return  t_deriv+conv+diff;
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double t)
    {
        double time=1.-t;
        return (1.-cos(2*Pi*p[0]))*sin(Pi/2.*(1.-p[1]))*cos(Pi/2.*time);
    }
    static DROPS::RegisterScalarFunction regscaq("Adjoint2_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Adjoint2_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Adjoint2_Solution",     DROPS::instat_scalar_fun_ptr(Solution));
    static DROPS::RegisterScalarFunction regscaa("Adjoint2_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
    static DROPS::RegisterVectorFunction regscav("Adjoint2_Flowfield",    DROPS::instat_vector_fun_ptr(Flowfield)   );
}//end of namespace

 namespace IA2Sensi{

    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double){
        double pi= 3.1415926;
        return -sin(2*pi*p[0]);
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double)
    {
        double pi= 3.1415926;
        return sin(2*pi*p[0]);
    }
    static DROPS::RegisterScalarFunction regscaq("Sensi_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Sensi_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Sensi_Solution",     DROPS::instat_scalar_fun_ptr(Solution));
    static DROPS::RegisterScalarFunction regscaa("Sensi_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
   }//end of namespace

 namespace IA2Grad{
    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL& p, double){
        return 2*exp(-p[0]);
    }
    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1.0;
    }
    /// \brief Solution
    double Solution( const DROPS::Point3DCL& p, double)
    {
        return p[0]*exp(-p[0]);
    }
    static DROPS::RegisterScalarFunction regscaq("Grad_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)    );
    static DROPS::RegisterScalarFunction regscaf("Grad_Source",       DROPS::instat_scalar_fun_ptr(Source)      );
    static DROPS::RegisterScalarFunction regscas("Grad_Solution",     DROPS::instat_scalar_fun_ptr(Solution));
    static DROPS::RegisterScalarFunction regscaa("Grad_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)   );
}//end of namespace

/****************************************
 *  Balakotaiah:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
 namespace Balakotaiah{
/***********************************************************************************************************************************/
//added
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

int getn(std::string dfile)
{
    int n=0;
    double tmp;
    std::ifstream fd;
    fd.open(dfile.c_str());
    if(fd.fail())
	std::cout << dfile << " not available" << std::endl;

    while(!fd.eof()){
        if(fd.good()){
	   fd >> tmp;
	   n++;
	}
    }

    fd.close();
    return n;
}

int countTimesteps(const int starttime, const int stoptime, const int timeinc)
{
    int tmp = 0;
    for(int i=starttime; i<=stoptime; i+=timeinc)
        tmp++;

    return tmp;
}

void getData(double *d, std::string dfile, std::string time)
{
    std::ifstream fd;
    std::string file;
    file = dfile + time + ".dat";
    fd.open(file.c_str());
    int i=0;

    while(!fd.eof()){
        if(fd.good()){
            fd >> d[i];
            i++;
        }
    }

    fd.close();
}

void uvh(double* res, int t, double x, double y)
{
    static bool     init=false;

    static int      n;

    static gsl_spline ** hspline,
                      ** qspline;

    static gsl_interp_accel * acc;

    double          hres[2],
                    qres[2];

    const std::string hfile="Data/longh",
                      qfile="Data/longq";

    const int       starttime=161,
                    stoptime=673,
                    timeinc=1;

    if(!init){
        std::string file;
        std::stringstream out;
        out << starttime;
        file=hfile+out.str()+".dat";
        n=getn(file);

        int tsteps = countTimesteps(starttime,stoptime,timeinc);
	acc = gsl_interp_accel_alloc();
        hspline = new gsl_spline*[tsteps];
        qspline = new gsl_spline*[tsteps];

        double * h,
               * q,
               * xa;

        h = new double[n];
        q = new double[n];
        xa = new double[n];

        for(int j=0; j<n; j++)
            xa[j] = j+1.;

	std::cout << "Getting Data... ";
        for(int i=0; i<tsteps; i++){
            std::stringstream time;
            time << starttime + i*timeinc;
            getData(h, hfile, time.str());
            hspline[i] = gsl_spline_alloc(gsl_interp_cspline, n);
            gsl_spline_init(hspline[i], xa, h, n);
            getData(q, qfile, time.str());
            qspline[i] = gsl_spline_alloc(gsl_interp_cspline, n);
            gsl_spline_init(qspline[i], xa, q, n);
        }
	std::cout << " done" << std::endl;


        init=true;

        delete[] h;
        delete[] q;
        delete[] xa;
    }

    hres[0] = gsl_spline_eval(hspline[(int)((t+0.1)/timeinc)], x, acc);
    hres[1] = gsl_spline_eval_deriv(hspline[(int)((t+0.1)/timeinc)], x, acc);
    qres[0] = gsl_spline_eval(qspline[(int)((t+0.1)/timeinc)], x, acc);
    qres[1] = gsl_spline_eval_deriv(qspline[(int)((t+0.1)/timeinc)], x, acc);

    /*
    res[0]=3*qres[0]/pow(hres[0],3)*(hres[0]*y-y*y/2);
    res[1]=-(1.5*qres[0]/pow(hres[0],4)*hres[1]-0.5*qres[1]/pow(hres[0],3))*pow(y,3)
            -(1.5*qres[1]/pow(hres[0],2)-3*qres[0]/pow(hres[0],3)*hres[1])*y*y;
    res[2]=hres[0];
    */
    double h1 = hres[0];
    double h2 = h1*h1;
    double h3 = h2*h1;
    double h4 = h2*h2;
    double y2 = y*y;
    double y3 = y2*y;
    res[0]=3*qres[0]/h3*(hres[0]*y-y2/2);
    res[1]=-(1.5*qres[0]/h4*hres[1]-0.5*qres[1]/h3)*y3
            -(1.5*qres[1]/h2-3*qres[0]/h3*hres[1])*y*y;
    res[2]=hres[0];
}
/*****************************************************************************************************************/

    double Interface( const DROPS::Point3DCL& p, double t)
    {
        double res[3];
        //std::cout << "t: " << t  << "x: " << p[0] << "y: " << p[1] << std::endl;
        uvh(res,(int)(t+0.1),p[0],p[1]);
        return res[2];
    }

    /// \brief Reaction: no reaction
    double Reaction(const DROPS::Point3DCL&, double){
        return 0.0;
    }

    /// \brief Diffusion
    double Diffusion(const DROPS::Point3DCL&, double){
        return 1e-9;
    }

    /// \brief Initial value
    double InitialValue( const DROPS::Point3DCL& , double){
        return 1e-5;
    }

    DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double t){
      double res[3];
//std::cout << "t: " << t << std::endl;
	uvh(res,(int)(t+0.1),p[0],p[1]);
        DROPS::Point3DCL v(0.);
        v[0] = res[0];
        v[1] = res[1];
	v[2] = 0.;
        return v;
    }

    /// \brief Solution
    double Solution( const DROPS::Point3DCL&, double)
    {
        return 1e-3;
    }

    /// \brief Right-hand side
    double Source(const DROPS::Point3DCL&, double){
        return 0.;
    }
    static DROPS::RegisterScalarFunction regscaq("Balakotaiah_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)     );
    static DROPS::RegisterScalarFunction regscaf("Balakotaiah_Source",       DROPS::instat_scalar_fun_ptr(Source)       );
    static DROPS::RegisterScalarFunction regscaa("Balakotaiah_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)    );
    static DROPS::RegisterScalarFunction regscaint("Balakotaiah_Interface",  DROPS::instat_scalar_fun_ptr(Interface)    );
    static DROPS::RegisterScalarFunction regscas("Balakotaiah_Solution",     DROPS::instat_scalar_fun_ptr(Solution)     );
    static DROPS::RegisterVectorFunction regscav("Balakotaiah_Velocity",     DROPS::instat_vector_fun_ptr(Flowfield)    );
    static DROPS::RegisterScalarFunction regscai("Balakotaiah_InitialValue", DROPS::instat_scalar_fun_ptr(InitialValue) );
}//end of namespace
