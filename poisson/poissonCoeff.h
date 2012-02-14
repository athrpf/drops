/// \file poissonCoeff.h
/// \brief PoissonCoeffCL the Coefficient-Function Container
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Liang Zhang, Thorolf Schulte; SC RWTH Aachen:
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

#ifndef POISSONCOEFF_H_
#define POISSONCOEFF_H_

#include "misc/container.h"
#include "misc/params.h"
#include "misc/bndmap.h"
#include <sstream>


namespace DROPS{


/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/

class PoissonCoeffCL
{
  private:
    static ParamCL C_;
    static double dx_, dy_;
    static int    nx_, ny_;
    static double dt_;

  public:
    //reaction
    static instat_scalar_fun_ptr q;
    //diffusion
    static instat_scalar_fun_ptr alpha;
    //source
    static instat_scalar_fun_ptr f;
    //initial condition
    static instat_scalar_fun_ptr InitialCondition;
    //solution
    static instat_scalar_fun_ptr Solution;
    //velocity
    static instat_vector_fun_ptr Vel;
    //Free interface function
    static instat_scalar_fun_ptr interface;  
    
    PoissonCoeffCL( ParamCL& P){
        C_=P;
        int nz_;
        double dz_;
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx_ >> dy_ >> dz_ >> nx_ >> ny_ >> nz_;
        dt_=P.get<double>("Time.StepSize");  //step size used in ALEVelocity
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        q = scamap[P.get<std::string>("PoissonCoeff.Reaction")];
        alpha = scamap[P.get<std::string>("PoissonCoeff.Diffusion")];
        f = scamap[P.get<std::string>("PoissonCoeff.Source")];
        Solution = scamap[P.get<std::string>("PoissonCoeff.Solution")];
        InitialCondition = scamap[P.get<std::string>("PoissonCoeff.InitialVal")];
        DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
        if(P.get<int>("PoissonCoeff.Convection")==0)
            Vel = vecmap["ZeroVel"];
        else
            Vel = vecmap[P.get<std::string>("PoissonCoeff.Flowfield")];
        interface = scamap[P.get<std::string>("ALE.Interface")];
    }
    PoissonCoeffCL( ParamCL& P, instat_scalar_fun_ptr diffusion, instat_scalar_fun_ptr source, instat_scalar_fun_ptr init){
        C_=P;
        int nz_;
        double dz_;
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx_ >> dy_ >> dz_ >> nx_ >> ny_ >> nz_;
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        q = scamap["Zero"];
        alpha = diffusion;
        f = source;
        //f = scamap[P.get<std::string>("PoissonCoeff.Source")];
        if (P.get<int >("Poisson.SolutionIsKnown")==0)
        Solution = scamap["Zero"];
        else
        Solution = scamap[P.get<std::string>("PoissonCoeff.Solution")];
        InitialCondition = init;
        DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
        if(P.get<int>("PoissonCoeff.Convection", 1)==0)
            Vel = vecmap["ZeroVel"];
        else
            Vel = vecmap[P.get<std::string>("PoissonCoeff.Flowfield")];
    }
    
    static Point3DCL ALEVelocity(const DROPS::Point3DCL& p, double t)
    {
        double eps =1.0e-7;
        DROPS::Point3DCL ret;
        ret  = Vel(p, t);
        ret[1] -= p[1]/interface(p, t)*(interface(p, t+eps)-interface(p, t))/eps;  //y/h(p,t)*h_p'(t)
        return ret;
    }
};

}//end of namespace
#endif
