/// \file poissonCoeff.h
/// \brief PoissonCoeffCL the Coefficient-Function Container
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Liang Zhang; SC RWTH Aachen:
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
#include <sstream>

namespace DROPS{

typedef double    (*scalar_fun_ptr)       (const Point3DCL&);
typedef Point3DCL (*vector_fun_ptr)       (const Point3DCL&);
typedef double    (*instat_scalar_fun_ptr)(const Point3DCL&, double);
typedef Point3DCL (*instat_vector_fun_ptr)(const Point3DCL&, double);  
  
  
/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/

template<class ParamsT>
class PoissonCoeffCL
{
  private:
    static ParamsT C_;
    static double dx_, dy_;

  public:
    //reaction
    static instat_scalar_fun_ptr q;
    
    static double alpha;

    //source
    static instat_scalar_fun_ptr f;
    //initial condition
    static instat_scalar_fun_ptr InitialCondition;
    //solution
    static instat_scalar_fun_ptr Solution;
    //velocity
    static instat_vector_fun_ptr Vel;
  
    PoissonCoeffCL( ParamCL& P){
        int nx_,ny_,nz_;
        double dz_;
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx_ >> dy_ >> dz_ >> nx_ >> ny_ >> nz_;
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        q = scamap[P.get<std::string>("PoissonCoeff.Reaction")];
        alpha = P.get<double>("PoissonCoeff.Diffusion");
        f = scamap[P.get<std::string>("PoissonCoeff.Source")];
        Solution = scamap[P.get<std::string>("PoissonCoeff.Solution")];
        InitialCondition = scamap[P.get<std::string>("PoissonCoeff.InitialVal")];
        DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
        Vel = vecmap[P.get<std::string>("PoissonCoeff.Flowfield")];
    }
    
};

template<class ParamsT>
ParamsT PoissonCoeffCL<ParamsT>::C_;

template<class ParamsT>
double PoissonCoeffCL<ParamsT>::dx_;

template<class ParamsT>
double PoissonCoeffCL<ParamsT>::dy_;

template<class ParamsT>
instat_scalar_fun_ptr PoissonCoeffCL<ParamsT>::q;
template<class ParamsT>
double PoissonCoeffCL<ParamsT>::alpha;
template<class ParamsT>
instat_scalar_fun_ptr PoissonCoeffCL<ParamsT>::f;
template<class ParamsT>
instat_scalar_fun_ptr PoissonCoeffCL<ParamsT>::Solution;
template<class ParamsT>
instat_scalar_fun_ptr PoissonCoeffCL<ParamsT>::InitialCondition;
template<class ParamsT>
instat_vector_fun_ptr PoissonCoeffCL<ParamsT>::Vel;

}//end of namespace
#endif
