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
#include <sstream>


namespace DROPS{

  //typedef double    (*scalar_fun_ptr)       (const Point3DCL&);
  //typedef Point3DCL (*vector_fun_ptr)       (const Point3DCL&);
//typedef double    (*instat_scalar_fun_ptr)(const Point3DCL&, double);
//typedef Point3DCL (*instat_vector_fun_ptr)(const Point3DCL&, double);


/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/

template<class ParamsT>
class PoissonCoeffCL
{
  private:
    static ParamCL C_;
    static double dx_, dy_;
    static int    nx_, ny_;
    static int    Ref_;   //Times of refinements

  public:
    //reaction
    static instat_scalar_fun_ptr q;

    static instat_scalar_fun_ptr alpha;

    //source
    static instat_scalar_fun_ptr f;
    //initial condition
    static instat_scalar_fun_ptr InitialCondition;
    //solution
    static instat_scalar_fun_ptr Solution;
    //velocity
    static instat_vector_fun_ptr Vel;

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
        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        q = scamap[P.get<std::string>("PoissonCoeff.Reaction")];
        alpha = scamap[P.get<std::string>("PoissonCoeff.Diffusion")];
        f = scamap[P.get<std::string>("PoissonCoeff.Source")];
        Solution = scamap[P.get<std::string>("PoissonCoeff.Solution")];
        InitialCondition = scamap[P.get<std::string>("PoissonCoeff.InitialVal")];
        DROPS::InVecMap & vecmap = DROPS::InVecMap::getInstance();
        Vel = vecmap[P.get<std::string>("PoissonCoeff.Flowfield")];
        Ref_=P.get<int>("DomainCond.RefineSteps");
    }
    //Only used for flat film case
    static double h_Value()
    {//mesh size in flow direction
        double h=dx_/(nx_*std::pow(2, Ref_));
        return h;
    }
    static double Sta_Coeff(const DROPS::Point3DCL& p, double t)
    {//Stabilization coefficient
        double h  =h_Value();
        double Pec=0.;
        Pec=Vel(p, t).norm()*h/(2.*alpha(p, t));  //compute mesh Peclet number
        if (Pec<=1)
            return 0.0;
        else
            return h/(2.*Vel(p, t).norm())*(1.-1./Pec);
    }
    static void Show_Pec()
    {
        double U=9.81*1.e3*dy_*dy_/(2*C_.get<double>("Exp.Mu"));

        const char line[] ="----------------------------------------------------------------------------------\n";
        std::cout<<line<<"The estimate of Peclet number is: "<<U*h_Value()/(2.*C_.get<double>("PoissonCoeff.Diffusion"))<<std::endl;
    }
};

template<class ParamsT>
ParamCL PoissonCoeffCL<ParamsT>::C_;

template<class ParamsT>
double PoissonCoeffCL<ParamsT>::dx_;

template<class ParamsT>
double PoissonCoeffCL<ParamsT>::dy_;
template<class ParamsT>
int PoissonCoeffCL<ParamsT>::nx_;
template<class ParamsT>
int PoissonCoeffCL<ParamsT>::ny_;
template<class ParamsT>
int PoissonCoeffCL<ParamsT>::Ref_;

template<class ParamsT>
instat_scalar_fun_ptr PoissonCoeffCL<ParamsT>::q;
template<class ParamsT>
instat_scalar_fun_ptr PoissonCoeffCL<ParamsT>::alpha;
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
