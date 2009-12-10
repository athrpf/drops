/// \file surfacetension.cpp
/// \brief compute the interfacial tension
/// \author LNM RWTH Aachen: Hieu Nguyen; SC RWTH Aachen:

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

#include "levelset/surfacetension.h"

namespace DROPS
{

double SurfaceTensionCL::sigma_c(double c) const
{
    if (c<0.) c=0.;
    if (c>1.) c=1.;
    double x= 1./(1.-C_[4]*c), y= c-cp_, z=y*y;
    return 0.001*(C_[0] + C_[1]*y+ C_[2]*z +C_[3]*y*z)*x;
}

double SurfaceTensionCL::dsigma_dc(double c) const
{
    if (c<0.) c=0.;
    if (c>1.) c=1.;
    double x= 1./(1.-C_[4]*c), y= c-cp_, z= y*y;
    return 0.001*x*((C_[1]+ 2*C_[2]*y +3* C_[3]*z) + C_[4]*x*(C_[0] + C_[1]*y+ C_[2]*z +C_[3]*y*z));
}

void SurfaceTensionCL::ComputeSF(const TetraCL& t, const BaryCoordCL * const p,
                          Quad5_2DCL<>& qsigma, Quad5_2DCL<Point3DCL>& qgradsigma) const
{
    switch (input_)
    {
        case Sigma_X: {
                          qsigma.assign(t, p, sigma_);
                          qgradsigma.assign(t, p, grad_sigma_);
                      }
                      break;
        case Sigma_C: {
                          LocalP1CL<> p1_c(t, *c_, cBnd_, time_);
                          LocalP2CL<> p2_c(p1_c);
                          Quad5_2DCL<> q5_c(p2_c, p);
                          Point3DCL gradc(0.), G[4];
                          double det;
                          P1DiscCL::GetGradients( G, det, t);
                          for (Uint i=0; i < 4; ++i)
                              gradc+= p1_c[i]*G[i];
                          for (Uint i= 0; i < Quad5_2DDataCL::NumNodesC; ++i) {
                              qsigma[i]= sigma_c(q5_c[i]);
                              qgradsigma[i]= dsigma_dc(q5_c[i])*gradc;
                          }
                      }
                      break;
        case Sigma_S: {
                          std::cerr<<"SurfaceTensionCL::ComputeSF: Sorry. Sigma_S not implemented yet."<<std::endl;
                      }
                      break;
         default:     throw DROPSErrCL("SurfaceTensionCL::ComputeSF: unknown method\n");
    }
    return ;
}

} // end of namespace DROPS
