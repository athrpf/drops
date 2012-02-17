/// \file  ale.h
/// \brief classes that move the grids for ale method
/// \author LNM RWTH Aachen: Liang Zhang; SC RWTH Aachen;

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

#ifndef DROPS_ALE_H
#define DROPS_ALE_H

#include "misc/container.h"
#include "misc/params.h"
#include "misc/bndmap.h"
#include "poisson/poissonCoeff.h"
#include <sstream>

namespace DROPS{
    
class ALECL{
///A class to handle ALE method for free surface one phase scalar problem;
    private:
    bool IfALE_;
    ParamCL Para_;
    MultiGridCL& mg_;  
    double dt_;
    double Ly_;   

  public:
    //free surface function
    instat_scalar_fun_ptr interface_;  

    ALECL(ParamCL P, MultiGridCL& mg):
    IfALE_(P.get<int>("ALE.wavy")),
    Para_(P), mg_(mg), dt_(Para_.get<double>("Time.StepSize"))
    {
        double Lx_, Lz_;
        std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
        size_t idx_;
        while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx_]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> Lx_ >> Ly_ >> Lz_ ;

        DROPS::InScaMap & scamap = DROPS::InScaMap::getInstance();
        interface_ = scamap[P.get<std::string>("ALE.Interface")];
    }
    bool GetALE() {return IfALE_;}
    //Initialize the grids
    void InitGrid();    
    //Scale the grids according to the free surface functions
    void MovGrid(double t);
};

} 
#endif