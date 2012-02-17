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
#include "num/discretize.h"
#include "misc/container.h"
#include "misc/params.h"
#include "poisson/ale.h"
#include <sstream>


namespace DROPS{
    
void ALECL::InitGrid()
{
    Point3DCL New_Coord(0.);
    int lvl= mg_.GetLastLevel();
    for (MultiGridCL::TriangVertexIteratorCL sit= mg_.GetTriangVertexBegin(lvl), send= mg_.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        Point3DCL Old_Coord(0.);
        Old_Coord = sit->GetCoord();
        New_Coord[0] = Old_Coord[0];
        New_Coord[1] = Old_Coord[1] * interface_(Old_Coord, 0.)/Ly_;
        New_Coord[2] = Old_Coord[2];
        sit->ChangeCoord(New_Coord);
    }
}

void ALECL::MovGrid(double t)
{
    Point3DCL New_Coord(0.);
    int lvl= mg_.GetLastLevel();  
    for (MultiGridCL::TriangVertexIteratorCL sit= mg_.GetTriangVertexBegin(lvl), send= mg_.GetTriangVertexEnd(lvl);
        sit != send; ++sit)
    {
        Point3DCL Old_Coord(0.);
        Old_Coord = sit->GetCoord();
        New_Coord[0] = Old_Coord[0];
        New_Coord[1] = Old_Coord[1] * interface_(Old_Coord, t + dt_)/interface_(Old_Coord, t);
        New_Coord[2] = Old_Coord[2];
        sit->ChangeCoord(New_Coord);
    }
}

} 