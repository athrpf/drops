/// \file splitboundary.cpp
/// \brief tests splitting of boundary elements
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "out/output.h"
#include <fstream>
#include <cmath>
#include <cstdio>

int main()
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    DROPS::MultiGridCL mg(brick);

    mg.SplitMultiBoundaryTetras();
    mg.MakeConsistentNumbering();
    std::cout << "Checking Sanity...\n" << std::flush;
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    std::cout << "...ok" << std::endl;

    std::ofstream geomview("splitboundary.geo");
    DROPS::GeomMGOutCL ggg( mg, -1, false, 3.0);
    geomview << ggg << std::flush;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
