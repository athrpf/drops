/// \file tetrabuildertest.cpp
/// \brief tests tetra builder with refinement rules
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

#include "geom/builder.h"

DROPS::Uint Rule(DROPS::Uint r)
{
    return r < 64 ? r : 127;
}

int TestReMark()
{
    int ret= 0;
    for (DROPS::Uint i= 0; i<=64; ++i) {
        for (DROPS::Uint j= 0; j<=64; ++j) {
            DROPS::TetraBuilderCL tet( Rule( i));
            DROPS::MultiGridCL mg( tet);
            tet.BogoReMark( mg, Rule( j));
            if ( mg.GetTetrasBegin( 0)->GetRefRule() != Rule( j)) {
                std::cout << "Aerger mit Regel " << Rule( i) << std::endl;
                ++ret;
            }
        }
    }
    return ret;
}

int main ()
{
  try {
    int ret= 0;
    std::cout << "Testing building single rules." << std::endl;
    for (DROPS::Uint i= 0; i<64; ++i) {
        DROPS::TetraBuilderCL brick( i);
        DROPS::MultiGridCL mg( brick);
        if ( mg.GetTetrasBegin( 0)->GetRefRule() != i) {
            std::cout << "Aerger mit Regel " << i << std::endl;
            ++ret;
        }
    }

    std::cout << "Testing remarking between different rules." << std::endl;
    ret+= TestReMark();
    return ret;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
