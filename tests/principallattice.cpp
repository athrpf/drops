/// \file combinatorialcut.cpp
/// \brief tests the PrincipalLattice-class
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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

#include "geom/principallattice.h"
#include <iostream>

int main()
{
    try {
        for (int i= 1; i <= 4; ++i) {
            const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( i);
            std::cout << "=======================================" << lat.num_intervals() << ' ' << lat.num_vertexes() << " " << lat.num_tetras() << std::endl;
            for (DROPS::PrincipalLatticeCL::const_vertex_iterator v= lat.vertex_begin(), end= lat.vertex_end(); v != end; ++v) {
                std::cout << lat.num_intervals()*(*v) << std::endl;
            }
            std:: cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
            for (DROPS::PrincipalLatticeCL::const_tetra_iterator v= lat.tetra_begin(), end= lat.tetra_end(); v != end; ++v) {
                std::cout << (*v)[0] << ' '  << (*v)[1] << ' ' << (*v)[2] << ' ' << (*v)[3] << ' ' << std::endl;
            }
        }
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}