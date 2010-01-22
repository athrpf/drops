/// \file serialization.cpp
/// \brief tests serialization of the grid
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/container.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include <fstream>

using namespace DROPS;
Uint rule = 0;

int main (int argc, char** argv)
{

    if (argc>1)
        rule = atoi(argv[1]);

    // Multigrid aufbauen
//    std::ifstream meshfile( "mesh.txt");
//    if (!meshfile) {
//        std::cout << "error while opening mesh file " << "mesh.txt" << "\n";
//        return 1;
//    }
//    DROPS::ReadMeshBuilderCL builder( meshfile);

    TetraBuilderCL builder(rule);
    MultiGridCL mg( builder);
//    MarkAll( mg);
//    mg.Refine();
//    MarkAll( mg);
//    mg.Refine();
    mg.SizeInfo( std::cout);
    std::cout << SanityMGOutCL(mg) << std::endl;

    MGSerializationCL serial (mg, "out-");
    serial.WriteMG();

//    std::ifstream meshfile2( "mesh.txt");
//    DROPS::ReadMeshBuilderCL builder2( meshfile2);
    TetraBuilderCL builder2(rule);
    FileBuilderCL buildmg("out-", &builder2);
    MultiGridCL mg2(buildmg);
    std::cout << "\n \n MultiGrid aus Datei gelesen\n \n";

    mg2.SizeInfo( std::cout);
    std::cout << SanityMGOutCL(mg2) << std::endl;

    MGSerializationCL serial2 (mg2, "neu-");
    serial2.WriteMG();
    return 0;
}
