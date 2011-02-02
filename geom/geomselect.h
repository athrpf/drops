/// \file geomselect.h
/// \brief offers build/create routines for some standard domains
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Christoph Lehrenfeld, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef GEOMSELECT_H_
#define GEOMSELECT_H_

#include "geom/multigrid.h"
#include "misc/bndmap.h"

namespace DROPS
{

/// \brief creates a standard domain
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param meshfile_name name of the gambit-meshfile or list of integers/double which determine the L-/B-/Cavity-/Brick domain, e.g. 1x1x1@2x3x2
 * @param GeomType type of the geometry: 0: GAMBIT mesh, 1: brick, 2: cavity, 3: l-brick, 4: b-brick
 * @param deserialization_file file that contains the (serialized) grid, "none": use standard builder
 * @param r_inlet inlet radius
*/
void BuildDomain( MultiGridCL* &mgp, const std::string& meshfile_name, int GeomType,
        const std::string& deserialization_file, double& r_inlet);


/// \brief determines boundary conditions
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param bnddata a pointer of the boundary data object
 * @param bnd_type function pointer which is used to construct the boundary description
 * @param bnd_type specify boundary types, e.g. 0!0!0!2!0!0, Dir0BC= 0, DirBC= 2
 * @param bnd_funcs specify boundary functions
 */
template< class BoundaryT>
void BuildBoundaryData( MultiGridCL* &mgp, BoundaryT* &bnddata,
        const std::string& bnd_type, const std::string& bnd_funcs);

} // end of namespace drops

#include "geom/geomselect.tpp"

#endif /* GEOMSELECT_H_ */
