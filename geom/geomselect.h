/// \file geomselect.h
/// \brief offers build/create routines for some standard domains
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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
#include "stokes/instatstokes2phase.h"
#include "poisson/poisson.h"

namespace DROPS
{

/// \brief creates a standard domain
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param meshfile_name name of the gambit-meshfile or list of integers/double which determine the L-/B-/Cavity-/Brick domain, e.g. 1x1x1@2x3x2
 * @param GeomType type of the geometry: 0: GAMBIT mesh, 1: brick, 2: cavity, 3: l-brick, 4: b-brick
 * @param deserialization_file file that contains the (serialized) grid, "none": use standard builder
 * @param r_inlet inlet radius
 * @param BC description of the boundary
*/
void BuildDomain( MultiGridCL* &mgp, const std::string& meshfile_name, int GeomType,
        const std::string& deserialization_file, double& r_inlet, std::vector<BndCondT>& BC);


/// \brief determines boundary conditions for some standard poisson problems
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param bnddata a pointer of the boundary data object
 * @param fun scalar function pointer which is used to construct the boundary description
 * @param GeomType type of the geometry: 0: GAMBIT mesh, 1: brick, 2: cavity, 3: l-brick, 4: b-brick
 * @param bnd_type type of the boundary: 0: measuring cell, 1: hom. Dirichlet, 2-5: different Neumann conditions, a.t.m. only valid for a brick domain
 * @param BC description of the boundary
 */
void BuildPoissonBoundaryData( MultiGridCL* &mgp, PoissonBndDataCL* &bnddata,
        instat_scalar_fun_ptr fun, int GeomType, int bnd_type, std::vector<BndCondT>& BC);

/// \brief determines boundary conditions for some standard stokes problems
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param bnddata a pointer of the boundary data object
 * @param inflow vector function pointer which is used to construct the boundary description, e.g. inflow velocity
 * @param GeomType type of the geometry: 0: GAMBIT mesh, 1: brick, 2: cavity, 3: l-brick, 4: b-brick
 * @param bnd_type type of the boundary: 0: measuring cell, 1: hom. Dirichlet, 2-5: different Neumann conditions, a.t.m. only valid for a brick domain
 * @param BC description of the boundary
 */
void BuildStokesBoundaryData( MultiGridCL* &mgp, StokesBndDataCL* &bnddata,
        instat_vector_fun_ptr inflow, int GeomType, int bnd_type, std::vector<BndCondT>& BC);

/// \brief Create geometry of a Mzelle or a brick
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param bnddata a pointer of the boundary data object
 * @param inflow vector function pointer which is used to construct the boundary description, e.g. inflow velocity
 * @param meshfile_name name of the gambit-meshfile or list of integers/double which determine the L-/B-/Cavity-/Brick domain, e.g. 1x1x1@2x3x2
 * @param GeomType type of the geometry: 0: GAMBIT mesh, 1: brick, 2: cavity, 3: l-brick, 4: b-brick
 * @param bnd_type type of the boundary: 0: measuring cell, 1: hom. Dirichlet, 2-5: different Neumann conditions, a.t.m. only valid for a brick domain
 * @param deserialization_file file that contains the (serialized) grid, "none": use standard builder
 * @param BC description of the boundary
 * @param r_inlet inlet radius
 */
void CreateGeom (MultiGridCL* &mgp, StokesBndDataCL* &bnddata,
                 instat_vector_fun_ptr inflow,
                 const std::string& meshfile_name,
                 int GeomType, int bnd_type,
                 const std::string& deserialization_file, double& r_inlet);

/// \brief Create geometry for a poisson problem
/**
 * @param mgp a pointer of the MultiGridCL object
 * @param bnddata a pointer of the boundary data object
 * @param inflow scalar function pointer which is used to construct the boundary description, e.g. known solution
 * @param meshfile_name name of the gambit-meshfile or list of integers/double which determine the L-/B-/Cavity-/Brick domain, e.g. 1x1x1@2x3x2
 * @param GeomType type of the geometry: 0: GAMBIT mesh, 1: brick, 2: cavity, 3: l-brick, 4: b-brick
 * @param bnd_type type of the boundary: 0: measuring cell, 1: hom. Dirichlet, 2-5: different Neumann conditions, a.t.m. only valid for a brick domain
 * @param deserialization_file file that contains the (serialized) grid, "none": use standard builder
 * @param BC description of the boundary
 * @param r_inlet inlet radius
 */
void CreateGeomPoisson (MultiGridCL* &mgp, PoissonBndDataCL* &bnddata,
                 instat_scalar_fun_ptr inflow,
                 const std::string& meshfile_name,
                 int GeomType, int bnd_type,
                 const std::string& deserialization_file, double& r_inlet);

} // end of namespace drops

#endif /* GEOMSELECT_H_ */
