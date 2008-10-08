/// \file
/// \brief Discretization for PDEs on an interface.

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/discretize.h"

#ifndef DROPS_IFACETRANSP_H
#define DROPS_IFACETRANSP_H

namespace DROPS {

/// \brief Routine to number unknowns on the vertices surrounding an
/// interface.
///
/// \todo This belongs, together with XFe, in a file in misc for interface dependent FE.
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are cut by the zero level of lset. If InterfacePatchCL::IntersectsInterior() ist true,
/// all vertices are numbered, if only InterfacePatchCL::Intersects() is true, only the
/// vertices with InterfacePatchCL::GetSign(vertex) == 0 are numbered. The other vertices in
/// such tetrahedra obtain NoIdx as number, but they are not counted as unknowns.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// A more user friendly interface is provided by CreateNumbOnInterface.
void CreateNumbOnInterfaceVertex (const Uint idx, IdxT& counter, Uint stride,
    const MultiGridCL::TriangTetraIteratorCL& begin,
    const MultiGridCL::TriangTetraIteratorCL& end,
    const VecDescCL& ls);

/// \brief Routine to number unknowns on the vertices surrounding an
/// interface.
///
/// \todo This belongs, together with XFe, in a file in misc for interface dependent FE.
///
/// Uses CreateNumbOnInterfaceVertex on the triangulation with level level on the multigrid mg.
/// One can only create P1-elements.
void CreateNumbOnInterface(Uint level, IdxDescCL& idx, MultiGridCL& mg,
    const VecDescCL& ls);

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of x
///        are copied to xext and the remaining values of xext are set to 0.
void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext);

void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls);

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls);

void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, instat_scalar_fun_ptr f);

} // end of namespace DROPS

#endif
