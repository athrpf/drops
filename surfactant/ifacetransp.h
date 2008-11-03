/// \file
/// \brief Discretization for PDEs on an interface.

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/discretize.h"
#include "levelset/mgobserve.h"

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
    const VecDescCL& ls, double omit_bound= -1./*default to using all dof*/);


/// \brief Routine to number unknowns on the vertices surrounding an
/// interface.
///
/// \todo This belongs, together with XFe, in a file in misc for interface dependent FE.
///
/// Uses CreateNumbOnInterfaceVertex on the triangulation with level level on the multigrid mg.
/// One can only create P1-elements.
void CreateNumbOnInterface(Uint level, IdxDescCL& idx, MultiGridCL& mg,
    const VecDescCL& ls, double omit_bound= -1./*default to using all dof*/);

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of x
///        are copied to xext and the remaining values of xext are set to 0.
void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext);

/// \brief Given a P1-vecdesc on the interface x and a P1-vecdesc xext on mg, the values of xext
///        are copied to x, where x has an unknown.
void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x);

/// \brief The routine sets up the mass-matrix in matM on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls);


/// \brief The routine sets up the Laplace-Beltrami-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// D is the diffusion-coefficient
void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, double D);


/// \brief The routine sets up the convection-matrix in mat on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
///
/// The template-parameter is only used to circumvent the exact type of the discrete
/// velocity solution in the Stokes classes.
template <class DiscVelSolT>
void SetupConvectionP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const DiscVelSolT& v);

/// \brief Helper of SetupConvectionP1
void SetupConvectionP1OnTriangle (const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> qp1[4],
    const LocalP2CL<Point3DCL>& u, Point3DCL grad[4], double coup[4][4]);

/// \brief The routine sets up the mass-matrix scaled with div_\Gamma v in mat on the interface
///        defined by ls. It belongs to the FE induced by standard P1-elements.
///
/// The template-parameter is only used to circumvent the exact type of the discrete
/// velocity solution in the Stokes classes.
template <class DiscVelSolT>
void SetupMassDivP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const DiscVelSolT& v);

/// \brief Helper of SetupMassDivP1
void SetupMassDivP1OnTriangle (const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> qp1[4],
    const LocalP2CL<Point3DCL>& u, LocalP1CL<Point3DCL> gradp2[10], const Point3DCL& n, double coup[4][4]);

/// \brief The routine sets up the load-vector in v on the interface defined by ls.
///        It belongs to the FE induced by standard P1-elements.
void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, instat_scalar_fun_ptr f);

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair an interface p1-function.
///
/// The function is extended to a P1-function on \Omega and then repaired as such. In
/// post_refine_sequence() the new interface-index is generated and the function is restricted
/// to it.
class InterfaceP1RepairCL : public MGObserverCL
{
  private:
    MultiGridCL& mg_;
    const VecDescCL&   lset_vd_;
    VecDescCL&         u_;

    IdxDescCL          fullp1idx_;
    VecDescCL          fullu_;

  public:
    InterfaceP1RepairCL (MultiGridCL& mg, const VecDescCL& lset_vd, VecDescCL& u)
        : mg_( mg), lset_vd_( lset_vd), u_( u), fullp1idx_( P1_FE), fullu_( &fullp1idx_) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  ();
    void post_refine_sequence ();
};

} // end of namespace DROPS

#include "surfactant/ifacetransp.tpp"

#endif
