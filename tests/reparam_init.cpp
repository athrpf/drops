/// \file reparam_init.cpp
/// \brief tests initialization phase of reparametrization
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

#include "misc/utils.h"
#include "geom/multigrid.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "num/discretize.h"
#include "num/fe.h"
#include "num/spmat.h"
#include "misc/problem.h"
#include "out/ensightOut.h"
#include "levelset/adaptriang.h"
#include <tr1/unordered_map>

using namespace DROPS;

Point3DCL orig( 0.);
const double r= 0.6;

double sphere2 (const Point3DCL& p, double)
{
    return (p - orig).norm_sq() - r*r;
}

double sphere2_stat (const Point3DCL& p)
{
    return sphere2( p, 0.);
}

double sphere_dist (const Point3DCL& p, double)
{
    return (p - orig).norm() - r;
}

double sphere_dist_stat (const Point3DCL& p)
{
    return sphere_dist( p, 0.);
}

Point3DCL sphere_dist_grad (const Point3DCL& p, double)
{
    return p.norm() == 0. ? Point3DCL() : (p - orig)/(p - orig).norm();
}


// Computes the maximum of |n - ng|. n is the normal on the patch, ng is the gradient of the exact distance function in in barycenter of the patch.
double facet_sup_norm(const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls)
{
    const DROPS::Uint lvl= ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;

    double maxn= 0.;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls);
	if (!triangle.Intersects()) continue; // We are not at the phase boundary.

	for (int ch= 0; ch < 8; ++ch) {
	    if (!triangle.ComputeForChild( ch)) continue; // Child ch has no intersection

	    for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                const Point3DCL n= triangle.GetNormal();
                const Point3DCL bary= 1./3.*(triangle.GetPoint( tri) + triangle.GetPoint( tri + 1) + triangle.GetPoint( tri + 2));
		const Point3DCL ng= sphere_dist_grad( bary, 0.);
                maxn= std::max( maxn, (n -ng).norm());
	    }
	}
    }
    return maxn;
}

double dist_to_sphere(const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls)
{
    const DROPS::Uint lvl= ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;

    double dd= 0.;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls);
	if (triangle.Intersects()) { // We are at the phase boundary.
	    for (int ch= 0; ch < 8; ++ch) {
	    	triangle.ComputeForChild( ch);
		for (Uint tript= 0; tript < triangle.GetNumPoints(); ++tript) {
		    dd= std::max( dd, std::fabs( sphere_dist( triangle.GetPoint( tript), 0.)));
		}
            }
	}
    }
    return dd;
}

template<class DiscP2FunType>
double vertex_sup_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls,
    const DiscP2FunType& f)
{
    const DROPS::Uint lvl= ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;

    double dd=  0.;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls);
	if (triangle.Intersects()) { // We are at the phase boundary.
	    for (int ch= 0; ch < 8; ++ch) {
	    	triangle.ComputeForChild( ch);
		for (Uint tript= 0; tript < triangle.GetNumPoints(); ++tript) {
		        dd= std::max( dd, std::fabs( f.val( *it, triangle.GetBary( tript))));
		}
            }
	}
    }
    return dd;
}

void CheckSigns (MultiGridCL& mg, const VecDescCL& Phi, const VecDescCL* Phi2= 0)
{
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    for (MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
        end= mg.GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        if ( Phi.Data[it->Unknowns(idx)]*sphere_dist_stat( it->GetCoord()) <= 0.) {
            std::cout << "Different signs on vertex." << std::endl;
         }
        if (Phi2) {
            if (Phi.Data[it->Unknowns(idx)]*Phi2->Data[it->Unknowns(idx)] <= 0.) {
                std::cout << "reparametrization and original: Different signs on vertex." << std::endl;
            }
        }
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(lvl),
        end= mg.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        if ( Phi.Data[it->Unknowns(idx)]*sphere_dist_stat( GetBaryCenter( *it)) <= 0.) {
            std::cout << "Different signs on edge." << std::endl;
        }
        if (Phi2) {
            if (Phi.Data[it->Unknowns(idx)]*Phi2->Data[it->Unknowns(idx)] <= 0.) {
                std::cout << "reparametrization and original: Different signs on edge." << std::endl;
            }
        }
    }
}


int main ()
{
  try {
    orig[0]= 0.005; orig[1]= 0.003; orig[2]= 0.001;
    // TestDist();
    // return 0;

    int numref;
    std::cout << " numref: ";
    std::cin >> numref;
    int f_lvl;
    std::cout << " f_lvl: ";
    std::cin >> f_lvl;

    DROPS::BrickBuilderCL brick( Point3DCL( -1.0),
                                  2.0*std_basis<3>(1),
                                  2.0*std_basis<3>(2),
                                  2.0*std_basis<3>(3),
                                  numref, numref, numref);
    MultiGridCL mg( brick);
    AdapTriangCL adap( mg, 0.1, 0, f_lvl);
    adap.MakeInitialTriang( sphere_dist_stat);

    // lset contains the quadratic levelset-function of the sphere
    // vd_dist ist the piecewise quadratic interpolation of the signed distance-function of the sphere.
    SurfaceTensionCL sf( /*surface tension*/ &sphere2);
    LevelsetP2CL lset( mg, sf);
    lset.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset.Phi.SetIdx( &lset.idx);

    lset.Init( &sphere2_stat);
    BndDataCL<> nobnd( 0);

    LevelsetP2CL lset_d( mg, sf);
    lset_d.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset_d.Phi.SetIdx( &lset_d.idx);
    lset_d.Init( &sphere_dist_stat);

    VecDescCL vd_dist( &lset.idx);
    vd_dist.Data= lset_d.Phi.Data;
    std::cout << "sup of (the P2-interpolant of) dist_\\Gamma  on \\Gamma_h: "
        << vertex_sup_norm( mg, lset.Phi, lset_d.GetSolution()) << std::endl;
    std::cout << "sup of dist_\\Gamma on \\Gamma_h: " << dist_to_sphere( mg, lset.Phi) << std::endl;
    std::cout << "sup of gradient-difference on \\Gamma_h: "
        << facet_sup_norm( mg, lset.Phi) << std::endl;
    // CheckSigns ( mg, lset.Phi);


    // Initialize Ensight6 output
    std::string ensf( "./reparam_init");
    Ensight6OutCL ensight( "reparam_init.case", 0);
    ensight.Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(),   "cube    ",     ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),      "P2_Levelset",  ensf + "_p2.scl"));
    ensight.Register( make_Ensight6Scalar( lset_d.GetSolution(),  "Dist",     ensf + "_dist.scl"));

    LevelsetP2CL lset_rep( mg, sf);
    lset_rep.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset_rep.Phi.SetIdx( &lset_rep.idx);
    lset_rep.Init( &sphere2_stat);
    lset_rep.ReparamFastMarching( true, false, false, false, 2);
    ensight.Register( make_Ensight6Scalar( lset_rep.GetSolution(),  "Reparam",     ensf + "_reparam.scl"));
    ensight.Write();

    std::cout << "after reparametrization: sup of (the P2-interpolant of) dist_\\Gamma  on \\Gamma_h: "
        << vertex_sup_norm( mg, lset_rep.Phi, make_P2Eval( mg, nobnd, vd_dist)) << std::endl;
    std::cout << "after reparametrization: sup of dist_\\Gamma on \\Gamma_h: "
        << dist_to_sphere( mg, lset_rep.Phi) << std::endl;

    std::cout << "sup of gradient-difference on \\Gamma_h: "
        << facet_sup_norm( mg, lset_rep.Phi) << std::endl;
    // CheckSigns ( mg, lset_rep.Phi, &lset.Phi);
  }
  catch( DROPSErrCL d) {
      d.handle();
  }
    return 0;
}
