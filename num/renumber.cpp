/// \file renumber.cpp
/// \brief implementation of bandwith reduction and downwind numbering
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

#include "num/renumber.h"
#include "num/discretize.h"
#include "misc/problem.h"

namespace DROPS
{

void DownwindAccu_P2CL::begin_accumulation ()
{
    std::cout << "entering DownwindAccu_P2CL::begin_accumulation";
    const size_t num_unks= idx.NumUnknowns()/idx.NumUnknownsVertex();
    mC= new SparseMatBuilderCL<double>( &C, num_unks, num_unks);
}

void DownwindAccu_P2CL::finalize_accumulation ()
{
    mC->Build();
    delete mC;
#ifndef _PAR
    std::cout << ": " << C.num_nonzeros() << " nonzeros in C. ";
#endif
    std::cout << '\n';
}

void DownwindAccu_P2CL::visit (const TetraCL& tet)
{
    SMatrixCL<10,10>     loc; ///< local matrix
    LocalNumbP2CL        n;   ///< global numbering of the P2-unknowns
    LocalP2CL<Point3DCL> vel_loc( tet, vel, velBndData); ///< local velocity

    n.assign_indices_only( tet, idx);
    local_setup( tet, vel_loc, loc);
    update_global_system( n, loc);
}

void DownwindAccu_P2CL::local_setup (const TetraCL& tet, const LocalP2CL<Point3DCL>& vel_loc,
    SMatrixCL<10,10>& loc)
{
    for (Uint i= 0; i < 10; ++i) {
        const BaryCoordCL& bi= i < 4 ? std_basis<4>( i + 1)
                                     : 0.5*(std_basis<4>( VertOfEdge(i - 4, 0) + 1) + std_basis<4>( VertOfEdge(i - 4, 1) + 1));
        const Point3DCL& pi=   i < 4 ? tet.GetVertex( i)->GetCoord() : GetBaryCenter( *tet.GetEdge( i - 4));
        for (Uint j= 0; j < 10; ++j) {
            const BaryCoordCL& bj= j < 4 ? std_basis<4>( j + 1)
                                         : 0.5*(std_basis<4>( VertOfEdge(j - 4, 0) + 1) + std_basis<4>( VertOfEdge(j - 4, 1) + 1));
            const Point3DCL& pj=   j < 4 ? tet.GetVertex( j)->GetCoord() : GetBaryCenter( *tet.GetEdge( j - 4));
            loc( i, j)= inner_prod( pi - pj,  vel_loc( 0.5*(bj + bi)));
         }
    }
}

void DownwindAccu_P2CL::update_global_system (const LocalNumbP2CL& n, const SMatrixCL<10,10>& loc)
{
    SparseMatBuilderCL<double>& C= *mC;
    const Uint num_components= idx.NumUnknownsVertex();

    for(int i= 0; i < 10; ++i)  // assemble row Numb[i]
        if (n.WithUnknowns( i)) // dof i is not on a Dirichlet boundary
            for(int j= 0; j < 10; ++j)
                if (n.WithUnknowns( j)) // dof j is not on a Dirichlet boundary
                    C( n.num[i]/num_components, n.num[j]/num_components)= loc( i, j);
}


void MassAccu_P2CL::begin_accumulation ()
{
    std::cout << "entering MassAccu_P2CL::begin_accumulation";
    const size_t num_unks= idx.NumUnknowns()/idx.NumUnknownsVertex();
    Mb= new SparseMatBuilderCL<double>( &M, num_unks, num_unks);
}

void MassAccu_P2CL::finalize_accumulation ()
{
    Mb->Build();
    delete Mb;
#ifndef _PAR
    std::cout << ": " << M.num_nonzeros() << " nonzeros in M. ";
#endif
    std::cout << '\n';
}

void MassAccu_P2CL::visit (const TetraCL& tet)
{
    SMatrixCL<10,10>     loc;     ///< local matrix
    LocalNumbP2CL        n;       ///< global numbering of the P2-unknowns

    n.assign_indices_only( tet, idx);
    local_setup( tet, loc);
    update_global_system( n, loc);
}

void MassAccu_P2CL::local_setup (const TetraCL& tet, SMatrixCL<10,10>& loc)
{
    const double absdet= tet.GetVolume()*6.;
    for (Uint i= 0; i < 10; ++i) {
        for (Uint j= 0; j < 10; ++j) {
            loc( i, j)= absdet*P2DiscCL::GetMass( i, j);
         }
    }
}

void MassAccu_P2CL::update_global_system (const LocalNumbP2CL& n, const SMatrixCL<10,10>& loc)
{
    SparseMatBuilderCL<double>& M= *Mb;
    const Uint num_components= idx.NumUnknownsVertex();

    for(int i= 0; i < 10; ++i)  // assemble row Numb[i]
        if (n.WithUnknowns( i)) // dof i is not on a Dirichlet boundary
            for(int j= 0; j < 10; ++j)
                if (n.WithUnknowns( j)) // dof j is not on a Dirichlet boundary
                    M( n.num[i]/num_components, n.num[j]/num_components)+= loc( i, j);
}

} // end of namspace DROPS
