/// \file subtriangulation.tpp
/// \brief Triangulation of a principal-lattice of a tetra adapted to a piecewise linear level-set function.
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

#include "geom/topo.h"

#include <iterator>

namespace DROPS
{

///\brief Write the sign of the levelset function src to the sequence dst.
/// \return end-iterator of the sequence of written signs
inline void
copy_levelset_sign ( const std::valarray<double>& src, std::valarray<byte>& dst)
{
    dst.resize( src.size());
    for (size_t i= 0; i < src.size(); ++i)
        dst[i]= sign( src[i]);
}

template <class VertexPartitionPolicyT,
          class VertexCutMergingPolicyT>
  void
  TetraPartitionCL::partition_principal_lattice (Uint num_intervals, const std::valarray<double>& ls)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( num_intervals);
    const Uint lattice_num_vertexes= lat.num_vertexes();

    tetras_.resize( 0);
    pos_tetra_begin_= 0;
    vertexes_.resize( 0);

    std::valarray<byte> ls_sign;
    copy_levelset_sign( ls, ls_sign);

    VertexCutMergingPolicyT edgecut( lat.vertex_begin()); // Stores the genuine cuts.

    TetraContT loc_tetras; // temporary container for the positive tetras.
    byte lset[4];
    Uint loc_vert_num;
    TetraT tet;
    for (PrincipalLatticeCL::const_tetra_iterator lattice_tet= lat.tetra_begin(), lattice_end= lat.tetra_end(); lattice_tet != lattice_end; ++lattice_tet) {
        for (Uint i= 0; i < 4; ++i)
            lset[i]= ls_sign[(*lattice_tet)[i]]; 
        const RefTetraPartitionCL& cut= RefTetraPartitionCL::instance( lset);
        for (RefTetraPartitionCL::const_tetra_iterator it= cut.tetra_begin(), end= cut.tetra_end(); it != end; ++it) {
            for (Uint j= 0; j < 4; ++j) {
                loc_vert_num= (*it)[j];
                if (loc_vert_num < 4)
                    tet[j]= (*lattice_tet)[loc_vert_num];
                else { // Cut vertex
                    const Ubyte v0= VertOfEdge( loc_vert_num - 4, 0),
                                v1= VertOfEdge( loc_vert_num - 4, 1);
                    tet[j]= lattice_num_vertexes + edgecut( (*lattice_tet)[v0], (*lattice_tet)[v1], lset[v0], lset[v1]);
                }
            }
            (cut.sign( it) == -1 ? tetras_ : loc_tetras).push_back( tet);
        }
    }
    pos_tetra_begin_= tetras_.size();
    std::copy( loc_tetras. begin(), loc_tetras.end(), std::back_inserter( tetras_));

    VertexPartitionPolicyT vertex_order_policy( lat, ls, tetras_.begin(), tetras_.end(), pos_tetra_begin_);
    vertex_order_policy.sort_vertexes( vertexes_, edgecut.cut_vertex_container(), pos_vertex_begin_, neg_vertex_end_);
}

} // end of namespace DROPS