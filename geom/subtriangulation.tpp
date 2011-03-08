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

namespace DROPS
{

template <class VertexPartitionPolicyT,
          class VertexCutMergingPolicyT>
  void
  TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>::partition_principal_lattice (Uint num_intervals, const GridFunctionCL<>& ls)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( num_intervals);
    tetras_.resize( 0);
    pos_tetra_begin_= 0;
    vertexes_.resize( 0);

    VertexPartitionPolicyT vertex_policy( vertexes_, lat);
    VertexCutMergingPolicyT edgecut( lat.vertex_begin(), vertex_policy.cut_vertex_container());

    TetraContT loc_tetras; // temporary container for the positive tetras.
    double lset[4];
    Uint loc_vert_num;
    TetraT tet;
    for (PrincipalLatticeCL::const_tetra_iterator lattice_tet= lat.tetra_begin(), lattice_end= lat.tetra_end(); lattice_tet != lattice_end; ++lattice_tet) {
        for (Uint i= 0; i < 4; ++i)
            lset[i]= ls[(*lattice_tet)[i]]; 
        const RefTetraPartitionCL& cut= RefTetraPartitionCL::instance( lset);
        for (RefTetraPartitionCL::const_tetra_iterator it= cut.tetra_begin(), end= cut.tetra_end(); it != end; ++it) {
            for (Uint j= 0; j < 4; ++j) {
                loc_vert_num= (*it)[j];
                if (loc_vert_num < 4)
                    tet[j]= (*lattice_tet)[loc_vert_num];
                else { // Cut vertex
                    const Ubyte v0= VertOfEdge( loc_vert_num - 4, 0),
                                v1= VertOfEdge( loc_vert_num - 4, 1);
                    tet[j]= vertex_policy.cut_index_offset() + edgecut( (*lattice_tet)[v0], (*lattice_tet)[v1], lset[v0], lset[v1]);
                }
            }
            (cut.sign( it) == -1 ? tetras_ : loc_tetras).push_back( tet);
        }
    }
    pos_tetra_begin_= tetras_.size();
    std::copy( loc_tetras. begin(), loc_tetras.end(), std::back_inserter( tetras_));
    vertex_policy.sort_vertexes( ls, tetras_.begin(), tetras_.end(), pos_tetra_begin_);
}

template <class VertexPartitionPolicyT,
          class VertexCutMergingPolicyT>
  void
  write_paraview_vtu (std::ostream& file_, const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& t, TetraSignEnum s)
{
    file_ << "<?xml version=\"1.0\"?>"  << '\n'
          << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << '\n'
          << "<UnstructuredGrid>"   << '\n';

    file_<< "<Piece NumberOfPoints=\""<< t.vertexes_.size() <<"\" NumberOfCells=\""<< t.tetra_size( s) << "\">";
    file_<< "\n\t<Points>"
         << "\n\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << "ascii\">\n\t\t";
    for(LatticePartitionTypesNS::const_vertex_iterator it= t.vertexes_.begin(), end= t.vertexes_.end(); it != end; ++it) {
        file_ << it[0][1] << ' ' << it[0][2] << ' ' << it[0][3] << ' ';
    }
    file_<< "\n\t\t</DataArray> \n"
         << "\t</Points>\n";

    file_   << "\t<Cells>\n"
            << "\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\""
            <<"ascii\">\n\t\t";
    std::copy( t.tetra_begin( s), t.tetra_end( s), std::ostream_iterator<LatticePartitionTypesNS::TetraT>( file_));
    file_ << "\n\t\t</DataArray>\n";
    file_ << "\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t";
    for(Uint i= 1; i <= t.tetra_size( s); ++i) {
        file_ << i*4 << " ";
    }
    file_ << "\n\t\t</DataArray>";
    file_ << "\n\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t";
    const int tetraType= 10;
    for(Uint i= 1; i <= t.tetra_size( s); ++i) {
            file_ << tetraType<<" ";
    }
    file_<<"\n\t\t</DataArray>"
         <<"\n\t</Cells>";

    file_ <<"\n</Piece>"
          <<"\n</UnstructuredGrid>"
          <<"\n</VTKFile>";
}

template <class VertexPartitionPolicyT,
          class VertexCutMergingPolicyT>
  std::ostream&
  operator<< (std::ostream& out, const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& t)
{
    out << t.tetras_.size() << ' ' << t.pos_tetra_begin_ << ' ' << t.pos_vertex_begin_ << ' ' << t.neg_vertex_end_ << '\n';
    std::copy( t.tetras_.begin(), t.tetras_.end(), std::ostream_iterator<LatticePartitionTypesNS::TetraT>( out));
    out << '\n';
    std::copy( t.vertexes_.begin(), t.vertexes_.end(), std::ostream_iterator<BaryCoordCL>( out)) ;
    return out;
}

} // end of namespace DROPS