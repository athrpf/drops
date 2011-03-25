/// \file subtriangulation.cpp
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

#include "geom/subtriangulation.h"

namespace DROPS
{

void
write_paraview_vtu (std::ostream& file_, const TetraPartitionCL& t, TetraSignEnum s)
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

std::ostream&
operator<< (std::ostream& out, const TetraPartitionCL& t)
{
    out << t.tetras_.size() << ' ' << t.pos_tetra_begin_ << ' ' << t.pos_vertex_begin_ << ' ' << t.neg_vertex_end_ << '\n';
    std::copy( t.tetras_.begin(), t.tetras_.end(), std::ostream_iterator<LatticePartitionTypesNS::TetraT>( out));
    out << '\n';
    std::copy( t.vertexes_.begin(), t.vertexes_.end(), std::ostream_iterator<BaryCoordCL>( out)) ;
    return out;
}


void
SurfacePatchCL::partition_principal_lattice (Uint num_intervals, const std::valarray<double>& ls)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( num_intervals);
    PrincipalLatticeCL::const_vertex_iterator lattice_verts= lat.vertex_begin();

    triangles_.resize( 0);
    is_boundary_triangle_.resize( 0);

    std::valarray<byte> ls_sign;
    copy_levelset_sign( ls, ls_sign);

    MergeCutPolicyCL edgecut( lat.vertex_begin()); // Fixme
    std::vector<Uint> copied_vertexes( lat.num_vertexes(), static_cast<Uint>( -1));

    byte lset[4];
    Uint loc_vert_num;
    TriangleT tri;
    for (PrincipalLatticeCL::const_tetra_iterator lattice_tet= lat.tetra_begin(), lattice_end= lat.tetra_end(); lattice_tet != lattice_end; ++lattice_tet) {
        for (Uint i= 0; i < 4; ++i)
            lset[i]= ls_sign[(*lattice_tet)[i]]; 
        const RefTetraSurfacePatchCL& cut= RefTetraSurfacePatchCL::instance( lset);
        if (cut.empty()) continue;
        for (RefTetraSurfacePatchCL::const_triangle_iterator it= cut.triangle_begin(), end= cut.triangle_end(); it != end; ++it) {
            for (Uint j= 0; j < 3; ++j) {
                loc_vert_num= (*it)[j];
                if (loc_vert_num < 4) {
                    const Uint lattice_vert_num= (*lattice_tet)[loc_vert_num];
                    if (copied_vertexes[lattice_vert_num] == static_cast<Uint>( -1)) {
                        vertexes_.push_back( lattice_verts[lattice_vert_num]);
                        copied_vertexes[lattice_vert_num]= vertexes_.size() - 1;
                    }
                    tri[j]= copied_vertexes[lattice_vert_num];
                }
                else { // Cut vertex
                    const Ubyte v0= VertOfEdge( loc_vert_num - 4, 0),
                                v1= VertOfEdge( loc_vert_num - 4, 1);
                    tri[j]= edgecut( (*lattice_tet)[v0], (*lattice_tet)[v1], lset[v0], lset[v1]);
                }
            }
            triangles_.push_back( tri);
            is_boundary_triangle_.push_back( cut.is_boundary_triangle());
        }
    }
}

void
write_paraview_vtu (std::ostream& file_, const SurfacePatchCL& t)
{
    file_ << "<?xml version=\"1.0\"?>"  << '\n'
          << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << '\n'
          << "<UnstructuredGrid>"   << '\n';

    file_<< "<Piece NumberOfPoints=\""<< t.vertexes_.size() <<"\" NumberOfCells=\""<< t.triangles_.size() << "\">";
    file_<< "\n\t<Points>"
         << "\n\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << "ascii\">\n\t\t";
    for(SurfacePatchCL::const_vertex_iterator it= t.vertexes_.begin(), end= t.vertexes_.end(); it != end; ++it) {
        file_ << it[0][1] << ' ' << it[0][2] << ' ' << it[0][3] << ' ';
    }
    file_<< "\n\t\t</DataArray> \n"
         << "\t</Points>\n";

    file_   << "\t<Cells>\n"
            << "\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\""
            <<"ascii\">\n\t\t";
    std::copy( t.triangles_.begin(), t.triangles_.end(), std::ostream_iterator<LatticePartitionTypesNS::TriangleT>( file_));
    file_ << "\n\t\t</DataArray>\n";
    file_ << "\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t";
    for(Uint i= 1; i <= t.triangles_.size(); ++i) {
        file_ << i*3 << " ";
    }
    file_ << "\n\t\t</DataArray>";
    file_ << "\n\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t";
    const int Type= 5; // Triangles
    for(Uint i= 1; i <= t.triangles_.size(); ++i) {
            file_ << Type<<" ";
    }
    file_<<"\n\t\t</DataArray>"
         <<"\n\t</Cells>";

    file_ <<"\n</Piece>"
          <<"\n</UnstructuredGrid>"
          <<"\n</VTKFile>";
}


void
UnorderedVertexPolicyCL::sort_vertexes (VertexContT& vertexes, VertexContT& cut_vertexes,
    size_t& pos_vertex_begin, size_t& neg_vertex_end)
{
    std::copy( cut_vertexes.begin(), cut_vertexes.end(), std::back_inserter( vertexes));
    pos_vertex_begin= neg_vertex_end= 0;
}

void
SortedVertexPolicyCL::sort_vertexes (VertexContT& vertexes, VertexContT& cut_vertexes,
    size_t& pos_vertex_begin, size_t& neg_vertex_end)
{
    const Uint lattice_num_vertexes= lat_.num_vertexes();
    const PrincipalLatticeCL::const_vertex_iterator lattice_vertex_begin= lat_.vertex_begin();

    std::valarray<byte> ls_sign;
    copy_levelset_sign( ls_, ls_sign);

    // Count signs
    Uint num_sign_arr[3]= { 0, 0, 0 };
    Uint* const num_sign= num_sign_arr + 1; // num_sign[i] == number of verts with sign i
    for (Uint i= 0; i < lattice_num_vertexes; ++i)
        ++num_sign[ls_sign[i]];
    const Uint num_zero_vertexes= num_sign[0] + cut_vertexes.size();

    vertexes.resize( num_sign[-1] + num_sign[1] + num_zero_vertexes);
    pos_vertex_begin= num_sign[-1];
    neg_vertex_end=   num_sign[-1] + num_zero_vertexes;

    std::vector<Uint> new_pos( num_sign[-1] + num_sign[1] + num_zero_vertexes); // maps old vertex-index to the new one
    size_t cursor_arr[3];
    size_t* const cursor= cursor_arr + 1; // Insertion cursors for the sorted-by-sign vertex numbers
    cursor[-1]= 0;
    cursor[0]= num_sign[-1];
    cursor[1]= num_sign[-1] + num_zero_vertexes;
    for (Uint i= 0; i < lattice_num_vertexes; ++i) {
        size_t& cur= cursor[ls_sign[i]];
        new_pos[i]= cur;
        vertexes[cur]= lattice_vertex_begin[i];
        ++cur;
    }
    size_t& cur= cursor[0];
    for (Uint i= 0; i < cut_vertexes.size(); ++i, ++cur) {
        new_pos[i + lattice_num_vertexes]= cur;
        vertexes[cur]= cut_vertexes[i];
    }
    // Reorder the indices in the tetras
    for (TetraContT::iterator it= tetra_begin_; it != tetra_end_; ++it)
        for (Uint i= 0; i < 4; ++i)
            (*it)[i]= new_pos[(*it)[i]];
}

void
PartitionedVertexPolicyCL::sort_vertexes (VertexContT& vertexes, VertexContT& cut_vertexes, size_t& pos_vertex_begin, size_t& neg_vertex_end)
{
    pol_.sort_vertexes( vertexes, cut_vertexes, pos_vertex_begin, neg_vertex_end);
    const Uint num_zero_vertexes= neg_vertex_end - pos_vertex_begin;

    vertexes.reserve( vertexes.size() + num_zero_vertexes); // Make sure that inserting the new zero vertexes will not invalidate the iterators indicating the sequence of zero-vertexes itself in the same container.
    vertexes.insert( vertexes.begin() + neg_vertex_end, vertexes.begin() + pos_vertex_begin, vertexes.begin() + neg_vertex_end);
    pos_vertex_begin= neg_vertex_end; // This is the partitioning of the sequence of vertexes.

    // Adjust the indices of all vertexes in the positive tetras
    for (TetraContT::iterator it= tetra_begin_ + pos_tetra_begin_; it < tetra_end_; ++it)
        for (Uint i= 0; i < 4; ++i)
            (*it)[i]+= num_zero_vertexes;
}

} // end of namespace DROPS