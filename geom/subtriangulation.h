/// \file subtriangulation.h
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

#ifndef DROPS_SUBTRIANGULATION_H
#define DROPS_SUBTRIANGULATION_H

#include "misc/container.h"
#include "geom/principallattice.h"

#include <vector>
#include <unordered_map>


namespace DROPS {

/// \brief Common types used by TetraPartitionCL, SurfacePatchCL and their helpers
namespace LatticePartitionTypesNS {

typedef SArrayCL<Uint, 4>          TetraT; ///< representation of a tetra of the partition via its vertices: index in the vertex_-array
typedef std::vector<TetraT>        TetraContT;
typedef TetraContT::const_iterator const_tetra_iterator;

typedef SArrayCL<Uint, 3>             TriangleT; ///< representation of a triangle of the interface via its vertices: index in the vertex_-array
typedef std::vector<TriangleT>        TriangleContT;
typedef TriangleContT::const_iterator const_triangle_iterator;

typedef std::vector<BaryCoordCL>    VertexContT;
typedef VertexContT::const_iterator const_vertex_iterator;

} // end of namespace LatticePartitionTypesNS


/// \brief Vertices are ordered with respect to the sign of the levelset function as follows: First the vertexes from the principal lattice, then all proper cut-vertexes.
class UnorderedVertexPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::TetraContT  TetraContT;
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    VertexContT& partition_vertexes_;

  public:
    UnorderedVertexPolicyCL (VertexContT& vertexes, const PrincipalLatticeCL& lat) : partition_vertexes_( vertexes) {
        partition_vertexes_.resize( 0);
        partition_vertexes_.reserve( lat.num_vertexes());
        std::copy( lat.vertex_begin(), lat.vertex_end(), std::back_inserter( partition_vertexes_));
    }

    ///\brief The container, in which TetraPartitionCL will store proper cuts of an edge: Simply the container for the vertexes of the principal lattice
    VertexContT& cut_vertex_container () { return partition_vertexes_; }
    ///\brief Offset of the proper cuts in relation to the vertexes of the principal lattice
    /// Returns zero, as the vertexes are stored in the same container as the partition vertexes.
    Uint cut_index_offset () { return 0; }

    /// \brief Sort the vertexes and update the vertex numbers in the tetras: a no-op for this policy
    void sort_vertexes (const GridFunctionCL<>&, TetraContT::iterator, TetraContT::iterator, size_t) {}
};

/// \brief Vertices are ordered with respect to the sign of the levelset function as follows: First the negative vertexes, then the zero vertexes, then the positive vertexes. The vertexes of the negative and of the positive tetras potentially overlap in the zero vertexes.
class SortedVertexPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::TetraContT  TetraContT;
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    VertexContT& vertexes_;
    VertexContT  cut_vertexes_;
    const PrincipalLatticeCL::const_vertex_iterator lattice_vertex_begin_;
    const Uint lattice_num_vertexes_;

  public:
    SortedVertexPolicyCL (VertexContT& vertexes, const PrincipalLatticeCL& lat)
        : vertexes_( vertexes),  lattice_vertex_begin_( lat.vertex_begin()),
          lattice_num_vertexes_( lat.num_vertexes()) {}

    ///\brief The container, in which TetraPartitionCL will store proper cuts of an edge
    VertexContT& cut_vertex_container () { return cut_vertexes_; }
    ///\brief Offset of the proper cuts in relation to the vertexes of the principal lattice
    Uint cut_index_offset () { return lattice_num_vertexes_; }

    /// \brief Sort the vertexes and update the vertex numbers in the tetras.
    void sort_vertexes (const GridFunctionCL<>& ls, TetraContT::iterator tetra_begin, TetraContT::iterator tetra_end, size_t pos_tetra_begin);
};

/// \brief Vertices are ordered with respect to the sign of the levelset function as follows: First the negative vertexes, then the zero vertexes, then the positive vertexes. Opposite to SortedVertexPolicyC, all zero vertexes are duplicated, such that the vertexes of the negative tetras and of the positive tetras are disjoint.
class PartitionedVertexPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::TetraContT  TetraContT;
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    VertexContT& partition_vertexes_;
    VertexContT  cut_vertexes_;
    const PrincipalLatticeCL::const_vertex_iterator lattice_vertex_begin_;
    const Uint lattice_num_vertexes_;

  public:
    PartitionedVertexPolicyCL (VertexContT& vertexes, const PrincipalLatticeCL& lat)
        : partition_vertexes_( vertexes),  lattice_vertex_begin_( lat.vertex_begin()),
          lattice_num_vertexes_( lat.num_vertexes()) {}

    ///\brief The container, in which TetraPartitionCL will store proper cuts of an edge
    VertexContT& cut_vertex_container () { return cut_vertexes_; }
    ///\brief Offset of the proper cuts in relation to the vertexes of the principal lattice
    Uint cut_index_offset () { return lattice_num_vertexes_; }

    /// \brief Sort the vertexes and update the vertex numbers in the tetras: Special care must be taken for the duplicated vertexes
    void sort_vertexes (const GridFunctionCL<>& ls, TetraContT::iterator tetra_begin,
        TetraContT::iterator tetra_end, size_t pos_tetra_begin);
};

///\brief A cut-vertex is added to the list of vertexes for each tetra, on which it is discovered; fast, but leads to more vertices, which in turn leads to more dof for quadrature rules that use the vertexes of the partition.
class DuplicateCutPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    const PrincipalLatticeCL::const_vertex_iterator partition_vertexes_;
    VertexContT& vertexes_;

  public:
    DuplicateCutPolicyCL (PrincipalLatticeCL::const_vertex_iterator partition_vertexes, VertexContT& vertexes) : partition_vertexes_( partition_vertexes), vertexes_( vertexes) {}

    ///\brief Add the cut vertex and return its number.
    Uint operator() (Uint v0, Uint v1, double ls0, double ls1) {
        const double edge_bary1_cut= ls0/(ls0 - ls1); // the root of the level set function on the edge
        vertexes_.push_back( ConvexComb( edge_bary1_cut, partition_vertexes_[v0], partition_vertexes_[v1]));
        return vertexes_.size() - 1;
    }
};

///\brief A cut-vertex is added to the list of vertexes only by the first tetra, on which it is discovered: cuts are memoized for each edge.
class MergeCutPolicyCL
{
  private:
    typedef LatticePartitionTypesNS::VertexContT VertexContT;

    struct UintPairHasherCL
    {
        size_t operator() (const std::pair<Uint, Uint>& p) const
            { return p.first << (4*sizeof(Uint)) ^ p.second; } // for less than 2^(sizeof(Uint)/2) vertices this is a bijection
    };

    typedef std::pair<Uint, Uint> EdgeT;
    typedef std::tr1::unordered_map<EdgeT, Uint, UintPairHasherCL> EdgeToCutMapT;

    const PrincipalLatticeCL::const_vertex_iterator partition_vertexes_;
    VertexContT& vertexes_;
    EdgeToCutMapT edge_to_cut_;

  public:
    MergeCutPolicyCL (PrincipalLatticeCL::const_vertex_iterator partition_vertexes, VertexContT& vertexes)
        : partition_vertexes_( partition_vertexes), vertexes_( vertexes) {}

    ///\brief Return the number of the cut vertex, if it is already memoized, otherwise add it and return its number.
    Uint operator() (Uint v0, Uint v1, double ls0, double ls1) {
        const EdgeT e= v0 < v1 ? std::make_pair( v0, v1) : std::make_pair( v1, v0);
        EdgeToCutMapT::const_iterator e_it= edge_to_cut_.find( e);
        if (e_it != edge_to_cut_.end())
            return e_it->second;
        else {
            const double edge_bary1_cut= ls0/(ls0 - ls1); // the root of the level set function on the edge
            vertexes_.push_back( ConvexComb( edge_bary1_cut, partition_vertexes_[v0], partition_vertexes_[v1]));
            return edge_to_cut_[e]= vertexes_.size() - 1;
        }
    }
};


///\brief forward declaration of TetraPartitionCL
template <class= SortedVertexPolicyCL, class= MergeCutPolicyCL> class TetraPartitionCL;

///\brief declaration of debug output (neccessary due to friend declaration in TetraPartitionCL)
template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  std::ostream&
  operator<< (std::ostream&, const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>&);

///\brief declaration of .vtu output (neccessary due to friend declaration in TetraPartitionCL)
template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  void
  write_paraview_vtu (std::ostream&, const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>&, TetraSignEnum= AllTetraC);


///\brief Partition the principal lattice of a tetra t (n intervals on each edge) according to the sign of a levelset function ls.
///
/// The sub-tetras, which are cut by the interface are triangulated to match the interface. The values of the level set function on the vertices of the principal lattice must be prescribed. The sequence of all tetrahedra contains the negative tetrahedra as initial subsequence.
///
/// For the vertexes, there are two properties, which can be selected by the policy template-parameters:
/// VertexPartitionPolicyT: Determines, how the vertices are ordered with respect to the sign of the levelset function. This is important for the generation of composite quadrature rules.
///     Unordered:   First the vertexes from the principal lattice, then all proper cut-vertexes.
///     Sorted:      First the negative vertexes, then the zero vertexes, then the positive vertexes. The vertexes of the negative and of the positive tetras potentially overlap in the zero vertexes.
///     Partitioned: Like sorted, and all zero vertexes are duplicated, such that the vertexes of the negative tetras and of the positive tetras are disjoint.
/// Unordered is fastest and appropriate for quadrature rules that use no vertexes of the triangulation. Quadrature rules that use vertexes need Sorted, if a continuous integrand is prescribed (and integrated on the negative or positive or both domains), and Partitioned, if the integrand is discontinuous.
///
/// VertexCutMergingPolicyT: Determines, how often cut vertexes are stored.
///     Duplicate: A cut-vertex is added for each tetra, on which it is discovered; fast, but leads to more vertices, which in turn leads to more dof for quadrature rules that use the vertexes of the partition.
///     Merge: The edge cuts are memoized for each edge and added only once -- default as it leads to the minimal amount of cut vertexes.
template <class VertexPartitionPolicyT,
          class VertexCutMergingPolicyT>
class TetraPartitionCL
{
  public:
    typedef LatticePartitionTypesNS::TetraT               TetraT;
    typedef LatticePartitionTypesNS::TetraContT           TetraContT;
    typedef LatticePartitionTypesNS::const_tetra_iterator const_tetra_iterator;

    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

  private:
    TetraContT tetras_;       ///< All tetras of the partition.
    size_t pos_tetra_begin_;  ///< begin of the subsequence of positive tetras

    VertexContT vertexes_;     ///< All vertices of the partition. 
    size_t pos_vertex_begin_;  ///< begin of the subsequence of vertexes of positive tetras
    size_t neg_vertex_end_;    ///< end of the subsequence of of vertexes of negative tetras

  public:
    TetraPartitionCL () : pos_tetra_begin_( 0), pos_vertex_begin_( 0), neg_vertex_end_( 0) {} ///< Empty default-cut

    ///\brief Computes the partition of the principal lattice with num_intervals on each edge of the reference-tetra given the level set values in ls.
    void partition_principal_lattice (Uint num_intervals, const GridFunctionCL<>& ls);

    size_t tetra_size  (TetraSignEnum s= AllTetraC) const ///< number of tetras with given sign
         { return tetra_end( s) - tetra_begin( s); }
    size_t vertex_size (TetraSignEnum s= AllTetraC) const ///< number of vertexes used by the tetras of the corresponding sign; depending on VertexCutMergingPolicyT, interface vertexes occur multiple times.
         { return vertex_end( s) - vertex_begin( s); }

    /// Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum.
    ///@{
    const_tetra_iterator tetra_begin (TetraSignEnum s= AllTetraC) const
        { return tetras_.begin() + (s == PosTetraC ? pos_tetra_begin_ : 0); }
    const_tetra_iterator tetra_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? tetras_.begin() + pos_tetra_begin_ : tetras_.end(); }
    ///@}
    /// Random-access to the vertices.
    ///@{
    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + (s == PosTetraC ? pos_vertex_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + neg_vertex_end_ : vertexes_.end(); }
    ///@}

    friend std::ostream& operator<< <> (std::ostream&, const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>&); ///< Debug-output to a stream (dumps all members)
    friend void write_paraview_vtu <>(std::ostream&, const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>&, TetraSignEnum);  ///< Debug-output to a stream: VTU-format for paraview.
};


///\brief Partition the principal lattice of a tetra t (n intervals on each edge) according to the sign of a levelset function ls. This class computes the triangles of the resulting piecewise trianglular interface.
class SurfacePatchCL
{
    typedef LatticePartitionTypesNS::TriangleT               TriangleT;
    typedef LatticePartitionTypesNS::TriangleContT           TriangleContT;
    typedef LatticePartitionTypesNS::const_triangle_iterator const_triangle_iterator;

    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

  private:
    TriangleContT triangles_; ///< All triangles of the interface.
    std::vector<bool> is_boundary_triangle_; ///< True, iff the triangle is a face of one of the tetras of the principal lattice.

    VertexContT vertexes_;

  public:
    /// Empty default-interface

    ///\brief Computes the piecewise triangular interface for the principal lattice with num_intervals on each edge of the reference-tetra given the level set values in ls.
    void partition_principal_lattice (Uint num_intervals, const GridFunctionCL<>& ls);

    /// True, iff the triangle is a face of one of the tetras of the principal lattice.
    ///@{
    bool is_boundary_triangle (Uint i) const { return is_boundary_triangle_[i]; }
    bool is_boundary_triangle (const_triangle_iterator it) const { return is_boundary_triangle_[it - triangles_.begin()]; }
    ///@}

    /// Random-access to the tetras and vertices.
    ///@{
    const_triangle_iterator triangle_begin () const { return triangles_.begin(); }
    const_triangle_iterator triangle_end   () const { return triangles_.end(); }
    const_vertex_iterator vertex_begin () const { return vertexes_.begin(); }
    const_vertex_iterator vertex_end   () const { return vertexes_.end(); }
    ///@}

    friend void write_paraview_vtu (std::ostream&, const SurfacePatchCL&); ///< Debug-output to a stream: VTU-format for paraview.
};

} // end of namespace DROPS

#include "geom/subtriangulation.tpp"

#endif