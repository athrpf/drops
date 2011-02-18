/// \file combinatorialcut.cpp
/// \brief tests the PrincipalLattice-class
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

#include "geom/principallattice.h"
#include "misc/container.h"
#include "num/discretize.h"

#include <iostream>
#include <sstream>

namespace DROPS {

///\brief Represents the reference tetra, which is cut by a linear level set function ls. The values of the latter are prescribed on the vertices.
class SignPatternTraitCL {
  private:
    Ubyte num_root_vert_;  ///< number of vertices, where the level set function is zero.
    Ubyte num_root_;       ///< number of roots of the level set function; invariant: num_root_vert <= num_root
    byte sign_[4];         ///< Sign of the level set function of the vertices; \f$\in\{-1,0,1\}\f$
    Ubyte cut_simplex_[4]; ///< local number with respect to the reference tetra of the object on the cut: [0,num_root_vert): vertex numbers in (0..3); [num_root_vert, num_root): edge numbers in (0..5). Both parts are sorted in increasing order.
    Ubyte cut_simplex_rep_[4]; ///< local number of the object on the cut: (0..9)

    bool is_3d () const { return num_root_vert_ == 4; }

  public:
    SignPatternTraitCL () : num_root_vert_( 4), num_root_( 4) {} ///< Uninitialized default state
    SignPatternTraitCL (const double ls[4]) { assign( ls); } ///< Assign a sign pattern on the vertices; throws if ls is identically 0.
    void assign (const double ls[4]); ///< Assign a sign pattern on the vertices; throws if ls is identically 0.

    bool empty () const { return num_root_ == 0; } ///< True, iff there is no intersection.
    bool is_2d () const { return num_root_ > 2; }  ///< True, iff the intersection has positive area.
    bool no_zero_vertex () const { return num_root_vert_ == 0; } ///< True, iff there is no vertex, in which ls vanishes.

    Ubyte num_cut_simplexes () const { return num_root_; } ///< Number of edges and vertices with a root of ls.
    Ubyte num_zero_vertexes () const { return num_root_vert_; } ///< Number of vertices of the tetra that are roots of ls.

    byte sign (int i) const { return sign_[i]; } ///< -1,0,1; sign of vertex i.

    /// Return local number of edges/verts with a root of ls. For edges, [] returns a edge number in 0..5, and () returns an extended vertex number in 4..9.
    ///@{
    Ubyte operator[] (int i) const { return cut_simplex_[i]; }
    Ubyte operator() (int i) const { return cut_simplex_rep_[i]; }
    ///@}

    friend std::ostream& operator<< (std::ostream&, const SignPatternTraitCL&); ///< Debug-output to a stream (dumps all members)
};

///\brief Helper to instance_idx below
inline Ubyte sign_plus_one (double d)
{
    return d > 0. ? 2 : (d < 0. ? 0 : 1);
}

///\brief Return the array-index for the possible 3^4 sign-patterns on the vertices of a tetra.
inline Ubyte instance_idx (const double ls[4])
{
    return  sign_plus_one( ls[0])*27 + sign_plus_one( ls[1])*9
          + sign_plus_one( ls[2])* 3 + sign_plus_one( ls[3]);
}

///\brief The triangles of the intersection of the reference-tetra with a linear levelset-function.
///
/// The class memoizes used sign-patterns if the triangulations are accessed via the instance( ls)-function. Individual instances may still be constructed (useful for debugging).
class RefTetraSurfacePatchCL {
  public:
    typedef SArrayCL<Ubyte, 3> TriangleT; ///< the vertices of a triangle of the cut: the tetra's vertices are denoted by 0..3, the edge-cuts by edge-num + 4, which is in 4..9.
    typedef const TriangleT* const_triangle_iterator;
    typedef       TriangleT*       triangle_iterator;

  private:
    TriangleT triangle_[2];      ///< at most two triangles
    Ubyte size_;                 ///< number of triangles
    Ubyte is_boundary_triangle_; ///< true if the triangle is one of the tetra's faces.

    Ubyte num_triangles (const SignPatternTraitCL& cut) const { return cut.is_2d() ? cut.num_cut_simplexes() - 2 : 0; }
    TriangleT MakeTriangle (Ubyte v0, Ubyte v1, Ubyte v2) const { return MakeSArray( v0, v1, v2); }

  public:
    RefTetraSurfacePatchCL () : size_( -1), is_boundary_triangle_( 0) {} ///< Uninitialized default state
    RefTetraSurfacePatchCL (const SignPatternTraitCL& cut) { assign( cut); } ///< Initialize with sign pattern on the vertices
    bool assign (const SignPatternTraitCL& cut); ///< Assign a sign pattern on the vertices; returns the value of empty()

    bool  is_initialized () const { return size_ <= 2; } ///< True after assign(...)

    static inline const RefTetraSurfacePatchCL&
    instance (const double ls[4]); ///< Recommended access to the triangles for a given sign-pattern; memoizes the result

    bool is_boundary_triangle () const { return is_boundary_triangle_ == 1; } ///< true, iff the triangle is one of the tetra's faces.

    bool  empty () const { return size_ == 0; } ///< true, iff the area of the intersection is 0.
    size_t size () const { return size_; }      ///< Number of triangles, 0, 1, or 2

    /// Random-access to the triangles
    ///@{
    const_triangle_iterator triangle_begin () const { return triangle_; }
    const_triangle_iterator triangle_end   () const { return triangle_ + size_; }
    ///@}
};

enum TetraSignEnum { AllTetraC, NegTetraC, PosTetraC }; ///< Designates the subset of tetras one is interested in RefTetraPartitionCL and TetraPartitionCL

///\brief The tetras partition the positive and negative part of the reference-tetra with respect to a linear levelset-function ls.
///
/// The class memoizes used sign-patterns if the triangulations are accessed via the instance( ls)-function. Individual instances may still be constructed (useful for debugging).
///
/// Layout of the tetra-sequence: tetras_ has at most 6 members, of which at most 3 are positive and 3 are negative.
/// [tetras_ ... <= ... begin_ ... <= ... tetras_ + 3 ... <= ... end_ ... <= ...tetras_ + 6]
/// [begin_..tetras_+3) are the negative tetras, [tetras_ +3..end_) are the positive tetras.
class RefTetraPartitionCL
{
  public:
    typedef SArrayCL<Ubyte, 4> TetraT; ///< representation of a tetra of the partition via its vertices: (0..3): original vertices; (4..9): original edges + 4
    typedef const TetraT* const_tetra_iterator;
    typedef       TetraT*       tetra_iterator;

  private:
    TetraT tetras_[6]; ///< at most six tetras
    tetra_iterator begin_;
    tetra_iterator end_;

    TetraT MakeTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3) const { return  MakeSArray( v0, v1, v2, v3); }

    inline void AddTetra (Ubyte v0, Ubyte v1, Ubyte v2, Ubyte v3, int sign) ///< The sequences grow away from tetras_+3
        { (sign == -1 ? *--begin_ : *end_++)= MakeTetra( v0, v1, v2, v3); }
     ///\brief e,f,g are assumed to be the equally-oriented edges of the quadrilateral faces.
    inline void AddPrism (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, Ubyte g0, Ubyte g1, int sign);

    ///\brief If the intersection is quadrilateral, this returns the first uncut edge.
    Ubyte first_uncut_edge (const SignPatternTraitCL& cut) const { return cut[0] == 1 ? 0 : (cut[1] == 2 ? 1 : 2); }
    Ubyte some_non_zero_vertex (const SignPatternTraitCL& cut) const;

  public:
    RefTetraPartitionCL () : begin_( tetras_ + 3), end_( tetras_) {} ///< Uninitialized default state
    RefTetraPartitionCL (const SignPatternTraitCL& cut) { assign( cut); } ///< Initialize with sign pattern on the vertices
    bool assign (const SignPatternTraitCL& cut); ///< Assign a sign pattern on the vertices; returns the value of is_uncut()
    bool is_initialized () const { return begin_ < end_; } ///< True after assign(...)

    static inline const RefTetraPartitionCL&
    instance (const double ls[4]); ///< Recommended access to the triangles for a given sign-pattern; memoizes the result

    bool is_uncut () const { return end_ == begin_ + 1; } ///< True, iff the partition has exactly one tetra
    int sign (const_tetra_iterator t) const { return t < tetras_ + 3 ? -1 : 1; } ///< Sign of the tetra, to which t points

    Ubyte tetra_size (TetraSignEnum s= AllTetraC) const ///< number of tetras with given sign
         { return tetra_end( s) - tetra_begin( s); }

     /// Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum
    ///@{
   const_tetra_iterator tetra_begin (TetraSignEnum s= AllTetraC) const
        { return s == PosTetraC ? tetras_ + 3 : begin_; }
    const_tetra_iterator tetra_end (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? tetras_ + 3 : end_; }
    ///@}

    friend std::ostream& operator<< (std::ostream&, const RefTetraPartitionCL&); ///< Debug-output to a stream (dumps all members)
};

void SignPatternTraitCL::assign (const double ls[4])
{
    num_root_vert_= num_root_= 0;

    int sum= 0;
    for (Ubyte i= 0; i < NumVertsC; ++i)
        sum+= (sign_[i]= ls[i] < 0. ? -1 : (ls[i] > 0. ? 1 : 0));
    if (sum == 4 || sum == -4) // optimize the case of uncut tetras
        return;

    for (Ubyte i= 0; i < NumVertsC; ++i)
        if (sign( i) == 0) cut_simplex_[num_root_vert_++]= i;
    num_root_= num_root_vert_;
    for (Ubyte i= 0; i < NumEdgesC; ++i)
        if (sign( VertOfEdge( i, 0))*sign( VertOfEdge( i, 1)) == -1)
            cut_simplex_[num_root_++]= i;
    std::memcpy( cut_simplex_rep_, cut_simplex_, 4*sizeof(byte));
    for (int i= num_root_vert_; i < num_root_; ++i)
        cut_simplex_rep_[i]+= NumVertsC;

    if (is_3d())
        throw DROPSErrCL( "InterfacePatchCL::assign: found 3-dim. zero level set, grid is too coarse!");
}

std::ostream& operator<< (std::ostream& out, const SignPatternTraitCL& c)
{
    out << static_cast<int>( c.num_root_vert_) << ' ' << static_cast<int>( c.num_root_) << '\n';
    for (int i= 0; i < 4; ++i)
        out << static_cast<int>( c.sign_[i]) << ' ';
    out << '\n';
    for (int i= 0; i < 4; ++i)
        out << static_cast<int>( c.cut_simplex_[i]) << ' ';
    out << '\n';
    for (int i= 0; i < 4; ++i)
        out << static_cast<int>( c.cut_simplex_rep_[i]) << ' ';
    return out << '\n';
}

inline bool
RefTetraSurfacePatchCL::assign (const SignPatternTraitCL& cut)
{
    for (size_= 0; size_ < num_triangles( cut); ++size_)
        triangle_[size_]= MakeTriangle( cut(size_), cut(size_ + 1), cut(size_ + 2));
    return empty();
}

inline const RefTetraSurfacePatchCL&
RefTetraSurfacePatchCL::instance (const double ls[4])
{
    static RefTetraSurfacePatchCL instances_[81]; // 81 = 3^4 = all possible sign-patterns on the vertices

    RefTetraSurfacePatchCL& instance= instances_[instance_idx ( ls)];
    if ( !instance.is_initialized())
        instance.assign( SignPatternTraitCL( ls));
    return instance;
}

inline void
RefTetraPartitionCL::AddPrism (Ubyte e0, Ubyte e1, Ubyte f0, Ubyte f1, Ubyte g0, Ubyte g1, int sign)
{
    AddTetra( e0, e1, f1, g1, sign);
    AddTetra( e0, f0, f1, g1, sign);
    AddTetra( e0, f0, g0, g1, sign);
}

inline Ubyte
RefTetraPartitionCL::some_non_zero_vertex (const SignPatternTraitCL& cut) const
{
    Ubyte v;
    for (v= 0; cut[v] == v && v < cut.num_zero_vertexes(); ++v)
        /*empty body*/;
    return v;
}

std::ostream& operator<< (std::ostream& out, const RefTetraPartitionCL& c)
{
    out << c.end_ - c.begin_ << ' ' << c.tetras_ - c.begin_ << '\n';
    for (Uint i= 0; i < c.end_ - c.begin_; ++i) {
        for (Uint j= 0; j < 4; ++j)
            out << static_cast<Uint>( c.begin_[i][j]) << ' ';
        out << '\n';
    }
    return out;
}

bool RefTetraPartitionCL::assign (const SignPatternTraitCL& cut)
{
    end_= begin_= tetras_ + 3;

    if (cut.empty()) { // Most common case: no cut.
        AddTetra( 0, 1, 2, 3, cut.sign( 0));
    }
    else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level
        if (cut.num_cut_simplexes() == 3) { // triangular cut: a tetra and a remaining prism
            const Ubyte v= VertByEdge( cut[0], cut[1]);
            AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
            AddPrism( OppVertOfEdge( cut[0], v), cut(0),
                      OppVertOfEdge( cut[1], v), cut(1),
                      OppVertOfEdge( cut[2], v), cut(2),
                      -cut.sign( v));
        }
        else if (cut.num_cut_simplexes() == 4) { // quadrilateral cut: two prisms
            const Ubyte e= first_uncut_edge( cut);
            const Ubyte f= OppEdge( e);
            AddPrism( VertOfEdge(e, 0), VertOfEdge( e, 1),
                      cut(0), cut(1),
                      cut(2), cut(3),
                      cut.sign(VertOfEdge( e, 0)));
            AddPrism( VertOfEdge(f, 0), VertOfEdge( f, 1),
                      cut(0), cut(2),
                      cut(1), cut(3),
                      cut.sign(VertOfEdge( f, 0)));
        }
    }
    else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next common case: there are cut edges, and also 1 or 2 vertices of the tetra with value 0 (the latter as we are in the else-part of cut.no_zero_vertex())
        if (cut.num_zero_vertexes() == 1) { // triangular cut through a vertex: a tetra and a remaining pyramid with quadrilateral base
            const Ubyte e= cut[1], f= cut[2];
            const Ubyte v= VertByEdge( e, f);
            AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
            const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1) : VertOfEdge( e, 0);
            const Ubyte opp_v_in_f= v == VertOfEdge( f, 0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
            // the pyramid
            AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
            AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
        }
        else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2 vertexes: two tetras
            const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
            const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
            AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
            AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
        }
    }
    else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
        AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
    return is_uncut();
}

inline const RefTetraPartitionCL&
RefTetraPartitionCL::instance (const double ls[4])
{
    static RefTetraPartitionCL instances_[81]; // 81 = 3^4 = all possible sign-patterns on the vertices

    RefTetraPartitionCL& instance= instances_[instance_idx (ls)];
    if ( !instance.is_initialized())
        instance.assign( SignPatternTraitCL( ls));
    return instance;
}



///\brief Partition the principal lattice of a tetra t (n intervals on each edge) according to the sign of a levelset function ls.
///
/// The sub-tetras, which are cut by the interface are triangulated to match the interface. The values of the level set function on the vertices of the principal lattice must be prescribed. The sequence of all tetrahedra contains the negative tetrahedra as initial subsequence.
class TetraPartitionCL
{
  public:
    typedef SArrayCL<Uint, 4> TetraT; ///< representation of a tetra of the partition via its vertices: index in the vertex_-array
    typedef std::vector<TetraT> TetraContT;
    typedef TetraContT::const_iterator const_tetra_iterator;
    typedef std::vector<BaryCoordCL> VertexContT;
    typedef VertexContT::const_iterator const_vertex_iterator;

  private:
    TetraContT tetras_; ///< All tetras of the partition.
    size_t pos_begin_;  ///< begin of the subsequence of positive tetras

    VertexContT vertexes_; ///< All vertices of the partition. 

  public:
    TetraPartitionCL () : pos_begin_( 0) {} ///< Empty default-cut

    ///\brief Computes the partition of the principal lattice with num_intervals on each edge of the reference-tetra given the level set values in ls.
    void partition_principal_lattice (Uint num_intervals, const GridFunctionCL<>& ls);

    size_t tetra_size (TetraSignEnum s= AllTetraC) const ///< number of tetras with given sign
         { return tetra_end( s) - tetra_begin( s); }

    /// Random-access to the tetras: all tetras, or negative and positive tetras separately, see TetraSignEnum.
    ///@{
    const_tetra_iterator tetra_begin (TetraSignEnum s= AllTetraC) const
        { return tetras_.begin() + ( s == PosTetraC ? pos_begin_ : 0); }
    const_tetra_iterator tetra_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? tetras_.begin() + pos_begin_ : tetras_.end(); }
    ///@}
    /// Random-access to the vertices.
    ///@{
    const_vertex_iterator vertex_begin () const { return vertexes_.begin(); }
    const_vertex_iterator vertex_end   () const { return vertexes_.end(); }
    ///@}

    friend std::ostream& operator<< (std::ostream&, const TetraPartitionCL&); ///< Debug-output to a stream (dumps all members)
    friend void write_paraview_vtu (std::ostream&, const TetraPartitionCL&, TetraSignEnum= AllTetraC);  ///< Debug-output to a stream: VTU-format for paraview.
};

///\brief Partition the principal lattice of a tetra t (n intervals on each edge) according to the sign of a levelset function ls. This class computes the triangles of the resultin piecewise trianglular interface.
class SurfacePatchCL
{
    typedef SArrayCL<Uint, 3> TriangleT; ///< representation of a triangle of the interface via its vertices: index in the vertex_-array
    typedef std::vector<TriangleT> TriangleContT;
    typedef TriangleContT::const_iterator const_triangle_iterator;
    typedef std::vector<BaryCoordCL> VertexContT;
    typedef VertexContT::const_iterator const_vertex_iterator;

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

// GridFunDomainT: const TetraCL* tetra () const (may return 0); iterator (ueber BaryCoordCL) vertex_begin(), vertex_end()
// TetraSubTriangCL + tetra_begin(), tetra_end()
// PartitionedTetraSubTriangCL + tetra_begin() und tetra_end() haben ein default enum-arg : {AllTetra, PosTetra, NegTetra}
// TetraSurfacePatchCL: GridFunDomainT + triangle_begin(), triangle_end()
// template <class GridFunT, class GridFunDomainT, class ResultContT> Evaluate (ResultContT& r, GridFunT& f, const GridFunDomainT& dom) { EvaluatorCL<GridFunT, GridFunDomainT>::evaluate( r, f, dom); }

///\todo The roots on the edges are inserted multiple times into vertexes_: use an associative map to cache them.
void TetraPartitionCL::partition_principal_lattice (Uint num_intervals, const GridFunctionCL<>& ls)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( num_intervals);

    tetras_.resize( 0);
    pos_begin_= 0;
    vertexes_.resize( 0);

    vertexes_.reserve( lat.num_vertexes());
    std::copy( lat.vertex_begin(), lat.vertex_end(), std::back_inserter( vertexes_));

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
                    const double edge_bary1_cut= lset[v0]/(lset[v0] - lset[v1]); // the root of the level set function on the edge
                    tet[j]= vertexes_.size();
                    vertexes_.push_back( ConvexComb( edge_bary1_cut, vertexes_[(*lattice_tet)[v0]], vertexes_[(*lattice_tet)[v1]]));
                }
            }
            (cut.sign( it) == -1 ? tetras_ : loc_tetras).push_back( tet);
        }
    }
    pos_begin_= tetras_.size();
    std::copy( loc_tetras. begin(), loc_tetras.end(), std::back_inserter( tetras_));
}

///\todo The roots on the edges are inserted multiple times into vertexes_: use an associative map to cache them.
void SurfacePatchCL::partition_principal_lattice (Uint num_intervals, const GridFunctionCL<>& ls)
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( num_intervals);
    PrincipalLatticeCL::const_vertex_iterator lattice_verts= lat.vertex_begin();

    triangles_.resize( 0);
    is_boundary_triangle_.resize( 0);

    double lset[4];
    Uint loc_vert_num;
    TriangleT tri;
    for (PrincipalLatticeCL::const_tetra_iterator lattice_tet= lat.tetra_begin(), lattice_end= lat.tetra_end(); lattice_tet != lattice_end; ++lattice_tet) {
        for (Uint i= 0; i < 4; ++i)
            lset[i]= ls[(*lattice_tet)[i]]; 
        const RefTetraSurfacePatchCL& cut= RefTetraSurfacePatchCL::instance( lset);
        if (cut.empty()) continue;
        for (RefTetraSurfacePatchCL::const_triangle_iterator it= cut.triangle_begin(), end= cut.triangle_end(); it != end; ++it) {
            for (Uint j= 0; j < 3; ++j) {
                tri[j]= vertexes_.size();
                loc_vert_num= (*it)[j];
                if (loc_vert_num < 4) {
                    const Uint lattice_vert_num= (*lattice_tet)[loc_vert_num];
                    vertexes_.push_back( lattice_verts[lattice_vert_num]);
                }
                else { // Cut vertex
                    const Ubyte v0= VertOfEdge( loc_vert_num - 4, 0),
                                v1= VertOfEdge( loc_vert_num - 4, 1);
                    const double edge_bary1_cut= lset[v0]/(lset[v0] - lset[v1]);
                    vertexes_.push_back( ConvexComb( edge_bary1_cut, lattice_verts[(*lattice_tet)[v0]], lattice_verts[(*lattice_tet)[v1]]));
                }
            }
            triangles_.push_back( tri);
            is_boundary_triangle_.push_back( cut.is_boundary_triangle());
        }
    }
}

std::ostream& operator<< (std::ostream& out, const TetraPartitionCL& t)
{
    out << t.tetras_.size() << ' ' << t.pos_begin_ << '\n';
    std::copy( t.tetras_.begin(), t.tetras_.end(), std::ostream_iterator<TetraPartitionCL::TetraT>( out));
    out << '\n';
    std::copy( t.vertexes_.begin(), t.vertexes_.end(), std::ostream_iterator<BaryCoordCL>( out)) ;
    return out;
}


void write_paraview_vtu (std::ostream& file_, const TetraPartitionCL& t, TetraSignEnum s)
{
    file_ << "<?xml version=\"1.0\"?>"  << '\n'
          << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << '\n'
          << "<UnstructuredGrid>"   << '\n';

    file_<< "<Piece NumberOfPoints=\""<< t.vertexes_.size() <<"\" NumberOfCells=\""<< t.tetra_size( s) << "\">";
    file_<< "\n\t<Points>"
         << "\n\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << "ascii\">\n\t\t";
    for(TetraPartitionCL::const_vertex_iterator it= t.vertexes_.begin(), end= t.vertexes_.end(); it != end; ++it) {
        file_ << it[0][1] << ' ' << it[0][2] << ' ' << it[0][3] << ' ';
    }
    file_<< "\n\t\t</DataArray> \n"
         << "\t</Points>\n";

    file_   << "\t<Cells>\n"
            << "\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\""
            <<"ascii\">\n\t\t";
    std::copy( t.tetra_begin( s), t.tetra_end( s), std::ostream_iterator<TetraPartitionCL::TetraT>( file_));
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

void write_paraview_vtu (std::ostream& file_, const SurfacePatchCL& t)
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
    std::copy( t.triangles_.begin(), t.triangles_.end(), std::ostream_iterator<SurfacePatchCL::TriangleT>( file_));
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

} // end of namespace DROPS


void test_tetra_cut ()
{
    DROPS::GridFunctionCL<> ls( 4);
    ls[0]= -1.; ls[1]= 0.; ls[2]= 0.; ls[3]= 0.;
    DROPS::TetraPartitionCL tet;
    // tet.partition_principal_lattice ( 1, ls);
    // std::cerr << tet;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              if (i == 0 && j == 0 && k == 0 && l == 0) continue;
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << ls[0] << ' ' << ls[1] << ' ' << ls[2] << ' ' << ls[3] << std::endl;
              DROPS::RefTetraPartitionCL cut( static_cast<double*>(&ls[0]));
              DROPS::SignPatternTraitCL comb_cut( static_cast<double*>(&ls[0]));
              tet.partition_principal_lattice ( 1, ls);
//              if (c == 5) {
//                  std::cerr << comb_cut << std::endl;
//                  std:: cerr << cut << std::endl;
//                  std::cerr << tet << std::endl << std::endl;
//              }
              std::ostringstream name;
              name << "hallo" << c << ".vtu";
              std::ofstream file( name.str().c_str());
              DROPS::write_paraview_vtu( file, tet);
          }
}

void test_cut_surface ()
{
    DROPS::GridFunctionCL<> ls( 4);
    ls[0]= -1.; ls[1]= 0.; ls[2]= 0.; ls[3]= 0.;
    DROPS::SurfacePatchCL tet;
    // tet.partition_principal_lattice ( 1, ls);
    // std::cerr << tet;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              if (i == 0 && j == 0 && k == 0 && l == 0) continue;
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << ls[0] << ' ' << ls[1] << ' ' << ls[2] << ' ' << ls[3] << std::endl;
              DROPS::RefTetraPartitionCL cut( static_cast<double*>(&ls[0]));
              DROPS::SignPatternTraitCL comb_cut( static_cast<double*>(&ls[0]));
              tet.partition_principal_lattice ( 1, ls);
              std::ostringstream name;
              name << "hallo_surf" << c << ".vtu";
              std::ofstream file( name.str().c_str());
              DROPS::write_paraview_vtu( file, tet);
          }
}

void test_principal_lattice ()
{
    for (int i= 1; i <= 4; ++i) {
        const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( i);
        std::cout << "=======================================" << lat.num_intervals() << ' ' << lat.num_vertexes() << " " << lat.num_tetras() << std::endl;
        for (DROPS::PrincipalLatticeCL::const_vertex_iterator v= lat.vertex_begin(), end= lat.vertex_end(); v != end; ++v) {
            std::cout << lat.num_intervals()*(*v) << std::endl;
        }
        std:: cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        for (DROPS::PrincipalLatticeCL::const_tetra_iterator v= lat.tetra_begin(), end= lat.tetra_end(); v != end; ++v) {
            std::cout << (*v)[0] << ' '  << (*v)[1] << ' ' << (*v)[2] << ' ' << (*v)[3] << ' ' << std::endl;
        }
    }
}

inline double sphere (const DROPS::Point3DCL& p)
{
    return p.norm() - 0.5;
}

void evaluate (DROPS::GridFunctionCL<>& dest, const DROPS::PrincipalLatticeCL& lat, double (*f)(const DROPS::Point3DCL& p))
{
    for (DROPS::PrincipalLatticeCL::const_vertex_iterator it= lat.vertex_begin(), end= lat.vertex_end(); it != end; ++it) {
        dest[it - lat.vertex_begin()]= f( DROPS::MakePoint3D( it[0][1], it[0][2], it[0][3]));
    }
}

void test_sphere_cut ()
{
    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( 10);
    DROPS::GridFunctionCL<> ls( lat.num_vertexes());
    evaluate( ls, lat, &sphere);
    DROPS::TetraPartitionCL tet;
    tet.partition_principal_lattice( 10, ls);
    std::ostringstream name;
    name << "sphere.vtu";
    std::ofstream file( name.str().c_str());
    DROPS::write_paraview_vtu( file, tet, DROPS::NegTetraC);
    file.close();

    DROPS::SurfacePatchCL surf;
    surf.partition_principal_lattice( 10, ls);
    name.str( "");
    name << "sphere_surf.vtu";
    file.open( name.str().c_str());
    DROPS::write_paraview_vtu( file, surf);
}

int main()
{
    try {
        // test_tetra_cut();
        // test_cut_surface();
        // test_principal_lattice();
        test_sphere_cut ();
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}