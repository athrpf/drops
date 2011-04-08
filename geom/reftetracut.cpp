/// \file reftetracut.cpp
/// \brief Triangulation of the reference tetraeder adapted to a linear level-set function.
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

#include "geom/reftetracut.h"
#include "geom/topo.h"

#include <iostream>
#include <cstring>

namespace DROPS {

void
SignPatternTraitCL::compute_cuts ()
{
    for (Ubyte i= 0; i < NumVertsC; ++i)
        if (sign( i) == 0)
            cut_simplex_[num_root_vert_++]= i;
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

void
SignPatternTraitCL::assign (const byte ls[4])
{
    num_root_vert_= num_root_= 0;

    byte sum= 0;
    for (Ubyte i= 0; i < NumVertsC; ++i)
        sum+= (sign_[i]= ls[i]);
    if (sum == 4 || sum == -4) // optimize the case of uncut tetras
        return;

    compute_cuts ();
}

void
SignPatternTraitCL::assign (const double ls[4])
{
    num_root_vert_= num_root_= 0;

    byte sum= 0;
    for (Ubyte i= 0; i < NumVertsC; ++i)
        sum+= (sign_[i]= sign( ls[i]));
    if (sum == 4 || sum == -4) // optimize the case of uncut tetras
        return;

    compute_cuts ();
}

std::ostream&
operator<< (std::ostream& out, const SignPatternTraitCL& c)
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


bool
RefTetraPatchCL::assign (const SignPatternTraitCL& cut)
{
    for (size_= 0; size_ < num_triangles( cut); ++size_)
        triangle_[size_]= MakeTriangle( cut(size_), cut(size_ + 1), cut(size_ + 2));
    return empty();
}

Ubyte
RefTetraPartitionCL::some_non_zero_vertex (const SignPatternTraitCL& cut) const
{
    Ubyte v;
    for (v= 0; v < cut.num_zero_vertexes() && cut[v] == v; ++v)
        /*empty body*/;
    return v;
}

bool
RefTetraPartitionCL::assign (const SignPatternTraitCL& cut)
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

std::ostream&
operator<< (std::ostream& out, const RefTetraPartitionCL& c)
{
    out << c.end_ - c.begin_ << ' ' << c.tetras_ - c.begin_ << '\n';
    for (Uint i= 0; i < c.end_ - c.begin_; ++i) {
        for (Uint j= 0; j < 4; ++j)
            out << static_cast<Uint>( c.begin_[i][j]) << ' ';
        out << '\n';
    }
    return out;
}

} // end of namespace DROPS