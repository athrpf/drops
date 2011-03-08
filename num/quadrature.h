/// \file quadrature.h
/// \brief numerical integration at the interface
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

#ifndef DROPS_QUADRATURE_H
#define DROPS_QUADRATURE_H

#include "misc/container.h"
#include "geom/subtriangulation.h"

#include <valarray>"

namespace DROPS {

class CompositeQuad2DomainCL
{
  public:
    typedef std::vector<BaryCoordCL> VertexContT;
    typedef VertexContT::const_iterator const_vertex_iterator;

    typedef std::valarrayL<double> WeightContT;

  private:
    VertexContT vertexes_;
    size_t pos_begin_;  ///< begin of the subsequence of positive tetras

    WeightContT weights_;

  public:
    CompositeQuad2DomainCL () : pos_begin_( 0), weights_( 0) {}
    void assign (const TetraPartitionCL& partition);

    size_t size (TetraSignEnum s= AllTetraC) const
        { return vertex_end( s) - vertex_begin( s); }

    const WeightContT& weights () const { return weights_; }

    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + ( s == PosTetraC ? pos_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + pos_begin_ : vertexes_.end(); }

    template <class GridFunT>
    typename ValueHelperCL<GridFunT>::value_type quad (const GridFunT<ValueT& f, double absdet, TetraSignEnum s= AllTetraC);
    template <class GridFunT>
    void quad (const GridFunT& f, double absdet, typename ValueHelperCL<GridFunT>::value_type& pos_int, typename ValueHelperCL<GridFunT>::value_type& neg_int);
};


class CompositeQuad5DomainCL
{
  public:
    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

    typedef std::valarrayL<double> WeightContT;

  private:
    VertexContT vertexes_;
    size_t pos_begin_;  ///< begin of the subsequence of vertices in positive tetras

    WeightContT weights_;

  public:
    CompositeQuad5DomainCL () : pos_begin_( 0), weights_( 0) {}
    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
    CompositeQuad5DomainCL (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition)
        : pos_begin_( 0), weights_( 0) { assign( partition); }
    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
    void assign (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition);

    size_t size (TetraSignEnum s= AllTetraC) const
        { return vertex_end( s) - vertex_begin( s); }

    const WeightContT& weights () const { return weights_; }

    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + ( s == PosTetraC ? pos_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + pos_begin_ : vertexes_.end(); }

    template <class GridFunT>
    typename ValueHelperCL<GridFunT>::value_type quad (const GridFunT& f, double absdet, TetraSignEnum s= AllTetraC);
    template <class GridFunT>
    void quad (const GridFunT& f, double absdet, typename ValueHelperCL<GridFunT>::value_type& pos_int, typename ValueHelperCL<GridFunT>::value_type& neg_int);
};

} // end of namespace DROPS

#include "num/quadrature.tpp"

#endif