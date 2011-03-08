/// \file quadrature.tpp
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

namespace DROPS {

template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
void CompositeQuad2DomainCL::assign (const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p)
{
    typedef TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT> TetraPartitionT;
    const Uint num_nodes= Quad2DataCL::NumNodesC;

    vertexes_.resize( 0);
    vertexes_.reserve( p.vertex_size() + p.tetra_size());
    std::copy( p.vertex_begin( NegTetraC), p.vertex_end( NegTetraC), std::back_inserter( vertexes_));

    weights_.resize(   p.vertex_size() + p.tetra_size());

    const TetraPartitionCL::const_vertex_iterator partition_vertexes= p.vertex_begin();
    const WeightContT tetra_weights( Quad5DataCL::Weight, num_nodes);
    Uint w_begin= 0;
    QRDecompCL<4,4> qr;
    SMatrixCL<4,4>& T= qr.GetMatrix();
    for (TetraPartitionCL::const_tetra_iterator it= p.tetra_begin(); it != p.tetra_end(); ++it, w_begin+= num_nodes) {
        for (int i= 0; i < 4; ++i)
            T.col( i)= partition_vertexes[(*it)[i]];
        for (Uint i= 0; i < num_nodes; ++i)
            vertexes_.push_back( T*Quad2DataCL::Node[i]);
        qr.prepare_solve();
        weights_[std::slice( w_begin, num_nodes, 1)]= std::fabs( qr.Determinant_R())*tetra_weights;
    }
}


template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
void CompositeQuad5DomainCL::assign (const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p)
{
    typedef TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT> TetraPartitionT;
    const Uint num_nodes= Quad5DataCL::NumNodesC;

    vertexes_.resize( 0);
    vertexes_.reserve( num_nodes*p.tetra_size());
    pos_begin_= num_nodes*p.tetra_size( NegTetraC);
    weights_.resize( num_nodes*p.tetra_size());

    const typename TetraPartitionT::const_vertex_iterator partition_vertexes= p.vertex_begin();
    const WeightContT tetra_weights( Quad5DataCL::Weight, num_nodes);
    Uint w_begin= 0;
    QRDecompCL<4,4> qr;
    SMatrixCL<4,4>& T= qr.GetMatrix();
    for (typename TetraPartitionT::const_tetra_iterator it= p.tetra_begin(); it != p.tetra_end(); ++it, w_begin+= num_nodes) {
        for (int i= 0; i < 4; ++i)
            T.col( i, partition_vertexes[(*it)[i]]);
        for (Uint i= 0; i < num_nodes; ++i)
            vertexes_.push_back( T*Quad5DataCL::Node[i]);
        qr.prepare_solve();
        weights_[std::slice( w_begin, num_nodes, 1)]= std::fabs( qr.Determinant_R())*tetra_weights;
    }
}


template <class GridFunT>
typename ValueHelperCL<GridFunT>::value_type CompositeQuad5DomainCL::quad (const GridFunT& f, double absdet, TetraSignEnum s)
{
    const WeightContT& theweights= weights();
          Uint begin= s == PosTetraC ? size( NegTetraC) : 0;
    const Uint end=   s == NegTetraC ? size( NegTetraC) : size( AllTetraC);

    typedef typename ValueHelperCL<GridFunT>::value_type value_type;
    value_type sum= value_type();
    for (; begin < end; ++begin)
       sum+= theweights[begin]*f[begin];
    return sum*absdet;
}

template <class GridFunT>
void CompositeQuad5DomainCL::quad (const GridFunT& f, double absdet, typename ValueHelperCL<GridFunT>::value_type& neg_int, typename ValueHelperCL<GridFunT>::value_type& pos_int)
{
    const WeightContT& theweights= weights();
    pos_int= neg_int= typename ValueHelperCL<GridFunT>::value_type();
    const Uint pos_begin= size( NegTetraC);

    for (Uint i= 0; i < pos_begin; ++i)
       neg_int+= theweights[i]*f[i];
    neg_int*= absdet;

    const Uint end= size( AllTetraC);
    for (Uint i= pos_begin; i < end; ++i)
       pos_int+= theweights[i]*f[i];
    pos_int*= absdet;
}


template <class GridFunT>
typename ValueHelperCL<GridFunT>::value_type CompositeQuad2DomainCL::quad (const GridFunT<ValueT>& f, double absdet, TetraSignEnum s)
{
    const Uint begin= s == PosTetraC ? size( NegTetraC) : 0;
    const Uint end=   s == NegTetraC ? size( NegTetraC) : size( AllTetraC);
    const WeightContT& weights= weights();

    ValueT sum= ValueT();
    for (; begin < end; ++begin)
       sum+=  weights[begin]*f[begin]
    return sum*absdet;
}

template <class GridFunT>
void CompositeQuad2DomainCL::quad (const GridFunT<ValueT>& f, double absdet, typename ValueHelperCL<GridFunT>::value_type& neg_int, typename ValueHelperCL<GridFunT>::value_type& pos_int)
{
    typedef typename ValueHelperCL<GridFunT>::value_type value_type;
    const WeightContT& weights= weights();
    pos_int= neg_int= value_type();
    const Uint pos_begin= size( NegTetraC);

    for (Uint i= 0; i < pos_begin; ++i)
       neg_int+=  weights[i]*f[i];
    neg_int*= absdet;

    const Uint end= size( AllTetraC);
    for (Uint i= pos_begin; i < end; ++i)
       pos_int+=  weights[i]*f[i];
    pos_int*= abs_det;
}

} // end of namespace DROPS