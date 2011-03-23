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

#include "num/discretize.h"

namespace DROPS {

template <class GridFunT, class DomainT>
  typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const DomainT& dom, TetraSignEnum s=AllTetraC)
{
          Uint begin= dom.dof_begin( s);
    const Uint end=   dom.dof_end(   s);
    typename DomainT::const_weight_iterator w_iter= dom.weight_begin( s);

    typedef typename ValueHelperCL<GridFunT>::value_type value_type;
    value_type sum= value_type();
    while (begin != end)
       sum+= (*w_iter++)*f[begin++];
    return sum*absdet;
}

template <class GridFunT, class DomainT>
  inline void
  quad (const GridFunT& f, double absdet, const DomainT& dom,
    typename ValueHelperCL<GridFunT>::value_type& neg_int,
    typename ValueHelperCL<GridFunT>::value_type& pos_int)
{
    neg_int= quad( f, absdet, dom, NegTetraC);
    pos_int= quad( f, absdet, dom, PosTetraC);
}

///\brief Helper to quad_{neg,pos}_integrand
/// Integrate a integrand, that is defined only on either the negative or the positive tetras. It does not work for standard integrands.
template <class GridFunT, class DomainT>
  typename ValueHelperCL<GridFunT>::value_type
  quad_single_domain_integrand (const GridFunT& f, double absdet, const DomainT& dom, TetraSignEnum s)
{
    typename DomainT::const_weight_iterator w_iter= dom.weight_begin( s);
          Uint begin= 0;
    const Uint end= dom.dof_end( s) - dom.dof_begin( s);

    typedef typename ValueHelperCL<GridFunT>::value_type value_type;
    value_type sum= value_type();
    while (begin != end)
       sum+= (*w_iter++)*f[begin++];
    return sum*absdet;
}

template <class GridFunT, class DomainT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_neg_integrand (const GridFunT& f, double absdet, const DomainT& dom)
{
    return quad_single_domain_integrand( f, absdet, dom, NegTetraC);
}

template <class GridFunT, class DomainT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_pos_integrand (const GridFunT& f, double absdet, const DomainT& dom)
{
    return quad_single_domain_integrand( f, absdet, dom, PosTetraC);
}

template <class QuadDataT, class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  const QuadDomainCL&
  make_CompositeQuadDomain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p)
{
    typedef TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT> TetraPartitionT;
    const Uint num_nodes= QuadDataT::NumNodesC;

    q.vertexes_.resize( 0);
    q.vertexes_.reserve( num_nodes*p.tetra_size());
    q.pos_begin_= q.neg_end_= num_nodes*p.tetra_size( NegTetraC);
    q.weights_.resize( num_nodes*p.tetra_size());
    q.all_weights_begin_= 0;
    q.pos_weights_begin_= q.pos_begin_;

    const typename TetraPartitionT::const_vertex_iterator partition_vertexes= p.vertex_begin();
    const typename QuadDomainCL::WeightContT tetra_weights( QuadDataT::Weight, num_nodes);
    Uint w_begin= 0;
    QRDecompCL<4,4> qr;
    SMatrixCL<4,4>& T= qr.GetMatrix();
    for (typename TetraPartitionT::const_tetra_iterator it= p.tetra_begin(); it != p.tetra_end();
        ++it, w_begin+= num_nodes) {
        for (int i= 0; i < 4; ++i)
            T.col( i, partition_vertexes[(*it)[i]]);
        for (Uint i= 0; i < num_nodes; ++i)
            q.vertexes_.push_back( T*QuadDataT::Node[i]);
        qr.prepare_solve();
        q.weights_[std::slice( w_begin, num_nodes, 1)]= std::fabs( qr.Determinant_R())*tetra_weights;
    }
    return q;
}

template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  const QuadDomainCL&
  make_CompositeQuad5Domain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p)
{
    return make_CompositeQuadDomain<Quad5DataCL>( q, p);
}

template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  const QuadDomainCL&
  make_CompositeQuad3Domain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p)
{
    return make_CompositeQuadDomain<Quad3DataCL>( q, p);
}


template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  const QuadDomainCL&
  make_CompositeQuad2Domain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p)
{
    typedef TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT> TetraPartitionT;

    q.neg_end_= p.vertex_size( NegTetraC) + p.tetra_size( NegTetraC);
    q.pos_begin_= p.vertex_begin( PosTetraC) - p.vertex_begin() + p.tetra_size( NegTetraC);

    q.weights_.resize( p.tetra_size( NegTetraC)  + p.vertex_size( NegTetraC) // weights for NegTetraC
                      +p.vertex_size( PosTetraC) + p.tetra_size( PosTetraC)  // weights for PosTetraC
                      +p.vertex_size() + p.tetra_size());                   // weights for AllTetraC
    q.pos_weights_begin_= p.vertex_size( NegTetraC) + p.tetra_size( NegTetraC);
    q.all_weights_begin_= q.pos_weights_begin_ + p.vertex_size( PosTetraC) + p.tetra_size( PosTetraC);

    typename QuadDomainCL::VertexContT neg_tetra_bary;
    neg_tetra_bary.reserve( p.tetra_size( NegTetraC));
    typename QuadDomainCL::VertexContT pos_tetra_bary;
    pos_tetra_bary.reserve( p.tetra_size( PosTetraC));

    const typename TetraPartitionT::const_vertex_iterator partition_vertexes= p.vertex_begin();
    QRDecompCL<4,4> qr;
    SMatrixCL<4,4>& T= qr.GetMatrix();
    for (typename TetraPartitionT::const_tetra_iterator it= p.tetra_begin(); it != p.tetra_end(); ++it) {
        for (Uint i= 0; i < NumVertsC; ++i)
            T.col( i, partition_vertexes[(*it)[i]]);
        const bool is_neg= p.sign( it) == -1;
        (is_neg ? neg_tetra_bary : pos_tetra_bary).push_back( T*Quad2DataCL::Node[4]);
        double* w= Addr( q.weights_) + (is_neg ? 0 : q.pos_weights_begin_);
        const Uint vertex_weight_begin= is_neg ? p.tetra_size( NegTetraC) : 0;
        const Uint vertex_beg= is_neg ? 0 : p.vertex_begin( PosTetraC) - p.vertex_begin();
        qr.prepare_solve();
        const double absdet= std::fabs( qr.Determinant_R());
        for (int i= 0; i < 4; ++i)
            w[(*it)[i] - vertex_beg + vertex_weight_begin]+= absdet*Quad2DataCL::Wght[0];
        const Uint tetra_weight_begin= is_neg ? 0 : p.vertex_size( PosTetraC);
        const typename TetraPartitionT::const_tetra_iterator tetra_beg=
            p.tetra_begin( is_neg ? NegTetraC : PosTetraC);
        w[it - tetra_beg + tetra_weight_begin]+= absdet*Quad2DataCL::Wght[1];
    }

    q.vertexes_.resize( 0);
    q.vertexes_.reserve( p.vertex_size() + p.tetra_size());
    std::copy( neg_tetra_bary.begin(), neg_tetra_bary.end(), std::back_inserter( q.vertexes_));
    std::copy( p.vertex_begin(), p.vertex_end(), std::back_inserter( q.vertexes_));
    std::copy( pos_tetra_bary.begin(), pos_tetra_bary.end(), std::back_inserter( q.vertexes_));

    q.weights_[std::slice( q.all_weights_begin_, q.size( NegTetraC), 1)]=
        q.weights_[std::slice( 0, q.size( NegTetraC), 1)];
    q.weights_[std::slice( q.all_weights_begin_ + q.dof_begin( PosTetraC), q.size( PosTetraC), 1)]+=
        q.weights_[std::slice( q.size( NegTetraC), q.size( PosTetraC), 1)];

    return q;
}

template <class SubdivisionT>
  ExtrapolationToZeroCL::ExtrapolationToZeroCL (Uint num_level, const SubdivisionT& s)
    : num_level_( num_level), f0_( num_level)
{
    if (num_level == 0)
        throw DROPSErrCL( "ExtrapolationToZeroCL: At least one level is needed.");

    std::vector<VecT> f( num_level);
    VecT x( num_level);
    for (Uint i= 0; i < num_level; ++i) {
        x[i]= 1./s( i);
        f[i].resize( num_level);
        f[i][i]= 1.;
    }
    compute_divided_differences( x, f);
    VecT der0( num_level);
    evaluate_newton_polynomial_and_derivative( x, f, 0., f0, der0);
    eliminate_linear_term( x, f0, der0);
}

/// \brief Multiply the weight for each level with the extrapolation factor and copy it to weights.
void
copy_weights (const std::vector<CompositeQuadratureTypesNS::WeightContT>& w_vec, const std::vector<Uint>& w_pos_begin,
    const std::valarray<double>& w_factor, CompositeQuadratureTypesNS::WeightContT& weights);

template <class QuadDataT, class LocalFET>
  const QuadDomainCL&
  make_ExtrapolatedQuadDomain (QuadDomainCL& q, Uint num_level, const LocalFET& ls, const ExtrapolationToZeroCL& extra)
{
    q.vertexes_.resize( 0);

    typedef TetraPartitionCL<SortedVertexPolicyCL, MergeCutPolicyCL> TetraPartitionT;
    typename QuadDomainCL::VertexContT pos_vertexes; // temporary container for the positive vertexes
    std::vector<QuadDomainCL::WeightContT> w_vec; // the weights for each level
    w_vec.reserve( num_level);
    std::vector<Uint> w_pos_begin; // begin of the positive weights on each level
    w_pos_begin.reserve( num_level);

    TetraPartitionT partition;
    QuadDomainCL qdom;
    std::valarray<double> ls_val; // values of the level-set function in the lattice-vertexes
    // Accumulate quadrature-points and weights for each level
    for (Uint i= 0; i < num_level; ++i) {
        const Uint num_intervals= sub( i);
        const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance( num_intervals);
        ls_val.resize( lat.num_vertexes());
        for (typename PrincipalLatticeCL::const_vertex_iterator it= lat.vertex_begin(), end= lat.vertex_end(); it != end; ++it)
            ls_val[it - lat.vertex_begin()]= ls( *it);
        partition.partition_principal_lattice( num_intervals, ls_val);
        make_CompositeQuadDomain<QuadDataT>( qdom, partition);
        std::copy( qdom.vertex_begin( NegTetraC), qdom.vertex_end( NegTetraC), std::back_inserter( q.vertexes_));
        std::copy( qdom.vertex_begin( PosTetraC), qdom.vertex_end( PosTetraC), std::back_inserter( pos_vertexes));
        w_vec.push_back( QuadDomainCL::WeightContT( qdom.weight_begin(), qdom.size()));
        w_pos_begin.push_back( qdom.size( NegTetraC));
    }
    // Setup the data for the quadrature points
    q.pos_begin_= q.neg_end_= q.vertexes_.size();
    q.vertexes_.resize( q.vertexes_.size() + pos_vertexes.size());
    std::copy( pos_vertexes.begin(), pos_vertexes.end(), q.vertexes_.begin() + q.pos_begin_);

    // Compute the extrapolated weights
    q.pos_weights_begin_= q.pos_begin_;
    q.all_weights_begin_= 0;
for (Uint i= 0; i < extra.weights().size(); ++i) {
    std::cerr << extra.weights()[i] << ' ';
}
std::cerr << std::endl;
    copy_weights( w_vec, w_pos_begin, extra.weights(), q.weights_);
    return q;
}

template <class LocalFET, class SubdivisionT>
  const QuadDomainCL&
  make_ExtrapolatedQuad5Domain (QuadDomainCL& q, Uint num_level, const LocalFET& ls, SubdivisionT sub)
{
    return make_ExtrapolatedQuadDomain<Quad5DataCL>( q, num_level, ls, sub);
}

} // end of namespace DROPS
