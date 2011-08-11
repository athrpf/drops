/// \file quadrature.cpp
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

#include "num/quadrature.h"

namespace DROPS {

QuadDomainCL&
QuadDomainCL::operator= (const QuadDomainCL& q)
{
    vertexes_=  q.vertexes_;
    pos_begin_= q.pos_begin_;
    neg_end_=   q.neg_end_;

    weights_.resize( q.weights_.size());
    weights_=           q.weights_;
    pos_weights_begin_= q.pos_weights_begin_;
    all_weights_begin_= q.all_weights_begin_;

    return *this;
}

const QuadDomainCL&
make_CompositeQuad2Domain (QuadDomainCL& q, const TetraPartitionCL& p)
{
    q.neg_end_= p.vertex_size( NegTetraC) + p.tetra_size( NegTetraC);
    q.pos_begin_= p.vertex_begin( PosTetraC) - p.vertex_begin() + p.tetra_size( NegTetraC);

    q.weights_.resize( p.tetra_size( NegTetraC)  + p.vertex_size( NegTetraC) // weights for NegTetraC
                      +p.vertex_size( PosTetraC) + p.tetra_size( PosTetraC)  // weights for PosTetraC
                      +p.vertex_size() + p.tetra_size());                    // weights for AllTetraC
    q.pos_weights_begin_= p.vertex_size( NegTetraC) + p.tetra_size( NegTetraC);
    q.all_weights_begin_= q.pos_weights_begin_ + p.vertex_size( PosTetraC) + p.tetra_size( PosTetraC);

    QuadDomainCL::VertexContT neg_tetra_bary; // positive barycenters
    neg_tetra_bary.reserve( p.tetra_size( NegTetraC));
    QuadDomainCL::VertexContT pos_tetra_bary; // negative barycenters
    pos_tetra_bary.reserve( p.tetra_size( PosTetraC));

    const TetraPartitionCL::const_vertex_iterator partition_vertexes= p.vertex_begin();
    SMatrixCL<4,4> T;
    for (TetraPartitionCL::const_tetra_iterator it= p.tetra_begin(); it != p.tetra_end(); ++it) {
        for (Uint i= 0; i < NumVertsC; ++i)
            T.col( i, partition_vertexes[(*it)[i]]);
        const bool is_neg= p.sign( it) == -1;
        (is_neg ? neg_tetra_bary : pos_tetra_bary).push_back( T*Quad2DataCL::Node[4]);
        double* w= Addr( q.weights_) + (is_neg ? 0 : q.pos_weights_begin_);
        const Uint vertex_weight_begin= is_neg ? p.tetra_size( NegTetraC) : 0;
        const Uint vertex_beg= is_neg ? 0 : p.vertex_begin( PosTetraC) - p.vertex_begin();
        const double absdet= std::fabs( VolFrac( T));
        for (int i= 0; i < 4; ++i)
            w[(*it)[i] - vertex_beg + vertex_weight_begin]+= absdet*Quad2DataCL::Wght[0];
        const Uint tetra_weight_begin= is_neg ? 0 : p.vertex_size( PosTetraC);
        const TetraPartitionCL::const_tetra_iterator tetra_beg=
            p.tetra_begin( is_neg ? NegTetraC : PosTetraC);
        w[it - tetra_beg + tetra_weight_begin]+= absdet*Quad2DataCL::Wght[1];
    }

    q.vertexes_.resize( 0);
    q.vertexes_.reserve( p.vertex_size() + p.tetra_size());
    std::copy( neg_tetra_bary.begin(), neg_tetra_bary.end(), std::back_inserter( q.vertexes_));
    std::copy( p.vertex_begin(), p.vertex_end(), std::back_inserter( q.vertexes_));
    std::copy( pos_tetra_bary.begin(), pos_tetra_bary.end(), std::back_inserter( q.vertexes_));

    q.weights_[std::slice( q.all_weights_begin_, q.vertex_size( NegTetraC), 1)]=
        q.weights_[std::slice( 0, q.vertex_size( NegTetraC), 1)];
    q.weights_[std::slice( q.all_weights_begin_ + q.dof_begin( PosTetraC), q.vertex_size( PosTetraC), 1)]+=
        q.weights_[std::slice( q.vertex_size( NegTetraC), q.vertex_size( PosTetraC), 1)];

    return q;
}

void
copy_weights (const std::vector<CompositeQuadratureTypesNS::WeightContT>& w_vec, const std::vector<Uint>& w_pos_begin,
    const std::valarray<double >& w_factor, CompositeQuadratureTypesNS::WeightContT& weights)
{
    Uint s= 0, s_neg= 0;
    for (Uint i= 0; i < w_vec.size(); ++i) {
        s+= w_vec[i].size();
        s_neg+= w_pos_begin[i];
    }
    weights.resize( s);

    Uint neg_it= 0, pos_it= s_neg;
    for (Uint i= 0; i < w_vec.size(); ++i) {
        const Uint j= w_vec.size() - 1 - i; // To access the extrapolation-weights, which are ordered from coarse to fine level
        weights[std::slice( neg_it, w_pos_begin[i], 1)]= w_factor[j]*w_vec[i][std::slice( 0, w_pos_begin[i], 1)];
        weights[std::slice( pos_it, w_vec[i].size() - w_pos_begin[i], 1)]
            = w_factor[j]*w_vec[i][std::slice( w_pos_begin[i], w_vec[i].size() - w_pos_begin[i], 1)];
        neg_it+= w_pos_begin[i];
        pos_it+= w_vec[i].size() - w_pos_begin[i];
    }
}

QuadDomain2DCL&
QuadDomain2DCL::operator= (const QuadDomain2DCL& q)
{
    vertexes_= q.vertexes_;

    weights_.resize( q.weights_.size());
    weights_= q.weights_;

    return *this;
}

void
copy_weights_surface (const std::vector<CompositeQuadratureTypesNS::WeightContT>& w_vec,
    const std::valarray<double >& w_factor, CompositeQuadratureTypesNS::WeightContT& weights)
{
    Uint s= 0;
    for (Uint i= 0; i < w_vec.size(); ++i)
        s+= w_vec[i].size();
    weights.resize( s);

    Uint it= 0;
    for (Uint i= 0; i < w_vec.size(); ++i) {
        const Uint j= w_vec.size() - 1 - i; // To access the extrapolation-weights, which are ordered from coarse to fine level
        weights[std::slice( it, w_vec[i].size(), 1)]= w_factor[j]*w_vec[i];
        it+= w_vec[i].size();
    }
}

void
ExtrapolationToZeroCL::compute_divided_differences (const std::valarray<double>& x, std::vector<std::valarray<double> >& w)
{
    for (Uint j= 1; j < x.size(); ++j)
        for (Uint i= x.size() - 1 ; i >= j; --i)
            w[i]= (w[i] - w[i - 1])/(x[i] - x[i - j]);
}

void
ExtrapolationToZeroCL::evaluate_newton_polynomial_and_derivative (const VecT& x, const std::vector<VecT>& c, double p, VecT& f, VecT& der)
{
    f= c[x.size() - 1];
    der= 0;

    for (Uint i= x.size() - 2; i < x.size(); --i) {
        der= f + (p - x[i])*der;
        f= c[i] + (p - x[i])*f;
    }
}

void
ExtrapolationToZeroCL::eliminate_linear_term (const std::valarray<double>& x,
    std::valarray<double>& val0, const std::valarray<double>& der0)
{
    // Evaluate the next Newton-basis-polynomial and its derivative at 0.
    double omega_n= 1.;
    double der_omega_n= 0;

    for (Uint i= x.size() - 1; i < x.size(); --i) {
        der_omega_n= omega_n - x[i]*der_omega_n;
        omega_n*= -x[i];
    }
    // Eliminate the linear term from the interpolation-polynomial
    val0-= (omega_n/der_omega_n)*der0;

}

} // end of namespace DROPS
