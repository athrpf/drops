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

#include <valarray>

namespace DROPS {

/// Integration of full-sized integrands, which have size( AllTetraC) components.
///@{
/// \brief Integrate on the negative, the positive or all tetras.
template <class GridFunT, class DomainT>
  typename ValueHelperCL<GridFunT>::value_type
  quad (const GridFunT& f, double absdet, const DomainT& dom, TetraSignEnum s=AllTetraC);

/// \brief Integrate on the negative and the positive tetras.
template <class GridFunT, class DomainT>
  inline void
  quad (const GridFunT& f, double absdet, const DomainT& dom,
    typename ValueHelperCL<GridFunT>::value_type& neg_int,
    typename ValueHelperCL<GridFunT>::value_type& pos_int);
///@}

/// Integration of small integrands, which have size( NegTetraC) or size( PosTetraC) components
///@{
/// \brief Integrate an integrand, that is defined only on the negative tetras. It does not work for full-sized integrands. Use quad for the latter.
template <class GridFunT, class DomainT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_neg_integrand (const GridFunT& f, double absdet, const DomainT& dom);

/// \brief Integrate an integrand, that is defined only on the positive tetras. It does not work for standard integrands. Use quad for the latter.
template <class GridFunT, class DomainT>
  inline typename ValueHelperCL<GridFunT>::value_type
  quad_pos_integrand (const GridFunT& f, double absdet, const DomainT& dom);
///@}

namespace CompositeQuadratureTypesNS {

typedef std::valarray<double> WeightContT;
typedef const double* const_weight_iterator;

} // end of namespace DROPS::CompositeQudratureTypesNS

/// \brief Prototype of a quadrature-domain
///
/// Reimplementation of the quadrature rule of Quad5CL
/// A quadrature rule is defined (and implemented) as a collection of quadrature points and a corresponding collection of weights.
// class QuadDomainCL
// {
//   public:
//     typedef LatticePartitionTypesNS::VertexContT VertexContT; ///\brief Container for barycentric coordinates of quadrature points.
//     typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;
// 
//     typedef CompositeQuadratureTypesNS::WeightContT WeightContT; ///< Container for the quadrature weights
//     typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;
// 
//   public:
//     QuadDomainCL ();
//     template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
//     QuadDomainCL (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition)
//         { assign( partition); }
//     template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
//     void assign (const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& partition);
// 
//     /// \brief sequence of the indices of the quadrature points for the given domain
//     ///@{
//     Uint dof_begin (TetraSignEnum s= AllTetraC) const;
//     Uint dof_end   (TetraSignEnum s= AllTetraC) const;
//     ///@}
// 
//     size_t size (TetraSignEnum s= AllTetraC) const; ///< Number of quadrature points in the given domain
// 
//     /// \brief Begin of the sequence of weights for integration on the given domain
//     const_weight_iterator weight_begin (TetraSignEnum s= AllTetraC) const; 
// 
//     /// \brief sequence of quadrature points in the given domain.
//     ///@{
//     const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const;
//     const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const;
//     ///@}
// };

class CompositeQuad2DomainCL
{
  public:
    typedef LatticePartitionTypesNS::VertexContT VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

    typedef CompositeQuadratureTypesNS::WeightContT WeightContT;
    typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;

  private:
    VertexContT vertexes_;
    size_t pos_begin_; ///< begin of the subsequence of vertexes of positive tetras
    size_t neg_end_;   ///< end of the subsequence of vertexes of negative tetras

    WeightContT neg_weights_;
    WeightContT pos_weights_;
    WeightContT all_weights_;

  public:
    CompositeQuad2DomainCL () : pos_begin_( 0), neg_end_( 0), neg_weights_( 0), pos_weights_( 0), all_weights_( 0){}
    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
    CompositeQuad2DomainCL (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition)
        { assign( partition); }
    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
    void assign (const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& partition);

    Uint dof_begin (TetraSignEnum s= AllTetraC) const
        { return s == PosTetraC ? pos_begin_ : 0; }
    Uint dof_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? neg_end_ : vertexes_.size(); }

    size_t size (TetraSignEnum s= AllTetraC) const
        { return dof_end( s) - dof_begin( s); }

    const_weight_iterator weight_begin (TetraSignEnum s= AllTetraC) const {
        return (s == NegTetraC ? Addr( neg_weights_)
                               : (s == PosTetraC ? Addr( pos_weights_)
                                                 : Addr( all_weights_)));
    }

    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + ( s == PosTetraC ? pos_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + neg_end_ : vertexes_.end(); }
};


class CompositeQuad5DomainCL
{
  public:
    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

    typedef CompositeQuadratureTypesNS::WeightContT           WeightContT;
    typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;

  private:
    VertexContT vertexes_;
    size_t pos_begin_;  ///< begin of the subsequence of vertices in positive tetras

    WeightContT weights_;

  public:
    CompositeQuad5DomainCL () : pos_begin_( 0), weights_( 0) {}
    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
    CompositeQuad5DomainCL (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition)
        { assign( partition); }
    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
    void assign (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition);

    Uint dof_begin (TetraSignEnum s= AllTetraC) const
        { return s == PosTetraC ? pos_begin_ : 0; }
    Uint dof_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? pos_begin_ : vertexes_.size(); }

    size_t size (TetraSignEnum s= AllTetraC) const
        { return dof_end( s) - dof_begin( s); }

    const_weight_iterator weight_begin (TetraSignEnum s= AllTetraC) const
        { return Addr( weights_) + (s == PosTetraC ? pos_begin_ : 0); }

    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + ( s == PosTetraC ? pos_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + pos_begin_ : vertexes_.end(); }
};


// class ExtrapolatedQuad5DomainCL
// {
//   public:
//     typedef LatticePartitionTypesNS::VertexContT           VertexContT;
//     typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;
// 
//     typedef CompositeQuadratureTypesNS::WeightContT           WeightContT;
//     typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;
// 
//   private:
//     VertexContT vertexes_;
//     size_t pos_begin_;  ///< begin of the subsequence of vertices in positive tetras
// 
//     WeightContT weights_;
// 
//   public:
//     ExtrapolatedQuad5DomainCL () : pos_begin_( 0), weights_( 0) {}
//     template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
//     ExtrapolatedQuad5DomainCL (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition)
//         { assign( partition); }
//     template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
//     void assign (const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& partition);
// 
//     Uint dof_begin (TetraSignEnum s= AllTetraC) const
//         { return s == PosTetraC ? pos_begin_ : 0; }
//     Uint dof_end   (TetraSignEnum s= AllTetraC) const
//         { return s == NegTetraC ? pos_begin_ : vertexes_.size(); }
// 
//     size_t size (TetraSignEnum s= AllTetraC) const
//         { return dof_end( s) - dof_begin( s); }
// 
//     const_weight_iterator weight_begin (TetraSignEnum s= AllTetraC) const
//         { return Addr( weights_) + (s == PosTetraC ? pos_begin_ : 0); }
// 
//     const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
//         { return vertexes_.begin() + ( s == PosTetraC ? pos_begin_ : 0); }
//     const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
//         { return s == NegTetraC ? vertexes_.begin() + pos_begin_ : vertexes_.end(); }
// };

///\brief Write the sign of the levelset function ls in the quadrature points [vert_begin, vert_end) to the sequence beginning at begin.
/// \return end-iterator of the sequence of written signs
// template<class sign_iterator>
//   sign_iterator
//   copy_levelset_sign (const LocalP2CL<>& ls,
//     LatticePartitionTypesNS::const_vertex_iterator vert_begin,
//     LatticePartitionTypesNS::const_vertex_iterator vert_end,
//     sign_iterator begin)
// {
//     while (vert_begin != vert_end) {
//         *begin++= sign( ls( *vert_begin++));
//     }
//     return begin;
// }



} // end of namespace DROPS

#include "num/quadrature.tpp"

#endif