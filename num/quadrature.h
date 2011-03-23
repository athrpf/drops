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

// Create a base class QuadDomainCL with dof_begin/dof_end, size, weight_begin, vertex_begin/vertex_end. Derive Quad2/Quad5DomainCL, CompositeQuad2/Quad5DomainCL, ExtrapolatedQuad2/Quad5DomainCL 

class QuadDomainCL; // forward declaration for the factory-methods

template <class QuadDataT, class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  const QuadDomainCL&
  make_CompositeQuadDomain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>& p);

template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  const QuadDomainCL&
  make_CompositeQuad2Domain (QuadDomainCL& qdom,
    const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& p);

template <class LocalFET>
  const QuadDomainCL&
  make_ExtrapolatedQuadDomain (QuadDomainCL& q, const LocalFET& ls, const ExtrapolationToZeroCL& extra);

/// \brief General quadrature-domain
/// A quadrature rule is defined (and implemented) as a collection of quadrature points and a corresponding collection of weights.
class QuadDomainCL
{
  public:
     /// \brief Container for barycentric coordinates of quadrature points.
    typedef LatticePartitionTypesNS::VertexContT           VertexContT;
    typedef LatticePartitionTypesNS::const_vertex_iterator const_vertex_iterator;

     ///\brief Container for the quadrature weights
    typedef CompositeQuadratureTypesNS::WeightContT           WeightContT;
    typedef CompositeQuadratureTypesNS::const_weight_iterator const_weight_iterator;

    /// Friend declaration for the factory methods; if their number becomes to big, a more elaborate factory-design is required.
    ///@{
    template <class QuadDataT, class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
      friend const QuadDomainCL&
      make_CompositeQuadDomain (QuadDomainCL&,
        const TetraPartitionCL<VertexPartitionPolicyT,VertexCutMergingPolicyT>&);

    template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
      friend const QuadDomainCL&
      make_CompositeQuad2Domain (QuadDomainCL&,
        const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>&);

    template <class QuadDataT, class LocalFET>
      friend const QuadDomainCL&
      make_ExtrapolatedQuadDomain (QuadDomainCL&, const LocalFET&, const ExtrapolationToZeroCL&);
    ///@}

  private:
    VertexContT vertexes_; ///< sequence of all vertexes; some may be used for both, the positive and the negative domain
    size_t pos_begin_; ///< begin of the subsequence of vertexes of positive tetras
    size_t neg_end_;   ///< end of the subsequence of vertexes of negative tetras

    WeightContT weights_; ///< sequence of all weights; if there are vertexes on the interface, which are used for both domains, the weights for the whole domain are appended and all_weights_begin_ > 0.
    size_t pos_weights_begin_;
    size_t all_weights_begin_;

  public:
    QuadDomainCL () ///< empty default constructor
        : pos_begin_( 0), neg_end_( 0), weights_( 0), pos_weights_begin_( 0), all_weights_begin_( 0) {}

    /// \brief sequence of the indices of the vertexes (quadrature points) for the given domain
    ///@{
    Uint dof_begin (TetraSignEnum s= AllTetraC) const
        { return s == PosTetraC ? pos_begin_ : 0; }
    Uint dof_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? neg_end_ : vertexes_.size(); }
    ///@}

    size_t size (TetraSignEnum s= AllTetraC) const ///< Number of quadrature points in the given domain
        { return dof_end( s) - dof_begin( s); }

    /// \brief Begin of the sequence of weights for integration on the given domain
    const_weight_iterator weight_begin (TetraSignEnum s= AllTetraC) const {
        return Addr( weights_) + (s == NegTetraC ? 0
            : (s == PosTetraC ? pos_weights_begin_
                              : all_weights_begin_));
    }

    /// \brief sequence of quadrature points in the given domain.
    ///@{
    const_vertex_iterator vertex_begin (TetraSignEnum s= AllTetraC) const
        { return vertexes_.begin() + (s == PosTetraC ? pos_begin_ : 0); }
    const_vertex_iterator vertex_end   (TetraSignEnum s= AllTetraC) const
        { return s == NegTetraC ? vertexes_.begin() + neg_end_ : vertexes_.end(); }
    ///@}
};

template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  inline const QuadDomainCL&
  make_CompositeQuad3Domain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& p);

template <class VertexPartitionPolicyT, class VertexCutMergingPolicyT>
  inline const QuadDomainCL&
  make_CompositeQuad5Domain (QuadDomainCL& q,
    const TetraPartitionCL<VertexPartitionPolicyT, VertexCutMergingPolicyT>& p);

template <class LocalFET>
  inline const QuadDomainCL&
  make_ExtrapolatedQuad5Domain (QuadDomainCL& q, const LocalFET& ls, const ExtrapolationToZeroCL& extra);

/// Determine, how many subdivisions of the tetra-edges are required for extrapolation on level i.
///@{
struct RombergSubdivisionCL {
    Uint operator() (Uint i) { return 1 << i; }
};

struct HarmonicSubdivisionCL {
    Uint operator() (Uint i) { return i + 1; }
};
///@}

class ExtrapolationToZeroCL
{
  private:
    typedef std::valarray<double> VecT;

    Uint num_level_; ///< number of extrapolation-levels
    VecT f0_; ///< Vector of weighting factors for the level-wise quadrature weights due to extrapolation and elimination of the linear error term.

    /// \brief Compute the coefficients of the polynomial interpolating (x_i, w_i) in the Newton-basis.
    /// The computation is performed in-place.
    void compute_divided_differences (const VecT& x, std::vector<VecT>& w);
    /// \brief Evaluate the first derivative of a polynomial and the polynomial itself given in the Newton-basis.
    /// This is nested multiplication combined with the product rule.
    void evaluate_newton_polynomial_and_derivative (const VecT& x, const std::vector<VecT>& c, double p, VecT& f, VecT& der);
    /// \brief Raise the extrapolation order by one by exploiting, that there is no linear term in the error.
    /// The equation 'linear-coefficient = 0' is used to compute the next Newton-basis-coefficient of the extrapolation-polynomial.
    void eliminate_linear_term (const VecT& x, VecT& val0, const VecT& der0);

  public:
    template <class SubdivisionT>
    ExtrapolationToZeroCL (Uint num_level, const SubdivisionT& s);
    const std::valarray<double>& weights () const { return f0_; }
    const Uint num_level () const { return num_level_; }
};

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