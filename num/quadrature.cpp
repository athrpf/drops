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
        weights[std::slice( neg_it, w_pos_begin[i], 1)]= w_factor[i]*w_vec[i][std::slice( 0, w_pos_begin[i], 1)];
        weights[std::slice( pos_it, w_vec[i].size() - w_pos_begin[i], 1)]
            = w_factor[i]*w_vec[i][std::slice( w_pos_begin[i], w_vec[i].size() - w_pos_begin[i], 1)];
        neg_it+= w_pos_begin[i];
        pos_it+= w_vec[i].size() - w_pos_begin[i];
    }
}

void
compute_divided_differences (const std::valarray<double>& x, std::vector<std::valarray<double> >& w)
{
    for (Uint j= 1; j < x.size(); ++j)
        for (Uint i= x.size() - 1 ; i >= j; --i)
            w[i]= (w[i] - w[i - 1])/(x[i] - x[i - j]);
}

void
evaluate_newton_polynomial_and_derivative (const VecT& x, const std::vector<VecT>& c, double p, VecT& f, VecT& der)
{
    f= c[x.size() - 1];
    der= 0;

    for (Uint i= x.size() - 2; i < x.size(); --i) {
        der= f + (p - x[i])*der;
        f= c[i] + (p - x[i])*f;
    }
}

void
eliminate_linear_term (const std::valarray<double>& x,
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
