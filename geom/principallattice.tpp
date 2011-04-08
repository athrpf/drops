/// \file principallattice.tpp
/// \brief The principal lattice on the reference tetra obtained by slicing with planes parallel to the faces.
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

#include "misc/utils.h"

namespace DROPS {

template <class LocalFET, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const LocalFET& ls, const DomainT& dom, ResultIterT result_iterator)
{
    return std::transform( dom.vertex_begin(), dom.vertex_end(), result_iterator, ls);
}

template <class LocalFET, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const LocalFET& ls, const DomainT& dom, TetraSignEnum s, ResultIterT result_iterator)
{
    return std::transform( dom.vertex_begin( s), dom.vertex_end( s), result_iterator, ls);
}

template <class LocalFET, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (const LocalFET& ls, const DomainT& dom, ResultContT& result_container)
{
    result_container.resize( dom.vertex_size());
    std::transform( dom.vertex_begin(), dom.vertex_end(), sequence_begin( result_container), ls);
    return result_container;
}

template <class LocalFET, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (const LocalFET& ls, const DomainT& dom, TetraSignEnum s, ResultContT& result_container)
{
    result_container.resize( dom.vertex_size( s));
    std::transform( dom.vertex_begin( s), dom.vertex_end( s), sequence_begin( result_container), ls);
    return result_container;
}

} // end of namespace DROPS