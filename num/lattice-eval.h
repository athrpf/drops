/// \file lattice-eval.h
/// \brief Evaluate (level-set-) functions on the domains used for two-phase quadrature.
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

#ifndef DROPS_LATTICE_EVAL_H
#define DROPS_LATTICE_EVAL_H

#include "misc/container.h"

namespace DROPS {

///\brief Evaluate the LocalP2CL-like function f on [dom.vertex_begin(), dom.vertex_end()).
/// The result is stored to the sequence starting at result_iterator.
/// LocalFET must provide double operator() (const BaryCoordCL&).
/// \return end-iterator of the sequence of written signs or reference to the container of results.
///
/// This works for PrincipalLatticeCL, TetraPartitionCL, SurfacePatchCL, QuadDomainCL and QuadDomain2DCL.
///@{
template <class LocalFET, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const LocalFET& f, const DomainT& dom, ResultIterT result_iterator);

template <class LocalFET, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const LocalFET& f, const DomainT& dom, TetraSignEnum s, ResultIterT result_iterator);
///@}

///\brief Evaluate the LocalP2CL-like function f on [dom.vertex_begin(), dom.vertex_end()).
/// The result is stored to result_container. The latter is resized to dom_vertex_size().
/// ResultContainerT may be std::valarray and its derivatives, e.g. GridFunctionCL.
/// LocalFET must provide double operator() (const BaryCoordCL&).
/// \return reference to the container of results.
///
/// This works for PrincipalLatticeCL, TetraPartitionCL, SurfacePatchCL, QuadDomainCL and QuadDomain2DCL.
///@{
template <class LocalFET, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (const LocalFET& f, const DomainT& dom, ResultContT& result_container);

template <class LocalFET, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (const LocalFET& f, const DomainT& dom, TetraSignEnum s, ResultContT& result_container);
///@}

///\brief Evaluate the plain function f on [dom.vertex_begin(), dom.vertex_end()).
/// The result is stored to the sequence starting at result_iterator.
/// \return end-iterator of the sequence of written signs or reference to the container of results.
///
/// This works for PrincipalLatticeCL, TetraPartitionCL, SurfacePatchCL, QuadDomainCL and QuadDomain2DCL.
///@{
template <class T, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraCL& tet, const DomainT& dom, double t, ResultIterT result_iterator);

template <class T, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraCL& tet, const DomainT& dom, TetraSignEnum s, double t, ResultIterT result_iterator);
///@}

///\brief Evaluate the plain function f on [dom.vertex_begin(), dom.vertex_end()).
/// The result is stored to result_container. The latter is resized to dom_vertex_size().
/// ResultContainerT may be std::valarray and its derivatives, e.g. GridFunctionCL.
/// \return reference to the container of results.
///
/// This works for PrincipalLatticeCL, TetraPartitionCL, SurfacePatchCL, QuadDomainCL and QuadDomain2DCL.
///@{
template <class T, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraCL& tet, const DomainT& dom, double t, ResultContT& result_container);

template <class T, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (T (*f)(const Point3DCL&, double), const TetraCL& tet, const DomainT& dom, TetraSignEnum s, double t, ResultContT& result_container);
///@}

///\brief Evaluate the P2EvalCL-like function f on [dom.vertex_begin(), dom.vertex_end()).
/// The result is stored to the sequence starting at result_iterator.
/// PEvalT must provide the typedef LocalFET for the corresponding LocalP2CL-like class.
/// \return end-iterator of the sequence of written signs or reference to the container of results.
///
/// This works for PrincipalLatticeCL, TetraPartitionCL, SurfacePatchCL, QuadDomainCL and QuadDomain2DCL.
///@{
template <class PEvalT, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const PEvalT& f, const TetraCL& tet, const DomainT& dom, ResultIterT result_iterator);

template <class PEvalT, class DomainT, class ResultIterT>
  inline ResultIterT
  evaluate_on_vertexes (const PEvalT& f, const TetraCL& tet, const DomainT& dom, TetraSignEnum s, ResultIterT result_iterator);
///@}
///

///\brief Evaluate the P2EvalCL-like function f on [dom.vertex_begin(), dom.vertex_end()).
/// The result is stored to result_container. The latter is resized to dom_vertex_size().
/// ResultContainerT may be std::valarray and its derivatives, e.g. GridFunctionCL.
/// \return reference to the container of results.
///
/// This works for PrincipalLatticeCL, TetraPartitionCL, SurfacePatchCL, QuadDomainCL and QuadDomain2DCL.
///@{
template <class PEvalT, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraCL& tet, const DomainT& dom, ResultContT& result_container);

template <class PEvalT, class DomainT, class ResultContT>
  inline const ResultContT&
  resize_and_evaluate_on_vertexes (const PEvalT& f, const TetraCL& tet, const DomainT& dom, TetraSignEnum s, ResultContT& result_container);
///@}


} // end of namespace DROPS

#include "num/lattice-eval.tpp"

#endif
