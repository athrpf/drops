/// \file mgobserve.h
/// \brief Observer-base-class for MultiGridCL-changes through AdaptriangCL
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_MGOBSERVER_H
#define DROPS_MGOBSERVER_H

#include "misc/utils.h"
#ifdef _PAR
#  include "parallel/parmultigrid.h"
#endif

namespace DROPS
{

/// \brief Observer-base-class for the observer-pattern
///
/// AdapTriangCL calls these methods around multigrid-changes. These can be used
/// to save information prior to grid changes or to repair functions on a multigrid
/// after such changes. Specific behavior must be added through subclassing.
/// \todo (merge) Put a reference on MultiGridCL (or ParMultiGridCL as member?) to this class?
/// \todo Number of index (velocity, pressure, levelset) for the ParMultiGridCL should be more comfortable.
class MGObserverCL
{
#ifdef _PAR
  protected:
    ParMultiGridCL* pmg_;       ///< All parallel Observers need the ParMultiGridCL
#endif

  public:
#ifdef _PAR
    MGObserverCL(ParMultiGridCL& pmg) : pmg_(&pmg) {}
    MGObserverCL() : pmg_(0) {}

    /// Get parallel multigrid.
    ParMultiGridCL& GetPMG() { return *pmg_; }
    /// Set parallel multigrid.
    void SetPMG( ParMultiGridCL& pmg) { pmg_= &pmg; }
#endif
    virtual ~MGObserverCL () {};

    /// Called immediately before MultiGridCL::Refine() in AdapTriangCL::ModifyGridStep.
    virtual void pre_refine  ()= 0;
    /// Called immediately after MultiGridCL::Refine() in AdapTriangCL::ModifyGridStep.
    virtual void post_refine ()= 0;

    /// Called at the beginning of AdapTriangCL::UpdateTriang().
    virtual void pre_refine_sequence  ()= 0;
    /// Called at the end of AdapTriangCL::UpdateTriang().
    virtual void post_refine_sequence ()= 0;
};

} // end of namespace DROPS

#endif
