/// \file
/// \brief Observer-base-class for MultiGridCL-changes through AdaptriangCL
/// \author Sven Gross, Joerg Grande, Patrick Esser IGPM RWTH Aachen
///         Oliver Fortmeier SC RWTH Aachen

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
