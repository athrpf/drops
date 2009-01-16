//**************************************************************************
// File:    adaptriang.h                                                   *
// Content: adaptive triangulation based on position of the interface      *
//          provided by the levelset function                              *
// Author:  Sven Gross, Joerg Grande, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************

#ifndef DROPS_ADAPTRIANG_H
#define DROPS_ADAPTRIANG_H

#include "levelset/mgobserve.h"
#include "levelset/levelset.h"
#include <vector>

#ifdef _PAR
# include "parallel/parallel.h"
# include "parallel/parmultigrid.h"
# include "parallel/loadbal.h"
# include "parallel/logger.h"
#endif

namespace DROPS
{

class AdapTriangCL
{
  private:
    MultiGridCL& mg_;
#ifdef _PAR
    ParMultiGridCL&   pmg_;
    LoadBalHandlerCL& lb_;
#endif
    double width_;
    int c_level_, f_level_;
    bool modified_;

    typedef std::vector<MGObserverCL*> ObserverContT;       // type for observing the multigird dependent FE functions
    ObserverContT observer_;                                // handler for manipulate FE-functions to due grid changes

    template <class DistFctT>
    double GetValue( const DistFctT& dist, const VertexCL& v) { return dist.val( v); }
    template <class DistFctT>
    double GetValue( const DistFctT& dist, const EdgeCL& e)   { return dist.val( e); }
    template <class DistFctT>
    double GetValue( const DistFctT& dist, const TetraCL& t)  { return dist.val( t, 0.25, 0.25, 0.25); }
    double GetValue( scalar_fun_ptr dist, const VertexCL& v)  { return dist( v.GetCoord() ); }
    double GetValue( scalar_fun_ptr dist, const EdgeCL& e)    { return dist( GetBaryCenter( e) ); }
    double GetValue( scalar_fun_ptr dist, const TetraCL& t)   { return dist( GetBaryCenter( t) ); }

#ifndef _PAR
    template <class DistFctT>
    bool ModifyGridStep( DistFctT&);
#else
    template <class DistFctT>
    bool ModifyGridStep( DistFctT&, bool lb=true);
#endif
    // One step of grid change; returns true if modifications were necessary,
    // false, if nothing changed.

    /// \name Call handlers (MGObserverCL) to manipulate FE-functions
    // @{
    /// \brief Tell Observer, that MG will be refined (and migration will be performed)
    void notify_pre_refine () {
#ifndef _PAR
        for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
            (*obs)->pre_refine();
#else
        if (!observer_.empty()){
            pmg_.DeleteVecDesc();
            for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
              (*obs)->pre_refine();
        }
#endif
    }

    /// \brief Tell Observer, that MG has been refined (and migration was performed)
    void notify_post_refine () {
        for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
            (*obs)->post_refine();
#ifdef _PAR
        if ( !observer_.empty() ){
            pmg_.DelAllUnkRecv();
            pmg_.DeleteRecvBuffer();
        }
#endif
    }

    /// \brief Tell Observer, that a sequence of refinements (and migrations) will take place
    void notify_pre_refine_sequence() {
        for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
            (*obs)->pre_refine_sequence();
    }

    /// \brief Tell Observer, that a sequence of refinements (and migrations) has taken place
    void notify_post_refine_sequence() {
        for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
            (*obs)->post_refine_sequence();
    }
    //@}

  public:
#ifndef _PAR
    AdapTriangCL( MultiGridCL& mg, double width, int c_level, int f_level)
      : mg_(mg), width_(width), c_level_(c_level), f_level_(f_level), modified_(false)
      { Assert( 0<=c_level && c_level<=f_level, "AdapTriangCL: Levels are cheesy.\n", ~0); }
#else
    AdapTriangCL( ParMultiGridCL& pmg, LoadBalHandlerCL& lb, double width, int c_level, int f_level)
      : mg_(pmg.GetMG()), pmg_(pmg), lb_(lb), width_(width), c_level_(c_level), f_level_(f_level), modified_(false)
      { Assert( 0<=c_level && c_level<=f_level, "AdapTriangCL: Levels are cheesy.\n", ~0); }

    /// \brief Get a reference onto the parallel MultiGrid
    ParMultiGridCL& GetPMG() { return pmg_; }

    /// \brief Get a constant reference onto the parallel MultiGrid
    const ParMultiGridCL& GetPMG() const { return pmg_; }
#endif

    /// \brief Make initial triangulation according to a distance function
    template <class DistFctT>
      void MakeInitialTriang (DistFctT&);

    /// \brief Coupling of all necessary steps update the triangulation according to a levelset function
    void UpdateTriang (const LevelsetP2CL& ls);

    /// \brief Check if the triangulation has been modified within last update
    bool WasModified () const { return modified_; }

    /// \brief Get a reference onto the MultiGrid
    const MultiGridCL& GetMG() const { return mg_; }

    /// \brief Push back a handler for FE-functions to apply changes due to grid modifications
    void push_back (MGObserverCL* o) { observer_.push_back( o); }
};

} // end of namespace DROPS

#include "levelset/adaptriang.tpp"

#endif
