/// \file adaptriang.h
/// \brief adaptive triangulation based on position of the interface provided by the levelset function
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
    ParMultiGridCL*  pmg_;
    LoadBalHandlerCL lb_;
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
            pmg_->DeleteVecDesc();
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
            pmg_->DelAllUnkRecv();
            pmg_->DeleteRecvBuffer();
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
    AdapTriangCL(  MultiGridCL& mg, double width, int c_level, int f_level, __UNUSED__ int refineStrategy = 1)
      : mg_( mg),
#ifdef _PAR
      pmg_( ParMultiGridCL::InstancePtr()), lb_( mg_),
#endif
      width_(width), c_level_(c_level), f_level_(f_level), modified_(false)
      {
        Assert( 0<=c_level && c_level<=f_level, "AdapTriangCL: Levels are cheesy.\n", ~0);
#ifdef _PAR
        pmg_->AttachTo( mg_);
        if (refineStrategy>=0)
            lb_.DoInitDistribution( ProcCL::Master());
        else
            refineStrategy*=-1;
        switch ( refineStrategy) {
            case 0 : lb_.SetStrategy( NoMig);     break;
            case 1 : lb_.SetStrategy( Adaptive);  break;
            case 2 : lb_.SetStrategy( Recursive); break;
        }
#endif
      }

#ifdef _PAR
    /// \brief Get a reference onto the parallel MultiGrid
    ParMultiGridCL& GetPMG() { return *pmg_; }

    /// \brief Get a constant reference onto the parallel MultiGrid
    const ParMultiGridCL& GetPMG() const { return *pmg_; }

    /// \brief Get a reference onto the LoadBalHandlerCL
    LoadBalHandlerCL& GetLb() { return lb_; }

    /// \brief Get a constant reference onto the LoadBalHandlerCL
    const LoadBalHandlerCL& GetLb() const { return lb_; }
#endif

    void SetWidth       (double width) { width_  = width; }
    void SetCoarseLevel (int c_level)  { c_level_= c_level; }
    void SetFineLevel   (int f_level)  { f_level_= f_level; }

    /// \brief Make initial triangulation according to a distance function
    template <class DistFctT>
      void MakeInitialTriang (DistFctT&);

    /// \brief Coupling of all necessary steps update the triangulation according to a levelset function
    void UpdateTriang (const LevelsetP2CL& ls);

    /// \brief Check if the triangulation has been modified within last update
    bool WasModified () const { return modified_; }

    /// \name Get a reference onto the MultiGrid
    //@{
    const MultiGridCL& GetMG() const { return mg_; }
    MultiGridCL&       GetMG()       { return mg_; }
    //@}

    /// \brief Push back a handler for FE-functions to apply changes due to grid modifications
    void push_back (MGObserverCL* o)
    {
        observer_.push_back( o);
#ifdef _PAR
        observer_.back()->SetPMG( *pmg_);
#endif
    }
};

} // end of namespace DROPS

#include "levelset/adaptriang.tpp"

#endif
