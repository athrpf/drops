/// \file accumulator.h
/// \brief Basic classes to define accumulators, e. g. for the setup of stiffness matrices
/// \author LNM RWTH Aachen: Nils Gerhard, Joerg Grande; SC RWTH Aachen:

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

#ifndef DROPS_ACCUMULATOR_H
#define DROPS_ACCUMULATOR_H

#include <vector>
#include <functional>
#include <algorithm>

namespace DROPS
{

/// \brief Base class of accumulators.
///
/// State and the action performed on each VisitedT& t is added through subclassing.
/// A common use case is obtained for VisitedT == TetraCL: the setup of global stiffness matrices and load vectors
template <class VisitedT>
class AccumulatorCL
{
  public:
    /// \brief Called to initiate a new accumulation before any call of visit.
    virtual void begin_accumulation   () {}
    /// \brief Called to finalize the accumulation after all calls of visit.
    virtual void finalize_accumulation() {}

    /// \brief Called exactly once for each element of the visited sequence.
    virtual void visit (const VisitedT& t)= 0;
    
    /// \brief Returns a pointer to a copy of the actual instantiation of this class
    virtual AccumulatorCL* clone ()= 0;

    virtual ~AccumulatorCL () {}
};

/// \brief Accumulation over sequences of TetraCL.
typedef AccumulatorCL<TetraCL> TetraAccumulatorCL;


/// \brief A tuple of accumulators plus the iteration logic.
///
/// The accumulators are stored via pointers to AccumulatorCL. The iteration logic comes from a pair of external iterators defining the sequence of VisitedT-objects to be visited.
///
/// For each visited  object t, the accumulators are called in the sequence of their registration.
template <class VisitedT>
class AccumulatorTupleCL
{
  private:
    typedef std::vector<AccumulatorCL<VisitedT>*> ContainerT;
    ContainerT accus_; ///\brief the individual accumulators

    /// \brief Calls begin_iteration for each accumulator before the iteration.
    inline void begin_iteration ();
    /// \brief Calls finalize_iteration for each accumulator after the iteration.
    inline void finalize_iteration ();

  public:
    /// \brief Registers a new accumulator.
    void push_back (AccumulatorCL<VisitedT>* p) { accus_.push_back( p); }

    /// \brief Clones the vector accus_ for every thread, please note the first thread (thread #0) gets no clone, but accus_ instead
    void clone_accus(std::vector<ContainerT> &cloned);
    
    /// \brief Deletes the Clones defined from clone_accus
    void delete_clones(std::vector<ContainerT> &cloned);

    /// \brief Calls the accumulators for each object in [begin, end).
    template <class ExternalIteratorCL>
    void operator() (ExternalIteratorCL begin, ExternalIteratorCL end);
    
    /// \brief Calls the accumulators for each object by using a multigridgraph.
    void operator() ( const MultiGridCL::IndependentTetraCT& graph);
};

template <class VisitedT>
inline void AccumulatorTupleCL<VisitedT>::begin_iteration ()
{
    std::for_each( accus_.begin(), accus_.end(), std::mem_fun( &AccumulatorCL<VisitedT>::begin_accumulation));
}

template <class VisitedT>
inline void AccumulatorTupleCL<VisitedT>::finalize_iteration ()
{
    std::for_each( accus_.begin(), accus_.end(), std::mem_fun( &AccumulatorCL<VisitedT>::finalize_accumulation));
}

template<class VisitedT>
template <class ExternalIteratorCL>
void AccumulatorTupleCL<VisitedT>::operator() (ExternalIteratorCL begin, ExternalIteratorCL end)
{
    begin_iteration();
    for ( ; begin != end; ++begin)
        std::for_each( accus_.begin(), accus_.end(), std::bind2nd( std::mem_fun( &AccumulatorCL<VisitedT>::visit), *begin));
    finalize_iteration();
}

template<class VisitedT>
void AccumulatorTupleCL<VisitedT>::operator() ( const MultiGridCL::IndependentTetraCT& graph)
{
    std::vector<ContainerT> clones(omp_get_max_threads());
    begin_iteration();
    clone_accus(clones);    
    
    for( size_t i = 0 ; i < graph.size() ; ++i)
        #pragma omp parallel for
        for( size_t j = 0 ; j < graph[i].size(); ++j)
            std::for_each( clones[omp_get_thread_num()].begin(), clones[omp_get_thread_num()].end(), std::bind2nd( std::mem_fun( &AccumulatorCL<VisitedT>::visit), *(graph[i][j])));

    finalize_iteration();
    delete_clones(clones);
}

template<class VisitedT>
void AccumulatorTupleCL<VisitedT>::clone_accus(std::vector<ContainerT>& cloned)
{
    cloned[0]= accus_;  
    for( size_t i = 1 ; i < cloned.size() ; ++i){ 
        cloned[i]= accus_;
        for( size_t j = 0 ; j < accus_.size() ; ++j){
            cloned[i][j] = accus_[j]->clone();
        }
    }
}

template<class VisitedT>   
void AccumulatorTupleCL<VisitedT>::delete_clones(std::vector<ContainerT>& cloned)
{
    for( size_t i = 1 ; i < cloned.size() ; ++i)  
        for( size_t j = 0 ; j < cloned[i].size() ; ++j)
            delete cloned[i][j]; 
}


/// \brief Accumulation over sequences of TetraCL.
typedef AccumulatorTupleCL<TetraCL> TetraAccumulatorTupleCL;

} // end of namespace DROPS

#endif
