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
    /// \param clone_id the thread-id, in which the clone will run. Useful to locate cloned helper objects.
    virtual AccumulatorCL* clone (int clone_id)= 0;

    virtual ~AccumulatorCL () {}
};

/// \brief Accumulation over sequences of TetraCL.
typedef AccumulatorCL<TetraCL> TetraAccumulatorCL;


/// \brief A tuple of accumulators plus the iteration logic.
///
/// The accumulators are stored via pointers to AccumulatorCL. The iteration logic comes from a pair of external iterators defining the sequence of VisitedT-objects to be visited.
/// For iteration over a MultiGridCL::IndependentTetraCT, the accumulators are cloned after begin_accumulation, visit is called OpenMP-parallel, and the clones are destroyed. finalize_accumulation is called only for the original accumulators.
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

    /// \brief Clones the vector accus_ for every thread; the first thread (thread 0) gets no clone, but a copy of accus_ instead
    void clone_accus(std::vector<ContainerT>& clones);

    /// \brief Deletes the clones defined from clone_accus; obviously, accus_ is not deleted
    void delete_clones(std::vector<ContainerT>& clones);

    /// \brief Calls the accumulators for each object in [begin, end).
    template <class ExternalIteratorCL>
    void operator() (ExternalIteratorCL begin, ExternalIteratorCL end);

    /// \brief Calls the accumulators for each object by using a ColorClassesCL.
    void operator() (const ColorClassesCL& colors);
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
void AccumulatorTupleCL<VisitedT>::operator() (const ColorClassesCL& colors)
{
    begin_iteration();

    std::vector<ContainerT> clones( omp_get_max_threads());
    clone_accus( clones);
    for (ColorClassesCL::const_iterator cit= colors.begin(); cit != colors.end() ;++cit) {
#       pragma omp parallel
        {
            const int t_id= omp_get_thread_num();
            const ColorClassesCL::ColorClassT& cc= *cit;
#           pragma omp for
            for (size_t j= 0; j < cc.size(); ++j)
                std::for_each( clones[t_id].begin(), clones[t_id].end(), std::bind2nd( std::mem_fun( &AccumulatorCL<VisitedT>::visit), *cc[j]));
        }
    }
    delete_clones(clones);

    finalize_iteration();
}

template<class VisitedT>
void AccumulatorTupleCL<VisitedT>::clone_accus(std::vector<ContainerT>& clones)
{
    clones[0]= accus_;
    for (size_t i= 1; i < clones.size(); ++i) {
        clones[i].resize( accus_.size());
        std::transform( accus_.begin(), accus_.end(), clones[i].begin(), std::bind2nd( std::mem_fun( &AccumulatorCL<VisitedT>::clone), i));
    }
}

template<class VisitedT>
void AccumulatorTupleCL<VisitedT>::delete_clones(std::vector<ContainerT>& clones)
{
    for (size_t i= 1; i < clones.size(); ++i)
        for (size_t j= 0; j < clones[i].size(); ++j)
            delete clones[i][j];
}


/// \brief Accumulation over sequences of TetraCL.
typedef AccumulatorTupleCL<TetraCL> TetraAccumulatorTupleCL;

} // end of namespace DROPS

#endif
