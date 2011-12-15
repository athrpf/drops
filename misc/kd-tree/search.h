/// \file   search.h
/// \brief  Search operations on kd-trees
/// \author Oliver Fortmeier, fortmeier@sc.rwth-aachen.de

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



#ifndef SEARCH_H
#define SEARCH_H

#include "misc/kd-tree/kd_tree_utils.h"
#include "misc/kd-tree/metric.h"
#include "misc/kd-tree/tree.h"
#include "misc/kd-tree/result.h"

namespace DROPS{
namespace KDTree{

    namespace internal {

        /* ******************************************************************** */
        /*  B A S E  S E A R C H  C L A S S                                     */
        /* ******************************************************************** */

        /// \brief Base class for performing search operations on kd-trees \anchor BaseSearchCL
        /** Search operations on kd-trees implementing the strategy (and template) pattern.
            \tparam metric      which metric should be used (2 = Euclidian, 0 = sup-metric, ...)
            \tparam T           type of coordinates, e.g., double of float
            \tparam K           dimension
            \tparam BucketSize  size of the bucket
        */
        template <usint metric, typename T, usint K, int BucketSize=12>
        class BaseSearchCL
        {
        public:
            typedef TreeCL<T, K, BucketSize> tree_type;                 ///< type of the tree

        protected:  // ------- type definitions ------- 
            typedef internal::NodeCL<T,K,BucketSize> node_type;         ///< type of a node

        protected:  // ------- member variables ------- 
            const tree_type& p_tree;                                    ///< reference to the kd-tree
            T const *        p_query;                                   ///< query point

        protected:  // ------- member functions ------- 
            virtual void initialize() = 0;                              ///< initialize the search
            virtual void do_search( node_type const * const) = 0;       ///< perform the search
            virtual void finalize() = 0;                                ///< finalize the search

            BaseSearchCL& operator=( const BaseSearchCL&);              ///< dummy for getting rid of the warning

        public:  // ------- member functions ------- 
            /// \brief Building the search class
            BaseSearchCL ( const tree_type& kdtree, T const * query) 
                : p_tree( kdtree), p_query( query) {}
            virtual ~BaseSearchCL() {}

            void search();                                              ///< perform the search by the template pattern
            inline const T* const& query() const { return p_query; }    ///< get the query point
        };



        /* ******************************************************************** */
        /*  B A S E  S E A R C H  R E C U R S I V E  C L A S S                  */
        /* ******************************************************************** */

        /// \brief Implements the recursive search in a kd-tree
        /** Searching recursively, efficiently in a kd-tree.
            \tparam metric      which metric should be used (2 = Euclidian, 0 = sup-metric, ...)
            \tparam T           type of coordinates, e.g., double of float
            \tparam K           dimension
            \tparam BucketSize  size of the bucket
        */
        template <usint metric, typename T, usint K, int BucketSize=12>
        class BaseSearchRecursiveCL : public internal::BaseSearchCL<metric,T,K,BucketSize>
        {
        protected:  // ------- type definitions ------- 
            typedef internal::BaseSearchCL<metric,T,K,BucketSize> base;     ///< base class
            typedef BaseSearchRecursiveCL<metric,T,K,BucketSize>  self;     ///< this class
            typedef typename base::tree_type                      tree_type;///< type of a tree
            typedef typename base::node_type                      node_type;///< type of a node

        protected:  // ------- member functions ------- 
            virtual void update( const size_t&, const T&) = 0;              ///< update the result
            virtual T distance() const = 0;                                 ///< get the recent distance
            void do_search( node_type const * const);                       ///< search recursively

        public:  // ------- member functions ------- 
            BaseSearchRecursiveCL( const tree_type& kdtree, T const * query)
                : base( kdtree, query) {}
            virtual ~BaseSearchRecursiveCL() {}
        };
    }   // end of namespace internal



    /* ******************************************************************** */
    /*  S E A R C H  N E A R E S T  N E I G H B O R S  C L A S S            */
    /* ******************************************************************** */

    /// \brief Search recursively for nnn nearest neighbors
    /** Search recursively for nn nearest neighbors where nn is the minimum of m 
        and the number of points in the tree. Therefore, the metric metric is used.
        \tparam metric      which metric should be used (2 = Euclidian, 0 = sup-metric, ...)
        \tparam T           type of coordinates, e.g., double of float
        \tparam K           dimension
        \tparam BucketSize  size of the bucket
    */        
    template <usint metric, typename T, usint K, int BucketSize=12>
    class SearchNearestNeighborsCL : public internal::BaseSearchRecursiveCL<metric,T,K,BucketSize>
    {
    public:     // ------- type definitions ------- 
        typedef NearestNeighResultCL<T>                                result_type; ///< type of the result

    protected:  // ------- type definitions ------- 
        typedef internal::BaseSearchRecursiveCL<metric,T,K,BucketSize> base;        ///< base class
        typedef typename base::tree_type                               tree_type;   ///< type of a tree
        typedef typename base::node_type                               node_type;   ///< type of a node

    protected:  // ------- member variables ------- 
        result_type p_res;                                                          ///< the result

    protected:  // ------- member functions ------- 
        void initialize();                                                          ///< put the first nnn points in the result
        void finalize();                                                            ///< for specific metrics, determine the roots
        inline void update( const size_t& idx, const T& d) { p_res.update(idx,d); } ///< update the list of nearest neighbors
        inline T distance() const { return p_res.distance(); }                      ///< get the distance to the point which is far away

    public:  // ------- member functions ------- 
        /// \brief Build a base nearest neighbor class
        SearchNearestNeighborsCL( const tree_type& kdtree, T const * query, const size_t m)
            : base( kdtree, query), p_res( std::min( m, base::p_tree.size())) {}
        virtual ~SearchNearestNeighborsCL() {}

        /// \name getters
        //@{
        inline size_t nnn() const            { return p_res.nnn();      }           ///< return number of nearest neighbors to be searched
        const result_type& result() const    { return p_res; }                      ///< return the result
        //@}
    };



    /* ******************************************************************************** */
    /*  S E A R C H  N E A R E S T  N E I G H B O R S  B R U T E  F O R C E  C L A S S  */
    /* ******************************************************************************** */

    /// \brief Search for nnn nearest neighbors by brute force
    /** Search for nn nearest neighbors where nn is the minimum of m 
        and the number of points in the tree. Therefore, the metric metric is used.
        \tparam metric      which metric should be used (2 = Euclidian, 0 = sup-metric, ...)
        \tparam T           type of coordinates, e.g., double of float
        \tparam K           dimension
        \tparam BucketSize  size of the bucket
    */
    template <usint metric, typename T, usint K, int BucketSize=12>
    class SearchNearestNeighborsBruteForceCL : public SearchNearestNeighborsCL<metric,T,K,BucketSize>
    {
    protected:  // ------- type definitions ------- 
        typedef SearchNearestNeighborsCL<metric,T,K,BucketSize> base;           ///< base class 
        typedef typename base::tree_type                        tree_type;      ///< type of a tree
        typedef typename base::node_type                        node_type;      ///< type of a node
         

    protected:  // ------- member functions ------- 
        void do_search( node_type const * const);                               ///< search by the brute force algorithm

    public:  // ------- member functions ------- 
        /// \brief Build the brute force neareast neighbor search class
        SearchNearestNeighborsBruteForceCL( const tree_type& kdtree, T const * query, const size_t m)
            : base( kdtree, query, m) {}
    };



    /* ******************************************************************** */
    /*  S E A R C H  N E A R E S T  N E I G H B O R   C L A S S             */
    /* ******************************************************************** */

    /// \brief Search recursively for 1 nearest neighbors
    /** Search recursively for the nearest neighbor in the tree. Therefore, 
        the metric metric is used.
        \tparam metric      which metric should be used (2 = Euclidian, 0 = sup-metric, ...)
        \tparam T           type of coordinates, e.g., double of float
        \tparam K           dimension
        \tparam BucketSize  size of the bucket
    */        
    template <usint metric, typename T, usint K, int BucketSize=12>
    class SearchNearestNeighborCL : public internal::BaseSearchRecursiveCL<metric,T,K,BucketSize>
    {
    public:     // ------- type definitions ------- 
        typedef internal::NearestNeighCL<T>                            result_type; ///< type of the result

    protected:  // ------- type definitions ------- 
        typedef internal::BaseSearchRecursiveCL<metric,T,K,BucketSize> base;        ///< base class
        typedef typename base::tree_type                               tree_type;   ///< type of a tree
        typedef typename base::node_type                               node_type;   ///< type of a node

    protected:  // ------- member variables ------- 
        result_type p_res;                                                          ///< the result

    protected:  // ------- member functions ------- 
        void initialize() {};                                                       ///< do nothing
        void finalize();                                                            ///< for specific metrics, determine the root
        inline void update( const size_t& idx, const T& d) { p_res.update(idx,d); } ///< update the nearest neighbor
        inline T distance() const { return p_res.distance(); }                      ///< get the distance to the recent point

    public:  // ------- member functions ------- 
        /// \brief Build a base nearest neighbor class
        SearchNearestNeighborCL( const tree_type&, T const *);
        virtual ~SearchNearestNeighborCL() {}

        /// \name getters
        //@{
        const result_type& result() const { return p_res; }                         ///< return the result
        //@}
    };



    /* ******************************************************************** */
    /*  S E A R C H  N E I G H B O R S  C L A S S                           */
    /* ******************************************************************** */

    /// \brief Search for neighbors in a given sphere around a query point
    /** Search for neighbors around a query point within a given radius.
        Therefore, the metric metric is used.
        \tparam metric      which metric should be used (2 = Euclidian, 0 = sup-metric, ...)
        \tparam T           type of coordinates, e.g., double of float
        \tparam K           dimension
        \tparam BucketSize  size of the bucket
    */
    template <usint metric, typename T, usint K, int BucketSize=12>
    class SearchNeighborsCL : public internal::BaseSearchRecursiveCL<metric, T, K, BucketSize>
    {
    public:
        typedef NeighResultCL<T>                                       result_type; ///< type of the result

    private:  // ------- type definitions ------- 
        typedef internal::BaseSearchRecursiveCL<metric,T,K,BucketSize> base;        ///< base class
        typedef typename base::tree_type                               tree_type;   ///< type of a tree
        typedef typename base::node_type                               node_type;   ///< type of a node

    private:  // ------- member variables ------- 
        result_type p_res;                                                          ///< the result

    private:  // ------- member functions ------- 
        void initialize();                                                          ///< for specific metric, build the power of the distance
        void finalize();                                                            ///< normalize the distance

    public:  // ------- member functions ------- 
        SearchNeighborsCL( const tree_type& kdtree, T const * query, const T d)
            : base( kdtree, query), p_res( d) {}

        /// \name getters
        //@{
        inline const T* const& query() const { return base::query();     }           ///< get the query point
        inline void update( const size_t& idx, const T& d){ p_res.update( idx,d);}   ///< update the result
        inline T distance() const            { return p_res.radius();    }           ///< return the search radius
        inline size_t nn() const             { return p_res.size();      }           ///< return number of neighbors of a query point
        inline const result_type& result() const    { return p_res;      }           ///< return the result
        //@}
    };



    /* ******************************************************************** */
    /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S     */
    /* ******************************************************************** */

    namespace internal {
        // ======================================================================
        // B A S E  S E A R C H  C L
        // ======================================================================

        /** This is the template for performing search operations. That is,
            first, the data structures are initialized, then, the search is 
            performed, and finally, the search is finalized.
        */
        template <usint metric, typename T, usint K, int BucketSize>
        void BaseSearchCL<metric,T,K,BucketSize>::search()
        {
            initialize();
            do_search( p_tree.root());
            finalize();
        }


        // ======================================================================
        // B A S E  S E A R C H  R E C U R S I V E  C L
        // ======================================================================
        /** Search efficiently by a recursive method. Therefore, only nodes are traversed
            whose bounding box intersects the given sphere.
            \param node the node to be traversed
        */
        template <usint metric, typename T, usint K, int BucketSize>
        void BaseSearchRecursiveCL<metric,T,K,BucketSize>::do_search( node_type const * const node)
        {
            // determine if the bounding box intersects the search sphere
            typedef typename MetricTraitsCL<metric>::metric_type metric_type;
            const internal::BoundingBoxCL<T,K>& bb= node->bounding_box();
            const bool intersects= metric_type::template intersects<T,K>( base::query(), distance(), bb);
            if ( intersects){                                                       // bounding box intersects the sphere around the nearest neighbor
                if ( !node->isLeaf()){                                              // internal node
                    if ( base::query()[node->sdim()]<=node->sval()){                // left is closer, so search left first
                        self::do_search( node->left());
                        self::do_search( node->right());
                    }
                    else{                                                           // right is closer, so search right first
                        self::do_search( node->right());
                        self::do_search( node->left());
                    }
                }
                else{                                                               // leaf node
                    const internal::BucketCL<BucketSize>& bucket= node->bucket();
                    size_t idx=NoIdx;
                    for ( int i=0; i<BucketSize && (idx=bucket[i])!=NoIdx; ++i){    // iterate over bucket and update the result
                        const T dist= metric_type::template distance<T,K>( base::query(), base::p_tree.addr(idx));
                        if ( dist<distance()) {                                     // this point is closer to the query point
                            update( idx, dist);
                        }
                    }
                }
            }
        }

    }   // end of namespace internal


    // ======================================================================
    // S E A R C H  N E A R E S T  N E I G H B O R S  C L
    // ======================================================================

    /** Put the first nnn points of the kd-tree into the result class. */
    template <usint metric, typename T, usint K, int BucketSize>
    void SearchNearestNeighborsCL<metric, T, K, BucketSize>::initialize() 
    {
        typedef typename internal::MetricTraitsCL<metric>::metric_type metric_type;
        for ( size_t i=0; i<nnn(); ++i){
            const T dist_to_point= metric_type::template distance<T,K>( base::query(), base::p_tree.addr(i));
            p_res.update( i, dist_to_point);
        }
    }


    /** For specific metric, the root of the distances has to be computed. */
    template <usint metric, typename T, usint K, int BucketSize>
    void SearchNearestNeighborsCL<metric, T, K, BucketSize>::finalize() 
    {
        if ( metric!=0 || metric!=1){
            const T divisor= metric!=0 ? metric : 1;
            for ( typename result_type::iterator it( p_res.begin()); it!=p_res.end(); ++it){
                it->distance()= metric==2 ? std::sqrt(it->distance()) : std::pow(it->distance(), static_cast<T>(1.0)/divisor);
            }
        }
    }


    // ======================================================================
    // S E A R C H  N E A R E S T  N E I G H B O R  B R U T E  F O R C E  C L
    // ======================================================================

    /** Search not-efficiently for the nearest neighbor by a brute force algorithm. */
    template <usint metric, typename T, usint K, int BucketSize>
    void SearchNearestNeighborsBruteForceCL<metric, T, K, BucketSize>::do_search( node_type const * const) 
    {
        typedef typename internal::MetricTraitsCL<metric>::metric_type metric_type;
        for ( size_t i=base::nnn(); i<base::p_tree.size(); ++i){
            const T dist_to_point= metric_type::template distance<T,K>( base::query(), base::p_tree.addr(i));
            if ( dist_to_point<base::p_res.distance())
                base::p_res.update( i, dist_to_point);
        }
    }



    /* ******************************************************************** */
    /*  S E A R C H  N E A R E S T  N E I G H B O R   C L A S S             */
    /* ******************************************************************** */

    /** Build the search structure, the first nearest neighbor is set to the 
        first point in the tree.
        \param kdtree   the tree to be searched
        \param query    the query point
    */
    template <usint metric, typename T, usint K, int BucketSize>
    SearchNearestNeighborCL<metric,T,K,BucketSize>::SearchNearestNeighborCL( const tree_type& kdtree, T const * query)
            : base( kdtree, query), 
              p_res( 0, internal::MetricTraitsCL<metric>::metric_type::template distance<T,K>( base::query(), base::p_tree.addr(0))) 
    {}


    /** Denormalize the distance. */
    template <usint metric, typename T, usint K, int BucketSize>
    void SearchNearestNeighborCL<metric,T,K,BucketSize>::finalize()
    {
        if ( metric!=0 || metric!=1){
            const T divisor= metric!=0 ? metric : 1;
            p_res.distance()= metric==2 ? std::sqrt( p_res.distance()) : std::pow(p_res.distance(), static_cast<T>(1.0)/divisor);
        }
    }


    // ======================================================================
    // S E A R C H  N E I G H B O R S  C L
    // ======================================================================

    /** For specific metric, normalize the distance */
    template <usint metric, typename T, usint K, int BucketSize>
    void SearchNeighborsCL<metric, T, K, BucketSize>::initialize()
    {
        if ( metric!=0 && metric!=1)
            p_res.radius()= std::pow( p_res.radius(), static_cast<T>(metric));
    }


    /** Renormalize the distances. */
    template <usint metric, typename T, usint K, int BucketSize>
    void SearchNeighborsCL<metric, T, K, BucketSize>::finalize()
    {
        if ( metric!=0 && metric!=1)
            p_res.radius()= std::pow( p_res.radius(), static_cast<T>(1)/static_cast<T>(metric));
    }

}
}
#endif
