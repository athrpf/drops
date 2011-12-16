/// \file   result.h
/// \brief  Result of a nearest-neighbor search using a kd-tree
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



#ifndef RESULT_H
#define RESULT_H

#include "misc/kd-tree/kd_tree_utils.h"
#include "misc/kd-tree/metric.h"

namespace DROPS{
namespace KDTree{

    namespace internal{
        /* ******************************************************************** */
        /*  N E A R E S T  N E I G H   C L A S S                                */
        /* ******************************************************************** */

        /// \brief Class for storing a single nearest neighbor
        /** The nearest neighbor is specified by its coordinates and the distance to 
            the query point.
            \tparam T   type of cooprinates, e.g., double of float
        */
        template <typename T>
        class NearestNeighCL
        {
        private:  // ------- member variables ------- 
            size_t p_index;                 ///< index of the nearest neighbor
            T      p_distance;              ///< distance to the nearest neighbor

        public:  // ------- member functions ------- 
            /// \brief constructor
            NearestNeighCL( const size_t& idx, const T& d) 
                : p_index( idx), p_distance(d) {}

            /// \name get distance
            //@{
            inline const T& distance() const { return p_distance; }
            inline       T& distance()       { return p_distance; }
            inline size_t   get_idx()  const { return p_index; }
            //@}

            /// \brief Update the point
            inline void update( const size_t& idx, const T& d);
        };
    }


    /* ******************************************************************** */
    /*  N E A R E S T  N E I G H  R E S U L T   C L A S S                   */
    /* ******************************************************************** */

    /// \brief Storing all nearest neighbors
    /** The nearest neighbors are stored in a list (std::vector), such that the 
        nearest neighbor is at the first position and the m-th nearest neighbor
        at the last position.
        \tparam T   type of cooprinates, e.g., double of float
        \tparam K   dimension
    */
    template <typename T>
    class NearestNeighResultCL : public std::vector< internal::NearestNeighCL<T> >
    {
    public:  // ------- type definitions ------- 
        typedef internal::NearestNeighCL<T>     value_type;     ///< value type, i.e., \ref internal::NearestNeighCL
        typedef std::vector< value_type >       base;           ///< base class
        typedef typename base::iterator         iterator;       ///< iterator of all nearest neighbors
        typedef typename base::const_iterator   const_iterator; ///< constant iterator of all nearest neighbors

    private:  // ------- member variables ------- 
        size_t p_nnn;                                           ///< number of nearest neighbors

    public:  // ------- member functions ------- 
        /// \name constructors (explicit, so no constructor of the base class is allowed)
        //@{
        explicit NearestNeighResultCL( T const * point, const T& d, const size_t nnn)
            : p_nnn( nnn) { update( point, d); }
        explicit NearestNeighResultCL( const size_t nnn) 
            : p_nnn( nnn) {}
        //@}

        /// \name using the following members from the base class
        //{@
        using base::size;
        using base::begin;
        using base::end;
        //@}

        inline void update( const size_t&, const T&);               ///< put a nearest neighbor into the list of m nearest neighbors
        inline T distance() const;                                  ///< return the largest distance
        inline size_t nnn() const { return p_nnn; }                 ///< return number of nearest neighbors to be searched
    };



    /* ******************************************************************** */
    /*  N E I G H  R E S U L T  C L A S S                                   */
    /* ******************************************************************** */

    /// \brief Class for storing a list of neighbors
    /** This list is, e.g., used to store all neighbors in a 
        given sphere around a query point.
        \tparam T   type of cooprinates, e.g., double of float
        \tparam K   dimension
    */
    template <typename T>
    class NeighResultCL : public std::vector< internal::NearestNeighCL<T> >
    {
    public:  // ------- type definitions ------- 
        typedef internal::NearestNeighCL<T>     value_type;     ///< value type, i.e., \ref internal::NearestNeighCL
        typedef std::vector< value_type >       base;           ///< base class
        typedef typename base::iterator         iterator;       ///< iterator of all nearest neighbors
        typedef typename base::const_iterator   const_iterator; ///< constant iterator of all nearest neighbors

    private:  // ------- member variables ------- 
        T p_radius;                                             ///< radius aropund the query point

    public:  // ------- member functions ------- 
        /// \name constructors (explicit, so no constructor of the base class is allowed)
        //@{
        explicit NeighResultCL( const T& d) 
            : p_radius(d) {}
        //@}

        inline T radius() const;                                    ///< return the radius
        inline T& radius();                                         ///< return the radius
        inline void update( const size_t&, const T&);               ///< update this list
    };



    /* ******************************************************************** */
    /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S     */
    /* ******************************************************************** */

    namespace internal{
        /// \brief Compare two nearest neighbors by their distance
        /** strikt ordering of the nearest neighbors by their distance.

            \return true iff distance(a,query)<distance(b,query)
        */
        template <typename T>
        inline bool operator< ( const NearestNeighCL<T>& a, const NearestNeighCL<T>& b)
        {
            return a.distance()<b.distance();
        }

        /** Store the new index and the new distance.
            \param idx  new index
            \param d    new distance
        */
        template <typename T>
        inline void NearestNeighCL<T>::update( const size_t& idx, const T& d)
        {
            if ( d<p_distance){
                p_index= idx;
                p_distance= d;
            }
        }
    }   // end of namespace internal



    /* ******************************************************************** */
    /*  N E A R E S T  N E I G H  R E S U L T   C L A S S                   */
    /* ******************************************************************** */

    /** The nearest neighbor \p point is added to the list of nearest neighbors. 
        Afterwards, the list is shrinked to the number of nearest neigbors. 

        \param idx  index of the coordinates of the nearest neighbor
        \param d    distance of the nearest neighbor to the query point
    */
    template <typename T>
    inline void NearestNeighResultCL<T>::update( const size_t& idx, const T& d)
    {
        if ( size()<p_nnn){
            base::push_back( value_type( idx, d));
            std::push_heap( this->begin(), this->end());
        }
        else {
            std::pop_heap( begin(), end());
            base::back().update( idx, d);
            std::push_heap( this->begin(), this->end());
        }
    }


    /** Get the distance of the m-th nearest neighbor to the query point. If not m
        nearest neighbors are already found, \f$ \infty \f$ is returned.

        \return distance the the m-th nearest neighbor
    */
    template <typename T>
    inline T NearestNeighResultCL<T>::distance() const
    {
        return size()<p_nnn ? std::numeric_limits<T>::max() : base::front().distance();
    }



    /* ******************************************************************** */
    /*  N E I G H  R E S U L T  C L A S S                                   */
    /* ******************************************************************** */

    /** Get the radius of the sphere for reading.
        \return radius of the sphere
    */
    template <typename T>
    inline T NeighResultCL<T>::radius() const
    {
        return p_radius;
    }


    /** Get the radius of the sphere for writing.
        \return radius of the sphere
    */
    template <typename T>
    inline T& NeighResultCL<T>::radius()
    {
        return p_radius;
    }


    /** Append the point p referenced by the index idx to the list of neighbors. 
        \param idx  index of the coordinates of the point to be added
        \param d    distance of the point to the query
    */
    template <typename T>
    inline void NeighResultCL<T>::update( const size_t& idx, const T& d) 
    { 
        base::push_back( value_type( idx, d)); 
    }
}
}
#endif
