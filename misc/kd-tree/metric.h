/// \file   metric.h
/// \brief  All methods used to implement a metric
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


#ifndef METRIC_H
#define METRIC_H

#include "misc/kd-tree/kd_tree_utils.h"
#include "misc/kd-tree/bounding_box.h"

namespace DROPS{
namespace KDTree{

    namespace internal{


        /* ********************************************************************* */
        /*  M E T R I C  C L A S S E S                                           */
        /* ********************************************************************* */
        /// \name functions that depends on the chosen metric
        //@{

        /// \brief Class for functions using the supremum metric
        struct SupMetricCL
        {
            /// \brief distance between two points
            template <typename T, usint K>
            static inline T distance( T const *, T const *);
            /// \brief does the bounding box and a sphere intersects
            template <typename T, usint K>
            static inline bool intersects( const T * const, const T&, const BoundingBoxCL<T,K>&);
        };


        /// \brief Class for functions using the ome metric
        struct OneMetricCL
        {
            /// \brief distance between two points
            template <typename T, usint K>
            static inline T distance( T const *, T const *);
            /// \brief does the bounding box and a sphere intersects
            template <typename T, usint K>
            static inline bool intersects( const T * const, const T&, const BoundingBoxCL<T,K>&);
        };


        /// \brief Class for functions using the Euclidean metric \anchor EuclideanMetricCL
        struct EuclideanMetricCL
        {
            /// \brief distance between two points
            template <typename T, usint K>
            static inline T distance( T const *, T const *);
            /// \brief does the bounding box and a sphere intersects
            template <typename T, usint K>
            static inline bool intersects( const T * const, const T&, const BoundingBoxCL<T,K>&);
        };


        /// \brief Class for functions using the Markowski metric
        template <usint metric>
        struct MarkowskiMetricCL
        {
            /// \brief distance between two points
            template <typename T, usint K>
            static inline T distance( T const *, T const *);
            /// \brief does the bounding box and a sphere intersects
            template <typename T, usint K>
            static inline bool intersects( const T * const, const T&, const BoundingBoxCL<T,K>&);
        };
        //@}


        /* ********************************************************************* */
        /*  T R A I T S  F O R  C H O S I N G  T H E  M E T R I C                */
        /* ********************************************************************* */

        /// \brief Trait for chosing the metric by a template parameter, i.e., Markowski metric
        template <usint metric>
        struct MetricTraitsCL { typedef internal::MarkowskiMetricCL<metric> metric_type; };

        /// \brief Trait for chosing the supremum metric by the template parameter 0
        template <>
        struct MetricTraitsCL<0> { typedef internal::SupMetricCL metric_type; };

        /// \brief Trait for chosing the one metric by the template parameter 1
        template <>
        struct MetricTraitsCL<1> { typedef internal::OneMetricCL metric_type; };

        /// \brief Trait for chosing the Euclidean metric by the template parameter 2
        template <>
        struct MetricTraitsCL<2> { typedef internal::EuclideanMetricCL metric_type; };



        /* ********************************************************************* */
        /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S      */
        /* ********************************************************************* */


        /* ******************************************************************** */
        /*  S U P  M E T R I C  C L A S S                                       */
        /* ******************************************************************** */

        /** Determine the distance of two points a and b by the supremum norm.
            \param a    first point
            \param b    second point
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \return \f$ \| a-b \|_\infty \f$
        */
        template <typename T, usint K>
        inline T SupMetricCL::distance( T const * a, T const * b)
        {
            T result= static_cast<T>(0);
            for ( usint i=0; i<K; ++i){
                result = std::max( static_cast<T>(std::abs(a[i]-b[i])), result);
            }
            return result;
        }


        /** Check if an axis-aligned bounding box intersects a sphere which is given in
            the supremum norm.
            \param C    center of the sphere
            \param r    radius of the sphere
            \param bb   the axis-aligned bounding box
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \return true iff the bounding box intersects the sphere
        */
        template <typename T, usint K>
        inline bool SupMetricCL::intersects( const T * const C, const T& r, const BoundingBoxCL<T,K>& bb)
        {
            // in each direction, there must be an intersection of the lines
            for ( usint i=0; i<K; ++i){
                if( C[i]-r < bb[2*i]) {
                    if ( C[i]+r < bb[2*i]){
                        return false;
                    }
                }
                else if( C[i]+r > bb[2*i+1] ){
                    if ( C[i]-r > bb[2*i+1]) {
                        return false;
                    }
                }
            }
            return true;
        }



        /* ******************************************************************** */
        /*  O N E  M E T R I C  C L A S S                                       */
        /* ******************************************************************** */

        /** Determine the distance of two points a and b by the one norm.
            \param a    first point
            \param b    second point
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \return \f$ \| a-b \|_1 \f$
        */
        template <typename T, usint K>
        inline T OneMetricCL::distance( T const * a, T const * b)
        {
            T result= static_cast<T>(0);
            for ( usint i=0; i<K; ++i){
                result += std::abs(a[i]-b[i]);
            }
            return result;
        }


        /** Check if an axis-aligned bounding box intersects a sphere which is given in
            the one norm.
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \todo Implement me
        */
        template <typename T, usint K>
        inline bool OneMetricCL::intersects( const T * const, const T&, const BoundingBoxCL<T,K>&)
        {
            return false;
        }



        /* ******************************************************************** */
        /*  E U C L I D E A N  M E T R I C  C L A S S                           */
        /* ******************************************************************** */

        /** Determine the distance of two points a and b by the Euclidean norm.
            \param a    first point
            \param b    second point
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \return \f$ \| a-b \|_2 \f$
        */
        template <typename T, usint K>
        inline T EuclideanMetricCL::distance( T const * a, T const * b)
        {
            T result= static_cast<T>(0);
            for ( usint i=0; i<K; ++i){
                result += std::pow( a[i]-b[i], 2);
            }
            return result;
        }


        /** Check if an axis-aligned bounding box intersects a sphere which is given in
            the Euclidean norm.
            \param C    center of the sphere
            \param r    radius of the sphere
            \param bb   the axis-aligned bounding box
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \return true iff the bounding box intersects the sphere
        */
        template <typename T, usint K>
        inline bool EuclideanMetricCL::intersects( const T * const C, const T& r, const BoundingBoxCL<T,K>& bb)
        {
            T dmin = static_cast<T>(0);
            for( int i=0; i<K; ++i) {
                if( C[i] < bb[2*i]) 
                    dmin += std::pow( C[i] - bb[2*i], 2); 
                else if( C[i] > bb[2*i+1] ) 
                    dmin += std::pow( C[i] - bb[2*i+1], 2);     
                }
            return (dmin <= r);
        }



        /* ******************************************************************** */
        /*  M A R K O W S K I  M E T R I C  C L A S S                           */
        /* ******************************************************************** */

        /** Determine the distance of two points a and b by a Markowski norm.
            \param a    first point
            \param b    second point
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \tparam metric Markowski metric
            \return \f$ \| a-b \|_p \f$ with \a p=metric
        */
        template <usint metric>
        template <typename T, usint K>
        inline T MarkowskiMetricCL<metric>::distance( T const * a, T const * b)
        {
            T result= static_cast<T>(0);
            const usint p= metric!=0 ? metric : std::numeric_limits<usint>::max();
            for ( usint i=0; i<K; ++i){
                result += static_cast<T>(std::pow( std::abs( a[i]-b[i]), p));
            }
            return std::pow( result, T(1.)/p);
        }


        /** Check if an axis-aligned bounding box intersects a sphere which is given in
            a Markowski norm.
            \tparam T   type of coordinates, e.g., double of float
            \tparam K   dimension
            \tparam metric Markowski metric
            \return true iff the bounding box intersects the sphere
            \todo Implement me
        */
        template <usint metric>
        template <typename T, usint K>
        inline bool MarkowskiMetricCL<metric>::intersects( const T * const, const T&, const BoundingBoxCL<T,K>&)
        {
            return false;
        }
    }    
}
}
#endif
