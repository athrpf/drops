/// \file   bounding_box.h
/// \brief  Axis-aligned bounding box of a node of a kd-tree
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


#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "misc/kd-tree/kd_tree_utils.h"

namespace DROPS{
namespace KDTree{

    namespace internal{

        /// \brief Axis-aligned bounding box of a node
        /** Each node contains a bounding box that encloses all points
            located in the subtree spanned by the node. This box is given
            by 2*K elements of type T.
            \tparam T   type of cooprinates, e.g., double of float
            \tparam K   dimension
        */
        template <typename T, usint K>
        class BoundingBoxCL
        {
        private:  // ------- member variables ------- 
            T p_bound[2*K];                                         ///< borders of the bounding box

        public:  // ------- member functions ------- 
            /// \brief empty bounding box (no values are assigned)
            BoundingBoxCL() { }

            /// \brief default destructor
            ~BoundingBoxCL() {}

            /// \brief Assignement operator
            inline BoundingBoxCL& operator=( const BoundingBoxCL&);

            /// \brief suggest a split dimension with the largest spread
            usint suggestSplitDim() const; 

            /// \brief accessing elements
            inline const T& operator[] (int i) const { return p_bound[i]; }

            /// \brief Assigning values
            inline T& operator[] (int i) { return p_bound[i]; }
            
            /// \brief memory for storing a bounding box
            static size_t memory() { return sizeof(T)*2*K; }
        };



        /* ********************************************************************* */
        /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S      */
        /* ********************************************************************* */

        template <typename T, usint K>
        inline BoundingBoxCL<T,K>& BoundingBoxCL<T,K>::operator=( const BoundingBoxCL<T,K>& bb)
        {
            std::copy( bb.p_bound, bb.p_bound+2*K, p_bound);
            return *this;
        }

        /** In common cases, it is advantageous to split the space in the dimension,
            where the largest spread occurs, i.e., the dimension i for which holds
            \f$ |x[i]-y[i]| \geq |x[j]-y[j]| \;\forall\,j \neq i \f$ where x and y
            are the left and right side of the bounding box. If the maximal spread 
            is smaller than the constant eps, than as suggested split dimension K is
            returned.
        */
        template <typename T, usint K>
        inline usint BoundingBoxCL<T,K>::suggestSplitDim() const 
        {
            usint sdim  = 0;
            T     spread= std::abs( p_bound[1]-p_bound[0]);
            for ( usint i=1; i<K; ++i){
                const T spread_dim= std::abs( p_bound[2*i+1]-p_bound[2*i]);
                if ( spread_dim>spread){
                    spread= spread_dim;
                    sdim  = i;
                }
            }
            return spread>std::numeric_limits<T>::epsilon() ? sdim : K;
        }
    }

    /// \brief write a bounding box on a stream
    template <typename T, usint K>
    inline std::ostream& operator << ( std::ostream& os, const internal::BoundingBoxCL<T,K>& bb)
    {
        os << '[';
        for (int i=0; i<2*K; i+=2)
            os << bb[i] << ' ';
        os << "] x [";
        for (int i=1; i<2*K; i+=2)
            os << bb[i] << ' ';
        os << ']';
        return os;
    }
}
}
#endif
