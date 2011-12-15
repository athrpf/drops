/// \file   node.h
/// \brief  Node of a kd-tree
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


#ifndef NODE_H
#define NODE_H

#include "misc/kd-tree/bounding_box.h"
#include "misc/kd-tree/bucket.h"
#include "misc/kd-tree/kd_tree_utils.h"

namespace DROPS{
namespace KDTree{

    namespace internal{

        /// \brief A node of the kd-tree.
        /** A node of a kd-tree is part of the binary search tree. If the node is
            a leaf, then a bucket is stored. Otherwise, this node only points to a left
            and a right node.
            \tparam T            type of values, e.g., double, float, ...
            \tparam K            dimension of a point, e.g., 2 or 3
            \tparam BucketSize   size of the buckets
        */
        template <typename T, usint K, int BucketSize>
        class NodeCL
        {
        public:   // ------- type definitions ------- 
            typedef NodeCL<T,K,BucketSize>  node_type;                      ///< type of a node
            typedef BoundingBoxCL<T,K>      bounding_box_type;              ///< type of the bounding box
            typedef BucketCL<BucketSize>    bucket_type;                    ///< type of the bucket

        private:  // ------- member variables ------- 
            node_type*        p_left;                                       ///< left child
            node_type*        p_right;                                      ///< right child
            bounding_box_type p_bbox;                                       ///< bounding box
            usint             p_sdim;                                       ///< split dimension
            T                 p_sval;                                       ///< split value
            bucket_type       p_bucket;                                     ///< bucket of points

        public:  // ------- member functions ------- 
            /// \brief Construct a node of a kd tree.
            NodeCL( const usint sdim=K)
                : p_left(0), p_right(0), p_sdim(sdim), p_sval()
            {}

            /// \brief Destructor recursively deletes its children and also deletes the bucket.
            ~NodeCL();

            /// \brief check if the node is a leaf (a node is leaf, it it has no left (and no right) child)
            inline bool isLeaf() const { return p_left==0; }

            /// \name getters
            //@{
            /// \name check for children
            //@{
            inline bool hasLeft() const { return p_left!=0; }
            inline bool hasRight() const { return p_right!=0; }
            inline bool hasChildren() const { return hasLeft() || hasRight(); }
            //@}

            inline const usint&             sdim()         const { return p_sdim; }     ///< get split dimension
            inline const T&                 sval()         const { return p_sval; }     ///< get split value
            inline const bucket_type&       bucket()       const { return p_bucket; }   ///< get bucket
            inline const bounding_box_type& bounding_box() const { return p_bbox; }     ///< get bounding box
            inline node_type* const&        left()         const { return p_left; }     ///< get pointer to left child
            inline node_type* const&        right()        const { return p_right; }    ///< get pointer to right child
            //@}

            /// \name setters for the builder
            //@{
            inline usint&             sdim()         { return p_sdim; }     ///< set the split dimension
            inline T&                 sval()         { return p_sval; }     ///< set the split value
            inline bucket_type&       bucket()       { return p_bucket; }   ///< set the bucket
            inline bounding_box_type& bounding_box() { return p_bbox; }     ///< set the bounding box
            inline node_type*&        left()         { return p_left; }     ///< set left child
            inline node_type*&        right()        { return p_right; }    ///< set right child
            //@}

            /// \brief Get the memory used to store this node.
            inline size_t memory() const;
        };



        /* ********************************************************************* */
        /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S      */
        /* ********************************************************************* */

        /** Delete the memory for left and right child as well as for the bucket.*/
        template <typename T, usint K, int BucketSize>
        inline NodeCL<T,K,BucketSize>::~NodeCL()
        {
            if ( hasLeft())  
                delete p_left;
            if ( hasRight()) 
                delete p_right;
        }


        /** Just adding up the memory that is allocated by this node, and 
            that it points to.
            \return memory in Bytes
        */
        template <typename T, usint K, int BucketSize>
        inline size_t NodeCL<T,K,BucketSize>::memory() const
        {
            return 2*sizeof(node_type*)                         // pointer to left and right child
                + bounding_box_type::memory()                   // size of storing the bounding box
                + sizeof(usint)                                 // split dimension
                + sizeof(T)                                     // split value
                + bucket_type::memory()                         // size of the bucket
                + (hasLeft() ? p_left->memory() : 0)            // size of left node
                + (hasRight() ? p_right->memory() : 0);         // size of right node
        }
    }
}
}
#endif
