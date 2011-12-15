/// \file   tree.h
/// \brief  kd-tree header and implementation
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



#ifndef TREE_H
#define TREE_H

#include "misc/kd-tree/node.h"
#include "misc/kd-tree/bucket.h"
#include "misc/kd-tree/kd_tree_utils.h"

namespace DROPS{
namespace KDTree{

    /// \brief Storing a kd-tree. \anchor TreeCL
    /** This class only represents a kd-tree and is capable of performing search
        operations on the tree. To construct a tree, use the class tree_builder.
        The tree is given as a pseudo, generalized, bucket tree. The type of
        each coordinate is given by the template parameter T, the dimension by K
        and the size of the bucket by BucketSize (default=12).

        To obtain information about the storage and the height of the tree, call the
        function \ref write_info. This function additionally provide you with a detailed
        view of the tree if the height of the tree is small.

        \tparam T           type of coordinates, e.g., double of float
        \tparam K           dimension
        \tparam BucketSize  size of the bucket
    */
    template <typename T, usint K, int BucketSize=12>
    class TreeCL
    {
    public:   // ------- type definitions ------- 
        typedef internal::NodeCL<T,K,BucketSize> node_type;                 ///< type of a node
        typedef std::vector<T>                   value_vec_type;            ///< type of a vector storing all elements
        typedef std::vector<size_t>              index_vec_type;            ///< tpye of a vector storing indices

    private:  // ------- member variables ------- 
        node_type*     p_root;                                              ///< the root node
        value_vec_type p_data;                                              ///< storing the data in an optimized fashion
        index_vec_type p_origidx;                                           ///< original index of a permuted index

    private:  // ------- member functions ------- 
        /// \name recursive functions
        //@{
        /// \brief Determine recursively the minimal and maximal level of a node.
        void level_info_rec( size_t&, size_t&, size_t&, const size_t, const node_type * const) const;
        //@}

    public:  // ------- member functions ------- 
        TreeCL() : p_root(0), p_data() {}                                   ///< constructor
        ~TreeCL() { delete p_root; }                                        ///< delete the tree recursively

        /// \brief get pointer to the root
        inline node_type* const&  root() const { return p_root; }
        /// \brief get reference to the data
        inline const value_vec_type& data() const { return p_data; }

        /// \name setters for the builder
        //@{
        inline node_type*&     root() { return p_root; }                    ///< access to the root node
        inline value_vec_type& data() { return p_data; }                    ///< access to the data array
        inline index_vec_type& origidx() { return p_origidx; }              ///< access to the index array
        //@}

        /// \brief get the number of points stored by the tree
        inline size_t size() const { return data().size()/K; }
        /// \brief get the memory used for storing the tree
        inline size_t memory() const { return p_root->memory()+sizeof(T)*p_data.size()+ sizeof(size_t)*p_origidx.size()+ sizeof(value_vec_type); }
        /// \brief Address of a point
        T const * addr( size_t i) const { return &p_data[0]+K*i;}
        /// \brief Get the point of an index
        inline std::vector<T> get_point( const size_t& idx) const;
        /// \brief Get the original index of a permuted index
        inline size_t get_orig( const size_t& idx) const { return p_origidx[idx]; }

        /// \brief write information of the tree on a stream
        void level_info( size_t& min, size_t& max, size_t& num_nodes) const;
    };



    /* ********************************************************************* */
    /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S      */
    /* ********************************************************************* */
    
    /** Traverse the tree in depth-first to determine recursively the minimal and 
        maximal level of a node. Furthermore, get the number of nodes stored by 
        the tree
        \param[out] min         minimal level of a node
        \param[out] max         maximal level of a node
        \param[ou] num_nodes    number of nodes
        \param[in] my_level     level of the node, e.g., 0 of the root node
        \param[in] node         node to be traversed
    */
    template <typename T, usint K, int BucketSize>
    void TreeCL<T,K,BucketSize>::level_info_rec( size_t& min, size_t& max, size_t& num_nodes, const size_t my_level, const node_type * const node) const
    {        
        ++num_nodes;                                                // increase number of nodes
        if ( node->hasLeft()){                                      // traverse left
            level_info_rec( min, max, num_nodes, my_level+1, node->left());
        }
        if ( node->hasRight()){                                     // traverse right
            level_info_rec( min, max, num_nodes, my_level+1, node->right());
        }
        if ( node->isLeaf()){                                       // leaf node, update max and min
            if ( my_level<min)
                min=my_level;
            if ( my_level>max)
                max=my_level;
        }
    }


    /** Get the point that belongs to the index idx.
        \param idx  index of the point in the rearranged data
        \return point belonging to idx
    */
    template <typename T, usint K, int BucketSize>
    inline std::vector<T> TreeCL<T,K,BucketSize>::get_point( const size_t& idx) const
    {
        return std::vector<T>( addr(idx), addr(idx)+K);
    }


    /** Determine the minimal and maximal level of leaf nodes
        \param min  minimal level of a leaf node
        \param max  maximal level of a leaf node
    */
    template <typename T, usint K, int BucketSize>
    void TreeCL<T,K,BucketSize>::level_info( size_t& min, size_t& max, size_t& num_nodes) const
    {
        min= std::numeric_limits<size_t>::max();
        max= std::numeric_limits<size_t>::min();
        num_nodes= 0;
        level_info_rec( min, max, num_nodes, 0, root());
    }

}
}
#endif
