/// \file   tree_builder.h
/// \brief  Builder pattern of a kd-tree
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


#ifndef TREE_BUILDER_H
#define TREE_BUILDER_H

#include "misc/kd-tree/tree.h"
#include "misc/kd-tree/node.h"
#include "misc/kd-tree/bounding_box.h"
#include "misc/kd-tree/bucket.h"
#include "misc/kd-tree/kd_tree_utils.h"

namespace DROPS{
namespace KDTree{

    /// \brief Builder pattern to generate a kd-tree \anchor TreeBuilderCL
    /** This class can construct a kd-tree in \f$ O(n \log(n)) \f$. It 
        constructs a generalied, pseudo, bucket kd-tree. It uses the median
        in the largest spread to bipartite the space.

        This builder is also capable of detecting equal points during the 
        construction. If such points are found, the equal points are represented 
        in the kd-tree by a single point.

        \tparam T            type of a scalar data entry, e.g., component of a vector
        \tparam K            dimension of the space
        \tparam BucketSize   size of buckets
    */
    template <typename T, usint K, int BucketSize=12>
    class TreeBuilderCL
    {
    public:    // ------- type definitions ------- 
        typedef TreeCL<T, K, BucketSize>            tree_type;              ///< type of the tree
        typedef internal::NodeCL<T,K,BucketSize>    node_type;              ///< type of a node
        typedef internal::BoundingBoxCL<T,K>        bounding_box_type;      ///< type of the bounding box
        typedef internal::BucketCL<BucketSize>      bucket_type;            ///< type of the bucket


    private:   // ------- classes ------- 
        class  Pred_less;                                                   ///< compare values in a given dimension
        class  Pred_isLeft;                                                 ///< compare value in a dimension with a given split value
        struct BuildTask;                                                   ///< task of a node to be built for the interval [first, last) in the index vector

    private:   // ------- type definitions ------- 
        typedef std::vector<size_t>               ivector_type;             ///< type of an index vector
        typedef ivector_type::iterator            ivector_iterator;         ///< type of iterator of an index vector
        typedef std::queue<BuildTask>             build_queue_type;         ///< type of the queue that stores all nodes to be built
        

    private:  // ------- member variables ------- 
        tree_type&          p_tree;                                         ///< the tree to be built
        T const *           p_data;                                         ///< points in the tree
        ivector_type        p_idxvec;                                       ///< index vector of all elements
        size_t              p_skipped;                                      ///< number of points merged into a single point (because they are indentical)


    private:  // ------- member functions ------- 
        /// \name Helper functions for constructing a tree
        //@{
        void determineBB( const size_t&, const size_t&, bounding_box_type&);///< determine the bounding box
        void estimateBB( const node_type*, const bool, bounding_box_type&); ///< estimate the bounding box
        void uniteBB(node_type* node);                                      ///< determine the boounding box
        T sort( const size_t&, const size_t&, const usint, size_t&);        ///< partition the index vector and return the median
        void buildLeaf( const size_t&, const size_t&, node_type*&);         ///< build the bucket of a leaf node and determine the exact bounding box
        void buildNode( BuildTask&, build_queue_type&);                     ///< build a single node
        void buildRoot( std::vector<build_queue_type>&, const size_t);      ///< build the "widened" root node
        void optimize( size_t&, node_type*&);                               ///< optimize the storage of the points
        void clear() { p_idxvec.clear(); }                                  ///< free the memory
        //@}

        /// \brief dummy declaration to get rid of warnings
        TreeBuilderCL<T,K,BucketSize>& operator= ( const TreeBuilderCL<T,K,BucketSize>&);


    public:  // ------- member functions ------- 
        TreeBuilderCL( tree_type& tree)                                     ///< Constructor
            : p_tree( tree), p_data(0), p_skipped(0) {}
        ~TreeBuilderCL() {}                                                 ///< Destructor
        /// \name build the tree
        //@{
        void build( T const *, const size_t);                               ///< Calling by a c-array
        void build( const std::vector<T>&);                                 ///< Calling by a std::vector
        void build( const std::valarray<T>&);                               ///< Calling by a std::valarray
        //@}        
    };



    /* ******************************************************************** */
    /*  D E F I N I T I O N   O F   I N T E R N A L   C L A S S E S         */
    /* ******************************************************************** */

    /// \brief Compare values in a given dimension
    /** This subclass is used to sort ranges of the index vector.*/
    template <typename T, usint K, int BucketSize>
    class TreeBuilderCL<T,K,BucketSize>::Pred_less 
    {
    protected:
        const usint p_sdim;                                             ///< dimension
        T const *   p_data;                                             ///< access to data
        Pred_less& operator=(const Pred_less&);                         ///< dummy declaration to get rid of warnings

    public:
        /// \brief Constructor
        /** This class compares two points in a given dimension 
            \param dim  the dimension in which the points are compared
            \param d    the data array
        */
        Pred_less( const usint dim, T const * d)
            : p_sdim(dim), p_data(d) {}
        /// \brief Compare two points
        /** The points are located at i and j.
            \param i    the position of the first point
            \param j    the position of the second point
            \return true iff the coordinate in \a dim dimension of point i is smaller than of point j
        */
        inline bool operator() (const size_t i, const size_t j) const
            { return p_data[K*i+p_sdim]<p_data[K*j+p_sdim]; }
    };


    /// \brief Compare value in a dimension with a given split value
    /** Derived from the class Pred_less because it makes use of common attributes. */
    template <typename T, usint K, int BucketSize>
    class TreeBuilderCL<T,K,BucketSize>::Pred_isLeft : public TreeBuilderCL<T,K,BucketSize>::Pred_less 
    {
    private:
        const T& p_val;                                                 ///< value for comparing
        Pred_isLeft& operator=(const Pred_isLeft&);                     ///< dummy declaration to get rid of warnings
    public:
        /// \brief Constructor
        /** This class checks if a point is located to the left w.r.t. to a given dimension and a split value
            \param dim  the dimension in which the points are compared
            \param val  split value
            \param d    the data array
        */
        Pred_isLeft( const usint dim, const T& val, T const * d)
            : Pred_less(dim, d), p_val(val) {}
        /// \brief Point i is located to the left of the split value
        /** Check if the value of the point i in dimension \a dim is left of the split value given 
            by \a val.
            \param i    index of the point
            \return true iff the coordinate of the point i in dim dimension is smaller that the split value
        */
        inline bool operator() (const size_t i) const
            { return Pred_less::p_data[K*i+Pred_less::p_sdim]<p_val; }
    };


    /// \brief Task of a node to be built for the interval [first, last) in the index vector
    /** Store all nodes that will span a subtree in a queue. This task specifies the node. */
    template <typename T, usint K, int BucketSize>
    struct TreeBuilderCL<T,K,BucketSize>::BuildTask {
        const size_t     first;                                 ///< first element in the index vector
        const size_t     last;                                  ///< last element in the index vector
        node_type*&      node;                                  ///< pointer to the node to be constructed
        const node_type* parent_node;                           ///< pointer to the parent bucket
        const bool       left;                                  ///< create a left or right node
        /// \brief Constructor
        /** The Build task is given by first and last index, as well as the node to be constructed.
            \param f        first index
            \param l        last index
            \param n        the node
            \param parnode  parent node
            \param le       left or right node to be built
        */
        BuildTask( size_t f, size_t l, node_type*& n, const node_type* parnode, const bool le) 
            : first(f), last(l), node(n), parent_node(parnode), left(le) {}
        BuildTask& operator=(const BuildTask&);
    };



    /* ******************************************************************** */
    /*  D E F I N I T I O N   O F   T E M P L A T E   F U N C T I O N S     */
    /* ******************************************************************** */

    /** The bounding box is a rectange determining the largest and 
        smallest value in each dimension in the data addressed by
        the index vector in the range [first,last)

        If \f$ n \f$ is the length of the interval [first,last), this function
        has a time complexity of \f$ O(n*K) \f$.

        \param first    first index in the index vector
        \param last     behind last index in the index vector
        \param bb       output, i.e., the bounding box
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::determineBB( const size_t& first, const size_t& last, bounding_box_type&bb)
    {
        // fill the bounding box with first element
        for ( int i=0; i<K; ++i){
            bb[ 2*i]  = p_data[ p_idxvec[first]*K+i];
            bb[ 2*i+1]= p_data[ p_idxvec[first]*K+i];
        }

        // determine the largest and smallest value in each dimension
        for ( size_t idx=first+1; idx<last; ++idx){
            for ( int i=0; i<K; ++i){
                bb[ 2*i]  = std::min( bb[ 2*i],   p_data[ p_idxvec[idx]*K+i]);
                bb[ 2*i+1]= std::max( bb[ 2*i+1], p_data[ p_idxvec[idx]*K+i]);
            }
        }
    }


    /** Estimate the bounding box of a node by using the bounding box of the parent 
        node and the median. If this function is called for the root node, the exact
        bounding box is determined.
        \param parnode  parent node
        \param left     true iff this is the left child of the parent node
        \param bb       bounding box of this node (i.e., the bounding box to be estimated)
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::estimateBB( const node_type* parnode, const bool left, bounding_box_type& bb)
    {
        if ( parnode!=0) {  // not the root node
            const usint sdim              = parnode->sdim();
            const T&    sval              = parnode->sval();
            const bounding_box_type& parbb= parnode->bounding_box();
            // Copy the bounding box of the parent and cut it by the median
            bb= parbb;
            bb[2*sdim + (left ? 1 : 0)]= sval;
        }
        else {              // this only happens for the root node
            determineBB( 0, p_idxvec.size(), bb);
        }
    }


    /** Build a bounding box of a node by uniting the bounding box of the children. 
        That is the bounding box includes the bounding boxes of both children.
        \param node     a node that has either two children or no child.
        \pre each node has either both children or no child!
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::uniteBB(node_type* node)
    {
        if ( node->hasChildren()){
            bounding_box_type&       bb      = node->bounding_box();
            const bounding_box_type& bb_left = node->left()->bounding_box();
            const bounding_box_type& bb_right= node->right()->bounding_box();

            for ( usint i=0; i<K; ++i){
                bb[2*i  ]= std::min( bb_left[2*i],   bb_right[2*i]);
                bb[2*i+1]= std::max( bb_left[2*i+1], bb_right[2*i+1]);
            }
        }
    }


    /** Split the indices in the index vector in the range [first,last), such that
        the indices [first,start_right) addresses points which have a smaller value
        in the dimension sdim as the median. Thus, the elements [start_right,last)
        points to data which are larger than median in the dimension sdim.

        On output, the parameter start_right contains the split index and the median
        is returned.

        The complexity of this function has on average a complexity of \f$ O(n) \f$
        where \f$ n \f$ is the length of the interval [first,last).
        
        \param[in]  first       first index in the index vector
        \param[in]  last        behind last index in the index vector
        \param[in]  sdim        split dimension
        \param[out] start_right index to first element which a larger value than the median
        \return median in the dimension sdim
    */
    template <typename T, usint K, int BucketSize>
    T TreeBuilderCL<T,K,BucketSize>::sort( const size_t& first, const size_t& last, const usint sdim, size_t& start_right)
    {
        //
        // determine the median as the element in the middle in the interval [first,last)
        //
        const size_t mid= (last+first)/2;           // index of the middle
        Pred_less less( sdim, p_data);              // comparing elements in the dimension sdim
        // get the element in the middle ...
        std::nth_element( p_idxvec.begin()+first, p_idxvec.begin()+mid, p_idxvec.begin()+last, less);
        // ... which is the median
        const T median= p_data[ (*(p_idxvec.begin()+mid)*K)+sdim];

        //
        // now, get the index of the first element that is larger or equal to the median.
        //
        Pred_isLeft isleft( sdim, median, p_data);  // check if an element is smaller than the median
        // use std::partition to partition the elements
        ivector_iterator first_right= std::partition( p_idxvec.begin()+first, p_idxvec.begin()+last, isleft);
        // pointer arithmetics leads to the first element in the right part
        ivector_iterator first_left= p_idxvec.begin();
        start_right= static_cast<size_t>( first_right-first_left);
        //start_right= static_cast<size_t>( std::distance( first_left, first_right));   // does not work for Oracle studio compiler

        // return the median
        return median;
    }


    /** If a leaf node is found, i.e. the node represents points less than BucketSize elements,
        the bucket is constructed. The indices are located in the range [first,last). Building 
        the bucket has a time complexity of \f$ O(BucketSize) \f$. Furthermore, the exact bounding
        box of this node is determined and the split dimension is set to K (i.e., made invalid).

        \param[in]  first       first index in the index vector
        \param[in]  last        behind last index in the index vector
        \param[in]  node        the node that contains the bucket
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::buildLeaf( const size_t& first, const size_t& last, node_type*& node)
    {
        // Build the bucket
        const int num_elems= last-first;
        for ( int i=0; i<BucketSize; ++i){
            (node->bucket())[i]= i<num_elems ? p_idxvec[first+i] : NoIdx;
        }
        // Build the bounding box
        determineBB( first, last, node->bounding_box());
        // Mark node as a leaf
        node->sdim()=K;
    }


    /** This is the main function to build a kd-tree. Therefore, the
        node \a node is generated that represents the points addressed by the
        indices in [first,last). If the range contains more elements than
        \a BucketSize, than an internal node is constructed, otherwise,
        a leaf node is constructed.

        In case of an internal node, the indices in the interval
        [task.first,task.last) are divided by the median of the largest 
        spread.

        This function also handles identical points by combining them in a 
        singe point.

        \param[in]  task    this task specifies the node that should be built
        \param      queue   the queue where to put the next tasks
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::buildNode( BuildTask& task, build_queue_type& queue)
    {
        // abbrevations
        const size_t first      = task.first;
        const size_t last       = task.last;
        node_type*&  node       = task.node;
        const node_type* parnode= task.parent_node;
        const bool side         = task.left;

        // allocate memory for the node
        node= new node_type();
   
        // build an internal node or a leaf
        if ( last-first>BucketSize) {                               // build an internal node
            // estimate the bounding box
            estimateBB( parnode, side, node->bounding_box());
            // first, try split dimension with the largest spread
            node->sdim()= node->bounding_box().suggestSplitDim();

            if ( node->sdim()==K) {     // this becomes a bucket node, since the largest spread is < eps
                // merge the points in the interval [first,last)
                buildLeaf( first, first+1, node);
                p_skipped += last-first-1;
            }
            else {                      // the bounding box has a volume>eps^K
                size_t start_right=last;
                bool success= false;
                for ( usint d=0; d<K && !success; ++d){
                    node->sval()= sort( first, last, node->sdim(), start_right);
                    if (start_right!=first && start_right!=last){       // we can split the interval
                        // this is the normal case!
                        queue.push( BuildTask( first,       start_right, node->left(),  node, true));
                        queue.push( BuildTask( start_right, last,        node->right(), node, false));
                        success= true;
                    }
                    else {                                              // try the next dimension
                        node->sdim() = (d+node->sdim()) % K;
                    }
                }
                if ( !success){             // merge the points in the interval [first,last)
                    determineBB( first, last, node->bounding_box());
                    node->sdim()= node->bounding_box().suggestSplitDim();
                    if ( node->sdim()!=K){  // warn if something strange has happened
                        std::cerr << "Cannot split the interval " << node->bounding_box() << " with points \n";
                        for ( size_t i=first; i<last; ++i){
                            std::cerr << "- (";
                            for ( usint j=0; j<K; ++j) std::cerr << p_data[ p_idxvec[i]*K+j] << ' ';
                            std::cerr << ")\n";
                        }
                        std::cerr << "Giving up and combine all points in a single node :-(" << std::endl;
                    }
                    buildLeaf( first, first+1, node);
                    p_skipped += last-first-1;
                }
            }
        }
        else {                                                      // build a leaf node
            buildLeaf( first, last, node);
        }
    }


    /** The "widened" root of the tree is constructed in the following approach. The root node 
        (on level 0) is constructed by a single thread. The nodes on level 1 are constructed by
        two threads, etc. That is, the nodes on level i are constructed by 2^i threads.

        \param queues   the queue in queues[i] is assigned to the thread i
        \param n        number of data points
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::buildRoot( std::vector<build_queue_type>& queues, const size_t n)
    {
        // if only a single thread is used, only the root is pushed into the queue
        if ( get_num_threads()==1){
            queues[0].push( BuildTask( 0, n, p_tree.root(), 0, true));
            return;
        }

        // use maximal num_threads/2 threads to generate the sub roots of the first log2(num_threads) levels.
        for ( int num_thr=1; num_thr<=get_num_threads()/2; num_thr*=2){
#pragma omp parallel for num_threads(num_thr)
            for ( int thread=0; thread<num_thr; ++thread){
                BuildTask* sub_root=0;                                  // this node should be constructed
                if ( num_thr==1){                                       // 0-th level, this is the root of the tree
                    sub_root= new BuildTask( 0, n, p_tree.root(), 0, true);      
                }
                else{                                                   // it is constructed by the previous level
                    if ( !queues[thread].empty()){
                        sub_root= new BuildTask( queues[thread].front()); 
                        queues[thread].pop();
                    }
                }
                if ( sub_root!=0){                                      // for very small trees it may happen, that the tree is already constructed
                    buildNode( *sub_root, queues[thread]);
                    if ( !queues[thread].empty()){
                        BuildTask left(  queues[thread].front()); queues[thread].pop();
                        queues[thread*2+0].push( left);
                    }
                    if ( !queues[thread].empty()){
                        BuildTask right( queues[thread].front()); queues[thread].pop();
                        queues[thread*2+1].push( right);
                    }
                    delete sub_root;
                }
            }
        }
    }


    /** The memory layout to store all the points is optimized in the following 
        way. All points addressed by a single bucket are aligned consecutively
        in the memory. 
        
        This causes an additionally computational work of \f$ O(n) \f$ while
        constructing the tree.

        Furthermore, the bounding boxes of all internal nodes are determined by
        uniting the bounding boxes of their children.

        \param[in,out] first_free   first index in the vector p_tree.data() that has not been assigned so far.
        \param[in]     node         the node that should be optimized
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::optimize( size_t& first_free, node_type*& node)
    {
        // go down in the tree until a leaf node is found
        if ( node->hasLeft()){
            optimize( first_free, node->left());
        }
        if ( node->hasRight()){
            optimize( first_free, node->right());
        }
        // a leaf stores a bucket. So optimize the data for that bucket.
        if ( node->isLeaf()){
            bucket_type& bucket= node->bucket();
            typename tree_type::value_vec_type& tree_data= p_tree.data();
            typename tree_type::index_vec_type& origidx  = p_tree.origidx();
            for ( int i=0; i<BucketSize && bucket[i]!=NoIdx; ++i){
                for ( usint j=0; j<K; ++j){
                    tree_data[first_free*K+j]= p_data[bucket[i]*K+j];                    
                }
                origidx[first_free]= bucket[i];
                bucket[i]= first_free++;
            }
        }
        // unite the bounding boxes of its children
        if ( node->hasChildren()){
            uniteBB( node); // do not do that for leaf nodes
        }
    }


    /** Build a kd-tree that covers the \a n points located in the field. The field
        data is not changed by this routine. However, the kd-tree does not
        actually store the data but references to the data. Therefore, it is necessary
        that the data is not changed by the user before all search operations
        are performed. The field data contains all points. The K-dimensional 
        point i is located at K*i,...,K*i+(K-1). Thus, the length of data should
        equal n*K.
        \param[in] data     data containing all points 
        \param[in] n        number of points in data
    */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::build( T const * data, const size_t n)
    {
        p_data= data; p_skipped= 0;

        // construct the index vector
        p_idxvec.resize( n);
        for ( size_t i=0; i<n; ++i) 
            p_idxvec[i]= i;

        // for each thread a single queue is used
        std::vector<build_queue_type> queues( get_num_threads());

        // build the hat (= "widend" root) of the tree
        buildRoot( queues, n);

#pragma omp parallel for
        for ( int thread=0; thread<get_num_threads(); ++thread){
            while ( !queues[thread].empty()){
                BuildTask task( queues[thread].front());
                queues[thread].pop();
                buildNode( task, queues[thread]);
            }
        }

        // optimize the data
        p_tree.data().resize( (n-p_skipped)*K);
        p_tree.origidx().resize( (n-p_skipped)*K);
        size_t first_free=0;
        optimize( first_free, p_tree.root());
        
        // free the memory
        clear();
    }


    /** For a detailed description, see build by a c-array. */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::build( const std::vector<T>& data)
    {
        build( &data[0], data.size()/K);
    }


    /** For a detailed description, see build by a c-array. */
    template <typename T, usint K, int BucketSize>
    void TreeBuilderCL<T,K,BucketSize>::build( const std::valarray<T>& data)
    {
        build( &data[0], data.size()/K);
    }
}
}
#endif
