/*============================================================================*/
/*       1         2         3         4         5         6         7        */
/*3456789012345678901234567890123456789012345678901234567890123456789012345678*/
/*============================================================================*/
/*                                             .                              */
/*                                               RRRR WW  WW   WTTTTTTHH  HH  */
/*                                               RR RR WW WWW  W  TT  HH  HH  */
/*      Header   : VveKdTree                     RRRR   WWWWWWWW  TT  HHHHHH  */
/*                                               RR RR   WWW WWW  TT  HH  HH  */
/*      Module   : VistaVisExt                   RR  R    WW  WW  TT  HH  HH  */
/*                                                                            */
/*      Project  : ViSTA                           Rheinisch-Westfaelische    */
/*                                               Technische Hochschule Aachen */
/*      Purpose  :  ...                                                       */
/*                                                                            */
/*                                                 Copyright (c)  1998-2000   */
/*                                                 by  RWTH-Aachen, Germany   */
/*                                                 All rights reserved.       */
/*                                             .                              */
/*============================================================================*/
/*                                                                            */
/*      CLASS DEFINITIONS:                                                    */
/*                                                                            */
/*        - KDTreeCL                                                          */
/*        - KDTreeNodeCL                                                      */
/*                                                                            */
/*============================================================================*/

#ifndef DROPS_KDTREE_H
#define DROPS_KDTREE_H

/**
 * (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
 *
 * Licensed under the Academic Free License version 1.1 found in file LICENSE
 * with additional provisions in that same file.
 *
 * The code was modified by Marc Wolter in order to remove boost data structures,
 * adapt the code to the majority of ViSTA style guides and improve readability.
 *
 * The code was once again modified by Oliver Fortmeier in order to adapt the code
 * to the majority of DROPS style guides and augment classes with serialization.
 **/

/*============================================================================*/
/* INCLUDES                                                                   */
/*============================================================================*/

#include "misc/utils.h"
#ifdef _PAR
#  include "parallel/parallel.h"
#endif
#include <vector>
#include <algorithm>

namespace DROPS {
/*============================================================================*/
/* DEFINES                                                                    */
/*============================================================================*/
template<typename T>
struct IntervalST
{
    T fLower, fUpper;
};

/*============================================================================*/
/* FORWARD DECLARATIONS                                                       */
/*============================================================================*/
template<typename T> class KDTreeNodeCL;
template<typename T> class CSearchRecord;
template<typename T> class KDTreeResultVectorCL;
template<typename T> struct KDTreeResultST;
template<typename T> class KDTreeSerialNode;

/*============================================================================*/
/* CLASS DEFINITIONS                                                          */
/*============================================================================*/

/**
 * class KDTreeCL
 *
 * Implement a kd tree for fast searching of points in a fixed data base
 * in k-dimensional Euclidean space.
 *
 * The main data structure, one for each k-d tree, pointing
 * to a tree of an indeterminate number of "KDTreeNodeCL"s.
 *
 **/

template<typename T>
class KDTreeCL
{
  public:
#ifdef _PAR
    /// \brief Types used for sending a KD-Tree
    struct MPITypesCL
    {
        ProcCL::DatatypeT IntervalSTType;        ///< Type for sending IntervalST
        ProcCL::DatatypeT KDTreeSerialNodeType;  ///< Type for sending KDTreeSerialNode

        MPITypesCL(int dim) {
            CreateTypes(dim);
        }

        ~MPITypesCL() {
            FreeTypes();
        }

        void CreateTypes(int dim);
        void FreeTypes();
    };
#endif

  private:
    KDTreeNodeCL<T>* Root_; // the root pointer
    /// pointing either to the_data or an internal
    /// rearranged data as necessary
    const T* InternalData_;
    /// the index for the tree leaves.  Data in a leaf with bounds [l,u] are
    /// in  'the_data[ind[l],*] to the_data[ind[u],*]
    std::vector<int> vecTreeLeaves_;
    /// if rearrange is true then this is the rearranged data storage.
    T* RearrangedData_;

    static const int BucketSize_ = 12; // global constant.

    /**
     * "the_data" is a reference to the underlying multi_array of the
     * data to be included in the tree.
     *
     * NOTE: this structure does *NOT* own the storage underlying this.
     * Hence, it would be a very bad idea to change the underlying data
     * during use of the search facilities of this tree.
     * Also, the user must deallocate the memory underlying it.
     **/
    const T* Data_;

    const int NumberOfPoints_; // number of data points
    int Dim_; //
    bool SortResults_; // USERS set to 'true'.
    const bool Rearrange_; // are we rearranging?

    KDTreeNodeCL<T>* readNodes(const std::vector<KDTreeSerialNode<T> > &nodes, int idx, int dim, std::vector<IntervalST<T> > &vecBounds);
    void SetData(T& din);
    void BuildTree(); // builds the tree.  Used upon construction.
    KDTreeNodeCL<T>* BuildTreeForRange(int l, int u, KDTreeNodeCL<T>* parent);
    void SelectOnCoordinate(int c, int k, int l, int u);
    int SelectOnCoordinateValue(int c, T alpha, int l, int u);
    void SpreadInCoordinate(int c, int l, int u, IntervalST<T>& interv);

public:
    /**
     * Constructor, passing in a 2D-float array containing the data, with the number of points iN
     * and for each point iDim values.
     *
     * Constructor has optional 'dim_in' feature, to use only
     * first 'dim_in' components for definition of nearest neighbors.
     */
    KDTreeCL(T* data_in, int iN, int iDim, bool rearrange_in = true, int dim_in = -1);

    /** Constructor, for MPI recieved KDTree
     *
     * Constructor is similar to the one originally used
     */
    KDTreeCL(const std::vector<KDTreeSerialNode<T> > &nodes, std::vector<T> &rearranged_data, int dim, std::vector<IntervalST<T> > vecBounds);

    // destructor
    ~KDTreeCL();

public:
    /// \name search routines

    /// \brief search for n nearest to a given query vector 'qv' usin
    /// exhaustive slow search.  For debugging, usually.
    void GetNNearest_BruteForce(std::vector<T>& qv, int nn, KDTreeResultVectorCL<T>& result);

    /// \brief search for n nearest to a given query vector 'qv'.
    void GetNNearest(std::vector<T>& qv, int nn, KDTreeResultVectorCL<T>& result);

    /// \brief search for 'nn' nearest to point [idxin] of the input data, excluding
    /// neighbors within correltime
    void GetNNearest_AroundPoint(int idxin, int correltime, int nn, KDTreeResultVectorCL<T>& result);

    /// \brief  search for all neighbors in ball of size (square Euclidean distance)
    /// r2.   Return number of neighbors in 'result.size()',
    void GetNearestByRadius(std::vector<T>& qv, T r2, KDTreeResultVectorCL<T>& result);

    /// \brief like 'r_nearest', but around existing point, with decorrelation
    /// interval.
    void GetNearestByRadius_AroundPoint(int idxin, int correltime, T r2, KDTreeResultVectorCL<T>& result);

    /// \brief count number of neighbors within square distance r2.
    int GetNumberOfNeighborsWithin(std::vector<T>& qv, T r2);

    /// \brief like r_count, c
    int GetNumberOfNeighbors_AroundPoint(int idxin, int correltime, T r2);

    /// \brief serialze the tree
    void SerializeTree(std::vector<KDTreeSerialNode<T> > &SerialNodes, std::vector<T> &SerialData, std::vector<IntervalST<T> > &VecBounds);

#ifdef _PAR
    /// \brief Gather tree on all processors, i.e., allgather
    void SerializeAndGather(std::vector<KDTreeCL<T> >&);
#endif

    friend class KDTreeNodeCL<T> ;
    friend class CSearchRecord<T> ;
};

/**
 * class KDTreeNodeCL
 *
 * a node in the tree.  Many are created per tree dynamically..
 *
 **/
template<typename T>
class KDTreeNodeCL
{
private:
    /// visible to self and KDTreeCL.
    friend class KDTreeCL<T> ; ///< allow KDTreeCL to access private
    friend class MPITypesCL;
    friend class KDTreeSerialNode<T> ;

    int CutDimension_; ///< dimension to cut;
    T CutValue_, CutValueLeft_, CutValueRight_; ///< cut value
    int l, u; ///< extents in index array for searching

    std::vector<IntervalST<T> > vecBounds_; ///< [min,max] of the box enclosing all points

    KDTreeNodeCL<T> *LeftNode_, *RightNode_; ///< pointers to left and right nodes.

    /// recursive innermost core routine for searching..
    void search(CSearchRecord<T>& sr);

    /// return true if the bounding box for this node is within the
    /// search range given by the searchvector and maximum ballsize in 'sr'.
    bool box_in_search_range(CSearchRecord<T>& sr);

    /// for processing final buckets.
    void process_terminal_node(CSearchRecord<T>& sr);
    void process_terminal_node_fixedball(CSearchRecord<T>& sr);

    /// for transmitting through MPI
    void SerializeNodes(std::vector<KDTreeSerialNode<T> > &nodes, int dim, std::vector<IntervalST<T> > &VecBounds);

public:
    /// constructor
    KDTreeNodeCL(int dim);
    KDTreeNodeCL(const KDTreeSerialNode<T> &iNode, int dim);

    /// destructor
    ~KDTreeNodeCL();
};

/// \brief The search routines return a (wrapped) vector
/// of these.
template<typename T>
struct KDTreeResultST
{
public:
    T Distance; // its square Euclidean distance
    int Index; // which neighbor was found
};

/// \brief inherit a std::vector<KDTreeResultST>
/// but, optionally maintain it in heap form as a priority
/// queue.
template<typename T>
class KDTreeResultVectorCL: public std::vector<KDTreeResultST<T> >
{
public:

    /// \brief add one new element to the list of results, and
    /// keep it in heap order.  To keep it in ordinary, as inserted,
    /// order, then simply use push_back() as inherited
    /// via vector<>
    void push_element_and_heapify(KDTreeResultST<T>&);
    T replace_maxpri_elt_return_new_maxpri(KDTreeResultST<T>&);

    /// \brief return the distance which has the maximum value of all on list,
    /// assuming that ALL insertions were made by
    /// push_element_and_heapify()
    T max_value();
};

/// \brief Data structure for serialization of a KD-Tree
template<typename T>
struct KDTreeSerialNode
{
    /// \brief  position of left and right child in a serialize tree.
    /** -1 if left or right Child does not exist */
    int leftChild, rightChild;

    ///information from the KDTreeNodeCL
    int CutDimension_;                          ///< dimension to cut;
    T CutValue_, CutValueLeft_, CutValueRight_; ///< cut value
    int l, u;                                   ///< extents in index array for searching

    KDTreeSerialNode(const KDTreeNodeCL<T> &iNode);
    KDTreeSerialNode();
    ~KDTreeSerialNode() { }
};

} // end of namespace
#include "misc/KDtree.tpp"

#endif
