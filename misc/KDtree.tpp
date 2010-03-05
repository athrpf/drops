//
// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.


#include <algorithm>
#include <iostream>
#include <stdio.h>

namespace DROPS {

// utility

template<typename T>
inline T squared(const T x)
{
    return (x * x);
}

//
//       KDTREE2_RESULT implementation
//
template<typename T>
inline bool operator<(const KDTreeResultST<T>& e1, const KDTreeResultST<T>& e2)
{
    return (e1.Distance < e2.Distance);
}

//
//       KDTREE2_RESULT_VECTOR implementation
//
template<typename T>
T KDTreeResultVectorCL<T>::max_value()
{
    return (this->begin()->Distance); // very first element
}

template<typename T>
void KDTreeResultVectorCL<T>::push_element_and_heapify(KDTreeResultST<T>& e)
{
    push_back(e); // what a vector does.
    std::push_heap(this->begin(), this->end()); // and now heapify it, with the new elt.
}

template<typename T>
T KDTreeResultVectorCL<T>::replace_maxpri_elt_return_new_maxpri(KDTreeResultST<T>& e)
{
    // remove the maximum priority element on the queue and replace it
    // with 'e', and return its priority.
    //
    // here, it means replacing the first element [0] with e, and re heapifying.

    std::pop_heap(this->begin(), this->end());
    this->pop_back();
    this->push_back(e); // insert new
    std::push_heap(this->begin(), this->end()); // and heapify.
    return ((*this)[0].Distance);
}

//
//        KDTREE2 implementation
//

// constructor
template<typename T>
KDTreeCL<T>::KDTreeCL(T* data_in, int iN, int iDim, bool rearrange_in, int dim_in) :
    Root_(NULL),
    InternalData_(NULL),
    vecTreeLeaves_(iN), RearrangedData_(NULL),
    Data_(data_in), NumberOfPoints_(iN), Dim_(iDim), SortResults_(false), Rearrange_(rearrange_in)
{
    //
    // initialize the constant references using this unusual C++
    // feature.
    //
    if (dim_in > 0)
        Dim_ = dim_in;

    BuildTree();

    if (Rearrange_) {
        // if we have a rearranged tree.
        // allocate the memory for it.
        // printf("rearranging\n");
        delete[] RearrangedData_;
        RearrangedData_ = new T[NumberOfPoints_ * Dim_];

        // permute the data for it.
        for (int i = 0; i < NumberOfPoints_; i++) {
            for (int j = 0; j < Dim_; j++) {
                RearrangedData_[i * Dim_ + j] = Data_[vecTreeLeaves_[i] * Dim_ + j];
            }
        }
        InternalData_ = RearrangedData_;
    }
    else {
        InternalData_ = Data_;
    }
}

// function for recursively reading nodes from a serialized Tree
template<typename T>
KDTreeNodeCL<T>* KDTreeCL<T>::readNodes(const std::vector<KDTreeSerialNode<T> > &nodes, int idx, int dim, std::vector<IntervalST<T> > &vecBounds)
{
    KDTreeNodeCL<T>* newNode = new KDTreeNodeCL<T> (nodes[idx], dim);
    for (int i = 0; i < dim; i++)
        newNode->vecBounds_.push_back(vecBounds[idx * dim + i]);
    if (nodes[idx].leftChild != -1)
        newNode->LeftNode_ = readNodes(nodes, nodes[idx].leftChild, dim, vecBounds);
    if (nodes[idx].rightChild != -1)
        newNode->RightNode_ = readNodes(nodes, nodes[idx].rightChild, dim, vecBounds);
    return newNode;
}

// second constructor for a serialized Tree
template<typename T>
KDTreeCL<T>::KDTreeCL(const std::vector<KDTreeSerialNode<T> > &nodes, std::vector<T> &rearranged_data, int dim, std::vector<IntervalST<T> > vecBounds) :
    Root_(NULL),
    InternalData_(NULL),
    vecTreeLeaves_(rearranged_data.size() / dim), RearrangedData_(NULL),
    Data_(NULL), NumberOfPoints_(rearranged_data.size() / dim), Dim_(dim), SortResults_(false), Rearrange_(true)
{
    delete[] RearrangedData_;
    RearrangedData_ = new T[rearranged_data.size()];
    for (uint i = 0; i < rearranged_data.size(); ++i) {
        RearrangedData_[i] = rearranged_data[i];
    }
    InternalData_ = RearrangedData_;
    Root_ = readNodes(nodes, 0, dim, vecBounds);
}

// destructor
template<typename T>
KDTreeCL<T>::~KDTreeCL()
{
    delete Root_;
    delete[] RearrangedData_;
}

// building routines
template<typename T>
void KDTreeCL<T>::BuildTree()
{
    for (int i = 0; i < NumberOfPoints_; i++)
        vecTreeLeaves_[i] = i;
    Root_ = BuildTreeForRange(0, NumberOfPoints_ - 1, NULL);
}

template<typename T>
KDTreeNodeCL<T>* KDTreeCL<T>::BuildTreeForRange(int l, int u, KDTreeNodeCL<T>* parent)
{
    // recursive function to build
    KDTreeNodeCL<T>* node = new KDTreeNodeCL<T> (Dim_);
    // the newly created node.

    if (u < l) {
        return (NULL); // no data in this node.
    }

    if ((u - l) <= BucketSize_) {
        // create a terminal node.

        // always compute true bounding box for terminal node.
        for (int i = 0; i < Dim_; i++) {
            SpreadInCoordinate(i, l, u, node->vecBounds_[i]);
        }

        node->CutDimension_ = 0;
        node->CutValue_ = 0.0;
        node->l = l;
        node->u = u;
        node->LeftNode_ = node->RightNode_ = NULL;

    }
    else {
        //
        // Compute an APPROXIMATE bounding box for this node.
        // if parent == NULL, then this is the root node, and
        // we compute for all dimensions.
        // Otherwise, we copy the bounding box from the parent for
        // all coordinates except for the parent's cut dimension.
        // That, we recompute ourself.
        //
        int c = -1;
        T maxspread = 0.0;
        int m;

        for (int i = 0; i < Dim_; i++) {
            if ((parent == NULL) || (parent->CutDimension_ == i)) {
                SpreadInCoordinate(i, l, u, node->vecBounds_[i]);
            }
            else {
                node->vecBounds_[i] = parent->vecBounds_[i];
            }
            T spread = node->vecBounds_[i].fUpper - node->vecBounds_[i].fLower;
            if (spread > maxspread) {
                maxspread = spread;
                c = i;
            }
        }

        //
        // now, c is the identity of which coordinate has the greatest spread
        //

        if (false) {
            m = (l + u) / 2;
            SelectOnCoordinate(c, m, l, u);
        }
        else {
            T sum;
            T average;

            if (true) {
                sum = 0.0;
                for (int k = l; k <= u; k++) {
                    sum += Data_[vecTreeLeaves_[k] * Dim_ + c];
                }
                average = sum / static_cast<T> (u - l + 1);
            }
            else {
                // average of top and bottom nodes.
                average = (node->vecBounds_[c].fUpper + node->vecBounds_[c].fLower) * 0.5f;
            }

            m = SelectOnCoordinateValue(c, average, l, u);
        }

        // move the indices around to cut on dim 'c'.
        node->CutDimension_ = c;
        node->l = l;
        node->u = u;

        node->LeftNode_ = BuildTreeForRange(l, m, node);
        node->RightNode_ = BuildTreeForRange(m + 1, u, node);

        if (node->RightNode_ == NULL) {
            for (int i = 0; i < Dim_; i++)
                node->vecBounds_[i] = node->LeftNode_->vecBounds_[i];
            node->CutValue_ = node->LeftNode_->vecBounds_[c].fUpper;
            node->CutValueLeft_ = node->CutValueRight_ = node->CutValue_;
        }
        else if (node->LeftNode_ == NULL) {
            for (int i = 0; i < Dim_; i++)
                node->vecBounds_[i] = node->RightNode_->vecBounds_[i];
            node->CutValue_ = node->RightNode_->vecBounds_[c].fUpper;
            node->CutValueLeft_ = node->CutValueRight_ = node->CutValue_;
        }
        else {
            node->CutValueRight_ = node->RightNode_->vecBounds_[c].fLower;
            node->CutValueLeft_ = node->LeftNode_->vecBounds_[c].fUpper;
            node->CutValue_ = (node->CutValueLeft_ + node->CutValueRight_) / 2.0f;
            //
            // now recompute true bounding box as union of subtree boxes.
            // This is now faster having built the tree, being logarithmic in
            // N, not linear as would be from naive method.
            //
            for (int i = 0; i < Dim_; i++) {
                node->vecBounds_[i].fUpper = std::max(node->LeftNode_->vecBounds_[i].fUpper, node->RightNode_->vecBounds_[i].fUpper);
                node->vecBounds_[i].fLower = std::min(node->LeftNode_->vecBounds_[i].fLower, node->RightNode_->vecBounds_[i].fLower);
            }
        }
    }
    return (node);
}

template<typename T>
void KDTreeCL<T>::SpreadInCoordinate(int c, int l, int u, IntervalST<T>& interv)
{
    // return the minimum and maximum of the indexed data between l and u in
    // smin_out and smax_out.
    T smin, smax;
    T lmin, lmax;
    int i;

    smin = Data_[vecTreeLeaves_[l] * Dim_ + c];
    smax = smin;

    // process two at a time.
    for (i = l + 2; i <= u; i += 2) {
        lmin = Data_[vecTreeLeaves_[i - 1] * Dim_ + c];
        lmax = Data_[vecTreeLeaves_[i] * Dim_ + c];

        if (lmin > lmax) {
            std::swap(lmin, lmax);
            //      float t = lmin;
            //      lmin = lmax;
            //      lmax = t;
        }

        if (smin > lmin)
            smin = lmin;
        if (smax < lmax)
            smax = lmax;
    }
    // is there one more element?
    if (i == u + 1) {
        T last = Data_[vecTreeLeaves_[u] * Dim_ + c];
        if (smin > last)
            smin = last;
        if (smax < last)
            smax = last;
    }
    interv.fLower = smin;
    interv.fUpper = smax;
    //  printf("Spread in coordinate %d=[%f,%f]\n",c,smin,smax);
}

template<typename T>
void KDTreeCL<T>::SelectOnCoordinate(int c, int k, int l, int u)
{
    //
    //  Move indices in ind[l..u] so that the elements in [l .. k]
    //  are less than the [k+1..u] elmeents, viewed across dimension 'c'.
    //
    while (l < u) {
        int t = vecTreeLeaves_[l];
        int m = l;

        for (int i = l + 1; i <= u; i++) {
            if (Data_[vecTreeLeaves_[i] * Dim_ + c] < Data_[t * Dim_ + c]) {
                m++;
                std::swap(vecTreeLeaves_[i], vecTreeLeaves_[m]);
            }
        } // for i
        std::swap(vecTreeLeaves_[l], vecTreeLeaves_[m]);

        if (m <= k)
            l = m + 1;
        if (m >= k)
            u = m - 1;
    } // while loop
}

template<typename T>
int KDTreeCL<T>::SelectOnCoordinateValue(int c, T alpha, int l, int u)
{
    //
    //  Move indices in ind[l..u] so that the elements in [l .. return]
    //  are <= alpha, and hence are less than the [return+1..u]
    //  elmeents, viewed across dimension 'c'.
    //
    int lb = l, ub = u;

    while (lb < ub) {
        if (Data_[vecTreeLeaves_[lb] * Dim_ + c] <= alpha) {
            lb++; // good where it is.
        }
        else {
            std::swap(vecTreeLeaves_[lb], vecTreeLeaves_[ub]);
            ub--;
        }
    }

    // here ub=lb
    if (Data_[vecTreeLeaves_[lb] * Dim_ + c] <= alpha)
        return (lb);
    else
        return (lb - 1);

}

// void kdtree2::dump_data() {
//   int upper1, upper2;
//
//   upper1 = N;
//   upper2 = dim;
//
//   printf("Rearrange=%d\n",rearrange);
//   printf("N=%d, dim=%d\n", upper1, upper2);
//   for (int i=0; i<upper1; i++) {
//     printf("the_data[%d][*]=",i);
//     for (int j=0; j<upper2; j++)
//       printf("%f,",the_data[i][j]);
//     printf("\n");
//   }
//   for (int i=0; i<upper1; i++)
//     printf("Indexes[%d]=%d\n",i,ind[i]);
//   for (int i=0; i<upper1; i++) {
//     printf("data[%d][*]=",i);
//     for (int j=0; j<upper2; j++)
//       printf("%f,",(*data)[i][j]);
//     printf("\n");
//   }
// }


//
/// \brief search record substructure
//
/// one of these is created for each search.
/// this holds useful information  to be used
/// during the search
template<typename T>
class CSearchRecord
{
private:
    friend class KDTreeCL<T> ;
    friend class KDTreeNodeCL<T> ;

    std::vector<T>& qv;
    int dim;
    bool rearrange;
    unsigned int nn; // , nfound;
    T ballsize;
    int centeridx, correltime;

    KDTreeResultVectorCL<T>& result; // results
    const T* data;
    const std::vector<int>& ind;
    T infinity;

public:
    CSearchRecord(std::vector<T>& qv_in, KDTreeCL<T>& tree_in, KDTreeResultVectorCL<T>& result_in) :
        qv(qv_in), result(result_in), data(tree_in.InternalData_), ind(tree_in.vecTreeLeaves_), infinity(std::numeric_limits<T>::max())
    {
        dim = tree_in.Dim_;
        rearrange = tree_in.Rearrange_;
        ballsize = infinity;
        nn = 0;
    }
};

template<typename T>
void KDTreeCL<T>::GetNNearest_BruteForce(std::vector<T>& qv, int nn, KDTreeResultVectorCL<T>& result)
{

    result.clear();

    for (int i = 0; i < NumberOfPoints_; i++) {
        T dis = 0.0;
        KDTreeResultST<T> e;
        for (int j = 0; j < Dim_; j++) {
            dis += squared(Data_[i * Dim_ + j] - qv[j]);
        }
        e.Distance = dis;
        e.Index = i;
        result.push_back(e);
    }
    std::sort(result.begin(), result.end());

}

template<typename T>
void KDTreeCL<T>::GetNNearest(std::vector<T>& qv, int nn, KDTreeResultVectorCL<T>& result)
{
    CSearchRecord<T> sr(qv, *this, result);
    std::vector<T> vdiff(Dim_, 0.0);

    result.clear();

    sr.centeridx = -1;
    sr.correltime = 0;
    sr.nn = nn;

    Root_->search(sr);

    if (SortResults_)
        std::sort(result.begin(), result.end());

}
// search for n nearest to a given query vector 'qv'.

template<typename T>
void KDTreeCL<T>::GetNNearest_AroundPoint(int idxin, int correltime, int nn, KDTreeResultVectorCL<T>& result)
{
    std::vector<T> qv(Dim_); //  query vector

    result.clear();

    for (int i = 0; i < Dim_; i++) {
        qv[i] = Data_[idxin * Dim_ + i];
    }
    // copy the query vector.

    {
        CSearchRecord<T> sr(qv, *this, result);
        // construct the search record.
        sr.centeridx = idxin;
        sr.correltime = correltime;
        sr.nn = nn;
        Root_->search(sr);
    }

    if (SortResults_)
        std::sort(result.begin(), result.end());

}

template<typename T>
void KDTreeCL<T>::GetNearestByRadius(std::vector<T>& qv, T r2, KDTreeResultVectorCL<T>& result)
{
    // search for all within a ball of a certain radius
    CSearchRecord<T> sr(qv, *this, result);
    std::vector<T> vdiff(Dim_, 0.0);

    result.clear();

    sr.centeridx = -1;
    sr.correltime = 0;
    sr.nn = 0;
    sr.ballsize = r2;

    Root_->search(sr);

    if (SortResults_)
        std::sort(result.begin(), result.end());

}

template<typename T>
int KDTreeCL<T>::GetNumberOfNeighborsWithin(std::vector<T>& qv, T r2)
{
    // search for all within a ball of a certain radius
    {
        KDTreeResultVectorCL<T> result;
        CSearchRecord<T> sr(qv, *this, result);

        sr.centeridx = -1;
        sr.correltime = 0;
        sr.nn = 0;
        sr.ballsize = r2;

        Root_->search(sr);
        return (result.size());
    }

}

template<typename T>
void KDTreeCL<T>::GetNearestByRadius_AroundPoint(int idxin, int correltime, T r2, KDTreeResultVectorCL<T>& result)
{
    std::vector<T> qv(Dim_); //  query vector

    result.clear();

    for (int i = 0; i < Dim_; i++) {
        qv[i] = Data_[idxin * Dim_ + i];
    }
    // copy the query vector.

    {
        CSearchRecord<T> sr(qv, *this, result);
        // construct the search record.
        sr.centeridx = idxin;
        sr.correltime = correltime;
        sr.ballsize = r2;
        sr.nn = 0;
        Root_->search(sr);
    }

    if (SortResults_)
        std::sort(result.begin(), result.end());

}

template<typename T>
int KDTreeCL<T>::GetNumberOfNeighbors_AroundPoint(int idxin, int correltime, T r2)
{
    std::vector<T> qv(Dim_); //  query vector


    for (int i = 0; i < Dim_; i++) {
        qv[i] = Data_[idxin * Dim_ + i];
    }
    // copy the query vector.

    {
        KDTreeResultVectorCL<T> result;
        CSearchRecord<T> sr(qv, *this, result);
        // construct the search record.
        sr.centeridx = idxin;
        sr.correltime = correltime;
        sr.ballsize = r2;
        sr.nn = 0;
        Root_->search(sr);
        return (result.size());
    }

}

template<typename T>
void KDTreeCL<T>::SerializeTree(std::vector<KDTreeSerialNode<T> > &nodes, std::vector<T> &SerialData, std::vector<IntervalST<T> > &VecBounds)
{
    for (int i = 0; i < NumberOfPoints_ * Dim_; ++i) {
        if (Rearrange_)
            SerialData.push_back(RearrangedData_[i]);
        else
            SerialData.push_back(Data_[i]);
    }
    Root_->SerializeNodes(nodes, Dim_, VecBounds);
}

//
//        KDTREE2_SERIALNODE
//

//constructor
template<typename T>
KDTreeSerialNode<T>::KDTreeSerialNode(const KDTreeNodeCL<T> &iNode)
{
    this->leftChild = -1;
    this->rightChild = -1;
    l = iNode.l;
    u = iNode.u;
    //copy by hand all information in KDTreeNodeCL except vecBounds_ and LeftNode and RightNode
    CutDimension_ = iNode.CutDimension_;
    CutValue_ = iNode.CutValue_;
    CutValueLeft_ = iNode.CutValueLeft_;
    CutValueRight_ = iNode.CutValueRight_;
}

template<typename T>
KDTreeSerialNode<T>::KDTreeSerialNode()
{
    this->leftChild = -1;
    this->rightChild = -1;
    l = 0;
    u = 0;
}

//
//        KDTREE2_NODE implementation
//

// constructor
template<typename T>
KDTreeNodeCL<T>::KDTreeNodeCL(const KDTreeSerialNode<T> &iNode, int dim)
{
    l = iNode.l;
    u = iNode.u;
    //copy by hand all information in KDTreeNodeCL except vecBounds_ and LeftNode and RightNode
    CutDimension_ = iNode.CutDimension_;
    CutValue_ = iNode.CutValue_;
    CutValueLeft_ = iNode.CutValueLeft_;
    CutValueRight_ = iNode.CutValueRight_;
    LeftNode_ = NULL;
    RightNode_ = NULL;
}

template<typename T>
KDTreeNodeCL<T>::KDTreeNodeCL(int dim) :
    vecBounds_(dim)
{
    LeftNode_ = RightNode_ = NULL;
    //
    // all other construction is handled for real in the
    // kdtree2 building operations.
    //
}

// destructor
template<typename T>
KDTreeNodeCL<T>::~KDTreeNodeCL()
{
    if (LeftNode_ != NULL)
        delete LeftNode_;
    if (RightNode_ != NULL)
        delete RightNode_;
    // maxbox and minbox
    // will be automatically deleted in their own destructors.
}

template<typename T>
void KDTreeNodeCL<T>::search(CSearchRecord<T>& sr)
{
    // the core search routine.
    // This uses true distance to bounding box as the
    // criterion to search the secondary node.
    //
    // This results in somewhat fewer searches of the secondary nodes
    // than 'search', which uses the vdiff vector,  but as this
    // takes more computational time, the overall performance may not
    // be improved in actual run time.
    //

    if ((LeftNode_ == NULL) && (RightNode_ == NULL)) {
        // we are on a terminal node
        if (sr.nn == 0) {
            process_terminal_node_fixedball(sr);
        }
        else {
            process_terminal_node(sr);
        }
    }
    else {
        KDTreeNodeCL<T> *ncloser, *nfarther;

        T extra;
        T qval = sr.qv[CutDimension_];
        // value of the wall boundary on the cut dimension.
        if (qval < CutValue_) {
            ncloser = LeftNode_;
            nfarther = RightNode_;
            extra = CutValueRight_ - qval;
        }
        else {
            ncloser = RightNode_;
            nfarther = LeftNode_;
            extra = qval - CutValueLeft_;
        };

        if (ncloser != NULL)
            ncloser->search(sr);

        if ((nfarther != NULL) && (squared(extra) < sr.ballsize)) {
            // first cut
            if (nfarther->box_in_search_range(sr)) {
                nfarther->search(sr);
            }
        }
    }
}

template<typename T>
void KDTreeNodeCL<T>::SerializeNodes(std::vector<KDTreeSerialNode<T> > &nodes, int dim, std::vector<IntervalST<T> > &vecBounds)
{
    int myIndex = nodes.size();

    nodes.push_back(KDTreeSerialNode<T> (*this));
    for (int i = 0; i < dim; i++)
        vecBounds.push_back(vecBounds_[i]);

    // the index of the Left and Right Child
    if (LeftNode_) {
        nodes[myIndex].leftChild = nodes.size();
        LeftNode_->SerializeNodes(nodes, dim, vecBounds);
    }
    if (RightNode_) {
        nodes[myIndex].rightChild = nodes.size();
        RightNode_->SerializeNodes(nodes, dim, vecBounds);
    }
}

template<typename T>
inline T dis_from_bnd(T x, T amin, T amax)
{
    if (x > amax) {
        return (x - amax);
    }
    else if (x < amin)
        return (amin - x);
    else
        return 0.0;

}

template<typename T>
inline bool KDTreeNodeCL<T>::box_in_search_range(CSearchRecord<T>& sr)
{
    //
    // does the bounding box, represented by minbox[*],maxbox[*]
    // have any point which is within 'sr.ballsize' to 'sr.qv'??
    //

    int dim = sr.dim;
    T dis2 = 0.0;
    T ballsize = sr.ballsize;
    for (int i = 0; i < dim; i++) {
        dis2 += squared(dis_from_bnd(sr.qv[i], vecBounds_[i].fLower, vecBounds_[i].fUpper));
        if (dis2 > ballsize)
            return (false);
    }
    return (true);
}

template<typename T>
void KDTreeNodeCL<T>::process_terminal_node(CSearchRecord<T>& sr)
{
    int centeridx = sr.centeridx;
    int correltime = sr.correltime;
    unsigned int nn = sr.nn;
    int dim = sr.dim;
    T ballsize = sr.ballsize;
    //
    bool rearrange = sr.rearrange;
    const T* data = sr.data;

    const bool debug = false;

    if (debug) {
        printf("Processing terminal node %d, %d\n", l, u);
        std::cout << "Query vector = [";
        for (int i = 0; i < dim; i++)
            std::cout << sr.qv[i] << ',';
        std::cout << "]\n";
        std::cout << "nn = " << nn << '\n';
    }

    for (int i = l; i <= u; i++) {
        int indexofi; // sr.ind[i];
        T dis;
        bool early_exit;

        if (rearrange) {
            early_exit = false;
            dis = 0.0;
            for (int k = 0; k < dim; k++) {
                dis += squared(data[i * dim + k] - sr.qv[k]);
                if (dis > ballsize) {
                    early_exit = true;
                    break;
                }
            }
            if (early_exit)
                continue; // next iteration of mainloop
            // why do we do things like this?  because if we take an early
            // exit (due to distance being too large) which is common, then
            // we need not read in the actual point index, thus saving main
            // memory bandwidth.  If the distance to point is less than the
            // ballsize, though, then we need the index.
            //
            indexofi = sr.ind[i];
        }
        else {
            //
            // but if we are not using the rearranged data, then
            // we must always
            indexofi = sr.ind[i];
            early_exit = false;
            dis = 0.0;
            for (int k = 0; k < dim; k++) {
                dis += squared(data[indexofi * dim + k] - sr.qv[k]);
                if (dis > ballsize) {
                    early_exit = true;
                    break;
                }
            }
            if (early_exit)
                continue; // next iteration of mainloop
        } // end if rearrange.

        if (centeridx > 0) {
            // we are doing decorrelation interval
            if (std::abs(indexofi - centeridx) < correltime)
                continue; // skip this point.
        }

        // here the point must be added to the list.
        //
        // two choices for any point.  The list so far is either
        // undersized, or it is not.
        //
        if (sr.result.size() < nn) {
            KDTreeResultST<T> e;
            e.Index = indexofi;
            e.Distance = dis;
            sr.result.push_element_and_heapify(e);
            if (debug)
                std::cout << "unilaterally pushed dis=" << dis;
            if (sr.result.size() == nn)
                ballsize = sr.result.max_value();
            // Set the ball radius to the largest on the list (maximum priority).
            if (debug) {
                std::cout << " ballsize = " << ballsize << "\n";
                std::cout << "sr.result.size() = " << sr.result.size() << '\n';
            }
        }
        else {
            //
            // if we get here then the current node, has a squared
            // distance smaller
            // than the last on the list, and belongs on the list.
            //
            KDTreeResultST<T> e;
            e.Index = indexofi;
            e.Distance = dis;
            ballsize = sr.result.replace_maxpri_elt_return_new_maxpri(e);
            if (debug) {
                std::cout << "Replaced maximum dis with dis=" << dis << " new ballsize =" << ballsize << '\n';
            }
        }
    } // main loop
    sr.ballsize = ballsize;
}

template<typename T>
void KDTreeNodeCL<T>::process_terminal_node_fixedball(CSearchRecord<T>& sr)
{
    int centeridx = sr.centeridx;
    int correltime = sr.correltime;
    int dim = sr.dim;
    T ballsize = sr.ballsize;
    //
    bool rearrange = sr.rearrange;
    const T* data = sr.data;

    for (int i = l; i <= u; i++) {
        int indexofi = sr.ind[i];
        T dis;
        bool early_exit;

        if (rearrange) {
            early_exit = false;
            dis = 0.0;
            for (int k = 0; k < dim; k++) {
                dis += squared(data[i * dim + k] - sr.qv[k]);
                if (dis > ballsize) {
                    early_exit = true;
                    break;
                }
            }
            if (early_exit)
                continue; // next iteration of mainloop
            // why do we do things like this?  because if we take an early
            // exit (due to distance being too large) which is common, then
            // we need not read in the actual point index, thus saving main
            // memory bandwidth.  If the distance to point is less than the
            // ballsize, though, then we need the index.
            //
            indexofi = sr.ind[i];
        }
        else {
            //
            // but if we are not using the rearranged data, then
            // we must always
            indexofi = sr.ind[i];
            early_exit = false;
            dis = 0.0;
            for (int k = 0; k < dim; k++) {
                dis += squared(data[indexofi * dim + k] - sr.qv[k]);
                if (dis > ballsize) {
                    early_exit = true;
                    break;
                }
            }
            if (early_exit)
                continue; // next iteration of mainloop
        } // end if rearrange.

        if (centeridx > 0) {
            // we are doing decorrelation interval
            if (std::abs(indexofi - centeridx) < correltime)
                continue; // skip this point.
        }

        {
            KDTreeResultST<T> e;
            e.Index = indexofi;
            e.Distance = dis;
            sr.result.push_back(e);
        }

    }
}

#ifdef _PAR
template<typename T>
void KDTreeCL<T>::MPITypesCL::CreateTypes(int dim)
{
    // IntervalSTType
    if ( IntervalSTType==ProcCL::NullDataType ){
        IntervalST<T> interval;
        ProcCL::DatatypeT type1[7] = { ProcCL::MPI_TT<T>::dtype, ProcCL::MPI_TT<T>::dtype };
        int blocklen1[9] = { 1, 1 };
        ProcCL::AintT disp1[2], orig1;
        orig1= ProcCL::Get_address(&interval);
        disp1[0]= ProcCL::Get_address(&interval.fLower) - orig1;
        disp1[1]= ProcCL::Get_address(&interval.fUpper) - orig1;
        IntervalSTType= ProcCL::CreateStruct(2, blocklen1, disp1, type1);
        ProcCL::Commit(IntervalSTType);
    }

    ///KDTreeSerialNodeType
    if ( KDTreeSerialNodeType==ProcCL::NullDataType ){
        KDTreeSerialNode<double> serialnode;
        ProcCL::DatatypeT type2[8] = {
                ProcCL::MPI_TT<int>::dtype, ProcCL::MPI_TT<int>::dtype, ProcCL::MPI_TT<int>::dtype,
                ProcCL::MPI_TT<T>::dtype,   ProcCL::MPI_TT<T>::dtype,   ProcCL::MPI_TT<T>::dtype,
                ProcCL::MPI_TT<int>::dtype, ProcCL::MPI_TT<int>::dtype
                };
        int blocklen2[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
        ProcCL::AintT disp2[8], orig2;
        orig2= ProcCL::Get_address(&serialnode);
        disp2[0]= ProcCL::Get_address( &serialnode.leftChild)     - orig2;
        disp2[1]= ProcCL::Get_address( &serialnode.rightChild)    - orig2;
        disp2[2]= ProcCL::Get_address( &serialnode.CutDimension_) - orig2;
        disp2[3]= ProcCL::Get_address( &serialnode.CutValue_)     - orig2;
        disp2[4]= ProcCL::Get_address( &serialnode.CutValueLeft_) - orig2;
        disp2[5]= ProcCL::Get_address( &serialnode.CutValueRight_)- orig2;
        disp2[6]= ProcCL::Get_address( &serialnode.l)             - orig2;
        disp2[7]= ProcCL::Get_address( &serialnode.u)             - orig2;
        KDTreeSerialNodeType= ProcCL::CreateStruct( 8, blocklen2, disp2, type2);
        ProcCL::Commit(KDTreeSerialNodeType);
    }
}

template<typename T>
void KDTreeCL<T>::MPITypesCL::FreeTypes()
{
    if ( IntervalSTType!=ProcCL::NullDataType )
        ProcCL::Free( IntervalSTType);
    if ( KDTreeSerialNodeType!=ProcCL::NullDataType )
        ProcCL::Free( KDTreeSerialNodeType);
}

/** Each processor owns its kd-tree. These trees are serialized and gathered on
    all processors. The result is returned in \a result
    \param result a vector of P kd-trees, i.e., result[p] stores kd-tree of processor p
*/
template<typename T>
void KDTreeCL<T>::SerializeAndGather(std::vector<KDTreeCL<T> >& result)
{
    result.clear();

    // Serialize the own tree ...
    std::vector< KDTreeSerialNode <T> > nodes;
    std::vector< T > data;
    std::vector< IntervalST<T> > vecBounds;
    SerializeTree(nodes, data, vecBounds);

    // Gather serialized data. Since Gather expects the same message length
    // (same number of nodes) this operation has to be done by Isend and Irecv :-(
    MPITypesCL types(Dim_);
    std::vector<ProcCL::RequestT> req(3*(ProcCL::Size()-1)); // Send 3 messages to all other processors
    std::vector< std::vector< KDTreeSerialNode <T> >* > all_nodes(ProcCL::Size(), 0);
    std::vector< std::vector< T >* >                    all_data(ProcCL::Size(), 0);
    std::vector< std::vector< IntervalST<T> >* >        all_vecBounds(ProcCL::Size(), 0);

    // Start communication
    int reqPos=0;
    for (int p=0; p<ProcCL::Size(); ++p){
        if ( p!=ProcCL::MyRank() ){
            req[reqPos++]= ProcCL::Isend(nodes,     types.KDTreeSerialNodeType, p, 4001);
            req[reqPos++]= ProcCL::Isend(data,      ProcCL::MPI_TT<T>::dtype,   p, 4002);
            req[reqPos++]= ProcCL::Isend(vecBounds, types.IntervalSTType,       p, 4003);
        }
    }

    // Receive data
    ProcCL::StatusT stat;
    int count;
    for (int p=0; p<ProcCL::Size(); ++p){
        if ( p!=ProcCL::MyRank() ){     // Receive
            ProcCL::Probe( p, 4001, stat);
            count= ProcCL::GetCount( stat, types.KDTreeSerialNodeType);
            all_nodes[p]= new std::vector< KDTreeSerialNode <T> >( count);
            ProcCL::Recv( Addr(*(all_nodes[p])), count, types.KDTreeSerialNodeType, p, 4001);

            ProcCL::Probe( p, 4002, stat);
            count= ProcCL::GetCount( stat, ProcCL::MPI_TT<T>::dtype);
            all_data[p]= new std::vector< T >( count);
            ProcCL::Recv( Addr(*(all_data[p])), count, ProcCL::MPI_TT<T>::dtype, p, 4002);

            ProcCL::Probe( p, 4003, stat);
            count= ProcCL::GetCount( stat, types.IntervalSTType);
            all_vecBounds[p]= new std::vector< IntervalST<T> >( count);
            ProcCL::Recv( Addr(*(all_vecBounds[p])), count, types.IntervalSTType, p, 4003);
        }
        else {  // Umbiegen der Pointer
            all_nodes[ProcCL::MyRank()]    = &nodes;
            all_data[ProcCL::MyRank()]     = &data;
            all_vecBounds[ProcCL::MyRank()]= &vecBounds;
        }
    }

    // Create trees
    for (int p=0; p<ProcCL::Size(); ++p){
        result.push_back( KDTreeCL<T>(all_nodes[p], all_data[p], Dim_, all_vecBounds[p]));
    }

    // Free memory
    for (int p=0; p<ProcCL::Size(); ++p){
           if ( p!=ProcCL::MyRank() ){
               delete (all_nodes[p]);
               delete (all_data[p]);
               delete (all_vecBounds[p]);
           }
    }
}

#endif
} // end of namespace
