/***************************************************************************
*  File:    partitioner.h                                                  *
*  Content: - Interface for partitioners                                   *
*           - Metis partitioner
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           18.09.2009                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file partitioner.h
/// \brief Interface for all partitioners

#ifndef _DROPS_PARTITIONER_H
#define _DROPS_PARTITIONER_H

#include <vector>
#include "misc/utils.h"
#include "parallel/parallel.h"

namespace DROPS
{

// fwd declaration
template <typename IdxT> class CSRPartitionerCL;

/// \brief Class for storing parallel graphs in CSR (compressed storage format)
template <typename IdxT>
class ParallelCSRGraphCL
{
  public:
    typedef std::vector<IdxT> VectorT;      ///< Vector-type to store data

  private:
    VectorT adjpointer_;                    ///< pointer to adjacency list
    VectorT adjacencies_;                   ///< adjacencies
    VectorT vtxdist_;                       ///< number of stored nodes on all processes
    VectorT vwgt_;                          ///< weight of nodes
    VectorT adjwgt_;                        ///< weight of adjacencies

  public:
    /// \brief Constructor
    ParallelCSRGraphCL() {}
    /// \brief Destructor
    ~ParallelCSRGraphCL() {}

    /// \brief Free memory
    void Clear(){
        adjpointer_.clear(); adjacencies_.clear();
        vtxdist_.clear();
        vwgt_.clear(); adjwgt_.clear;
    }

    /// PArtitioner may access all data
    friend class CSRPartitionerCL<IdxT>;
};

/// \enum PartMethod Various methods to be used to partition a graph
enum PartMethod{
    KWay,       ///< Multilevel method
    Recursive,  ///< bisection
    Adaptive,   ///< adaptive recomputation of graph partitioning
    Identity,   ///< Take partition as it is
    NoMig       ///< No migration!
};

/// \brief Base (abstract) class for all derived partitioner
/** This class stores (parallel) CSR graphs and provides an interface
    for all partitioners.
 */
template <typename IdxT>
class CSRPartitionerCL
{
  public:
    typedef std::vector<int>                  PartitionT;   ///< return type for distribution of a node
    typedef ParallelCSRGraphCL<IdxT>::VectorT VectorT;      ///< Vector-type to store data

  protected:
    ParallelCSRGraphCL<IdxT> graph_;          ///< the graph
    PartitionT               part_;           ///< partitioning of the graph
    PartMethod               method_;         ///< used method to partition a graph

  public:
    /// \brief Constructor
    CSRPartitionerCL() {}
    /// \brief Destructor
    virtual ~CSRPartitionerCL() { Clear(); }

    /// \brief Get number of adjacencies
    size_t GetNumEdges() { return adjacencies_.size(); }
    /// \brief Get number of vertices
    size_t GetNumVertices() { return adjpointer_.size(); }

    /// \brief Return used method
    PartMethod GetMethod() const{ return method_; }
    /// \brief Set Method, overload this method to perform error checking
    virtual void SetMethod( const PartMethod&) =0;

    /// \brief Partition a local stored graph
    virtual const PartitionT& PartSerial() = 0;
    /// \brief Partition a parallel stored graph
    virtual const PartitionT& PartParallel()=0;
    /// \brief Get number of cutted edges
    virtual size_t GetEdgeCut() const =0;

    /// \name Get constant references on graph data
    //@{
    const VectorT& GetAdjPtr()    const { return graph_.adjpointer_; }
    const VectorT& GetAdj()       const { return graph_.adjacencies_; }
    const VectorT& GetVtxDist()   const { return graph_.vtxdist_; }
    const VectorT& GetVtxWeight() const { return graph_.vwgt_; }
    const VectorT& GetAdjWeight() const { return graph_.adjwgt_; }
    const PartitionT& GetPart()   const { return part_; }
    //@}

    /// \name Get references on graph data
    //@{
    VectorT& GetAdjPtr()    { return graph_.adjpointer_; }
    VectorT& GetAdj()       { return graph_.adjacencies_; }
    VectorT& GetVtxDist()   { return graph_.vtxdist_; }
    VectorT& GetVtxWeight() { return graph_.vwgt_; }
    VectorT& GetAdjWeight() { return graph_.adjwgt_; }
    //@}
};


/// \brief Class that uses Metis and ParMetis to partition a graph
template <typename IdxT=idxtype>
class MetisPartitionerCL : public CSRPartitionerCL<IdxT>
{
  public:
    typedef MetisPartotionerCL<IdxT> self;  ///< this class
    typedef CSRPartitionerCL<IdxT>   base;  ///< base class

  private:

  public:
    /// \brief Set used method and check, if the method can be performed by Metis
    void SetMethod(const PartMethod& method) {
        if (method!=KWay || method!=Recursive || method!=Adaptive || method!=Identity || method!=NoMig){
            throw DROPSErrCL("MetisPartitionersCL::GetMethod: No such partitioning method known for Metis");
        }
        method_= method;
    }

    const PartitionT& PartSerial();
};

}   // end of namespace

#include "parallel/partitioner.tpp"
#endif
