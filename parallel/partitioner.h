/// \file partitionerclass.h
/// \brief Header file for the class that implements different types of 
/// partitioners. Basically each new partitioner used needs a modification on
/// the input vectors to suit their corresponding interface. (metis is an exception)
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Alin Bastea
/// Begin: 18.01.2010


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

#ifndef DROPS_PARTITIONER_H
#define DROPS_PARTITIONER_H

#include "parallel/parallel.h"
#include "parallel/partime.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "misc/utils.h"
#ifdef _ZOLTAN
#  include <zoltan.h>
#endif
#include <parmetis.h>
#include <metis.h>

namespace DROPS
{
typedef idxtype* IndexArray;        ///< idxtype is defined as Integer

/// \enum PartOption tells which partitioner should pe used to compute graph partition problem
enum PartOption{
                metis,          ///< Parmetis
                zoltan,         ///< Zoltan
                scotch          ///< Scotch
};

/// \enum PartMethod tells which method should be used to compute graph partition problem
enum PartMethod{
    KWay,       ///< Multilevel method
    Recursive,  ///< bisection
    Adaptive,   ///< adaptive re-computation of graph partitioning
    Identity,   ///< Take partition as it is
    NoMig       ///< No migration!
};

/// \brief structure of the inputs and the outputs of the partitioners
struct GraphST
{
    IndexArray     xadj,                                        ///< Starting index in the array adjncy, where node[i] writes its neighbors
                   adjncy,                                      ///< adjacencies of the nodes
                   vtxdist,                                     ///< number of nodes, that is stored by all procs
                   vwgt,                                        ///< weight of the Nodes
                   adjwgt,                                      ///< weight of the edges
                   part;                                        ///< resulting array, where each node should be send to
    float          ubvec;                                       ///< quality of graph partitioning
    float*         xyz;                                         ///< geometric information about the nodes (==barymetric Center of the Tetras)
    int            myVerts;                                     ///< number of vertices on this proc
    int            myAdjs;                                      ///< number of Adjacencies
    int            edgecut;                                     ///< number of edges, that are cut by ParMETIS
    int            movedMultiNodes;                             ///< number of multinodes, that are moved by last migration
    bool           geom;                                        ///< flag, that indicates, if the geometrical datas should be used
    idxtype        myfirstVert;                                 ///< first vertex # on the current proc
    GraphST():xadj(0), adjncy(0), vtxdist(0), vwgt(0), adjwgt(0), part(0), xyz(0), myVerts(0), myAdjs(0),edgecut(0), movedMultiNodes(0),geom(0), myfirstVert(0){} ///< constructor
    ~GraphST();                                                 ///< destructor
    void Resize (int numadj, int myVerts, bool geom);           ///< The next method will allocate memory for most of the arrays in the struct
    void ResizeVtxDist();
    void Clear();                                               ///< Liberating the memory used for storing the arrays in the struct

    int GetProc(int globalid)
    {
       int i;
       for (i=0;i< ProcCL::Size();i++)
               if ( (globalid>=vtxdist[i]) && (globalid<vtxdist[i+1]) )
                   return i;
       throw DROPSErrCL("GraphST::GetProc: Wrong node id");
       return -1;
    }
};

/// \brief Abstract partitioner class. Through the factory design pattern it will be inherited by different
///  partitioners
/** All derived class must implement the functions CreateGraph, PartGraphPar, and PartGraphSer */
class PartitionerCL
{
  private:
    GraphST* allArrays_;                        ///< the input and output arrays
  public:
    PartitionerCL();                            ///< Constructor
    virtual ~PartitionerCL();                   ///< Destructor
    virtual void CreateGraph() =0;              ///< Create the graph using the input data
    virtual void PartGraphPar() =0;             ///< Parallel partitioning of the graph (highly used in applications)
    virtual void PartGraphSer( int master) =0;  ///< First serial partitioning of the graph needed to distribute the
                                                /// different partitioning components to the processors. (usually used once/application)
    static PartitionerCL* newPartitioner(PartOption partOption,float quality,PartMethod method); ///< Factory function for the factory design pattern implementation
    GraphST& GetGraph();                        ///< Getter for the graph structure
};//end of the abstract class PartitionerCL

/// \brief Derived partitioner class from the PartitionerCL class ===> implements the ParMetis graph partitioner
class ParMetisCL : public PartitionerCL
{
  private:
    PartMethod meth_;                   ///< attribute used to partition the graph with ParMetis
    float quality_;                     ///< quality of the Adaptive method partitioning
  public:
    ParMetisCL(PartMethod meth,float ubvec, float quality=1000.0);  //< Constructor
    ~ParMetisCL();                      ///< Implicit Destructor
    void CreateGraph();                 ///< Method that helps with the protocol communication between DROPS and the partitioner
    void PartGraphPar();                ///< Implemented virtual method of the abstract parent class PartitionerCL
    void PartGraphSer( int master);     ///< Implemented virtual method of the abstract parent class PartitionerCL
};//end of the derived class ParMetisCL

#ifdef _ZOLTAN
/// \brief Derived partitioner from the Partitioner class ===> implements the Zoltan graph partitioner
class ZoltanCL : public PartitionerCL
{
  private:
    Zoltan_Struct *zz;
    PartMethod meth_;// For Zoltan the only possibility is a parallel or a serial partitioning for now
  public:
    ZoltanCL(PartMethod meth);               ///< Constructor
    ~ZoltanCL();                             ///< Destructor
    /// \brief Returns the total number of vertices on a processor
    static int get_number_of_vertices(void *, int *);
    /// \brief Returns the vertices list
    static void get_vertex_list(void *, int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int *);
    /// \ brief Returns the total number of edges
    static void get_num_edges_list(void *, int , int, int , ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *);
    /// \ brief Returns the list of edges of each processor
    static void get_edge_list(void *, int , int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *, ZOLTAN_ID_PTR , int *, int , float *, int *);
    void CreateGraph();                     ///< Method that helps with the protocol communication between DROPS and the partitioner
    void PartGraphPar();                    ///< Implemented virtual method of the abstract parent class PartitionerCL
    void PartGraphSer( int master);         ///< Implemented virtual method of the abstract parent class PartitionerCL
};//end of the derived class ZoltanCL
#endif

/// \brief Derived partitioner from the Partitioner class ===> implements the Scotch graph partitioner
class ScotchCL : public PartitionerCL
{
  public:
    ScotchCL() : PartitionerCL(){};                 ///< Constructor
    ~ScotchCL();                                    ///< Destructor
};//end of the derived class ScotchCL
}// end of namespace
#endif //DROPS_PARTITIONER_H
