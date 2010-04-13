/// \file partitionerclass.cpp
/// \brief Implementation of the methods of the classes inside partitioner.h file
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

#include "parallel/partitioner.h"

namespace DROPS
{

//Implementation of the methods of the GraphST structure
//========================================================
/// \brief Destructor of the GraphST structure
GraphST::~GraphST()
{
    Clear();
}

/// \brief Method used to allocate memory for the vtxdist vector in the struct
void GraphST::ResizeVtxDist()
{
    if (vtxdist!=0) delete[] vtxdist;
    vtxdist = new int[ProcCL::Size()+1];
}

/**This method will allocate memory for most of the arrays in the GraphST structure
 \param numadj  number of neighbors
        myVerts number of verticies
        geom    flag for geometry
 */
void GraphST::Resize(int numadj, int myVerts, bool geom)
{
    if (xadj!=0) delete[] xadj;
    xadj = new idxtype[myVerts+1];
    if (adjncy!=0) delete[] adjncy;
    adjncy = new idxtype[numadj];
    if (vwgt!=0) delete[] vwgt;
    vwgt = new idxtype[myVerts];
    if (adjwgt!=0) delete[] adjwgt;
    adjwgt = new idxtype[numadj];
    if (geom)
    {
        if (xyz!=0) delete[] xyz;
        xyz = new float[3*myVerts];
    }
    if (part!=0) delete[] part;
    part = new idxtype[myVerts];
}
/// \brief This method will liberate the memory used for storing the arrays in the structure
void GraphST::Clear()
{
    /** Free all allocated arrays */
    if (xadj!=0)    delete[] xadj;
    if (adjncy!=0)  delete[] adjncy;
    if (vtxdist!=0) delete[] vtxdist;
    if (part!=0)    delete[] part;
    if (vwgt!=0)    delete[] vwgt;
    if (adjwgt!=0)  delete[] adjwgt;
    if (xyz!=0)     delete[] xyz;
    xadj=0; adjncy=0; vtxdist=0; part=0;
    vwgt=0; adjwgt=0; xyz=0;
}
//End of the implementations of the methods in the structure GraphST
//============================================================================================

//Implementation of the methods of the parent class partitioner
//===========================================================================================

/// \brief Constructor of the base partitioner class (allocates memory for the structure attribute)
PartitionerCL::PartitionerCL()
{
    allArrays_ = new GraphST();
}

/// \brief Destructor of the base partitioner class (empties memory)
PartitionerCL::~PartitionerCL()
{
    delete allArrays_;
}

/// \brief Getter/setter for the structure member inside the base partitioner class
GraphST& PartitionerCL::GetGraph()
{
    return *allArrays_;
}

/**Factory method that based on the partOption parameter allocates memory for the chosen type of object
  (it implements the Factory design pattern)
 \param partOption  type of partitioner used
        quality     quality of the partitioning
        method      partition method invoked by the partitioner
 */
PartitionerCL* PartitionerCL::newPartitioner(PartOption partOption, float quality, PartMethod method)
{
    PartitionerCL* pa = NULL;
     switch (partOption)
     {
         case metis:
         {
             pa = new ParMetisCL(method, quality);
             break;
         }
         case zoltan:
         {
#ifdef _ZOLTAN
             pa = new ZoltanCL(method);
#else
             throw DROPSErrCL("PartitionerCL::newPartitioner: Zoltan library is not included");
#endif
             break;
         }
         case scotch:
         {
             /*pa = new ScotchCL();*/
             break;
         }
     }
     return pa;
}

//End of the implementations of the methods in the abstract class PartitionerCL
//============================================================================================

//Implementation of the methods in the derived class ParMetisCL
//=================================================================================================

/** Constructor of the derived class ParMetisCL (parent PartitionerCL)
 \param ubvec       type of partitioner used
        quality     quality of the partitioning
        meth        quality parameter
 */
ParMetisCL::ParMetisCL(PartMethod meth, float ubvec, float quality) : PartitionerCL()
{   meth_=meth;
    GetGraph().ubvec = ubvec;
    quality_ = quality;
}

/// \brief Destructor of the derived class ParMetisCL (parent PartitionerCL)
ParMetisCL::~ParMetisCL()
{

}

/// \brief For ParMetis partitioner the next function has no real value but it must be implemented for
/// other types of partitioners (as a manipulation of the inputs)
void ParMetisCL::CreateGraph()
{

}

/** Serial partitioning
 \param master      rank of the master thread
 */
void ParMetisCL::PartGraphSer(int master)
{
    if ( ProcCL::MyRank()!=master)
        return;
    
    Comment("- Start calculate LoadBalanace with METIS-"<<(meth_ == KWay ? "Kway" : "Recursive")<<std::endl, DebugLoadBalC);

    if (GetGraph().myVerts != GetGraph().vtxdist[DynamicDataInterfaceCL::InfoProcs()])
        Comment("LoadBalCL: PartGraphSer: This procedure is called by a proc, that does not have all nodes!"<<std::endl, DebugLoadBalC);

    int    wgtflag    = 3,                                      ///< Weights on vertices and adjacencies are given
           numflag    = 0,                                      ///< Numbering of verts starts by 0 (C-Style)
           nparts     = DynamicDataInterfaceCL::InfoProcs(),    ///< Number of subdomains (per proc one)
           n          = GetGraph().myVerts,
           options[5] = {0,0,0,0,0};                            ///< Default options
    //Depending on the method chosen KWay or Recursive the appropriate parmetis function is called
    if (meth_ == KWay)
        METIS_PartGraphKway(      &n, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt, &wgtflag, &numflag,  &nparts, options,&GetGraph().edgecut, GetGraph().part);
    else if (meth_==Recursive)
        METIS_PartGraphRecursive( &n, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt, &wgtflag, &numflag,  &nparts, options,&GetGraph().edgecut,GetGraph().part);

    Comment("  * Number of Edgecut: "<<edgecut_<<std::endl, DebugLoadBalC);
}

/// \brief Parallel partitioning
void ParMetisCL::PartGraphPar()
    /** This procedure uses ParMetis in order to compute a partitioning of the dual
        reduced graph, that has been set up with the member function
        CreateDualRedGraph. The previous distribution of the vertices among the
        processors is used.
        \pre CreateDualRedGraph() must have been called before
    */
    {
        
        Comment("- Start calculate LoadBalanace with ParMETIS-"<<(meth_ == Adaptive ? "AdaptiveRepart")<<(meth_ == KWay ? "KWay")<<(meth_ == Identity ? "Identity")<<std::endl,DebugLoadBalC);

        Assert( GetGraph().xadj && GetGraph().adjncy_ && GetGraph().vwgt && GetGraph().adjwgt,
                DROPSErrCL("LoadBalCL::PartGraphPar: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
                DebugLoadBalC);

        int    wgtflag = 3,                 // Weights on vertices and adjacencies are given
               numflag = 0,                 // numbering of verts starts by 0 (C-Style)
               ncon    = 1,                 // one weight per vertex
               nparts  = DynamicDataInterfaceCL::InfoProcs();   // number of subdomains (per proc one)
        float *tpwgts  = new float[nparts], // weight of partition
               itr     = quality_,          // how much an exchange costs
               ubvec   = GetGraph().ubvec;  // imbalace tolerance for eacht vertex weight
        int* options = new int[4];
    //     options[0]=1; options[1]=3; options[2]=15, options[3]=1;    // display times within parmetis
        options[0]=0;                                               // no options for parmetis

        if ((meth_ == KWay) || (meth_ == Adaptive))
            std::fill(tpwgts, tpwgts+nparts, 1./(float)nparts);
        if (GetGraph().part == 0)
            GetGraph().part = new idxtype[GetGraph().myVerts];
        if (meth_ == Identity)
            std::fill(GetGraph().part, GetGraph().part+GetGraph().myVerts, ProcCL::MyRank());

        MPI_Comm comm = MPI_COMM_WORLD;
        //Depending on the method choosen Adaptive, KWay or Identity the apropiate parmetis function is called
        switch (meth_)
        {       
            case Adaptive:
            ParMETIS_V3_AdaptiveRepart(
                    GetGraph().vtxdist, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().vwgt,
                    GetGraph().adjwgt, &wgtflag, &numflag, &ncon, &nparts,tpwgts,
                    &ubvec, &itr, options, &GetGraph().edgecut, GetGraph().part, &comm);break;
            case KWay:
            ParMETIS_V3_PartKway(
            GetGraph().vtxdist, GetGraph().xadj, GetGraph().adjncy, GetGraph().vwgt, GetGraph().adjwgt,
            &wgtflag, &numflag, &ncon, &nparts,
            tpwgts, &ubvec,0 /*options*/,
            &GetGraph().edgecut, GetGraph().part, &comm);break;
            case Identity:break;
            case Recursive:break;
            case NoMig:break;
        }

        if ((meth_ == KWay) || (meth_ == Adaptive))
        delete[] tpwgts;
        delete[] options;
}
//end of the implementations of the methods in the derived class ParMetisCL (parent PartitionerCL)
//================================================================================================

#ifdef _ZOLTAN
/// \brief Constructor of the derived class ZoltanCL (parent class PartitionerCL)
ZoltanCL::ZoltanCL(PartMethod meth):PartitionerCL() //Constructor of the Zoltan partitioner
    {
        meth_ = meth;
        zz = Zoltan_Create(ProcCL::GetComm());
    }

/// \brief DEstructor of the derived class ZoltanCL (parent class PartitionerCL)
ZoltanCL::~ZoltanCL()
{
    Zoltan_Destroy(&zz);
}

/** Querry function required by Zoltan
 \param data      pointer to the GraphST attribute
        ierr      zoltan error type
 */
int ZoltanCL::get_number_of_vertices(void *data, int *ierr)
{
    GraphST *graph = (GraphST *)data;
    *ierr = ZOLTAN_OK;
     return graph->myVerts;
}

/** Querry function required by Zoltan
 \param data        pointer to the GraphST attribute
        globalID    global id of the actual vertex
        localID     local id of the actual vertex
        obj_wgts    weights of the vertices
        ierr        zoltan error type
 */
void ZoltanCL::get_vertex_list(void *data, int, int, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int, float *obj_wgts, int *ierr)
{
    int i;
    GraphST *graph = (GraphST *)data;
    *ierr = ZOLTAN_OK;
    // Return the IDs of our verticies, and the weights.
    for (i=0; i<graph->myVerts; i++)
    {
      globalID[i] = i+ graph->myfirstVert;
      localID[i] = i;
      obj_wgts[i] = graph->vwgt[i];
    }
}

/** Querry function required by Zoltan
 \param data        pointer to the GraphST attribute
        sizeGID     dimension size of the global id of the actual vertex
        sizeLID     dimension size of the local id of the actual vertex
        num_obj     number of vertices
        localID     local id of the actual vertex
        num_edges   number of edges in the graph
        nbor_GID    global id of the neighboring vertices
        nbor_Proc   processor ot which the neighbor is assigned
        wgt_dim     dimension size of the weights of the edges
        ewgts       weights of the edges
        ierr        Zoltan error type
 */
void ZoltanCL::get_edge_list(void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR , ZOLTAN_ID_PTR localID, int *num_edges, ZOLTAN_ID_PTR nborGID, int *nborProc, int wgt_dim, float *ewgts, int *ierr)//returns the list of edges of each processor
{
    int i, j, from, to;
    int *nextNbor, *nextProc;
    float *nextWgt;

    GraphST *graph = (GraphST *)data;
    *ierr = ZOLTAN_OK;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->myVerts)|| (wgt_dim != 1))
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

      nextNbor = (int *)nborGID;
      nextProc = nborProc;
      nextWgt = ewgts;

      for (i=0; i < num_obj; i++)
      {
        to = graph->xadj[localID[i]+1];
        from = graph->xadj[localID[i]];
        if ((to - from) != num_edges[i])
        {
            *ierr = ZOLTAN_FATAL;
            return;
        }
        for (j=from; j < to; j++)
        {
          *nextNbor++ = graph->adjncy[j];//?????
          *nextProc++ = graph->GetProc(graph->adjncy[j]);
          *nextWgt++ = graph->adjwgt[j];
        }
      }
}

/** Querry function required by Zoltan
 \param data        pointer to the GraphST attribute
        sizeGID     dimension size of the global id of the actual vertex
        sizeLID     dimension size of the local id of the actual vertex
        num_obj     number of vertices
        num_Edges   number of edges
        ierr        Zoltan error type
 */
void ZoltanCL::get_num_edges_list(void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *numEdges, int *ierr)// returns the total number of edges
{
    int i;
      GraphST *graph = (GraphST *)data;
      if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->myVerts))
      {
        *ierr = ZOLTAN_FATAL;
        return;
      }
      for (i=0;  i < num_obj ; i++)
        numEdges[i] = graph->xadj[i+1] - graph->xadj[i];
      *ierr = ZOLTAN_OK;
}

/// \brief Method that sets some parameters inside Zoltan so that the communication between drops and the partitioner is correct
void ZoltanCL::CreateGraph()
{

    Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");


    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &GetGraph());
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list,  &GetGraph());
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list,  &GetGraph());
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list,  &GetGraph());
}

/// \brief Serial partitioning
void ZoltanCL::PartGraphSer( int)
{
    Comment("- Start calculate LoadBalanace with Zoltan-Partition"<<std::endl, DebugLoadBalC);

    if (GetGraph().myVerts != GetGraph().vtxdist[DynamicDataInterfaceCL::InfoProcs()])
        Comment("LoadBalCL: PartGraphSer: This procedure is called by a proc, that does not have all nodes!"<<std::endl, DebugLoadBalC);
    int i;
    int changes ; // Set to 1 or .TRUE. if the decomposition was changed by the load-balancing method; 0 or .FALSE. otherwise.
    int numGidEntries; //the number of array entries used to describe a single global ID
    int numLidEntries; //the number of array entries used to describe a single local ID
    int numImport;
    ZOLTAN_ID_PTR importGlobalGids;
    ZOLTAN_ID_PTR importLocalGids;
    int *importProcs;
    int *importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalGids;
    ZOLTAN_ID_PTR exportLocalGids;
    int *exportProcs;
    int *exportToPart;
    int rc;

    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");//Partition "from scratch," not taking into account the current data distribution;
                                                //this option is recommended for static load balancing
    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
          &changes,     /* 1 if partitioning was changed, 0 otherwise */
          &numGidEntries,  /* Number of integers used for a global ID */
          &numLidEntries,  /* Number of integers used for a local ID */
          &numImport,      /* Number of vertices to be sent to me */
          &importGlobalGids,  /* Global IDs of vertices to be sent to me */
          &importLocalGids,   /* Local IDs of vertices to be sent to me */
          &importProcs,    /* Process rank for source of each incoming vertex */
          &importToPart,   /* New partition for each incoming vertex */
          &numExport,      /* Number of vertices I must send to other processes*/
          &exportGlobalGids,  /* Global IDs of the vertices I must send */
          &exportLocalGids,   /* Local IDs of the vertices I must send */
          &exportProcs,    /* Process to which I send each of the vertices */
          &exportToPart);


    if (rc != ZOLTAN_OK)
    {
      printf("sorry...\n");
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
    }

    for (i=0; i < GetGraph().myVerts; i++)
        GetGraph().part[i] = ProcCL::MyRank();

    for (i=0; i < numExport; i++)
        GetGraph().part[exportLocalGids[i]] = exportToPart[i];

}

/// \brief Parallel partitioning
void ZoltanCL::PartGraphPar()
{
    Comment("- Start calculate LoadBalanace with Zoltan-Repartition"<<std::endl,DebugLoadBalC);

    Assert( GetGraph().xadj && GetGraph().adjncy_ && GetGraph().vwgt && GetGraph().adjwgt,
            DROPSErrCL("LoadBalCL::PartGraphPar: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
            DebugLoadBalC);
    int i;
    int changes; // Set to 1 or .TRUE. if the decomposition was changed by the load-balancing method; 0 or .FALSE. otherwise.
    int numGidEntries; //the number of array entries used to describe a single global ID
    int numLidEntries; //the number of array entries used to describe a single local ID
    int numImport;
    ZOLTAN_ID_PTR importGlobalGids;
    ZOLTAN_ID_PTR importLocalGids;
    int *importProcs;
    int *importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalGids;
    ZOLTAN_ID_PTR exportLocalGids;
    int *exportProcs;
    int *exportToPart;
    int rc;

    Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");//Partition but take into account current data distribution to keep data migration low;
                                                 //this option is recommended for dynamic load balancing

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
          &changes,        /* 1 if partitioning was changed, 0 otherwise */
          &numGidEntries,  /* Number of integers used for a global ID */
          &numLidEntries,  /* Number of integers used for a local ID */
          &numImport,      /* Number of vertices to be sent to me */
          &importGlobalGids,  /* Global IDs of vertices to be sent to me */
          &importLocalGids,   /* Local IDs of vertices to be sent to me */
          &importProcs,    /* Process rank for source of each incoming vertex */
          &importToPart,   /* New partition for each incoming vertex */
          &numExport,      /* Number of vertices I must send to other processes*/
          &exportGlobalGids,  /* Global IDs of the vertices I must send */
          &exportLocalGids,   /* Local IDs of the vertices I must send */
          &exportProcs,    /* Process to which I send each of the vertices */
          &exportToPart);

    if (rc != ZOLTAN_OK)
    {
      printf("sorry...\n");
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
    }
    for (i=0; i < GetGraph().myVerts; i++)
        GetGraph().part[i] = ProcCL::MyRank();

    for (i=0; i < numExport; i++)
        GetGraph().part[exportLocalGids[i]] = exportToPart[i];
}
# endif
//end of the implementations of the methods in the derived class ZoltanCL (parent PartitionerCL)
//================================================================================================
}//end of namespace
