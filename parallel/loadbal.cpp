/// \file loadbal.cpp
/// \brief Loadbalancing of tetrahedal multigrids
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier, Timo Henrich

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

#include "parallel/loadbal.h"
#include "parallel/parallel.h"
#include "num/interfacePatch.h"
#include <iomanip>

namespace DROPS{

/****************************************************************************
* L O A D  B A L  C L A S S                                                 *
****************************************************************************/
// static initialisation
Uint     LoadBalCL::TriangLevel_ = 0;
IFT      LoadBalCL::FaceIF_      = 0;
idxtype* LoadBalCL::myfirstVert_ = 0;


/****************************************************************************
* W R A P P E R                                                             *
*****************************************************************************
* Converts C++ functions to C functions, so DDD can use them correctly      *
****************************************************************************/
extern "C" int HandlerScatterC( OBJT o, void *d) {return LoadBalCL::HandlerScatter(o,d);}
extern "C" int HandlerGatherC ( OBJT o, void *d) {return LoadBalCL::HandlerGather(o,d);}


/// \brief Constructor
LoadBalCL::LoadBalCL(MultiGridCL& mg, int partitioner, int TriLevel, PartMethod meth, int weightFct) : idx_(), lset_(0), weightFct_(weightFct)
/// \param mg          Reference on the multigrid
/// \param partitioner Choose a partitioner: 1 - Metis, 2 - Zoltan, 3 - Scotch
/// \param TriLevel    level that should be balanced
/// \param meth        Type of method used for partitioning
/// \param weightFct   Which information for weighting the dual reduced graph should be used. This parameter
///                    is a binary representation (b2 b1 b0) of the three methods:
///                    b0 : use information about children
///                    b1 : use information about unknowns
///                    b2 : use information about intersected subs
/// \todo (of) estimate good parameter ubvec for parmetis
{
    partitioner_ = PartitionerCL::newPartitioner( Partitioner(partitioner), 1.05, meth);   // Pointer to the partitioner class
    mg_=&mg;
    if (FaceIF_==0)
        InitIF();

    if (TriLevel<0)
        TriangLevel_=mg_->GetLastLevel();
    else
        TriangLevel_=TriLevel;
    partitioner_->GetGraph().movedMultiNodes=0;
}


/// \brief Destructor
LoadBalCL::~LoadBalCL()
{
    DeleteGraph();
}


/// \brief Init the DDD interface for communicate adjacenies between processors
void LoadBalCL::InitIF()
/** \todo (of) Falls nicht nur das letzte Triangulierungslevel balanciert werden
    soll, mussen auch LB-Nummern auf Ghosts versendet werden!
*/
{
    TypeT  O[8];
    PrioT  A[3], B[3];
    /* O is an array of object types, A and B are arrays of priorities */

    /* interface of faces */
    O[0] = FaceCL::GetType();
    A[0]= PrioHasUnk;   B[0]= PrioHasUnk;
    A[1]= PrioMaster;   B[1]= PrioMaster;
    A[2]= PrioGhost;    B[2]= PrioGhost;

    // Da auf der letzten Triangulierung keine Ghosts vorkommen, sind die Ghost-
    // Subsimplices fuer die Kommunikation der Graph-Verts bedeutungslos. Fuer
    // die Lastverteilung von anderen Triangulierungen oder die Lastverteilung
    // fuer Mehrgitterverfahren muessen die Ghosts evtl. beruecksichtigt werden...

    FaceIF_ = DynamicDataInterfaceCL::IFDefine(1, O, 2, A, 2, B);
    DynamicDataInterfaceCL::IFSetName( FaceIF_, (char*)"Face-Interface for LB");
}


/// \brief Communicate adjacencies between processors
/// \pre pointer to myfirstvert_ has to be set!
void LoadBalCL::CommunicateAdjacency()
{
    DynamicDataInterfaceCL::IFExchange( FaceIF_, sizeof(idxtype), &HandlerGatherC, &HandlerScatterC);
    myfirstVert_=0;
}


/// \brief Gather adjacency information at processor-boundary
int LoadBalCL::HandlerGather( OBJT obj, void *buf)
/** This function is called by DDD within the CommunicateAdjacency procedure
    \param obj a pointer to a face on a processor boundary
    \param buf buffer, where to place the load balancing number
*/
{
    FaceCL* const fp = ddd_cast<FaceCL*>( obj);            // transform DDD-Object to a FaceCL-pointer
    idxtype* sendbuf = static_cast<idxtype*>(buf);         // transform buf into a storage for idxtype

    if (!fp->IsInTriang( TriangLevel_) )                   // this is not correct
        return 1;

    const TetraCL* tp= fp->GetSomeTetra();                 // get a arbitrary tetraeder

    if (tp->HasGhost() )                                   // if this tetra has ghost, this is not the right one
        tp= fp->GetNeighborTetra( tp);                     // get the next one

    Assert(tp->HasLbNr(),
           DROPSErrCL("LoadBalCL::HandlerGather: Tetra without loadbalance-number"),
           DebugLoadBalC);

    if (tp->HasLbNr())                                     // if the neighbor tetra has a LoadBalance Nr (children are numbered too!)
        *sendbuf= tp->GetLbNr() + *myfirstVert_ ;// put the number into the send buffer

    return 0;
}


/// \brief Gather adjacency information at processor-boundary
int LoadBalCL::HandlerScatter( OBJT obj, void *buf)
/**  Store the recieved load balance number in the FaceCL
    \param obj a pointer to a face on a processor boundary
    \param buf recieved load balancing number
*/
{
    FaceCL* const fp= ddd_cast<FaceCL*>( obj);                 // transform DDD-Object to a FaceCL-pointer
    const idxtype* const recvbuf= static_cast<idxtype*>(buf);  // transform recieved data to a idxtype-pointe
    fp->SetLbNeigh(*recvbuf);                                  // store the number correct
    return 0;
}


/// \brief Iterator to the first tetra in LoadBalSet
LbIteratorCL LoadBalCL::GetLbTetraBegin(int TriLevel) const
/** Take tetra out of the coarsest level of the multigrid an check if it is in
    the loadbalancing set
    \param TriLevel if greater than zero the level to be partitioned, else
     the finest level
*/
{
    // if function is called with default TriLevel (-1), set TriLevel to last level
    if ( TriLevel<0)
        TriLevel = mg_->GetLastLevel();

    // Compute smalles number of non-empty level
    Uint i=0;
    while (mg_->GetTetras().IsLevelEmpty(i) && i<mg_->GetTetras().GetNumLevel()-1)
        ++i;

    // Get first iterator of this level
    LbIteratorCL ret(mg_, mg_->GetTetrasBegin(i), i, TriLevel);

    // if the tetra, this iterator points to isn't in the LbSet, get the first
    // one, that is in the LbSet.
    return ( !mg_->GetTetras().IsLevelEmpty(i) && ret.IsInLbSet())? ret : ++ret;
}


/// \brief iterator to the last tetra in LoadBalSet
LbIteratorCL LoadBalCL::GetLbTetraEnd( int Level) const
{
     if (Level<0)
         Level= mg_->GetLastLevel();                                // if last TriangLevel is wished, get it
     return LbIteratorCL( mg_, mg_->GetTetrasEnd(Level), Level, Level);
}


/// \brief Create the numbering of the multinodes in the LoadBalSet
void LoadBalCL::CreateNumbering(int TriLevel)
/** Assign each tetrahedra in the loadbalancing set a loadbalancing number
    \param TriLevel level of the triangulation that should be load balanced
*/
{
    idxtype counter=0;                      // Nodes are first numbered on all procs starting by 0. Offset is added later

    LbIteratorCL  it (GetLbTetraBegin()),   // first tetra in LbSet
                 end ( GetLbTetraEnd() );   // last tetra in LbSet

    // for all tetras in the LbSet
    for (; it!=end; ++it)
    {
        it->SetLbNr(counter);               // set number to counter

        // assign the same vertex number (=loadbalnr) to all children
        for (TetraCL::ChildPIterator ch(it->GetChildBegin()), chend (it->GetChildEnd()); ch!=chend; ++ch)
            if ((*ch)->IsInTriang( TriLevel) )
                (*ch)->SetLbNr(counter);

        // increase the counter, to number the next tetra(s)
        counter++;
    }

    // number of vertices is determined by the counter
    partitioner_->GetGraph().myVerts=counter;
}


/// \brief Remove all loadbalancing-numbers from the tetrahedra
void LoadBalCL::RemoveLbNr()
{
    MultiGridCL::TetraIterator  sit(mg_->GetAllTetraBegin()), end(mg_->GetAllTetraEnd());
    for (; sit!=end; ++sit)
    {
        sit->DelLbNr();
    }
}

/** Put the weight of a vertex in the array vwgt of the graph structure at
    the index wgtpos. The weight is computed by the number of children or
    set to 1 if no children exists. Afterwards, this wgtpos is increased by 1.
    \param t parent tetrahedron
    \param wgtpos  in: where to put the weight in the array vwgt,
                  out: where to put the next weight in the array vwgt
*/
void LoadBalCL::GetWeightRef( const TetraCL& t, size_t& wgtpos)
{
    if (t.IsUnrefined()){
        GetPartitioner()->GetGraph().vwgt[wgtpos++]= 1;
    }
    else{
        GetPartitioner()->GetGraph().vwgt[wgtpos++]= t.GetRefData().ChildNum;
    }
}

/** Helper function for determining the weights by unknowns
    \param t tetrahedron of the finest level
    \param list list of unknowns
*/
void LoadBalCL::UnkOnSingleTetra( const TetraCL& t, LoadBalCL::UnkWghtListT& list) const
{
    for ( size_t i=0; i<idx_.size(); ++i){
        for ( TetraCL::const_VertexPIterator it=t.GetVertBegin(); it!=t.GetVertEnd(); ++it)
            UnkOnSimplex( **it, list, idx_[i]);
        for ( TetraCL::const_EdgePIterator it=t.GetEdgesBegin(); it!=t.GetEdgesEnd(); ++it)
            UnkOnSimplex( **it, list, idx_[i]);
    }
}

/** Determine weight of a vertx by counting the number of unknowns
    \param t parent tetrahedron
    \param wgtpos  in: where to put the weight in the array vwgt,
                  out: where to put the next weight in the array vwgt
*/
void LoadBalCL::GetWeightUnk( const TetraCL& t, size_t& wgtpos)
{
    UnkWghtListT unklist;
    // Collect information
    if ( t.IsUnrefined()){
        UnkOnSingleTetra(t, unklist);
    }
    else{
        for ( TetraCL::const_ChildPIterator it=t.GetChildBegin(); it!=t.GetChildEnd(); ++it)
            UnkOnSingleTetra( **it, unklist);
    }
    // Count number of dof on tetra
    int numUnks=0;
    for ( UnkWghtListT::const_iterator it= unklist.begin(); it!=unklist.end(); ++it){
        numUnks += it->numUnk;
    }
    GetPartitioner()->GetGraph().vwgt[wgtpos++]= numUnks;
}


/** Put the weight of a vertex in the array vwgt of the graph structure at
    the index wgtpos. The weight consists of two constrains. The first one is
    computed by GetWeightRef and the second one by the number of subtetrahedra
    intersected by the interface of the level set function. Afterwards, this
    wgtpos is increased by 2.
    \param t parent tetrahedron
    \param wgtpos  in: where to put the weight in the array vwgt,
                  out: where to put the next weight in the array vwgt
*/
void LoadBalCL::GetWeightLset( const TetraCL& t, size_t& wgtpos)
{
    InterfacePatchCL patch;                              // check for intersection
    if ( t.IsUnrefined()){
        if (CheckForLsetUnk(t)){
            patch.Init( t, *lset_, *lsetbnd_);
            GetPartitioner()->GetGraph().vwgt[wgtpos++]= patch.GetNumIntersectedSubTetras();
        }
        else{
            GetPartitioner()->GetGraph().vwgt[wgtpos++]= 0;
        }
    }
    else{
        int numItersectedSub=0;
        for ( TetraCL::const_ChildPIterator it=t.GetChildBegin(); it!=t.GetChildEnd(); ++it){
            if (CheckForLsetUnk(**it)){
                patch.Init( **it, *lset_, *lsetbnd_);
                numItersectedSub += patch.GetNumIntersectedSubTetras();
            }
        }
        GetPartitioner()->GetGraph().vwgt[wgtpos++]= numItersectedSub;
    }
}

/// \brief Estimate adjacencies of an unrefined tetra
void LoadBalCL::AdjUnrefined( TetraCL& t, int& edgecount, size_t& vwgpos)
/** Put the adjacenzies into adjncy_ and estimate the weight of the node, assoziated with the tetra
    if the tetraeder is unrefined
    \param t the tetraeder
    \param edgecount IN/OUT: smallest unused edgenumber
    \param vwgpos  IN: where to put the vertex weight in the array vwgt of the graph structure
                  OUT: next position where to put the weight of a vertex
*/
{
    for (Uint face= 0; face<NumFacesC; ++face){                             // Adjacences are faces
        if (!t.IsBndSeg(face)){                                             // If this face belongs to a domain boundary make nothing
            if (t.GetFace(face)->IsOnProcBnd() ){                           // If this face belongs to a boundary between procs
                partitioner_->GetGraph().adjwgt[edgecount]= 1;                                      // set edgeweight to one
                partitioner_->GetGraph().adjncy[edgecount++]= t.GetFace(face)->GetLbNeigh();        // put neighbor into adjacenz list
            }
            else{                                                           // Neighbor tetra is on the same proc
                const TetraCL* const neigh= t.GetNeighbor( face);           // get a pointer to this tetra
                Assert(neigh->HasLbNr(),
                       DROPSErrCL("LoadBalCL::AdjUnrefined: Tetra without loadbalance-number"),
                       DebugLoadBalC);
                if (neigh->HasLbNr()){                                      // test if the neighbor is in the LoadBalSet
                    partitioner_->GetGraph().adjwgt[edgecount]= 1;                                  // set edgeweight to one
                    partitioner_->GetGraph().adjncy[edgecount++]= neigh->GetLbNr()+ partitioner_->GetGraph().myfirstVert;    // put neighbor into the adjacenz list
                }
            }
        }
    }

    // Compute weight of the vertex
    if ( weightFct_&1)
        GetWeightRef( t, vwgpos);
    if ( weightFct_&2){
        if ( idx_.empty())
            std::cout << "No information about unknowns is known to LoadBalCL.\nNo information of unknowns is considered for load balancing" << std::endl;
        else
            GetWeightUnk(t, vwgpos);
    }
    if ( weightFct_&4){
        if ( lset_==0)
            std::cout << "No information about level set is known to LoadBalCL.\nNo information of level set is considered for load balancing" << std::endl;
        else
            GetWeightLset( t, vwgpos);
    }
}


/// \brief Estimate adjacencies of an refined tetra
void LoadBalCL::AdjRefined( TetraCL& t, int& edgecount, size_t& vwgpos)
/** Put the adjacenzies into adjncy_ and estimate the weight of the node, assoziated with the tetra
    if the tetraeder is unrefined
    \param t the tetraeder
    \param edgecount IN/OUT: smallest unused edgenumber
    \param vwgpos  IN: where to put the vertex weight in the array vwgt of the graph structure
                  OUT: next position where to put the weight of a vertex
*/
{
    GraphEdgeCT thisAdj;                                                // store adjacencies of this tetraeder-node

    TetraCL::ChildPIterator ch    (t.GetChildBegin()),                  // iterator through the children
                            chend ( t.GetChildEnd() );                  // last child

    for (;ch!=chend; ++ch){
        if ((*ch)->IsInTriang(TriangLevel_)){                           // if this child is the triangulation, that should be balanced
            for (Uint face= 0; face<NumFacesC; ++face){                 // iterate over all faces
                if (!(*ch)->IsBndSeg(face)){                            // if this face belong to the domain boundary do nothing
                    if ((*ch)->GetFace(face)->IsOnProcBnd() )           // neighbor is found on another proc
                        ++thisAdj[(*ch)->GetFace(face)->GetLbNeigh()];  // put neighbor into the adjacency map and increase the weight by one
                    else {                                               // tetra has neighbor tetra on this proc
                        const TetraCL* const neigh= (*ch)->GetNeighInTriang( face, TriangLevel_);
                        Assert(neigh->HasLbNr(),
                               DROPSErrCL("LoadBalCL::AdjRefined: Tetra without loadbalance-number"),
                               DebugLoadBalC);

                        if (neigh->HasLbNr() )                          // test if this neighbor has a number
                            ++thisAdj[neigh->GetLbNr() +  partitioner_->GetGraph().myfirstVert]; // put neighbor into the adjazenz map and increase the weight by one
                    }
                }
            }
        }
    }

    thisAdj.erase( t.GetLbNr() + partitioner_->GetGraph().myfirstVert);                         // delete my one number (set into this adjacenclist by children of t, which have the same number as t)

    GraphEdgeCT::iterator  it = thisAdj.begin(),                        // start at the first adjacenz
                          end = thisAdj.end();                          // end at the last adjacenz
    for (; it!=end; ++it, ++edgecount){                                 // and put all adjacencies into the list for parMETIS
        partitioner_->GetGraph().adjncy[edgecount]= it->first;
        partitioner_->GetGraph().adjwgt[edgecount]= it->second;
    }

    // Compute weight of the vertex
    if ( weightFct_&1)
        GetWeightRef( t, vwgpos);
    if ( weightFct_&2){
        if ( idx_.empty())
            std::cout << "No information about unknowns is known to LoadBalCL.\nNo information of unknowns is considered for load balancing" << std::endl;
        else
            GetWeightUnk(t, vwgpos);
    }
    if ( weightFct_&4){
        if ( lset_==0)
            std::cout << "No information about level set is known to LoadBalCL.\nNo information of level set is considered for load balancing" << std::endl;
        else
            GetWeightLset( t, vwgpos);
    }
}


/// \brief Compute the number of adjacencies on this proc
Uint LoadBalCL::EstimateAdj()
/// \return number of adjacencies
{
    GraphNeighborCT adj;                                                // Set of neighbors of actual multinode
    Uint adjCounter = 0;                                                // counter of all adjacencies
    LbIteratorCL begin= GetLbTetraBegin(), end= GetLbTetraEnd();        // iterate through the LoadBalSet

    for (LbIteratorCL it=begin; it!=end; ++it)
    {
        if (it->IsUnrefined() ){                                        // Unrefined => (number of faces - faces of domain boundary)
            for (Uint face= 0; face<NumFacesC; ++face)
                if (!it->IsBndSeg(face))
                    ++adjCounter;
        }
        else{

            adj.clear();                                                // clear the List

            for (TetraCL::ChildPIterator ch (it->GetChildBegin()), chend (it->GetChildEnd()); ch!=chend; ++ch){
                if ((*ch)->IsInTriang(TriangLevel_) ){
                    for (Uint face= 0; face<NumFacesC; ++face){
                        const FaceCL* const fp= (*ch)->GetFace(face);
                        if (fp->IsOnBoundary())
                            continue;
                        if (fp->IsOnProcBnd() )
                            adj.insert( fp->GetLbNeigh() );
                        else
                        { // tetra has neighbor tetra
                            const TetraCL* const neigh ( (*ch)->GetNeighInTriang( face, TriangLevel_) );

                            Assert(neigh,
                                   DROPSErrCL("LoadBalCL::EstimateAdj: No neighbor found!"),
                                   DebugLoadBalC);

                            Assert(neigh->HasLbNr(),
                                   DROPSErrCL("LoadBalCL::EstimateAdj:  Tetra without loadbalance-number"),
                                   DebugLoadBalC);

                            if (neigh->HasLbNr() )
                                adj.insert( neigh->GetLbNr() +  partitioner_->GetGraph().myfirstVert);
                        }
                    }
                }
            }
            // no adjacencies with members of my multinode...
            adj.erase( it->GetLbNr() +  partitioner_->GetGraph().myfirstVert);
            adjCounter+= adj.size();
        }
    }

    partitioner_->GetGraph().myAdjs = adjCounter;
    return adjCounter;
}


/// \brief Create the dual reduced Graph for ParMETIS on the last triangulation level
void LoadBalCL::CreateDualRedGraph(bool geom)
/** Set up the graph for (Par)Metis.
    \param geom use of geometric information
    \todo (of) Create Graph of arbitrary triangulation level
*/
{
    Comment("- Setting up dual, reduced graph with " << partitioner_->GetGraph().myAdjs<< " edges " << std::endl, DebugLoadBalC);

    partitioner_->GetGraph().geom = geom;                                       // sets the attribute geom of the partitioner_
    TriangLevel_= mg_->GetLastLevel();

    CreateNumbering();                                                          // Create the numbers of the multinodes and set the number of my verts

    // calculate the vtxdist-Array
    partitioner_->GetGraph().ResizeVtxDist();
    IndexArray vtx_rcv= new int[ProcCL::Size()];                                // receive buffer for number of nodes on other procs
    ProcCL::Gather( (int)partitioner_->GetGraph().myVerts, vtx_rcv, -1);        // communicate the number of nodes

    partitioner_->GetGraph().vtxdist[0]= 0;                                     // proc 0 starts with node 0
    for (int i=0; i<ProcCL::Size(); ++i)
        partitioner_->GetGraph().vtxdist[i+1]= partitioner_->GetGraph().vtxdist[i] + vtx_rcv[i];    // proc i+1 starts with node \sum_{j=0}^i nodesOnProc(i)
    partitioner_->GetGraph().myfirstVert = partitioner_->GetGraph().vtxdist[ProcCL::MyRank()];     // vtxdist_[i] is starting node of proc i

    delete[] vtx_rcv;


     // Compute Adjacencies
    myfirstVert_=&(partitioner_->GetGraph().myfirstVert);
    CommunicateAdjacency();                                     // create adjacencies over proc boundaries
    Uint numadj= EstimateAdj();                                 // compute the number of adjacencies on this proc

     // Number of conditions per vertex
    int ncon=0;
    if ( weightFct_&1 || weightFct_==0) ++ncon;
    if ( weightFct_&2 && !idx_.empty()) ++ncon;
    if ( weightFct_&4 && lset_!=0)      ++ncon;
     // Allocate space for the Arrays
    partitioner_->GetGraph().Resize(numadj, partitioner_->GetGraph().myVerts, partitioner_->GetGraph().geom,  ncon);

     // put all nodes and adjacencies into the lists
    LbIteratorCL begin = GetLbTetraBegin(),
                   end = GetLbTetraEnd();
    idxtype vertCount=0,                                        // number of actual node
            edgeCount=0;                                        // number of actual edge
    size_t vwgtcount=0;                                         // counter where to put the weight of a vertex
    for ( LbIteratorCL it= begin; it!=end; ++it)                // iterate through the multinodes of the LoadBalSet
    {
        partitioner_->GetGraph().xadj[vertCount] = edgeCount;   // remember, where node vertCount writes its neighbors
        if (it->IsUnrefined() )
            AdjUnrefined( *it, edgeCount, vwgtcount);           // compute adjacencies, number of adjacencies and node-weight for unrefined tetra
        else
            AdjRefined( *it, edgeCount, vwgtcount);             // compute adjacencies, number of adjacencies and node-weight for refined tetra
        if (partitioner_->GetGraph().geom)                      // if with geom, compute the barycenter of the tetra
        {
            const Point3DCL coord= GetBaryCenter( *it);
            partitioner_->GetGraph().xyz[3*vertCount]  = coord[0];
            partitioner_->GetGraph().xyz[3*vertCount+1]= coord[1];
            partitioner_->GetGraph().xyz[3*vertCount+2]= coord[2];
        }

        Assert( edgeCount<=(int)numadj,
                DROPSErrCL( "LoadBalanceCL: CreateDualRedGraph: Too much adjacencies!"),
                DebugLoadBalC);

        ++vertCount;                                            // go to the next node
    }

    partitioner_->GetGraph().xadj[vertCount]= edgeCount;
    Assert( vertCount==partitioner_->GetGraph().myVerts,
            DROPSErrCL("LoadBalaceCL: CreateDualRedGraph: number of vertices is not correct!"),
            DebugLoadBalC);
    Assert( (partitioner_->GetGraph().ncon)*vertCount== (int)vwgtcount,
            DROPSErrCL("LoadBalaceCL: CreateDualRedGraph: Number of vertex weights does not match"),
            DebugLoadBalC);

    partitioner_->CreateGraph(); //Implement for the other partitioners(not necessary for metis)
}


/// \brief Remove information about the graph
void LoadBalCL::DeleteGraph()
/** Free all allocated arrays */
{
    partitioner_->GetGraph().Clear();
    partitioner_->GetGraph().movedMultiNodes=0;
}

/// \brief Do migration of the tetrahedra
void LoadBalCL::Migrate()
/** Iteratate over the tetrahedra that are in the loadbalancing set and tell
    DDD on which processer these tetras belong. For a detailed description see
    Diploma thesis of Sven Gross.
*/
{
#if DROPSDebugC
	GIDT observe1 = 207360, observe2=0, observe3=0;
#endif
    if (ProcCL::MyRank()==0)
        Comment("- Start Migrating"<<std::endl, DebugLoadBalC);

    Uint me = ProcCL::MyRank();
    partitioner_->GetGraph().movedMultiNodes=0;

    PROCT dest;
    for (LbIteratorCL it= GetLbTetraBegin(), end= GetLbTetraEnd(); it!=end; ++it)
    {
        dest =  static_cast<PROCT>(partitioner_->GetGraph().part[it->GetLbNr()]);

        if (dest==me) continue;
        partitioner_->GetGraph().movedMultiNodes++;
#if DROPSDebugC
        if ( it->GetGID()==observe1 || it->GetGID()==observe2 || it->GetGID()==observe3)
            std::cout << "["<<ProcCL::MyRank()<<"] ===> Transfer des Tetras mit GID "<<it->GetGID() << " nach " << dest << " als ";
#endif

        if (it->IsUnrefined() )
        { // E1-Xfer
            ParMultiGridCL::TXfer( *it, dest, PrioMaster, true);
#if DROPSDebugC
            if ( it->GetGID()==observe1 || it->GetGID()==observe2 || it->GetGID()==observe3)
                std::cout << "E1-Xfer mit delete =1 und PrioMaster" << std::endl;
#endif
        }
        else
        { // E2-Xfer
        	PrioT asPrio=PrioGhost;

            for( int* proclist= DynamicDataInterfaceExtraCL::InfoProcList( it->GetHdr() ); *proclist!=-1; proclist+= 2){
                if (*proclist==dest){
                    if (proclist[1]>=PrioMaster){
                        asPrio=PrioMaster;
                    }
                }
            }
#if DROPSDebugC
            if ( it->GetGID()==observe1 || it->GetGID()==observe2 || it->GetGID()==observe3)
                std::cout << "E2-Xfer mit delete ="<< (it->GetPrio()==PrioGhost)
                        << " und Prio"<<(asPrio==PrioMaster?"Master":"Ghost")<<" Unrefined=" << it->IsUnrefined()
                        << std::endl;
#endif

            ParMultiGridCL::TXfer( *it, dest, asPrio, it->GetPrio()==PrioGhost);
        }

        it->DelLbNr();

        if (!it->IsUnrefined() )
        {
            for (TetraCL::ChildPIterator ch(it->GetChildBegin()), chend(it->GetChildEnd()); ch!=chend; ++ch)
            {
                if ((*ch)->IsUnrefined() || (*ch)->HasGhost() )
                { // M1-Xfer
                    ParMultiGridCL::TXfer( **ch, dest, PrioMaster, true);
#if DROPSDebugC
                    if ( it->GetGID()==observe1 || (*ch)->GetGID()==observe1 || (*ch)->GetGID()==observe2 || (*ch)->GetGID()==observe3)
                        std::cout << "["<<ProcCL::MyRank()<<"]===> Transfer des Tetras mit GID "<< (*ch)->GetGID()
                                << " als Kind von " << it->GetGID() << " nach " << dest
                                << " als M1-Xfer mit delete =1 und PrioMaster" << std::endl;
#endif
                }
                else
                { // M2-Xfer
                    const bool E2Xfer= it.IsInLbSet( **ch) && partitioner_->GetGraph().part[(*ch)->GetLbNr()]!=static_cast<idxtype>(me);

                    ParMultiGridCL::TXfer( **ch, dest, PrioMaster, E2Xfer);
                    if (!E2Xfer)
                    {
                        ParMultiGridCL::PrioChange(*ch,PrioGhost);
                    }
#if DROPSDebugC
                    if ( (*ch)->GetGID()==observe1 || (*ch)->GetGID()==observe2 || (*ch)->GetGID()==observe3)
                        std::cout << "["<<ProcCL::MyRank()<<"]===> Transfer des Tetras mit GID "<< (*ch)->GetGID()
                                << " als Kind von " << it->GetGID() << " nach " << dest
                                << " als M2-Xfer mit delete =" << E2Xfer
                                << " und PrioMaster und ChangePrio to Prio"<< (E2Xfer?"Master":"Ghost") << std::endl;
#endif
                }
            }
        }
    }

    partitioner_->GetGraph().movedMultiNodes= ProcCL::GlobalSum(partitioner_->GetGraph().movedMultiNodes);
}


/// \brief Print the Graph on the ostream
void LoadBalCL::ShowGraph(std::ostream& os)
{
    IndexArray part     = partitioner_->GetGraph().part;
    IndexArray vwgt     = partitioner_->GetGraph().vwgt;
    IndexArray xadj     = partitioner_->GetGraph().xadj;
    IndexArray adjncy   = partitioner_->GetGraph().adjncy;
    int myVerts         = partitioner_->GetGraph().myVerts;
    int myAdjs          = partitioner_->GetGraph().myAdjs;
    int counter=0;
    if (xadj==0)
        os << "No Graph created!"<<std::endl;
    else
    {
        os << "Proc "<<ProcCL::MyRank()<<" stores "<<myVerts<<" of "<<GetNumAllVerts()
           << " with "<<myAdjs<<" adjacencies"<<std::endl
           << "The graph is:" <<std::endl
           << " Node | Weight | Dest | Neighbors"<<std::endl
           << "------+--------+------+--------------------------"<<std::endl;

        for (int i=0; i<myVerts; ++i)
        {
            os << std::setw(5) << i << " |" << std::setw(7) << vwgt[i] << " |";
            if (part!=0)
                os << std::setw(5) << part[i];
            else
                os << "  X  ";
            os <<" | ";
            for (int j=xadj[i]; j<xadj[i+1]; ++j)
            {
                os << std::setw(4) <<adjncy[j] << "  ";
                ++counter;
            }
            os<<std::endl;
        }
    }
}


/// \brief writes a serial graph onto the ostream
std::ostream& operator << (std::ostream &os ,const LoadBalCL &lb)
/** The format of the output is the inputfile for METIS. So the vertices starts by 1
    and the edges are not count twice, like in the CSR-format!
*/
{
    os << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    os << "%  Graphfile for METIS               %\n";
    os << "%  created by DROPS of the LoadBalCL %\n";
    os << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    os << lb.partitioner_->GetGraph().myVerts                                               // number of nodes
       << " " << (lb.partitioner_->GetGraph().myAdjs/2)                                     // number of "real" adjacencies
       << " 11" <<std::endl;                                        // there are weights onto the nodes and the adjacencies

    for (int i=0; i<lb.partitioner_->GetGraph().myVerts; ++i)                               // go over all nodes and write the information onto the stream
    {
        os << lb.partitioner_->GetGraph().vwgt[i] << " ";                                   // first number in the line is the weight of the node i
        for (int j=lb.partitioner_->GetGraph().xadj[i]; j<lb.partitioner_->GetGraph().xadj[i+1]; ++j)
            os << lb.partitioner_->GetGraph().adjncy[j]+1 << " " << lb.partitioner_->GetGraph().adjwgt[j] << " ";   // then the adjacence with adjacenc weight
        os << std::endl;
    }
    return os;
}


/// \brief Write graph, so that it can be read easily
void LoadBalCL::PrintGraphInfo(std::ostream& os) const
{   IndexArray vtxdist = partitioner_->GetGraph().vtxdist;
    IndexArray xadj    = partitioner_->GetGraph().xadj;
    IndexArray adjncy  = partitioner_->GetGraph().adjncy;
    IndexArray adjwgt  = partitioner_->GetGraph().adjwgt;
    IndexArray vwgt    = partitioner_->GetGraph().vwgt;

    const int me=ProcCL::MyRank(), size = ProcCL::Size();
    const int loc_num_verts=vtxdist[me+1]-vtxdist[me],
              loc_num_edges=xadj[loc_num_verts];

    // print general information
    os << size << ' ' << loc_num_verts << ' ' << loc_num_edges << std::endl;

    // print vtxdist
    for (int i=0; i<=size; ++i){
        os << vtxdist[i] << ' ';
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print xadj
    for (int i=0; i<=loc_num_verts; ++i){
        os << xadj[i] << ' ';
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print adjncy_
    for (int i=0; i<loc_num_edges; ++i){
        os << adjncy[i] << ' ';
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print vwgt
    for (int i=0; i<loc_num_verts; ++i){
        os << vwgt[i];
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print adjwgt
    for (int i=0; i<loc_num_edges; ++i){
         os << adjwgt[i];
        if (i%10==0) os << std::endl;
    }
    os << std::endl;
}


/****************************************************************************
* L O A D  B A L  H A N D L E R  C L A S S                                  *
****************************************************************************/

LoadBalHandlerCL::LoadBalHandlerCL(MultiGridCL& mg, Partitioner part, float /*ub*/) : mg_(&mg)
{
    lb_ = new LoadBalCL(mg, part);
    strategy_ = Adaptive;
    xferUnknowns_ = false;
    debugMode_    = false;
}

LoadBalHandlerCL::~LoadBalHandlerCL()
{
    if (lb_) delete lb_; lb_=0;
}

LoadBalHandlerCL::LoadBalHandlerCL(const MGBuilderCL &builder, int master, Partitioner part, PartMethod meth, bool geom, bool debug)
/// \param[in] builder Builder that creates a multigrid
/// \param[in] master  The master creates the whole multigrid, the other procs creates an empty multigrid
/// \param[in] meth    Methode that is used to compute the graph-partitioning
/// \param[in] geom    Should geometric information be used to compute the graph-partitioning
/// \param[in] debug   Print information about moved multinodes, edgecut and time
{
    strategy_ = Adaptive;
    xferUnknowns_ = false;
    debugMode_    = debug;

    // Create a multigrid
    mg_ = new MultiGridCL(builder);
    // tell ParMultiGridCL about the multigrid
    ParMultiGridCL::AttachTo(*mg_);
    // Create new LoadBalancingCL
    lb_ = new LoadBalCL(*mg_, part, -1, meth);

    ParTimerCL timer;
    double duration;

    lb_->DeleteGraph();
    if (debugMode_ && ProcCL::IamMaster())
        std::cout << "  - Create dual reduced graph ...\n";

    // Create the graph for partitioning the initial grid
    if (debugMode_) timer.Reset();
    lb_->CreateDualRedGraph(geom);

    if (debugMode_){
        timer.Stop();
        duration = timer.GetMaxTime();
        if (ProcCL::IamMaster()) std::cout << "       --> "<<duration<<" sec\n";
    }

    if (debugMode_ && ProcCL::IamMaster())
        std::cout << "  - Compute Graphpartitioning ...\n";

    if (debugMode_) timer.Reset();
//    if (ProcCL::MyRank()==master)
        lb_->PartitionSer(master);

    if (debugMode_){
        timer.Stop();
        duration = timer.GetMaxTime();
        if (ProcCL::IamMaster()) std::cout << "       --> "<<duration<<" sec\n";
    }

    // distribute the multigrid
    if (debugMode_ && ProcCL::IamMaster())
        std::cout << "  - Migration ...\n";

    if (debugMode_) timer.Reset();
    ParMultiGridCL::XferStart();
    lb_->Migrate();
    ParMultiGridCL::XferEnd();
    if (debugMode_){
        timer.Stop();
        duration = timer.GetMaxTime();
        if (ProcCL::IamMaster()) std::cout << "       --> "<<duration<<" sec\n";
    }
    movedNodes_ = lb_->GetMovedMultiNodes();
    edgeCut_    = lb_->GetEdgeCut();
    lb_->DeleteGraph();
}


void LoadBalHandlerCL::DoMigration()
/** This function encapsulate all necessary steps to perform a loadbalancing
    step. So it create the graph, call ParMetis to compute the partitioning
    and finally do the migration.
*/
{
    // Just do a migration if this is wished
    if (strategy_ == NoMig) return;
    if (ProcCL::Size() == 1){
        std::cout << "Skip migration, because only one proc is involved!\n"; return;
    }

    // Time measurement
    ParTimerCL timer;
    double duration;

    Assert(ProcCL::Size()>1, DROPSErrCL("LoadBalHandlerCL::DoMigration: Only one proc found"), DebugLoadBalC);

    if (debugMode_ && ProcCL::IamMaster()) std::cout << "  - Create dual reduced graph ... \n";

    if (debugMode_) timer.Reset();
    lb_->CreateDualRedGraph();

    if (debugMode_){
        timer.Stop();
        duration = timer.GetMaxTime();
        if (ProcCL::IamMaster()) std::cout << "       --> "<<duration<<" sec\n";
    }

    int allAdjacencies= ProcCL::GlobalSum(lb_->GetNumLocalAdjacencies());

    if (debugMode_ && ProcCL::IamMaster()) std::cout << "  - Compute graph partitioning ... \n";

    if (debugMode_) timer.Reset();

    lb_->PartitionPar();

    if (debugMode_){
        timer.Stop();
        duration = timer.GetMaxTime();
        if (ProcCL::IamMaster()) std::cout << "       --> "<<duration<<" sec for "<<lb_->GetNumAllVerts()<<" vertices and "<<allAdjacencies<<" adjacencies\n";
    }


    if (debugMode_ && ProcCL::IamMaster()) std::cout << "  - Migration ... \n";
    if (debugMode_) timer.Reset();

    ParMultiGridCL::XferStart();
    lb_->Migrate();
    ParMultiGridCL::XferEnd();

    // After the transfer there may be some unknowns on ghost. These can be destroyed now.
//     if (xferUnknowns_){
//         ParMultiGridCL::DeleteUnksOnGhosts();
//     }
    ParMultiGridCL::MarkSimplicesForUnknowns();

    movedNodes_ = lb_->GetMovedMultiNodes();
    edgeCut_    = lb_->GetEdgeCut();

    if (debugMode_){
        timer.Stop();
        duration = timer.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cout << "       --> "<<duration<<" sec\n";
            std::cout << "       --> "<<GetMovedMultiNodes()<<" moved multinodes\n";
            std::cout << "       --> "<<GetEdgeCut()<<" edge cut\n";
        }
    }
    lb_->DeleteGraph();
    lb_->RemoveLbNr();
}


void LoadBalHandlerCL::DoInitDistribution(int)
/** This function encapsulate all necessary steps to distribute the initial
    grid, that is stored on the master processor. So it create the graph, call
    Metis to compute the partitioning and finally do the migration.
*/
{
    if (ProcCL::Size()==1){
        movedNodes_=0;
        edgeCut_=0;
        std::cout << "Skip migration, because only one proc is involved!\n"; return;
    }

    if (debugMode_ && ProcCL::IamMaster())
        std::cout << "  - Create dual reduced graph ...\n";
    lb_->CreateDualRedGraph(true);

    if (debugMode_)
        std::cout << "  - Create graph partition ...\n";
    lb_->PartitionSer( ProcCL::Master());

    if (debugMode_ && ProcCL::IamMaster())
        std::cout << "  - Migration ...\n";
    ParMultiGridCL::XferStart();
    lb_->Migrate();
    ParMultiGridCL::XferEnd();
    movedNodes_ = lb_->GetMovedMultiNodes();
    edgeCut_    = lb_->GetEdgeCut();
    ParMultiGridCL::MarkSimplicesForUnknowns();
    lb_->DeleteGraph();
}
}   // end of namespace DROPS
