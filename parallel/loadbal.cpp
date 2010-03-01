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
#include <iomanip>

namespace DROPS{

/****************************************************************************
* L O A D  B A L  C L A S S                                                 *
****************************************************************************/
// static initialisation
Uint    LoadBalCL::TriangLevel_ = 0;
IFT LoadBalCL::FaceIF_      = 0;
idxtype LoadBalCL::myfirstVert_ = 0;


/****************************************************************************
* W R A P P E R                                                             *
*****************************************************************************
* Converts C++ functions to C functions, so DDD can use them correctly      *
****************************************************************************/
extern "C" int HandlerScatterC( OBJT o, void *d) {return LoadBalCL::HandlerScatter(o,d);}
extern "C" int HandlerGatherC ( OBJT o, void *d) {return LoadBalCL::HandlerGather(o,d);}


/// \brief Constructor
LoadBalCL::LoadBalCL(MultiGridCL& mg, float ub, int TriLevel) : idx_(), ubvec_(ub), useUnkInfo_(false)
/// \param mg       Reference on the multigrid
/// \param ub       imbalace tolerance for eacht vertex weight, suggestion from parMETIS is 1.05
/// \param TriLevel level that should be balanced
/// \todo (of) estimate good parameter ubvec for parmetis
{
    mg_=&mg;
    xadj_=0;
    adjncy_=0;
    vtxdist_=0;
    part_=0;
    vwgt_=0;
    adjwgt_=0;
    xyz_=0;

    if (FaceIF_==0)
        InitIF();

    if (TriLevel<0)
        TriangLevel_=mg_->GetLastLevel();
    else
        TriangLevel_=TriLevel;
    movedMultiNodes_=0;
}


/// \brief Destructor
LoadBalCL::~LoadBalCL()
{
    DeleteGraph();
}


/// \brief Init the DDD interface for communicate adjacenies between processors
void LoadBalCL::InitIF()
/** \todo (of) Falls nicht nur das letzte Triangulierungslevel balanciert werden
    soll, m�ssen auch LB-Nummern auf Ghosts versendet werden!
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
void LoadBalCL::CommunicateAdjacency()
{
	DynamicDataInterfaceCL::IFExchange( FaceIF_, sizeof(idxtype), &HandlerGatherC, &HandlerScatterC);
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
        *sendbuf= tp->GetLbNr() + myfirstVert_;            // put the number into the send buffer

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
    myVerts_=counter;
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


/// \brief Inserts the indices of unknowns into a set and returns this set.
std::set<IdxT> LoadBalCL::UnkOnTetra(const TetraCL& t) const
/** This is used for alternative weights on tetras
    \param t the tetraeder
*/
{
    std::set<IdxT> unknowns;
    if (idx_.empty())
        return unknowns;

    IdxT maxUnk=1;
    for (Uint i=0; i<idx_.size(); ++i)
        maxUnk= std::max( idx_[i]->NumUnknowns(), maxUnk);

    for (Uint i=0; i<idx_.size(); ++i)
    {
        Uint index=idx_[i]->GetIdx();
        if (idx_[i]->NumUnknownsVertex()>0){
            for (TetraCL::const_VertexPIterator it(t.GetVertBegin()), end(t.GetVertEnd()); it!=end; ++it){
                if ((*it)->Unknowns.Exist() && (*it)->Unknowns.Exist(index)){
                    const IdxT dof= (*it)->Unknowns(index);
                    for (Uint j=0; j<idx_[i]->NumUnknownsVertex(); ++j)
                        unknowns.insert(dof+j+i*maxUnk);
                    if ( idx_[i]->IsExtended() )
                        for (Uint j=0; j<idx_[i]->NumUnknownsVertex(); ++j)
                            unknowns.insert(idx_[i]->GetXidx()[dof]+j+i*maxUnk);
                }
            }
        }
        if (idx_[i]->NumUnknownsEdge()>0){
            for (TetraCL::const_EdgePIterator it(t.GetEdgesBegin()), end(t.GetEdgesEnd()); it!=end; ++it){
                if ((*it)->Unknowns.Exist() && (*it)->Unknowns.Exist(index)){
                    const IdxT dof= (*it)->Unknowns(index);
                    for (Uint j=0; j<idx_[i]->NumUnknownsEdge(); ++j)
                        unknowns.insert(dof+j+i*maxUnk);
                    if ( idx_[i]->IsExtended() )
                        for (Uint j=0; j<idx_[i]->NumUnknownsEdge(); ++j)
                            unknowns.insert(idx_[i]->GetXidx()[dof]+j+i*maxUnk);
                }
            }
        }
    }
    return unknowns;
}


/// \brief Compute weight of a tetra
Uint LoadBalCL::GetWeight(const TetraCL& t) const
/** This function computes weights of a tetraeder that is a multinode in the
    graph. Therefore we distinguish between two cases: information about
    unknowns are given or not.
    \param t the tetraeder
*/
{
    // No information about unknowns are given
    if (idx_.empty() || !useUnkInfo_)
    {
        if (t.IsUnrefined())
            return 1;
        else
            return t.GetRefData().ChildNum;
    }
    else
    {
        std::set<IdxT> all_unknowns;
        if (t.IsUnrefined())
            all_unknowns= UnkOnTetra(t);
        else {
            TetraCL::const_ChildPIterator ch    (t.GetChildBegin()),            // iterator through the children
                                          chend ( t.GetChildEnd() );            // last child

            for (; /*ch && */ch!=chend; ++ch){
                if ((*ch)->IsInTriang(TriangLevel_)){                           // if this child is the triangulation, that should be balanced
                    std::set<IdxT> tetra_unks= UnkOnTetra(t);
                    std::insert_iterator< std::set<IdxT> > inserter(all_unknowns, all_unknowns.begin());
                    std::set_union(all_unknowns.begin(), all_unknowns.end(),
                                   tetra_unks.begin(), tetra_unks.end(),
                                   inserter);
                }
            }
        }
        return all_unknowns.size();
    }
}


/// \brief Estimate adjacencies of an unrefined tetra
Uint LoadBalCL::AdjUnrefined( TetraCL& t, int& edgecount)
/** Put the adjacenzies into adjncy_ and estimate the weight of the node, assoziated with the tetra
    if the tetraeder is unrefined
    \param t the tetraeder
    \param edgecount IN/OUT: smallest unused edgenumber
*/
{
    for (Uint face= 0; face<NumFacesC; ++face){                             // Adjacences are faces
        if (!t.IsBndSeg(face)){                                             // If this face belongs to a domain boundary make nothing
            if (t.GetFace(face)->IsOnProcBnd() ){                           // If this face belongs to a boundary between procs
                adjwgt_[edgecount]= 1;                                      // set edgeweight to one
                adjncy_[edgecount++]= t.GetFace(face)->GetLbNeigh();        // put neighbor into adjacenz list
            }
            else{                                                           // Neighbor tetra is on the same proc
                const TetraCL* const neigh= t.GetNeighbor( face);           // get a pointer to this tetra
                Assert(neigh->HasLbNr(),
                       DROPSErrCL("LoadBalCL::AdjUnrefined: Tetra without loadbalance-number"),
                       DebugLoadBalC);
                if (neigh->HasLbNr()){                                      // test if the neighbor is in the LoadBalSet
                    adjwgt_[edgecount]= 1;                                  // set edgeweight to one
                    adjncy_[edgecount++]= neigh->GetLbNr()+myfirstVert_;    // put neighbor into the adjacenz list
                }
            }
        }
    }

    return GetWeight(t);
}


/// \brief Estimate adjacencies of an refined tetra
Uint LoadBalCL::AdjRefined( TetraCL& t, int& edgecount)
/** Put the adjacenzies into adjncy_ and estimate the weight of the node, assoziated with the tetra
    if the tetraeder is unrefined
    \param t the tetraeder
    \param edgecount IN/OUT: smallest unused edgenumber
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
                            ++thisAdj[neigh->GetLbNr() + myfirstVert_]; // put neighbor into the adjazenz map and increase the weight by one
                    }
                }
            }
        }
    }

    thisAdj.erase( t.GetLbNr() + myfirstVert_);                         // delete my one number (set into this adjacenclist by children of t, which have the same number as t)

    GraphEdgeCT::iterator  it = thisAdj.begin(),                        // start at the first adjacenz
                          end = thisAdj.end();                          // end at the last adjacenz
    for (; it!=end; ++it, ++edgecount){                                 // and put all adjacencies into the list for parMETIS
        adjncy_[edgecount]= it->first;
        adjwgt_[edgecount]= it->second;
    }

    return GetWeight(t);
}


/// \brief Compute the number of adjacencies on this proc
Uint LoadBalCL::EstimateAdj()
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
                                adj.insert( neigh->GetLbNr() + myfirstVert_);
                        }
                    }
                }
            }
            // no adjacencies with members of my multinode...
            adj.erase( it->GetLbNr() + myfirstVert_);
            adjCounter+= adj.size();
        }
    }

    myAdjs_ = adjCounter;
    return adjCounter;
}


/// \brief Create the dual reduced Graph for ParMETIS on the last triangulation level
void LoadBalCL::CreateDualRedGraph(bool geom)
/** Set up the graph for (Par)Metis.
    \param geom use of geometric information
    \todo (of) Create Graph of arbitrary triangulation level
*/
{
    Comment("- Setting up dual, reduced graph"<<std::endl,DebugLoadBalC);

    geom_=geom;
    TriangLevel_= mg_->GetLastLevel();

    CreateNumbering();                                          // Create the numbers of the multinodes and set the number of my verts

    // calculate the vtxdist-Array
    vtxdist_ = new int[ProcCL::Size()+1];                       // parameter for parmetis
    IndexArray vtx_rcv= new int[ProcCL::Size()];                // recieve buffer for number of nodes on other procs
    ProcCL::Gather( (int)myVerts_, vtx_rcv, -1);                // communicate the number of nodes

    vtxdist_[0]= 0;                                             // proc 0 starts with node 0
    for (int i=0; i<DynamicDataInterfaceCL::InfoProcs(); ++i)
        vtxdist_[i+1]= vtxdist_[i] + vtx_rcv[i];                // proc i+1 starts with node \sum_{j=0}^i nodesOnProc(i)
    myfirstVert_= vtxdist_[ProcCL::MyRank()];                   // vtxdist_[i] is starting node of proc i

    delete[] vtx_rcv;


     // Compute Adjacencies
    CommunicateAdjacency();                                     // create adjacencies over proc boundaries
    Uint numadj= EstimateAdj();                                 // compute the number of adjacencies on this proc

     // Allocate space for the Arrays
    xadj_   = new idxtype[myVerts_+1];
    adjncy_ = new idxtype[numadj];
    vwgt_   = new idxtype[myVerts_];
    adjwgt_ = new idxtype[numadj];
    if (geom_)
        xyz_ = new float[3*myVerts_];
    part_   = new idxtype[myVerts_];

     // put all nodes and adjacencies into the lists
    LbIteratorCL begin = GetLbTetraBegin(),
                   end = GetLbTetraEnd();
    idxtype vertCount=0,                                        // number of actual node
            edgeCount=0;                                        // number of actual edge

    for ( LbIteratorCL it= begin; it!=end; ++it)                // iterate through the multinodes of the LoadBalSet
    {
        xadj_[vertCount] = edgeCount;                           // remember, where node vertCount writes its neighbors
        if (it->IsUnrefined() )
            vwgt_[vertCount] = AdjUnrefined( *it, edgeCount);   // compute adjacencies, number of adjacencies and node-weight for unrefined tetra
        else
            vwgt_[vertCount] = AdjRefined( *it, edgeCount);     // compute adjacencies, number of adjacencies and node-weight for refined tetra
        if (geom_)                                              // if with geom, compute the barycenter of the tetra
        {
            const Point3DCL coord= GetBaryCenter( *it);
            xyz_[3*vertCount]  = coord[0];
            xyz_[3*vertCount+1]= coord[1];
            xyz_[3*vertCount+2]= coord[2];
        }

        Assert( edgeCount<=(int)numadj,
                DROPSErrCL( "LoadBalanceCL: CreateDualRedGraph: Too much adjacencies!"),
                DebugLoadBalC);

        ++vertCount;                                            // go to the next node
    }

    xadj_[vertCount]= edgeCount;
    Assert( vertCount==myVerts_,
            DROPSErrCL("LoadBalaceCL: CreateDualRedGraph: number of vertices is not correct!"),
            DebugLoadBalC);
}


/// \brief Remove information about the graph
void LoadBalCL::DeleteGraph()
/** Free all allocated arrays */
{
    if (xadj_!=0)    delete[] xadj_;
    if (adjncy_!=0)  delete[] adjncy_;
    if (vtxdist_!=0) delete[] vtxdist_;
    if (part_!=0)    delete[] part_;
    if (vwgt_!=0)    delete[] vwgt_;
    if (adjwgt_!=0)  delete[] adjwgt_;
    if (xyz_!=0)     delete[] xyz_;

    xadj_=0; adjncy_=0; vtxdist_=0; part_=0;
    vwgt_=0; adjwgt_=0; xyz_=0;
    movedMultiNodes_=0;
}


/// \brief Compute a partitioning of the dual reduced graph with ParMetis
void LoadBalCL::ParPartKWay()
/** This procedure uses ParMetis in order to compute a partitioning of the dual
    reduced graph, that has been set up with the member function
    CreateDualRedGraph. No information about a previous distribution of the
    vertices among the processors is used.
    \pre CreateDualRedGraph() must have been called before
*/
{
    Comment("- Start calculate LoadBalanace with ParMETIS"<<std::endl, DebugLoadBalC);

    Assert( xadj_ && adjncy_ && vwgt_ && adjwgt_,
            DROPSErrCL("LoadBalCL::ParPartKWay: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
            DebugLoadBalC);

    int    wgtflag = 3,                 // Weights on vertices and adjacencies are given
           numflag = 0,                 // numbering of verts starts by 0 (C-Style)
           ncon    = 1,                 // one weight per vertex
           nparts  = DynamicDataInterfaceCL::InfoProcs();   // number of subdomains (per proc one)
    float *tpwgts  = new float[nparts], // weight of partion
           ubvec   = ubvec_;            // imbalace tolerance for eacht vertex weight
    int   *options = new int[1];        // default options

    options[0]=0;

    std::fill(tpwgts, tpwgts+nparts, (float)1./(float)nparts);
    if (part_==0)
        part_ = new idxtype[myVerts_];

    MPI_Comm comm = MPI_COMM_WORLD;

    ParMETIS_V3_PartKway(
            vtxdist_, xadj_, adjncy_, vwgt_, adjwgt_,
            &wgtflag, &numflag, &ncon, &nparts,
            tpwgts, &ubvec,0 /*options*/,
            &edgecut_, part_, &comm);

    delete[] tpwgts;
    delete[] options;

}

/// \brief Computes a re-partitioning of a parallel distributed graph
void LoadBalCL::AdaptRepart(float quality)
/** This procedure uses ParMetis in order to compute a partitioning of the dual
    reduced graph, that has been set up with the member function
    CreateDualRedGraph. The previous distribution of the vertices among the
    processors is used.
    \param quality Ratio between communication time to data redistribution time,
      ParMetis recommends 1000
    \pre CreateDualRedGraph() must have been called before
*/
{
    Comment("- Start calculate LoadBalanace with ParMETIS-AdaptiveRepart"<<std::endl,DebugLoadBalC);

    Assert( xadj_ && adjncy_ && vwgt_ && adjwgt_,
            DROPSErrCL("LoadBalCL::AdaptRepart: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
            DebugLoadBalC);

    int    wgtflag = 3,                 // Weights on vertices and adjacencies are given
           numflag = 0,                 // numbering of verts starts by 0 (C-Style)
           ncon    = 1,                 // one weight per vertex
           nparts  = DynamicDataInterfaceCL::InfoProcs();   // number of subdomains (per proc one)
    float *tpwgts  = new float[nparts], // weight of partion
           itr     = quality,           // how much an exchange costs
           ubvec   = ubvec_;            // imbalace tolerance for eacht vertex weight
    int   *options = new int[4];        // default options
//     options[0]=1; options[1]=3; options[2]=15, options[3]=1;    // display times within parmetis
    options[0]=0;                                               // no options for parmetis

    std::fill(tpwgts, tpwgts+nparts, (float)1/(float)nparts);
    if (part_==0)
        part_ = new idxtype[myVerts_];

    MPI_Comm comm = MPI_COMM_WORLD;

    ParMETIS_V3_AdaptiveRepart(
            vtxdist_, xadj_, adjncy_, vwgt_, vwgt_,
            adjwgt_, &wgtflag, &numflag, &ncon, &nparts,tpwgts,
            &ubvec, &itr, options, &edgecut_, part_, &comm);

    delete[] tpwgts;
    delete[] options;
}


/// \brief Compute serial the partitioning of a given Graph in the CSR-format with metis
void LoadBalCL::SerPartKWay(PartMethod meth)
/** Compute the distribution of a graph, that is stored on a single processor in
    order to get a partitioning of the vertices among the processors.
    \param meth specifiy the method that should be used
*/
{
    Comment("- Start calculate LoadBalanace with METIS-"<<(meth==KWay ? "Kway" : "Recursive")<<std::endl, DebugLoadBalC);

    if (myVerts_!=GetNumAllVerts())
        Comment("LoadBalCL: SerPartKWay: This procedure is called by a proc, that does not have all nodes!"<<std::endl, DebugLoadBalC);

    int    wgtflag    = 3,                  // Weights on vertices and adjacencies are given
           numflag    = 0,                  // numbering of verts starts by 0 (C-Style)
           nparts     = DynamicDataInterfaceCL::InfoProcs(),    // number of subdomains (per proc one)
           n          = myVerts_,
           options[5] = {0,0,0,0,0};        // default options

    if (meth==KWay)
        METIS_PartGraphKway(      &n, xadj_, adjncy_, vwgt_, adjwgt_, &wgtflag, &numflag,  &nparts, options,&edgecut_, part_);
    else if (meth==Recursive)
        METIS_PartGraphRecursive( &n, xadj_, adjncy_, vwgt_, adjwgt_, &wgtflag, &numflag,  &nparts, options,&edgecut_, part_);

    Comment("  * Number of Edgecut: "<<edgecut_<<std::endl, DebugLoadBalC);
}

/// \brief Leave partition as it is
void LoadBalCL::IdentityPart()
{
    Assert( xadj_ && adjncy_ && vwgt_ && adjwgt_,
            DROPSErrCL("LoadBalCL::IdentityPart: Graph has not been set up. Maybe use CreateDualRedGraph before calling this routine"),
            DebugLoadBalC);

    if (part_==0)
        part_=part_ = new idxtype[myVerts_];
    std::fill(part_, part_+myVerts_, ProcCL::MyRank());
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

    movedMultiNodes_=0;

    PROCT dest;
    for (LbIteratorCL it= GetLbTetraBegin(), end= GetLbTetraEnd(); it!=end; ++it)
    {
        dest =  static_cast<PROCT>(part_[it->GetLbNr()]);

        if (dest==me) continue;
        movedMultiNodes_++;
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
                    const bool E2Xfer= it.IsInLbSet( **ch) && part_[(*ch)->GetLbNr()]!=static_cast<idxtype>(me);

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

    movedMultiNodes_= ProcCL::GlobalSum(movedMultiNodes_);
}


/// \brief Print the Graph on the ostream
void LoadBalCL::ShowGraph(std::ostream& os)
{
    int counter=0;
    if (xadj_==0)
        os << "No Graph created!"<<std::endl;
    else
    {

        os << "Proc "<<ProcCL::MyRank()<<" stores "<<myVerts_<<" of "<<GetNumAllVerts()
           << " with "<<myAdjs_<<" adjacencies"<<std::endl
           << "The graph is:" <<std::endl
           << " Node | Weight | Dest | Neighbors"<<std::endl
           << "------+--------+------+--------------------------"<<std::endl;

        for (int i=0; i<myVerts_; ++i)
        {
            os << std::setw(5) << i << " |" << std::setw(7) << vwgt_[i] << " |";
            if (part_!=0)
                os << std::setw(5) << part_[i];
            else
                os << "  X  ";
            os <<" | ";
            for (int j=xadj_[i]; j<xadj_[i+1]; ++j)
            {
                os << std::setw(4) <<adjncy_[j] << "  ";
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
    os << lb.myVerts_                                               // number of nodes
       << " " << (lb.myAdjs_/2)                                     // number of "real" adjacencies
       << " 11" <<std::endl;                                        // there are weights onto the nodes and the adjacencies

    for (int i=0; i<lb.myVerts_; ++i)                               // go over all nodes and write the information onto the stream
    {
        os << lb.vwgt_[i] << " ";                                   // first number in the line is the weight of the node i
        for (int j=lb.xadj_[i]; j<lb.xadj_[i+1]; ++j)
            os << lb.adjncy_[j]+1 << " " << lb.adjwgt_[j] << " ";   // then the adjacence with adjacenc weight
        os << std::endl;
    }
    return os;
}


/// \brief Write graph, so that it can be read easily
void LoadBalCL::PrintGraphInfo(std::ostream& os) const
{
    const int me=ProcCL::MyRank(), size = ProcCL::Size();
    const int loc_num_verts=vtxdist_[me+1]-vtxdist_[me],
              loc_num_edges=xadj_[loc_num_verts];

    // print general information
    os << size << ' ' << loc_num_verts << ' ' << loc_num_edges << std::endl;

    // print vtxdist
    for (int i=0; i<=size; ++i){
        os << vtxdist_[i] << ' ';
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print xadj
    for (int i=0; i<=loc_num_verts; ++i){
        os << xadj_[i] << ' ';
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print adjncy_
    for (int i=0; i<loc_num_edges; ++i){
        os << adjncy_[i] << ' ';
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print vwgt
    for (int i=0; i<loc_num_verts; ++i){
        os << vwgt_[i];
        if (i%10==0) os << std::endl;
    }
    os << std::endl;

    // print adjwgt
    for (int i=0; i<loc_num_edges; ++i){
         os << adjwgt_[i];
        if (i%10==0) os << std::endl;
    }
    os << std::endl;
}


/****************************************************************************
* L O A D  B A L  H A N D L E R  C L A S S                                  *
****************************************************************************/

LoadBalHandlerCL::LoadBalHandlerCL(MultiGridCL& mg, float ub) : mg_(&mg)
{
    lb_ = new LoadBalCL(mg, ub);
    strategy_ = Adaptive;
    xferUnknowns_ = false;
    debugMode_    = false;
}

LoadBalHandlerCL::~LoadBalHandlerCL()
{
    if (lb_) delete lb_; lb_=0;
}

LoadBalHandlerCL::LoadBalHandlerCL(const MGBuilderCL &builder, int master, PartMethod meth, bool geom, bool debug)
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
    lb_ = new LoadBalCL(*mg_);

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
    if (ProcCL::MyRank()==master)
        lb_->SerPartKWay(meth);

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
    if (strategy_==NoMig) return;
    if (ProcCL::Size()==1){
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

    // std::cout<< " *** Partioner " << GetStrategy() << std::endl;

    switch (GetStrategy())
    {
        case NoMig     : break;
        case Adaptive  : lb_->AdaptRepart(); break;
        case Recursive : lb_->ParPartKWay(); break;
        case Identity  : lb_->IdentityPart(); break;
        default        : std::cout << "No such method known for migration strategy ERROR\nEXIT"; std::exit(0);
    }

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

    if (debugMode_ && ProcCL::IamMaster())
        std::cout << "  - Create graph partition ...\n";
    if (ProcCL::IamMaster()){
        lb_->SerPartKWay();
    }

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
