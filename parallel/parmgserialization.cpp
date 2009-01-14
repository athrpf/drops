//**************************************************************************
// File:    parmgserializationcl.cpp                                       *
// Content: Declarations of Class ParMGSerializationCL                     *
// Author:  Joerg Peters, Volker Reichelt, Patrick Esser, IGPM RWTH Aachen *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// History: begin - November, 12 2007                                      *
//**************************************************************************

/*******************************************************************
*   P A R  M G  S E R I A L I Z A T I O N   C L                    *
*******************************************************************/

#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "parallel/parmgserialization.h"

namespace DROPS
{

typedef double (*scal_fun)( const DROPS::Point3DCL&, double);
typedef Point3DCL (*vec_fun)( const DROPS::Point3DCL&, double);

inline void ParMGSerializationCL::VertexInfoCL::Init(const VertexCL &src)
{
    Point3DCL tmp = src.GetCoord();
    gid=src.GetGID();
    choords[0]=tmp[0];
    choords[1]=tmp[1];
    choords[2]=tmp[2];
    level = src.GetLevel();
    remove = src.IsMarkedForRemovement();
    isOnBoundary = src.IsOnBoundary();
}

inline void ParMGSerializationCL::BndVertexInfoCL::Init(const VertexCL &srcVertex, const BndPointCL &srcBndVertex)
{
    gid = srcVertex.GetGID();
    bndIdx = srcBndVertex.GetBndIdx();
    choords[0] = srcBndVertex.GetCoord2D()[0];
    choords[1] = srcBndVertex.GetCoord2D()[1];
}

inline void ParMGSerializationCL::EdgeInfoCL::Init(const EdgeCL &src)
{
    gid = src.GetGID();
    level = src.GetLevel();
    vertices[0]=0;
    vertices[1]=0;
    accMFR = src.GetAccMFR();
    if(src.GetVertex(0)!= NULL)
        vertices[0] = (src.GetVertex(0)->GetGID());
    if(src.GetVertex(1)!= NULL)
        vertices[1] = (src.GetVertex(1)->GetGID());

    if(src.GetMidVertex()!=NULL)
        midVertex =src.GetMidVertex()->GetGID();
    else
        midVertex =0;

    remove = src.IsMarkedForRemovement();
// collect boundary information
    bndIdxStart = *(src.GetBndIdxBegin());
    bndIdxEnd = *(src.GetBndIdxBegin() +1);
}

inline void ParMGSerializationCL::FaceInfoCL::Init(const FaceCL &src)
{
    gid=src.GetGID();
    for (int neigh=0; neigh<4; ++neigh){
        neighbor[neigh] = src.GetNeighbor(neigh) ? src.GetNeighbor(neigh)->GetGID() : 0;
    }
    level = src.GetLevel();
    remove = src.IsMarkedForRemovement();
    bndIdx =  src.GetBndIdx();
}

inline void ParMGSerializationCL::TetraInfoCL::Init(const TetraCL &src)
{
    gid=src.GetGID();
    level=src.GetLevel();

    for (Uint vert=0; vert<NumVertsC; ++vert){
        vertices[vert] = src.GetVertex(vert)->GetGID();
    }
    for (Uint edge=0; edge<NumEdgesC; ++edge){
        edges[edge] = src.GetEdge(edge)->GetGID();
    }
    for (Uint face=0; face<NumFacesC; ++face){
        faces[face] = src.GetFace(face)->GetGID();
    }

    refRule    = src.GetRefRule();
    refMark    = src.GetRefMark();
    parentTetra= src.GetParent() ? src.GetParent()->GetGID() : 0;
}

inline void ParMGSerializationCL::TetraChildsInfoCL::Init(const TetraCL &src)
{
    Uint num_children = src.GetRefData().ChildNum;
    gid = src.GetGID();
    for (Uint child=0; child<MaxChildrenC; ++child)
        childs[child] = child<num_children ? src.GetChild(child)->GetGID() : 0;
}

template <>
bool ParMGSerializationCL::AmIResponsibleProc<EdgeCL>(const EdgeCL& s)
/** For edges the priority has to be greater or equal to PrioGhost. An edge with
    PrioVGhost must not store a reference to a midvertex even though the edge
    is refined.
    \param s edge
*/
{
    return s.IsExclusive(PrioGhost);
}


/*!
 * This function is called by all processes but the master-proc.
 * It sends all edges from the current grid back to the master-proc.
 */
void ParMGSerializationCL::SendEdges()
{
    // Move all edges to local buffer
    MoveEdges();
    // Fill Buffer
    ProcCL::RequestT req;
    EdgeInfoCL * pEdges = &edgeBuffer_[0];
    req= ProcCL::Isend(pEdges, edgeBuffer_.size(), edgeStMPI_, masterProc_, edgeTag_);
    ProcCL::Wait(req);
}


/*!
 * Send all vertices back to master-processor
 *
 */
void ParMGSerializationCL::SendVertices()
{
    // Move all vertices to local buffer
    MoveVertices();
    // Send all Vertex and BndVertex data back to master-proc.
    std::valarray<ProcCL::RequestT> req(2);
    VertexInfoCL * pVertices= &vertexBuffer_[0];
    req[0]= ProcCL::Isend(pVertices, vertexBuffer_.size(), vertexStMPI_, masterProc_, vertexTag_);
    BndVertexInfoCL * pBndVertex = &bndVertexBuffer_[0];
    req[1]= ProcCL::Isend(pBndVertex, bndVertexBuffer_.size(), bndVertexStMPI_, masterProc_, bndVertexTag_);
    ProcCL::WaitAll(req);
}

/*!
 * This function is called by all processes but the master-proc.
 * It sends all faces from the current grid back to the master-proc.
 */
void ParMGSerializationCL::SendFaces()
{
    // Move all faces to local buffer
    MoveFaces();
    // Send face-data back to master proc
    ProcCL::RequestT req;
    FaceInfoCL * pFace = &faceBuffer_[0];
    req= ProcCL::Isend(pFace, faceBuffer_.size(), faceStMPI_, masterProc_, faceTag_);
    ProcCL::Wait(req);
}


/*!
 * This function is called by all processes but the master-proc.
 * It sends all tetras from the current grid back to the master-proc.
 */
void ParMGSerializationCL::SendTetras()
{
    // Move all Tetras and their childs to local buffer.
    MoveTetras();
    // Send all tettras and their childs back to master proc.
    std::valarray<ProcCL::RequestT> req(2);
    TetraInfoCL * pTetra = &tetraBuffer_[0];
    req[0]= ProcCL::Isend(pTetra, tetraBuffer_.size(), tetraStMPI_, masterProc_, tetraTag_);
    TetraChildsInfoCL * pTetraChilds = &tetraChildsBuffer_[0];
    req[1]= ProcCL::Isend(pTetraChilds, tetraChildsBuffer_.size(), tetraChildsStMPI_, masterProc_, tetraChildTag_);
    ProcCL::WaitAll(req);
}


/*!
 * Moves all edges from the current grid are moved into the object's buffer.
 */
void ParMGSerializationCL::MoveEdges()
{
    // Count all edges
    Uint numEdges  = 0;
    for (MultiGridCL::const_EdgeIterator sit=mg_.GetAllEdgeBegin(); sit!=mg_.GetAllEdgeEnd(); ++sit)
        if (AmIResponsibleProc(*sit))
            numEdges++;

    // Resize edge buffer
    edgeBuffer_.resize(numEdges);

    Uint cntEdges = 0 ;
    for (MultiGridCL::const_EdgeIterator sit=mg_.GetAllEdgeBegin(); sit!=mg_.GetAllEdgeEnd(); ++sit){
        // only move exclusive edges
        if (!AmIResponsibleProc(*sit))
            continue;

        edgeBuffer_[cntEdges].Init(*sit);
        cntEdges++;
    }
}


/*!
 * Moves all faces from the current grid into the object's buffer.
 */
void ParMGSerializationCL::MoveFaces()
{
    // Count faces
    Uint numFaces = 0 ;
    for (MultiGridCL::const_FaceIterator sit=mg_.GetAllFaceBegin(); sit!=mg_.GetAllFaceEnd(); ++sit)
          numFaces ++ ;

    Uint cntFaces = 0;

    faceBuffer_.resize(numFaces);

    // Fill Buffer with information about each face
    for (MultiGridCL::const_FaceIterator sit=mg_.GetAllFaceBegin(); sit!=mg_.GetAllFaceEnd(); ++sit){
        faceBuffer_[cntFaces].Init(*sit);
        cntFaces ++ ;
    }
}


/*!
 * Moves all vertices from the current grid into the object's buffer.
 */
void ParMGSerializationCL::MoveVertices()
{
    // Count all Vertices and BoundaryVertices
    Uint numVertices = 0 , numBndVertices = 0 ;
    for (MultiGridCL::const_VertexIterator sit=mg_.GetAllVertexBegin(); sit!=mg_.GetAllVertexEnd(); ++sit){
        if (AmIResponsibleProc(*sit)) {
            numVertices ++ ;
            if(sit->IsOnBoundary()){
                for (VertexCL::const_BndVertIt it= sit->GetBndVertBegin(); it != sit->GetBndVertEnd(); ++it)
                    numBndVertices ++ ;
            }
        }
    }

    // Resize all used buffers
    vertexBuffer_.resize(numVertices);
    bndVertexBuffer_.resize(numBndVertices);


    // Fill Buffer
    Uint cntVertex = 0, cntVertexBnd = 0 ;
    for (MultiGridCL::const_VertexIterator sit=mg_.GetAllVertexBegin(); sit!=mg_.GetAllVertexEnd(); ++sit) {
        // only move exclusive vertices
        if(!AmIResponsibleProc(*sit))
            continue;

        // Collect information about vertex
        vertexBuffer_[cntVertex].Init(*sit);

        // Collect all boundary vertices
        if (sit->IsOnBoundary()) {
            for (VertexCL::const_BndVertIt it= sit->GetBndVertBegin(); it != sit->GetBndVertEnd(); ++it) {
                   bndVertexBuffer_[cntVertexBnd].Init(*sit,*it);
                   cntVertexBnd++;
            }
        }
        cntVertex ++ ;
    }
}


/*!
 * Moves all tetras from the current grid into the object's buffer.
 */
void ParMGSerializationCL::MoveTetras()
{
    // Count number of master tetrahedra and parent tetrahedra on this processor
    Uint numMasterTetras=0, numParentTetras=0;
    for (MultiGridCL::const_TetraIterator sit=mg_.GetAllTetraBegin(); sit!=mg_.GetAllTetraEnd(); ++sit)
    {
        if (sit->IsMaster())
            ++numMasterTetras;
        if (sit->GetRefData().ChildNum>0 && !sit->HasGhost())
            ++numParentTetras;
    }

    // Allocate memory for these tetrahedra and parent tetrahedra
    tetraBuffer_.resize(numMasterTetras);
    tetraChildsBuffer_.resize(numParentTetras);

    // Gather all information
    Uint masterTetra=0, parentTetra=0;
    for (MultiGridCL::const_TetraIterator sit=mg_.GetAllTetraBegin(); sit!=mg_.GetAllTetraEnd(); ++sit)
    {
        // Check wheather the tetrahedron is a master tetrahedra. Only these
        // tetras are used for the serialization
        if (sit->IsMaster())
            tetraBuffer_[masterTetra++].Init(*sit);

        // If tetrahedron *sit is refined, collect information about children.
        // Collect this information via tetrahedra that do not have a ghost
        // copy. (Otherwise only ghosts have access to all children)
        if (sit->GetRefData().ChildNum>0 && !sit->HasGhost())
            tetraChildsBuffer_[parentTetra++].Init(*sit);
    }
}


/**!
 * Writes all vertices to disk.
 */
void ParMGSerializationCL::WriteVertices(const std::string &path)
{

    std::ofstream vertex_file((path+"Vertices").c_str());
    std::ofstream bndvtx_file((path+"BoundaryVertices").c_str());

    if (!vertex_file)
        throw ParSerializeErrCL("Cannot open file for writing vertices", 1);
    if (!bndvtx_file)
        throw ParSerializeErrCL("Cannot open file for writing boundary vertices", 1);

    int i=0, j=0;
    // Write all vertices back into output-file
    for ( std::vector<VertexInfoCL>::iterator it=vertexBuffer_.begin() ; it != vertexBuffer_.end(); it++ , i++)
    {
        if (i!=0) vertex_file << '\n';

        vertex_file << AddrMap_getLocalId(addrMapVertex_,it->gid) << " " << std::scientific << std::setprecision(16)
                    << it->choords[0]<<" "<<it->choords[1]<<" "<<it->choords[2]<< " "
                    << it->level         << " " << it->remove;// <<'\n';


        std::vector<BndVertexInfoCL> tmpVector=bndVertexMap_[it->gid];

        if(tmpVector.size()>0)
        {
            for (std::vector<BndVertexInfoCL>::iterator ix=tmpVector.begin() ; ix != tmpVector.end();ix++, ++j) {
                if (j!=0) bndvtx_file << '\n';
                bndvtx_file << AddrMap_getLocalId(addrMapVertex_,ix->gid) << " " << ix->bndIdx << " "
                        << std::scientific << std::setprecision(16) << ix->choords[0] << " " << ix->choords[1];
            }
        }

    }

}


/**!
 * Writes all edges to disk.
 */
void ParMGSerializationCL::WriteEdges(const std::string &path)
{
    std::ofstream edge_file((path+"Edges").c_str());
    if (!edge_file)
        throw ParSerializeErrCL("Cannot open file for writing edges", 1);

    int i=0;
    // Write all edges back into output-file
    for (std::vector<EdgeInfoCL>::iterator it=edgeBuffer_.begin() ; it != edgeBuffer_.end(); it++ , i++) {

        if (i!=0) edge_file << '\n';
        edge_file << AddrMap_getLocalId(addrMapVertex_,it->vertices[0]) << " "
                  << AddrMap_getLocalId(addrMapVertex_,it->vertices[1]) << " "
                  << AddrMap_getLocalId(addrMapVertex_,it->midVertex) << " "
                  << it->bndIdxStart<< " " << it->bndIdxEnd << " " << it->accMFR << " "
                  << it->level << " " << it->remove;
    }
}


/**!
 * Writes all faces to disk.
 */
void ParMGSerializationCL::WriteFaces(const std::string &path)
{

    std::ofstream face_file((path+"Faces").c_str());
    if (!face_file)
        throw ParSerializeErrCL("Cannot open file for writing faces", 1);

    int i=0;
    // Write all faces  back into output-file
    for (std::vector<FaceInfoCL>::iterator it=faceBuffer_.begin() ; it != faceBuffer_.end(); it++ , i++) {

        // Skip non exclusive faces
        if(!it->isExclusive) continue;

        if (i!=0) face_file << '\n';
        face_file << AddrMap_getLocalId(addrMapTetra_,it->neighbor[0]) << " "
                  << AddrMap_getLocalId(addrMapTetra_,it->neighbor[1]) << " "
                  << AddrMap_getLocalId(addrMapTetra_,it->neighbor[2]) << " "
                  << AddrMap_getLocalId(addrMapTetra_,it->neighbor[3]) << " "
                  << it->bndIdx << " " << it->level << " " << it->remove;
    }
}


/**!
 * Writes all tetras the child relations to disk.
 */
void ParMGSerializationCL::WriteTetras(const std::string &path)
{

    std::ofstream tetra_file ((path+"Tetras").c_str());
    std::ofstream child_file ((path+"Children").c_str());
    if (!tetra_file)
        throw ParSerializeErrCL("Cannot open file for writing tetrahedra", 1);
    if (!child_file)
        throw ParSerializeErrCL("Cannot open file for writing children", 1);


    Uint idx=0;

    bool start=true, child_start=true;
    // Write all tetras back into output-file
    for (std::vector<TetraInfoCL>::iterator it=tetraBuffer_.begin() ; it != tetraBuffer_.end(); it++) {

        ++idx;

        if (!start) tetra_file << '\n';

        tetra_file << AddrMap_getLocalId(addrMapTetra_,it->gid) << " "
                   << it->level << " " << it->refRule << " " << it->refMark << " "
                   << AddrMap_getLocalId(addrMapVertex_,it->vertices[0]) << " " << AddrMap_getLocalId(addrMapVertex_,it->vertices[1]) << " "
                   << AddrMap_getLocalId(addrMapVertex_,it->vertices[2]) << " " << AddrMap_getLocalId(addrMapVertex_,it->vertices[3]) << " "
                   << AddrMap_getLocalId(addrMapEdge_,it->edges[0]) << " " << AddrMap_getLocalId(addrMapEdge_,it->edges[1]) << " "
                   << AddrMap_getLocalId(addrMapEdge_,it->edges[2]) << " " << AddrMap_getLocalId(addrMapEdge_,it->edges[3]) << " "
                   << AddrMap_getLocalId(addrMapEdge_,it->edges[4]) << " " << AddrMap_getLocalId(addrMapEdge_,it->edges[5]) << " "
                   << AddrMap_getLocalId(addrMapFace_,it->faces[0]) << " " << AddrMap_getLocalId(addrMapFace_,it->faces[1]) << " "
                   << AddrMap_getLocalId(addrMapFace_,it->faces[2]) << " " << AddrMap_getLocalId(addrMapFace_,it->faces[3]) << " "
                   << AddrMap_getLocalId(addrMapTetra_,it->parentTetra);
        start=false;
    } // end for: write all tetras

    // Write all tetra childs back into output-file
    for (std::vector<TetraChildsInfoCL>::iterator it=tetraChildsBuffer_.begin() ; it != tetraChildsBuffer_.end(); it++) {
        if (!child_start)
            child_file << '\n';
        else
            child_start=false;
        child_file << AddrMap_getLocalId(addrMapTetra_,it->gid) << " ";
        for (Uint i=0; i<8; ++i) {
            child_file << AddrMap_getLocalId(addrMapTetra_,it->childs[i]) << " ";
        }
    }
}


/**
 * Receives the data from all running processes. This function is only called
 * by the master process
 */
void ParMGSerializationCL::FetchData()
{
    Assert(ProcCL::MyRank()==masterProc_,
           ParSerializeErrCL("ParMGSerializationCL::FetchData: Called by non master processor", 1), DebugOutPutC);


    // collect information from all processes
    for(int proc= 0 ; proc < ProcCL::Size(); proc ++ )
    {
        // skip master
        if(proc == masterProc_)
            continue;

        // Status for checking number of received elements
        ProcCL::StatusT status;

        Comment("Receive data from Proc:  "<<proc<<std::endl, DebugOutPutC);

        //---------------------------
        // Receive vertices
        //---------------------------
        ProcCL::Probe(proc, vertexTag_, status);
        Uint numRecieveVertices= ProcCL::GetCount(status, vertexStMPI_);
        // Enlarge local buffer
        vertexBuffer_.reserve(vertexBuffer_.size() + numRecieveVertices);
        VertexInfoCL * pRecvBuffer1 = new VertexInfoCL[numRecieveVertices];
        ProcCL::Recv(pRecvBuffer1, numRecieveVertices, vertexStMPI_, proc, vertexTag_);
        for(Uint recCnt = 0 ; recCnt < numRecieveVertices ; recCnt ++)
            vertexBuffer_.push_back(pRecvBuffer1[recCnt]);
        delete[] pRecvBuffer1;
        Comment("\t"<<numRecieveVertices<<" vertices received!\n", DebugOutPutC);

        //---------------------------
        // Receive boundary vertices
        //---------------------------
        ProcCL::Probe(proc, bndVertexTag_, status);
        Uint numRecieveBndVertices= (Uint)ProcCL::GetCount(status, bndVertexStMPI_);
        // Enlarge local buffer
        bndVertexBuffer_.reserve(bndVertexBuffer_.size() + numRecieveBndVertices);
        BndVertexInfoCL * pRecvBuffer2 = new BndVertexInfoCL[numRecieveBndVertices];
        ProcCL::Recv(pRecvBuffer2, numRecieveBndVertices, bndVertexStMPI_, proc, bndVertexTag_);
        for(Uint recCnt = 0 ; recCnt < numRecieveBndVertices ; recCnt ++)
            bndVertexBuffer_.push_back(pRecvBuffer2[recCnt]);
        delete[] pRecvBuffer2;
        Comment("\t"<<numRecieveBndVertices<<" boundary vertices received!\n", DebugOutPutC);

        //---------------------------
        // Receive edges
        //---------------------------
        ProcCL::Probe(proc, edgeTag_, status);
        Uint numRecieveEdges= ProcCL::GetCount(status, edgeStMPI_);
        // Enlarge local buffer
        edgeBuffer_.reserve(edgeBuffer_.size() + numRecieveEdges);
        EdgeInfoCL * pRecvBuffer3 = new EdgeInfoCL[numRecieveEdges];
        ProcCL::Recv(pRecvBuffer3, numRecieveEdges, edgeStMPI_, proc, edgeTag_);
        for(Uint recCnt = 0 ; recCnt < numRecieveEdges ; recCnt ++)
            edgeBuffer_.push_back(pRecvBuffer3[recCnt]);
        delete[] pRecvBuffer3;
        Comment("\t"<<numRecieveEdges<<" edges received!\n", DebugOutPutC);

        //---------------------------
        // Receive faces
        //---------------------------
        ProcCL::Probe(proc, faceTag_, status);
        Uint numRecieveFaces= (Uint)ProcCL::GetCount(status, faceStMPI_);
        // Enlarge local buffer
        faceBuffer_.reserve(faceBuffer_.size() + numRecieveFaces);
        FaceInfoCL * pRecvBuffer4 = new FaceInfoCL[numRecieveFaces];
        ProcCL::Recv(pRecvBuffer4, numRecieveFaces, faceStMPI_, proc, faceTag_);
        for(Uint recCnt = 0 ; recCnt < numRecieveFaces ; recCnt ++)
            faceBuffer_.push_back(pRecvBuffer4[recCnt]);
        delete[] pRecvBuffer4;
        Comment("\t"<<numRecieveFaces<<" faces received!\n", DebugOutPutC);

        //---------------------------
        // Receive tetras
        //---------------------------
        ProcCL::Probe(proc, tetraTag_, status);
        Uint numRecieveTetras= (Uint) ProcCL::GetCount(status, tetraStMPI_);
        // Enlarge local buffer
        tetraBuffer_.reserve(tetraBuffer_.size() + numRecieveTetras);
        TetraInfoCL * pRecvBuffer5 = new TetraInfoCL[numRecieveTetras];
        ProcCL::Recv(pRecvBuffer5, numRecieveTetras, tetraStMPI_, proc, tetraTag_);
        for(Uint recCnt = 0 ; recCnt < numRecieveTetras ; recCnt ++)
            tetraBuffer_.push_back(pRecvBuffer5[recCnt]);
        delete[] pRecvBuffer5;
        Comment("\t"<<numRecieveTetras<<" tetras received!\n", DebugOutPutC);

        //---------------------------
        // Receive  tetra childs
        //---------------------------
        ProcCL::Probe(proc, tetraChildTag_, status);
        Uint numRecieveTetraChilds= (Uint)ProcCL::GetCount(status, tetraChildsStMPI_);
        // Enlarge local buffer
        tetraChildsBuffer_.reserve(tetraChildsBuffer_.size() + numRecieveTetraChilds);
        TetraChildsInfoCL * pRecvBuffer6 = new TetraChildsInfoCL[numRecieveTetraChilds];
        ProcCL::Recv(pRecvBuffer6, numRecieveTetraChilds, tetraChildsStMPI_, proc, tetraChildTag_);
        for(Uint recCnt = 0 ; recCnt < numRecieveTetraChilds ; recCnt ++)
            tetraChildsBuffer_.push_back(pRecvBuffer6[recCnt]);
        delete[] pRecvBuffer6;
        Comment("\t"<<numRecieveTetraChilds<<" children information received!\n", DebugOutPutC);
    } // end for : fetch data from all procs


} // END OF FUNCTION


/// \brief Registers the mpi-datatype used to transmit vertex-data
void ParMGSerializationCL::CreateVertexMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (vertexStMPI_ != ProcCL::NullDataType)
        return;

    VertexInfoCL tmp;
    const int count               = 4;
    const int block[4]            = {2,3,2,2};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[4]= { ProcCL::MPI_TT<Uint>::dtype, ProcCL::MPI_TT<double>::dtype,
                                      ProcCL::MPI_TT<int>::dtype, ProcCL::MPI_TT<Usint>::dtype };
    ProcCL::AintT displace[4];

    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.choords[0]) -sAddr;
    displace[2]= ProcCL::Get_address(&tmp.remove) -sAddr;
    displace[3]= ProcCL::Get_address(&tmp.bndIdxStart) -sAddr;

    vertexStMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(vertexStMPI_);
}


/// \brief Registers the mpi-datatype used to transmit boundary vertex-data
void ParMGSerializationCL::CreateBndVertexMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (bndVertexStMPI_ != ProcCL::NullDataType)
        return;

    BndVertexInfoCL tmp;
    const int count               = 2;
    const int block[2]            = {2,2};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[2]= { ProcCL::MPI_TT<Uint>::dtype, ProcCL::MPI_TT<double>::dtype };
    ProcCL::AintT displace[2];

    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.choords[0]) -sAddr;

    bndVertexStMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(bndVertexStMPI_);
}


/// \brief Registers the mpi-datatype used to transmit edge-data
void ParMGSerializationCL::CreateEdgeMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (edgeStMPI_ != ProcCL::NullDataType)
        return;

    EdgeInfoCL tmp;
    const int count               = 3;
    const int block[3]            = {5,3,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[3]= { ProcCL::MPI_TT<Uint>::dtype, ProcCL::MPI_TT<Usint>::dtype,
                                      ProcCL::MPI_TT<int>::dtype };
    ProcCL::AintT displace[3];

    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.accMFR) -sAddr;
    displace[2]= ProcCL::Get_address(&tmp.remove) -sAddr;

    edgeStMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(edgeStMPI_);
}


/// \brief Registers the mpi-datatype used to transmit face-data
void ParMGSerializationCL::CreateFaceMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (faceStMPI_ != ProcCL::NullDataType)
        return;

    FaceInfoCL tmp;
    const int count               = 3;
    const int block[3]            = {6,1,2};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[3]= { ProcCL::MPI_TT<Uint>::dtype, ProcCL::MPI_TT<Usint>::dtype,
                                      ProcCL::MPI_TT<int>::dtype};
    ProcCL::AintT displace[3];

    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.bndIdx) -sAddr;
    displace[2]= ProcCL::Get_address(&tmp.remove) -sAddr;

    faceStMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(faceStMPI_);
}


/// \brief Registers the mpi-datatype used to transmit tetra-data
void ParMGSerializationCL::CreateTetraMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (tetraStMPI_ != ProcCL::NullDataType)
        return;

    TetraInfoCL tmp;
    const int count               = 1;
    const int block[1]            = {19};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[1]= { ProcCL::MPI_TT<Uint>::dtype};
    ProcCL::AintT displace[1];
    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;

    tetraStMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(tetraStMPI_);
}


/// \brief Registers the mpi-datatype used to transmit tetra-child-data
void ParMGSerializationCL::CreateTetraChildsMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (tetraChildsStMPI_ != ProcCL::NullDataType)
        return;

    TetraChildsInfoCL tmp;
    const int count               = 1;
    const int block[1]            = {9};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[1]= {ProcCL::MPI_TT<Uint>::dtype};
    ProcCL::AintT displace[1];
    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;

    tetraChildsStMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(tetraChildsStMPI_);
}

/// \brief Registers the mpi-datatype used to transmit a scalar dof value
void ParMGSerializationCL::CreateDofScalarMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (dofScalarMPI_ != ProcCL::NullDataType)
        return;

    DofScalarCL tmp;
    const int count               = 2;
    const int block[2]            = {1,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[2]= {ProcCL::MPI_TT<Uint>::dtype, ProcCL::MPI_TT<double>::dtype };

    ProcCL::AintT displace[2];
    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.data)-sAddr;

    dofScalarMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(dofScalarMPI_);
}

/// \brief Registers the mpi-datatype used to transmit a vectorial dof value
void ParMGSerializationCL::CreateDofVecMPI()
{
    // if MPI-Datatype has been already created, do nothing
    if (dofVecMPI_ != ProcCL::NullDataType)
        return;

    DofVecCL tmp;
    const int count               = 2;
    const int block[2]            = {1,3};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp.gid);
    const ProcCL::DatatypeT typ[2]= {ProcCL::MPI_TT<Uint>::dtype, ProcCL::MPI_TT<double>::dtype};
    ProcCL::AintT displace[2];

    displace[0]= ProcCL::Get_address(&tmp.gid)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.data)-sAddr;

    dofVecMPI_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(dofVecMPI_);
}

/// \brief Assign each Face a natural number and mark copies of faces
void ParMGSerializationCL::RegisterFace(FaceInfoCL& f)
/** Register faces only once and store neighbors at position 0,1 for neighbors
    lying on the same level as the face and on position 2,3 for green children
    of parental neighbors.
*/
{
    if ( AddrMap_storeGID(addrMapFace_,f.gid) ){
        // Face is not known, so mark this face as exclusive
        f.isExclusive= true;
    }
    else{
        // remeber, that this face has not to be written into the file
        f.isExclusive= false;

        // search for stored face
        std::vector<FaceInfoCL>::iterator storedFace= faceBuffer_.begin();
        for ( ; storedFace != faceBuffer_.end(); storedFace++){
            if (storedFace->gid==f.gid && storedFace->isExclusive)
                break;
        }
        Assert(storedFace!=faceBuffer_.end() && storedFace->gid==f.gid,
               ParSerializeErrCL("ParMGSerializationCL::RegisterFace: No or wrong face found", 1), DebugOutPutC);

        for (Uint neigh=0; neigh<4; ++neigh){
            // get number of neighbor tetrahedron
            const Uint neighNo=f.neighbor[neigh];

            // If neighbor is not known, go on
            if (neighNo==0)
                continue;

            // Check whether neighbor is on same level as the face or on next finer level
            const Uint offset= neigh<2 ? 0 : 2;

            // If neighbor is already known, go on
            if (storedFace->neighbor[offset+0]==neighNo || storedFace->neighbor[offset+1]==neighNo)
                continue;

            // right here this neighbor is not known and should be set, so check this
            Assert(storedFace->neighbor[offset+0]==0 || storedFace->neighbor[offset+1]==0,
                   ParSerializeErrCL("ParMGSerializationCL::CreateLocalAddrMap: To many neighbors for a face", 1), DebugOutPutC);

            if (storedFace->neighbor[offset+0]==0)
                storedFace->neighbor[offset+0]= neighNo;
            else
                storedFace->neighbor[offset+1]= neighNo;
        }
    }
}


void ParMGSerializationCL::MakeNeighborsConsistent(FaceInfoCL& f)
// As declared in the assumptions of FaceCL, a green child and its parent with a
// common face are both stored on the same side. Therefore the neighbors 2 and 3
// (if at least one of them exists) have to be sorted. This has to be done only
// for exclusive faces
{
    if (!f.neighbor[2])
        return;

    // Check if parent of neighbor 2 is stored at position 0
    const Uint neighIdx= AddrMap_getLocalId(addrMapTetra_, f.neighbor[2])-1;
    const Uint parentNo= tetraBuffer_[neighIdx].parentTetra;
    Assert(parentNo!=0,
           ParSerializeErrCL("ParMGSerializationCL::MakeNeighborsConsistent: No parent found", 1), DebugOutPutC);
    if (parentNo!=f.neighbor[0])
        std::swap(f.neighbor[0], f.neighbor[1]);

    // Check for errors
    if (f.neighbor[2]){
        Assert(f.neighbor[0]==tetraBuffer_[AddrMap_getLocalId(addrMapTetra_, f.neighbor[2])-1].parentTetra,
               ParSerializeErrCL("ParMGSerializationCL::MakeNeighborsConsistent: Wrong neighbors in face at position 2", 1), DebugOutPutC);
    }
    if (f.neighbor[3]){
        Assert(f.neighbor[1]==tetraBuffer_[AddrMap_getLocalId(addrMapTetra_, f.neighbor[3])-1].parentTetra,
               ParSerializeErrCL("ParMGSerializationCL::MakeNeighborsConsistent: Wrong neighbors in face at position 3", 1), DebugOutPutC);
    }

}


/// \brief Prepare data for writing
void ParMGSerializationCL::PrepareData()
/** Sorts all received data and initializes the pAddrMap that represent a memory
    map from the natural-numbers into theGIDs of all elements.*/
{
    // Sort all buffers
    std::sort(vertexBuffer_.begin(),vertexBuffer_.end(),BasicGeomCL::compareLevel);
    std::sort(faceBuffer_.begin(),faceBuffer_.end(),BasicGeomCL::compareLevel);
    std::sort(edgeBuffer_.begin(),edgeBuffer_.end(),BasicGeomCL::compareLevel);
    std::sort(tetraBuffer_.begin(),tetraBuffer_.end(),BasicGeomCL::compareLevel);

    // Register GIDS for vertices, edges and tetras
    for ( std::vector<VertexInfoCL>::iterator it=vertexBuffer_.begin() ; it != vertexBuffer_.end(); it++ )
        AddrMap_storeGID(addrMapVertex_,it->gid);

    for ( std::vector<EdgeInfoCL>::iterator it=edgeBuffer_.begin() ; it != edgeBuffer_.end(); it++ )
        AddrMap_storeGID(addrMapEdge_,it->gid);

    for ( std::vector<TetraInfoCL>::iterator it=tetraBuffer_.begin() ; it != tetraBuffer_.end(); it++ )
        AddrMap_storeGID(addrMapTetra_,it->gid);

    // Register face GIDs in two steps
    for ( std::vector<FaceInfoCL>::iterator it=faceBuffer_.begin() ; it != faceBuffer_.end(); it++ )
        RegisterFace(*it);                       // register faces
    for ( std::vector<FaceInfoCL>::iterator it=faceBuffer_.begin() ; it != faceBuffer_.end(); it++ )
        if (it->isExclusive)
            MakeNeighborsConsistent(*it);       // register order neighbors

    // move bndVertices from the receive-vector in an associative-map

    for (std::vector<BndVertexInfoCL>::iterator it=bndVertexBuffer_.begin() ; it != bndVertexBuffer_.end();it++)
        bndVertexMap_[it->gid].push_back(*it);


}


/**
 *  Write out the grid data into files in the given folder with the desired prefix.
 *  This function internally allocates some memory to store the edges,vertices etc. This memory is normally cleaned up
 *  with the next call of WriteMG().
 *  If you are not going to use WriteMG() anymore and want free the memory, just call ::Clear() afterwards.
 */
void ParMGSerializationCL::WriteMG(const std::string path)
{
    // Clean up memory allocated by former use of WriteMG
    Clear();

    std::string outputPath= path_;
    if(path.length()>0)
        outputPath = path ;

    if(masterProc_ == ProcCL::MyRank())
    {
        // check path
        std::string testfile_name = outputPath+"test.txt";
        std::ofstream testfile(testfile_name.c_str());
        if (!testfile)
            throw ParSerializeErrCL("Cannot create files in given directory!", 0);
        testfile.close();
        DeleteFile(testfile_name);

        // Move all geom elements
        // from MultiGrid to local object buffer.
        MoveVertices();
        MoveEdges();
        MoveFaces();
        MoveTetras();

        Comment("Start transfer\n", DebugOutPutC);
        FetchData();
        Comment("Transfer finished\n", DebugOutPutC);

        // Create adress map
        PrepareData();
#if DROPSDebugC&DebugOutPutC
        printStat();
#endif

        // Write all buffered geom-data to disk
        WriteVertices (outputPath);
        WriteEdges(outputPath);
        WriteFaces(outputPath);
        WriteTetras(outputPath);
    }
    else
    {
        // Send all simplices to master processor
        SendVertices();
        SendEdges();
        SendFaces();
        SendTetras();
    }
}

/**  This procedure collects information about dofs on vertices and edges (if
    there are dofs). The dofs are put into the vector under the same conditions
    as within the WriteMG() procedure, i.e. if the calling processor has
    collected infomation about a simplex the same processor collect the dof
    on that simplex.
    \param vec          describer of the dof
    \param localDofVert vector of all dof of the local processor on vertices
    \param localDofEdge vector of all dof of the local processor on edges
*/
void ParMGSerializationCL::CollectDOF(const VecDescCL* vec)
{
    IdxT idx=vec->RowIdx->GetIdx();

    // Collect dofs from all vertices
    for (MultiGridCL::const_VertexIterator sit=mg_.GetAllVertexBegin(); sit!=mg_.GetAllVertexEnd(); ++sit)
    {
        // Check if the proc is responsible for the current vertex
        if (sit->IsExclusive(PrioHasUnk)){
            if(sit->Unknowns.Exist(idx)){
                switch(vec->RowIdx->NumUnknownsVertex())
                {
                    // Scalar value
                    case 1:
                        DofScalarCL tmpS;
                        tmpS.gid = sit->GetGID();
                        tmpS.data= vec->Data[sit->Unknowns(idx)];
                        dofScalBuffer_.push_back(tmpS);
                        break;

                    // Vectorial value
                    case 3:
                        DofVecCL tmpV;
                        tmpV.gid = sit->GetGID();
                        for (Uint i=0; i<vec->RowIdx->NumUnknownsVertex(); ++i){
                            tmpV.data[i]=(vec->Data[sit->Unknowns(idx)+i]);
                        }
                        dofVecBuffer_.push_back(tmpV);
                        break;

                    // Unexpected type of data
                    default:
                        Assert(false, ParSerializeErrCL(" Invalid number of unknowns!", 1),DebugOutPutC);
                        break;
                }   // END SWITCH
            }   // END IF (Unknowns Exist)
        }   // END IF (IsExclusive)
    }   // END FOR

    // If there are dofs on the edges, collect them too.
    if (vec->RowIdx->NumUnknownsEdge()>0)
    {
        // Resize Buffer
        for (MultiGridCL::const_EdgeIterator sit=mg_.GetAllEdgeBegin(); sit!=mg_.GetAllEdgeEnd(); ++sit)
        {
            // Check if the proc is responsible for the current edge
            if (sit->IsExclusive(PrioHasUnk)){
                if(sit->Unknowns.Exist(idx)){
                    switch(vec->RowIdx->NumUnknownsEdge())
                    {
                        // Scalar value
                        case 1:
                            DofScalarCL tmpS;
                            tmpS.gid = sit->GetGID();
                            tmpS.data= vec->Data[sit->Unknowns(idx)];
                            dofScalBuffer_.push_back(tmpS);
                            break;

                        // Vectorial value
                        case 3:
                            DofVecCL tmpV;
                            tmpV.gid = sit->GetGID();

                            // Collect all data for the vector
                            for (Uint i=0; i<vec->RowIdx->NumUnknownsEdge(); ++i){
                                tmpV.data[i]=(vec->Data[sit->Unknowns(idx)+i]);
                            }
                            dofVecBuffer_.push_back(tmpV);
                            break;

                        // Unexpected type of data
                        default:
                            Assert(false, ParSerializeErrCL("Invalid number of unknowns.", 1), DebugOutPutC);
                            break;
                    }   // END OF SWITCH
                }   // END IF (Unknowns Exist)
            }   // END IF (IsExclusive)
        }   // END FOR
    }   // END IF (Unknowns on edges)
}

/**
 * Free the local buffers for edge and vertex dofs.
 */
void ParMGSerializationCL::FreeLocalDOFBuffer()
{
    // Buffer wieder freigeben
    dofScalBuffer_.clear();
    dofScalBuffer_.resize(0);
    dofVecBuffer_.clear();
    dofVecBuffer_.resize(0);

    dofScalMap_.clear();

    for (std::map<Uint,double*>::iterator it=dofVecMap_.begin() ; it != dofVecMap_.end(); it++ )
        delete (*it).second;

    dofVecMap_.clear();
}

/**  Send all dof on vertices and all dof on edges in two messages to the master
    processor. This processor writes in a proceeding step all dof sorted into
    a file.
    \param localDofVert local dofs on vertices
    \param localDofEdge local dofs on edges
 */
void ParMGSerializationCL::SendDOF()
{
    std::valarray<ProcCL::RequestT> req(2);
    // Send scalar data
    if(dofScalBuffer_.size()>0)
        req[0]= ProcCL::Isend(dofScalBuffer_, dofScalarMPI_, masterProc_, dofScalarTag_);

    // Send vectorial data
    if(dofVecBuffer_.size()>0)
        req[1]= ProcCL::Isend(dofVecBuffer_, dofVecMPI_, masterProc_, dofVecTag_);

    // Wait for finishing all send operations
    if(dofScalBuffer_.size()>0)
        ProcCL::Wait(req[0]);
    if(dofVecBuffer_.size()>0)
        ProcCL::Wait(req[1]);

}
/**
 * Receive all vertex and edge dofs from the designated processor
 * \param recieveBufferVerts All received vertex dofs will be stored in this buffer.
 * \param recieveBufferEdges All received edge dofs will be stored in this buffer.
 * \param p Number of the processor.
 */
void ParMGSerializationCL::RecieveDOF(const VecDescCL* vec, int p)
{
    // Status for checking number of recieved elements
    ProcCL::StatusT status;

    // Function can not be called from the master
    if (p==masterProc_ )
      return ;

    DofScalarCL *pRecvBufferScal = NULL ;
    DofVecCL    *pRecvBufferVec  = NULL ;

    Uint numRecvDof= 0;

    // Receive
    switch(vec->RowIdx->NumUnknownsVertex())
    {
        // recieve scalar data
        case 1:
            ProcCL::Probe(p, dofScalarTag_, status);
            numRecvDof= (Uint)ProcCL::GetCount(status, dofScalarMPI_);
            pRecvBufferScal = new DofScalarCL[numRecvDof*sizeof(DofScalarCL)];
            ProcCL::Recv(pRecvBufferScal, numRecvDof, dofScalarMPI_, p, dofScalarTag_);

            // map local gids to their new gid on the master
            ExpandDOFMap(pRecvBufferScal,numRecvDof);

            break;

        // recieve vectorial data
        case 3:
            ProcCL::Probe(p, dofVecTag_, status);
            numRecvDof= (Uint)ProcCL::GetCount(status, dofVecMPI_);
            pRecvBufferVec = new DofVecCL[numRecvDof*sizeof(DofVecCL)];
            ProcCL::Recv(pRecvBufferVec, numRecvDof, dofVecMPI_, p, dofVecTag_);

            // map local gids to their new gid on the master
            ExpandDOFMap(pRecvBufferVec,numRecvDof);

            break;

        // recieve unknown data
        default:
            break;

    }; // END OF SWITCH

    delete pRecvBufferScal;
	delete pRecvBufferVec;
}
/**
 * \brief Expand the current DOF-Map with dof supplied in buffer.
 * \param pBuffer A buffer containing scalar dof values
 * \param bufferSize Number of elements in pBuffer
 **/
void ParMGSerializationCL::ExpandDOFMap(const DofScalarCL* pBuffer,Uint bufferSize )
{
    if(pBuffer==NULL) return ;

    for(Uint i = 0 ; i < bufferSize ; i++)
        dofScalMap_[pBuffer[i].gid] = pBuffer[i].data;
}


/**
 * \brief Expand the current DOF-Map with dof supplied in buffer.
 * \param pBuffer A buffer containing vectorial dof values
 * \param bufferSize Number of elements in pBuffer
 **/
void ParMGSerializationCL::ExpandDOFMap(const DofVecCL* pBuffer,Uint bufferSize )
{
    if(pBuffer==NULL) return ;

    for(Uint i = 0 ; i < bufferSize ; i++)
    {
        // determine global gid
        //Uint gid = std::max(AddrMap_getLocalId(addrMapVertex_,pBuffer[i].gid),AddrMap_getLocalId(addrMapEdge_,pBuffer[i].gid));

        double * tmp = new double[3];
        tmp[0] = pBuffer[i].data[0];
        tmp[1] = pBuffer[i].data[1];
        tmp[2] = pBuffer[i].data[2];

        dofVecMap_[pBuffer[i].gid] = tmp;
    }
}


void ParMGSerializationCL::WriteDOF(const VecDescCL* vec, const std::string& name,const std::string path)
/** This procedure writes out the values of the DOF given by the vec into a file
    specified by the name that is also given as a parameter. I.e. all information
    are gathered by all processes and written out by the master process. For the
    sake of simplicity it is assumed that only unkwons on vertices and edges are
    present. This is consistent to the rest of the DROPS package.
    \param vec  describer of the dof
    \param name name of the unknowns and the file where to put the unknowns.
*/
{
    // Create used MPI datatypes
    CreateDofScalarMPI();
    CreateDofVecMPI();

    //
    //  Master-Proc
    //
    if (ProcCL::MyRank()==masterProc_){

        // Check path
        std::string testfile_name = path_+"test.txt";
        std::ofstream testfile(testfile_name.c_str());
        if (!testfile)
            throw ParSerializeErrCL("Cannot create files in given directory!", 0);
        testfile.close();
        DeleteFile(testfile_name);


        // Collect DOF from Master-Proc
        CollectDOF(vec);

        // Move local DOFs to GID map

        if(dofScalBuffer_.size()>0)
            ExpandDOFMap(&dofScalBuffer_[0],dofScalBuffer_.size());
        else if (dofVecBuffer_.size()>0)
            ExpandDOFMap(&dofVecBuffer_[0],dofVecBuffer_.size());


        // Receive Edges/Vertice DOFs

        for (int p=0; p<ProcCL::Size(); ++p){
            if( p == ProcCL::MyRank() )
                continue;
            RecieveDOF(vec, p);
        }

       // Open files

        const std::string commentStr("# Written Dofs: ");
        std::string placeHolder ="#";

        placeHolder.resize(72 ,' ');

        std::string outputPath = path_;
        if(path.length()>0)
                outputPath = path;

        std::ofstream vertexDofFile((outputPath+name+"_Vertices").c_str());
        // Reserve some bytes at the beginning of file
        vertexDofFile<<placeHolder<< std::endl;
        std::ofstream edgeDofFile((outputPath+name+"_Edges").c_str());
        edgeDofFile<<placeHolder<< std::endl;

        Uint writtenVertexDofCnt = 0,writtenEdgeDofCnt = 0 ;

        // Write out scalar or vectorial data
        switch(vec->RowIdx->NumUnknownsVertex())
        {
            // scalar data
            case 1:
                // vertices
                for ( std::vector<VertexInfoCL>::iterator it=vertexBuffer_.begin() ; it != vertexBuffer_.end(); it++)
                {
                    vertexDofFile << std::scientific << std::setprecision(16) << dofScalMap_[it->gid] << std::endl << std::flush;
                    writtenVertexDofCnt++;
                }

                // edges
                for (std::vector<EdgeInfoCL>::iterator it=edgeBuffer_.begin() ; it != edgeBuffer_.end(); it++)
                {

                    Uint gid = it->gid; //  AddrMap_getLocalId(addrMapEdge_,it->gid) ;

                    // skip edge is no unknowns on it
                    if(dofScalMap_[gid]==0)
                        continue;

                    edgeDofFile << std::scientific << std::setprecision(16) << dofScalMap_[gid]<< std::endl << std::flush ;
                    writtenEdgeDofCnt++;
                }


                break;

            // vectorial data
            case 3:

                // vertices
                for ( std::vector<VertexInfoCL>::iterator it=vertexBuffer_.begin() ; it != vertexBuffer_.end(); it++ )
                {

                    Uint gid = it->gid ; //  AddrMap_getLocalId(addrMapVertex_,it->gid) ;

                    if(dofVecMap_[gid] == 0 )
                        continue;

                    vertexDofFile << std::scientific << std::setprecision(16) << dofVecMap_[gid][0] << " "
                             <<dofVecMap_[gid][1]<<" "
                             <<dofVecMap_[gid][2] <<std::endl;

                    writtenVertexDofCnt++;
                }


                // edges
                for (std::vector<EdgeInfoCL>::iterator it=edgeBuffer_.begin() ; it != edgeBuffer_.end(); it++)
                {

                    Uint gid = it->gid;

                    // skip edge is no unknowns on it
                    if(dofVecMap_[gid]==0)
                        continue;


                    edgeDofFile << std::scientific << std::setprecision(16) << dofVecMap_[gid][0] << " "
                             << dofVecMap_[gid][1]<<" "
                             << dofVecMap_[gid][2] <<std::endl;

                    writtenEdgeDofCnt++;
                }



                break;

        }; // END OF SWITCH

        // Move put-pointer back to the file beginning
        vertexDofFile.seekp(0);
        vertexDofFile <<commentStr << writtenVertexDofCnt << std::flush;

        edgeDofFile.seekp(0);
        edgeDofFile << commentStr<< writtenEdgeDofCnt << std::flush;

    }
    //
    //  Worker-Procs
    //
    else
    {
        // collect all dof data
        CollectDOF(vec);
        // send data back to master
        SendDOF();
    }

    // free local dof buffers
    FreeLocalDOFBuffer();
}

/*!
    * Returns the new ID of an registered GID.
    * If the searchGID wasn't registered before 0 is returned and an error message
    * is printed.
*/
Uint ParMGSerializationCL::AddrMap_getLocalId(MapIdToGIDT &tmpMap,Uint searchGID)
{
    if(searchGID == 0 )
        return 0 ;

    MapIdToGIDT::iterator it= tmpMap.find(searchGID);
    if (it != tmpMap.end() )
        return it->second;
    else
        return 0;
}


/*!
    * Stores a new GID into the adress-map and returns it's
    * new id.
    * If the GID was already registered before true is returned.
*/
bool ParMGSerializationCL::AddrMap_storeGID(MapIdToGIDT &tmpMap,Uint newGID)
{
    // check if gid is already registered
    if(tmpMap.find(newGID) == tmpMap.end()){
        Uint newKey = tmpMap.size() +1 ;
        tmpMap[newGID]=newKey;
        return true;
    }
    else{
        return false;
    }
}

ParMGSerializationCL::~ParMGSerializationCL()
/// Free all MPI-datatypes
{
    if ( vertexStMPI_!=ProcCL::NullDataType )
        ProcCL::Free(vertexStMPI_);
    if ( bndVertexStMPI_!=ProcCL::NullDataType )
        ProcCL::Free(bndVertexStMPI_);
    if ( edgeStMPI_!=ProcCL::NullDataType )
        ProcCL::Free(edgeStMPI_);
    if ( faceStMPI_!=ProcCL::NullDataType )
        ProcCL::Free(faceStMPI_);
    if ( tetraStMPI_!=ProcCL::NullDataType )
        ProcCL::Free(tetraStMPI_);
    if ( tetraChildsStMPI_!=ProcCL::NullDataType )
        ProcCL::Free(tetraChildsStMPI_);
    if ( dofScalarMPI_!=ProcCL::NullDataType )
        ProcCL::Free(dofScalarMPI_);
    if ( dofVecMPI_!=ProcCL::NullDataType )
        ProcCL::Free(dofVecMPI_);


}

void ParMGSerializationCL::Clear()
{
    vertexBuffer_.clear();
    bndVertexBuffer_.clear();
    edgeBuffer_.clear();
    faceBuffer_.clear();
    tetraBuffer_.clear();
    tetraChildsBuffer_.clear();
    addrMapVertex_.clear();
    addrMapEdge_.clear();
    addrMapFace_.clear();
    addrMapTetra_.clear();

    // Clear Bnd-Vertex Map and it associated vectors
    bndVertexMap_.clear();
    /// \todo Is it necessary to iterater over all stored vectors and clear them,too?
}


void ParMGSerializationCL::printStat()
{
    std::cout << "+++++++ Queues +++++++\n";
    std::cout << "Vertices:" << vertexBuffer_.size() <<"\n";
    std::cout << "Vertices Bnd:" << bndVertexBuffer_.size() <<"\n";
    std::cout << "Faces:" << faceBuffer_.size() <<"\n";
    std::cout << "Edges:" << edgeBuffer_.size() <<"\n";
    std::cout << "Tetras:" << tetraBuffer_.size() <<"\n";
    std::cout << "Tetra Childs:" << tetraChildsBuffer_.size() <<"\n";
    std::cout << "----------------------\n\tTotal:";
    std::cout << (vertexBuffer_.size() + faceBuffer_.size() + edgeBuffer_.size() +tetraBuffer_.size())<< "\n";
    std::cout << "**** Registerd GIDs:: \n";
    std::cout << "\t" << (addrMapTetra_.size() +addrMapVertex_.size() + addrMapEdge_.size() + addrMapFace_.size()) << "\n";
}

void ReadDOF(const MultiGridCL& mg, VecDescCL* vec, const std::string& file)
/** In order to read dof out of a file, make sure, that the row index within
    vec has been created and numbered correctly before calling this function.
    \param mg   corresponding multigrid
    \param vec  vector describer of the dof
    \param file name of the file, where the dof has been stored
    \pre   row index of vec has to be created
*/
{
    // Open corresponding files for reading
    std::ifstream srcFileVert((file + "_Vertices").c_str());
    std::ifstream srcFileEdge((file + "_Edges").c_str());

    if (!srcFileVert)
        throw ParSerializeErrCL("Error while reading file '"  + file + "'", 1);
    if (!srcFileEdge)
        throw ParSerializeErrCL("Error while reading file '"  + file + "'", 1);

    // Remove comment line at begining of the file
    char buffer[1024];
    srcFileVert.getline(buffer, sizeof(buffer));
    srcFileEdge.getline(buffer, sizeof(buffer));

    const IdxT idx=vec->RowIdx->GetIdx();

    IdxT counter=0;
    // Iterate over all vertices
    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
         sit!=mg.GetTriangVertexEnd() && !srcFileVert.eof();
         ++sit)
    {


        if (sit->Unknowns.Exist(idx)){
            if (srcFileVert.eof()){
                throw ParSerializeErrCL("EndOfFile reached!", 1);
            }
            for (Uint i=0; i<vec->RowIdx->NumUnknownsVertex(); ++i)
            {
                counter++;
                srcFileVert >> vec->Data[ sit->Unknowns(idx)+i ];
            }
        }
    }

    counter=0;

    // If data exists, also iterate over the edges
    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
         (sit!=mg.GetTriangEdgeEnd() && !srcFileEdge.eof());
         ++sit){
        if (sit->Unknowns.Exist(idx)){
            for (Uint i=0; i<vec->RowIdx->NumUnknownsVertex(); ++i)
            {
                counter++;
                srcFileEdge >> vec->Data[ sit->Unknowns(idx)+i ];
            }
        }
    }

    counter=0;


    // Close files
    srcFileVert.close();
    srcFileEdge.close();

} // END OF FUNCTIONs

} //end of namespace DROPS
