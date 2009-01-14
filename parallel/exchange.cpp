//**************************************************************************
// File:    exchange.cpp                                                   *
// Content: Classes that handle DDD-Interfaces for accumulations           *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:    May, 23th 2006                                                 *
//          - ExchangeBlockCL added                                        *
//          August, 25th 2006                                              *
//          - Sysnums on ohter procs are computed                          *
// Begin:   March, 08th 2006                                               *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file exchange.cpp

#include "parallel/exchange.h"
#include "parallel/parmultigrid.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <limits>

namespace DROPS{
// ------------------------------------
// E X C H A N G E  D A T A  C L A S S
// ------------------------------------
ExchangeDataCL::ExchangeDataCL(int proc) : toProc_(proc), SendType_(MPI_DATATYPE_NULL) { }

ExchangeDataCL::ExchangeDataCL(const ExchangeDataCL &ex)
{
    toProc_   = ex.toProc_;
    SendType_ = ex.SendType_;
    Sysnums_ = ex.Sysnums_;
#ifdef DebugParallelNumC
    SendTypeSize_=NoIdx;
#endif
}

/// \brief Destructor
ExchangeDataCL::~ExchangeDataCL()
{
    Sysnums_.resize(0);
    if (SendType_!=ProcCL::NullDataType)
        ProcCL::Free(SendType_);
}

/// \brief Set the corresponding proc
void ExchangeDataCL::SetToProc(const int Proc)
{
    toProc_ = Proc;
}

/// \brief Create Datatype for sending
void ExchangeDataCL::CreateDataType(const int count, const int blocklength[], const int array_of_displacements[])
{
    if (SendType_!=ProcCL::NullDataType)
        ProcCL::Free(SendType_);
    SendType_ = ProcCL::CreateIndexed<double>(count, blocklength, array_of_displacements);
    ProcCL::Commit(SendType_);
#ifdef DebugParallelNumC
    SendTypeSize_= array_of_displacements[count-1] + blocklength[count-1];
#endif
}

/// \brief Create describtion of the positions, where to store the recieved unknowns
void ExchangeDataCL::CreateSysnums(const SysnumListCT &sys)
{
    Sysnums_=sys;
}

void ExchangeDataCL::DebugInfo(std::ostream &os) const
{
    os << "ExchangeDataCL between " << ProcCL::MyRank() << " and " << toProc_ << std::endl
            << " Positions, where to store recieved numbers: " << std::endl;
    for (size_t i=0; i<Sysnums_.size(); ++i)
        os << Sysnums_[i] << " ";
    os << std::endl;
}

// --------------------------
// E X C H A N G E  C L A S S
// --------------------------

// Init of static members
ExchangeCL::SendList2ProcT ExchangeCL::SendList_     = ExchangeCL::SendList2ProcT();
ExchangeCL::RecvSysnumCT   ExchangeCL::RecvSysnums_  = ExchangeCL::RecvSysnumCT();
ExchangeCL::MappingIdxCT   ExchangeCL::tmpMappingIdx_= ExchangeCL::MappingIdxCT();
ExchangeCL::SysnumProcCT   ExchangeCL::tmpSysProc_   = ExchangeCL::SysnumProcCT();
IdxDescCL*               ExchangeCL::RowIdx_         = 0;
int                        ExchangeCL::maxNeighs_    = 0;

ExchangeCL::ExchangeCL()
{
    created_   = false;
    mapCreated_= false;
    tag_       = 1001;
}

ExchangeCL::~ExchangeCL()
{
    clear();
}

/// \brief remove all information
void ExchangeCL::clear()
{
    ExList_.clear();
    LocalIndex.resize(0);
    DistrIndex.resize(0);
    numLocalIdx_=0;
    numDistrIdx_=0;
    vecSize_=0;
    numNeighs_=0;
    SendRecvReq_.resize(0);
    created_=false;
    recvBuf_.resize(0);
    MappingIdx_.clear();
    SysProc_.clear();
    mapCreated_=false;
    accIdxCreated_=false;
}

/// \name Definition of the wrappers for DDD
//@{
extern "C" int HandlerGatherSysnumsVertexC(DDD_OBJ objp, void* buf){
    return ExchangeCL::HandlerGatherSysnums<VertexCL>(objp, buf);
}
extern "C" int HandlerScatterSysnumsVertexC(DDD_OBJ objp, void* buf){
    return ExchangeCL::HandlerScatterSysnums<VertexCL>(objp, buf);
}
extern "C" int HandlerGatherSysnumsEdgeC(DDD_OBJ objp, void* buf){
    return ExchangeCL::HandlerGatherSysnums<EdgeCL>(objp, buf);
}
extern "C" int HandlerScatterSysnumsEdgeC(DDD_OBJ objp, void* buf){
    return ExchangeCL::HandlerScatterSysnums<EdgeCL>(objp, buf);
}
//@}


void ExchangeCL::CreateExchangeDataMPIType()
/// Iterate over the SendList_ and create for each neighbor processor a single ExchangeDataCL
/// with information of elements in a vector to be sent.
/// Information about vertices and edges are handled in the same manner, hence the number of DOF on vertices
/// and edges must be the same.
{
    if (RowIdx_->NumUnknownsEdge()>0 && RowIdx_->NumUnknownsEdge()!=RowIdx_->NumUnknownsVertex()){
        throw DROPSErrCL("ExchangeCL::CreateExchangeDataClasses: Unknowns on edges (if exist) and vertices must be the same");
    }

    std::valarray<int> blocklength, displacements;
    for (const_SendList2ProcIter it(SendList_.begin()), end(SendList_.end()); it!=end; ++it){
        const size_t numSendUnks=RowIdx_->NumUnknownsVertex()*it->second.size();
        blocklength.resize(numSendUnks, RowIdx_->NumUnknownsVertex());
        displacements.resize(it->second.size());
        std::copy(it->second.begin(), it->second.end(), Addr(displacements));
        ExList_.push_back(ExchangeDataCL(it->first));
        ExList_.back().CreateDataType(it->second.size(), Addr(blocklength), Addr(displacements));
    }
}

void ExchangeCL::TransferSendOrder(bool CreateMap)
/// This function calls the DDD-Interface for vertices and edges to transfere
/// the send order. The interface of the InterfaceCL<VertexCL> and InterfaceCL<EdgeCL>
/// are used. The messages consits of a list of four numbers:
///  (1) the number of the sending processor
///  (2) the number of the receiving processor
///  (3) the position of the sysnum within the send process
///  (4) and the local sysnum
/// Unfortunately the DDD-Interface does not offer some information about the processor,
/// which gather or scatters. So the massage has a size of maxNeighs (over all processors) *4.
{
    numNeighs_= SendList_.size();
    maxNeighs_=GlobalMax((int)numNeighs_);

    // if there are unknowns on vertices
    if (RowIdx_->NumUnknownsVertex()>0){
        DDD_IFExchange(InterfaceCL<VertexCL>::GetIF(),  // exchange datas over distributed vertices
                       4*maxNeighs_*sizeof(IdxT),       // number of datas to be exchanged
                       HandlerGatherSysnumsVertexC,     // how to gather datas
                       HandlerScatterSysnumsVertexC     // how to scatter datas
                      );
    }
    if (RowIdx_->NumUnknownsEdge()>0){
        DDD_IFExchange(InterfaceCL<EdgeCL>::GetIF(),    // exchange datas over distributed edges
                       4*maxNeighs_*sizeof(IdxT),       // number of datas to be exchanged
                       HandlerGatherSysnumsEdgeC,       // how to gather datas
                       HandlerScatterSysnumsEdgeC       // how to scatter datas
                      );
    }

    // Static members are just used for the creation via DDD-Interface
    if (CreateMap){
        MappingIdx_ = tmpMappingIdx_;
        SysProc_.resize(tmpSysProc_.size());
        std::copy(tmpSysProc_.begin(), tmpSysProc_.end(), SysProc_.begin());
    }

#if DROPSDebugC&DebugParallelNumC
    // Check if all values of RecvSysnums_ are set
    for( RecvSysnumCT::const_iterator it(RecvSysnums_.begin()); it!=RecvSysnums_.end(); ++it){
        for (ExchangeDataCL::SysnumListCT::const_iterator sysit(it->second.begin()); sysit!=it->second.end(); ++sysit){
            if ( (*sysit)>=(int)RowIdx_->NumUnknowns()){
                std::cerr << "["<<ProcCL::MyRank()<<"] Found a strange sysnum "<<(*sysit)<<std::endl;
                throw DROPSErrCL("ExchangeCL::TransferSendOrder: Recieve sequence has not been set correct");
            }
        }
    }
#endif
}

// Create Indices, if forAccParDot is set, indices for performing inner products
// of two accumulated vectors are created too
void ExchangeCL::CreateIndices(IdxDescCL *RowIdx, const VectorBaseCL<bool>& DistSysnums, bool forAccParDot)
{
    // Count the local and distributed unknowns
    numLocalIdx_=0; numDistrIdx_=0;
    for (Ulint i=0; i<RowIdx->NumUnknowns(); ++i)
    {
        if (DistSysnums[i])
            ++numDistrIdx_;
        else
            ++numLocalIdx_;
    }
    // Allocate memory
    LocalIndex.resize(numLocalIdx_); DistrIndex.resize(numDistrIdx_);
    //  Create lists
    Ulint loc_pos=0, dist_pos=0;
    for (Ulint i=0; i<RowIdx->NumUnknowns(); ++i)
    {
        if (DistSysnums[i])
            DistrIndex[dist_pos++]=i;
        else
            LocalIndex[loc_pos++]=i;
    }

    // Falls gewuenscht, werden jetzt noch die Indices erstellt, so dass akkumulierte Vektoren behandelt werden
    // koennen. Dies ist im Prinzip die Umkehrung der akkumulierten in die verteilte Form.
    //
    // Die geteilten DoFs werden so aufgeteilt, dass für eine DoF genau ein Prozessor zustaendig ist. Dazu
    // sammelt jeder Prozessor die verteilten DoFs ein, falls er der Prozessor mit der kleinsten Id ist.
    /// \todo (of) mehr Dokumentation
    if (forAccParDot)
    {
        Assert(mapCreated_,
               DROPSErrCL("ExchangeCL::CreateIndices: If inner products for accumulated Vectors should be performed, in ExchangeCL::CreateList thet flag CreateMap must be set!"),
               DebugParallelNumC);
        AccDistIndex.clear();
        typedef std::vector<IdxT>           IdxVecT;
        typedef std::map<ProcNumT, IdxVecT> CoupProcIdxT;
        CoupProcIdxT coup;

        // Nachbar Prozessoren
        ProcNumCT neighs=GetNeighbors();
        for (ProcNumCT::const_iterator proc(neighs.begin()), pend(neighs.end()); proc!=pend; ++proc)
            coup[*proc];


        const ProcNumT me= ProcCL::MyRank();

        // Sortiere die zu sendenden DoF nach Prozessor in coup ein
        for (IdxT i=0; i<numDistrIdx_; ++i)
        {
            IdxT dof= DistrIndex[i];
            ProcNumCT onProc= GetProcs(dof);
            ProcNumT minProc= *std::min_element(onProc.begin(), onProc.end());
            if (minProc>me)
                coup[minProc].push_back(dof);
        }

        // Anlegen der Sendefelder
        const Uint numSends=coup.size();
        std::vector<ProcCL::RequestT> req(numSends);
        std::vector<VectorBaseCL<IdxT> > sendBuf(numSends, VectorBaseCL<IdxT>(1));
        int sendpos=0;

        // Speichere eigene DoFs und sende fremde DoFs
        for (CoupProcIdxT::const_iterator it(coup.begin()), end(coup.end()); it!=end; ++it){
            // Die erste Hälfte der Indices, die mit einem anderen Proc geteilt werden, übernimmt der Proc selber, die anderen werden gesendet
            const Uint size=it->second.size();
            const Uint count=size/2;
            sendBuf[sendpos].resize(count);
            int pos=0;

            // Indices, für die dieser Proc zuständig ist
            for (Uint i=0; i<count+(size%2); ++i){
                AccDistIndex.push_back(it->second[i]);
            }

            // Indices, für die der andere Proc zuständig ist
            for (Uint i=count+(size%2); i<size; ++i){
                sendBuf[sendpos][pos++] = GetExternalIdxFromProc(it->second[i],it->first);
            }

            // Senden
            req[sendpos] = ProcCL::Isend(Addr(sendBuf[sendpos]), pos, it->first, tag_+4);

            // nächster Prozessor
            ++sendpos;
        }


        // Empfangen der DoFs, für die dieser Proc zuständig ist
        VectorBaseCL<IdxT> recvBuf;
        for (ProcNumCT::const_iterator proc(neighs.begin()), pend(neighs.end()); proc!=pend; ++proc){
            if (*proc<me){
                ProcCL::StatusT stat;
                ProcCL::Probe(*proc, tag_+4, stat);
                const int count = ProcCL::GetCount<Ulint>(stat);
                recvBuf.resize(count);
                ProcCL::Recv(recvBuf, *proc, tag_+4);

                for (int i=0; i<count; ++i){
                    AccDistIndex.push_back(recvBuf[i]);
                }
            }
        }

        // Warten bis alle MPI-Sends abgeschlossen ist, damit die Felder gelöscht werden können
        ProcCL::WaitAll(req);

        // Setzen der Flags
        accIdxCreated_=true;
    }
}

/// \brief Create all information, to perform exchange of numerical data
void ExchangeCL::CreateList(const MultiGridCL& mg, IdxDescCL *RowIdx, bool CreateMap, bool CreateAccDist)
/// In order to handle numerical data, all processors must know something about
/// the distributed sysnums. I. e. which sysnums are stored local and which are
/// distributed. And if they do have a distributed sysnum, they have to know,
/// with wich other processors, they share this sysnum. All neccessary lists
/// and data structures are set up by this function.
/// \param mg        the multigrid
/// \param RowIdx    a pointer to the indix of the dof (e.g. velocity, pressure, ...)
/// \param CreateMap This function can create a mapping (localsysnum, proc)->external sysnum. Set this
///    this flag, to create such a mapping. This map is also used to perform inner products on accumulated
///    vectors.
/// \param CreateAccDist Tell this function to create also lists, so that the ExchangeCL is able to
///    handle accumulated vectors
{
    if (created_)
        clear();
    SendList_.clear();

    if (ProcCL::Size()==1){
        std::cerr << "Skipping CreateExchangeCL, because only one proc is involved!\n";
        created_    = true;
        accIdxCreated_=true;
        mapCreated_=false;
        vecSize_    = RowIdx->NumUnknowns();
        numLocalIdx_= vecSize_;
        numDistrIdx_= 0;
        numAllRecvUnk_=0;
        LocalIndex.resize(numLocalIdx_);
        AccDistIndex.resize(0);
        for (Uint i=0; i<numLocalIdx_; ++i)
            LocalIndex[i]=i;
        numNeighs_=0;
        return;
    }

    RowIdx_=RowIdx;
    Uint numUnkVert=RowIdx_->NumUnknownsVertex(), numUnkEdge=RowIdx_->NumUnknownsEdge();
    Uint numUnk=RowIdx_->NumUnknowns();

      // Allocate memory
    VectorBaseCL<bool> DistSysnums(false, numUnk);
    tmpSysProc_.resize(numUnk);

      // Collect distributed sysnums and put them into the lists
    if (numUnkVert>0){
        CollectSendSysNums(mg.GetTriangVertexBegin(), mg.GetTriangVertexEnd(), DistSysnums);
    }
    if (numUnkEdge>0){
        CollectSendSysNums(mg.GetTriangEdgeBegin(), mg.GetTriangEdgeEnd(), DistSysnums);
    }

      // sort the sysnums, to create the ExchangeDataCLs
    for (SendList2ProcIter it(SendList_.begin()), end(SendList_.end()); it!=end; ++it){
        std::sort( it->second.begin(), it->second.end() );
    }
      // Create the ExchangeDataCLs
    CreateExchangeDataMPIType();

      // Since number of send DOF to a processor is the same as the number of receiving
      // elements, resize the RecvSysnums_
    for (SendList2ProcIter it(SendList_.begin()), end(SendList_.end()); it!=end; ++it){
        const size_t numSendUnks=numUnkVert*it->second.size();
        RecvSysnums_.insert( RecvSysnumElemT(it->first, ExchangeDataCL::SysnumListCT(numSendUnks, NoIdx)) );
    }

      // Tell neighbor processors about the send order of the unknowns via a
      // DDD-interface.
    TransferSendOrder(CreateMap);
    if (!CreateMap){
        MappingIdx_.clear();
        SysProc_.clear();
        mapCreated_=false;
    }
    else{
        mapCreated_=true;
    }

      // Tell DataExchangeCL about the recieve sequence
    for (CommListCT::iterator it(ExList_.begin()), end(ExList_.end()); it!=end; ++it){
        it->CreateSysnums( RecvSysnums_[it->GetProc()] );
    }

      // Create the local and distributed indices
    CreateIndices(RowIdx, DistSysnums, CreateAccDist);

      // Set flags for error preventing
    created_ = true;
    accIdxCreated_=CreateAccDist;
    vecSize_ = RowIdx->NumUnknowns();
    numNeighs_= ExList_.size();

      // Free memory and reset static members
    SendList_.clear();
    RecvSysnums_.clear();
    tmpMappingIdx_.clear();
    tmpSysProc_.clear();
    RowIdx_=0;
    maxNeighs_=-1;

      // Allocate memory for receiving and set offsets
    SendRecvReq_.resize(2*numNeighs_);
    recvOffsets_.resize(numNeighs_);
    recvOffsets_[0]=0;
    Uint i=1;
    for (CommListCT::iterator it(ExList_.begin()), end(ExList_.end()); it!=end; ++it, ++i){
        if (i<numNeighs_){
            recvOffsets_[i] = recvOffsets_[i-1] + it->GetNumRecvEntries();
        }
        else{
           numAllRecvUnk_ =  recvOffsets_[i-1] + it->GetNumRecvEntries();
        }
    }
    recvBuf_.resize(numAllRecvUnk_);
}

/// \brief Debug information about the distributed and shared indices of a vector
void ExchangeCL::DebugInfo(std::ostream &os) const
{
    os << " Vectors, that can be handled by this class must have dimension: " << vecSize_ << std::endl;
    os << " Indices of local sysnums ("<<LocalIndex.size()<<"):\n";
    for (size_t i=0; i<GetNumLocIdx(); ++i)
        os << LocalIndex[i] << "  ";
    os << "\n Indices of distributed sysnums ("<<DistrIndex.size()<<"):\n";
    for (size_t i=0; i<GetNumDistIdx(); ++i)
        os << DistrIndex[i] << " ";
    os << "\n Exchnages between Procs:" <<std::endl;
    for (CommListCT::const_iterator it(ExList_.begin()), end(ExList_.end()); it!=end; ++it)
    {
        it->DebugInfo(os);
        os << std::endl;
    }
}

/// \brief Information about the size of sended data between proc
/** This function requieres that all procs call this function*/
void ExchangeCL::SizeInfo(std::ostream &os, const int Proc) const
/// \param[in] Proc proc, that print the information
/// \param os Where to print the information
{
    const int me   = ProcCL::MyRank(),          // who I am
              size = ProcCL::Size();            // how many are out there
    int * SendSize = new int [size-me-1];       // Information are symmetric, so store only message size of procs with greater proc-number than I have


    for (int i=0; i<size-me-1; ++i)             // reset to zero
        SendSize[i]=0;

    // Collect message size to all procs of greater proc-number that "me"
    for (CommListCT::const_iterator it(ExList_.begin()), end (ExList_.end()); it!=end; ++it)
        if (it->toProc_>me)
            SendSize[it->toProc_-me-1] = it->Sysnums_.size();

    if (me!=Proc){       // Send my information to the proc, that makes the output
        ProcCL::RequestT req=ProcCL::Isend(SendSize, size-me-1, Proc, tag_);
        ProcCL::Wait(req);
    }
    else
    {
        DMatrixCL<double> MsgSize(size,size);   // matrix of the message-sizes
        for (int i=0; i<size; ++i)
            for (int j=0; j<size;++j)
                MsgSize(i,j)=0;                 // set size to zero
        int *recvBuf = new int[size-1];
        for (int p=0; p<size; ++p)
        {
            if (p!=Proc)                        // recieve information from the other procs
                ProcCL::Recv(recvBuf, size-1-p, p, tag_);
            for (int i=0; i<size-1-p; ++i)      // put information into the matrix
                MsgSize(p+1+i,p) = MsgSize(p,p+1+i) = (p!=Proc) ? recvBuf[i] : SendSize[i];
        }
        // display informatrion
        os.precision(3);
        os.setf(std::ios::fixed);
        os << "Message Size between procs (in kb)\n\t";
        for (int j=0; j<size; ++j)
            os << std::setw(6) << j << ' ';
        os << std::endl;
        for (int i=0; i<size; ++i)
        {
            os << std::setw(7) << i;
            for (int j=0; j<size; ++j)
                os << std::setw(7) << (double(MsgSize(i,j))/128.);
            os << std::endl;
        }
        delete[] recvBuf;               // free memory
        os.precision(6);                // reset precision to standard
    }
    delete[] SendSize;                  // free memory
}


/// \brief Checks if a vector is accumulated
bool ExchangeCL::IsAcc(const VectorCL& vec) const
/// All vectors that have the same DoF must have the same value
/// \pre  The Mapping (my index,toProc) -> external index have to be created
{
    if (ProcCL::Size()==1)
        return true;
    Assert(mapCreated_, DROPSErrCL("ExchangeCL::IsAcc: Mapping has to be created before!"),DebugParallelNumC);
    bool ret=true;
    typedef std::pair<IdxT,double>                   CoupT;         // Coupling of DoF and value
    typedef std::map<ProcNumT, std::vector<CoupT> >  CoupCT;        // container for procs
    CoupCT ToSend;
    ProcNumCT neigh=GetNeighbors();

      // Collect information of distributed DoFs
    for (IdxT i=0; i<vec.size(); ++i)
        if (IsDist(i)){
            ProcNumCT onProcs=GetProcs(i);
            for (ProcNumCT::iterator it(onProcs.begin()), end(onProcs.end()); it!=end; ++it)
                ToSend[*it].push_back(CoupT(i,vec[i]));
        }


     // buffer for sending
    IdxT   **senddof_buf = new IdxT*[neigh.size()];
    double **sendval_buf = new double*[neigh.size()];
    std::valarray<ProcCL::RequestT> req(2*neigh.size());
    int proc_pos=0;
    for (CoupCT::iterator it(ToSend.begin()), end(ToSend.end()); it!=end; ++it, ++proc_pos){
        // allocate mem for sending
        const int count = it->second.size();
        senddof_buf[proc_pos] = new IdxT[count];
        sendval_buf[proc_pos] = new double[count];
        for (int i=0; i<count; ++i){
            senddof_buf[proc_pos][i]=GetExternalIdxFromProc(it->second[i].first,it->first);
            sendval_buf[proc_pos][i]=it->second[i].second;
        }
        // send
        req[proc_pos+0]= ProcCL::Isend(senddof_buf[proc_pos], count, it->first, tag_);
        req[proc_pos+1]= ProcCL::Isend(sendval_buf[proc_pos], count, it->first, tag_+1);
    }

    // recieve and check
    for (ProcNumCT::iterator proc(neigh.begin()), end(neigh.end()); proc!=end; ++proc){
        ProcCL::StatusT stat;
        ProcCL::Probe(*proc, tag_, stat);
        const int count = ProcCL::GetCount<IdxT>(stat);
        IdxT *recvdof_buf  = new IdxT[count];
        double *recvval_buf= new double[count];
        ProcCL::Recv(recvdof_buf, count, *proc, tag_);
        ProcCL::Recv(recvval_buf, count, *proc, tag_+1);
        for (int i=0; i<count; ++i){
            if (std::fabs(recvval_buf[i]-vec[recvdof_buf[i]])>1e-10){
                ret=false;
                  std::cerr << "["<<ProcCL::MyRank()<<"] Differ at Pos " << recvdof_buf[i]
                            << " by loc " <<vec[recvdof_buf[i]]<<" ext "<<recvval_buf[i]
                            << " on proc "<<*proc<<std::endl;

            }
        }
        delete[] recvdof_buf;
        delete[] recvval_buf;
    }

    ProcCL::WaitAll(req);
    for (Uint i=0; i<neigh.size(); ++i){
        delete[] senddof_buf[i];
        delete[] sendval_buf[i];
    }
    delete[] senddof_buf;
    delete[] sendval_buf;
    return Check(ret);
}

/// \brief Check if two ExchangeCL's seems to be the same
/** Check equality by: number of local and distributed indices, size of handleable vectors,
    same local and distributed indices.<br>
    Then check included ExchangeDataCL's for equality by: proc number, number and position of recieving sysnums*/
bool ExchangeCL::IsEqual(const ExchangeCL &ex, std::ostream* os) const
{
    if (numLocalIdx_ != ex.numLocalIdx_){
        if (os)
            (*os) << "["<<ProcCL::MyRank()<<"] number of local indices is not equal("<<numLocalIdx_<<"!="<<ex.numLocalIdx_<<")!" << std::endl;
        return false;
    }
    if (numDistrIdx_ != ex.numDistrIdx_){
        if (os)
            (*os) << "["<<ProcCL::MyRank()<<"] number of local distributed is not equal("<<numDistrIdx_<<"!="<<ex.numDistrIdx_<<")!" << std::endl;
        return false;
    }
    if (vecSize_!= ex.vecSize_){
        std::cerr << "["<<ProcCL::MyRank()<<"] vecSize is not equal!" << std::endl;
        return false;
    }

    for (size_t i=0; i<numLocalIdx_; ++i)
        if (LocalIndex[i] != ex.LocalIndex[i]){
            if (os)
                (*os) << "["<<ProcCL::MyRank()<<"] LocalIndex["<<i<<"] is not equal!"<<std::endl;
            return false;
        }

    for (size_t i=0; i<numDistrIdx_; ++i)
        if (DistrIndex[i]!=ex.DistrIndex[i]){
            if (os)
                (*os)  << "["<<ProcCL::MyRank()<<"] DistrIndex["<<i<<"] is not equal!"<<std::endl;
            return false;
        }

    for (CommListCT::const_iterator ex_it(ex.ExList_.begin()), ex_end(ex.ExList_.end()),
         my_it(ExList_.begin()),    my_end(ExList_.end());
         ex_it!=ex_end && my_it!=my_end; ++my_it, ++ex_it){
        // Check, ob ExchangeDataCL gleich ist!
        if (my_it->GetProc() != ex_it->GetProc()){
            if (os)
                (*os)  <<"["<<ProcCL::MyRank()<<"] Exchanging Proc is not equal!"<<std::endl;
            return false;
        }
        if (my_it->Sysnums_.size() != ex_it->Sysnums_.size()){
            if (os)
                (*os)  << "["<<ProcCL::MyRank()<<"] Size of Synum="<<my_it->Sysnums_.size()<<"!="<<ex_it->Sysnums_.size()<<std::endl;
            return false;
        }
        for (size_t i=0; i<my_it->Sysnums_.size(); ++i)
            if (my_it->Sysnums_[i] != ex_it->Sysnums_[i])
            {
                if (os)
                    (*os)  << "["<<ProcCL::MyRank()<<"] Sysnum["<<i<<"]="<<my_it->Sysnums_[i]<<"!="<<ex_it->Sysnums_[i]<<std::endl;
                return false;
            }
    }
    return true;
}

// -------------------------------------
// E X C H A N G E  B L O C K  C L A S S
// -------------------------------------

///\brief Construct a class that can store exchange information for m blocks
ExchangeBlockCL::ExchangeBlockCL(size_t m) : m_(m), ExchgList_(m), Block_(m+1)
{
    created_      = false;
    blockCreated_ = false;
    start_tag_    = 1001;
}

/// \brief Remove all information
void ExchangeBlockCL::clear()
{
    created_      = false;
    blockCreated_ = false;
    for (ExchangeCT::iterator it(ExchgList_.begin()), end(ExchgList_.end()); it!=end; ++it)
        it->clear();
}

/// \brief Calculate and store indices, where the blocks starts
void ExchangeBlockCL::CreateBlockIdx()
{
    Assert(created_ && !blockCreated_, DROPSErrCL("ExchangeBlockCL::CreateBlockIdx: not all index-describer classes recieved, or block_idx has been created allready"), DebugNumericC);
    Block_[0] = 0;
    for (size_t i=0; i<m_; ++i)
        Block_[i+1] = Block_[i] + ExchgList_[i].GetNum();
    blockCreated_ = true;
   // Allocate memory for requests
    SendRecvReq_.resize(GetNumBlocks());
    if (ProcCL::Size()>1)
        for (size_t i=0; i<GetNumBlocks(); ++i)
            SendRecvReq_[i].resize(2*ExchgList_[i].GetNumNeighs());

}

/// \brief Create for a single index-describer an exchange structure
void ExchangeBlockCL::CreateList(const MultiGridCL& mg, size_t i, MLIdxDescCL *idx, bool CreateMap, bool CreateAccDist)
/// \param[in] mg Multigrid of that the ExchangeCL should be build
/// \param[in] i   index of the block to create an ExchangeCL for that block
/// \param[in] idx row-index describer, that describes the dependences between num and geom
/// \param[in} CreateMap Create Mapping of indices
/// \param[in] CreateAccDist Create indices for AccParDot
/** \todo (of): Bis jetzt wird für jeden Index eine Nachricht geschickt. Es wäre besser, wenn
        zwischen zwei Prozessoren genau eine Nachricht geschickt werden würde! */
{
    Assert(i<m_, DROPSErrCL("ExchangeBlockCL::CreateList: block does not exist"), DebugNumericC);
    ExchgList_[i].CreateList(mg,idx, CreateMap, CreateAccDist);
    Created();              // Check if all indices have been submitted. If all indices have been recieved, create  block-indices!
}

/// \brief Create for all index-describers an single exchange structure
void ExchangeBlockCL::CreateList(const MultiGridCL& mg, const MLIdxDescCT &idxs, bool CreateMap, bool CreateAccDist)
/// \param[in] mg Multigrid of that the ExchangeCL should be build
/// \param[in] idxs      std::vector of pointers to index-describer classes. i-th element of this vector should describe the i-th block
/// \param[in] CreateMap Create Mapping of indices
/// \param[in] CreateAccDist Create indices for AccParDot
{
    Assert(idxs.size()==m_, DROPSErrCL("ExchangeBlockCL::CreateList: given idx-describers does not fit"), DebugNumericC);
    for (size_t i=0; i<m_; ++i)
        CreateList(mg,i,idxs[i],CreateMap,CreateAccDist);
    Created();
    Assert(created_ && blockCreated_, DROPSErrCL("ExchangeBlockCL::CreateList: Internal Error, lists not created!"), DebugNumericC);
}
} // end of namespace DROPS
