//**************************************************************************
// File:    exchange.tpp                                                   *
// Content: Classes that handle DDD-Interfaces for accumulations           *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:    May, 23th 2006                                                 *
//          - ExchangeBlockCL added
// Begin:   March, 08th 2006                                               *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file exchange.tpp

namespace DROPS{

// ------------------------------------
// E X C H A N G E  D A T A  C L A S S
// ------------------------------------
/// \brief Get the corresponding proc
int ExchangeDataCL::GetProc() const
{
    return toProc_;
}

/// \brief Get number of received elements
size_t ExchangeDataCL::GetNumRecvEntries() const
{
    return Sysnums_.size();
}

/// \brief Send data
ProcCL::RequestT ExchangeDataCL::Isend(const VectorCL& v, int tag, Ulint offset) const
/** This procedure send the datas with an non-blocking non-synchronous MPI Send.
    \param v      vector, that contains elements for sending (local data)
    \param tag    used tag for communication
    \param offset start element of the vector
    \pre v has to be big enough
*/
{
    Assert(v.size()>=SendTypeSize_+offset, DROPSErrCL("ExchangeDataCL::Isend: Vector is not long enough for transfer!"), DebugParallelNumC);
    return ProcCL::Isend(Addr(v)+offset, 1, SendType_, toProc_, tag);
}

/// \brief Recieve data
ProcCL::RequestT ExchangeDataCL::Irecv(int tag, VectorCL& recvBuf, Ulint offset) const
/** This procedure receives the datas with an non-blocking non-synchronous MPI
    Receive.
    \param tag     used tag for communication
    \param recvBuf buffer for storing received unknowns
    \param offset  first position, where to store unknowns
    \pre recvBuf has to be big enough to store all data
    \pre no other procedures is allowed to work on the memory (in particular no other Irecv!)
*/
{
    Assert(recvBuf.size()>=GetNumRecvEntries(), DROPSErrCL("ExchangeDataCL::Irecv: Receive buffer is not long enough!"), DebugParallelNumC);
    return ProcCL::Irecv(Addr(recvBuf)+offset, GetNumRecvEntries(), toProc_, tag);
}

/// \brief Add data
void ExchangeDataCL::Accumulate(VectorCL& v, Ulint offsetV, VectorCL& recvBuf, Ulint offsetRecv) const
/** This procedure accumulates received data. It assumes, that the data has been
    received.
    \param v          original value, that contains all local unknowns
    \param offsetV    start element of the vector, that contains local unknowns
    \param recvBuf    vector of all received elements
    \param offsetRecv first element in receive buffer
    \pre Communication has to be done before entering this procedure
*/
{
    // add the datas to the positions described by Sysnums_
    for (Uint i=0; i<GetNumRecvEntries(); ++i)
        v[Sysnums_[i]+offsetV] += recvBuf[i+offsetRecv];
}


// --------------------------
// E X C H A N G E  C L A S S
// --------------------------

template <typename SimplexIterT>
void ExchangeCL::CollectSendSysNums(const SimplexIterT& begin, const SimplexIterT& end,
                                    VectorBaseCL<bool>& DistSysnums)
/// Iterate over all simplices between begin and end and put the sysnum (if distributed
/// sysnums exits on the simplex) into the SendList_ and set flag in DistSysnum, that
/// the sysnum is distributed.
{
    Uint idx(RowIdx_->GetIdx());                                                    // index of the unknowns
    for (SimplexIterT sit(begin); sit!=end; ++sit){                                 // for all simplices
        if (!sit->IsLocal() && sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){  // check for distributed sysnum
            for (int *procList=sit->GetProcList(); *procList!=-1; procList+=2){     // for all processors, that owns the simplex too
                if ( *procList!=ProcCL::MyRank() && *(procList+1)==PrioHasUnk ){    // if other processor stores unknowns
                    SendList_[*procList].push_back(sit->Unknowns(idx));
                    for (Uint i=0; i<RowIdx_->NumUnknownsVertex(); ++i){             // remember, which sysnums are distributed
                        DistSysnums[sit->Unknowns(idx)+i]= true;
                    }
                }
            }
        }
    }
}

/// \brief Find the position of an element in a sorted vector.
template <typename T>
  IdxT ExchangeCL::findPos(const std::vector<T>& a, const T& elem)
/// Find the position of an element in a sorted vector. If the element is
/// not found, return NoIdx
{
    typename std::vector<T>::const_iterator it= std::lower_bound(a.begin(), a.end(), elem);
    if (it==a.end() || *it!=elem)
        return NoIdx;
    return (IdxT)std::distance(a.begin(), it);
}


template <typename SimplexT>
  int ExchangeCL::HandlerGatherSysnums(DDD_OBJ objp, void* buf)
/// Gather sendposition on sender-side. Therefore iterate over the list of
/// sending sysnums and put them into the buffer. The content of the buffer
/// is described detailed in the documentation of TransferSendOrder.
/// (Should be private! Just public, so DDD can call this function)
/// \param objp pointer to the simplex
/// \param buf buffer
{
    SimplexT* const sp = ddd_cast<SimplexT*>(objp);
    IdxT* buffer= static_cast<IdxT*>(buf);
    Uint idx=RowIdx_->GetIdx();
    int pos=0;

    for (const_SendList2ProcIter it(SendList_.begin()), end(SendList_.end()); it!=end; ++it, ++pos){
        *(buffer++)= (IdxT)ProcCL::MyRank();                        // from processor
        *(buffer++)= (IdxT)it->first;                               // to processor
        // unknown exist on simplex and processor
        if (sp->Unknowns.Exist() && sp->Unknowns.Exist(idx)){
            *(buffer++)= findPos(it->second, sp->Unknowns(idx));    // send position
            *(buffer++)= sp->Unknowns(idx);                         // local sysnum
        }
        else{
            *(buffer++)=NoIdx;
            *(buffer++)=NoIdx;
        }
        Assert(pos<maxNeighs_, DROPSErrCL("ExchangeCL::HandlerGatherSysnums: To many neigh processors"), DebugParallelNumC);
    }
    // Fill the rest of the buffer with dummy-values:
    for (; pos<maxNeighs_; ++pos){
        *(buffer++)= (IdxT)ProcCL::Size();
        *(buffer++)= (IdxT)ProcCL::Size();
        *(buffer++)= NoIdx;
        *(buffer++)= NoIdx;
    }
    return 0;
}

template <typename SimplexT>
  int ExchangeCL::HandlerScatterSysnums(DDD_OBJ objp, void* buf)
/// Scatter sendposition on reciever-side. Therefore iterate over the
/// content of the buffer and check for the right content for this processor.
/// The content of the message is described detailed in the documentation of
/// TransferSendOrder. (Should be private! Just public, so DDD can call this function)
/// \param objp pointer to the simplex
/// \param buf buffer
{
    SimplexT* const sp = ddd_cast<SimplexT*>(objp);
    IdxT* buffer= static_cast<IdxT*>(buf);
    Uint idx=RowIdx_->GetIdx();
    int  fromProc=-1;
    IdxT sendPos=NoIdx;
    IdxT localSysnum, remoteSysnum;
    IdxT myRank=(IdxT)ProcCL::MyRank();

    // Check if there are unknowns on the simplex and processor
    if (sp->Unknowns.Exist() && sp->Unknowns.Exist(idx)){
        for (int i=0; i<maxNeighs_; ++i){
            // Check if this is the correct processor
            if ( buffer[4*i+1]==myRank){
                fromProc= buffer[4*i];                              // sending processor
                if (buffer[4*i+2]!=NoIdx){
                    sendPos = buffer[4*i+2]*RowIdx_->NumUnknownsVertex(); // position of DOF in messages
                    remoteSysnum=buffer[4*i+3];                           // sysnum on sending processor
                    localSysnum=sp->Unknowns(idx);                        // local sysnum
                    // put these information into the lists
                    for (Uint j=0; j<RowIdx_->NumUnknownsVertex(); ++j){
                        RecvSysnums_[fromProc][sendPos+j]=localSysnum + j;          // where to add received dof
                        tmpMappingIdx_[fromProc][localSysnum+j]= remoteSysnum +j;   // mapping (proc,localsysnum)->remote sysnum
                        tmpSysProc_[localSysnum+j].push_back(fromProc);             // mapping localsysnum->[processors, owning sysnum]
                    }
                }
                else{
                    throw DROPSErrCL("ExchangeCL::HandlerScatterSysnums: Recieved NoIdx as sendposition!");
                }
            }
        }
    }
    return 0;
}

/// \brief Get the number of local sysnums
Ulint ExchangeCL::GetNumLocIdx()  const
{
    Assert(created_, DROPSErrCL("ExchangeCL::GetNumLocIdx: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    return numLocalIdx_;
}

/// \brief Get the number of ditributed sysnums
Ulint ExchangeCL::GetNumDistIdx() const
{
    Assert(created_, DROPSErrCL("ExchangeCL::GetNumDistIdx: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    return numDistrIdx_;
}

/// \brief get the size of vector, that can be accumulated
Ulint ExchangeCL::GetNum() const
{
    Assert(created_, DROPSErrCL("ExchangeCL::GetNum: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    return vecSize_;
}

/// \brief Get number of elements, that should be received (i.e. size of the receive buffer)
Ulint ExchangeCL::GetNumReceiveElements() const
{
    return numAllRecvUnk_;
}

/// \brief check if the list has been created
bool  ExchangeCL::Created() const{
    return created_;
}
/// \brief check if the mapping has been created
bool  ExchangeCL::MapCreated() const
{
    return mapCreated_;
}

/// \brief check if index for accumulated inner products are set
bool  ExchangeCL::AccIdxCreated() const
{
    return accIdxCreated_;
}

/// \brief get number of neighbor processors
Uint ExchangeCL::GetNumNeighs() const
{
    return numNeighs_;
}

/// \brief Start communication
void ExchangeCL::InitCommunication(const VectorCL &vec, RequestCT& req, int tag, Ulint offset, VectorCL* recvBuf) const
/** This procedure initializes the communication with neighbor processors, i.e.
    it calls Isend and Irecv to all neighbor processors.
    \param vec     all entries of vec, that are shared, are send to other procs, that owns this unknown too. Watch out, not to change these entries bevore recieving!
    \param req     Requests that should be used for sending and receiving
    \param tag     default -1: the tagused by this function. If tag=-1, then use the default tag=1001
    \param offset  default 0: For blocked vectors this offset is used to enter a special block
    \param recvBuf default 0: Buffer for receiving unknowns (default 0: using the standard receive buffer of this class)
    \pre List has to be created
    \pre request container has to be the size of (2*number of neighbors)
    \pre recvBuf must be big enough (if given) or standard receive buffer must be big enough (if not given, default)
*/
{
    Assert(created_,
           DROPSErrCL("ExchangeCL::InitCommunication: Lists have not been created (Maybe use CreateList before!\n"),
           DebugParallelNumC);
    Assert(req.size()==2*GetNumNeighs(),
           DROPSErrCL("ExchangeCL::InitCommunication: Request container has wrong length"),
           DebugParallelNumC);
    Assert((recvBuf==0 && recvBuf_.size()>=numDistrIdx_) || (recvBuf!=0 && recvBuf->size()>=numDistrIdx_),
           DROPSErrCL("ExchangeCL::InitCommunication: Receive Buffer is not big enough"),
           DebugParallelNumC);

    // if no receive buffer is explicitly given, use the default one
    if (!recvBuf)
        recvBuf= &recvBuf_;

    // set tag for sending and receiving
    const int mytag= (tag==-1) ? tag_ : tag;

    // iterate over all neighbors and init communication
    CommListCT::const_iterator lit=ExList_.begin(), end=ExList_.end();
    const size_t num_neigh= GetNumNeighs();
    Uint i=0;
    for (; lit!=end; ++lit, ++i){
        req[i]          = lit->Isend(vec, mytag, offset);
        req[i+num_neigh]= lit->Irecv(mytag, *recvBuf, recvOffsets_[i]);
    }
}

/// \brief Accumulate the Vector
void ExchangeCL::AccFromAllProc(VectorCL &vec, RequestCT& req, Ulint offset, VectorCL* recvBuf) const
/** This procedure waits until all (I)sends from this proc into the wild
    proc-world and all (I)receive operations are finished. Then the accumulation
    can be performed.

    \param vec     recieve from other procs the values of the shared entries and add them
    \param tag     default -1: the tag used by this function. If tag=-1, then use the default tag=1001
    \param offset  default 0: For blocked vectors this offset is used to enter a special block
    \param recvBuf default 0: Buffer of received unknowns (default 0: using the standard receive buffer of this class)
    \pre List has to be created
    \pre request container has to be the size of (2*number of neighbors)
    \pre recvBuf must be big enough (if given) or standard receive buffer must be big enough (if not given, default)
*/
{
    Assert(created_,
           DROPSErrCL("ExchangeCL::AccFromAllProc: Lists have not been created (Maybe use CreateList before!\n"),
           DebugParallelNumC);
    Assert(req.size()==2*GetNumNeighs(),
           DROPSErrCL("ExchangeCL::InitCommunication: Request container has wrong length"),
           DebugParallelNumC);
    Assert((recvBuf==0 && recvBuf_.size()>=numDistrIdx_) || (recvBuf!=0 && recvBuf->size()>=numDistrIdx_),
           DROPSErrCL("ExchangeCL::InitCommunication: Receive Buffer is not big enough"),
           DebugParallelNumC);

    // if no receive buffer is explicitly given, use the default one
    if (!recvBuf)
        recvBuf= &recvBuf_;

    // Wait untill all sends and revceives are completed
    ProcCL::WaitAll(req);

    // start accumulation
    CommListCT::const_iterator lit=ExList_.begin(), end=ExList_.end();
    Uint i=0;
    for (; lit!=end; ++lit, ++i)
        lit->Accumulate(vec, offset, *recvBuf, recvOffsets_[i]);
}

/****************************************************
*   L O C  D O T _                                  *
*****************************************************
*  inner product of two distributed vectors x and   *
*  y without performing an global reduce. The       *
*  vector x will be accumulated and stored in x_acc *
****************************************************/
double ExchangeCL::LocDot_(const VectorCL& x, const VectorCL& y, VectorCL* x_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeCL::LocDot_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::LocDot_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::LocDot_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

    bool newx = (x_acc==0);

    double loc_sum=0;                       // sum of local entries
    double dist_sum=0;                      // sum of distributed entries

    InitCommunication(x, SendRecvReq_);        // send shared entries of x to all neighbor procs

    // assign all values of x to x_acc
    if (newx)
        x_acc= new VectorCL(x);
    else{
        Assert(x_acc->size()==vecSize_, DROPSErrCL("ExchangeCL::LocDot_: vector x_acc has not the right length"),DebugParallelNumC);
        *x_acc=x;
    }

    for (Ulint i=0; i<numLocalIdx_; ++i)    // do local summation
        loc_sum += x[LocalIndex[i]] * y[LocalIndex[i]];

    AccFromAllProc(*x_acc, SendRecvReq_);   // recieve values from neighbors and sum them up

    for (Ulint i=0; i<numDistrIdx_; ++i)    // do summation of distributed entries
        dist_sum += (*x_acc)[DistrIndex[i]] * y[DistrIndex[i]];


    if (newx){                              // if new x_acc where created give memory free
        delete x_acc;
        x_acc=0;
    }

    return loc_sum+dist_sum;                // return result
}

/****************************************************
*   A C C U R  L O C  D O T  N O  A C C _           *
*****************************************************
*  inner product of two accumulated vectors x and   *
*  y without performing an global reduce.           *
*****************************************************
*  pre: both vectors are accumulated                *
****************************************************/
double ExchangeCL::AccurLocDotNoAcc_(const VectorCL& x, const VectorCL& y) const
{
    Assert(created_, DROPSErrCL("ExchangeCL::AccurLocDotNoAcc_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::AccurLocDotNoAcc_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::AccurLocDotNoAcc_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::AccurLocDotNoAcc_: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    double loc_sum=0,       // sum of local entries
           acc_sum=0;       // sum of distributed entries

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x[LocalIndex[i]] * y[LocalIndex[i]];


    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += x[AccDistIndex[i]] * y[AccDistIndex[i]];

    return loc_sum+acc_sum;
}

/****************************************************
*   A C C U R  L O C  D O T  O N E  A C C _         *
*****************************************************
*  Accurate inner product of one accumulated vector *
*  (x_acc) and one distributed vecor (y) without    *
*  performing an global reduce. If y_acc not equal  *
*  zero, the acumulated form of y will be given     *
*  by return in y_acc.                              *
*****************************************************
*  pre: x_acc is accumulated, y is distributed      *
****************************************************/
double ExchangeCL::AccurLocDotOneAcc_(const VectorCL& x_acc, const VectorCL& y, VectorCL* y_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeCL::AccurLocDotOneAcc_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x_acc.size()==y.size(), DROPSErrCL("ExchangeCL::AccurLocDotOneAcc_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x_acc.size()==vecSize_, DROPSErrCL("ExchangeCL::AccurLocDotOneAcc_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::AccurLocDotOneAcc_: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    bool newy = (y_acc==0);

    double loc_sum=0,                       // sum of local entries
           acc_sum=0;                       // sum of distributed entries

    InitCommunication(y, SendRecvReq_);     // send shared entries of x to all neighbor procs

    if (newy)
        y_acc = new VectorCL(y);
    else{
        Assert(y_acc->size()==x_acc.size(), DROPSErrCL("ExchangeCL::AccurLocDotOneAcc_: y_acc has not the right size"), DebugParallelNumC);
        *y_acc=y;
    }

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x_acc[LocalIndex[i]] * y[LocalIndex[i]];

    AccFromAllProc(*y_acc, SendRecvReq_);   // recieve values from neighbors and sum them up

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += x_acc[AccDistIndex[i]] * (*y_acc)[AccDistIndex[i]];

    // free memory
    if (newy){
        delete y_acc;
        y_acc=0;
    }

    return loc_sum+acc_sum;
}

/****************************************************
*   A C C U R  L O C  D O T  B O T H   A C C _      *
*****************************************************
*  Accurate inner product of two distributed        *
*  vectors (x and y) without performing an global   *
*  reduce. If x_acc or y_acc not equal zero, the    *
*  acumulated form of x or y will be given by return*
*  in x_acc or y_acc.                               *
*****************************************************
*  pre: x and y have distributed form               *
****************************************************/
double ExchangeCL::AccurLocDotBothAcc_(const VectorCL& x, const VectorCL& y, VectorCL* x_acc, VectorCL* y_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeCL::AccurLocDotBothAcc_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::AccurLocDotBothAcc_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::AccurLocDotBothAcc_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::AccurLocDotBothAcc_: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    bool newx = (x_acc==0);
    bool newy = (y_acc==0);

    double loc_sum=0,                           // sum of local entries
           acc_sum=0;                           // sum of distributed entries

    RequestCT req_y(2*ExList_.size());
    VectorCL  secondRecvBuf(recvBuf_.size());

    InitCommunication(x, SendRecvReq_, tag_);               // send shared entries of x to all neighbor procs
    InitCommunication(y, req_y, tag_+1, 0, &secondRecvBuf); // send shared entries of y to all neighbor procs

    if (newx)
        x_acc = new VectorCL(x);
    else{
        Assert(x_acc->size()==x.size(), DROPSErrCL("ExchangeCL::AccurLocDotBothAcc_: x_acc has not the right size"), DebugParallelNumC);
        *x_acc=x;
    }

    if (newy)
        y_acc = new VectorCL(y);
    else{
        Assert(y_acc->size()==x.size(), DROPSErrCL("ExchangeCL::AccurLocDotBothAcc_: y_acc has not the right size"), DebugParallelNumC);
        *y_acc=y;
    }

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x[LocalIndex[i]] * y[LocalIndex[i]];

    AccFromAllProc(*x_acc, SendRecvReq_);                   // recieve values from neighbors and sum them up
    AccFromAllProc(*y_acc, req_y, 0, &secondRecvBuf);       // recieve values from neighbors and sum them up

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += (*x_acc)[AccDistIndex[i]] * (*y_acc)[AccDistIndex[i]];

    // free memory
    if (newx){
        delete x_acc;
        x_acc=0;
    }

    if (newy){
        delete y_acc;
        y_acc=0;
    }

    return loc_sum+acc_sum;
}

/// \brief Local inner product without global reduce
/** This function computes the inner product of two vectors without global reduce. The vectors may be available
    in different forms. Therefore the \a acc_x and \a acc_y flag exists. If an accumulation is performed, the
    result can be given to the caller.*/
double ExchangeCL::LocDot (const VectorCL& x, bool acc_x,
                           const VectorCL& y, bool acc_y,
                           bool useAccur,
                           VectorCL* x_acc, VectorCL *y_acc) const
    /// \param[in]  x        first vector
    /// \param[in]  acc_x    is vector \a x given in accumulated form
    /// \param[in]  y        second vector
    /// \param[in]  acc_y    is vector \a y given in accumulated form
    /// \param[in]  useAccur should the accurater but slower version be used
    /// \param[out] x_acc    if not equal zero, pointer on the accumulated form of \a x (also non distributed values will be copied)
    /// \param[out] y_acc    if not equal zero, pointer on the accumulated form of \a y (also non distributed values will be copied)
    /// \return              inner product of \a x and \a y without global reduce
{
    Assert(created_, DROPSErrCL("ExchangeCL::LocDot: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::LocDot: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::LocDot: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

    if (acc_x && x_acc!=0)
        *x_acc=x;
    if (acc_y && y_acc!=0)
        *y_acc=y;

    if (useAccur){
        if (acc_x && acc_y)
            return AccurLocDotNoAcc_(x,y);
        if (acc_x && !acc_y)
            return AccurLocDotOneAcc_(x,y,y_acc);
        if (acc_y && !acc_x)
            return AccurLocDotOneAcc_(y,x,x_acc);
        if (!acc_x && !acc_y)
            return AccurLocDotBothAcc_(x,y,x_acc,y_acc);
    }
    else{
        if (!acc_x && !acc_y)
        {   // x and y are not accumulated, so accumulate one
            if (x_acc && y_acc==0)              // accumulate x
                return LocDot_(x,y,x_acc);
            if (y_acc && x_acc==0)              // accumulate y
                return LocDot_(y,x,y_acc);
            if (x_acc==0 && y_acc==0)           // no accumulated form is wished
                return LocDot_(x,y,0);
            // x and y should be accumulated. That makes no sence
            throw DROPSErrCL("ExchangeCL::LocDot: It makes no sence to do LocDot to accumulate both vectors. Better use accurate version");
        }
        // form here on at least one vector is accumulated
        if (acc_x && !acc_y)
            return dot(x,y);
        if (acc_y && !acc_x)
            return dot(x,y);
        // if both vectors are accumulated, you are not allowed to call this function with useAccur==false
        if (acc_x && acc_y)
            throw DROPSErrCL("ExchangeCL::LocDot: Cannot perform a normal inner product on two accumulated vectors, set useAccur=true");
    }
    throw DROPSErrCL("ExchangeCL::LocDot: Internal error, no matching found");
}

/// \brief Local inner product with global reduce
double ExchangeCL::ParDot (const VectorCL& x, bool acc_x,
                           const VectorCL& y, bool acc_y,
                           bool useAccur,
                           VectorCL* x_acc, VectorCL *y_acc) const
/** This function computes the inner product of two vectors. The vectors may be available in different forms.
    Therefore the \a acc_x and \a acc_y flag exists. If an accumulation is performed, the
    result can be given to the caller.*/
    /// \param[in]  x        first vector
    /// \param[in]  acc_x    is vector \a x given in accumulated form
    /// \param[in]  y        second vector
    /// \param[in]  acc_y    is vector \a y given in accumulated form
    /// \param[in]  useAccur should the accurater but slower version be used
    /// \param[out] x_acc    if not equal zero, pointer on the accumulated form of \a x
    /// \param[out] y_acc    if not equal zero, pointer on the accumulated form of \a y
    /// \return              inner product of \a x and \a y without global reduce
{
    return GlobalSum(LocDot(x, acc_x, y, acc_y, useAccur, x_acc, y_acc));
}

/// \brief Squared norm of a vector without global reduce
/** This function computes the squared norm of a Vector without global reduce. The Vector \a r can available in
    accumulated or distributed form (according to switch \a acc_r). If wished the accumulated form of \a r will
    is given on return within \a r_acc.*/
double ExchangeCL::LocNorm_sq(const VectorCL &r, bool acc_r, bool useAccur, VectorCL* r_acc) const
    /// \param[in] r        Vector of which the norm should be computed
    /// \param[in] acc_r    is vector \a r given in accumulated form
    /// \param[in] useAccur should the accurater but slower version be used
    /// \param[out] r_acc   if not equal zero, pointer on the accumulated form of \a r
    /// \return             norm of the vector \a r
{
    double loc_norm_sq;
    if (acc_r && useAccur){
        loc_norm_sq=LocAccNorm_sq(r);
        if (r_acc){
            if (r_acc->size()!=r.size())
                r_acc->resize(r.size());
            *r_acc =r;
        }
    }
    else if (!acc_r && useAccur){
        if (r_acc && r_acc->size()!=r.size())
            r_acc->resize(r.size());

        bool new_tmp= (r_acc==0);
        VectorCL *tmp=0;

        if (new_tmp)
            tmp = new VectorCL(r.size());
        else
            tmp = r_acc;

        loc_norm_sq= LocAccNorm_sq(r,*tmp);

        if (new_tmp)
            delete tmp;
    }
    else if (!acc_r && !useAccur){
        bool newr= (r_acc==0);
        if (newr)
            r_acc = new VectorCL(r);
        else{
            if (r_acc->size()!=r.size()) r_acc->resize(r.size());
            *r_acc=r;
        }
        loc_norm_sq= DotAcc(*r_acc,r);
        if (newr){
            delete r_acc;
            r_acc=0;
        }

    }
    else
        throw DROPSErrCL("ExchangeCL::LocNorm_sq: Cannot perform norm on accumulated vector without flag useAccur\n");

    return loc_norm_sq;
}


/// \brief Norm of a Vector
/** This function computes the norm of a Vector. The Vector \a r can available in accumulated or
    distributed form (according to switch \a acc_r). If wished the accumulated form of \a r will is given
    on return within \a r_acc.*/
double ExchangeCL::Norm(const VectorCL &r, bool acc_r, bool useAccur, VectorCL* r_acc) const
    /// \param[in] r        Vector of which the norm should be computed
    /// \param[in] acc_r    is vector \a r given in accumulated form
    /// \param[in] useAccur should the accurater but slower version be used
    /// \param[out] r_acc   if not equal zero, pointer on the accumulated form of \a r
    /// \return             norm of the vector \a r
{
    double norm_sq=GlobalSum(LocNorm_sq(r, acc_r, useAccur, r_acc));
    if (norm_sq<0.){
        std::cerr << "["<<ProcCL::MyRank()<<"] In function ExchangeCL::Norm:\n Norm of vector smaller than zero: "
                  << "acc_r "<<acc_r<<", use Accur "<<useAccur<<", (bool)r_acc "<<((bool)r_acc)
                  << ", squared value "<<norm_sq<<std::endl;
    }
    DROPS_Check_Norm(norm_sq, "ExchangeCL::Norm: negative squared norm because of accumulation!");

    return std::sqrt(norm_sq);
}

/// \brief Squared norm of a Vector
/** This function computes the square of a norm of a Vector. The Vector \a r can available in accumulated or
    distributed form (according to switch \a acc_r). If wished the accumulated form of \a r will is given
    on return within \a r_acc.*/
double ExchangeCL::Norm_sq(const VectorCL& r, bool acc_r, bool useAccur, VectorCL* r_acc) const
    /// \param[in] r        Vector of which the norm should be computed
    /// \param[in] acc_r    is vector \a r given in accumulated form
    /// \param[in] useAccur should the accurater but slower version be used
    /// \param[out] r_acc   if not equal zero, pointer on the accumulated form of \a r
    /// \return             square norm of the vector \a r
{
    double norm_sq=GlobalSum(LocNorm_sq(r, acc_r, useAccur, r_acc));
    DROPS_Check_Norm(norm_sq, "ExchangeCL::Norm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}

/// \brief Get index of a distributed index on another proc
IdxT ExchangeCL::GetExternalIdxFromProc(IdxT myIdx, ProcNumT proc) const
/** To use this function in CreateList the parameter CreateMapIdx must be set on true!<br>
    This function computes the index of a local index on the processor proc. If the unknown does not
    exist on both procs or the unknown is not distributed NoIdx is returned.
    \pre CreateList must be calld with CreateMap=true
*/
{
    Assert (mapCreated_,
            DROPSErrCL("ExchangeCL::GetExternalIdxFromProc: The mapping: (external idx) -> (my idx) is not created. \nMaybe set flag CreateMapIdx for CreateList()!"),
            DebugParallelNumC);

    MappingIdxCT::const_iterator       proc_it  (MappingIdx_.lower_bound(proc));
    const MappingIdxCT::const_iterator proc_end (MappingIdx_.end());

    if (proc_it==proc_end)
        return NoIdx;

    ProcMappingIdxCT::const_iterator       idx_it  ( (proc_it->second).lower_bound(myIdx) );
    const ProcMappingIdxCT::const_iterator idx_end ( (proc_it->second).end() );

    if (idx_it==idx_end)
        return NoIdx;

    return idx_it->second;
    // return MappingIdx_[proc][myIdx];
}

/// \brief Check if a sysnum is distributed
bool ExchangeCL::IsDist(IdxT i) const
{
    Assert(mapCreated_,
           DROPSErrCL("ExchangeCL::IsDist: The mapping: (external idx) -> (my idx) is not created. \nMaybe set flag CreateMapIdx for CreateList()!"),
           DebugParallelNumC);
    return SysProc_[i].size()>0;
}

/// \brief Check if a sysnum can be found on another proc
bool ExchangeCL::IsOnProc(IdxT i, ProcNumT p)
{
    Assert(mapCreated_,
           DROPSErrCL("ExchangeCL::IsOnProc: The mapping: (external idx) -> (my idx) is not created. \nMaybe set flag CreateMapIdx for CreateList()!"),
           DebugParallelNumC);
    return std::find(SysProc_[i].begin(),SysProc_[i].end(), p) != SysProc_[i].end();

}

/// \brief Get list of procs (except local proc) that owns a sysnum
ExchangeCL::ProcNumCT ExchangeCL::GetProcs(IdxT i) const
{
    Assert(mapCreated_,
           DROPSErrCL("ExchangeCL::GetProcs: The mapping: (external idx) -> (my idx) is not created. \nMaybe set flag CreateMapIdx for CreateList()!"),
           DebugParallelNumC);
    return SysProc_[i];
}

/// \brief Get number of procs, that owns a sysnum
Uint ExchangeCL::GetNumProcs(IdxT i) const
{
    if (ProcCL::Size()==1)
        return 0;
    Assert(mapCreated_,
           DROPSErrCL("ExchangeCL::GetNumProcs: The mapping: (external idx) -> (my idx) is not created. \nMaybe set flag CreateMapIdx for CreateList()!"),
           DebugParallelNumC);
    return SysProc_[i].size()+1;
}

/// \brief Get procs that shares at least one unknown with this proc
ExchangeCL::ProcNumCT ExchangeCL::GetNeighbors() const
{
    ProcNumCT ret;
    if (ProcCL::Size()==1)
        return ret;
    Assert(mapCreated_,
           DROPSErrCL("ExchangeCL::GetNumProcs: The mapping: (external idx) -> (my idx) is not created. \nMaybe set flag CreateMapIdx for CreateList()!"),
           DebugParallelNumC);
    for (CommListCT::const_iterator it(ExList_.begin()), end(ExList_.end()); it!=end; ++it)
        ret.push_back(it->toProc_);
    return ret;
}

/// \brief Parallel InnerProduct of two distributed vectors and accumulate of one vector
/** Parallel InnerProduct with two distributed (unaccumulated) vectors x and y. The result is
    Accumulate(x)^T * y. Watch out: The vector x will be accumulated after the procedure!
 */
double ExchangeCL::ParDotAcc(VectorCL& x, const VectorCL& y)const
    /// \param[in,out] x in local distributed form, after the routine, this vector will be accumulated
    /// \param[in]     y in local distributed form, not changed
    /// \return          Accumulate(x)^T * y
{
    Assert(created_, DROPSErrCL("ExchangeCL::ParDotAcc: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::ParDotAcc: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::ParDotAcc: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

    double loc_sum=0;   // local sum of all elements
    double acc_sum=0;   // local sum of all accumulated elements

    // first send all shared entries (nonblocking)
    InitCommunication(x, SendRecvReq_);

    // while sending, we can sum up unshared entries
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x[LocalIndex[i]] * y[LocalIndex[i]];

    // now we need the accumulated values
    AccFromAllProc(x, SendRecvReq_);

    // these has to be sumed up too
    for (Ulint i=0; i<numDistrIdx_; ++i)
        acc_sum += x[DistrIndex[i]] * y[DistrIndex[i]];

    // now reduce the local_sum's over all procs
    return GlobalSum(loc_sum+acc_sum);
}

/// \brief Parallel InnerProduct of two distributed vectors and accumulate of one vector
/** Parallel InnerProduct with two distributed (unaccumulated) vectors x and y. The result is
    Accumulate(x)^T * y. Watch out: The vector x will be accumulated after the procedure and no global
    reduce operation will be performed!
 */
inline double ExchangeCL::DotAcc(VectorCL& x, const VectorCL& y) const
/// \param[in,out] x in local distributed form, after the routine, this vector will be accumulated
/// \param[in]     y in local distributed form, not changed
/// \return          Accumulate(x)^T * y (without reduce operation!)
{
    Assert(created_, DROPSErrCL("ExchangeCL::DotAcc: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::DotAcc: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::DotAcc: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

    double loc_sum=0;
    RequestCT req(ExList_.size());
    InitCommunication(x, SendRecvReq_);

    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x[LocalIndex[i]] * y[LocalIndex[i]];

    AccFromAllProc(x, SendRecvReq_);

    for (Ulint i=0; i<numDistrIdx_; ++i)
        loc_sum += x[DistrIndex[i]] * y[DistrIndex[i]];

    return loc_sum;
}

/// \brief Parallel InnerProduc of two distributed vectors
/** This function is slower as the function ExchangeCL::ParDotAcc, because a local
    copie of one of the inputvectors have to be made
*/
double ExchangeCL::ParDot(const VectorCL &x, const VectorCL &y) const
    /// \param[in] x in distributed form
    /// \param[in] y in distributed form
    /// \return      InnerProduct of x and y
{
    VectorCL x_acc(x);          // temp vector for accumulation
    return ParDotAcc(x_acc,y);
}


/// \brief Parallel InnerProduct of two distributed vectors, get the accumulated vector
/** This function calculates the InnerProduct of two distributed vectors stores the
    accumulated form of the second parameter in the first parameter
*/
double ExchangeCL::ParDot(VectorCL &x_acc, const VectorCL &x, const VectorCL &y) const
    /// \param[out] x_acc accumulated form of x
    /// \param[in]  x     distributed form.
    /// \param[in]  y     distributed form
    /// \return           InnerProduct of x and y
{
    x_acc=x;
    return ParDotAcc(x_acc,y);
}

/// \brief Parallel InnerProduct of two distributed vectors. Computes both accumulated vectors (maybe accurater)
/** This function calculates the InnerProduct of two distributed vectors. Therefore both vectors will
    accumulate at first and then the InnerProduct will be performed. This function is slower as ParDotAcc(),
    because two accumulations have to be performed. But it should be more accurate.
*/
double ExchangeCL::AccParDot(const VectorCL& x, const VectorCL& y, VectorCL& x_acc, VectorCL& y_acc) const
{
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeCL::AccParDot: Vectors do not have the same size"), DebugParallelNumC);
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::AccParDot: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    if (x_acc.size()!=x.size())
        x_acc.resize(x.size());
    x_acc = x;

    if (y_acc.size()!=y.size())
        y_acc.resize(y.size());
    y_acc = y;

    // Send to all other procs
    RequestCT req_y(2*ExList_.size());
    VectorCL secondRecvBuf(recvBuf_.size());

    InitCommunication(x, SendRecvReq_, tag_, 0);
    InitCommunication(y, req_y, tag_+1, 0, &secondRecvBuf);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x[LocalIndex[i]] * y[LocalIndex[i]];

    // now we need the accumulated values
    AccFromAllProc(x_acc, SendRecvReq_, 0);
    AccFromAllProc(y_acc, req_y, 0, &secondRecvBuf);

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += x_acc[AccDistIndex[i]] * y_acc[AccDistIndex[i]];

    return GlobalSum(loc_sum+acc_sum);
}

/// \brief InnerProduct of two accumulated vectors without global reduce
double ExchangeCL::LocAccDot(const VectorCL& x_acc, const VectorCL& y_acc) const
/// \param[in] x_acc accumulated vector
/// \param[in] y_acc accumulated vector
/// \return          Inner product of x^T y, but without global reduce
{
    Assert(x_acc.size()==y_acc.size(), DROPSErrCL("ExchangeCL::LocAccDot: Vectors do not have the same size"), DebugParallelNumC);
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::LocAccDot: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    double loc_sum=0,
    acc_sum=0;

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x_acc[LocalIndex[i]] * y_acc[LocalIndex[i]];


    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += x_acc[AccDistIndex[i]] * y_acc[AccDistIndex[i]];

    return loc_sum+acc_sum;
}


/// \brief Norm of a distributed unaccumulated vector with higher accuracy
double ExchangeCL::AccNorm_sq(const VectorCL &r, VectorCL& r_acc) const
/// \param[in]  r     unaccumulated vector
/// \param[out] r_acc accumulated vektor r (has not be initialized)
{
    double norm_sq= GlobalSum(LocAccNorm_sq(r,r_acc));
    DROPS_Check_Norm(norm_sq,"ExchangeCL::AccNorm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}

/// \brief Norm of a distributed unaccumulated vector without global reduce (higher accuracy)
double ExchangeCL::LocAccNorm_sq(const VectorCL& r, VectorCL& r_acc) const
/// \param[in]  r     unaccumulated vector
/// \param[out] r_acc accumulated vektor r (has not be initialized)
{
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::LocAccNorm_sq: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    if (r.size()!=r_acc.size())
        r_acc.resize(r.size());
    r_acc = r;

    InitCommunication(r, SendRecvReq_, tag_, 0);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += r[LocalIndex[i]] * r[LocalIndex[i]];

    // now we need the accumulated values
    AccFromAllProc(r_acc, SendRecvReq_, 0);

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += r_acc[AccDistIndex[i]] * r_acc[AccDistIndex[i]];

    return loc_sum+acc_sum;
}

/// \brief Norm of an accumulated vector without global reduce (higher accuracy)
double ExchangeCL::LocAccNorm_sq(const VectorCL& r_acc) const
{
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::LocAccNorm_sq: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    double loc_sum=0,
           acc_sum=0;

    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += r_acc[LocalIndex[i]] * r_acc[LocalIndex[i]];

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += r_acc[AccDistIndex[i]] * r_acc[AccDistIndex[i]];

    return loc_sum+acc_sum;
}

/// \brief InnerProduct with two accumulated vectors
double ExchangeCL::AccParDot(const VectorCL &x, const VectorCL& y_acc, VectorCL& x_acc) const
/// \param[in] x     unaccumulated vector.
/// \param[in] y_acc accumulated vector
/// \param[out] x_acc accumulated form of vector x
/// \return x^T y
{
    Assert(x.size()==y_acc.size(), DROPSErrCL("ExchangeCL::AccParDot: Vectors do not have the same size"), DebugParallelNumC);
    Assert(accIdxCreated_, DROPSErrCL("ExchangeCL::AccParDot: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)"), DebugParallelNumC);

    if (x_acc.size()!=x.size())
        x_acc.resize(x.size());
    x_acc=x;

    InitCommunication(x, SendRecvReq_, tag_, 0);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += x[LocalIndex[i]] * y_acc[LocalIndex[i]];

    // now we need the accumulated values
    AccFromAllProc(x_acc, SendRecvReq_, 0);

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += x_acc[AccDistIndex[i]] * y_acc[AccDistIndex[i]];

    return GlobalSum(loc_sum+acc_sum);
}

/// \brief Norm of an accumulated vector
double ExchangeCL::AccNorm_sq(const VectorCL& r_acc) const
/// \param[in] r_acc accumulated vector
/// \return squared norm of the vector r_acc
{
    Assert(created_, DROPSErrCL("ExchangeCL::AccNrom_sq: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
    Assert(r_acc.size()==vecSize_, DROPSErrCL("ExchangeCL::AccNrom_sq: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (Ulint i=0; i<numLocalIdx_; ++i)
        loc_sum += r_acc[LocalIndex[i]] * r_acc[LocalIndex[i]];

    // now do the global summation
    for (Ulint i=0; i<AccDistIndex.size(); ++i)
        acc_sum += r_acc[AccDistIndex[i]] * r_acc[AccDistIndex[i]];

    double norm_sq = GlobalSum(loc_sum+acc_sum);
    DROPS_Check_Norm(norm_sq, "ExchangeCL::AccNorm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}

/// \brief Accumulation of a Vector
void ExchangeCL::Accumulate(VectorCL &x) const
{
    Assert(created_, DROPSErrCL("ExchangeCL::Accumulate: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
    if (x.size()!=vecSize_)
        std::cerr << "Vector has size "<<x.size()<< ", but should be "<<vecSize_<<std::endl;
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::Accumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);

    RequestCT req(ExList_.size());
    InitCommunication(x, SendRecvReq_);
    AccFromAllProc(x, SendRecvReq_);
}

/// \brief Returns the accumulated form of a vector
VectorCL ExchangeCL::GetAccumulate (const VectorCL &x) const
    /// \param[in] x vector in distributed form
    /// \return      accumulated form of x
{
    Assert(created_, DROPSErrCL("ExchangeCL::GetAccumulate: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
    Assert(x.size()==vecSize_, DROPSErrCL("ExchangeCL::GetAccumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);

    RequestCT req(ExList_.size());
    InitCommunication(x, SendRecvReq_);
    VectorCL x_acc(x);
    AccFromAllProc(x_acc, SendRecvReq_);
    return x_acc;
}

/// \brief Return accumulated vectors
std::vector<VectorCL> ExchangeCL::GetAccumulate (const std::vector<VectorCL>& x) const
{
    Assert(created_, DROPSErrCL("ExchangeCL::GetAccumulate: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<x.size(); ++i)
        Assert(x[i].size()==vecSize_, DROPSErrCL("ExchangeCL::GetAccumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);
#endif

    // allocate memory for requests and receive buffers
    std::valarray<RequestCT> req(x.size());
    std::valarray<VectorCL>  recvBufs(x.size());
    for (size_t i=0; i<x.size(); ++i){
        req[i].resize(ExList_.size());
        recvBufs[i].resize(recvBuf_.size());
    }

    // send and receive entries of x
    for (size_t i=0; i<x.size(); ++i)
        InitCommunication(x[i], req[i], tag_+i, 0, &recvBufs[i]);

    // allocate mem for x_acc and init with x
    std::vector<VectorCL> x_acc(x.size());
    for (size_t i=0; i<x.size(); ++i){
        x_acc[i].resize(x[i].size());
        x_acc[i]=x[i];
    }

    // do accumulation
    for (size_t i=0; i<x.size(); ++i)
        AccFromAllProc(x_acc[i], req[i], 0, &recvBufs[i]);

    return x_acc;
}

/// \brief Calculate the square of the euclidian-norm and accumulates the vector
double ExchangeCL::Norm_sq_Acc(VectorCL &r_acc, const VectorCL &r) const
    /// \param[out] r_acc accumulated form of r (can be uninitialized)
    /// \param[in]  r local distributed vecot
    /// \return       r^T * r
{
    r_acc=r;

    double norm_sq = ParDotAcc(r_acc,r);
    DROPS_Check_Norm(norm_sq,"ExchangeCL::Norm_sq_Acc: negative squared norm because of accumulation!");
    return norm_sq;
}

/// \brief Calculate the square of the euclidian-norm of a vector
double ExchangeCL::Norm_sq(const VectorCL &r) const
    /// \param[in] r distributed form of a vector
    /// \return      squared euclidian norm of the vector r
{
    VectorCL r_acc(r);

    double norm_sq = ParDotAcc(r_acc,r);
    DROPS_Check_Norm(norm_sq,"ExchangeCL::Norm_sq_Acc: negative squared norm because of accumulation!");
    return norm_sq;
}

/// \brief Returns the euclidian-norm of a vector
double ExchangeCL::Norm(const VectorCL &r) const
    /// \param[in] r distributed form of a vector
    /// \return      euclidian norm of the vector r
{
    VectorCL r_acc(r);

    double norm_sq = ParDotAcc(r_acc,r);

    if (norm_sq<0.){
        std::cerr << "["<<ProcCL::MyRank()<<"] In function ExchangeCL::Norm (1 arg):\n Norm of vector smaller than zero: "
                  << "squared value "<<norm_sq<<std::endl;
        throw DROPSErrCL("ExchangeCL::Norm: negative squared norm because of accumulation!");
    }
    DROPS_Check_Norm(norm_sq,"ExchangeCL::Norm: negative squared norm because of accumulation!");

    return std::sqrt(norm_sq);
}

// -------------------------------------
// E X C H A N G E  B L O C K  C L A S S
// -------------------------------------

// Information functions
//----------------------

/// \brief returns the ExchangeCL corresponding to the i-th block
const ExchangeCL& ExchangeBlockCL::Get(size_t i) const{
    Assert(i<m_, DROPSErrCL("ExchangeBlockCL::Get: block does not exist"), DebugNumericC);
    return ExchgList_[i];
}

/// \brief returns the ExchangeCL corresponding to the i-th block
ExchangeCL& ExchangeBlockCL::Get(size_t i){
    Assert(i<m_, DROPSErrCL("ExchangeBlockCL::Get: block does not exist"), DebugNumericC);
    return ExchgList_[i];
}

/// \brief returns the number of blocks handled by this class
size_t ExchangeBlockCL::GetNumBlocks() const{
    return m_;
}

/// \brief check if all exchange-lists have been created
bool ExchangeBlockCL::Created()
{
    if (created_) return true;
    bool check=true;
    for (size_t i=0; i<m_ && check; ++i)
        check = check && ExchgList_[i].Created();
    if (check){
        created_=true;
        CreateBlockIdx();
    }
    return check;
}

/// \brief get the size of vector, that can be accumulated
Ulint ExchangeBlockCL::GetNum() const
{
    Assert(created_ && blockCreated_, DROPSErrCL("ExchangeBlockCL::GetNum: Lists have not bee created yet"), DebugNumericC);
    return Block_[m_];
}

/// \brief get the size of i-th block
Ulint ExchangeBlockCL::GetNum(size_t i) const
{
    Assert(ExchgList_[i].Created(), DROPSErrCL("ExchangeBlockCL::GetNum(size_t i): i-th ExchangeCL not created yet"), DebugNumericC);
    return ExchgList_[i].GetNum();
}

// Communication functions
//------------------------

void ExchangeBlockCL::InitCommunication(const VectorCL &vec, VecRequestCT& req, int tag, std::vector<VectorCL>* recvBuf) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::InitCommunication: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);

    const int mytag= tag<0 ? start_tag_ : tag;
    for (size_t i=0; i<m_; ++i)
        ExchgList_[i].InitCommunication(vec, req[i], mytag+i, Block_[i], ((recvBuf) ? &((*recvBuf)[i]) : 0));
}

void ExchangeBlockCL::AccFromAllProc(VectorCL &vec, VecRequestCT& req, std::vector<VectorCL>* recvBuf) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::AccFromAllProc: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    for (size_t i=0; i<m_; ++i)
        ExchgList_[i].AccFromAllProc(vec, req[i], Block_[i], ((recvBuf) ? &((*recvBuf)[i]) : 0));
}

double ExchangeBlockCL::ParDotAcc(VectorCL &x, const VectorCL &y) const
{
    return GlobalSum(DotAcc(x,y));
}

double ExchangeBlockCL::DotAcc(VectorCL& x, const VectorCL& y) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::DotAcc: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeBlockCL::DotAcc: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::DotAcc: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

    // local sum of all elements
    double loc_sum=0;
    // Send shared entries
    InitCommunication(x,SendRecvReq_);
    // while sending, compute local parts
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x[ExchgList_[i].LocalIndex[j]+Block_[i]] * y[ExchgList_[i].LocalIndex[j]+Block_[i]]; // use offset
    // recieve shared parts
    AccFromAllProc(x,SendRecvReq_);
    // compute part of shared verts
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumDistIdx(); ++j)                                               // for all local indices of this block
            loc_sum += x[ExchgList_[i].DistrIndex[j]+Block_[i]] * y[ExchgList_[i].DistrIndex[j]+Block_[i]]; // use offset
    // return local sum
    return loc_sum;
}

// InnerProduct: Both vectors will be accumulated, this is more accurate
double ExchangeBlockCL::AccParDot(const VectorCL& x, const VectorCL& y, VectorCL& x_acc, VectorCL& y_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::AccParDot: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeBlockCL::AccParDot: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::AccParDot: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeCL::AccParDot: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    if (x_acc.size()!=x.size())
        x_acc.resize(x.size());
    x_acc = x;

    if (y_acc.size()!=y.size())
        y_acc.resize(y.size());
    y_acc = y;

    // Allocate memory for requests and a recieve buffer
    VecRequestCT req_y(GetNumBlocks());
    std::vector<VectorCL> secondRecvBuf(GetNumBlocks());
    for (size_t i=0; i<GetNumBlocks(); ++i){
        req_y[i].resize(2*ExchgList_[i].GetNumNeighs());
        secondRecvBuf[i].resize(ExchgList_[i].GetNumReceiveElements());
    }

    // Send to all other procs
    InitCommunication(x, SendRecvReq_, start_tag_);
    InitCommunication(y, req_y, start_tag_+m_, &secondRecvBuf);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x[ExchgList_[i].LocalIndex[j]+Block_[i]] * y[ExchgList_[i].LocalIndex[j]+Block_[i]];

    // now we need the accumulated values
    AccFromAllProc(x_acc, SendRecvReq_);
    AccFromAllProc(y_acc, req_y, &secondRecvBuf);

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += x_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]] * y_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    return GlobalSum(loc_sum+acc_sum);
}


// InnerProduct: with two accumulated vectors. The first given vec should be accumulated
double ExchangeBlockCL::AccParDot(const VectorCL& x, const VectorCL& y_acc, VectorCL& x_acc) const
{
    Assert(x.size()==y_acc.size(), DROPSErrCL("ExchangeCL::AccParDot: Vectors do not have the same size"), DebugParallelNumC);
    Assert(created_, DROPSErrCL("ExchangeBlockCL::DotAcc: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::DotAcc: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeCL::AccParDot: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    if (x_acc.size()!=x.size())
        x_acc.resize(x.size());
    x_acc=x;

    InitCommunication(x, SendRecvReq_);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x[ExchgList_[i].LocalIndex[j]+Block_[i]] * y_acc[ExchgList_[i].LocalIndex[j]+Block_[i]];

    // now we need the accumulated values
    AccFromAllProc(x_acc, SendRecvReq_);

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += x_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]] * y_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    return GlobalSum(loc_sum+acc_sum);
}

double ExchangeBlockCL::LocAccDot(const VectorCL& x_acc, const VectorCL& y_acc) const
{
    Assert(x_acc.size()==y_acc.size(), DROPSErrCL("ExchangeCL::LocAccDot: Vectors do not have the same size"), DebugParallelNumC);
    Assert(created_, DROPSErrCL("ExchangeBlockCL::LocAccDot: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x_acc.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::LocAccDot: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeCL::LocAccDot: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x_acc[ExchgList_[i].LocalIndex[j]+Block_[i]] * y_acc[ExchgList_[i].LocalIndex[j]+Block_[i]];

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += x_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]] * y_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    return loc_sum+acc_sum;
}

double ExchangeBlockCL::AccNorm_sq(const VectorCL& r, VectorCL& r_acc) const
{
    double norm_sq = GlobalSum(LocAccNorm_sq(r,r_acc));
    DROPS_Check_Norm(norm_sq,"ExchangeCL::Norm: negative squared norm because of accumulation!");
    return norm_sq;
}

double ExchangeBlockCL::LocAccNorm_sq(const VectorCL& r, VectorCL& r_acc) const
{
#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeCL::locAccNorm_sq: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif
    if (r.size()!=r_acc.size())
        r_acc.resize(r.size());
    r_acc = r;

    InitCommunication(r, SendRecvReq_);

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += r[ExchgList_[i].LocalIndex[j]+Block_[i]] * r[ExchgList_[i].LocalIndex[j]+Block_[i]];

    // now we need the accumulated values
    AccFromAllProc(r_acc, SendRecvReq_);

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += r_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]] * r_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    return loc_sum+acc_sum;
}

/// \brief Squared norm of an accumulated vector without a global reduce
double ExchangeBlockCL::LocAccNorm_sq(const VectorCL& r_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::LocAccNorm_sq: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(r_acc.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::LocAccNorm_sq: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);
#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeCL::LocAccNorm_sq: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    double loc_sum=0,
           acc_sum=0;

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += r_acc[ExchgList_[i].LocalIndex[j]+Block_[i]] * r_acc[ExchgList_[i].LocalIndex[j]+Block_[i]];


    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += r_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]] * r_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    return loc_sum+acc_sum;
}

/// \brief Squared norm of an accumulated vector
double ExchangeBlockCL::AccNorm_sq(const VectorCL& r_acc) const
{
    double norm_sq = GlobalSum(LocAccNorm_sq(r_acc));
    DROPS_Check_Norm(norm_sq,"ExchangeBlockCL::AccNorm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}


void ExchangeBlockCL::Accumulate(VectorCL &x) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::Accumulate: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::Accumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);

    InitCommunication(x, SendRecvReq_);
    AccFromAllProc(x, SendRecvReq_);
}

VectorCL ExchangeBlockCL::GetAccumulate (const VectorCL &x) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::GetAccumulate: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::GetAccumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);

    InitCommunication(x, SendRecvReq_);
    VectorCL x_acc(x);
    AccFromAllProc(x_acc, SendRecvReq_);
    return x_acc;
}

std::vector<VectorCL> ExchangeBlockCL::GetAccumulate (const std::vector<VectorCL>& x) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::GetAccumulate: Lists have not been created (Maybe use CreateList before!)\n"), DebugParallelNumC);
#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<x.size(); ++i)
        Assert(x[i].size()==GetNum(), DROPSErrCL("ExchangeBlockCL::GetAccumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);
#endif
    // Allocate memory for requests and receive buffers
    std::valarray<VecRequestCT> req(x.size());
    std::vector<std::vector<VectorCL> > RecvBuffers(x.size());
    for (size_t i=0; i<x.size(); ++i){
        req[i].resize(GetNumBlocks());
        RecvBuffers[i].resize(GetNumBlocks());
    }
    for (size_t i=0; i<x.size(); ++i){
         for (size_t j=0; j<GetNumBlocks(); ++j){
            req[i][j].resize(2*ExchgList_[i].GetNumNeighs());
            RecvBuffers[i][j].resize(ExchgList_[i].GetNumReceiveElements());
         }
    }

    // send and receive entires of x
    for (size_t i=0; i<x.size(); ++i)
        InitCommunication(x[i], req[i], start_tag_+i*m_, &(RecvBuffers[i]));

    // allocate mem for x_acc and init with x
    std::vector<VectorCL> x_acc(x.size());
    for (size_t i=0; i<x.size(); ++i){
        x_acc[i].resize(x[i].size());
        x_acc[i]=x[i];
    }

    // do accumulation
    for (size_t i=0; i<x.size(); ++i)
        AccFromAllProc(x_acc[i], req[i], &(RecvBuffers[i]));

    return x_acc;
}


double ExchangeBlockCL::Norm_sq_Acc(VectorCL &r_acc, const VectorCL &r) const
{
    if (r_acc.size()!=r.size())
        r_acc.resize(r.size());
    r_acc=r;
    double norm_sq = ParDotAcc(r_acc,r);
    DROPS_Check_Norm(norm_sq,"ExchangeBlockCL::AccNorm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}

double ExchangeBlockCL::Norm(const VectorCL &r) const
{
    VectorCL r_acc(r);
    double norm_sq = ParDotAcc(r_acc,r);
    DROPS_Check_Norm(norm_sq,"ExchangeBlockCL::Norm old: negative squared norm because of accumulation!");
    return std::sqrt(norm_sq);
}

double ExchangeBlockCL::Norm_sq(const VectorCL &r) const
{
    VectorCL r_acc(r);
    double norm_sq = ParDotAcc(r_acc,r);
    DROPS_Check_Norm(norm_sq,"ExchangeBlockCL::Norm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}

double ExchangeBlockCL::ParDot(const VectorCL &x, const VectorCL &y) const{
    VectorCL x_acc(x);          // temp vector for accumulation
    return ParDotAcc(x_acc,y);
}

double ExchangeBlockCL::ParDot(VectorCL &x_acc, const VectorCL &x, const VectorCL &y) const{
    if (x_acc.size()!=x.size())
        x_acc.resize(x.size());
    x_acc=x;
    return ParDotAcc(x_acc,y);
}

/****************************************************
*   L O C  D O T _                                  *
*****************************************************
*  inner product of two distributed vectors x and   *
*  y without performing an global reduce. The       *
*  vector x will be accumulated and stored in x_acc *
****************************************************/
double ExchangeBlockCL::LocDot_(const VectorCL& x, const VectorCL& y, VectorCL* x_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::LocDot_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::Accumulate: vector length does not fit to the created lists. (Maybe used a wrong IdxDescCL?)"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeBlockCL::LocDot_: Vectors do not have the same length"), DebugParallelNumC);

    bool newx = (x_acc==0);

    double loc_sum=0;                       // sum of local entries
    double dist_sum=0;                      // sum of distributed entries

    InitCommunication(x, SendRecvReq_);        // send shared entries of x to all neighbor procs

    // assign all values of x to x_acc
    if (newx)
        x_acc= new VectorCL(x);
    else{
        Assert(x_acc->size()==GetNum(), DROPSErrCL("ExchangeBlockCL::LocDot_: vector x_acc has not the right length"),DebugParallelNumC);
        *x_acc=x;
    }

    // local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x[ExchgList_[i].LocalIndex[j]+Block_[i]] * y[ExchgList_[i].LocalIndex[j]+Block_[i]]; // use offset

    AccFromAllProc(*x_acc, SendRecvReq_);   // recieve values from neighbors and sum them up

    // summation over shared indices
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumDistIdx(); ++j)                                               // for all local indices of this block
            dist_sum += (*x_acc)[ExchgList_[i].DistrIndex[j]+Block_[i]] * y[ExchgList_[i].DistrIndex[j]+Block_[i]]; // use offset

    if (newx){                              // if new x_acc where created give memory free
        delete x_acc;
        x_acc=0;
    }

    return loc_sum+dist_sum;                // return result
}

/****************************************************
*   A C C U R  L O C  D O T  N O  A C C _           *
*****************************************************
*  inner product of two accumulated vectors x and   *
*  y without performing an global reduce.           *
*****************************************************
*  pre: both vectors are accumulated                *
****************************************************/
double ExchangeBlockCL::AccurLocDotNoAcc_(const VectorCL& x, const VectorCL& y) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::AccurLocDotNoAcc_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeBlockCL::AccurLocDotNoAcc_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::AccurLocDotNoAcc_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeBlockCL::AccurLocDotNoAcc_: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    double loc_sum=0,       // sum of local entries
           acc_sum=0;       // sum of distributed entries

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x[ExchgList_[i].LocalIndex[j]+Block_[i]] * y[ExchgList_[i].LocalIndex[j]+Block_[i]];

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += x[ExchgList_[i].AccDistIndex[j]+Block_[i]] * y[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    return loc_sum+acc_sum;
}

/****************************************************
*   A C C U R  L O C  D O T  O N E  A C C _         *
*****************************************************
*  Accurate inner product of one accumulated vector *
*  (x_acc) and one distributed vecor (y) without    *
*  performing an global reduce. If y_acc not equal  *
*  zero, the acumulated form of y will be given     *
*  by return in y_acc.                              *
*****************************************************
*  pre: x_acc is accumulated, y is distributed      *
****************************************************/
double ExchangeBlockCL::AccurLocDotOneAcc_(const VectorCL& x_acc, const VectorCL& y, VectorCL* y_acc) const
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::AccurLocDotOneAcc_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x_acc.size()==y.size(), DROPSErrCL("ExchangeBlockCL::AccurLocDotOneAcc_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x_acc.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::AccurLocDotOneAcc_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeBlockCL::AccurLocDotOneAcc_: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    bool newy = (y_acc==0);

    double loc_sum=0,                       // sum of local entries
           acc_sum=0;                       // sum of distributed entries

    InitCommunication(y, SendRecvReq_);        // send shared entries of x to all neighbor procs

    if (newy)
        y_acc = new VectorCL(y);
    else{
        Assert(y_acc->size()==x_acc.size(), DROPSErrCL("ExchangeBlockCL::AccurLocDotOneAcc_: y_acc has not the right size"), DebugParallelNumC);
        *y_acc=y;
    }

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x_acc[ExchgList_[i].LocalIndex[j]+Block_[i]] * y[ExchgList_[i].LocalIndex[j]+Block_[i]];

    AccFromAllProc(*y_acc, SendRecvReq_);   // recieve values from neighbors and sum them up

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += x_acc[ExchgList_[i].AccDistIndex[j]+Block_[i]] * (*y_acc)[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    // free memory
    if (newy){
        delete y_acc;
        y_acc=0;
    }

    return loc_sum+acc_sum;
}

/****************************************************
*   A C C U R  L O C  D O T  B O T H   A C C _      *
*****************************************************
*  Accurate inner product of two distributed        *
*  vectors (x and y) without performing an global   *
*  reduce. If x_acc or y_acc not equal zero, the    *
*  acumulated form of x or y will be given by return*
*  in x_acc or y_acc.                               *
*****************************************************
*  pre: x and y have distributed form               *
****************************************************/
double ExchangeBlockCL::AccurLocDotBothAcc_(const VectorCL& x, const VectorCL& y, VectorCL* x_acc, VectorCL* y_acc) const
{
    //throw DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: This function does not work ... ");

    Assert(created_, DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

#if DROPSDebugC&DebugParallelNumC
    for (size_t i=0; i<m_; ++i)
        if (!ExchgList_[i].AccIdxCreated())
            throw DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: Indices for accumulated distributed indices has not been created. (Maybe use CreateAccDist for CreateList)");
#endif

    bool newx = (x_acc==0);
    bool newy = (y_acc==0);

    double loc_sum=0,                       // sum of local entries
           acc_sum=0;                       // sum of distributed entries

    // Allocate mem for requests for second vector
    VecRequestCT req_y(GetNumBlocks());
    std::vector<VectorCL> secondRecvBuf(GetNumBlocks());
    for (size_t i=0; i<GetNumBlocks(); ++i){
        req_y[i].resize(2*ExchgList_[i].GetNumNeighs());
        secondRecvBuf[i].resize(ExchgList_[i].GetNumReceiveElements());
    }

    InitCommunication(x, SendRecvReq_, start_tag_);             // send and receive shared entries of x to all neighbor procs
    InitCommunication(y, req_y, start_tag_+m_, &secondRecvBuf);   // send and receive shared entries of x to all neighbor procs

    if (newx)
        x_acc = new VectorCL(x);
    else{
        Assert(x_acc->size()==x.size(), DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: x_acc has not the right size"), DebugParallelNumC);
        *x_acc=x;
    }

    if (newy)
        y_acc = new VectorCL(y);
    else{
        Assert(y_acc->size()==x.size(), DROPSErrCL("ExchangeBlockCL::AccurLocDotBothAcc_: y_acc has not the right size"), DebugParallelNumC);
        *y_acc=y;
    }

    // do local summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].GetNumLocIdx(); ++j)                                                // for all local indices of this block
            loc_sum += x[ExchgList_[i].LocalIndex[j]+Block_[i]] * y[ExchgList_[i].LocalIndex[j]+Block_[i]];

    AccFromAllProc(*x_acc, SendRecvReq_);               // sum elemnts up
    AccFromAllProc(*y_acc, req_y, &secondRecvBuf);      // sum elemnts up

    // now do the global summation
    for (size_t i=0; i<m_; ++i)                                                                             // in every block
        for (Ulint j=0; j<ExchgList_[i].AccDistIndex.size(); ++j)                                           // for all local indices of this block
            acc_sum += (*x_acc)[ExchgList_[i].AccDistIndex[j]+Block_[i]] * (*y_acc)[ExchgList_[i].AccDistIndex[j]+Block_[i]];

    // free memory
    if (newx){
        delete x_acc;
        x_acc=0;
    }

    if (newy){
        delete y_acc;
        y_acc=0;
    }

    return loc_sum+acc_sum;
}

/// \brief Local inner product without global reduce
/** This function computes the inner product of two vectors without global reduce. The vectors may be available
    in different forms. Therefore the \a acc_x and \a acc_y flag exists. If an accumulation is performed, the
    result can be given to the caller.*/
double ExchangeBlockCL::LocDot (const VectorCL& x, bool acc_x,
                                const VectorCL& y, bool acc_y,
                                bool useAccur,
                                VectorCL* x_acc, VectorCL *y_acc) const
    /// \param[in]  x        first vector
    /// \param[in]  acc_x    is vector \a x given in accumulated form
    /// \param[in]  y        second vector
    /// \param[in]  acc_y    is vector \a y given in accumulated form
    /// \param[in]  useAccur should the accurater but slower version be used
    /// \param[out] x_acc    if not equal zero, pointer on the accumulated form of \a x (also non distributed values will be copied)
    /// \param[out] y_acc    if not equal zero, pointer on the accumulated form of \a y (also non distributed values will be copied)
    /// \return              inner product of \a x and \a y without global reduce
{
    Assert(created_, DROPSErrCL("ExchangeBlockCL::LocDot: Lists have not been created (Maybe use CreateList before!\n"), DebugParallelNumC);
    Assert(x.size()==y.size(), DROPSErrCL("ExchangeBlockCL::LocDot: Vectors do not have the same length"), DebugParallelNumC);
    Assert(x.size()==GetNum(), DROPSErrCL("ExchangeBlockCL::LocDot: vector length does not fit to the created lists. Maybe used a wrong IdxDescCL?"), DebugParallelNumC);

    if (acc_x && x_acc!=0)
        *x_acc=x;
    if (acc_y && y_acc!=0)
        *y_acc=y;

    if (useAccur){
        if (acc_x && acc_y)
            return AccurLocDotNoAcc_(x,y);
        if (acc_x && !acc_y)
            return AccurLocDotOneAcc_(x,y,y_acc);
        if (acc_y && !acc_x)
            return AccurLocDotOneAcc_(y,x,x_acc);
        if (!acc_x && !acc_y)
            return AccurLocDotBothAcc_(x,y,x_acc,y_acc);
    }
    else{
        if (!acc_x && !acc_y)
        {   // x and y are not accumulated, so accumulate one
            if (x_acc && y_acc==0)              // accumulate x
                return LocDot_(x,y,x_acc);
            if (y_acc && x_acc==0)              // accumulate y
                return LocDot_(y,x,y_acc);
            if (x_acc==0 && y_acc==0)           // no accumulated form is wished
                return LocDot_(x,y,0);
            // x and y should be accumulated. That makes no sence
            throw DROPSErrCL("ExchangeBlockCL::LocDot: It makes no sence to do LocDot to accumulate both vectors. Better use accurate version");
        }
        // form here on at least one vector is accumulated
        if (acc_x && !acc_y)
            return dot(x,y);
        if (acc_y && !acc_x)
            return dot(x,y);
        // if both vectors are accumulated, you are not allowed to call this function with useAccur==false
        if (acc_x && acc_y)
            throw DROPSErrCL("ExchangeBlockCL::LocDot: Cannot perform a normal inner product on two accumulated vectors, set useAccur=true");
    }
    throw DROPSErrCL("ExchangeBlockCL::LocDot: Internal error, no matching found");
}

/// \brief Local inner product with global reduce
double ExchangeBlockCL::ParDot (const VectorCL& x, bool acc_x,
                                const VectorCL& y, bool acc_y,
                                bool useAccur,
                                VectorCL* x_acc, VectorCL *y_acc) const
/** This function computes the inner product of two vectors. The vectors may be available in different forms.
    Therefore the \a acc_x and \a acc_y flag exists. If an accumulation is performed, the
    result can be given to the caller.*/
    /// \param[in]  x        first vector
    /// \param[in]  acc_x    is vector \a x given in accumulated form
    /// \param[in]  y        second vector
    /// \param[in]  acc_y    is vector \a y given in accumulated form
    /// \param[in]  useAccur should the accurater but slower version be used
    /// \param[out] x_acc    if not equal zero, pointer on the accumulated form of \a x
    /// \param[out] y_acc    if not equal zero, pointer on the accumulated form of \a y
    /// \return              inner product of \a x and \a y without global reduce
{
    return GlobalSum(LocDot(x, acc_x, y, acc_y, useAccur, x_acc, y_acc));
}

/// \brief Squared norm of a vector without global reduce
/** This function computes the squared norm of a Vector without global reduce. The Vector \a r can available in
    accumulated or distributed form (according to switch \a acc_r). If wished the accumulated form of \a r will
    is given on return within \a r_acc.*/
double ExchangeBlockCL::LocNorm_sq(const VectorCL &r, bool acc_r, bool useAccur, VectorCL* r_acc) const
    /// \param[in] r        Vector of which the norm should be computed
    /// \param[in] acc_r    is vector \a r given in accumulated form
    /// \param[in] useAccur should the accurater but slower version be used
    /// \param[out] r_acc   if not equal zero, pointer on the accumulated form of \a r
    /// \return             norm of the vector \a r
{
    double loc_norm_sq;//=LocDot(r, acc_r, r, acc_r, useAccur, r_acc, 0);
    if (acc_r && useAccur){
        loc_norm_sq= LocAccNorm_sq(r);
    }
    else if (!acc_r && useAccur){
        if (r_acc)
            loc_norm_sq= LocAccNorm_sq(r,*r_acc);
        else{
            VectorCL r_tmp(r.size());
            loc_norm_sq= LocAccNorm_sq(r,r_tmp);
        }
    }
    else
        loc_norm_sq=LocDot(r, acc_r, r, acc_r, useAccur, r_acc, 0);
    if (loc_norm_sq<0.){
        std::cerr << "EXBLOCKCL: local norm is "<<loc_norm_sq<<" computed by "<<(useAccur?"accur":"not accur")<<" variant"<<std::endl;
    }

    return loc_norm_sq;
}

/// \brief Norm of a Vector
/** This function computes the norm of a Vector. The Vector \a r can available in accumulated or
    distributed form (according to switch \a acc_r). If wished the accumulated form of \a r will is given
    on return within \a r_acc.*/
double ExchangeBlockCL::Norm(const VectorCL &r, bool acc_r, bool useAccur, VectorCL* r_acc) const
    /// \param[in] r        Vector of which the norm should be computed
    /// \param[in] acc_r    is vector \a r given in accumulated form
    /// \param[in] useAccur should the accurater but slower version be used
    /// \param[out] r_acc   if not equal zero, pointer on the accumulated form of \a r
    /// \return             norm of the vector \a r
{
    double norm_sq= GlobalSum(LocNorm_sq(r, acc_r, useAccur, r_acc));
    DROPS_Check_Norm(norm_sq,"ExchangeBlockCL::Norm new: negative squared norm because of accumulation!");
    return std::sqrt(norm_sq);
}

/// \brief Squared norm of a Vector
/** This function computes the square of a norm of a Vector. The Vector \a r can available in accumulated or
    distributed form (according to switch \a acc_r). If wished the accumulated form of \a r will is given
    on return within \a r_acc.*/
double ExchangeBlockCL::Norm_sq(const VectorCL& r, bool acc_r, bool useAccur, VectorCL* r_acc) const
    /// \param[in] r        Vector of which the norm should be computed
    /// \param[in] acc_r    is vector \a r given in accumulated form
    /// \param[in] useAccur should the accurater but slower version be used
    /// \param[out] r_acc   if not equal zero, pointer on the accumulated form of \a r
    /// \return             square norm of the vector \a r
{
    double norm_sq = GlobalSum(LocNorm_sq(r, acc_r, useAccur, r_acc));
    DROPS_Check_Norm(norm_sq,"ExchangeBlockCL::Norm_sq: negative squared norm because of accumulation!");
    return norm_sq;
}
} // end of namespace DROPS
