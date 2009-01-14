// **************************************************************************
// File:    parfastmarch.tpp                                                *
// Content: Classes for performing a parallel fastmarching algorithm        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen     *
//          Oliver Fortmeier, RZ RWTH Aachen                                *
// Version: 0.1                                                             *
// Date:                                                                    *
// Begin:   August  21th, 2006                                              *
// **************************************************************************
/// \author Oliver Fortmeier
/// \file parfastmarch.tpp

#include <map>
namespace DROPS{

template<typename ConT>
  bool IsIn(const ConT &T, const typename ConT::value_type val)
{
    return std::find(T.begin(), T.end(), val)!=T.end();
}

//------------------------
// R E P R  T E T R A  C L
//------------------------

IdxT ReprTetraCL::operator[] (Uint i) const
/** Get DoF
    \param[in] i number of the vert within the tetra
*/
{
    Assert(i<4, DROPSErrCL("ReprTetraCL::[]: Only four vertices of tetra exists"),DebugParallelNumC);
    return data_[i].second;
}

bool ReprTetraCL::IsGhost(Uint i) const
/** Check if a vertex is ghost
    \param[in] i number of the vert within the tetra
*/
{
    Assert(i<4, DROPSErrCL("ReprTetraCL::IsDist: Only four vertices of tetra exists"),DebugParallelNumC);
    return data_[i].first;
}

void ReprTetraCL::Set(Uint i, const FMIdxT& idx)
/** Set the DoF and a ghost flag
    \param[in] i number of the vert within the tetra
    \param[in] idx the DoF and ghost-flag of the vert i
*/
{
    Assert(i<4, DROPSErrCL("ReprTetraCL::Set: Only four vertices of tetra exists"),DebugParallelNumC);
    data_[i] = idx;
}

FMIdxT ReprTetraCL::Get(Uint i) const
/** Get DoF and ghost-flag
    \param[in] i number of the vert within the tetra
*/
{
    Assert(i<4, DROPSErrCL("ReprTetraCL::Get: Only four vertices of tetra exists"),DebugParallelNumC);
    return data_[i];
}


//--------------------------
// F M  T R A N S F E R  C L
//--------------------------

template<typename ExCL>
  FMTransferCL<ExCL>::FMTransferCL(VecDescCL& v, VectorBaseCL<Point3DCL>& coord,                    // local DoF
                                   VectorBaseCL<VertexNeighT>& neigh, VectorBaseCL<byte>& typ,      // neighbor tetras of DoF
                                   std::vector<double>& v_ext, std::vector<Point3DCL>& coord_ext,   // external DoF
                                   std::vector<byte>& typ_ext,                                      // types on external DoF
                                   VectorBaseCL<GhostsList>& hasghost,                              // hasGhost DoF
                                   ExCL& ex)                                                        // ExchangeCL
    : FlagIdxMPIT_(ProcCL::NullDataType), CoordValMPIT_(ProcCL::NullDataType),
      IdxValMPIT_(ProcCL::NullDataType), IdxTypeMPIT_(ProcCL::NullDataType),
      StopIdx(NoIdx-1), v_(v), Coord_(coord), neigh_(neigh), Typ_(typ), DataExt_(v_ext),
      CoordExt_(coord_ext), TypExt_(typ_ext), HasGhost_(hasghost), ex_(ex), neighs_(ex.GetNeighbors()),
      tag_(3001)
{
    CreateFlagIdxMPIT();
    CreateCoordValMPIT();
    CreateIdxValMPIT();
    CreateIdxTypeMPIT();
}

template<typename ExCL>
  FMTransferCL<ExCL>::~FMTransferCL()
{
    if (FlagIdxMPIT_==ProcCL::NullDataType)
        ProcCL::Free(FlagIdxMPIT_);
    if (CoordValMPIT_==ProcCL::NullDataType)
        ProcCL::Free(CoordValMPIT_);
    if (IdxValMPIT_==ProcCL::NullDataType)
        ProcCL::Free(IdxValMPIT_);
}

// Creating MPI::Datatypes
//------------------------

template<typename ExCL>
  void FMTransferCL<ExCL>::CreateFlagIdxMPIT()
{
    if (FlagIdxMPIT_==ProcCL::NullDataType)
        return;

    CoupFlagIdxS tmp;

    const int count               = 2;
    const int block[2]            = {1,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp);
    const ProcCL::DatatypeT typ[2]= { ProcCL::MPI_TT<char>::dtype, ProcCL::MPI_TT<IdxT>::dtype };

    ProcCL::AintT displace[2];
    displace[0]= ProcCL::Get_address(&tmp.flag)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.dof) -sAddr;

    FlagIdxMPIT_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(FlagIdxMPIT_);
}

template<typename ExCL>
  void FMTransferCL<ExCL>::CreateCoordValMPIT()
{
    CoupCoordValS tmp;

    const int count               = 3;
    const int block[3]            = {3,1,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp);
    const ProcCL::DatatypeT typ[3]= { ProcCL::MPI_TT<double>::dtype, ProcCL::MPI_TT<double>::dtype,
                                      ProcCL::MPI_TT<char>::dtype };

    ProcCL::AintT displace[3];
    displace[0]= ProcCL::Get_address(&tmp.coord)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.val)  -sAddr;
    displace[2]= ProcCL::Get_address(&tmp.typ)  -sAddr;

    CoordValMPIT_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(CoordValMPIT_);
}

template<typename ExCL>
  void FMTransferCL<ExCL>::CreateIdxValMPIT()
{
    CoupIdxValS tmp;

    const int count               = 3;
    const int block[3]            = {1,1,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp);
    const ProcCL::DatatypeT typ[3]= { ProcCL::MPI_TT<char>::dtype, ProcCL::MPI_TT<IdxT>::dtype,
                                      ProcCL::MPI_TT<double>::dtype };

    ProcCL::AintT displace[3];
    displace[0]= ProcCL::Get_address(&tmp.flag)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.dof) -sAddr;
    displace[2]= ProcCL::Get_address(&tmp.val) -sAddr;

    IdxValMPIT_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(IdxValMPIT_);
}

template<typename ExCL>
  void FMTransferCL<ExCL>::CreateIdxTypeMPIT()
{
    CoupIdxTypeS tmp;

    const int count               = 3;
    const int block[3]            = {1,1,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp);
    const ProcCL::DatatypeT typ[3]= { ProcCL::MPI_TT<char>::dtype, ProcCL::MPI_TT<IdxT>::dtype,
                                      ProcCL::MPI_TT<char>::dtype };

    ProcCL::AintT displace[3];
    displace[0]= ProcCL::Get_address(&tmp.flag)-sAddr;
    displace[1]= ProcCL::Get_address(&tmp.dof) -sAddr;
    displace[2]= ProcCL::Get_address(&tmp.val) -sAddr;

    IdxTypeMPIT_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(IdxTypeMPIT_);
}

template<typename ExCL>
  bool FMTransferCL<ExCL>::HasDist(const ReprTetraCL& tetra) const
{
    for (Uint i=0; i<4; ++i)
        if (ex_.IsDist(tetra[i]))
            return true;
    return false;
}

template<typename ExCL>
  typename FMTransferCL<ExCL>::ProcSet FMTransferCL<ExCL>::GetProcs(const ReprTetraCL& tetra) const
{
    ProcSet procs;
    for (Uint i=0; i<4; ++i){
        typename ExCL::ProcNumCT onProc=ex_.GetProcs(tetra[i]);
        for (typename ExCL::ProcNumCT::const_iterator it(onProc.begin()), end(onProc.end()); it!=end; ++it)
            procs.insert(*it);
    }
    return procs;
}

template<typename ExCL>
  void FMTransferCL<ExCL>::MarkTetraForXfer(const ReprTetraCL& tetra)
/// This function put the tetra into a list, if the tetra is not stored in the loist before. The check for
/// equality only checks the DoF.
/// The list will be send, if XferTetras is called.
/// \param[in] tetra the tetra, that should be transfered
/// \pre The tetra must not contain any ghost DoF
{
    Assert(!tetra.HasGhost(),
            DROPSErrCL("FMTransferCL::MarkTetraForXfer: Not allowed to transfer Tetra that contains ghosts"),
            DebugParallelNumC);
    //if (std::find(ToXferTetra_.begin(), ToXferTetra_.end(), tetra)!=ToXferTetra_.end())
//     if (!IsIn(ToXferTetra_,tetra))
        ToXferTetra_.push_back(tetra);
}

template<typename ExCL>
  void FMTransferCL<ExCL>::XferTetras()
/// All Tetras that are marked by the function MarkTetraForXfer are now transfered.
{
    // init neighs_, if has been changed
    neighs_ = ex_.GetNeighbors();
    // Collect DoF that should be send and sort them by proc number
    // We also creates a mapping from DoF to the index of the DoF in the local vector
    typedef std::vector<IdxT>        DoFVec;            // set of DoF that should be transfered
    typedef std::vector<ReprTetraCL> TetraVec;          // vector of tetras that should be send
    typedef std::map<int, DoFVec>    ToSendDoFCT;       // container of DoF that should be send
    typedef std::map<int, TetraVec>  ToSendTetraCT;     // container of tetras that should be send
    typedef std::map<IdxT,int>       MapDoFSetT;        // mapping DoF -> index in DoFVec
    typedef std::map<int,MapDoFSetT> MapDoFT;           // mapping for a proc

    ToSendDoFCT   toProcSendDoF;
    ToSendTetraCT toProcSendTetra;
    MapDoFT       DoFMap;

    for (typename XferTetraList::const_iterator it(ToXferTetra_.begin()), end(ToXferTetra_.end()); it!=end; ++it){
        ProcSet onProc=GetProcs(*it);
        for (ProcSet::const_iterator proc(onProc.begin()), pend(onProc.end()); proc!=pend; ++proc){
              // insert Tetra in the SendList for Tetra
            toProcSendTetra[*proc].push_back(*it);
            for (Uint i=0; i<4; ++i){
                  // check is DoF is allready stored (if DoFVec is empty, this works too!)
                IdxT Nr = (*it)[i];
                bool notStored= std::find(toProcSendDoF[*proc].begin(), toProcSendDoF[*proc].end(), Nr)==toProcSendDoF[*proc].end();
                if (notStored){
                    DoFMap[*proc][Nr] = toProcSendDoF[*proc].size();
                    toProcSendDoF[*proc].push_back(Nr);
                }
            }
        }
    }

    int NumNeighProcs = neighs_.size();
    CoupFlagIdxS  **sendIdx     = new CoupFlagIdxS*[NumNeighProcs];
    CoupCoordValS **sendCoordVal= new CoupCoordValS*[NumNeighProcs];
    ProcCL::RequestT  *regIdx   = new ProcCL::RequestT[NumNeighProcs];
    ProcCL::RequestT  *regCoord = new ProcCL::RequestT[NumNeighProcs];

    int procPos=0;

      // for all neighbor procs send tetra and dofs
    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc, ++procPos){
        const int countDoF  = toProcSendDoF[*proc].size();
        const int countTetra= toProcSendTetra[*proc].size();
        Assert(countTetra>0, DROPSErrCL("FMTransferCL::XferTetras: Each neighbor procs should gets an message"),
               DebugParallelNumC);

        // allocate mem for sending
        sendIdx[procPos]     = new CoupFlagIdxS[4*countTetra];
        sendCoordVal[procPos]= new CoupCoordValS[countDoF];

          // Fill sendIdx
        for (int tetra=0; tetra<countTetra; ++tetra){
            for (int dof=0; dof<4; ++dof){
                CoupFlagIdxS tmp;
                IdxT Nr=toProcSendTetra[*proc][tetra][dof];   // proc, number of tetra, DoF in tetra
                tmp.flag= (ex_.IsOnProc(Nr,*proc) ? 0 : 1);
                tmp.dof = (tmp.flag ? DoFMap[*proc][Nr] : ex_.GetExternalIdxFromProc(Nr,*proc));
                sendIdx[procPos][tetra*4+dof] = tmp;
            }
        }

          // Fill sendCoordVal
        for (int i=0; i<countDoF; ++i){
            CoupCoordValS tmp;
            tmp.coord[0]= Coord_[toProcSendDoF[*proc][i]][0];
            tmp.coord[1]= Coord_[toProcSendDoF[*proc][i]][1];
            tmp.coord[2]= Coord_[toProcSendDoF[*proc][i]][2];
            tmp.val     = v_.Data[toProcSendDoF[*proc][i]];
            tmp.typ     = Typ_[toProcSendDoF[*proc][i]];
            sendCoordVal[procPos][i]=tmp;
        }

          // real send
        regIdx[procPos]  = ProcCL::Isend(sendIdx[procPos], 4*countTetra, FlagIdxMPIT_, *proc, tag_);
        regCoord[procPos]= ProcCL::Isend(sendCoordVal[procPos], countDoF, CoordValMPIT_, *proc, tag_+1);
    }

     // recieve tetras and dofs
    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc){
          // recieve indices of tetras
        ProcCL::StatusT stat;
        ProcCL::Probe(*proc, tag_, stat);
        int countTetra = ProcCL::GetCount(stat, FlagIdxMPIT_);
        Assert(countTetra%4==0, DROPSErrCL("FMTransferCL::XferTetras: Wrong number of recieved indices"), DebugParallelNumC);
        countTetra = countTetra/4;
        CoupFlagIdxS  *recvIdx     = new CoupFlagIdxS[4*countTetra];
        ProcCL::Recv(recvIdx, 4*countTetra, FlagIdxMPIT_, *proc, tag_);

          // recieve coords and vals
        ProcCL::Probe(*proc, tag_+1, stat);
        const int countDoF = ProcCL::GetCount(stat, CoordValMPIT_);
        CoupCoordValS *recvCoordVal= new CoupCoordValS[countDoF];
        ProcCL::Recv(recvCoordVal, countDoF, CoordValMPIT_, *proc, tag_+1);

        const Uint startGhost=DataExt_.size();
        Assert(startGhost==CoordExt_.size(), DROPSErrCL("FMTransferCL::XferTetras: Ghost Vals and ghost coords do not have same length"), DebugParallelNumC);

          // Put Ghost values and coodinates into vector
        for (int i=0; i<countDoF; ++i){
            DataExt_.push_back(recvCoordVal[i].val);
            CoordExt_.push_back(MakePoint3D(recvCoordVal[i].coord[0], recvCoordVal[i].coord[1], recvCoordVal[i].coord[2]));
            TypExt_.push_back(recvCoordVal[i].typ);
        }

          // put Tetras to neighbors
        for (int tetra=0; tetra<countTetra; ++tetra){
              // create tetra
            ReprTetraCL t;
            for (Uint j=0; j<4; ++j){
                bool ghost=recvIdx[tetra*4+j].flag;
                IdxT Nr   = (ghost ? recvIdx[tetra*4+j].dof+startGhost : recvIdx[tetra*4+j].dof);
                t.Set(j,FMIdxT(ghost,Nr));
            }
              // put tetra to neighborhood of DoF
            for (Uint j=0; j<4; ++j)
                if (!t.IsGhost(j))
                    neigh_[t[j]].push_back(t);
        }

          // send start index for ghosts back
        ProcCL::Send(&startGhost, 1, *proc, tag_+2);

          // free recieve buffer
        delete[] recvIdx;
        delete[] recvCoordVal;
    }

    // free send buffers
    ProcCL::WaitAll(NumNeighProcs, regIdx);
    ProcCL::WaitAll(NumNeighProcs, regCoord);

    delete[] regIdx;
    delete[] regCoord;

    for (int i=0; i<NumNeighProcs; ++i){
        delete[] sendIdx[i];
        delete[] sendCoordVal[i];
    }

    delete[] sendIdx;
    delete[] sendCoordVal;


      // recieve Ghost-Indices
    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc){
        Uint startGhost=0;
        ProcCL::Recv(&startGhost, 1, *proc, tag_+2);
        IdxT pos=startGhost;
        for(Uint i=0; i<toProcSendDoF[*proc].size(); ++i){
            IdxT Nr=toProcSendDoF[*proc][i];
            if (!ex_.IsOnProc(Nr,*proc))
                HasGhost_[Nr].push_back(GhostIdxT(*proc, pos++));
        }
    }
}

template<typename ExCL>
  void FMTransferCL<ExCL>::ChangedHasGhost(IdxT dof)
/// If a DoF that has at least on ghost value on another proc tell these procs about the change.
/// Therefore put the representative tetras into a list. The transfer will be performed if the
/// function UpdateGhosts() will be called.
{
    Assert(!HasGhost_[dof].empty() || ex_.IsDist(dof), DROPSErrCL("FMTransferCL::ChangedHasGhost: ghostlist is empty"),DebugParallelNumC);
    if (!IsIn(ToUpdateGhost_,dof)){
//         if (ProcCL::MyRank()==1){
//             std::cout << "Send Change of dist dof (loc " <<dof<<", ext "<<(ex_.IsOnProc(dof,0)?ex_.GetExternalIdxFromProc(dof,0):-1)
//                       << ") mit Wert " << v_.Data[dof]<< std::endl;
//         }
        ToUpdateGhost_.push_back(dof);
    }
}

template<typename ExCL>
  void FMTransferCL<ExCL>::UpdateGhosts()
/// Send changed values of the DoF to the procs, that has these DoFs and update values
{
    typedef std::pair<IdxT,IdxT>                    CoupLocExtT; // Coupling of local dof extern ghost index
    typedef std::map <int, std::list<CoupLocExtT> > ProcGhCT;
    ProcGhCT sendGhosts;

    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc){
        sendGhosts[*proc];
    }

//     int tmpcount1=0, tmpcount2=0, tmpcount3=0, tmpcount4=0;
      // Collect DoFs that should be send
    for (typename IdxListT::const_iterator dof(ToUpdateGhost_.begin()), end(ToUpdateGhost_.end()); dof!=end; ++dof){
          // Collect distributed DoFs

        if ( ex_.IsDist(*dof) ){
            typename ExCL::ProcNumCT onProc=ex_.GetProcs(*dof);
            for (typename ExCL::ProcNumCT::iterator it(onProc.begin()), end(onProc.end()); it!=end; ++it){
                sendGhosts[*it].push_back(CoupLocExtT(*dof,ex_.GetExternalIdxFromProc(*dof,*it)));
//                 ++tmpcount1;
            }
        }
          // Collect HasGhosts DoFs for sending
        for (GhostsList::const_iterator it(HasGhost_[*dof].begin()), end(HasGhost_[*dof].end()); it!=end; ++it){
            sendGhosts[it->first].push_back(CoupLocExtT(*dof,it->second));
//             ++tmpcount2;
        }
    }

      // send updated ghosts
    for (ProcGhCT::const_iterator send(sendGhosts.begin()), end(sendGhosts.end()); send!=end; ++send){
          // Create send buf
        const int count= send->second.size();
        CoupIdxValS *sendBuf = new CoupIdxValS[count];

          // put vals into sendbuf
        int pos=0;
        for (std::list<CoupLocExtT>::const_iterator dof(send->second.begin()), dend(send->second.end()); dof!=dend; ++dof){
            sendBuf[pos].flag = ( ex_.IsOnProc(dof->first,send->first) ? 0 : 1);
//             if (!sendBuf[pos].flag)
//                 ++tmpcount3;
//             else
//                 ++tmpcount4;
            sendBuf[pos].dof  = dof->second;
            sendBuf[pos++].val= v_.Data[dof->first];
        }
//         if (ProcCL::MyRank()==1)
//             std::cout << ">>> Die Counter lauten: "<<tmpcount1<<", "<<tmpcount2<<", "<<tmpcount3<<", "<<tmpcount4<<"\n";
//         tmpcount1=0;
//         tmpcount2=0;

          // Send
        ProcCL::Send(sendBuf, count, IdxValMPIT_, send->first, tag_+3);
        delete[] sendBuf;
    }

     // Ghostinformation have been send, so list is not needed any more
    ToUpdateGhost_.clear();

      // recieve messages
    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc){
          // check how many ghosts have to be updated
        ProcCL::StatusT stat;
        ProcCL::Probe(*proc, tag_+3, stat);
        const int count = ProcCL::GetCount(stat, IdxValMPIT_);

          // recieving
        CoupIdxValS *recvBuf = new CoupIdxValS[count];
        ProcCL::Recv(recvBuf, count, IdxValMPIT_, *proc, tag_+3);

          // handle recieved data
        for (int i=0; i<count; ++i){
            if (recvBuf[i].flag){
                DataExt_[recvBuf[i].dof] = recvBuf[i].val;
//                 ++tmpcount1;
            }
            else{
//                 IF_MASTER
//                         std::cout << "Updating value ("<<recvBuf[i].dof<<") from "<<v_.Data[recvBuf[i].dof]
//                                   <<" to ("<<ex_.GetExternalIdxFromProc(recvBuf[i].dof,*proc)<<") "
//                                   <<recvBuf[i].val<<std::endl;
                v_.Data[recvBuf[i].dof] = std::min(recvBuf[i].val,v_.Data[recvBuf[i].dof]);
//                 ++tmpcount2;
            }
        }

        delete[] recvBuf;
    }
//     IF_MASTER
//             std::cout << "Proc 0: Counter: "<<tmpcount1<<", "<<tmpcount2<<"\n";
}


template<typename ExCL>
  void FMTransferCL<ExCL>::MarkChangeType(IdxT dof)
{
    Assert(!HasGhost_[dof].empty() || ex_.IsDist(dof), DROPSErrCL("FMTransferCL::MarkChangeType DoF is not distributed"),DebugParallelNumC);
    if (!IsIn(ToUpdateType_,dof))
        ToUpdateType_.push_back(dof);
}

template<typename ExCL>
  void FMTransferCL<ExCL>::UpdateType()
{
    typedef std::pair<IdxT,IdxT>                    CoupLocExtT; // Coupling of local dof extern ghost index
    typedef std::map <int, std::list<CoupLocExtT> > ProcTypeCT;
    ProcTypeCT sendTypes;

    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc){
        sendTypes[*proc];
    }

      // Collect updated DoFs for sending
    for (typename IdxListT::const_iterator dof(ToUpdateType_.begin()), end(ToUpdateType_.end()); dof!=end; ++dof){
        if ( ex_.IsDist(*dof) ){
            typename ExCL::ProcNumCT onProc=ex_.GetProcs(*dof);
            for (typename ExCL::ProcNumCT::iterator it(onProc.begin()), end(onProc.end()); it!=end; ++it)
                sendTypes[*it].push_back(CoupLocExtT(*dof,ex_.GetExternalIdxFromProc(*dof,*it)));
        }

        for (GhostsList::const_iterator it(HasGhost_[*dof].begin()), end(HasGhost_[*dof].end()); it!=end; ++it)
            sendTypes[it->first].push_back(CoupLocExtT(*dof,it->second));
    }

      // send updated DoF
    for (ProcTypeCT::const_iterator send(sendTypes.begin()), end(sendTypes.end()); send!=end; ++send){
          // Create send buf
        const int count= send->second.size();
        CoupIdxTypeS *sendBuf = new CoupIdxTypeS[count];

          // put vals into sendbuf
        int pos=0;
        for (std::list<CoupLocExtT>::const_iterator dof(send->second.begin()), dend(send->second.end()); dof!=dend; ++dof){
            sendBuf[pos].flag = ( ex_.IsOnProc(dof->first,send->first) ? 0 : 1);
            sendBuf[pos].dof= dof->second;
            sendBuf[pos++].val = Typ_[dof->first];
        }

          // Send
        ProcCL::Send(sendBuf, count, IdxTypeMPIT_, send->first, tag_+4);
        delete[] sendBuf;
    }

     // Updateinformation have been send, so list is not needed any more
    ToUpdateType_.clear();

      // recieve messages
    for (typename ExCL::ProcNumCT::const_iterator proc(neighs_.begin()), pend(neighs_.end()); proc!=pend; ++proc){
          // check how many ghosts have to be updated
        ProcCL::StatusT stat;
        ProcCL::Probe(*proc, tag_+4, stat);
        const int count = ProcCL::GetCount(stat, FlagIdxMPIT_);

          // recieving
        CoupIdxTypeS *recvBuf = new CoupIdxTypeS[count];
        ProcCL::Recv(recvBuf, count, IdxTypeMPIT_, *proc, tag_+4);

          // handle recieved data
        for (int i=0; i<count; ++i){
            if (recvBuf[i].flag){
                Assert(TypExt_[recvBuf[i].dof]>=recvBuf[i].val, DROPSErrCL("FMTransferCL::UpdateType: Making Type smaller!"),DebugParallelNumC);
                TypExt_[recvBuf[i].dof] = recvBuf[i].val;
            }
            else{
                Assert(Typ_[recvBuf[i].dof]>=recvBuf[i].val, DROPSErrCL("FMTransferCL::UpdateType: Making Type smaller!"),DebugParallelNumC);
                Typ_[recvBuf[i].dof] = recvBuf[i].val;
            }
        }
        delete[] recvBuf;
    }
}

template<typename ExCL>
  IdxT FMTransferCL<ExCL>::FindTrial(IdxT locDoF) const
  /// Search over all procs for DoF minimal value.
{
    double myVal=(locDoF!=NoIdx ? v_.Data[locDoF] : 1e99);
    double *globVal = new double[ProcCL::Size()];
    Gather(myVal, globVal,-1);

    int minproc=-1;
    double minVal=1e99;
    IdxT ret=NoIdx;

    for (int p=0; p<ProcCL::Size(); ++p){
        if (globVal[p]<minVal){
            minVal=globVal[p];
            minproc=p;
        }
    }

    if (minVal==1e99)      // No DoF is left
        ret = StopIdx;
    else
        ret = (minproc==ProcCL::MyRank() ? locDoF : NoIdx);

    return ret;
}


//------------------------------
// P A R  F A S T M A R C H  C L
//------------------------------

template<typename ExCL>
  ParFastMarchCL<ExCL>::ParFastMarchCL( MultiGridCL& mg, VecDescCL& v, ExCL& ex)
    : MG_(mg), ex_(ex), size_(v.RowIdx->NumUnknowns), v_(v),
      Typ_(Far, size_), Coord_(0), HasGhost_(0), CoordExt_(0), DataExt_(0),
      FMTransfer_(v_,Coord_,neigh_,Typ_,DataExt_,CoordExt_,TypExt_,HasGhost_,ex),
      initpar_(false)
/// \param mg The multigrid on whitch the fastmarching algorithm should be performed
/// \param v  Values of the levelset function
/// \param ex An ExchangeCL for v
{
}


template<typename ExCL>
  Point3DCL ParFastMarchCL<ExCL>::GetCoord(const FMIdxT &dof) const
{
    return (dof.first ? CoordExt_[dof.second] : Coord_[dof.second]);
}


template<typename ExCL>
  double ParFastMarchCL<ExCL>::GetData(const FMIdxT& dof) const
{
    return (dof.first ? DataExt_[dof.second] : v_.Data[dof.second]);
}


template<typename ExCL>
  void ParFastMarchCL<ExCL>::SetData(double val, IdxT dof)
{
    v_.Data[dof]= std::min(std::abs(val),std::abs(v_.Data[dof]));
    // if a value on a DoF that is marked as HasGhost changes, tell other procs about this change
    if (HasGhost(dof) || ex_.IsDist(dof))
        FMTransfer_.ChangedHasGhost(dof);
}

template<typename ExCL>
  byte ParFastMarchCL<ExCL>::GetType(const FMIdxT& dof) const
{
    return (dof.first ? TypExt_[dof.second] : Typ_[dof.second]);
}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::SetType(FMMark mark, IdxT dof)
{
    Typ_[dof]= mark;
    if (Typ_[dof]!=mark && HasGhost(dof))
        FMTransfer_.MarkChangeType(dof);
}


template<typename ExCL>
  void ParFastMarchCL<ExCL>::UpdateHasGhost(IdxT dof)
// Check if the given DoF has ghost the put this dof into update list
{
    Assert(HasGhost(dof),DROPSErrCL("ParFastMarchCL::UpdateHasGhost: Given DoF has no ghost"),DebugParallelNumC);
    if (HasGhost(dof) || ex_.IsDist(dof))
        FMTransfer_.ChangedHasGhost(dof);
}


template<typename ExCL>
        bool ParFastMarchCL<ExCL>::HasGhost(IdxT i) const
{
    return !HasGhost_[i].empty();
}


template<typename ExCL>
  void ParFastMarchCL<ExCL>::InitCoord()
/// Iteratate over all DoF and put coordinates into local list
{
    Comment("Init Coord local, backup old values and make all values greater zero\n", DebugParallelNumC);
    const Uint idx= v_.RowIdx->GetIdx();                // number of index
    Coord_.resize( size_);
    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(), end=MG_.GetTriangVertexEnd();
         it!=end; ++it)
        Coord_[it->Unknowns(idx)]= it->GetCoord();
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(), end=MG_.GetTriangEdgeEnd();
         it!=end; ++it)
        Coord_[it->Unknowns(idx)]= GetBaryCenter( *it);

    Old_.resize( size_);
    Old_= v_.Data;
    for (IdxT i=0; i<v_.Data.size(); ++i)
        v_.Data[i] = std::abs(v_.Data[i]);
}


template<typename ExCL>
  void ParFastMarchCL<ExCL>::InitNeighAndPar()
{
    Comment("Init Neigh and Par\n", DebugParallelNumC);
    const Uint idx= v_.RowIdx->GetIdx();                // number of index
    IdxT Numb[10];                                      // DoF for a tetra
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);   // regular ref of a tetra

    neigh_.resize( size_);                              // allocate mem for neighborship
    HasGhost_.resize(size_);                            // allocate mem for hasGhost marks

    for (MultiGridCL::TriangTetraIteratorCL it=MG_.GetTriangTetraBegin(), end=MG_.GetTriangTetraEnd();
         it!=end; ++it)
    {
          // collect data on all DoF
        for (int v=0; v<10; ++v){
            if (v<4)
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            else
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);
        }

          // for all regular children do
        for (int ch=0; ch<8; ++ch){
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            ReprTetraCL t;
            bool containsDist=false;

              // Fill tetra and search for distributed DoF
            for (int vert=0; vert<4; ++vert){
                const IdxT Nr= Numb[ data.Vertices[vert]];
                t.Set(vert,FMIdxT(false,Nr));
                containsDist = containsDist || ex_.IsDist(Nr);
            }

              // Mark tetras for transfer (if they a vert was detected to be distributed)
            if (containsDist)
                FMTransfer_.MarkTetraForXfer(t);

              // init neigh_
            for (int vert= 0; vert<4; ++vert)
                neigh_[t[vert]].push_back( t);
        }
    }
    FMTransfer_.XferTetras();
    initpar_= true;
}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::InitZero( bool ModifyZero)
/// This is the same function as in serial FastMarchCL.
{
    Comment("Init Zero\n", DebugParallelNumC);
    // Knoten an der Phasengrenze als Finished markieren
    // und Distanz zur Phasengrenze bestimmen (falls ModifyZero)
    const Uint idx= v_.RowIdx->GetIdx();
    int        sign[10];
    int        num_sign[3]; // - 0 +
    IdxT       Numb[10];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);
      // for all tetras
    for (MultiGridCL::TriangTetraIteratorCL it(MG_.GetTriangTetraBegin()), end(MG_.GetTriangTetraEnd()); it!=end; ++it)
    {
        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            Numb[v]= v<4 ? it->GetVertex(v)->Unknowns(idx)
                : it->GetEdge(v-4)->Unknowns(idx);
            sign[v]= std::abs(Old_[Numb[v]])<1e-8 ? 0 : (Old_[Numb[v]]>0 ? 1 : -1);
            if (sign[v]==0)
                SetType(Finished,Typ_[Numb[v]]);
        }

        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];

            const bool intersec= (num_sign[0]*num_sign[2]!=0); // Vorzeichenwechsel

            if (!intersec) continue;

            if (!ModifyZero)
            {
                for (int vert= 0; vert<4; ++vert)
                {
                    const IdxT Nr= Numb[data.Vertices[vert]];
                    //Typ_[Nr]= Finished;
                    //v_.Data[Nr]= std::abs( Old_[Nr]);
                    SetType(Finished, Nr);
                    SetData(std::abs(Old_[Nr]), Nr);
                }
                continue;
            }

            // from here on intersec and ModifyZero are true
            Point3DCL Schnitt[4];
            int num= 0;

            // Compute intersections with edges of tetra
            for (int vert= 0; vert<4; ++vert)
                if (sign[data.Vertices[vert]]==0)
                    Schnitt[num++]= Coord_[Numb[data.Vertices[vert]]];

            for (int edge= 0; edge<6 && num<4; ++edge)
            {
                const Ubyte v1= data.Vertices[ VertOfEdge( edge, 0)],
                v2= data.Vertices[ VertOfEdge( edge, 1)];
                if (sign[v1]*sign[v2] == -1) // Vorzeichenwechsel auf edge
                {
                    const IdxT Nr1= Numb[v1],
                    Nr2= Numb[v2];
                    const double bary= Old_[Nr1]/(Old_[Nr1]-Old_[Nr2]);
                    Schnitt[num++]= (1-bary)*Coord_[Nr1] + bary*Coord_[Nr2];
                }
            }
/*
            fil << "geom {OFF " << num << " 1 0\n";
            for (int i=0; i<num; ++i)
            {
            for (int j=0; j<3; ++j)
            fil << Schnitt[i][j] << ' ';
            fil << '\n';
        }
            if (num==3)
            fil << "3 0 1 2";
            else
            fil << "4 0 1 3 2";
            fil << "\n}\n";
*/
            if (num<3) throw DROPSErrCL("FastMarchCL::InitZero: intersection missing");

            for (int repeat=0; repeat<num-2; ++repeat)
            { // fuer num==4 (Schnitt ABDC ist viereckig)
              // zwei Dreiecke ABC + DBC betrachten
                if (repeat) Schnitt[0]= Schnitt[3];

                const Point3DCL a= Schnitt[1] - Schnitt[0],
                b= Schnitt[2] - Schnitt[0];

                for (int vert=0; vert<4; ++vert)
                {
                    if (sign[data.Vertices[vert]]==0) continue;

                    const IdxT Nr= Numb[data.Vertices[vert]];
                    const Point3DCL Crd= Coord_[Nr],
                    c=   Crd - Schnitt[0];
                    double dist= std::min( c.norm(), (Crd-Schnitt[1]).norm());
                    dist= std::min( dist, (Crd-Schnitt[2]).norm());

                    const double bary1= inner_prod(a,c)/a.norm_sq(),
                    bary2= inner_prod(b,c)/b.norm_sq();
                    if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
                    {
                        const Point3DCL lotfuss= (1-bary1-bary2)*Schnitt[0] + bary1*Schnitt[1] + bary2*Schnitt[2];
                        dist= std::min( dist, (lotfuss - Crd).norm());
                    }

                    if (Typ_[Nr] != Finished)
                    {
                        //Typ_[Nr]= Finished;
                        //v_.Data[Nr]= dist;
                        SetType(Finished, Nr);
                        SetData(dist,Nr);
                    }
                    else{
                        SetData(std::min( dist, v_.Data[Nr]), Nr);
//                         v_.Data[Nr]= std::min( dist, v_.Data[Nr]);
                    }
                }
            }
        }
    }
    UpdateGhost();
//fil << "}\n";
}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::InitClose()
  /// Mark DoF that are close to Finished marked DoF as Close and update value of levelset-function
{
    Comment("Init Close\n", DebugParallelNumC);
    // an Finished angrenzende Knoten mit Close markieren und dort v_ updaten
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    for (Uint Nr=0; Nr<size_; ++Nr)
        if (Typ_[Nr]==Finished)
            for (Uint n=0; n<neigh_[Nr].size(); ++n)
                for (int j=0; j<4; ++j)
                    Update(neigh_[Nr][n].Get(j));
    UpdateGhost();
}

template<typename ExCL>
  double ParFastMarchCL<ExCL>::CompValueProj( IdxT Nr, int num, const FMIdxT upd[3]) const
// Nr is local, upd may contain ghosts
{
    double val= 1e99;
    switch (num)
    {
        case 2: // Projektion auf Edge
        {
            const Point3DCL a= GetCoord(upd[1]) - GetCoord(upd[0]); //Coord_[upd[1]] - Coord_[upd[0]];
            const Point3DCL b= Coord_[Nr] - GetCoord(upd[0]);        //Coord_[  Nr  ] - Coord_[upd[0]];
            const double bary= inner_prod(a,b)/a.norm_sq();
            if (bary>=0 && bary<=1){
//              const Point3DCL lotfuss= (1-bary)*Coord_[upd[0]] + bary*Coord_[upd[1]];
                const Point3DCL lotfuss= (1-bary)*GetCoord(upd[0]) + bary*GetCoord(upd[1]);
//              const double y= (1-bary)*v_.Data[upd[0]] + bary*v_.Data[upd[1]];
                const double y= (1-bary)*GetData(upd[0]) + bary*GetData(upd[1]);
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
        break;

        case 3: // Projektion auf Face
        {
            const Point3DCL a= GetCoord(upd[1]) - GetCoord(upd[0]); // Coord_[upd[1]] - Coord_[upd[0]];
            const Point3DCL b= GetCoord(upd[2]) - GetCoord(upd[0]); // Coord_[upd[2]] - Coord_[upd[0]];
            const Point3DCL c= Coord_[  Nr  ]   - GetCoord(upd[0]); // Coord_[  Nr  ] - Coord_[upd[0]];
            const double bary1= inner_prod(a,c)/a.norm_sq(),
            bary2= inner_prod(b,c)/b.norm_sq();
            if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
            {
//              const Point3DCL lotfuss= (1-bary1-bary2)*Coord_[upd[0]] + bary1*Coord_[upd[1]] + bary2*Coord_[upd[2]];
                const Point3DCL lotfuss= (1-bary1-bary2)*GetCoord(upd[0]) + bary1*GetCoord(upd[1]) + bary2*GetCoord(upd[2]);
//              const double y= (1-bary1-bary2)*v_.Data[upd[0]] + bary1*v_.Data[upd[1]] + bary2*v_.Data[upd[2]];
                const double y= (1-bary1-bary2)*GetData(upd[0]) + bary1*GetData(upd[1]) + bary2*GetData(upd[2]);
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
    }

    return val;
}

template<typename ExCL>
  IdxT ParFastMarchCL<ExCL>::FindTrial() const
// Returns DoF in Close set with minimal value if this DoF is located on the calling proc. Otherwise
// return NoIdx.
{
    double min= 1e99;
    IdxT min_idx= NoIdx;

    for (typename std::set<IdxT>::const_iterator it(Close_.begin()), end(Close_.end()); it!=end; ++it)
    {
        if (v_.Data[*it]<=min)
        {
            min= v_.Data[*it];
            min_idx= *it;
        }
    }

    return FMTransfer_.FindTrial(min_idx);
}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::Update(const FMIdxT& dof)
{
    if (dof.first)          // No Updates on ghosts
        return;
    if (Typ_[dof.second]==Finished)
        return;             // No Update on finished DoF

    IdxT NrI= dof.second;   // NrI is local
    FMIdxT upd[3];
    double minval= (Typ_[NrI]==Close ? v_.Data[NrI] : 1e99);

    for (Uint n=0; n<neigh_[NrI].size(); ++n){
        int num= 0;
        for (int j=0; j<4; ++j){
            const FMIdxT dofJ= neigh_[NrI][n].Get(j);       // j may be on ghost

            if (GetType(dofJ) == Finished){
                upd[num++]= dofJ;
                minval= std::min( minval, GetData(dofJ) + (GetCoord(dofJ)-Coord_[NrI]).norm());
            }
        }
        minval= std::min( minval, CompValueProj( NrI, num, upd));
    }

    SetData(minval, NrI); // v_.Data[NrI]= minval;
    if (Typ_[NrI] != Close){
        Close_.insert( NrI);
        SetType(Close, NrI); //Typ_[NrI]= Close;
    }

}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::RestoreSigns()
{ // restore signs of v_
    for (IdxT i=0, N= Old_.size(); i<N; ++i)
        if (Old_[i]<0){
//             IF_MASTER
//               std::cout << "Mache Eintrag "<<i<<" negativ\n";
            v_.Data[i]*= -1;
        }
}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::UpdateGhost()
{
    FMTransfer_.UpdateType();
     FMTransfer_.UpdateGhosts();
}

template<typename ExCL>
  void ParFastMarchCL<ExCL>::Reparam( bool ModifyZero)
{
    InitCoord();
    InitNeighAndPar();
    InitZero( ModifyZero);
    InitClose();

    Comment("Reparametrization\n", DebugParallelNumC);
    IdxT next, step=0;

    while ((next= FindTrial()) != FMTransfer_.StopIdx)
    {
        if (next!=NoIdx){                                       // minimal DoF is on calling proc
            Close_.erase( next);
            Typ_[next]= Finished;

            std::set<FMIdxT> neighVerts;
            for (Uint n=0; n<neigh_[next].size(); ++n)
            { // collect all neighboring verts in neighVerts
                for (Uint i=0; i<4; ++i)
                    neighVerts.insert( neigh_[next][n].Get(i));
            }
            for (std::set<FMIdxT>::const_iterator it(neighVerts.begin()), end(neighVerts.end()); it!=end; ++it)
            { // update all neighboring verts, mark as Close
                Update( *it);
            }
            neigh_[next].clear(); // will not be needed anymore
        }
        UpdateGhost();
        ++step;
    }

    Assert(AllPos(), DROPSErrCL("ParFastMarchCL::Reparam: Not all values positive!"),DebugParallelNumC);
    RestoreSigns();
    Assert(ex_.IsAcc(v_.Data),DROPSErrCL("ParFastMarchCL::Reparam: Levelset is not accumulated any more!"), DebugParallelNumC);

    Comment("Reparametrization needs "<<step<<" steps\n",DebugParallelNumC);
}

template<typename ExCL>
  bool ParFastMarchCL<ExCL>::CheckHasGhosts()
{
    Comment("CheckHasGhosts()\n", DebugParallelNumC);
    for (IdxT i=0; i<HasGhost_.size(); ++i){
        for (typename GhostsList::iterator it(HasGhost_[i].begin()), end(HasGhost_[i].end()); it!=end; ++it){
            if (ex_.IsOnProc(i,it->first)){
                std::cerr << "["<<ProcCL::MyRank()<<"] Index "<<i<<" ist als Ghost auf Proc "<<it->first<<" vorhanden, dort aber auch als dist vorliegend\n";
                throw DROPSErrCL("Distributed Ghost!");
            }
        }
    }
    return true;
}

template<typename ExCL>
  bool ParFastMarchCL<ExCL>::AllPos()
{
    for (IdxT i=0; i<v_.Data.size(); ++i){
        if (v_.Data[i]<0){
            std::cout << "Data auf DoF " <<i<<" ist kleiner als 0\n";
            throw DROPSErrCL("Found negative value");
        }
    }
    return true;
}

} // end of namespace DROPS
