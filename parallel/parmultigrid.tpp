//**************************************************************************
// File:    parmultigrid.h                                                 *
// Content: Class that constitute the parallel multigrid                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, RZ RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   November, 14th, 2005                                           *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file parmultigrid.tpp

namespace DROPS{

template<class IterT>
void ShowSimplex( IterT begin, IterT end, DDD_GID gid, char *mesg, int proc= -1)
{
    if (proc!=DDD_InfoMe() && proc!=-1)
        return;
    for (IterT it= begin; it!=end; ++it)
        if (it->GetGID()==gid)
    {
        std::cout << "Show " << mesg << ":\n";
        it->DebugInfo( std::cout);
    }
}


/// \brief Check if a simplex lies on a boundary between procs
template<class SimplexT>
  bool ParMultiGridCL::IsOnProcBnd( const SimplexT* s)
{
    const DDD_HDR hdr= s->GetHdr();

    if (DDD_InfoIsLocal( hdr ) )
        return false;
    Uint counter= 0;
    for( int* proclist= DDD_InfoProcList( hdr); *proclist!=-1; proclist+= 2)
        if (*(proclist+1)!=PrioVGhost)
            ++counter;
    return counter>=2;
}

/// \brief Change the priority of a simplex
template<class SimplexT>
  void ParMultiGridCL::PrioChange(SimplexT* const Tp, Uint Prio)
{
//  Assert(PrioChangeMode || TransferMode, DROPSErrCL("ParMultiGridCL: PrioChange: There must be an active Xfer- or PrioChange-Mode, to run this procedure"), DebugParallelC);
    if (Tp->GetPrio()!=Prio)
    {
        DDD_XferPrioChange( Tp->GetHdr(), Prio);
        Tp->SetPrio(Prio);
    }
}

/// \brief Recieve Information, if all Vector-Describers are recieved
bool ParMultiGridCL::VecDescRecv()
{
    return !_VecDesc.empty();
}

/// \brief Assign a VecDescCL to ParMultiGridCL
template<typename BndT>
  void ParMultiGridCL::AttachTo(VecDescCL *x, const BndT* bndp)
/// If \a x is already known (checked by sysnum), this function overwrite the
/// old VecDescCL, else the pointer is "pushed_back" to the known VecDescCLs.
/// If the pointer is already known, this routine assumes, that the boundary
/// conditions did not changed, so there is no need to call
/// AttachTo( const IdxDescCL*, const BndDataCL*) again.
{
    // determine position, where to store the pointer, and store it
    size_t pos= GetStorePos(x->RowIdx);
    if ( pos==_VecDesc.size() ){
        _VecDesc.push_back( x);
        _ScalBnd.push_back( 0);
        _VecBnd.push_back( 0);
    }
    else {
        _VecDesc[pos]= x;
    }
    AttachTo( x->RowIdx, bndp);

    // set flags if unknowns exits on vertices, edges or tetras
    if (x->RowIdx->NumUnknownsVertex()>0) _UnkOnSimplex[0]=true;
    if (x->RowIdx->NumUnknownsEdge()>0)   _UnkOnSimplex[1]=true;
    if (x->RowIdx->NumUnknownsTetra()>0)  _UnkOnSimplex[2]=true;
}



/// \brief Get position where the IdxDesc is internally stored
Uint ParMultiGridCL::GetStorePos(const IdxDescCL* idxDesc)
/// Check (by comparing the sysnum), where the VecDescCL correspoding to
/// \a idxDesc is internally stored.
/// \return position within _VecDesc if \a idxDesc is known, else return the
///         size of _VecDesc
{
    Uint idx= idxDesc->GetIdx();
    Uint pos=0;
    while ( pos<_VecDesc.size() && _VecDesc[pos]->RowIdx->GetIdx()!=idx )
        ++pos;
    return pos;
}

/// \brief Check if unknowns on simplices are available
bool ParMultiGridCL::UnknownsOnSimplices()
{
    return _UnkOnSimplex[0]||_UnkOnSimplex[1]||_UnkOnSimplex[2];
}

/// \brief Gather unknowns on a simplex, that cannot be set by another proc
template<typename SimplexT>
int ParMultiGridCL::GatherInterpolValues (DDD_OBJ obj, void* buf)
/** For a detailt situation describtion look at ScatterInterpolValues*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    Uint idx   = _actualVec->RowIdx->GetIdx();
    Uint numUnk= _actualVec->RowIdx->NumUnknownsVertex();

    // if there are unknowns to the new index on the vertex, and this is an old
    // value, then put it into the buffer. Else mark the the unknown as invalide
    if (sp->Unknowns.Exist() && sp->Unknowns.Exist(idx) && sp->Unknowns.UnkRecieved(idx)){
        for (Uint i=0; i<numUnk; ++i){
            buffer[i].mark= true;
            buffer[i].idx = idx;
            buffer[i].val = _actualVec->Data[sp->Unknowns(idx)];
        }
    }
    else{
        buffer->mark= false;
        buffer->idx = idx;
    }

    return 0;
}

/// \brief Scatter unknowns on a simplex, that cannot be set on this proc
template<typename SimplexT>
int ParMultiGridCL::ScatterInterpolValues(DDD_OBJ obj, void* buf)
/** If a unknown on a distributed simplex has not been interpolated by a proc
    due to a lack of information, the information are exchanged here. */
{
    SimplexT* const sp= ddd_cast<SimplexT*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    Uint idx   = _actualVec->RowIdx->GetIdx();
    Uint numUnk= _actualVec->RowIdx->NumUnknownsVertex();         // unknowns on edges and vertices are the same

    // if there is an unknown to the new index on this vertex, that do not store
    // an old value and the recieved value is valide, then set this value and
    // mark the value as an old value
    if (sp->Unknowns.Exist() && sp->Unknowns.Exist(idx) && !sp->Unknowns.UnkRecieved(idx) && buffer->mark && buffer->idx==idx){
        for (Uint i=0; i<numUnk; ++i)
            _actualVec->Data[sp->Unknowns(idx)]= buffer[i].val;
        sp->Unknowns.SetUnkRecieved(idx);
    }

    return 0;
}

/****************************************************************************
* S E N D  U N K N O W N S                                                  *
*****************************************************************************
*   This procedure puts all unknowns of a vectoriel or scalar type into the *
*   the class AddVecCL or AddScalCL. The constructed array of Add<Type>CL   *
*   will be given to DDD and send to the recieving proc.                    *
*                                                                           *
*   The unknowns can be stored at two different places. If the unknown has  *
*   been recieved due to a killed ghost teraeder, the unknown can be found  *
*   in the buffer _RecvBuf. Otherwise it is stored in the vector, which has *
*   been given to the ParMultiGridCL before the refinement has been started.*
*                                                                           *
*   Parameter:                                                              *
*   - sp:   simplex, where to get the data                                  *
*   - type: scalar or vectoriel data                                        *
*   - buf:  buffer of the type Added<Type>CL                                *
*   - cnt:  number of unknowns (in terms of double)                         *
****************************************************************************/
template<class SimplexT>
  void ParMultiGridCL::SendUnknowns(SimplexT* sp, DDD_TYPE type, void *buf, int cnt)
{
    if (cnt==0) return;

    Uint counter=0;
    if (type == AddedScalCL::GetType()){                                    // scalar values are demanded
        AddedScalCL* buffer = static_cast<AddedScalCL*>(buf);               // cast buffer to AddedScalCL*
        if (sp->Unknowns.Exist()){                                          // check for unknowns on this simplex
            for (size_t i=0; i<_VecDesc.size(); ++i)                        // search through all VecDesc for scalar unknowns
            {
                // Error checking
                Assert(_VecDesc[i], DROPSErrCL("ParMultiGridCL::SendUnknowns: No Information about Unknowns"),
                       DebugParallelNumC);
                Assert(_VecDesc[i]->RowIdx, DROPSErrCL("ParMultiGridCL::SendUnknowns: RowIdx not set"),
                       DebugParallelNumC);
                const Uint idx=(_VecDesc[i])->RowIdx->GetIdx();             // index number

                // Check if there are the right unknowns on this simplex to this index
                if ( (_VecDesc[i])->RowIdx->GetNumUnknownsOnSimplex<SimplexT>()!=1   // only scalar unknowns are handled here
                     || !sp->Unknowns.Exist(idx) )                          // unknown must exist
                    continue;

                const IdxT sysnum= sp->Unknowns(idx);                       // systemnumber
                const double data= !sp->Unknowns.UnkRecieved(idx) ?         // value of the unknown
                                    (_VecDesc[i])->Data[ sysnum] :
                                    _RecvBuf[sysnum];
                buffer[counter].SetData(data);                              // put data into message
                buffer[counter].SetIdx(i);                                  // Tell other proc, where to put this unknown
                ++counter;                                                  // next element
            }
        }
    }
    else if (type == AddedVecCL::GetType())
    { // see above for documentation
        AddedVecCL* buffer = static_cast<AddedVecCL*>(buf);
        if (sp->Unknowns.Exist()){
            for (size_t i=0; i<_VecDesc.size(); ++i)
            {
                Assert(_VecDesc[i], DROPSErrCL("ParMultiGridCL::SendUnknowns: No Information about Unknowns"),
                       DebugParallelNumC);
                Assert(_VecDesc[i]->RowIdx, DROPSErrCL("ParMultiGridCL::SendUnknowns: RowIdx not set"),
                       DebugParallelNumC);
                const Uint idx=(_VecDesc[i])->RowIdx->GetIdx();

                if ( (_VecDesc[i])->RowIdx->GetNumUnknownsOnSimplex<SimplexT>()!=3   // vectorial unknowns
                     || !sp->Unknowns.Exist(idx) )
                    continue;

                const IdxT sysnum= sp->Unknowns(idx);
                Point3DCL data;
                if (!sp->Unknowns.UnkRecieved(idx))
                {
                    data[0]=(_VecDesc[i])->Data[sysnum+0];
                    data[1]=(_VecDesc[i])->Data[sysnum+1];
                    data[2]=(_VecDesc[i])->Data[sysnum+2];
                }
                else
                {
                    data[0]=_RecvBuf[sysnum+0];
                    data[1]=_RecvBuf[sysnum+1];
                    data[2]=_RecvBuf[sysnum+2];
                }
                buffer[counter].SetData(data);
                buffer[counter].SetIdx(i);
                ++counter;
            }
        }
    }
    Assert((int)counter==cnt,DROPSErrCL("ParMultiGridCL::SendUnknowns: to much or less datas added"), DebugParallelC);
}


/****************************************************************************
* R E C V  U N K N O W N S                                                  *
*****************************************************************************
*   If a simplex is recieved, this procedure is called in order to put the  *
*   recieved data into the recieve buffer, if the unknowns have not been    *
*   exists before this simplex is recieved.                                 *
****************************************************************************/
template<class SimplexT>
  void ParMultiGridCL::RecvUnknowns(SimplexT *sp, DDD_TYPE type, void *buf, int cnt)
{
    if (cnt==0) return;                                                     // nothing to recieved so skipp
    if (type == AddedScalCL::GetType())                                     // recieve scalar values
    {
        AddedScalCL* buffer = static_cast<AddedScalCL*>(buf);               // cast to AddedScalCL*

        // put all new recieved data into the _RecvBuf
        for (int i=0; i<cnt; ++i)                                           // go over all incomming messages
        {
            const Uint idxVecDesc = buffer[i].GetIdx(),                     // under which index from _VecDesc the datas should be stored
                  idx= (_VecDesc[idxVecDesc])->RowIdx->GetIdx();            // get the index of the unknown (e.g. pressure, levelset,...)

            if (!sp->Unknowns.Exist(idx) && !sp->Unknowns.UnkRecieved(idx)) // if this unknown is not known on this simplex
            {
                if(_RecvBufPos>=_RecvBuf.size())                            // check if recieve buffer is big enough
                    EnlargeRecieveBuffer();
                sp->Unknowns.Prepare(idx);                                  // create the UnknownIdxCL
                sp->Unknowns(idx)=_RecvBufPos;                              // remeber, where unknown has been put
                sp->Unknowns.SetUnkRecieved(idx);                           // remeber, this unknown is recieved
                _RecvBuf[_RecvBufPos]=buffer[i].GetData();                  // put data into the _RecvBuf
                _RecvBufPos+=1;
            }
        }
    }
    else if (type == AddedVecCL::GetType() )
    {   // the same as above!
        AddedVecCL* buffer = static_cast<AddedVecCL*>(buf);

        for (int i=0; i<cnt; ++i)
        {
            const Uint idxVecDesc = buffer[i].GetIdx(),
                  idx= (_VecDesc[idxVecDesc])->RowIdx->GetIdx();
            Point3DCL recv( buffer[i].GetData() );

            if (!sp->Unknowns.Exist(idx) && !sp->Unknowns.UnkRecieved(idx))
            {
                if(_RecvBufPos+3>_RecvBuf.size())
                    EnlargeRecieveBuffer();
                sp->Unknowns.Prepare(idx);
                sp->Unknowns(idx)=_RecvBufPos;
                sp->Unknowns.SetUnkRecieved(idx);
                _RecvBuf[_RecvBufPos+0]= recv[0];
                _RecvBuf[_RecvBufPos+1]= recv[1];
                _RecvBuf[_RecvBufPos+2]= recv[2];
                _RecvBufPos+=3;
            }
        }
    }
}

/****************************************************************************
* L I N E A R   I N T E R P O L A T I O N                                   *
*****************************************************************************
*   do a linear interpolation on an edge. This function returns true, iff   *
*   all needed values are available, else false.                            *
****************************************************************************/
template<typename BndT>
bool ParMultiGridCL::LinearInterpolation(const EdgeCL& e, Uint idx, const BndT* bnd, const VectorCL& data, typename BndT::bnd_type& new_dof)
{
    const VertexCL *vp0=e.GetVertex(0), *vp1=e.GetVertex(1);

    if (    ( !bnd->IsOnDirBnd(*vp0) && !vp0->Unknowns.Exist(idx) )
         || ( !bnd->IsOnDirBnd(*vp1) && !vp1->Unknowns.Exist(idx) ) )
        return false;

    typedef typename BndT::bnd_type DataT;
    DataT dof0, dof1;                                   // unknowns on vertices of edge

    // gather dof0 on vertex 0
    if (bnd->IsOnDirBnd(*vp0))                          // vertex lies on a dirichlet boundary
        dof0= bnd->GetDirBndValue(*vp0);
    else if (vp0->Unknowns.UnkRecieved(idx))            // unknowns has been recieved, so we can find them in the recieve buffer
        dof0= GetDofOutOfVector<VertexCL, BufferCT, DataT>()(*vp0, idx, _RecvBuf);
    else                                                // value cann be get out of the given source
        dof0= GetDofOutOfVector<VertexCL, VectorCL, DataT>()(*vp0, idx, data);

    // gather dof1 on vertex 1
    if (bnd->IsOnDirBnd(*vp1))                          // vertex lies on a dirichlet boundary
        dof1= bnd->GetDirBndValue(*vp1);
    else if (vp1->Unknowns.UnkRecieved(idx))            // unknowns has been recieved, so we can find them in the recieve buffer
        dof1= GetDofOutOfVector<VertexCL, BufferCT, DataT>()(*vp1, idx, _RecvBuf);
    else                                                // value cann be get out of the given source
        dof1= GetDofOutOfVector<VertexCL, VectorCL, DataT>()(*vp1, idx, data);

    new_dof=0.5*(dof0+dof1);
    return true;
}


/// \brief Get dof of a simplex (vectorial)
template<typename SimplexT>
  Point3DCL ParMultiGridCL::GetDof<SimplexT,Point3DCL>::operator() (const SimplexT& s, Uint idx, int pos)
/** Get unknowns on a simplex, if the index is known to the VecDesc.
    \param s the simplex
    \param idx index type
    \param pos position within the VecDesc. If -1, then this routine searches for right position
    \return value of the dof*/
{
    Assert(s.Unknowns.Exist() && s.Unknowns.Exist(idx),
           DROPSErrCL("ParMultiGridCL::GetDof<Point3DCL>: No Unknowns on the simplex"),
           DebugParallelNumC);
    Assert(pos<0 || _VecDesc[pos]->RowIdx->GetIdx()==idx,
           DROPSErrCL("ParMultiGridCL::GetDof<Point3DCL>: Wrong position is given"),
           DebugParallelNumC);

    Point3DCL dof;
    if (s.Unknowns.UnkRecieved(idx)){               // unknowns are stored within recieve buffer
        for (int i=0; i<3; ++i)
            dof[i]= _RecvBuf[s.Unknowns(idx)+i];
    }
    else{                                           // unknowns are stored within VecDesc
        // find index in _VecDesc
        if (pos<0){
            pos=std::numeric_limits<int>::max();
            for (int i=0; i<(int)_VecDesc.size() && pos>i; ++i)
                if (_VecDesc[i]->RowIdx->GetIdx()==idx)
                    pos=i;
        }

        Assert(pos!=std::numeric_limits<int>::max(),
               DROPSErrCL("ParMultiGridCL::GetDof<Point3DCL>: No appropriate index known"),
               DebugParallelNumC);

        for (int i=0; i<3; ++i)
            dof[i]= _VecDesc[pos]->Data[s.Unknowns(idx)+i];
    }
    return dof;
}

/// \brief Get dof of a simplex (scalar)
template<typename SimplexT>
  double ParMultiGridCL::GetDof<SimplexT,double>::operator() (const SimplexT& s, Uint idx, int pos)
/** Get unknowns on a simplex, if the index is known to the _VecDesc.
    \param s the simplex
    \param idx index type
    \param pos position within _VecDesc
    \return value of the dof*/
{
    Assert(s.Unknowns.Exist() && s.Unknowns.Exist(idx),
           DROPSErrCL("ParMultiGridCL::GetDof<double>: No Unknowns on the simplex"),
           DebugParallelNumC);

    Assert(pos<0 || _VecDesc[pos]->RowIdx->GetIdx()==idx,
           DROPSErrCL("ParMultiGridCL::GetDof<double>: Wrong position is given"),
           DebugParallelNumC);

    double dof;
    if (s.Unknowns.UnkRecieved(idx))    // unknown can be found within recieve buffer
        dof= _RecvBuf[s.Unknowns(idx)];
    else{                               // unknown can be found in VecDesc
        // find index in _VecDesc
        if (pos<0){
            pos=std::numeric_limits<int>::max();
            for (int i=0; i<(int)_VecDesc.size() && pos>i; ++i)
                if (_VecDesc[i]->RowIdx->GetIdx()==idx)
                    pos=i;
        }

        // check if found
        Assert(pos!=std::numeric_limits<int>::max(),
            DROPSErrCL("ParMultiGridCL::GetDof<double>: No appropriate index known"),
            DebugParallelNumC);

        dof= _VecDesc[pos]->Data[s.Unknowns(idx)];
    }
    return dof;
}

/// \brief Put vectorial data into a given vector
template<typename SimplexT>
  void ParMultiGridCL::SetDof(const SimplexT& s, Uint idx, VectorCL& data, const Point3DCL& dof)
/** In order to easyly put data into a given vector, this function can be used
    \param s    Vertex or edge, the unknown is living on
    \param idx  index type
    \param data vector of the unknowns
    \param dof  the value*/
{
    Assert(s.Unknowns.Exist() && s.Unknowns.Exist(idx),
           DROPSErrCL("ParMultiGridCL::SetDof<Point3DCL>: No Unknowns on the simplex"),
           DebugParallelNumC);
    data[s.Unknowns(idx)+0]=dof[0];
    data[s.Unknowns(idx)+1]=dof[1];
    data[s.Unknowns(idx)+2]=dof[2];
}

/// \brief Put scalar data into a given vector
template<typename SimplexT>
  void ParMultiGridCL::SetDof(const SimplexT& s, Uint idx, VectorCL& data, const double& dof)
/** In order to easyly put data into a given vector, this function can be used
    \param s    Vertex or edge, the unknown is living on
    \param idx  index type
    \param data vector of the unknowns
    \param dof  the value*/
{
    Assert(s.Unknowns.Exist() && s.Unknowns.Exist(idx),
           DROPSErrCL("ParMultiGridCL::SetDof<double>: No Unknowns on the simplex"),
           DebugParallelNumC);
    data[s.Unknowns(idx)]=dof;
}

/// \brief Get a vectorial unknown out of the recieve buffer
template<typename SimplexT, typename ContainerT>
  Point3DCL ParMultiGridCL::GetDofOutOfVector<SimplexT, ContainerT, Point3DCL>::operator()
          (const SimplexT& s, Uint idx, const ContainerT& data)
/** Easyly access of values stored in a vector
    \param s   simplex on which the unknown lies
    \param idx index of the unknown
    \return dof out of the given vector*/
{
    Assert(s.Unknowns.Exist() && s.Unknowns.Exist(idx),
           DROPSErrCL("ParMultiGridCL::GetDofOutOfRecvBuffer<Point3DCL>: No Unknowns on the simplex"),
           DebugParallelNumC);
    return MakePoint3D(data[s.Unknowns(idx)+0], data[s.Unknowns(idx)+1], data[s.Unknowns(idx)+2]);
}

/// \brief Get a scalar unknown out of the recieve buffer
template<typename SimplexT, typename ContainerT>
  double ParMultiGridCL::GetDofOutOfVector<SimplexT, ContainerT, double>::operator()
              (const SimplexT& s, Uint idx, const ContainerT& data)
/** Within transfer and refinement operation, some unknowns are put into a
    recieve buffer, which can be accessed by this funktion easily
    \param s   simplex on which the unknown lies
    \param idx index of the unknown
    \param dof at exit value of the unknown*/
{
    Assert(s.Unknowns.Exist() && s.Unknowns.Exist(idx),
           DROPSErrCL("ParMultiGridCL::GetDofOutOfRecvBuffer<double>: No Unknowns on the simplex"),
           DebugParallelNumC);
    return data[s.Unknowns(idx)];
}

/// \brief Put a value of vectorial unknown into the Recieve Buffer
template<typename SimplexT>
  void ParMultiGridCL::PutDofIntoRecvBuffer(SimplexT&s, Uint idx, const Point3DCL& dof)
{
    if (!s.Unknowns.Exist(idx)){
        if (_RecvBufPos+3>_RecvBuf.size())
            EnlargeRecieveBuffer();
        s.Unknowns.Prepare(idx);
        s.Unknowns(idx)=_RecvBufPos;
        _RecvBuf[_RecvBufPos+0]= dof[0];
        _RecvBuf[_RecvBufPos+1]= dof[1];
        _RecvBuf[_RecvBufPos+2]= dof[2];
        _RecvBufPos+=3;
    }
}

/// \brief Put a value of scalar unknown into the Recieve Buffer
template<typename SimplexT>
  void ParMultiGridCL::PutDofIntoRecvBuffer(SimplexT& s, Uint idx, const double& dof)
{
    if (!s.Unknowns.Exist(idx)){
        if (_RecvBufPos >= _RecvBuf.size() )
            EnlargeRecieveBuffer();
        s.Unknowns.Prepare(idx);
        s.Unknowns(idx)=_RecvBufPos;
        _RecvBuf[_RecvBufPos]= dof;
        _RecvBufPos+=1;
    }
}


/// \brief Put datas on an edae from the old vector onto a vertex into a new vector
template<typename BndT>
  void ParMultiGridCL::PutData(MultiGridCL::const_EdgeIterator& sit,
                               const VectorCL* const old_data, VectorCL* new_data,
                               const Uint old_idx, const Uint new_idx,
                               const IdxDescCL* idxDesc, const BndT* bnd)
/** This routine puts the unknowns on a simplex according to an index into
    a new data vector. Therefore datas are taken from the recieve buffer, if the
    data has been stored on another proc before the refinement and migration, or
    out of the "old" vector, if the calling proc has already owned these data
    before the refinement. <br>
    For P1 finite elements it may happen, that the function RepairAfterRefineP1
    cannot interpolate a new midvertex, because the values on the parent edge
    has been deleted. So this interpolation is done here, too.<br>
    For P2 finite elements it may happen, that unknowns on an edge has been
    moved from the edge to the new created midvertex. This is done by this
    procedure as well.
    \param sit         Iterator onto a simplex
    \param old_data    Pointer to the old data vector
    \param new_data    Pointer to the new data vector
    \param old_idx     old index
    \param new_idx     new index
    \param idxDesc     describer of index, just used for getting information
                       about number of unknowns on simplices
    \param bnd         boundary conditions for the index*/
{
    typedef typename BndT::bnd_type DataT;
    DataT new_dof;
    const Uint numUnknowns=idxDesc->NumUnknownsEdge();
    if (numUnknowns>0)
    {
        if (sit->Unknowns.Exist()
            && sit->Unknowns.Exist(new_idx)
            && sit->Unknowns.Exist(old_idx) )
        {
            if (sit->Unknowns.UnkRecieved(old_idx))                         // this is a recieved unknown
                new_dof= GetDofOutOfVector<EdgeCL, BufferCT, DataT>()(*sit, old_idx, _RecvBuf);
            else
                new_dof= GetDofOutOfVector<EdgeCL, VectorCL, DataT>()(*sit, old_idx, *old_data);
            SetDof(*sit, new_idx, *new_data, new_dof);
            sit->Unknowns.SetUnkRecieved(new_idx);
        }

        // This is the special case for edges.
        if (sit->Unknowns.Exist()
            && sit->Unknowns.Exist(old_idx)
            && !sit->Unknowns.Exist(new_idx)
            && sit->IsRefined()
            && sit->GetMidVertex()->Unknowns.Exist()
            && sit->GetMidVertex()->Unknowns.Exist(new_idx)
        )
        {
            if (sit->Unknowns.UnkRecieved(old_idx))                         // this is a recieved unknown
                new_dof= GetDofOutOfVector<EdgeCL, BufferCT, DataT>()(*sit, old_idx, _RecvBuf);
            else
                new_dof= GetDofOutOfVector<EdgeCL, VectorCL, DataT>()(*sit, old_idx, *old_data);
            SetDof(*sit->GetMidVertex(), new_idx, *new_data, new_dof);
            sit->GetMidVertex()->Unknowns.SetUnkRecieved(new_idx);          // set recieved flag as a flag, that this unknown is allready set
        // Unfortunally we cannot adapt the unknown of a midvertex, that should be deleted,
        // because the reference to the midvertex has been deleted
        }
    }
    else
    {   // Another special case for edges an P1 (yes, P1) functions. It may happen,
        // that a midvertex on an edge cannot be interpolated, because the function
        // RepairAfterRefineP1 has no access to the old values of the both vertices
        // of an edge, because the priority chances and hence the unknowns are deleted.
        // Therefore this linear interpolation is done here. This is not pretty, but it
        // works.
        if ( sit->IsRefined()
             && sit->GetMidVertex()->Unknowns.Exist(new_idx)
             && !sit->GetMidVertex()->Unknowns.UnkRecieved(new_idx)
             && ( !sit->GetVertex(0)->Unknowns.Exist(new_idx) || !sit->GetVertex(1)->Unknowns.Exist(new_idx) )
           )
        {
            LinearInterpolation(*sit, old_idx, bnd, *old_data, new_dof);
            SetDof(*sit->GetMidVertex(), new_idx, *new_data, new_dof);
            sit->GetMidVertex()->Unknowns.SetUnkRecieved(new_idx);
        }
    }
}
} // end of namespace DROPS
