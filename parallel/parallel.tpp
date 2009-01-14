/***************************************************************************
*  File:    parallel.tpp                                                   *
*  Content: Interface for parallel support                                 *
*           ProcCL - Management of the procs                               *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           10.11.2005                                             *
*  last modified:   10.11.2005                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file parallel.tpp

namespace DROPS
{

#ifdef _MPICXX_INTERFACE

/// \name MPI-Operations
//@{
# define MPI_MAX_Operation  MPI::MAX
# define MPI_MIN_Operation  MPI::MIN
# define MPI_SUM_Operation  MPI::SUM
# define MPI_LAND_Operation MPI::LAND
# define MPI_LOR_Operation  MPI::LOR
//@}

/// \name MPI-Calls with C++-Interface of MPI
//@{
inline int ProcCL::MyRank()
  { return Communicator_.Get_rank();}

inline int ProcCL::Size()
  { return Communicator_.Get_size(); }

template <typename T>
  inline void ProcCL::Reduce(const T* myData, T* globalData, int size, const OperationT& op, int root)
  { Communicator_.Reduce(myData, globalData, size, ProcCL::MPI_TT<T>::dtype, op, root); }

template <typename T>
  inline void ProcCL::AllReduce(const T* myData, T* globalData, int size, const ProcCL::OperationT& op)
  { Communicator_.Allreduce(myData, globalData, size, ProcCL::MPI_TT<T>::dtype, op); }

template <typename T>
  inline void ProcCL::Gather(const T* myData, T* globalData, int size, int root)
  { Communicator_.Gather(myData, size, ProcCL::MPI_TT<T>::dtype, globalData, size, ProcCL::MPI_TT<T>::dtype, root); }

template <typename T>
  inline void ProcCL::AllGather(const T* myData, T* globalData, int size)
  { Communicator_.Allgather(myData, size, ProcCL::MPI_TT<T>::dtype, globalData, size, ProcCL::MPI_TT<T>::dtype); }

inline void ProcCL::Probe(int source, int tag, ProcCL::StatusT& status)
  { Communicator_.Probe(source, tag, status); }

inline int ProcCL::GetCount(ProcCL::StatusT& status, const ProcCL::DatatypeT& type)
  { return status.Get_count(type); }

inline void ProcCL::Wait(RequestT& req)
  { req.Wait(); }

inline void ProcCL::WaitAll(int count, RequestT* req)
  { RequestT::Waitall(count, req); }

template <typename T>
  inline void ProcCL::Recv(T* data, int count, const ProcCL::DatatypeT& type, int source, int tag)
  { Communicator_.Recv(data, count, type, source, tag); }

template <typename T>
  inline ProcCL::RequestT ProcCL::Irecv(T* data, int count, int source, int tag)
  { return Communicator_.Irecv(data, count, ProcCL::MPI_TT<T>::dtype, source, tag); }

template <typename T>
  inline void ProcCL::Send(const T* data, int count, const DatatypeT& type, int dest, int tag)
  { Communicator_.Send(data, count, type, dest, tag); }

template <typename T>
  inline ProcCL::RequestT ProcCL::Isend(const T* data, int count, const ProcCL::DatatypeT& datatype, int dest, int tag)
  { return Communicator_.Isend(data, count, datatype, dest, tag); }

template <typename T>
  inline ProcCL::AintT ProcCL::Get_address(T* data)
  { return MPI::Get_address(data); }

template <typename T>
  inline ProcCL::DatatypeT ProcCL::CreateIndexed(int count, const int array_of_blocklengths[], const int array_of_displacements[])
  { return MPI_TT<T>::dtype.Create_indexed(count, array_of_blocklengths, array_of_displacements); }

inline ProcCL::DatatypeT ProcCL::CreateStruct(int count, const int* block, const ProcCL::AintT* distplace, const ProcCL::DatatypeT* tpye)
  { return MPI::Datatype::Create_struct(count, block, distplace, tpye); }

inline void ProcCL::Commit(DatatypeT& type)
  { type.Commit(); }

inline void ProcCL::Free(ProcCL::DatatypeT& type)
  { type.Free(); }


template<typename T>
  inline void ProcCL::Bcast(T* data, int size, int proc)
  { Communicator_.Bcast(data, size, ProcCL::MPI_TT<T>::dtype, proc); }

inline void ProcCL::Barrier()
  { Communicator_.Barrier(); }

inline void ProcCL::Abort(int code)
  { Communicator_.Abort(code); }
//@}

#else   // _MPICXX_INTERFACE

/// \name MPI-Operations
//@{
# define MPI_MAX_Operation  MPI_MAX
# define MPI_MIN_Operation  MPI_MIN
# define MPI_SUM_Operation  MPI_SUM
# define MPI_LAND_Operation MPI_LAND
# define MPI_LOR_Operation  MPI_LOR
//@}

/// \name MPI-Calls with C-Interface of MPI
//@{
inline int ProcCL::MyRank(){
    int rank;
    MPI_Comm_rank(Communicator_,&rank);
    return rank;
}

inline int ProcCL::Size(){
    int size;
    MPI_Comm_size(Communicator_, &size);
    return size;
}

template <typename T>
  inline void ProcCL::Reduce(const T* myData, T* globalData, int size, const OperationT& op, int root)
  { MPI_Reduce(const_cast<T*>(myData), globalData, size, ProcCL::MPI_TT<T>::dtype, op, root, Communicator_); }

template <typename T>
  inline void ProcCL::AllReduce(const T* myData, T* globalData, int size, const ProcCL::OperationT& op)
  { MPI_Allreduce(const_cast<T*>(myData), globalData, size, ProcCL::MPI_TT<T>::dtype, op, Communicator_); }

template <typename T>
  inline void ProcCL::Gather(const T* myData, T* globalData, int size, int root)
  { MPI_Gather(const_cast<T*>(myData), size, ProcCL::MPI_TT<T>::dtype, globalData, size, ProcCL::MPI_TT<T>::dtype, root,Communicator_); }

template <typename T>
  inline void ProcCL::AllGather(const T* myData, T* globalData, int size)
  { MPI_Allgather (const_cast<T*>(myData), size, ProcCL::MPI_TT<T>::dtype, globalData, size, ProcCL::MPI_TT<T>::dtype, Communicator_); }

inline void ProcCL::Probe(int source, int tag, ProcCL::StatusT& status)
  { MPI_Probe(source, tag, Communicator_, &status);}

inline int ProcCL::GetCount(ProcCL::StatusT& status, const ProcCL::DatatypeT& type){
    int count;
    MPI_Get_count(&status, type, &count);
    return count;
}

inline void ProcCL::Wait(RequestT& req){
    StatusT tmpStat;
    MPI_Wait(&req, &tmpStat);
}

inline void ProcCL::WaitAll(int count, RequestT* req){
    std::valarray<StatusT> tmpStat(count);
    MPI_Waitall(count, req, Addr(tmpStat));
}

template <typename T>
inline void ProcCL::Recv(T* data, int count, const ProcCL::DatatypeT& type, int source, int tag){
    StatusT tmpStat;
    MPI_Recv(data, count, type, source, tag, Communicator_, &tmpStat);
}

template <typename T>
  inline ProcCL::RequestT ProcCL::Irecv(T* data, int count, int source, int tag){
    RequestT req;
    MPI_Irecv(data, count, ProcCL::MPI_TT<T>::dtype, source, tag, Communicator_, &req);
    return req;
}

template <typename T>
  inline void ProcCL::Send(const T* data, int count, const DatatypeT& type, int dest, int tag)
  { MPI_Send(const_cast<T*>(data), count, type, dest, tag, Communicator_); }

template <typename T>
  inline ProcCL::RequestT ProcCL::Isend(const T* data, int count, const ProcCL::DatatypeT& datatype, int dest, int tag){
    RequestT req;
    MPI_Isend(const_cast<T*>(data), count, datatype, dest, tag, Communicator_, &req);
    return req;
}

template <typename T>
  inline ProcCL::AintT ProcCL::Get_address(T* data){
    AintT addr;
    MPI_Address(data, &addr);
    return addr;
}

template <typename T>
  inline ProcCL::DatatypeT ProcCL::CreateIndexed(int count, const int array_of_blocklengths[], const int array_of_displacements[]){
    DatatypeT newtype;
    MPI_Type_indexed(count, const_cast<int*>(array_of_blocklengths), const_cast<int*>(array_of_displacements), MPI_TT<T>::dtype, &newtype);
    return newtype;
}

inline ProcCL::DatatypeT ProcCL::CreateStruct(int count, const int* block, const ProcCL::AintT* distplace, const ProcCL::DatatypeT* type){
    DatatypeT newtype;
    MPI_Type_struct(count, const_cast<int*>(block), const_cast<AintT*>(distplace), const_cast<DatatypeT*>(type), &newtype);
    return newtype;
}

inline void ProcCL::Commit(DatatypeT& type)
  { MPI_Type_commit (&type); }

inline void ProcCL::Free(ProcCL::DatatypeT& type)
  { MPI_Type_free(&type); }

template<typename T>
  inline void ProcCL::Bcast(T* data, int size, int proc)
  { MPI_Bcast(data, size, ProcCL::MPI_TT<T>::dtype, proc, Communicator_); }

inline void ProcCL::Barrier()
  { MPI_Barrier(Communicator_); }

inline void ProcCL::Abort(int code)
  { MPI_Abort(Communicator_, code); }
//@}

#endif  // _MPICXX_INTERFACE

template <typename T>
 inline int ProcCL::GetCount(ProcCL::StatusT& status)
 { return GetCount(status, ProcCL::MPI_TT<T>::dtype); }

inline void ProcCL::WaitAll(std::valarray<ProcCL::RequestT>& reqs)
  { WaitAll((int)reqs.size(), Addr(reqs)); }

inline void ProcCL::WaitAll(std::vector<ProcCL::RequestT>& reqs)
  { WaitAll((int)reqs.size(), Addr(reqs)); }

template <typename T>
  inline void ProcCL::Recv(T* data, int count, int source, int tag)
  /// as MPI-datatype "MPI_TT<T>::dtype" is used
  { Recv(data, count, ProcCL::MPI_TT<T>::dtype, source, tag); }

template <typename T>
  inline void ProcCL::Send(const T* data, int count, int dest, int tag)
  /// as MPI-datatype "MPI_TT<T>::dtype" is used
  { ProcCL::Send(data, count, ProcCL::MPI_TT<T>::dtype, dest, tag); }

template <typename T>
  inline ProcCL::RequestT ProcCL::Isend(const T* data, int count, int dest, int tag)
  /// as MPI-datatype "MPI_TT<T>::dtype" is used
  { return Isend(data, count, ProcCL::MPI_TT<T>::dtype, dest, tag);}

template <typename T>
  inline ProcCL::RequestT ProcCL::Isend(const std::valarray<T>& data, int dest, int tag)
  /// as MPI-datatype "MPI_TT<T>::dtype" is used
  { return Isend(Addr(data), data.size(), MPI_TT<T>::dtype, dest, tag); }

template <typename T>
  inline ProcCL::RequestT ProcCL::Isend(const std::vector<T>& data, const DatatypeT& type, int dest, int tag)
  { return Isend(Addr(data), data.size(), type, dest, tag); }

template <typename T>
inline void ProcCL::Recv(std::valarray<T>& data, int source, int tag)
  { Recv(Addr(data), data.size(), source, tag); }

template <typename T>
  inline void ProcCL::Bcast(std::valarray<T>& vals, int root)
  { Bcast(Addr(vals), vals.size(), root); }


/// Perform a reduction operation over all processors with MPI.
/// \param[in] myData local copy of the variable
/// \param[in] proc   number of processor, that gets the result (if proc<0 then
///                   the result is distributed to all processors)
/// \param[in] op     operation (e.g. MPI::SUM, MPI::MAX, ...)
/// \return           result of the reduction
template<typename T>
  inline T ProcCL::GlobalOp(const T& myData, int proc, const ProcCL::OperationT& op)
{
    Assert(proc<ProcCL::Size(), DROPSErrCL("GlobalOp: proc does not exists"), DebugParallelC);
    T res=T();
    if (proc<0)
        AllReduce(&myData, &res, 1, op);
    else
        Reduce(&myData, &res, 1, op, proc);
    return res;
}

template<typename T> inline T GlobalMax(T myData, int proc)
  { return ProcCL::GlobalOp(myData, proc, MPI_MAX_Operation); }

template<typename T> inline T GlobalMin(T myData, int proc)
  { return ProcCL::GlobalOp(myData, proc, MPI_MIN_Operation); }

template<typename T> inline T GlobalSum(T myData, int proc)
  { return ProcCL::GlobalOp(myData, proc, MPI_SUM_Operation); }

inline bool Check(int myVal, int proc)
  { return ProcCL::GlobalOp(myVal, proc, MPI_LAND_Operation); }

inline bool GlobalOr(int myVal, int proc)
  { return ProcCL::GlobalOp(myVal, proc, MPI_LOR_Operation); }


/// Perform a reduction operation over all processors with a vector of values
/// with MPI.
/// \param[in]  myData  pointer to local copy of the variables
/// \param[out] allData pointer, where to store the result
/// \param[in]  cnt     number of elements in the array
/// \param[in]  proc    number of processor, that gets the result (if proc<0
///                     then the result is distributed to all processors)
/// \param[in]  op      operation (e.g. MPI::SUM, MPI::MAX, ...)
/// \pre memory for allData must be allocated before calling this function
template<typename T>
  inline void ProcCL::GlobalOp(const T* myData, T* allData, int cnt, int proc, const ProcCL::OperationT& op)
{
    Assert(proc<ProcCL::Size(), DROPSErrCL("GlobalOp: proc does not exists"), DebugParallelC);
    if (proc<0)
        ProcCL::AllReduce(myData, allData, cnt, op);
    else
        ProcCL::Reduce(myData, allData, cnt, op, proc);
}

template<typename T>
  inline void GlobalMax(const T* myData, T* allData, int cnt, int proc)
  { return ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_MAX_Operation); }

template<typename T>
  inline void GlobalMin(const T* myData, T* allData, int cnt, int proc)
  { return ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_MIN_Operation); }

template<typename T>
  inline void GlobalSum(const T* myData, T* allData, int cnt, int proc)
  { return ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_SUM_Operation); }

template<typename T>
  inline void Gather(T myData, T* allData, int proc)
/** Gather a single datum to a single or to all processor(s)
    \param myData  local copy of data
    \param allData result of gather operation (must be of size ProcCL::Size())
    \param proc    rank of proc, which recievs the result (<0 -> allgather)
*/
{
    Assert(proc<ProcCL::Size(), DROPSErrCL("Gather: proc does not exists"), DebugParallelC);
    if (proc<0)
        ProcCL::AllGather(&myData, allData, 1);
    else
        ProcCL::Gather(&myData, allData, 1, proc);
}

template<typename T>
  inline void Gather(const T* myData, T* allData, int cnt, int proc)
/** Gather data to a single or to all processor(s)
    \param myData  local copy of data
    \param allData result of gather operation (must be of size ProcCL::Size()*cnt)
    \param proc    rank of proc, which recievs the result (<0 -> allgather)
*/
{
    Assert(proc<ProcCL::Size(), DROPSErrCL("Gather: proc does not exists"), DebugParallelC);
    if (proc<0)
        ProcCL::AllGather(myData, allData, cnt);
    else
        ProcCL::Gather(myData, allData, cnt, proc);
}

} // namespace DROPS
