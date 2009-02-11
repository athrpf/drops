/***************************************************************************
*  File:    parallel.h                                                     *
*  Content: Interface for parallel support                                 *
*           ProcCL - Management of the procs                               *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           10.11.2005                                             *
*  last modified:   10.11.2005                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file parallel.h
/// \brief Management of processors and interface for MPI-Routines

#ifndef _DROPS_PARALLEL_H
#define _DROPS_PARALLEL_H

#include <iostream>
#include <valarray>
#include "misc/utils.h"
#include <ddd.h>

namespace DROPS
{

/***************************************************************************
*   C A S T I N G   O F   D D D - P O I N T E R S                          *
***************************************************************************/
template<typename T> T ddd_cast (DDD_OBJ p)
  { return reinterpret_cast<T>(p); }

template<typename T> DDD_OBJ ddd_cast (T* p)
  { return reinterpret_cast<DDD_OBJ>(p); }


/***************************************************************************
*   G L O B A L - O P E R A T I O N S                                      *
***************************************************************************/

/// \name MPI-Operations
//@{
#ifdef _MPICXX_INTERFACE
#  define MPI_MAX_Operation  MPI::MAX
#  define MPI_MIN_Operation  MPI::MIN
#  define MPI_SUM_Operation  MPI::SUM
#  define MPI_LAND_Operation MPI::LAND
#  define MPI_LOR_Operation  MPI::LOR
#else
#  define MPI_MAX_Operation  MPI_MAX
#  define MPI_MIN_Operation  MPI_MIN
#  define MPI_SUM_Operation  MPI_SUM
#  define MPI_LAND_Operation MPI_LAND
#  define MPI_LOR_Operation  MPI_LOR
#endif
//@}

/***************************************************************************
*   P R O C - C L A S S                                                    *
***************************************************************************/
/// \brief Manage several procs
class ProcCL
/** This class acts as an interface to MPI. */
{
  public:

#ifdef _MPICXX_INTERFACE
    typedef ::MPI::Op       OperationT;         ///< type of operations
    typedef ::MPI::Status   StatusT;            ///< type of stati
    typedef ::MPI::Request  RequestT;           ///< type of requests
    typedef ::MPI::Datatype DatatypeT;          ///< type of data-types
    typedef ::MPI::Comm     CommunicatorT;      ///< type of communicator
    typedef ::MPI::Aint     AintT;              ///< type of addresses
#else
    typedef MPI_Op          OperationT;         ///< type of operations
    typedef MPI_Status      StatusT;            ///< type of stati
    typedef MPI_Request     RequestT;           ///< type of requests
    typedef MPI_Datatype    DatatypeT;          ///< type of data-types
    typedef MPI_Comm        CommunicatorT;      ///< type of communicator
    typedef MPI_Aint        AintT;              ///< type of addresses
#endif

    template<typename> struct MPI_TT;           ///< Traits to determine the corresponding MPI_Datatype
                                                /// constant for a given type.
    static const DatatypeT  NullDataType;       ///< MPI-Datatype, which is not set

  private:
    static Uint my_rank_;                       // Which Id do I have?
    static Uint size_;                          // How many are out there?
    static const CommunicatorT& Communicator_;  // communicator (=MPI_COMM_WORLD, MPI::COMM_WORLD)

    /// \name helper functions for global operations
    //@{
      /// \brief Global operation with one argument
    template <typename T>
    static inline T GlobalOp(const T&, int, const ProcCL::OperationT&);
      /// \brief Global operation with multiple argument
    template <typename T>
    static inline void GlobalOp(const T*, T*, int, int, const ProcCL::OperationT&);
      /// \brief Global operation a valarray as argument
    template <typename T>
    static inline std::valarray<T> GlobalOp(const std::valarray<T>&, int, const ProcCL::OperationT&);
    //@}

    static ProcCL* instance_;                   ///< only one instance of ProcCL may exist (Singleton-Pattern)
    ProcCL(int*, char***);                      ///< constructor
    ~ProcCL();                                  ///< destructor

  public:
      /// \brief Get a pointer to the ProcCL (Singleton-Pattern)
    static ProcCL* InstancePtr(int* argc, char*** argv) { return instance_ ? instance_ : (instance_= new ProcCL(argc, argv)); }
      /// \brief Get a reference to the ProcCL (Singleton-Pattern)
    static ProcCL& Instance(int* argc, char*** argv)    { return *InstancePtr(argc, argv); }
      /// \brief Wait for an input of a proc
    static void Prompt(int);
      /// \brief Check if I am Master
    static bool IamMaster() { return my_rank_==Drops_MasterC; }     // Drops_MasterC defined in utils.h
      /// \brief Get rank of master processor
    static int  Master()    { return Drops_MasterC; }
      /// \brief Get used MPI communicator
    static const CommunicatorT& GetComm() { return Communicator_; }
      /// \brief check the rank, MPI has given to the calling proc
    static int MyRank()     { return my_rank_; }
      /// \brief check how many procs are used by this program
    static int Size()       { return size_; }

    /// \name plain MPI-Calls with C++- or C-Interface of MPI
    //@{
      /// \brief MPI-Reduce-wrapper
    template <typename T>
    static inline void Reduce(const T*, T*, int, const OperationT&, int);
      /// \brief MPI-Allreduce-wrapper
    template <typename T>
    static inline void AllReduce(const T*, T*, int, const OperationT&);
      /// \brief MPI-Gather-wrapper or MPI-Allgather-wrapper if root<0 (both data-types are the same)
    template <typename T>
    static inline void Gather(const T*, T*, int, int root);
      /// \brief MPI-Probe-wrapper
    static inline void Probe(int, int, StatusT&);
      /// \brief MPI-Get_count-wrapper
    static inline int GetCount(StatusT&, const DatatypeT&);
      /// \brief MPI-Wait-wrapper
    static inline void Wait(RequestT&);
      /// \brief MPI-Waitall-wrapper
    static inline void WaitAll(int, RequestT*);
      /// \brief MPI-Recv-wrapper
    template <typename T>
    static inline void Recv(T*, int, const DatatypeT&, int, int);
      /// \brief MPI-Irecv-wrapper
    template <typename T>
    static inline RequestT Irecv(T*, int, int, int);
      /// \brief MPI-Send-wrapper with given datatype
    template <typename T>
    static inline void Send(const T*, int, const DatatypeT&, int, int);
      /// \brief MPI-Isend-wrapper with given datatype
    template <typename T>
    static inline RequestT Isend(const T*, int, const DatatypeT&, int, int);
      /// \brief MPI-Get_address-wrapper
    template <typename T>
    static inline AintT Get_address(T*);
      /// \brief MPI-wrapper for creating an indexed datatype
    template <typename T>
    static inline DatatypeT CreateIndexed(int count, const int*, const int*);
      /// \brief MPI-wrapper for creating a structured datatype
    static inline DatatypeT CreateStruct(int, const int*, const AintT*, const DatatypeT*);
      /// \brief MPI-Wrapper for commiting a datatype
    static inline void Commit(DatatypeT&);
      /// \brief MPI-wrapper for freeing a datatype
    static inline void Free(ProcCL::DatatypeT& type);
      /// \brief MPI-Bcast-wrapper
    template <typename T>
    static inline void Bcast(T*, int, int);
      /// \brief MPI-Barrier-wrapper
    static inline void Barrier();
      /// \brief Abort MPI
    static inline void Abort(int);
    //@}

    /// \name Wrapper for MPI commands (no direct call of MPI-functions)
    //@{
      /// \brief MPI-Get_count-wrapper
    template <typename T>
    static inline int GetCount(StatusT&);
      /// \brief MPI-Waitall-wrapper for all requests in a valarray
    static inline void WaitAll(std::valarray<RequestT>&);
      /// \brief MPI-Waitall-wrapper for all requests in a vector
    static inline void WaitAll(std::vector<RequestT>&);
      /// \brief MPI-Recv-wrapper
    template <typename T>
    static inline void Recv(T*, int, int, int);
      /// \brief MPI-Send-wrapper with automatic generated datatype
    template <typename T>
    static inline void Send(const T*, int, int, int);
      /// \brief MPI-Isend-wrapper with automatic generated datatype
    template <typename T>
    static inline RequestT Isend(const T*, int, int, int);
    /// \name Specialized operations on vector-type classes
    //@{
    template <typename T>
    static inline RequestT Isend(const std::valarray<T>&, int, int);
    template <typename T>
    static inline RequestT Isend(const std::valarray<T>&, const DatatypeT&, int, int);
    template <typename T>
    static inline RequestT Isend(const std::vector<T>&, const DatatypeT&, int, int);
    template <typename T>
    static inline void Recv(std::valarray<T>&, int, int);
    template <typename T>
    static inline RequestT Irecv(std::valarray<T>&, int, int);
    template <typename T>
    static inline void Bcast(std::valarray<T>&, int);
    //@}
    //@}

    /// \name Global operations with synchronization
    //@{
    /// \name Global sum
    //@{
    template<typename T>
    static T GlobalSum(const T& myData, int proc=-1)
        { return ProcCL::GlobalOp(myData, proc, MPI_SUM_Operation); }
    template<typename T>
    static void GlobalSum(const T* myData, T* allData, int cnt, int proc=-1)
        { return ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_SUM_Operation); }
    template<typename T>
    static  std::valarray<T> GlobalSum(const std::valarray<T>& myData, int proc=-1)
        { return ProcCL::GlobalOp(myData, proc, MPI_SUM_Operation); }
    //@}
    /// \name Global maximum
    //@{
    template<typename T>
    static T GlobalMax(const T& myData, int proc=-1)
        { return ProcCL::GlobalOp(myData, proc, MPI_MAX_Operation); }
    template<typename T>
    static void GlobalMax(const T* myData, T* allData, int cnt, int proc=-1)
        { ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_MAX_Operation); }
    template<typename T>
    static std::valarray<T> GlobalMax(const std::valarray<T>& myData, int proc=-1)
        { return ProcCL::GlobalOp(myData, proc, MPI_MAX_Operation); }
    //@}
    /// \name Global minimum
    //@{
    template<typename T>
    static T GlobalMin(const T& myData, int proc=-1)
        { return ProcCL::GlobalOp(myData, proc, MPI_MIN_Operation); }
    template<typename T>
    static void GlobalMin(const T* myData, T* allData, int cnt, int proc=-1)
        { return ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_MIN_Operation); }
    template<typename T>
    static  std::valarray<T> GlobalMin(const std::valarray<T>& myData, int proc=-1)
        { return ProcCL::GlobalOp(myData, proc, MPI_MIN_Operation); }
    //@}
    /// \brief Check if a boolean value is true on at least one proc
    static bool GlobalOr(bool myVal, int proc=-1)
        { return (bool)ProcCL::GlobalOp((int)myVal, proc, MPI_LOR_Operation); }
    /// \brief Check if a boolean value is true on all procs
    static bool Check(bool myVal, int proc=-1)
        { return (bool)ProcCL::GlobalOp((int)myVal, proc, MPI_LAND_Operation); }
    //@}

    /// \name Gather operations
    //@{
    template<typename T>
    static void Gather(T myData, T* allData, int proc)
        { Gather(&myData, allData, 1, proc); }

    template<typename T>
    static std::valarray<T> Gather(const std::valarray<T>& myData, int proc) {
        std::valarray<T> allData(Size()*myData.size());
        Gather(Addr(myData), Addr(allData), myData.size(), proc);;
        return allData;
    }
    //@}

};

template<> struct ProcCL::MPI_TT<int>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<Uint>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<Ulint>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<Usint>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<double>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<char>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<byte>
  { static const ProcCL::DatatypeT& dtype; };

template<> struct ProcCL::MPI_TT<float>
  { static const ProcCL::DatatypeT& dtype; };

} // namespace DROPS

#include "parallel/parallel.tpp"        // for inline and/or template functions!
#endif // _DROPS_PARALLEL_H
