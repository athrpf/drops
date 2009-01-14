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

template<typename T> T ddd_cast (DDD_OBJ p)
  { return reinterpret_cast<T>(p); }

template<typename T> DDD_OBJ ddd_cast (T* p)
  { return reinterpret_cast<DDD_OBJ>(p); }

/***************************************************************************
*   D E C L A R A C T I O N   O F   F U N C T I O N S                      *
***************************************************************************/

/// \brief Get global sum
template<typename T>
  inline T GlobalSum(T, int proc=-1);

/// \brief Get global sum over a vector
template<typename T>
  inline void GlobalSum(const T*, T*, int cnt, int proc=-1);

/// \brief Get global maximum
template<typename T>
  inline T GlobalMax(T, int proc=-1);

/// \brief Get global maximum over a vector
template<typename T>
  inline void GlobalMax(const T*, T*, int cnt, int proc=-1);

/// \brief Get global minimum
template<typename T>
  inline T GlobalMin(T, int proc=-1);

/// \brief Get global minimum over a vector
template<typename T>
  inline void GlobalMin(const T*, T*, int cnt, int proc=-1);

/// \brief Collect one number from all procs
template<typename T>
  inline void Gather(T, T*, int proc);

/// \brief Collect multiple numbers from all procs
template<typename T>
  inline void Gather(const T*, T*, int cnt, int proc);

/// \brief Check if a boolean value is true on all procs
inline bool Check(int, int proc=-1);

/// \brief Check if a boolean value is true on at least one proc
inline bool GlobalOr(int, int proc=-1);


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
    static Uint _my_rank;                       // Which Id do I have?
    static Uint _size;                          // How many are out there?
    static const CommunicatorT& Communicator_;  // communicator (=MPI_COMM_WORLD, MPI::COMM_WORLD)

  public:
    ProcCL(int*, char***);                      ///< constructor
    ~ProcCL();                                  ///< destructor
      /// \brief Wait for an input of a proc
    static void Prompt(int);
      /// \brief Check if I am Master
    static bool IamMaster() { return MyRank()==Drops_MasterC; }     // Drops_MasterC defined in utils.h
      /// \brief Get rank of master processor
    static int  Master()    { return Drops_MasterC; }
      /// \brief Get used MPI communicator
    static const CommunicatorT& GetComm() { return Communicator_; }
    
    /// \name plain MPI-Calls with C++- or C-Interface of MPI
    //@{
      /// \brief check the rank, MPI has given to the calling proc
    static inline int MyRank();
      /// \brief check how many procs are used by this program
    static inline int Size();
      /// \brief MPI-Reduce-wrapper
    template <typename T>
    static inline void Reduce(const T*, T*, int, const OperationT&, int);
      /// \brief MPI-Allreduce-wrapper
    template <typename T>
    static inline void AllReduce(const T*, T*, int, const OperationT&);
      /// \brief MPI-Gather-wrapper (both data-types are the same)
    template <typename T>
    static inline void Gather(const T*, T*, int, int);
      /// \brief MPI-Allgather-wrapper (both data-types are the same)
    template <typename T>
    static inline void AllGather(const T*, T*, int);
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
    /// \name Specialized operations on vector classes
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
    static inline void Bcast(std::valarray<T>&, int);
    //@}
    //@}

    /// \name helper functions for global operations
    //@{
      /// \brief Global operation with one argument
    template <typename T>
    static inline T GlobalOp(const T&, int, const ProcCL::OperationT&);
      /// \brief Global operation with multiple argument
    template <typename T>
    static inline void GlobalOp(const T*, T*, int, int, const ProcCL::OperationT&);
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
