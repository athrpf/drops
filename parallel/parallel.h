/// \file parallel.h
/// \brief interface to MPI and helper functions to easily use communication
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_PARALLEL_H
#define DROPS_PARALLEL_H

#include <iostream>
#include <valarray>
#include <string>
#include <numeric>
#include "misc/utils.h"
#include "parallel/distributeddatatypes.h"
#ifdef HAVE_ZOLTAN
#include <zoltan.h>
#endif


namespace DROPS
{

/***************************************************************************
*   C A S T I N G   O F   D D D - P O I N T E R S                          *
***************************************************************************/
template<typename T> T ddd_cast (OBJT p)
  { return reinterpret_cast<T>(p); }

template<typename T> OBJT ddd_cast (T* p)
  { return reinterpret_cast<OBJT>(p); }


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

class ProcInitCL; //forward declaration
/***************************************************************************
*   P R O C - C L A S S                                                    *
***************************************************************************/
/// \brief Manage several procs
class ProcCL
/** This class acts as an interface to MPI. */
{
    friend class ProcInitCL;

  public:

#ifdef _MPICXX_INTERFACE
    typedef ::MPI::Op       OperationT;         ///< type of operations
    typedef ::MPI::Status   StatusT;            ///< type of stati
    typedef ::MPI::Request  RequestT;           ///< type of requests
    typedef ::MPI::Datatype DatatypeT;          ///< type of data-types
    typedef ::MPI::Comm     CommunicatorT;      ///< type of communicator
    typedef ::MPI::Aint     AintT;              ///< type of addresses
    typedef ::MPI::User_function FunctionT;     ///< type of user defined functions
#else
    typedef MPI_Op          OperationT;         ///< type of operations
    typedef MPI_Status      StatusT;            ///< type of stati
    typedef MPI_Request     RequestT;           ///< type of requests
    typedef MPI_Datatype    DatatypeT;          ///< type of data-types
    typedef MPI_Comm        CommunicatorT;      ///< type of communicator
    typedef MPI_Aint        AintT;              ///< type of addresses
    typedef MPI_User_function FunctionT;        ///< type of user defined functions
#endif

    template<typename> struct MPI_TT;           ///< Traits to determine the corresponding MPI_Datatype
                                                /// constant for a given type.
    static const DatatypeT  NullDataType;       ///< MPI-Datatype, which is not set

  private:
    static Uint my_rank_;                       // Which Id do I have?
    static Uint size_;                          // How many are out there?
    static int  procDigits_;                    // How many digits are necessary to decode rank of process?
    static const CommunicatorT& Communicator_;  // communicator (=MPI_COMM_WORLD, MPI::COMM_WORLD)
    static MuteStdOstreamCL* mute_;             // for muting std::cout, std::cout, std::clog

    static ProcCL* instance_;                   ///< only one instance of ProcCL may exist (Singleton-Pattern)
    ProcCL(int*, char***);                      ///< constructor, mutes all non-master standard output streams
    ~ProcCL();                                  ///< destructor

  public:
      /// \brief Get a pointer to the ProcCL (Singleton-Pattern)
    static ProcCL* InstancePtr(int* argc = 0, char*** argv = 0) { return instance_ ? instance_ : (instance_= new ProcCL(argc, argv)); }
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
    /// \name Parallel output
    //@{
      // \brief Mute output of standard output streams for all non-master procs
    static void MuteStdOstreams()    { if (!IamMaster()) mute_->Mute(); }
      // \brief Recover behavior of standard output streams
    static void RecoverStdOstreams() { mute_->Recover(); }
    //@}

    /// \brief Append rank of processor to an string
    static void AppendProcNum( std::string&);

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
      /// \brief MPI-Gatherv-wrapper or MPI-Allgatherv wrapper if root<0 (both data-types are the same)
    template<typename T>
    static inline void Gatherv( const T*, int, T*, const int*, const int*, int root);
      /// \brief MPI-Gatherv-wrapper or MPI-Allgatherv wrapper if root<0 (both data-types are the same)
    template<typename T>
    static inline void Scatterv( const T*, const int*, const int*, T*, int, int root=Drops_MasterC);
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
      /// \brief MPI-wrapper for creating an operation
    static inline void InitOp(OperationT&, FunctionT*, bool);
      /// \brief MPI-wrapper for freeing an operation
    static inline void FreeOp(OperationT&);
      /// \brief MPI-Bcast-wrapper
    template <typename T>
    static inline void Bcast(T*, int, int);
      /// \brief MPI get time
    static inline double Wtime();
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
    template <typename T>
    static inline int GetMessageLength(int, int);
    static inline int GetMessageLength(int, int, DatatypeT&);
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
    static inline RequestT Isend(const std::vector<T>&, int, int);
    template <typename T>
    static inline void Recv(std::valarray<T>&, int, int);
    template <typename T>
    static inline RequestT Irecv(std::valarray<T>&, int, int);
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
      /// \brief Global operation a valarray as argument
    template <typename T>
    static inline std::valarray<T> GlobalOp(const std::valarray<T>&, int, const ProcCL::OperationT&);
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
        { ProcCL::GlobalOp(myData, allData, cnt, proc, MPI_SUM_Operation); }
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
    static void Gather(T myData, T* allData, int proc=-1)
        { Gather(&myData, allData, 1, proc); }

    template<typename T>
    static std::valarray<T> Gather(const std::valarray<T>& myData, int proc) {
        std::valarray<T> allData(Size()*myData.size());
        Gather(Addr(myData), Addr(allData), myData.size(), proc);
        return allData;
    }
    template<typename T>
    static std::vector<T> Gather(const T& myData, int proc) {
        std::vector<T> allData( Size());
        Gather( &myData, Addr(allData), 1, proc);
        return allData;
    }
    template<typename T>
    static std::valarray<T> Gatherv( const std::valarray<T>& myData, int proc=-1)
    {
        const std::vector<int> recvc= ProcCL::Gather( (int)myData.size(), -1);      // number of elements per proc
        std::valarray<T> allData( std::accumulate( recvc.begin(), recvc.end(), 0)); // receive buffer
        std::vector<int> displ( recvc.size());                                    // displacements
        displ[0]=0; for ( size_t i=1; i<displ.size(); ++i) displ[i]= displ[i-1]+recvc[i-1];
        Gatherv( Addr(myData), (int)myData.size(), Addr(allData), Addr( recvc), Addr( displ), proc);
        return allData;
    }
    template<typename T>
    static std::valarray<T> Gatherv( const std::vector<T>& myData, int proc=-1)
    {
        const std::vector<int> recvc= ProcCL::Gather( (int)myData.size(), -1);      // number of elements per proc
        std::valarray<T> allData( std::accumulate( recvc.begin(), recvc.end(), 0)); // receive buffer
        std::vector<int> displ( recvc.size());                                    // displacements
        displ[0]=0; for ( size_t i=1; i<displ.size(); ++i) displ[i]= displ[i-1]+recvc[i-1];
        Gatherv( Addr(myData), (int)myData.size(), Addr(allData), Addr( recvc), Addr( displ), proc);
        return allData;
    }
    //@}
};

/// \brief Manage construction and destruction of ProcCL
class ProcInitCL
{
  public:
    ProcInitCL(int* argc, char*** argv) {
        ProcCL::Instance( argc, argv);
#ifdef HAVE_ZOLTAN
        float ver;
        if (Zoltan_Initialize( *argc, *argv, &ver)!=ZOLTAN_OK)
            throw DROPSErrCL("Cannot initialize Zoltan");
#endif
    }
    ~ProcInitCL() { if (ProcCL::InstancePtr()) delete ProcCL::InstancePtr(); }
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

template<> struct ProcCL::MPI_TT<short int>
  { static const ProcCL::DatatypeT& dtype; };

#ifdef DROPS_WIN
# ifdef WIN64
template<> struct ProcCL::MPI_TT<size_t>
  { static const ProcCL::DatatypeT& dtype; };
# endif
#endif



} // namespace DROPS

#include "parallel/parallel.tpp"        // for inline and/or template functions!
#endif // _DROPS_PARALLEL_H
