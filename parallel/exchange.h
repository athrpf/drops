/// \file exchange.h
/// \brief handling of a parallel distributed vectors and distributed matrices
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

/// These classes do not use the DDD-Interfaces. After the lists
/// are created no geometric datas are needed to do the
/// accumulation in opposite to the DDD-Interface. And this class
/// split the send and the recieve, so other work can be done
/// between these commands.

#ifndef DROPS_EXCHANGE_H
#define DROPS_EXCHANGE_H

#include "parallel/parallel.h"
#include <list>
#include <vector>
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "misc/problem.h"

namespace DROPS{

// fwd declaration
class ExchangeCL;
class ExchangeBlockCL;
class AccumulateMatrixCL;

/****************************************************************************
* E X C H A N G E  D A T A  S E N D  C L A S S                              *
****************************************************************************/
/// \brief Handle sending of numerical data among two processes
/** This class handles sending data among this processor and another
    processor.
    It uses a MPI Datatype for gathering.
 */
/****************************************************************************
* E X C H A N G E  D A T A  S E N D  C L A S S                              *
****************************************************************************/
class ExchangeDataSendCL
{
  public:
    friend class ExchangeMatrixCL;

  protected:
    int                     toProc_;                                            // handle communication among "me" and "toProc_"
    ProcCL::DatatypeT       SendType_;                                          // type for sending
    int                     count_;                                             // number of elements to send
#ifdef DebugParallelNumC
    Ulint                   SendTypeSize_;                                      // how big the transfered vector must be at least
#endif

    // Set processor, that receive data
    void SetToProc(const int proc) { toProc_= proc; }
    // Create MPI-Datatype for sending
    void CreateDataType(const int count, const int blocklength[], const int array_of_displacements[]);

  public:
    ExchangeDataSendCL(int proc);
    ExchangeDataSendCL();
    ExchangeDataSendCL(const ExchangeDataSendCL&);
    ~ExchangeDataSendCL();

    /// \brief Get rank of neighbor processor
    inline int GetProc() const { return toProc_; }
    // Send data to "toProc_" (nonblocking, asynchronous)
    template <typename VectorT>
    inline ProcCL::RequestT Isend(const VectorT&, int tag, Ulint offset) const;
};

/****************************************************************************
* E X C H A N G E  D A T A  C L A S S                                       *
****************************************************************************/
/// \brief Handle exchange of numerical data with one proc of one index
/** This class handels the exchange of data between this proc and another
    proc for one IdxDescCL (like pressure, levelset or velocity).
    It uses a MPI Datatype for gathering and a list of sysnums for
    scattering.
 */
/****************************************************************************
* E X C H A N G E  D A T A  C L A S S                                       *
****************************************************************************/
class ExchangeDataCL : public ExchangeDataSendCL
{
  public:
    friend class ExchangeCL;
    typedef ExchangeDataSendCL base;                 ///< base class
    typedef std::vector<int> SysnumListCT;           ///< Sequence, where to store the recieved unknowns

  private:
    SysnumListCT            Sysnums_;                ///< sysnums of the received data
    void CreateSysnums(const SysnumListCT&);         // Create SysnumListCT for receiving

  public:
    ExchangeDataCL(int proc);
    ExchangeDataCL();
    ExchangeDataCL(const ExchangeDataCL&);
    ~ExchangeDataCL();

    // Get number of received elements
    inline size_t GetNumRecvEntries() const;
    // Receive data (nonblocking)
    inline ProcCL::RequestT Irecv(int tag, VectorCL& recvBuf, Ulint offset) const;
    // add data from "toProc_"
    inline void Accumulate(VectorCL&, Ulint offsetV, VectorCL& recvBuf, Ulint offsetRecv) const;

    // print, where to store received unknowns
    void DebugInfo(std::ostream&) const;
};

/****************************************************************************
* E X C H A N G E  C L A S S                                                *
****************************************************************************/
/// \brief Handle exchange of all numerical data (for every index one class is needed!!!)
/** This class is the main class for handling the exchange of numerical datas.<br>
    On the one hand the class can accumulate a vector. Because this can be done
    in two steps (sending data and recieving datas) between neighboring procs,
    an inner product of a distributed and a accumulated vector can be performed
    effective. Between sending and recieving the numerical data, the local
    inner product can be done. <br>
    On the other hand the class can give information about all sysnums. So the
    class can decide whether a sysnum is just local on one proc or the sysnum can
    be found on different procs (and how many). It can also compute the external
    sysnum on a neighbor proc of a local sysnum. <br>
    \todo Handle Indices of various levels
*/
/****************************************************************************
* E X C H A N G E  C L A S S                                                *
****************************************************************************/
class ExchangeCL
{
    friend class IdxDescCL;

  public:
    typedef VectorBaseCL<Ulint>            IndexT;              ///< Type for storage for local and distributed sysnums
    typedef VectorBaseCL<ProcCL::RequestT> RequestCT;           ///< Type for storage Request for all neighbor procs
    typedef int ProcNumT;                                       ///< Type for number of procs
    typedef std::list<ProcNumT>            ProcNumCT;           ///< List of procs
    typedef std::vector<IdxT>              IdxVecT;             ///< Vector of indices

    IndexT  LocalIndex;                                         ///< Indices of local sysnums
    IndexT  DistrIndex;                                         ///< Indices of distributed sysnums
    IdxVecT AccDistIndex;                                       ///< Indices of distributed sysnums, this proc is responsible for

  private:
      // types for internal handling of exchanging numerical data and mapping of sysnums of proc-boundary
    typedef std::list< ExchangeDataCL >          CommListCT;        // Store information about one index
    typedef std::map<IdxT, IdxT>                 ProcMappingIdxCT;  // all local idx to external idx according to one proc
    typedef std::map<ProcNumT, ProcMappingIdxCT> MappingIdxCT;      // all idx mappings for all neighbor procs
    typedef std::vector<ProcNumCT>               SysnumProcCT;      // procs of a sysnum

      // internal handling of exchanging numerical data and mapping of sysnums of proc-boundary
    CommListCT        ExList_;          // Storage for all ExchangeData-Classes
    MappingIdxCT      MappingIdx_;      // maps local idx to external idx according to proc
    SysnumProcCT      SysProc_;         // procs that owns a sysnum
    mutable RequestCT SendRecvReq_;     // standard request handle for non-blocking sending and receiving

      // types for creating the ExchangeCL
    typedef std::vector<IdxT>                                 SendListSingleProcT;
    typedef std::map<ProcNumT, SendListSingleProcT>           SendList2ProcT;
    typedef SendList2ProcT::const_iterator                    const_SendList2ProcIter;
    typedef SendList2ProcT::iterator                          SendList2ProcIter;
    typedef std::map<ProcNumT, ExchangeDataCL::SysnumListCT>  RecvSysnumCT;
    typedef std::pair<ProcNumT, ExchangeDataCL::SysnumListCT> RecvSysnumElemT;

      // members for creating the ExchangeCL (static for DDD)
    static SendList2ProcT SendList_;
    static RecvSysnumCT   RecvSysnums_;
    static MappingIdxCT   tmpMappingIdx_;
    static SysnumProcCT   tmpSysProc_;
    static IdxDescCL*     RowIdx_;
    static int            maxNeighs_;

      // helper functions for creating the ExchangeCL
    template <typename SimplexIterT>
      void CollectSendSysNums(const SimplexIterT& begin, const SimplexIterT& end, VectorBaseCL<bool>& DistSysnums);
    template <typename T>
      static IdxT findPos(const std::vector<T>& a, const T& elem);
    void CreateExchangeDataMPIType();
    void TransferSendOrder(bool CreateMap);
    void CreateIndices(IdxDescCL*, const VectorBaseCL<bool>&, bool forAccParDot);

      // flags and sizes
    Ulint numLocalIdx_;    // number of local sysnums
    Ulint numDistrIdx_;    // number of distributed sysnums
    Ulint vecSize_;        // check, how long the vector for accumulation must be (cannot prevent all errors)
    Uint  numNeighs_;      // number of neighbors
    Ulint numAllRecvUnk_;  // size, that the receive buffer must has at least
    bool  created_;        // Flag for checking if the lists are created
    bool  mapCreated_;     // Flag if the mapping: (external idx) -> (my idx) is created
    bool  accIdxCreated_;  // Flag if AccDistIdx has been created
    int   tag_;            // internal used tag for MPI to create exchange lists

    mutable VectorCL recvBuf_;  // Standard buffer for receiving unknowns
    IdxVecT recvOffsets_;       // offsets for receiving (so each neighbor stores in different positions of the same receive buffer)

    inline double LocDot_(const VectorCL&, const VectorCL&, VectorCL* x_acc) const;                   // Inner Product of one accumulated and one distributed vector

    inline double AccurLocDotNoAcc_(const VectorCL&, const VectorCL&) const;                          // Accure inner product with no accumulation
    inline double AccurLocDotOneAcc_(const VectorCL&, const VectorCL&, VectorCL*) const;              // Accure inner product with already one accumulated vector
    inline double AccurLocDotBothAcc_(const VectorCL&, const VectorCL&, VectorCL*, VectorCL*) const;  // Accure inner product with two unaccumulated vectors

    void CreateList(const MultiGridCL& mg, IdxDescCL *RowIdx,
                    bool CreateMap=true, bool CreateAccDist=true);                      // create communication lists to an index

    /// \brief Create ExchangeCL due to a MLIdxDescCL. Only data on finest level can be exchanged!
    void CreateList(const MultiGridCL& mg, MLIdxDescCL *RowIdx, bool CreateMap=true, bool CreateAccDist=true){
        CreateList(mg, RowIdx->GetFinestPtr(), CreateMap, CreateAccDist);
    }

  public:
    ExchangeCL();                                                                       // just set the created-flag to false;
    ~ExchangeCL();                                                                      // delete ExList
    void clear();                                                                       // remove all information

    /// \name Helper function for DDD, should be private ...
    //@{
    template <typename SimplexT>
      static int HandlerGatherSysnums(OBJT objp, void* buf);
    template <typename SimplexT>
      static int HandlerScatterSysnums(OBJT objp, void* buf);
    //@}
    inline Ulint GetNumLocIdx()  const;                                                 // get the number of local sysnums
    inline Ulint GetNumDistIdx() const;                                                 // get the number of distributed sysnums
    inline Ulint GetNumDistAccIdx() const;                                              // get number of distributed sysnums, this proc is exclusively responsible for
    inline Ulint GetNum() const;                                                        // get the size of vector, that can be accumulated
    inline Ulint GetNumReceiveElements() const;                                         // get number of elements, that should be received (i.e. size of the receive buffer)
    inline bool  Created() const;                                                       // check if the list has been created
    inline bool  MapCreated() const;                                                    // check if the mapping has been created
    inline bool  AccIdxCreated() const;                                                 // check if index for accumulated inner products are set

    // start sending and receiving
    inline void   InitCommunication(const VectorCL&, RequestCT&, int tag=-1, Ulint offset=0, VectorCL* recvBuf=0) const;
    // finish communication and accumulate vector
    inline void   AccFromAllProc(VectorCL&, RequestCT&, Ulint offset=0, VectorCL* recvBuf=0) const;

    inline void                  Accumulate(VectorCL&) const;                           // Accumulate the Vector
    inline VectorCL              GetAccumulate (const VectorCL&) const;                 // Return accumulated Vector
    inline std::vector<VectorCL> GetAccumulate (const std::vector<VectorCL>&) const;    // Return accumulated vectors

    // Perform inner products (without and with global reduce)
    inline double LocDot    (const VectorCL&, bool, const VectorCL&, bool, bool useAccur=true, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    inline double ParDot    (const VectorCL&, bool, const VectorCL&, bool, bool useAccur=true, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    // Perform norms (without and with global reduce)
    inline double LocNorm_sq(const VectorCL&, bool, bool useAccur=true, VectorCL* r_acc=0) const;
    inline double Norm      (const VectorCL&, bool, bool useAccur=true, VectorCL* r_acc=0) const;
    inline double Norm_sq   (const VectorCL&, bool, bool useAccur=true, VectorCL* r_acc=0) const;

    // old interface should be removed
    inline double Norm_sq_Acc(VectorCL&, const VectorCL&) const;                        // x_acc^T * x (accumulates second to first parameter)
    inline double Norm_sq(const VectorCL&) const;                                       // \|x\|_2^2
    inline double Norm(const VectorCL&) const;                                          // returns the euclidian-norm of a vector

    inline double ParDotAcc(VectorCL&, const VectorCL&) const;                          // InnerProduct: first Vector will be accumulated after the procedure!
    inline double DotAcc(VectorCL&, const VectorCL&) const;                             // InnerProduct without global reduce. first Vector will be accumulated
    inline double ParDot(const VectorCL&, const VectorCL&) const;                       // InnerProduct: no accumulation of the input vectors but slower as function above
    inline double ParDot(VectorCL&, const VectorCL&, const VectorCL&) const;            // InnerProduct: store the accumulated second vector in the first parameter

    inline double AccParDot(const VectorCL&, const VectorCL&, VectorCL&, VectorCL&) const;// InnerProduct: Both vectors will be accumulated, this is more accurate
    inline double AccParDot(const VectorCL&, const VectorCL&, VectorCL&) const;         // InnerProduct: with two accumulated vectors. The first given vec should be accumulated
    inline double LocAccDot(const VectorCL&, const VectorCL&) const;                    // InnerProduct of two accumulated vectors without global reduce
    inline double AccNorm_sq(const VectorCL&, VectorCL&) const;                         // Norm of a distributed unaccumulated vector
    inline double AccNorm_sq(const VectorCL&) const;                                    // Norm of an accumulated vector
    inline double LocAccNorm_sq(const VectorCL&, VectorCL&) const;                      // Norm of a distributed unaccumulated vector without global reduce
    inline double LocAccNorm_sq(const VectorCL&) const;                                 // Norm of an accumulated vector without global reduce
    // end of old interface

    inline IdxT      GetExternalIdxFromProc(IdxT, ProcNumT) const;                      // Get index of a distributed index on another proc
    inline bool      IsDist(IdxT) const;                                                // Check if a sysnum is distributed
    inline bool      IsOnProc(IdxT,ProcNumT);                                           // Check if a sysnum can be found on another proc
    inline ProcNumCT GetProcs(IdxT) const;                                              // Get list of procs that owns a sysnum (except local proc)
    inline Uint      GetNumProcs(IdxT) const;                                           // Get number of procs, that owns a sysnum
    inline bool      IsExclusive(IdxT) const;                                           // Is a sysnum on the calling processor exclusive (i.e. this proc has the smallest proc id)

    inline ProcNumCT GetNeighbors() const;                                              // Get procs that shares at least one unknown with this proc
    inline Uint      GetNumNeighs() const;                                              // Get number of neighbor processes

      // Debugging and information
    void DebugInfo(std::ostream&) const;                                                // Debug Info
    void SizeInfo(std::ostream&, const int Proc=0) const;
    bool IsAcc(const VectorCL&) const;                                                  // Check if a vector is accumulated

    bool IsEqual(const ExchangeCL&, std::ostream*os=0) const;                           // for degubbing, check if to ExchangeCL'es seems to be equal
};

/// \name Wrapper for gathering and scattering data to create the ExchangeCL
//@{
extern "C" int HandlerGatherSysnumsVertexC(OBJT, void*);
extern "C" int HandlerScatterSysnumsVertexC(OBJT, void*);
extern "C" int HandlerGatherSysnumsEdgeC(OBJT, void*);
extern "C" int HandlerScatterSysnumsEdgeC(OBJT, void*);
//@}


/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/
/// \brief Handle exchange of all numerical data for a blocked vector, i.e.
///    vectors, that have multiple IdxDescCL
/** This class handles the exchange of a blocked vector containing multiple
    describers. This is used to perform a blocked version of an iterative
    solver. For example GCR can be used to solve the Oseen problem
 */
/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/

class ExchangeBlockCL
{
  public:
    typedef std::vector<const IdxDescCL*>      IdxDescCT;       ///< Container for IdxDescCL
    typedef std::vector<IdxT>                  BlockOffsetCT;   ///< Container of starting index of block elements
    typedef std::vector<ExchangeCL::RequestCT> VecRequestCT;    ///< Container for requests

  private:
    IdxDescCT            idxDesc_;      ///< store all index describers to access ExchangeCLs
    BlockOffsetCT        blockOffset_;  ///< store the length of vectors
    mutable VecRequestCT SendRecvReq_;  ///< standard requests for sending and receiving
    int                  startTag_;     ///< first Tag to be used for sending and receiving

    /// \brief start sending and receiving
    void InitCommunication(const VectorCL&, VecRequestCT&, int tag=-1, std::vector<VectorCL>* recvBuf =0) const;
    /// \brief finish communication and accumulate vector
    void AccFromAllProc(VectorCL&, VecRequestCT&, std::vector<VectorCL>* recvBuf=0) const;

    /// \brief Sum up local elements
    inline double SumUpLocal(const VectorCL&, const VectorCL&) const;
    /// \brief Sum up distributed elements
    inline double SumUpDist(const VectorCL&, const VectorCL&) const;

    /// \brief Accurate version of a local inner product with two given accumulated vector
    inline double AccurLocDotNoAcc(const VectorCL&, const VectorCL&) const;
    /// \brief Accurate version of a local inner product with one given accumulated vector
    inline double AccurLocDotOneAcc(const VectorCL&, const VectorCL&, VectorCL*) const;
    /// \brief Accurate version of a local inner product with no given accumulated vector
    inline double AccurLocDotBothAcc(const VectorCL&, const VectorCL&, VectorCL*, VectorCL*) const;
    /// \brief Inner product of one accumulated and one distributed vector
    inline double LocDot(const VectorCL&, const VectorCL&, VectorCL* x_acc) const;

  public:
    ExchangeBlockCL()
      : idxDesc_(), blockOffset_(), SendRecvReq_(), startTag_(1001) {}

    /// \brief Attach an index describer
    void AttachTo(const IdxDescCL&);
    /// \brief Ask for number of handled blocks
    size_t GetNumBlocks() const { return idxDesc_.size(); }
    /// \brief Ask for length of vectors, that can be accumulated
    IdxT GetNum() const { return blockOffset_.back(); }
    /// \brief Ask for an ExchangeCL
    const ExchangeCL& GetEx( size_t i) const { return idxDesc_[i]->GetEx(); }

    /// \brief Update of datastructure, i.e. blockoffset_
    void Update();

    /// \brief Perform an inner product without global reduction of the sum
    inline double LocDot    (const VectorCL&, bool, const VectorCL&, bool, bool useAccur=true, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    /// \brief Perform an inner product with global reduction of the sum
    inline double ParDot    (const VectorCL&, bool, const VectorCL&, bool, bool useAccur=true, VectorCL* x_acc=0, VectorCL* y_acc=0) const;
    /// \brief Perform squared Euklidian norm without global reduction of the sum
    inline double LocNorm_sq(const VectorCL&, bool, bool useAccur=true, VectorCL* r_acc=0) const;
    /// \brief Perform squared Euklidian norm with global reduction of the sum
    inline double Norm_sq   (const VectorCL&, bool, bool useAccur=true, VectorCL* r_acc=0) const;
    /// \brief Perform Euklidian norm with global reduction of the sum
    inline double Norm      (const VectorCL&, bool, bool useAccur=true, VectorCL* r_acc=0) const;
    /// \brief Accumulate the given vector
    inline void                  Accumulate(VectorCL&) const;
    /// \brief Return an accumulated vector
    inline VectorCL              GetAccumulate (const VectorCL&) const;
};


/****************************************************************************
* E X C H A N G E  M A T R I X  C L A S S                                   *
****************************************************************************/
/// \brief Handle the accumulation of a sparse matrix (MatrixCL)
/** This class is capable of determining the communication pattern for
    accumulating a sparse matrix, and performing the accumulation.
    \todo(par) Develope an "accure" version of accumulation
 */
/****************************************************************************
* E X C H A N G E  M A T R I X  C L A S S                                   *
****************************************************************************/
class ExchangeMatrixCL
{
  public:
    typedef ExchangeCL::ProcNumCT  ProcNumCT;       ///< Container for storing neighbor processes
    typedef ProcNumCT::iterator    ProcNum_iter;    ///< iterator of ProcNumCT
    typedef std::vector<size_t>    CouplingCT;      ///< Container of distributed matrix elements

  private:
    /// each element of ExList_ handles the send-process with a single neighbor processor
    std::vector<ExchangeDataSendCL>    ExList_;
    /// Buffer for receiving elements
    std::vector<VectorCL>              RecvBuf_;
    /// Where to add/store received non-zeroes
    std::vector<CouplingCT>            Coupl_;
    /// flag, if non-zero is not stored on local processor
    static size_t NoIdx_;

    /// Determine the intersection of two processor lists
    inline ProcNum_iter Intersect(ProcNumCT& a, ProcNumCT& b, ProcNumCT& result)
    /// Sort both lists and use the standard intersection algorithm
    {
        a.sort();
        b.sort();
        return std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), result.begin());
    }

    /// Determine the position, where a nonzero is stored
    inline size_t GetPosInVal(const size_t row, const size_t col, const MatrixCL& mat)
    /// if the non-zero (row,col) is not stored by the local processor, this function
    /// returns NoIdx_
    {
        Assert( row<mat.num_rows() && col<mat.num_cols(), DROPSErrCL("ExchangeMatrixCL::GetPosInVal: Row or col out of bounds"), DebugParallelNumC);
        const size_t *pos= std::lower_bound( mat.GetFirstCol(row), mat.GetFirstCol(row+1), col);
        return (pos != mat.GetFirstCol(row+1) && *pos==col) ? pos-mat.GetFirstCol(0) : NoIdx_;
    }

  public:
    // default constructors and destructors

    /// \brief Reset
    void Clear() { ExList_.clear(); RecvBuf_.clear(); Coupl_.clear(); }

    /// \brief Determine the communication pattern for accumulating a matrix
    void BuildCommPattern(const MatrixCL& mat, const IdxDescCL& RowIdx, const IdxDescCL& ColIdx){
        BuildCommPattern(mat, RowIdx.GetEx(), ColIdx.GetEx());
    }

    /// \brief Determine the communication pattern for accumulating a matrix
    void BuildCommPattern(const MatrixCL&, const ExchangeCL& RowEx, const ExchangeCL& ColEx);

    /// \brief Accumulate a matrix
    MatrixCL Accumulate(const MatrixCL&);
};

} // end of namespace DROPS

// File, where the inline an template-functions are declared
#include "parallel/exchange.tpp"

#endif
