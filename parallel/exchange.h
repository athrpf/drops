//**************************************************************************
// File:    exchange.h                                                     *
// Content: Class that handles the accumulation of the unknowns            *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:    May, 23th 2006                                                 *
//          - ExchangeBlockCL added                                        *
//          August, 25th 2006                                              *
//          - Sysnums on ohter procs are computed                          *
//          September, 7th 2006                                            *
//          - Inner products and norms can be computed on accumulated vecs *
// Begin:   March, 08th 2006                                               *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file exchange.h
/// \brief Exchange numerical Data and perform inner products
///
///          These classes do not use the DDD-Interfaces. After the lists   *
///          are created no geometric datas are needed to do the            *
///          accumulation in opposite to the DDD-Interface. And this class  *
///          split the send and the recieve, so other work can be done      *
///          between these commands.                                        *

#ifndef DROPS_EXCHANGE_H
#define DROPS_EXCHANGE_H

#include "parallel/parallel.h"
#include <list>
#include <vector>
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "misc/problem.h"

namespace DROPS{

class ExchangeCL;

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
class ExchangeDataCL
{
  friend class ExchangeCL;

  private:
    typedef std::vector<int> SysnumListCT;                                      // Sequence, where to store the recieved unknowns

    int                     toProc_;                                            // handle comm between "me" and "toProc_"
    ProcCL::DatatypeT       SendType_;                                          // Type for sending
    SysnumListCT            Sysnums_;                                           // sysnums of the recieved data
#ifdef DebugParallelNumC
    Ulint                   SendTypeSize_;                                      // how big the transfered vector must be at least
#endif

    void SetToProc(const int Proc);
    void CreateDataType(const int count, const int blocklength[], const int array_of_displacements[]);
    void CreateSysnums(const SysnumListCT&);                                    // Create SysnumListCT for recieving

  public:
    ExchangeDataCL(int proc);
    ExchangeDataCL();
    ExchangeDataCL(const ExchangeDataCL&);
    ~ExchangeDataCL();

    // Get rank of neighbor proc
    inline int GetProc() const;
    // Get number of received elements
    inline size_t GetNumRecvEntries() const;
    // Send data to "toProc_" (nonblocking, asynchronous)
    inline ProcCL::RequestT Isend(const VectorCL&, int tag, Ulint offset) const;
    // Recieve data (nonblocking)
    inline ProcCL::RequestT Irecv(int tag, VectorCL& recvBuf, Ulint offset) const;
    // add data from "toProc_"
    inline void Accumulate(VectorCL&, Ulint offsetV, VectorCL& recvBuf, Ulint offsetRecv) const;

    // print, where to store recieved unknowns
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
    Ulint vecSize_;        // check, how big the vector for accumulation mus be (cannot prevent all errors)
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

  public:
    ExchangeCL();                                                                       // just set the created-flag to false;
    ~ExchangeCL();                                                                      // delete ExList
    void clear();                                                                       // remove all information

    void CreateList(const MultiGridCL& mg, IdxDescCL *RowIdx,
                    bool CreateMap=true, bool CreateAccDist=true);                      // create communication lists to an index

    /// \brief Create ExchangeCL due to a MLIdxDescCL. Only data on finest level can be exchanged!
    void CreateList(const MultiGridCL& mg, MLIdxDescCL *RowIdx, bool CreateMap=true, bool CreateAccDist=true){
        CreateList(mg, RowIdx->GetFinestPtr(), CreateMap, CreateAccDist);
    }

    /// \name Helper function for DDD, should be private ...
    //@{
    template <typename SimplexT>
      static int HandlerGatherSysnums(DDD_OBJ objp, void* buf);
    template <typename SimplexT>
      static int HandlerScatterSysnums(DDD_OBJ objp, void* buf);
    //@}
    inline Ulint GetNumLocIdx()  const;                                                 // get the number of local sysnums
    inline Ulint GetNumDistIdx() const;                                                 // get the number of ditributed sysnums
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
extern "C" int HandlerGatherSysnumsVertexC(DDD_OBJ, void*);
extern "C" int HandlerScatterSysnumsVertexC(DDD_OBJ, void*);
extern "C" int HandlerGatherSysnumsEdgeC(DDD_OBJ, void*);
extern "C" int HandlerScatterSysnumsEdgeC(DDD_OBJ, void*);
//@}

/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/
/// \brief Handle exchange of all numerical data for a blocked vector
/** This class handles the exchange of a blocked vector containing m blocks.
    For each block an index-describer class is given.
    \todo (of) AccurLocDotBothAcc_ has to be fixed
 */
/****************************************************************************
* E X C H A N G E  B L O C K  C L A S S                                     *
****************************************************************************/
class ExchangeBlockCL
{
  public:
    typedef std::vector<ExchangeCL>   ExchangeCT;
    typedef std::vector<MLIdxDescCL*> MLIdxDescCT;
    typedef std::vector<Ulint>        BlockSizeCT;
    typedef VectorBaseCL<ExchangeCL::RequestCT> VecRequestCT;

  private:
    size_t       m_;                                                        // number of blocks
    bool         created_;                                                  // flag if all lists have been created
    bool         blockCreated_;                                             // flag if the blocks have been created
    ExchangeCT   ExchgList_;                                                // list of exchange classes
    BlockSizeCT  Block_;                                                    // (m+1) vector, that stores the first index of each block. Last index stores the size of handleable vectors
    int          start_tag_;                                                // standard tag for the first block. blocks get consecutive tags (=1001)
    mutable VecRequestCT SendRecvReq_;                                      // standard request handle for non-blocking sending and receiving

    void CreateBlockIdx();                                                  // calc and store the beginning indices of a blocked vector (all ExchangeCLs must have been created!)

    inline double LocDot_(const VectorCL&, const VectorCL&, VectorCL* x_acc) const;                   // Inner Product of one accumulated and one distributed vector

    inline double AccurLocDotNoAcc_(const VectorCL&, const VectorCL&) const;                          // Accure inner product with no accumulation
    inline double AccurLocDotOneAcc_(const VectorCL&, const VectorCL&, VectorCL*) const;              // Accure inner product with already one accumulated vector
    inline double AccurLocDotBothAcc_(const VectorCL&, const VectorCL&, VectorCL*, VectorCL*) const;  // Accure inner product with two unaccumulated vectors

  public:
    ExchangeBlockCL(size_t m);                                              // Construct a class that can store m blocks
    void clear();                                                           // remove all information
//  ~ExchangeBlockCL();

    void CreateList(const MultiGridCL&, const MLIdxDescCT&,                 // Create exchange-lists for all blocks
                    bool CreateMap=true, bool CreateAccDist=true);
    void CreateList(const MultiGridCL&, size_t i, MLIdxDescCL*,             // Create exchange-list for i-th block
                    bool CreateMap=true, bool CreateAccDist=true);

    inline const ExchangeCL& Get(size_t i) const;                           // Get the i-th exchange block corresponding to the i-th block
    inline       ExchangeCL& Get(size_t i);
    inline size_t      GetNumBlocks() const;                                // Get number of blocks
    inline bool        Created();                                           // check if all exchange-lists have been created;
    inline Ulint       GetNum() const;                                      // get the size of vector, that can be accumulated
    inline Ulint       GetNum(size_t i) const;                              // get the size of the i-th block

    /// \name neighborhood communication on blocked vectors
    /// \brief for detailed describtions see ExchangeCL
    // @{
    // start sending and receiving
    inline void   InitCommunication(const VectorCL&, VecRequestCT&, int tag=-1, std::vector<VectorCL>* recvBuf=0) const;
    // finish communication and accumulate vector
    inline void   AccFromAllProc(VectorCL&, VecRequestCT&, std::vector<VectorCL>* recvBuf=0) const;

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

    inline double Norm_sq_Acc(VectorCL&, const VectorCL&) const;            // x_acc^T * x (accumulates second to first parameter)
    inline double Norm(const VectorCL&) const;                              // returns the euclidian-norm of a vector
    inline double Norm_sq(const VectorCL&) const;                           // \|x\|_2^2

    inline double ParDotAcc(VectorCL&, const VectorCL&) const;              // InnerProduct: first Vector will be accumulated after the procedure!
    inline double DotAcc(VectorCL&, const VectorCL&) const;                 // InnerProduct without global reduce. first Vector will be accumulated
    inline double ParDot(const VectorCL&, const VectorCL&) const;           // InnerProduct: no accumulation of the input vectors but slower as function above
    inline double ParDot(VectorCL&, const VectorCL&, const VectorCL&) const;// InnerProduct: store the accumulated second vector in the first parameter

    inline double AccParDot(const VectorCL&, const VectorCL&, VectorCL&, VectorCL&) const;// InnerProduct: Both vectors will be accumulated, this is more accurate
    inline double AccParDot(const VectorCL&, const VectorCL&, VectorCL&) const;         // InnerProduct: with two accumulated vectors. The first given vec should be accumulated
    inline double AccNorm_sq(const VectorCL&, VectorCL&) const;                         // Norm of a distributed unaccumulated vector
    inline double AccNorm_sq(const VectorCL&) const;                                    // Norm of an accumulated vector
    inline double LocAccDot(const VectorCL&, const VectorCL&) const;                    // InnerProduct of two accumulated vectors without global reduce
    inline double LocAccNorm_sq(const VectorCL&, VectorCL&) const;                      // Norm of a distributed unaccumulated vector without global reduce
    inline double LocAccNorm_sq(const VectorCL&) const;                                 // Norm of an accumulated vector without global reduce
    // @}
};

} // end of namespace DROPS

// File, where the inline an template-functions are declared
#include "parallel/exchange.tpp"

#endif
