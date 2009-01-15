//**************************************************************************
// File:    parmultigrid.h                                                 *
// Content: Class   that constitute the parallel multigrid                 *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, RZ RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   November, 14th, 2005                                           *
//**************************************************************************
// remarks:                                                                *
//  + no unknown on faces can be transfered. See doxygen!                  *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file parmultigrid.h
/// \brief Handling of parallel MultiGird data structures
#ifndef DROPS_PAR_MULTIGRID_H
#define DROPS_PAR_MULTIGRID_H

#include <ddd.h>

#include "parallel/loadbal.h"
#include "parallel/addeddata.h"
#include "parallel/interface.h"

#include "misc/utils.h"
#include "misc/problem.h"

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "geom/boundary.h"
#include "geom/topo.h"
#include "num/unknowns.h"

namespace DROPS
{
/****************************************************************************
* P A R  M U L T I  G R I D  C L A S S                                      *
****************************************************************************/
/// \brief Constitute a multigrid on several procs
///
/// This is the major class to constitute a multigrid on several procs
/// It handles the transfers and the accumulation of numerical and geometric
/// values. This class makes intensiv use of the DDD-Libary, which is
/// developed for UG.
/// See [],[],[] for more information and how to use this class within DROPS
///
/// This class cannot handle unknwons on faces. If you want to use them,
/// you have to change the MultiGridCL (adding an UnknownHandleCL to the FaceCL)
/// and this class by writing an HandlerFGather and HandlerFScatter.
/// \todo(of) AttachTo: für VecDesc und Bnd in eine Funktion packen. Geht das?
/****************************************************************************
* P A R  M U L T I  G R I D  C L A S S                                      *
****************************************************************************/
class ParMultiGridCL
{
  public: // typedefs
      /// \brief Vector of pointers to  Vector-Describer-Classes.
      ///
      /// Used for storing the recieved numerical data and creating new numbering on the vertices, edges and tetras
      /// Each type of unknown has one VecDescCL* (like pressure or velocity)
    typedef std::vector<VecDescCL*>             VecDescPCT;
     /// \brief Buffer for double values
     ///
     /// Buffer type for storing double values. Used within refinement and migration
    typedef std::vector< double >               BufferCT;
     /// \brief Vector of TetraCL pointer for remebering special tetras
     ///
     /// Used for remebering Tetras, that have changed prio from ghost to master, to make MFR
     /// on edges consistent
    typedef std::vector<TetraCL*>               TetraPCT;
    /// \brief Vector of scalar boundary conditions
    typedef std::vector< const BndDataCL<double>* >    ScalBndCT;
    /// \brief Vector of vectorial boundary conditions
    typedef std::vector< const BndDataCL<Point3DCL>* > VecBndCT;

  private:    // variables
    static DDD_IF _EdgeIF;                  // Interface for accumulate the "Marked for Refinement"
    static DDD_IF _TetraIF;                 // Interface for communicateRefMarks() and TreatGhosts()

    static DDD_IF        _FaceIF;           // |
    static std::ostream* _os;               //  >  for IsSane()
    static bool          _sane;             // |

    static DDD_TYPE _BndPtT,                // Boundary-Point-Type for Xfer_AddData                 (no header!)
                    _ChildPtrT;             // Childer-Pointer-Type for Xfer Children of a Tetra    (no header!)

    static DDD_IF        NotMasterIF_;      // Vertices, edges and tetras that are not masters (for destroying unknowns after migration)

    static DDD_IF        _GhToMaTetra_IF;   // Ghost tetras to master-tetras

    static MultiGridCL*     _mg;            // Pointer to the stored multigrid
    static MG_VertexContT*  _VertCont;      // Simplices stored on this proc
    static MG_EdgeContT*    _EdgeCont;
    static MG_FaceContT*    _FaceCont;
    static MG_TetraContT*   _TetraCont;

    static TetraPCT         ToHandleTetra_; // Remeber tetras that have changed prio from ghost to master within transfer.

    static VecDescPCT     _VecDesc;         // Vector of Pointers to the vector describers, where the unknowns are stored
    static BufferCT       _RecvBuf;         // Buffer for the recieved numerical data (used for refinement and migration)!
    static ScalBndCT      _ScalBnd;         // Store scalar boundary conditions
    static VecBndCT       _VecBnd;          // Store vectorial boundary conditions

    static bool TransferMode,               // If an active DDD-Xfer-Mode is activ
                PrioChangeMode;             // if we are in a DDD-PrioChange enviroment
    static IdxT _RecvBufPos;                // last postion in recieve buffer that is not used;

    static int _level;                      // all level to _level in XferStart() and XferEnd()
    static bool _UnkOnSimplex[3];           // Remeber if there are unknowns on (Vertex,Edge,Tetra)
    /// \todo _UnkOnSimplex wird nicht mehr unbedingt gebraucht. Besser direkt über _VecDesc abfragen.

    static VecDescCL* _actualVec;           // actual index within HandleNewIdx (used within DDD-Gather and DDD-Scatter operation)


  public:
    /// \name Constructor and Destructor
    // @{
    ParMultiGridCL();
    ParMultiGridCL(const MGBuilderCL&);
    ParMultiGridCL(const ParMultiGridCL&);                                  // Copyconstructor is not implemented!
    ~ParMultiGridCL();
    // @}


    /// \name Functions concerning the MultiGrid
    // @{
    static void AttachTo(MultiGridCL&);                                     // attach the MultiGrid to the ParMultiGridCL
    static MultiGridCL& GetMG();                                            // Get a reference on the MultiGrid
    void Refine();                                                          // Refine the MultiGrid
    static void AdjustLevel();                                              // Apply all the same number of levels to all procs
    void MarkAll();                                                         // All Tetras of last level are marked
    static void MarkSimplicesForUnknowns(int Level=-1);                     // Set prio PrioHasUnk on all simplices that are able to store unknowns
    // @}


    /// \name Functions for transfers
    // @{
    static void XferStart(int Level=-1);                                    // Call this everytime before using an transfer command!
    static void XferEnd();                                                  // Call this everytime after all tetras are marked for xfer
    static void TXfer(TetraCL&, DDD_PROC, DDD_PRIO, bool del);              // Transfer a tetra to another proc
    // @}


    /// \name Functions concerning the handling of unknowns
    // @{
    template<typename BndT>
    void AttachTo(VecDescCL*, const BndT*);                                 // attach VecDescCL and boundary conditions
    void DeleteVecDesc();
    static inline bool UnknownsOnSimplices();                               // are there unknowns on simplices
    void HandleUnknownsAfterRefine();                                       // Handle Unknowns after a refine operation
    void HandleNewIdx(IdxDescCL* old, VecDescCL* newVec);                   // Handle unknowns after migration and refinement
    /// \brief Handle unknowns after migration and refinement with a MLIdxDescCL. This is just a wrapper for the IdxDescCL
    void HandleNewIdx(MLIdxDescCL* old, VecDescCL* newVec) {
        HandleNewIdx( old->GetFinestPtr(), newVec);
    }
    void CompleteRepair(VecDescCL* newVec);                                 // After RepairAfterRefine[P1|P2] call this function
    void DelAllUnkRecv();                                                   // Clear all RecieveUnkowns-Marks
    void DeleteRecvBuffer();                                                // Delete recieve Buffer
    static void DeleteUnksOnGhosts(int Level=-1);                           // Delete Unknowns on ghosts and vertical ghosts
    // @}


    /// \name Access to unknowns that are stored in a recieve buffer or an old index
    // @{
    template<typename SimplexT, typename DofT>
    struct GetDof{
      inline DofT operator() (const SimplexT&, Uint idx, int pos=-1);       // Get dof on a simplex (not implemented, see specializations down)
    };
    template<typename SimplexT>
    struct GetDof<SimplexT, Point3DCL>{
      inline Point3DCL operator() (const SimplexT&, Uint idx, int pos=-1);  // Spezialisation for vectorial data
    };
    template<typename SimplexT>
    struct GetDof<SimplexT, double>{
      inline double operator() (const SimplexT&, Uint idx, int pos=-1);     // Spezialisation for scalar data
    };
    // @}


    /// \name Checking and debug functions
    // @{
    template<class SimplexT> static bool IsOnProcBnd( const SimplexT* s);   // If a simplex lies on a proc-boundary
    static bool IsSane(std::ostream&, int Level= -1);                       // Check if distributed edges and faces have the same subsimlices
    static void DebugInfo(std::ostream&);                                   // writes usefull infos onto outputstream
    static void Show(DDD_GID gid, char *mesg, int proc= -1);                // Show the simplex with a given GID
    static void ShowTetraIF();                                              // Show the Tetra-Interface
    static void ShowEdgeIF();                                               // Show the Edge-Interface
    void ShowTypes() const;                                                 // Display all types that are defined and declared
    void ShowInterfaces() const;                                            // Display all interfaces used by ParMultiGridCL
    static void ConsCheck();                                                // DDD-Consisty-Check
    static double GetBalance();                                             // Calculate Imbalance of triangulation over procs on last triangulation-level --> Communication
    static void ShowVecDesc();
    static size_t GetRecvBufferSize();                                      // Get the size of _RecvBuf
    // @}


    /// \name Handlers for DDD-Interfaces, Priority-Enviroment and Identification-Enviroment
    // @{
    static void AccumulateMFR(int Level=-1);                                    // accumulate mark for refinements on Edges on Level (-1==all levels!)
    static void CommunicateRefMarks( Uint Level);                               // Tell everybody, if an tetra is marked for refinement
    static void TreatGhosts (int Level= -1);                                    // Rescue all ghosts-subsimplices
    static void RescueGhostVerts(Uint Level);                                   // Vertices are special
    static void TreatHasGhosts (int Level= -1);                                 // Rescue all subsimplices of tetra that has ghosts
    static void AdaptPrioOnSubs();                                              // Set prios of all subsimplices right
    static void RescueSubs(TetraCL&);                                           // All subsimplices of tetra are rescued and get prio PrioMaster
    static void AdaptMidVertex (int Level= -1);                                 // Set MidVertices with vertical overlapping
    template<class SimplexT>
    inline static void PrioChange(SimplexT* const, Uint Prio);                  // Change prio of a simplex
    static void IdentifyVertex( const EdgeCL* Parent);                                              // Identify an vertex by parant edge
    static void IdentifyEdge( EdgeCL* Me, const EdgeCL* Parent, Uint);                              // for subedge Me of parent edge Parent
    static void IdentifyEdge( EdgeCL* Me, const FaceCL* Parent, const VertexCL*, const VertexCL*);  // for subedge Me in parent face Parent
    static void IdentifyFace( FaceCL* Me, const FaceCL* Parent, Uint);
    // @}


  private:
    // functions concerning the interfaces and handlers
    // ------------------------------------------------
    void DeclareAll();                                              // Declare all types
    void DefineAll();                                               // Define all types
    void SetAllHandler();                                           // Set all Handlers to the wrapped handlers from this class
    void InitIF();                                                  // tell DDD about the interfaces

    static void VXfer(VertexCL&, DDD_PROC, DDD_PRIO, bool del);     // transfer a vertex, this can damage the MultiGrid structure
    static void EXfer(EdgeCL&, DDD_PROC, DDD_PRIO, bool del);       // transfer an edge, this can damage the MultiGrid structure
    static void FXfer(FaceCL&, DDD_PROC, DDD_PRIO, bool del);       // not implemented

    // declare and definetypes from this class
    static void DeclareBndPtT();                                    // Declare BoundaryPointer-Type
    static void DeclareChildPtrT();                                 // Declare ChildPointer-Type
    static void DefineBndPtT();                                     // Define BoundaryPointer-Type
    static void DefineChildPtrT();                                  // Define ChildPointer-Type


    // functions concerning the internal handling of unknowns
    // ------------------------------------------------------

    // checking and size-estimating functions
    //@{
    static inline bool VecDescRecv();                               // VecDesc revieved?
    static Uint NumberOfUnknownsOnTetra();                          // Get overall number of unknowns on a tetrahedron
    static Uint NumberOfUnknownsOnVertex();                         // Get number of unknowns on a vertex
    static Uint NumberOfUnknownsOnEdge();                           // Get number of unknowns on an edge
    //@}


    // Copy values and doing linear interpolation after refine/migrate
    //@{
    static inline Uint GetStorePos(const IdxDescCL*);                       // Get position where the IdxDesc is internally stored
    template<typename BndT>
    void AttachTo(const IdxDescCL*, const BndT*);                           // attach boundary conditions
    template<typename BndT>
    static inline bool LinearInterpolation(const EdgeCL&, Uint, const BndT*, const VectorCL&, typename BndT::bnd_type& new_dof);
    static void RescueUnknownsOnEdges();                                    // Rescue unknowns on edges, that are deleted
    void PutData(MultiGridCL::const_VertexIterator&,
                 const VectorCL* const, VectorCL*, const Uint,
                 const Uint, const IdxDescCL*);                             // Copy unknowns on a vertex into a new datafield
    template<typename BndT>
    void PutData(MultiGridCL::const_EdgeIterator&,
                 const VectorCL* const, VectorCL*, const Uint,
                 const Uint, const IdxDescCL*, const BndT*);                // Copy unknowns on an edge into a new datafield
    template<typename BndT>
    inline static const BndT* GetBndCond(const IdxDescCL*);                 // Get boundary condition to store VecDesCL
    //@}

    // Send and recieve unknowns
    //@{
    template<class SimplexT>
    static inline void SendUnknowns(SimplexT*, DDD_TYPE, void*, int);       // Send Unknwons within the Handler<SimplexT>Gather
    template<class SimplexT>
    static inline void RecvUnknowns(SimplexT*, DDD_TYPE, void*, int);       // Recieve Unknwons within the Handler<SimplexT>Gather
    static void EnlargeRecieveBuffer();                                     // Enlarge the _RecvBuf
    //@}


    // set and get functions for vectors
    //@{
    template<typename SimplexT>
    static void SetDof(const SimplexT&, Uint, VectorCL&, const Point3DCL&); // Put data into a given vector
    template<typename SimplexT>
    static void SetDof(const SimplexT&, Uint, VectorCL&, const double&);    // Put data into a given vector

    template<typename SimplexT>
    static void PutDofIntoRecvBuffer(SimplexT&, Uint, const Point3DCL&);    // Put a value of unknown into the Recieve Buffer
    template<typename SimplexT>
    static void PutDofIntoRecvBuffer(SimplexT&, Uint, const double&);       // Put values of unknown into the Recieve Buffer

    template<typename SimplexT, typename ContainerT, typename DofT>
    struct GetDofOutOfVector{
      inline DofT operator() (const SimplexT&, Uint, const ContainerT&);    // Get values of unknown out of a given vector
    };
    template<typename SimplexT, typename ContainerT>
    struct GetDofOutOfVector<SimplexT, ContainerT, Point3DCL>{
      inline Point3DCL operator() (const SimplexT&, Uint, const ContainerT&);// Get values of unknown out of a given vector
    };
    template<typename SimplexT, typename ContainerT>
    struct GetDofOutOfVector<SimplexT, ContainerT, double>{
      inline double operator() (const SimplexT&, Uint, const ContainerT&);  // Get values of unknown out of a given vector
    };
    //@}


  public:
    /// \name DDD-Handler (should be private)
    //@{
    static void DeleteObj(void *buffer, size_t size, int ddd_typ);          // see cpp-File for further information
    template<class SimplexT> static void HandlerDelete( DDD_OBJ);           // DDD-Handler for DDD_XferDeleteObj

    static DDD_OBJ HandlerVConstructor( size_t, DDD_PRIO, DDD_ATTR level);  // how the reciever can construct a VertexCL
    static void HandlerVXfer( DDD_OBJ , DDD_PROC, DDD_PRIO);                // how the sender can send a VertexCL
    static void HandlerVGather( DDD_OBJ, int, DDD_TYPE, void*);             // how to send the boundary-information and numerical Data
    static void HandlerVScatter( DDD_OBJ, int, DDD_TYPE, void*, int);       // how to recieve the boundary-information and numerical Data

    static DDD_OBJ HandlerEConstructor( size_t, DDD_PRIO, DDD_ATTR level);  // How to construct an edge,
    static void HandlerEXfer(DDD_OBJ, DDD_PROC, DDD_PRIO);                  // to send one,
    static void HandlerEGather( DDD_OBJ, int, DDD_TYPE, void*);             // to pack numerical data to message,
    static void HandlerEScatter( DDD_OBJ, int, DDD_TYPE, void*, int);       // and to unpack this data

    static DDD_OBJ HandlerFConstructor( size_t, DDD_PRIO, DDD_ATTR level);  // How to construct a face
    static void HandlerFXfer(DDD_OBJ, DDD_PROC, DDD_PRIO);                  // not implemented

    static DDD_OBJ HandlerTConstructor( size_t, DDD_PRIO, DDD_ATTR level);  // How to construct a tetraeder,
    static void HandlerTXfer(DDD_OBJ, DDD_PROC, DDD_PRIO);                  // to send one
    static void HandlerTGather( DDD_OBJ, int, DDD_TYPE, void*);             // to pack numerical data to message,
    static void HandlerTScatter( DDD_OBJ, int, DDD_TYPE, void*, int);       // to unpack this data
    static void HandlerTUpdate(DDD_OBJ);                                    // to link tetra to faces
    static void HandlerTObjMkCons( DDD_OBJ, int);                           // to set right MFR values onto the edges
    static void HandlerTSetPrio(DDD_OBJ, DDD_PRIO);                         // to set the priority

    static int GatherEdgeMFR( DDD_OBJ, void*);                              // send MFR from Edge
    static int ScatterEdgeMFR( DDD_OBJ, void*);                             // collect MFR on Edge

    static int GatherTetraRestrictMarks( DDD_OBJ, void*);                   // send if a tetra is marked for refinement
    static int ScatterTetraRestrictMarks( DDD_OBJ, void*);                  // recieve, if a tetra is marked for refinement

    static int ExecGhostRescue( DDD_OBJ);                                   // _TetraIF calls this functions with TreatGhosts()
    static int ExecGhVertRescue( DDD_OBJ);                                  // Vertices has to be treaten special!
    static int ExecHasGhost( DDD_OBJ);                                      // _TetraIF calls this functions with TreatHasGhosts()
    static int ExecAdaptVGhostMidVertex( DDD_OBJ);                          // _EdgeIF calls this function with AdaptMidVertex()

    static int GatherEdgeSane( DDD_OBJ, void*, DDD_PROC, DDD_ATTR);         // Check if Edges are sane
    static int ScatterEdgeSane( DDD_OBJ, void*, DDD_PROC, DDD_ATTR);
    static int GatherFaceSane( DDD_OBJ, void*, DDD_PROC, DDD_ATTR);         // Check if Faces are sane
    static int ScatterFaceSane( DDD_OBJ, void*, DDD_PROC, DDD_ATTR);

    template<class SimplexT> static int DestroyUnksOnSimplex(DDD_OBJ);      // Destroy Unks on a simplex
    static int GatherUnknownsRef (DDD_OBJ, void*);                          // Gather  all unknowns on tetra with subsimplices for RepairAfterRefine
    static int ScatterUnknownsRef(DDD_OBJ, void*);                          // Scatter all unknowns on tetra with subsimplices for RepairAfterRefine

    static int GatherUnknownsMigV (DDD_OBJ, void*);                         // Gather  all unknowns on a vertex for sending unknowns to a vertex that has changed prio from vghos/ghost to PrioHasUnk
    static int ScatterUnknownsMigV(DDD_OBJ, void*);                         // Scatter all unknowns on a vertex for sending unknowns on a vertex that has changed prio from vghos/ghost to PrioHasUnk
    static int GatherUnknownsMigE (DDD_OBJ, void*);                         // Gather  all unknowns on a edge for sending unknowns to a edge that has changed prio from vghos/ghost to PrioHasUnk
    static int ScatterUnknownsMigE(DDD_OBJ, void*);                         // Scatter all unknowns on a egde for sending unknowns on a edge that has changed prio from vghos/ghost to PrioHasUnk

    template<typename SimplexT>
    static int GatherInterpolValues (DDD_OBJ, void*);                       // Gather  unknowns of an interpolated simplex
    template<typename SimplexT>
    static int ScatterInterpolValues(DDD_OBJ, void*);                       // Scatter unknowns of an interpolated simplex
    //@}
};

// Declaration of spezialized template functions
//----------------------------------------------
template<>
void ParMultiGridCL::AttachTo<BndDataCL<double> >(const IdxDescCL*, const  BndDataCL<double>*);
template<>
void ParMultiGridCL::AttachTo<BndDataCL<Point3DCL> >(const IdxDescCL*, const  BndDataCL<Point3DCL>*);
template<>
  const BndDataCL<double>* ParMultiGridCL::GetBndCond<BndDataCL<double> >( const IdxDescCL*);
template<>
  const BndDataCL<Point3DCL>* ParMultiGridCL::GetBndCond<BndDataCL<Point3DCL> >( const IdxDescCL*);


// Declaration of DDD-Handlers used within the template functions
//---------------------------------------------------------------
extern "C" int GatherInterpolValuesVC (DDD_OBJ o, void* b);
extern "C" int ScatterInterpolValuesVC(DDD_OBJ o, void* b);
extern "C" int GatherUnknownsMigVC (DDD_OBJ o, void* b);
extern "C" int ScatterUnknownsMigVC(DDD_OBJ o, void* b);
extern "C" int GatherUnknownsMigEC (DDD_OBJ o, void* b);
extern "C" int ScatterUnknownsMigEC(DDD_OBJ o, void* b);


//  Helper classes and functions
//------------------------------

/// \brief For transfering unknowns of a killed ghost a marker has to be sent as well
///        to distinguish if this is a valid unknown
struct TransferUnkT
{
    double val;         ///< value
    Uint   idx;         ///< index
    bool   mark;        ///< marker
};

const int MIG=0, REF=1;         // constants for choosing the the filename in PrintMG

/// \brief Write Debug information into the file output/(proc_id)_MG_(REF|MIG)_(step).mg
void PrintMG(const ParMultiGridCL&, int type=MIG);

/// \brief Check parallel multigrid and write output into the file output/sane_(proc_id).chk
bool CheckParMultiGrid(const DROPS::ParMultiGridCL&);

} // end of namespace DROPS

#include "parallel/parmultigrid.tpp"     // implementation of the template or/and inline functions
#endif
