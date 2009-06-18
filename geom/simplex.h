/// \file
/// \brief simplex.h
/// \author Sven Gross, Joerg Grande, Patrick Esser, Eva IGPM, Oliver Fortmeier, SC

#ifndef DROPS_SIMPLEX_H
#define DROPS_SIMPLEX_H

#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include "misc/utils.h"
#include "geom/boundary.h"
#include "geom/topo.h"
#include "num/unknowns.h"

#ifdef _PAR
#  include "parallel/parallel.h"
#  include <ddd.h>
#  include <compiler.h>
#  include <map>
#  include <parmetis.h>   // for indextype
#endif

namespace DROPS
{

#ifdef _PAR
/// \enum Priority Priorities for distributed objects
enum Priority
{
    PrioNeutral     = 0,     ///< no priotity
    PrioKilledGhost = 1,     ///< prio of ghost tetras marked for removement
    PrioVGhost      = 2,     ///< for subs of overlapping master tetra (don't constitute a proc bnd)
    PrioGhost       = 3,     ///< for ghost tetras and subs only owned by ghost tetras (skipped by the public iterators)
    PrioMaster      = 4,     ///< master copy of tetras and their subs
    PrioHasUnk      = 5      ///< simplices having potentionally unknowns
};
#endif

// fwd decl
class BoundaryCL;
class PeriodicEdgesCL;

// Classes for simplices in the triangulation
class VertexCL;
class EdgeCL;
class FaceCL;
class TetraCL;

// Containers for storage of the simplices
typedef GlobalListCL<VertexCL> MG_VertexContT;
typedef GlobalListCL<EdgeCL>   MG_EdgeContT;
typedef GlobalListCL<FaceCL>   MG_FaceContT;
typedef GlobalListCL<TetraCL>  MG_TetraContT;


//**************************************************************************
// Classes: FaceWrapperCL, RecycleBinCL                                    *
// Purpose: When the triangulation is modified, some simplices will have   *
//          to be removed in one step of the algorithm, but will be needed *
//          again some steps later. In order not to lose the information   *
//          in the 'Unknowns', we store those simplices in a 'RecycleBin'. *
//          To identify a simplex in the 'RecycleBin' we use its vertices. *
//          'Edges' and 'Tetras' know their vertices, but 'faces' do not.  *
//          Therefore, we provide a 'FaceWrapper' that adds its 'Vertices' *
//          to a 'Face'.                                                   *
//**************************************************************************

struct FaceWrapperCL
{
  public:
    const FaceCL*   face;
    const VertexCL* vert1;
    const VertexCL* vert2;

    FaceWrapperCL(const FaceCL* f, const VertexCL* v1, const VertexCL* v2)
        : face(f), vert1(v1), vert2(v2) {}
};


class RecycleBinCL
{
  public:
    typedef std::list<const EdgeCL*>  EdgeContT;
    typedef std::list<FaceWrapperCL>  FaceWrapperContT;
    typedef std::list<const TetraCL*> TetraContT;

  private:
    EdgeContT        _Edges;
    FaceWrapperContT _Faces;
    TetraContT       _Tetras;

  public:
    inline const EdgeCL*  FindEdge  (const VertexCL* v) const;
    inline const FaceCL*  FindFace  (const VertexCL* v1, const VertexCL* v2) const;
    inline const TetraCL* FindTetra (const VertexCL* v1, const VertexCL* v2, const VertexCL* v3) const;

    void Recycle (const EdgeCL* ep)  { _Edges.push_back(ep); } // ??? nicht alles doppelt speichern ???
    void Recycle (const FaceCL* fp, const VertexCL* v1, const VertexCL* v2)
        { _Faces.push_back( FaceWrapperCL(fp,v1,v2) ); }
    void Recycle (const TetraCL* tp) { _Tetras.push_back(tp); }

    void DebugInfo(std::ostream&) const;
};


/*******************************************************************
*   V E R T E X  C L                                               *
*******************************************************************/
/// \brief Represents vertices in the multigrid
/** Contains the geometric part ('_Coord', '_BndVerts') of a
    point in the multigrid as well as some topological ('_Level')
    information.
    It also stores some algorithmic information like '_RemoveMark'
    and the RecycleBins.                                          */
/*******************************************************************
*   V E R T E X  C L                                               *
*******************************************************************/
class VertexCL
{
    friend class MultiGridCL;
#ifdef _PAR
    friend class ParMultiGridCL;
#endif

  public:
    typedef std::vector<BndPointCL>::iterator       BndVertIt;
    typedef std::vector<BndPointCL>::const_iterator const_BndVertIt;

  private:
    IdCL<VertexCL>           _Id;                                               // id of the vertex on this proc
    Point3DCL                _Coord;                                            // global coordinates of the vertex
    std::vector<BndPointCL>* _BndVerts;                                         // Parameterdarstellung dieses Knotens auf evtl. mehreren Randsegmenten
    RecycleBinCL*            _Bin;                                              // recycle-bin
    bool                     _RemoveMark;                                       // flag, if this vertex should be removed
#ifndef _PAR
    Uint                     _Level : 8;                                        // in parallel mode, the level is stored in _dddH.attr
#else
    // parallel
    static DDD_TYPE          _dddT;                                             // DDD-Type for Verts, for every Vertex the same
    DDD_HEADER               _dddH;                                             // VertexCL is a distributed DDD-Object
#if DROPSDebugC&DebugSubscribeC
    mutable bool subscribed_;
#endif
    VertexCL()                                                                  // standard constructor used for DDD to create a vertex
        : _Id(), _Coord(0.), _BndVerts(0), _Bin(0), _RemoveMark(false)
    {
        DDD_MarkHdrInvalid(&_dddH);                                             // The header will be constructed by DDD
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
    }
    static void Declare();                                                      // Declare the Vertex-Type by DDD ("parallel/parmultigrid.cpp")
    static void Define();                                                       // Define the Vertex Type by DDD  ("parallel/parmultigrid.cpp")
#endif
    mutable bool _needed;                                                       // Check within IsSane, if the Simples is needed by tetra or is forgotten to delete

  public:
    UnknownHandleCL          Unknowns;                                          ///< access to the unknowns on the vertex

  private:
    // RecycleBin
    bool                HasRecycleBin       () const { return _Bin; }           ///< check if recycle-bin is not empty
    const RecycleBinCL* GetRecycleBin       () const { return _Bin; }
          RecycleBinCL* GetCreateRecycleBin ()       { return HasRecycleBin() ? _Bin : (_Bin= new RecycleBinCL); }

// ===== Interface for refinement algorithm =====

  public:
    inline  VertexCL (const Point3DCL& Point3D, Uint FirstLevel, IdCL<VertexCL> id= IdCL<VertexCL>()); ///< create a vertex by coordinate and level; FileBuilderCL has to construct the _Id, too, thus it can optionally be set.
    inline  VertexCL (const VertexCL&);                                         ///< Danger!!! Copying simplices might corrupt the multigrid structure!!!
    inline ~VertexCL ();                                                        ///< also deletes the recycle-bin and boundary information

    // Boundary
    inline void AddBnd  (const BndPointCL& BndVert);                            ///< add boundary-information
    inline void BndSort ();                                                     ///< sort boundary-segments

    // RemovementMarks
    bool IsMarkedForRemovement () const { return _RemoveMark; }                 ///< check if vertex is marked for removement
    void SetRemoveMark         ()       { _RemoveMark = true; }                 ///< set mark for removement
    void ClearRemoveMark       ()       { _RemoveMark = false; }                ///< clear mark for removement

    // RecycleBin
    void DestroyRecycleBin () { delete _Bin; _Bin= 0; }                         ///< empty recycle-bin

    void Recycle (const EdgeCL* ep)                                             ///< put a pointer to an edge into the recycle-bin of this vertex
      { GetCreateRecycleBin()->Recycle(ep); }
    void Recycle (const FaceCL* fp, const VertexCL* vp1, const VertexCL* vp2)   ///< put a pointer to a face into the recycle-bin of this vertex
      { GetCreateRecycleBin()->Recycle(fp,vp1,vp2); }
    void Recycle (const TetraCL* tp)                                            ///< put a pointer to a tetra into the recycle-bin of this vertex
      { GetCreateRecycleBin()->Recycle(tp); }

    /// Find an edge in the recycle-bin by the opposite vertex. Returns 0 if no edge is found.
    EdgeCL*  FindEdge  (const VertexCL* v) const
      { return HasRecycleBin() ? const_cast<EdgeCL*>(GetRecycleBin()->FindEdge(v))          : 0; }
    /// Find a face in the recycle-bin by the other vertices. Returns 0 if no face is found.
    FaceCL*  FindFace  (const VertexCL* v1, const VertexCL* v2) const
      { return HasRecycleBin() ? const_cast<FaceCL*>(GetRecycleBin()->FindFace(v1,v2))      : 0; }
    /// Find a tetra in the recycle-bin by the other vertices. Returns 0 if no tetra is found.
    TetraCL* FindTetra (const VertexCL* v1, const VertexCL* v2, const VertexCL* v3) const
      { return HasRecycleBin() ? const_cast<TetraCL*>(GetRecycleBin()->FindTetra(v1,v2,v3)) : 0; }

// ===== Public interface =====

    const IdCL<VertexCL>& GetId           () const { return _Id; }                          ///< get id of this vertex (numbered locally on this proc)
#ifndef _PAR
    Uint                  GetLevel        () const { return _Level; }                       ///< get level of vertex (=first appearance in the multigrid)
#else
    Uint                  GetLevel        () const { return _dddH.attr;}                    ///< get level of vertex (=first appearance in the multigrid)
    // parallel functions
    static DDD_TYPE       GetType         ()       { return _dddT;}                         ///< get DDD-Vertex-Type
    DDD_PRIO              GetPrio         () const { return _dddH.prio;}                    ///< get priority
    DDD_GID               GetGID          () const { return _dddH.gid;}                     ///< get global id
    DDD_HDR               GetHdr          () const { return const_cast<DDD_HDR>(&_dddH);}   ///< get DDD-Hdr of this vertex
    bool                  IsMaster        () const { return GetPrio()>=PrioMaster; }        ///< check if vertex is master
    bool                  MayStoreUnk     () const { return GetPrio()==PrioHasUnk; }        ///< check for ability of storing unknowns due to priority
    bool                  IsLocal         () const { return DDD_InfoIsLocal(GetHdr());}     ///< check if vertex is local
    int                   GetNumDist      () const;                                         ///< get number of procs on which the vertex is stored
    bool                  IsExclusive     ( Priority prio=PrioMaster ) const;               ///< check if vertex is exclusive
    int *                 GetProcList     () const { return DDD_InfoProcList(GetHdr());}    ///< get list of procs and prios of this vertex
    void                  XferDelete      ()                                                ///< tell DDD vertex will be deleted
        { DDD_XferDeleteObj(&_dddH);
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
        }
    void                  SetPrio(DDD_PRIO p)      { _dddH.prio=p;}                         ///< set priority of this vertex (danger: no notification to DDD, use PrioChange!)
#endif
    const Point3DCL&      GetCoord        () const { return _Coord; }                       ///< get coordinate of this vertex
    bool                  IsOnBoundary    () const { return _BndVerts; }                    ///< check if this vertex lies on domain boundary
    const_BndVertIt       GetBndVertBegin () const { return _BndVerts->begin(); }
    const_BndVertIt       GetBndVertEnd   () const { return _BndVerts->end(); }
    bool                  HasBnd          (const BndPointCL& BndVert) const;                ///< check if Vertex belongs on a special boundary index given by a boundary vertex
    bool                  IsInTriang      (Uint TriLevel) const                             ///< check if the vertex can be found in a triangulation level
      { return  GetLevel() <= TriLevel; }

    // Debugging
    bool IsSane    (std::ostream&, const BoundaryCL& ) const;                               ///< check for sanity of this vertex
    void DebugInfo (std::ostream&) const;                                                   ///< get debug-information
    void SetNeeded (bool n) const {_needed=n;}                                              ///< Set if the simplex is needed
    bool GetNeeded () const {return _needed;}                                               ///< Get needed
};


/*******************************************************************
*   E D G E  C L                                                   *
*******************************************************************/
/// \brief Represents an edge in the multigrid
/** The refinement algorithm works by manipulating the marks '_MFR'
    of the edges. This is possible, because the refinement pattern
    of a face/tetrahedron is determined by the pattern of its edges.
    Marking the edges ensures consistency between the neighbors.  */
/*******************************************************************
*   E D G E  C L                                                   *
*******************************************************************/
class EdgeCL
{
  friend class PeriodicEdgesCL;
#ifdef _PAR
    friend class ParMultiGridCL;
#endif

  public:
    typedef MG_VertexContT::LevelCont VertContT;                                    ///< container for vertices
    typedef MG_EdgeContT::LevelCont   EdgeContT;                                    ///< container for subedges

  private:
    SArrayCL<VertexCL*, 2> _Vertices;                                               // "left" and "right" vertex of the edge
    VertexCL*              _MidVertex;                                              // midvertex, if the edge is refined
    SArrayCL<BndIdxT, 2>   _Bnd;                                                    // an edge can be found on (max) two boundary-segments
    mutable short int      _MFR;                                                    // mark, if the edge should be/is refined (set by refinement-algo)
    short int              _localMFR;                                               // MFR!=localMFR iff edge is on periodic boundary

    bool                   _RemoveMark;                                             // mark for removement
#ifndef _PAR
    Uint                   _Level : 8;                                              // level of the edge (according to owning tetras)
#else
    // parallel
    static DDD_TYPE        _dddT;                                                   // DDD-Type of edges
    DDD_HEADER             _dddH;                                                   // DDD-Header
    mutable short int      _AccMFR;                                                 // accumulated MFR over all procs
#if DROPSDebugC
    mutable bool subscribed_;
#endif

    EdgeCL() : _Vertices(static_cast<VertexCL*>(0)), _MidVertex(0), _Bnd(NoBndC),   // standard constructor used for DDD to create an edge
               _MFR(0), _RemoveMark(false), _AccMFR(0)
    {
        DDD_MarkHdrInvalid(&_dddH);                                                 // The header will be constructed by DDD
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
    }
    static void Declare();                                                          // declare DDD-Edge-Type (implemented in "parallel/parmultigrid.cpp")
    static void Define();                                                           // define DDD-Edge-Type (implemented in "parallel/parmultigrid.cpp")
#endif
    mutable bool _needed;                                                       // Check within IsSane, if the Simples is needed by tetra or is forgotten to delete

  public:
    UnknownHandleCL Unknowns;                                                   ///< access to unknowns on this edge

// ===== Interface for refinement =====
    ///< Create an edge
    inline EdgeCL (VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0= NoBndC, BndIdxT bnd1= NoBndC, short int MFR=0);
    EdgeCL (const EdgeCL&);                                                     ///< Danger!!! Copying simplices might corrupt the multigrid structure!!!
    ~EdgeCL() {
#ifdef _PAR
        // This exception may be thrown, by exiting the program. If the MultiGridCL is deleted, this can happen.
        Assert(!subscribed_, DROPSErrCL("Deleteting not unsubscribed Edge"), DebugSubscribeC);
#endif
    }
    // default dtor

    void AddBndIdx(BndIdxT idx)                                                 ///< add boundary-information
      { if (_Bnd[0]==NoBndC) _Bnd[0]= idx; else _Bnd[1]= idx; }

    // Midvertex
    VertexCL*   GetMidVertex   ()             { return _MidVertex; }            ///< get pointer to midvertex
    void        SetMidVertex   (VertexCL* vp) { _MidVertex= vp; }               ///< set pointer to midvertex
    void        RemoveMidVertex()             { _MidVertex= 0; }                ///< remove pointer to midvertex without deleting it
    void        BuildMidVertex (VertContT&, const BoundaryCL&);                 ///< create midvertex
    inline void BuildSubEdges  (EdgeContT&, VertContT&, const BoundaryCL&);     ///< build both subedges

    // Marks
#ifndef _PAR
    bool IsMarkedForRef       () const { return _MFR; }                         ///< check if this edge is marked for refinement
    void IncMarkForRef        ()       { ++_MFR; ++_localMFR;}                  ///< increase mark for refinement count
    void DecMarkForRef        ()       { --_MFR; --_localMFR;}                  ///< decrease mark for refinement count
    void ResetMarkForRef      ()       { _MFR= 0; _localMFR= 0; }               ///< remove mark for refinement
#else
    // parallel Marks (also accumulated marks over all procs are introduced)
    bool IsMarkedForRef       () const { return _AccMFR; }                      ///< check if this edge is marked for refinement
    void IncMarkForRef        () const { ++_MFR; ++_AccMFR; }                   ///< increase the local and accumulated mark for refinement count
    void DecMarkForRef        () const { --_MFR; --_AccMFR; }                   ///< decrease the local and accumulated mark for refinement count
    void ResetMarkForRef      ()       { _MFR= 0; _AccMFR= 0; }                 ///< remove the local and accumulated mark for refinement
    short int GetAccMFR       () const {return _AccMFR;}                        ///< get accumulated mark for refinement
    void SetAccMFR ( short int MFR)    { _AccMFR= MFR; }                        ///< Set accumulated MFR
#endif
    bool IsMarkedForRemovement() const { return _RemoveMark; }                  ///< check if edge is marked for removement
    void SetRemoveMark        ()       { _RemoveMark= true; }                   ///< set mark for removement
    void ClearRemoveMark      ()       { _RemoveMark= false; }                  ///< clear mark for removement

    // etc.
    void RecycleMe   () const { _Vertices[0]->Recycle(this); }                  ///< put a pointer to this edge into the recycle-bin of the "left" vertex
    void SortVertices()                                                         ///< sort vertices by id
      { if (_Vertices[1]->GetId() < _Vertices[0]->GetId()) std::swap(_Vertices[0],_Vertices[1]); }

// ===== Public interface =====
#ifndef _PAR
    Uint                  GetLevel        () const { return _Level; }                       ///< return level (stored within the class)
#else
    Uint                  GetLevel        () const { return _dddH.attr;}                    ///< return level (stored within the DDD-Header)
    static DDD_TYPE       GetType         ()       { return _dddT;}                         ///< return DDD-type of edges
    DDD_PRIO              GetPrio         () const { return _dddH.prio;}                    ///< return priority of this edge
    DDD_GID               GetGID          () const { return _dddH.gid;}                     ///< return global id of this edge
    DDD_HDR               GetHdr          () const { return const_cast<DDD_HDR>(&_dddH);}   ///< return DDD-Hdr of this edge
    bool                  IsMaster        () const { return GetPrio()>=PrioMaster; }        ///< check if edge is master
    bool                  MayStoreUnk     () const { return GetPrio()==PrioHasUnk; }        ///< check for ability of storing unknowns due to priority
    bool                  IsLocal         () const { return DDD_InfoIsLocal(GetHdr());}     ///< check if this edge is just stored on this proc
    int                   GetNumDist      () const;                                         ///< get number of procs on which the edge is stored
    int *                 GetProcList     () const { return DDD_InfoProcList(GetHdr());}    ///< get list of procs and prios of this edge
    bool                  IsExclusive     ( Priority prio=PrioMaster ) const;               ///< check if edge is exclusive
    void                  XferDelete      ()                                                ///< tell DDD that this edge will be deleted
    { DDD_XferDeleteObj(&_dddH);
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
        }
    void SetPrio                (DDD_PRIO p)       { _dddH.prio=p;}                         ///< set priority of this edge (no notification to DDD, use PrioChange instead!)

#endif
    const VertexCL* GetVertex     (Uint i)             const { return _Vertices[i]; }       ///< get pointer to the "left" or the "right" vertex
    const VertexCL* GetMidVertex  ()                   const { return _MidVertex; }         ///< get midvertex
    const VertexCL* GetNeighbor   (const VertexCL* vp) const { return vp==_Vertices[0] ? _Vertices[1] : _Vertices[0]; } ///< get opposite vertex
    bool            HasVertex     (const VertexCL* vp) const { return vp==_Vertices[0] || vp==_Vertices[1]; }           ///< check if the edge has a vertex
    bool            IsRefined     ()                   const { return _MidVertex; }         ///< check if edge is refined
    bool            IsOnBoundary  ()                   const { return _Bnd[0] != NoBndC; }  ///< check if edge lies on the domain boundary
    const BndIdxT*  GetBndIdxBegin()                   const { return _Bnd.begin(); }
    const BndIdxT*  GetBndIdxEnd  ()                   const
      { return IsOnBoundary() ? (_Bnd[1] == NoBndC ? _Bnd.begin()+1 : _Bnd.end() ) : _Bnd.begin(); }
    bool            IsInTriang    (Uint TriLevel)      const                                ///< check if edge can be found in a triangulation level
      { return GetLevel() == TriLevel || ( GetLevel() < TriLevel && !IsRefined() ); }
    short int GetMFR          () const {return _MFR;}                                       ///< get mark for refinement of this proc

    // Debugging
    bool IsSane    (std::ostream&) const;                                                   ///< check for sanity
    void DebugInfo (std::ostream&) const;                                                   ///< get debug-information
    void SetNeeded (bool n) const {_needed=n;}                                              ///< Set if the simplex is needed
    bool GetNeeded () const {return _needed;}                                               ///< Get needed
};


/*******************************************************************
*   F A C E  C L                                                   *
*******************************************************************/
/// \brief Represents a face in the multigrid
/** It can have neighbors on two levels: two regular ones (one per
    side) on level '_Level' and two green ones on the next.       */
/*******************************************************************
*   F A C E  C L                                                   *
*******************************************************************/
class FaceCL
{
#ifdef _PAR
    friend class ParMultiGridCL;
#endif

  private:
    SArrayCL<const TetraCL*,4> _Neighbors;                               // neighbor tetras of the face
    const BndIdxT              _Bnd;                                     // boundary-index of this face
    bool                       _RemoveMark;                              // mark for removement
#ifndef _PAR
    Uint                       _Level : 8;                               // level of the face (=level according to tetras) (in parallel stored in _dddH.attr)
#else
    static DDD_TYPE            _dddT;                                    // DDD-type of faces
    DDD_HEADER                 _dddH;                                    // DDD-Header of the distributed face
#if DROPSDebugC&DebugSubscribeC
    mutable bool subscribed_;
#endif
    idxtype                    _lbNoNeigh;                               // If this face has a neighbor tetra on a different proc, this number is stored here!

    FaceCL() : _Neighbors(static_cast<const TetraCL*>(0)), _Bnd(NoBndC), // standard constructor used for DDD to create an edge
               _RemoveMark(false), _lbNoNeigh(-1)
    {
        DDD_MarkHdrInvalid(&_dddH);                                      // the header will be constructed by DDD
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
    }
    static void Declare();                                               // declare DDD-type (implemented in "parallel/parmultigrid.cpp")
    static void Define();                                                // define DDD-type (implemented in "parallel/parmultigrid.cpp")
#endif
    mutable bool _needed;                                                // Check within IsSane, if the Simples is needed by tetra or is forgotten to delete

  public:
    UnknownHandleCL Unknowns;                                                      ///< access to unknowns on a face (not yet used)

// ===== Interface for refinement =====
    inline FaceCL (Uint Level, BndIdxT bnd= NoBndC);                               ///< create a face
    FaceCL (const FaceCL&);                                                        ///< Danger!!! Copying simplices might corrupt the multigrid structure!!!
    ~FaceCL() {
#ifdef _PAR
# if DROPSDebugC&DebugSubscribeC
        // This exception may be thrown, by exiting the program. If the MultiGridCL is deleted, this can happen.
        Assert(!subscribed_, DROPSErrCL("Deleteting not unsubscribed Face"), DebugSubscribeC);
# endif
#endif
    }
    // default dtor

    // RemovementMarks
    bool IsMarkedForRemovement() const { return _RemoveMark; }                     ///< check if marked for removement
    void SetRemoveMark        ()       { _RemoveMark= true; }                      ///< set mark for removement
    void ClearRemoveMark      ()       { _RemoveMark= false; }                     ///< clear mark for removement

    // Neighbors
    void LinkTetra  (const TetraCL*);                                              ///< link a tetra to face
    void UnlinkTetra(const TetraCL*);                                              ///< unlink a tetra of face
    void SetNeighbor(Uint i, TetraCL* tp) { _Neighbors[i]=tp; }                    ///< Set tetra as neighbor
#ifdef _PAR
    bool IsLinkedTo (const TetraCL*) const;                                        ///< check if a tetra is linked to this face
#endif

    // Recycling
    void RecycleMe(VertexCL* vp0, const VertexCL* vp1, const VertexCL* vp2) const  ///< put a pointer to this face into the recycle-bin of the first vertex
      { vp0->Recycle(this,vp1,vp2); }

// ===== Public Interface =====
#ifndef _PAR
    Uint            GetLevel       () const { return _Level; }                          ///< get level of the face (stored within the class)
#else
    Uint            GetLevel       () const { return _dddH.attr;}                       ///< get level of the face (stored within the DDD-Header)
    static DDD_TYPE GetType        ()       { return _dddT;}                            ///< get DDD-type of the faces
    DDD_PRIO        GetPrio        () const { return _dddH.prio;}                       ///< get priority of the face
    DDD_GID         GetGID         () const { return _dddH.gid;}                        ///< get global id of the face
    DDD_HDR         GetHdr         () const { return const_cast<DDD_HDR>(&_dddH);}      ///< get DDD-Hdr of the face
    bool            IsGhost        () const { return GetPrio()<PrioMaster; }            ///< check if face is ghost
    bool            IsMaster       () const { return GetPrio()>=PrioMaster; }           ///< check if face is master
    bool            MayStoreUnk    () const { return GetPrio()==PrioHasUnk; }           ///< check for ability of storing unknowns due to priority
    bool            IsLocal        () const { return DDD_InfoIsLocal(GetHdr());}        ///< check if this face is just stored on this proc
    int             GetNumDist     () const;                                            ///< get number of procs on which the face is stored
    int *           GetProcList    () const { return DDD_InfoProcList(GetHdr());}       ///< get list of procs and prios of this face
    bool            IsOnProcBnd    () const;                                            ///< check if face lies between two procs
    DDD_PROC        GetNeighborProc() const                                             ///< get neighbor proc over the face
        {return*(DDD_InfoProcList(const_cast<DDD_HDR>(&_dddH))+2);}
    void            XferDelete     ()                                                   ///< tell DDD that this face will be deleted
        {DDD_XferDeleteObj(&_dddH);
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
        }
    void            SetPrio  (DDD_PRIO p)   { _dddH.prio=p;}                            ///< set priority of this face (DANGER: no notification to DDD, use PrioChange instead!)
    idxtype         GetLbNeigh     () const { return _lbNoNeigh;}                       ///< get loadbalance-number of neighbor tetra, if this tetra is stored on a different proc
    void            SetLbNeigh (idxtype no) { _lbNoNeigh=no;}                           ///< set loadbalance-number
#endif

    bool        IsOnNextLevel() const { return _Neighbors[2] || _Neighbors[3]; }        ///< check if face can be found in the next level
    bool        IsOnBoundary () const { return _Bnd != NoBndC; }                        ///< check if face lies on the domain-boundary
    BndIdxT     GetBndIdx    () const { return _Bnd; }                                  ///< get index of the boundary-segment
    inline bool IsRefined    () const;                                                  ///< check if face is refined
    bool        IsInTriang   (Uint TriLevel) const                                      ///< check if face can be found in a triangulation level
      { return GetLevel() == TriLevel || ( GetLevel() < TriLevel && !IsRefined() ); }

    // Get simplex
    inline const VertexCL* GetVertex(Uint)       const;                                 ///< get i'th vertex of the face
    inline const EdgeCL*   GetEdge  (Uint)       const;                                 ///< get i'th edge of the face
    inline const TetraCL*  GetTetra (Uint, Uint) const;                                 ///< get tetra of level and number

    // Neighboring tetras
           const TetraCL* GetSomeTetra     () const { return _Neighbors[0]; }           ///< return pointer to first neighbor
    inline bool           HasNeighborTetra (const TetraCL*)       const;                ///< check if a tetra is neighbor
    inline const TetraCL* GetNeighborTetra (const TetraCL*)       const;                ///< get neighbor tetra of another tetra
           const TetraCL* GetNeighInTriang (const TetraCL*, Uint) const;                ///< get neighbor tetra in triangulation level
    inline Uint           GetFaceNumInTetra(const TetraCL*)       const;                ///< get number of face within a tetra
           const TetraCL* GetNeighbor      (Uint i) const { return _Neighbors[i];}      ///< get raw tetra-pointer from the array
    // Debugging
    bool IsSane   (std::ostream&) const;                                                ///< check for sanity
    void DebugInfo(std::ostream&) const;                                                ///< get debug-information
    void SetNeeded (bool n) const {_needed=n;}                                          ///< Set if the simplex is needed
    bool GetNeeded () const {return _needed;}                                           ///< Get needed
};


/*******************************************************************
*   T E T R A  C L                                                 *
*******************************************************************/
/// \brief Represents a tetrahedron in multigrid
/** This is probably the most important data structure of DROPS.
    All major routines that work on the grid (i.e. the refinement
    algorithm, the routines to set up a discretized system, the
    error estimator) "do it tetra by tetra".
    \todo (merge) _RefRule and _RefMark as single Usint OK?       */
/*******************************************************************
*   T E T R A  C L                                                 *
*******************************************************************/
class TetraCL
{
    friend class MultiGridCL;
#ifdef _PAR
    friend class ParMultiGridCL;
#endif

  public:
    typedef SArrayCL<VertexCL*,NumVertsC>::iterator         VertexPIterator;            ///< iterator of pointers to vertices of this tetra
    typedef SArrayCL<VertexCL*,NumVertsC>::const_iterator   const_VertexPIterator;      ///< const version
    typedef SArrayCL<EdgeCL*,NumEdgesC>::iterator           EdgePIterator;              ///< iterator of pointers to edges of this tetra
    typedef SArrayCL<EdgeCL*,NumEdgesC>::const_iterator     const_EdgePIterator;        ///< const version
    typedef SArrayCL<FaceCL*,NumFacesC>::iterator           FacePIterator;              ///< iterator of pointers to faces of this tetra
    typedef SArrayCL<FaceCL*,NumFacesC>::const_iterator     const_FacePIterator;        ///< const version
    typedef SArrayCL<TetraCL*,MaxChildrenC>::iterator       ChildPIterator;             ///< iterator of pointers to children of this tetra
    typedef SArrayCL<TetraCL*,MaxChildrenC>::const_iterator const_ChildPIterator;       ///< const version
    typedef MG_VertexContT::LevelCont                       VertContT;                  ///< container for verts for linking purpose
    typedef MG_EdgeContT::LevelCont                         EdgeContT;                  ///< container for edges for linking purpose
    typedef MG_FaceContT::LevelCont                         FaceContT;                  ///< container for faces for linking purpose
    typedef MG_TetraContT::LevelCont                        TetraContT;                 ///< container for children for linking purpose

  private:
    // static arrays for computations
    static SArrayCL<EdgeCL*, NumAllEdgesC> _ePtrs;                      // EdgePointers for linking edges within refinement
    static SArrayCL<FaceCL*, NumAllFacesC> _fPtrs;                      // FacePointers for linking faces within refinement

    IdCL<TetraCL> _Id;                                                  // id-number (locally numbered on one proc)
    Usint          _RefRule;                                            // actual refinement of the tetrahedron
    mutable Usint  _RefMark;                                           // refinement-mark (e.g. set by the error estimator)
#ifndef _PAR
    Uint          _Level : 8;
#else
    idxtype                  _lbNr;                                     // number of this Tetra used for LoadBalance with parmetis (typedef in parmetis.h), initialized with -1 as no number
    static DDD_TYPE          _dddT;                                     // DDD-Type for Tetras, for every Tetra the same
    DDD_HEADER               _dddH;                                     // DDD-Header
# if DROPSDebugC&DebugSubscribeC
    mutable bool subscribed_;
# endif
#endif

    // subsimplices, parent, children
    SArrayCL<VertexCL*,NumVertsC>    _Vertices;                         // container for verts of tetra
    SArrayCL<EdgeCL*,NumEdgesC>      _Edges;                            // container for edges of tetra
    SArrayCL<FaceCL*,NumFacesC>      _Faces;                            // container for faces of tetra
    TetraCL*                         _Parent;                           // container for parent of tetra
    SArrayCL<TetraCL*,MaxChildrenC>* _Children;                         // container for children of tetra, for leaves: null-pointer
// ??? TODO: Kann man das ohne Speicherluecken und ohne Fragmentieren hinkriegen ???

#ifdef _PAR
    TetraCL() : _Id(), _RefRule(UnRefRuleC),_RefMark(NoRefMarkC),       // standard constructor used by DDD to create an edge
                _lbNr(-1), _Vertices(static_cast<VertexCL*>(0)),
                _Edges(static_cast<EdgeCL*>(0)),_Faces(static_cast<FaceCL*>(0)),
                _Parent(0), _Children(0)
    {
        DDD_MarkHdrInvalid(&_dddH);                                     // the header will be constructed by DDD
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
    }

    static void Declare();                                              // declare Tetra type (implemented in "parallel/parmultigrid.cpp")
    static void Define();                                               // define Tetra type (implemented in "parallel/parmultigrid.cpp")
#endif

  public:
    UnknownHandleCL Unknowns;                                                   ///< access to unknowns on tetras

// ===== Interface for refinement =====
    inline  TetraCL (VertexCL*, VertexCL*, VertexCL*, VertexCL*, TetraCL*, IdCL<TetraCL> id= IdCL<TetraCL>());     ///< constructor of verts and parent; FileBuilderCL has to construct the _Id, too, thus it can optionally be set.
#ifdef _PAR
    inline  TetraCL (VertexCL*, VertexCL*, VertexCL*, VertexCL*, TetraCL*, Uint, IdCL<TetraCL> id= IdCL<TetraCL>());     ///< constructor of verts and level, if no parent is available; FileBuilderCL has to construct the _Id, too, thus it can optionally be set.
#endif
    TetraCL (const TetraCL&);                                                   ///< Danger!!! Copying simplices might corrupt the multigrid structure!!!
    inline ~TetraCL ();

    // access to children, vertices
    ChildPIterator GetChildBegin  ()                                            ///< "Pointer-Iterator" to first child
      { return _Children ? _Children->begin() : 0; }
    ChildPIterator GetChildEnd    ()                                            ///< "Pointer-Iterator" to end of children
      { return _Children ? _Children->begin() + GetRefData().ChildNum : 0; }
    VertexCL*      GetVertMidVert (Uint i)                                      ///< return pointer to midvertex of edge or vertex of tetra
      { return IsMidVert(i) ? _Edges[EdgeOfMidVert(i)]->GetMidVertex() : _Vertices[i]; }

    // rules and marks
    void        SetRefRule          (Uint RefRule) { _RefRule = RefRule; }      ///< set refinement rule for this tetra
    void        RestrictMark        ();                                         ///< set regular refinement if marked and manipulate MFR on edges
    inline void CommitRegRefMark    () const;                                   ///< increase MFR on edges
    inline void UnCommitRegRefMark  () const;                                   ///< decrease MFR on edges
    inline void Close               ();                                         ///< calculate green closure
    void        ClearAllRemoveMarks ();                                         ///< save all subsimplices of this tetra and all child tetras from removement

    // recycling
    void RecycleMe        () const { _Vertices[0]->Recycle(this); }             ///< put pointer to the tetra into recycle-bin of vertex(0)
    void RecycleReusables ();                                                   ///< put all subsimplices of the tetra and children that are needed in next refinement step into recycle-bins

    // building children
    void        CollectEdges           (const RefRuleCL&, VertContT&, EdgeContT&, const BoundaryCL&);   ///< build or unrecycle edges that are needed for refinement
    void        CollectFaces           (const RefRuleCL&, FaceContT&);                                  ///< build or unrecycle faces that are needed for refinement
    inline void LinkEdges              (const ChildDataCL&);                                            ///< link edges from "_ePtrs" to the tetra according to the ChildDataCL
    inline void LinkFaces              (const ChildDataCL&);                                            ///< link faces from "_fPtrs" to the tetra according to the ChildDataCL
    void CollectAndLinkChildren (const RefRuleCL&, TetraContT&);                                 ///< build, unrecycle and link children
    void        UnlinkFromFaces        ()                                                               ///< remove link from faces to the tetra
      { for (Uint face=0; face<NumFacesC; ++face) _Faces[face]->UnlinkTetra(this); }

    // used by builder
    void BuildEdges        (EdgeContT&);                                         ///< build edges
    void BuildAndLinkFaces (FaceContT&);                                         ///< build and link faces
    void SetFace           (Uint f, FaceCL* fp) { _Faces[f]= fp; }               ///< set face
    void SetEdge           (Uint e, EdgeCL* ep) { _Edges[e]= ep; }               ///< set edge
    void SetRefMark        (Uint refmark)       { _RefMark= refmark; }           ///< set RefMark
    void SetChild          (Uint, TetraCL*);                                     ///< set a child-pointer; if neccessary the _Children-Array is allocated first

//
// Public Interface
//
#ifndef _PAR
    Uint GetLevel() const { return _Level; }                                     ///< return level of tetra
#else
    static DDD_TYPE GetType() { return _dddT;}
    bool IsProcBnd (Uint face) const { return _Faces[face]->IsOnProcBnd(); }     ///< check if face of tetra is on processor-boundary

    Uint    GetLevel   () const { return _dddH.attr;}                            ///< return level of tetra (stored within DDD-Header)
    DDD_PRIO GetPrio    () const { return _dddH.prio;}                            ///< get priority of the tetra
    DDD_GID GetGID     () const { return _dddH.gid;}                             ///< get global id of the tetra
    DDD_HDR GetHdr     () const { return const_cast<DDD_HDR>( &_dddH); }         ///< get Hdr of the tetra (const)
    DDD_HDR GetHdr     ()       { return &_dddH; }                               ///< get Hdr of the tetra
    idxtype GetLbNr    () const { return _lbNr;}                                 ///< get number for load balance
    bool    IsGhost    () const { return GetPrio()<PrioMaster; }                 ///< check if tetra is ghost
    bool    IsMaster   () const { return GetPrio()>=PrioMaster; }                ///< check if tetra is master
    bool    MayStoreUnk() const { return GetPrio()==PrioMaster; }                ///< check for ability of storing unknowns due to priority
    bool    HasGhost   () const;                                                 ///< check if tetra has a ghost-copy somewhere
    bool    HasLbNr    () const { return _lbNr>=0;}                              ///< check if this tetra got a number for load balance
    bool    IsLocal    () const { return DDD_InfoIsLocal(GetHdr());}             ///< check if this tetra is just stored on this proc
    int     GetNumDist () const;                                                 ///< get number of procs on which the tetra is stored
    int *   GetProcList() const { return DDD_InfoProcList(GetHdr());}            ///< get list of procs and prios of this tetra
    bool    IsExclusive( Priority prio=PrioMaster ) const;                       ///< check if tetra is exclusive

    void SetPrio(DDD_PRIO p)     { _dddH.prio=p;}                                ///< set priority of this face (no notification to DDD, use PrioChange instead!)
    void SetLbNr(idxtype nr)     { _lbNr=nr;}                                    ///< set the number for load balance
    void DelLbNr()               { _lbNr=-1;}                                    ///< clear the number for load balance

    void XferDelete()                                                            ///< tell DDD that this tetra will be deleted
        { DDD_XferDeleteObj(&_dddH);
#if DROPSDebugC&DebugSubscribeC
        subscribed_=false;
#endif
        }
#endif

    const IdCL<TetraCL>& GetId () const { return _Id; }                          ///< get local id
    Uint GetRefMark            () const { return _RefMark; }                     ///< get refinement mark
    Uint GetRefRule            () const { return _RefRule; }                     ///< get refinement rule
    inline const RefRuleCL& GetRefData () const;                                 ///< get information about refinement data

// _RefMark is mutable
    void SetRegRefMark () const { _RefMark= RegRefMarkC; }                       ///< mark tetra for regular refinement
    void SetRemoveMark () const { _RefMark= RemoveMarkC; }                       ///< mark tetra for removement
    void SetNoRefMark  () const { _RefMark= NoRefMarkC; }                        ///< mark tetra for no refinement

    bool IsMarkEqRule   () const { return _RefMark == _RefRule; }                ///< check if tetra is refined as the mark says
    bool IsUnrefined    () const { return _RefRule == UnRefRuleC; }              ///< check if tetra is unrefined
    bool IsRegularlyRef () const { return _RefRule == RegRefRuleC; }             ///< check if tetra is regular refined
    bool IsRegular      () const                                                 ///< check if the tetra is regular
#ifndef _PAR
      { return GetLevel()!=0 ? _Parent->GetRefRule() == RegRefRuleC : true; }
#else
      { return GetLevel()!=0 ? (IsGhost() ? true : _Parent->GetRefRule() == RegRefRuleC ) : true; }
#endif

    bool IsMarkedForRef        () const                                          ///< check if tetra is marked for refinement
      { return _RefMark != NoRefMarkC && _RefMark != RemoveMarkC; }
    bool IsMarkedForRegRef     () const { return _RefMark == RegRefMarkC; }      ///< check if tetra is marked for regular refinement
    bool IsMarkedForRemovement () const { return _RefMark == RemoveMarkC; }      ///< check if tetra is marked for removement
    bool IsMarkedForNoRef      () const { return _RefMark == NoRefMarkC; }       ///< check if tetra is marked for no refinement

    /// \name access to subsimplices
    //@{
    const_VertexPIterator GetVertBegin ()   const { return _Vertices.begin(); }
    const_VertexPIterator GetVertEnd   ()   const { return _Vertices.end(); }
    const_EdgePIterator   GetEdgesBegin()   const { return _Edges.begin(); }
    const_EdgePIterator   GetEdgesEnd  ()   const { return _Edges.end(); }
    const_FacePIterator   GetFacesBegin()   const { return _Faces.begin(); }
    const_FacePIterator   GetFacesEnd  ()   const { return _Faces.end(); }
    const VertexCL*       GetVertex(Uint i) const { return _Vertices[i]; }
    const EdgeCL*         GetEdge  (Uint i) const { return _Edges[i]; }
    const FaceCL*         GetFace  (Uint i) const { return _Faces[i]; }
    //@}

    bool           IsBndSeg        (Uint face) const { return _Faces[face]->IsOnBoundary(); }           ///< check if face lies on domain-boundary
    bool           IsNeighbor      (Uint face) const { return _Faces[face]->HasNeighborTetra(this); }   ///< check if this tetra is neigbor to a face
    BndIdxT        GetBndIdx       (Uint face) const { return _Faces[face]->GetBndIdx(); }              ///< get boundary index of a face
    const TetraCL* GetNeighbor     (Uint face) const { return _Faces[face]->GetNeighborTetra(this); }   ///< get pointer to neighbor tetra over face
    const TetraCL* GetNeighInTriang(Uint face, Uint trilevel) const
      { return _Faces[face]->GetNeighInTriang( this, trilevel); }

    /// \name access to parent and children
    //@{
    const TetraCL*       GetParent     ()       const { return _Parent; }
    const_ChildPIterator GetChildBegin ()       const { return _Children ? _Children->begin() : 0; }
    const_ChildPIterator GetChildEnd   ()       const { return _Children ? _Children->begin() + GetRefData().ChildNum : 0; }
    const TetraCL*       GetChild      (Uint i) const { return (*_Children)[i]; }
    //@}

    double               GetVolume     () const;                                           ///< get volume of tetra
    double               GetNormal     (Uint face, Point3DCL& normal, double& dir) const;  ///< get normal onto face with direction
    double               GetOuterNormal(Uint face, Point3DCL& normal)              const;  ///< get outer normal onto face
    bool                 IsInTriang    (Uint TriLevel) const                               ///< check if tetra is in triangulation level
      { return GetLevel() == TriLevel || ( GetLevel() < TriLevel && IsUnrefined() ); }

    bool IsSane    (std::ostream&) const;                                                  ///< check for sanity
    void DebugInfo (std::ostream&) const;                                                  ///< get debug-information
};


/// \name barycentic center of simplices
//@{
inline
Point3DCL GetBaryCenter(const VertexCL& v) { return v.GetCoord(); }
Point3DCL GetBaryCenter(const EdgeCL&);
Point3DCL GetBaryCenter(const FaceCL&);
Point3DCL GetBaryCenter(const TetraCL&);
Point3DCL GetBaryCenter(const TetraCL& t, Uint face);
//@}

Point3DCL GetWorldCoord(const TetraCL&, const SVectorCL<3>&);
Point3DCL GetWorldCoord(const TetraCL&, Uint face, const SVectorCL<2>&);
// barycentric coordinates:
Point3DCL GetWorldCoord(const TetraCL&, const SVectorCL<4>&);

SVectorCL<3> FaceToTetraCoord(const TetraCL& t, Uint f, SVectorCL<2> c);

/// \brief Maps world-coordinates p to barycentric coordinates of the tetra t.
class World2BaryCoordCL
{
  private:
    QRDecompCL<4> qr_;

  public:
    World2BaryCoordCL (const TetraCL& t);
    BaryCoordCL operator() (const Point3DCL& p) const;
};



//**************************************************************************
//  I n l i n e   f u n c t i o n s  of the classes above                  *
//**************************************************************************


// ********** RecycleBinCL **********

inline const EdgeCL* RecycleBinCL::FindEdge (const VertexCL* v) const
{
    for (EdgeContT::const_iterator it= _Edges.begin(), end= _Edges.end(); it!=end; ++it)
        if ( (*it)->GetVertex(1) == v ) return *it;
    return 0;
}


inline const FaceCL* RecycleBinCL::FindFace (const VertexCL* v1, const VertexCL* v2) const
{
    for(FaceWrapperContT::const_iterator it= _Faces.begin(), end= _Faces.end(); it!=end; ++it)
        if ( it->vert1 == v1 && it->vert2 == v2 ) return it->face;
    return 0;
}


inline const TetraCL* RecycleBinCL::FindTetra (const VertexCL* v1, const VertexCL* v2, const VertexCL* v3) const
{
    for(TetraContT::const_iterator it= _Tetras.begin(), end= _Tetras.end(); it!=end; ++it)
        if ( (*it)->GetVertex(1) == v1 && (*it)->GetVertex(2) == v2 && (*it)->GetVertex(3) == v3 )
            return *it;
    return 0;
}


// ********** VertexCL **********

// Constructor for the VertexCL. The differences between parallel an serial mode are:
// the level is stored at different places and in the parallel mode, the Header of the vertex
// has to be created!

#ifndef _PAR
inline VertexCL::VertexCL (const Point3DCL& Coord, Uint FirstLevel, IdCL<VertexCL> id)
    : _Id(id), _Coord(Coord), _BndVerts(0), _Bin(0),_RemoveMark(false), _Level(FirstLevel)
{}
#else
inline VertexCL::VertexCL (const Point3DCL& Coord, Uint FirstLevel, IdCL<VertexCL> id)
    : _Id(id), _Coord(Coord), _BndVerts(0), _Bin(0),_RemoveMark(false)/*, ProcSysNum(0)*/
{
    DDD_HdrConstructor(&_dddH, _dddT, PrioMaster, FirstLevel);
#if DROPSDebugC&DebugSubscribeC
    subscribed_=true;
#endif
}
#endif

#ifndef _PAR
inline VertexCL::VertexCL (const VertexCL& v)
    : _Id(v._Id), _Coord(v._Coord),
      _BndVerts(v._BndVerts ? new std::vector<BndPointCL>(*v._BndVerts) : 0),
      _Bin(v._Bin ? new RecycleBinCL(*v._Bin) : 0),_RemoveMark(v._RemoveMark),
      _Level (v._Level), Unknowns(v.Unknowns)
{}
#else
inline VertexCL::VertexCL (const VertexCL& v)
    : _Id(v._Id), _Coord(v._Coord),
      _BndVerts(v._BndVerts ? new std::vector<BndPointCL>(*v._BndVerts) : 0),
      _Bin(v._Bin ? new RecycleBinCL(*v._Bin) : 0),
      _RemoveMark(v._RemoveMark), Unknowns(v.Unknowns)/*, ProcSysNum(0)*/
{
    DDD_HdrConstructorMove(&_dddH, const_cast<DDD_HEADER*>( &v._dddH));     // The Header of the Objekt has to be moved
#if DROPSDebugC&DebugSubscribeC
    v.subscribed_=false;
    subscribed_=true;
#endif
}
#endif


inline VertexCL::~VertexCL ()
{
  delete _BndVerts;
  DestroyRecycleBin();
#ifdef _PAR
# if DROPSDebugC&DebugSubscribeC
  // This exception may be thrown, by exiting the program. If the MultiGridCL is deleted, this can happen.
  Assert(!subscribed_, DROPSErrCL("Deleteting not unsubscribed Vertex"), DebugSubscribeC);
# endif
#endif
}

inline void VertexCL::AddBnd (const BndPointCL& BndVert)
{
    if (_BndVerts)
        _BndVerts->push_back(BndVert);
    else
        _BndVerts= new std::vector<BndPointCL>(1,BndVert);
}


inline void VertexCL::BndSort ()
    { std::sort(_BndVerts->begin(), _BndVerts->end(), BndPointSegLessCL()); }

#ifdef _PAR
inline bool VertexCL::IsExclusive( Priority prio) const
/** A vertex is exclusive, if the vertex is local stored or the vertex is master and owned by
    by the proc with the smallest id*/
{
    if (IsLocal())
        return true;
    int *procs = GetProcList();
    int minproc=ProcCL::Size();
    const int me = ProcCL::MyRank();
    int i=0;
    while (procs[i]!=-1)
    {
        if (procs[i]<minproc && procs[i+1]>=prio)
            minproc = procs[i];
        i+=2;
    }
    Assert(minproc != ProcCL::Size(), DROPSErrCL("VertexCL::IsExclusive: no exclusive proc for vertex!"), DebugParallelC);
    if (me==minproc)
        return true;
    else
        return false;
}
#endif

inline bool VertexCL::HasBnd(const BndPointCL& BndVert) const
{
    if (!IsOnBoundary()) return false;
        for (const_BndVertIt it(GetBndVertBegin()), end(GetBndVertEnd()); it!=end; ++it)
            if (it->GetBndIdx()==BndVert.GetBndIdx()) return true;
    return false;
}


// ********** EdgeCL **********

#ifndef _PAR
inline EdgeCL::EdgeCL (VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0, BndIdxT bnd1, short int MFR)
    : _MidVertex(0), _MFR(MFR), _localMFR(MFR), _RemoveMark(false), _Level(Level)
{
    _Vertices[0]= vp0; _Vertices[1]= vp1;
    _Bnd[0]= bnd0; _Bnd[1]= bnd1;
}
#else
inline EdgeCL::EdgeCL (VertexCL* vp0, VertexCL* vp1, Uint Level, BndIdxT bnd0, BndIdxT bnd1, short int MFR)
    : _MidVertex(0), _MFR(MFR), _localMFR(MFR), _RemoveMark(false), _AccMFR(MFR)/*, ProcSysNum(0)*/
{
    _Vertices[0]= vp0; _Vertices[1]= vp1;
    _Bnd[0]= bnd0; _Bnd[1]= bnd1;
    DDD_HdrConstructor(&_dddH, _dddT, PrioMaster, Level);
#if DROPSDebugC&DebugSubscribeC
    subscribed_=true;
#endif
}
#endif

#ifndef _PAR
inline EdgeCL::EdgeCL (const EdgeCL& e)
    : _Vertices(e._Vertices), _MidVertex(e._MidVertex), _Bnd(e._Bnd),
      _MFR(e._MFR), _localMFR(e._localMFR), _RemoveMark(e._RemoveMark), _Level(e._Level),
      Unknowns(e.Unknowns)
{}
#else
inline EdgeCL::EdgeCL (const EdgeCL& e)
    : _Vertices(e._Vertices), _MidVertex(e._MidVertex), _Bnd(e._Bnd),
      _MFR(e._MFR), _localMFR(e._localMFR), _RemoveMark(e._RemoveMark), _AccMFR(e._AccMFR),
      Unknowns(e.Unknowns)/*, ProcSysNum(0)*/
{
    DDD_HdrConstructorMove(&_dddH, const_cast<DDD_HEADER*>( &e._dddH));
#if DROPSDebugC&DebugSubscribeC
    e.subscribed_=false;
    subscribed_=true;
#endif
}
#endif


inline void EdgeCL::BuildSubEdges(EdgeContT& edgecont, VertContT& vertcont, const BoundaryCL& Bnd)
{
    BuildMidVertex(vertcont, Bnd);
    edgecont.push_back( EdgeCL(_Vertices[0], GetMidVertex(), GetLevel()+1, _Bnd[0], _Bnd[1]) );
    edgecont.push_back( EdgeCL(GetMidVertex(), _Vertices[1], GetLevel()+1, _Bnd[0], _Bnd[1]) );
}

#ifdef _PAR
inline bool EdgeCL::IsExclusive( Priority prio) const
/** An edge is exclusive, if the edge is local stored or the edge is master and owned by
    by the proc with the smallest id*/
{
    if (IsLocal())
        return true;
    int *procs = GetProcList();
    int minproc=ProcCL::Size();
    const int me = ProcCL::MyRank();
    int i=0;
    while (procs[i]!=-1)
    {
        if (procs[i]<minproc && procs[i+1]>=prio)
            minproc = procs[i];
        i+=2;
    }
    Assert(minproc != ProcCL::Size(), DROPSErrCL("EdegCL::IsExclusive: no exclusive proc for edge!"), DebugParallelC);
    if (me==minproc)
        return true;
    else
        return false;
}
#endif

// ********** FaceCL **********
#ifndef _PAR
inline FaceCL::FaceCL (Uint Level, BndIdxT bnd)
    : _Bnd(bnd), _RemoveMark(false), _Level(Level) {}
#else
inline FaceCL::FaceCL (Uint Level, BndIdxT bnd)
    : _Bnd(bnd), _RemoveMark(false), _lbNoNeigh(-1)
{
    DDD_HdrConstructor(&_dddH, _dddT, PrioMaster, Level);
#if DROPSDebugC&DebugSubscribeC
    subscribed_=true;
#endif
}
#endif

#ifndef _PAR
inline FaceCL::FaceCL (const FaceCL& f)
    : _Neighbors(f._Neighbors), _Bnd(f._Bnd),_RemoveMark(f._RemoveMark),
       _Level(f._Level), Unknowns(f.Unknowns) {}
#else
inline FaceCL::FaceCL (const FaceCL& f)
    : _Neighbors(f._Neighbors), _Bnd(f._Bnd),
      _RemoveMark(f._RemoveMark), _lbNoNeigh(-1), Unknowns(f.Unknowns)
{
    DDD_HdrConstructorMove(&_dddH, const_cast<DDD_HEADER*>( &f._dddH));
#if DROPSDebugC&DebugSubscribeC
    f.subscribed_=false;
    subscribed_=true;
#endif
}
#endif


inline bool FaceCL::IsRefined () const
{
    const TetraCL* const tp= GetSomeTetra();
    return RefinesFace( tp->GetRefRule(), GetFaceNumInTetra(tp) );
}


inline bool FaceCL::HasNeighborTetra (const TetraCL* tp) const
{
    if ( IsOnBoundary() ) return false;
    if ( tp->GetLevel() == GetLevel() ) // sequential version: always true, parallel: neighbor may be stored on a different proc
        return _Neighbors[0]==tp ? _Neighbors[1] : _Neighbors[0];
    Assert( IsOnNextLevel() && tp->GetLevel() == GetLevel()+1,
        DROPSErrCL("FaceCL::HasNeighborTetra: Illegal Level."), DebugRefineEasyC);
    return _Neighbors[2]==tp ? _Neighbors[3] : _Neighbors[2];
}

inline const TetraCL* FaceCL::GetNeighborTetra (const TetraCL* tp) const
{
    if ( tp->GetLevel()==GetLevel() )
        return _Neighbors[0]==tp ? _Neighbors[1] : _Neighbors[0];
    Assert( IsOnNextLevel() && tp->GetLevel() == GetLevel()+1,
        DROPSErrCL("FaceCL::GetNeighborTetra: Illegal Level."), DebugRefineEasyC);
    return _Neighbors[2]==tp ? _Neighbors[3] : _Neighbors[2];
}

inline Uint FaceCL::GetFaceNumInTetra (const TetraCL* tp) const
{
    for (Uint face=0; face<NumFacesC; ++face)
        if (tp->GetFace(face) == this) return face;
    throw DROPSErrCL("FaceCL::GetFaceNumInTetra: I'm not face of given tetra!");
}

inline const VertexCL* FaceCL::GetVertex (Uint v) const
{
    const TetraCL* const tp= GetSomeTetra();
    return tp->GetVertex( VertOfFace(GetFaceNumInTetra(tp), v) );
}

inline const EdgeCL* FaceCL::GetEdge (Uint e) const
{
    const TetraCL* const tp= GetSomeTetra();
    return tp->GetEdge( EdgeOfFace(GetFaceNumInTetra(tp), e) );
}

inline const TetraCL* FaceCL::GetTetra (Uint Level, Uint side) const
{
    if ( Level == GetLevel() ) return _Neighbors[side];
    else
    {
        Assert(IsOnNextLevel(), DROPSErrCL("FaceCL::GetTetra: No such level."), DebugRefineEasyC);
        return _Neighbors[side+2];
    }
}

// T e t r a C L

#ifndef _PAR
inline TetraCL::TetraCL (VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, IdCL<TetraCL> id)
    : _Id(id), _RefRule(UnRefRuleC), _RefMark(NoRefMarkC),
      _Level(Parent==0 ? 0 : Parent->GetLevel()+1),
      _Parent(Parent), _Children(0)
{
    _Vertices[0] = vp0; _Vertices[1] = vp1;
    _Vertices[2] = vp2; _Vertices[3] = vp3;
}
#else
inline TetraCL::TetraCL (VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, IdCL<TetraCL> id)
    : _Id(id), _RefRule(UnRefRuleC), _RefMark(NoRefMarkC), _lbNr(-1),
      _Parent(Parent), _Children(0)/*, ProcSysNum(0)*/
{
    _Vertices[0] = vp0; _Vertices[1] = vp1;
    _Vertices[2] = vp2; _Vertices[3] = vp3;
    DDD_HdrConstructor(&_dddH, _dddT, PrioMaster, Parent==0 ? 0 : Parent->GetLevel()+1);
#if DROPSDebugC&DebugSubscribeC
    subscribed_=true;
#endif
}

inline TetraCL::TetraCL (VertexCL* vp0, VertexCL* vp1, VertexCL* vp2, VertexCL* vp3, TetraCL* Parent, Uint lvl, IdCL<TetraCL> id)
    : _Id(id), _RefRule(UnRefRuleC), _RefMark(NoRefMarkC), _lbNr(-1),
      _Parent( Parent), _Children(0)/*, ProcSysNum(0)*/
{
    Assert(!Parent && Parent->GetLevel()!=lvl-1, DROPSErrCL("TetraCL::TetraCL: Parent and given level does not match"), DebugRefineEasyC);
    _Vertices[0] = vp0; _Vertices[1] = vp1;
    _Vertices[2] = vp2; _Vertices[3] = vp3;
    DDD_HdrConstructor(&_dddH, _dddT, PrioMaster, lvl);
#if DROPSDebugC&DebugSubscribeC
    subscribed_=true;
#endif
}

#endif


#ifndef _PAR
inline TetraCL::TetraCL (const TetraCL& T)
    : _Id(T._Id), _RefRule(T._RefRule), _RefMark(T._RefMark),
       _Level(T._Level), _Vertices(T._Vertices), _Edges(T._Edges),
      _Faces(T._Faces), _Parent(T._Parent),
      _Children(T._Children ? new SArrayCL<TetraCL*,MaxChildrenC> (*T._Children) : 0),
      Unknowns(T.Unknowns) {}
#else
inline TetraCL::TetraCL (const TetraCL& T)
    : _Id(T._Id), _RefRule(T._RefRule),
      _RefMark(T._RefMark), _lbNr(-1), _Vertices(T._Vertices), _Edges(T._Edges),
      _Faces(T._Faces), _Parent(T._Parent),
      _Children(T._Children ? new SArrayCL<TetraCL*,MaxChildrenC> (*T._Children) : 0),
      Unknowns(T.Unknowns)/*, ProcSysNum(0)*/
{
    DDD_HdrConstructorMove(&_dddH, const_cast<DDD_HEADER*>( &T._dddH));
#if DROPSDebugC&DebugSubscribeC
    T.subscribed_=false;
    subscribed_=true;
#endif
}
#endif


inline TetraCL::~TetraCL()
{
    if (_Children) delete _Children;
#ifdef _PAR
# if DROPSDebugC&DebugSubscribeC
    // This exception may be thrown, by exiting the program. If the MultiGridCL is deleted, this can happen.
    Assert(!subscribed_, DROPSErrCL("Deleteting not unsubscribed Tetra"), DebugSubscribeC);
# endif
#endif
}


inline void TetraCL::CommitRegRefMark() const
{
    for (const_EdgePIterator epiter(_Edges.begin()); epiter!=_Edges.end(); ++epiter)
        (*epiter)->IncMarkForRef();
}


inline void TetraCL::UnCommitRegRefMark() const
{
    for (const_EdgePIterator epiter(_Edges.begin()); epiter!=_Edges.end(); ++epiter)
        (*epiter)->DecMarkForRef();
}


inline void TetraCL::Close()
{
    Uint mask= 1;
    Uint newrefmark= 0;

    for (const_EdgePIterator edgep( _Edges.begin() ); edgep!=_Edges.end(); ++edgep, mask<<=1)
        if ( (*edgep)->IsMarkedForRef() ) newrefmark|= mask;
    switch (newrefmark)
    {
      case RegRefMarkC:
        _RefMark= GreenRegRefMarkC;
        return;

      case NoRefMarkC:
      // if this tetra was marked for removement, it can be removed,
      // i. e. we leave the mark for removement, so that RestrictMarks
      // will catch it in the next coarser level
        if (_RefMark!=RemoveMarkC) _RefMark= newrefmark;
        return;

      default:
        _RefMark= newrefmark;
        return;
    }
}

inline void TetraCL::LinkEdges (const ChildDataCL& childdat)
{
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
        Assert(!_ePtrs[childdat.Edges[edge]]->IsMarkedForRemovement(),"TetraCL::LinkEdge, link edge that is marked for removement", ~0);
        _Edges[edge]= _ePtrs[childdat.Edges[edge]];
    }
}

inline void TetraCL::LinkFaces (const ChildDataCL& childdat)
{
    for (Uint face=0; face<NumFacesC; ++face)
    {
        _Faces[face]= _fPtrs[childdat.Faces[face]];
        _Faces[face]->LinkTetra(this);
    }
}



inline const RefRuleCL& TetraCL::GetRefData () const
{
    return DROPS::GetRefRule(this->GetRefRule() & 63);
}

#ifdef _PAR
inline bool TetraCL::HasGhost () const
/// In the case of refinement with unknowns this procedure takes care of deleted ghost copies.
/// \return fase if there this tetra is a ghost, is local or this is master and the ghost copy on
///         another proc is marked for removement
{
    if (IsGhost())
        return false;
    if (DDD_InfoIsLocal(GetHdr()))
        return false;
    else{
        // Check if prio on other proc is ghost
        for (int *procList=GetProcList(); *procList!=-1; procList+=2)
            if (*procList!=ProcCL::MyRank() && procList[1]==PrioKilledGhost)
                return false;
    }
    return true;
}
#endif


} // end of namespace DROPS

#endif
