//**************************************************************************
// File:    multigrid.h                                                    *
// Content: Classes that constitute the multigrid                          *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Eva IGPM RWTH Aachen*
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.7                                                            *
// Date:    August, 1st, 2001                                              *
// Begin:   August, 3rd, 2000                                              *
//**************************************************************************

// TODO: Use information hiding, access control and const-qualification more
//       extensively to avoid accidental changes of the multigrid structure.

#ifndef DROPS_MULTIGRID_H
#define DROPS_MULTIGRID_H

#include "geom/simplex.h"


namespace DROPS
{

//**************************************************************************
// Classes that constitute a multigrid and helpers                         *
//**************************************************************************

class MultiGridCL;
class MGBuilderCL;

template <class SimplexT>
struct TriangFillCL;

template <class SimplexT>
class TriangCL
{
  public:
    typedef std::vector<SimplexT*> LevelCont;

    typedef SimplexT**                     ptr_iterator;
    typedef const SimplexT**         const_ptr_iterator;

    typedef ptr_iter<SimplexT>             iterator;
    typedef ptr_iter<const SimplexT> const_iterator;

  private:
    mutable std::vector<LevelCont> triang_;
    MultiGridCL&                   mg_;

    inline void MaybeCreate (int lvl) const;

  public:
    TriangCL (MultiGridCL& mg);

    void clear () { triang_.clear(); }
    Uint size  (int lvl= -1) const
        { MaybeCreate( lvl); return triang_[StdIndex( lvl)].size() - 1; }

    ptr_iterator begin (int lvl= -1)
        { MaybeCreate( lvl); return &*triang_[StdIndex( lvl)].begin(); }
    ptr_iterator end   (int lvl= -1)
        {
            MaybeCreate( lvl);
            return &*(triang_[StdIndex( lvl)].end() - 1);
        }

    const_ptr_iterator begin (int lvl= -1) const
        { MaybeCreate( lvl); return const_cast<const_ptr_iterator>( &*triang_[StdIndex( lvl)].begin()); }
    const_ptr_iterator end   (int lvl= -1) const
        {
            MaybeCreate( lvl);
            return const_cast<const_ptr_iterator>( &*( triang_[StdIndex( lvl)].end() - 1));
        }

    ///@{ Cave: The returned level-container contains a zero-pointer as last element, which serves as end-iterator for the end()-functions in this class.
    LevelCont&       operator[] (int lvl)
        { MaybeCreate( lvl); return triang_[StdIndex( lvl)]; }
    const LevelCont& operator[] (int lvl) const
        { MaybeCreate( lvl); return triang_[StdIndex( lvl)]; }
    ///@}

    inline int  StdIndex    (int lvl) const;

};

typedef  TriangCL<VertexCL> TriangVertexCL;
typedef  TriangCL<EdgeCL>   TriangEdgeCL;
typedef  TriangCL<FaceCL>   TriangFaceCL;
typedef  TriangCL<TetraCL>  TriangTetraCL;

/// \brief Type of functions used to identify points on periodic boundaries,
///     that share the same dof.
typedef bool (*match_fun)(const Point3DCL&, const Point3DCL&);


class BoundaryCL
/// \brief stores boundary segments and information on periodic boundaries (if some exist)
{
  friend class MGBuilderCL;

  public:
    enum BndType {
        Per1Bnd= 1,    ///< periodic boundary 1
        Per2Bnd= 2,    ///< periodic boundary 2
        OtherBnd= 0    ///< non-periodic boundary
    };

    typedef std::vector<BndSegCL*> SegPtrCont;
    typedef std::vector<BndType>   BndTypeCont;

  private:
    SegPtrCont           Bnd_;
    mutable BndTypeCont* BndType_;
    mutable match_fun    match_;

  public:
    BoundaryCL() : BndType_(0), match_(0) {}
    /// deletes the objects pointed to in Bnd_ and BndType_
    ~BoundaryCL();

    const BndSegCL* GetBndSeg(BndIdxT idx)  const { return Bnd_[idx]; }
    BndIdxT         GetNumBndSeg()          const { return Bnd_.size(); }
    BndType         GetBndType(BndIdxT idx) const { return BndType_? (*BndType_)[idx] : OtherBnd; }

    void      SetPeriodicBnd( const BndTypeCont& type, match_fun) const;
    match_fun GetMatchFun() const { return match_; }
    bool      Matching ( const Point3DCL& p, const Point3DCL& q) const { return match_(p,q); }
    bool      HasPeriodicBnd() const { return match_; }
};


#ifdef _PAR
class LbIteratorCL;
#endif

class MultiGridCL
{

  friend class MGBuilderCL;
#ifdef _PAR
  friend class ParMultiGridCL;
  friend class LbIteratorCL;
#endif

  public:
    typedef MG_VertexContT VertexCont;
    typedef MG_EdgeContT   EdgeCont;
    typedef MG_FaceContT   FaceCont;
    typedef MG_TetraContT  TetraCont;

    typedef VertexCont::LevelCont VertexLevelCont;
    typedef EdgeCont::LevelCont   EdgeLevelCont;
    typedef FaceCont::LevelCont   FaceLevelCont;
    typedef TetraCont::LevelCont  TetraLevelCont;

    typedef VertexCont::LevelIterator             VertexIterator;
    typedef EdgeCont::LevelIterator               EdgeIterator;
    typedef FaceCont::LevelIterator               FaceIterator;
    typedef TetraCont::LevelIterator              TetraIterator;
    typedef VertexCont::const_LevelIterator const_VertexIterator;
    typedef EdgeCont::const_LevelIterator   const_EdgeIterator;
    typedef FaceCont::const_LevelIterator   const_FaceIterator;
    typedef TetraCont::const_LevelIterator  const_TetraIterator;

    typedef TriangVertexCL::iterator             TriangVertexIteratorCL;
    typedef TriangEdgeCL::iterator               TriangEdgeIteratorCL;
    typedef TriangFaceCL::iterator               TriangFaceIteratorCL;
    typedef TriangTetraCL::iterator              TriangTetraIteratorCL;
    typedef TriangVertexCL::const_iterator const_TriangVertexIteratorCL;
    typedef TriangEdgeCL::const_iterator   const_TriangEdgeIteratorCL;
    typedef TriangFaceCL::const_iterator   const_TriangFaceIteratorCL;
    typedef TriangTetraCL::const_iterator  const_TriangTetraIteratorCL;

  private:
    BoundaryCL _Bnd;
    VertexCont _Vertices;
    EdgeCont   _Edges;
    FaceCont   _Faces;
    TetraCont  _Tetras;

    TriangVertexCL _TriangVertex;
    TriangEdgeCL   _TriangEdge;
    TriangFaceCL   _TriangFace;
    TriangTetraCL  _TriangTetra;

    size_t     _version;                            // each mudification of the multigrid increments this number

#ifdef _PAR
    bool killedGhostTetra_;                         // are there ghost tetras, that are marked for removement, but has not been removed so far
    bool withUnknowns_;                             // are the unknowns on simplices
    std::list<TetraIterator> toDelGhosts_;
    bool EmptyLevel(Uint lvl)
        { return _Vertices[lvl].empty()&&_Edges[lvl].empty()&&_Faces[lvl].empty()&&_Tetras[lvl].empty(); }
#endif

    void PrepareModify   () { _Vertices.PrepareModify(); _Edges.PrepareModify(); _Faces.PrepareModify(); _Tetras.PrepareModify(); }
    void FinalizeModify  () { _Vertices.FinalizeModify(); _Edges.FinalizeModify(); _Faces.FinalizeModify(); _Tetras.FinalizeModify(); }
    void AppendLevel     () { _Vertices.AppendLevel(); _Edges.AppendLevel(); _Faces.AppendLevel(); _Tetras.AppendLevel(); }
    void RemoveLastLevel () { _Vertices.RemoveLastLevel(); _Edges.RemoveLastLevel(); _Faces.RemoveLastLevel(); _Tetras.RemoveLastLevel(); }

    void ClearTriangCache () { _TriangVertex.clear(); _TriangEdge.clear(); _TriangFace.clear(); _TriangTetra.clear(); }

    void RestrictMarks (Uint Level) { std::for_each( _Tetras[Level].begin(), _Tetras[Level].end(), std::mem_fun_ref(&TetraCL::RestrictMark)); }
    void CloseGrid     (Uint);
    void UnrefineGrid  (Uint);
    void RefineGrid    (Uint);

  public:
    MultiGridCL (const MGBuilderCL& Builder);
    MultiGridCL (const MultiGridCL&); // Dummy
    // default ctor

#ifdef _PAR
    bool KilledGhosts()      const              /// Check if there are ghost tetras, that are marked for removement, but has not been removed so far
        { return killedGhostTetra_; }
    bool UnknownsForRefine() const              /// Check if there are unknowns on simplices
        { return withUnknowns_; }
#endif

    const BoundaryCL& GetBnd     () const { return _Bnd; }
    const VertexCont& GetVertices() const { return _Vertices; }
    const EdgeCont&   GetEdges   () const { return _Edges; }
    const FaceCont&   GetFaces   () const { return _Faces; }
    const TetraCont&  GetTetras  () const { return _Tetras; }

    const TriangVertexCL& GetTriangVertex () const { return _TriangVertex; }
    const TriangEdgeCL&   GetTriangEdge   () const { return _TriangEdge; }
    const TriangFaceCL&   GetTriangFace   () const { return _TriangFace; }
    const TriangTetraCL&  GetTriangTetra  () const { return _TriangTetra; }

    VertexIterator GetVerticesBegin (int Level=-1) { return _Vertices.level_begin( Level); }
    VertexIterator GetVerticesEnd   (int Level=-1) { return _Vertices.level_end( Level); }
    EdgeIterator   GetEdgesBegin    (int Level=-1)  { return _Edges.level_begin( Level); }
    EdgeIterator   GetEdgesEnd      (int Level=-1)  { return _Edges.level_end( Level); }
    FaceIterator   GetFacesBegin    (int Level=-1) { return _Faces.level_begin( Level); }
    FaceIterator   GetFacesEnd      (int Level=-1) { return _Faces.level_end( Level); }
    TetraIterator  GetTetrasBegin   (int Level=-1) { return _Tetras.level_begin( Level); }
    TetraIterator  GetTetrasEnd     (int Level=-1) { return _Tetras.level_end( Level); }
    const_VertexIterator GetVerticesBegin (int Level=-1) const { return _Vertices.level_begin( Level); }
    const_VertexIterator GetVerticesEnd   (int Level=-1) const { return _Vertices.level_end( Level); }
    const_EdgeIterator   GetEdgesBegin    (int Level=-1) const { return _Edges.level_begin( Level); }
    const_EdgeIterator   GetEdgesEnd      (int Level=-1) const { return _Edges.level_end( Level); }
    const_FaceIterator   GetFacesBegin    (int Level=-1) const { return _Faces.level_begin( Level); }
    const_FaceIterator   GetFacesEnd      (int Level=-1) const { return _Faces.level_end( Level); }
    const_TetraIterator  GetTetrasBegin   (int Level=-1) const { return _Tetras.level_begin( Level); }
    const_TetraIterator  GetTetrasEnd     (int Level=-1) const { return _Tetras.level_end( Level); }

    VertexIterator GetAllVertexBegin (int= -1     ) { return _Vertices.begin(); }
    VertexIterator GetAllVertexEnd   (int Level=-1) { return _Vertices.level_end( Level); }
    EdgeIterator   GetAllEdgeBegin   (int= -1     ) { return _Edges.begin(); }
    EdgeIterator   GetAllEdgeEnd     (int Level=-1) { return _Edges.level_end( Level); }
    FaceIterator   GetAllFaceBegin   (int= -1     ) { return _Faces.begin(); }
    FaceIterator   GetAllFaceEnd     (int Level=-1) { return _Faces.level_end( Level); }
    TetraIterator  GetAllTetraBegin  (int= -1     ) { return _Tetras.begin(); }
    TetraIterator  GetAllTetraEnd    (int Level=-1) { return _Tetras.level_end( Level); }
    const_VertexIterator GetAllVertexBegin (int= -1     ) const { return _Vertices.begin(); }
    const_VertexIterator GetAllVertexEnd   (int Level=-1) const { return _Vertices.level_end( Level); }
    const_EdgeIterator   GetAllEdgeBegin   (int= -1     ) const  { return _Edges.begin(); }
    const_EdgeIterator   GetAllEdgeEnd     (int Level=-1) const  { return _Edges.level_end( Level); }
    const_FaceIterator   GetAllFaceBegin   (int= -1     ) const { return _Faces.begin(); }
    const_FaceIterator   GetAllFaceEnd     (int Level=-1) const { return _Faces.level_end( Level); }
    const_TetraIterator  GetAllTetraBegin  (int= -1     ) const { return _Tetras.begin(); }
    const_TetraIterator  GetAllTetraEnd    (int Level=-1) const { return _Tetras.level_end( Level); }

    TriangVertexIteratorCL GetTriangVertexBegin (int Level=-1) { return _TriangVertex.begin( Level); }
    TriangVertexIteratorCL GetTriangVertexEnd   (int Level=-1) { return _TriangVertex.end( Level); }
    TriangEdgeIteratorCL   GetTriangEdgeBegin   (int Level=-1) { return _TriangEdge.begin( Level); }
    TriangEdgeIteratorCL   GetTriangEdgeEnd     (int Level=-1) { return _TriangEdge.end( Level); }
    TriangFaceIteratorCL   GetTriangFaceBegin   (int Level=-1) { return _TriangFace.begin( Level); }
    TriangFaceIteratorCL   GetTriangFaceEnd     (int Level=-1) { return _TriangFace.end( Level); }
    TriangTetraIteratorCL  GetTriangTetraBegin  (int Level=-1) { return _TriangTetra.begin( Level); }
    TriangTetraIteratorCL  GetTriangTetraEnd    (int Level=-1) { return _TriangTetra.end( Level); }
    const_TriangVertexIteratorCL GetTriangVertexBegin (int Level=-1) const { return _TriangVertex.begin( Level); }
    const_TriangVertexIteratorCL GetTriangVertexEnd   (int Level=-1) const { return _TriangVertex.end( Level); }
    const_TriangEdgeIteratorCL   GetTriangEdgeBegin   (int Level=-1) const { return _TriangEdge.begin( Level); }
    const_TriangEdgeIteratorCL   GetTriangEdgeEnd     (int Level=-1) const { return _TriangEdge.end( Level); }
    const_TriangFaceIteratorCL   GetTriangFaceBegin   (int Level=-1) const { return _TriangFace.begin( Level); }
    const_TriangFaceIteratorCL   GetTriangFaceEnd     (int Level=-1) const { return _TriangFace.end( Level); }
    const_TriangTetraIteratorCL  GetTriangTetraBegin  (int Level=-1) const { return _TriangTetra.begin( Level); }
    const_TriangTetraIteratorCL  GetTriangTetraEnd    (int Level=-1) const { return _TriangTetra.end( Level); }

    Uint GetLastLevel() const { return _Tetras.GetNumLevel()-1; }
    Uint GetNumLevel () const { return _Tetras.GetNumLevel(); }

    void   IncrementVersion() {++_version; }                    ///< Increment version of the multigrid
    size_t GetVersion() const { return _version; }              ///< Get version of the multigrid

    void Refine();                                              // in parallel mode, this function uses a parallel version for refinement!

    void Scale( double);
    void Transform( Point3DCL (*mapping)(const Point3DCL&));
    void MakeConsistentNumbering();
    void SizeInfo(std::ostream&);                               // all procs have to call this function in parallel mode!
    void ElemInfo(std::ostream&, int Level= -1);                // all procs have to call this function in parallel mode
#ifdef _PAR
    Uint GetNumDistributedObjects() const;                      // get number of distributed objects
    Uint GetNumTriangTetra(int Level=-1);                       // get number of tetraeder of a given level
    Uint GetNumTriangFace(int Level=-1);                        // get number of faces of a given level
    Uint GetNumDistributedFaces(int Level=-1);                  // get number of faces on processor boundary
#endif

    bool IsSane (std::ostream&, int Level=-1) const;
};


class PeriodicEdgesCL
/// \brief handles edges on periodic boundaries.
///
/// This class is used by the refinement algorithm in MultiGridCL to accumulate the MFR counters on linked periodic edges.
/// This assures that periodic boundaries are matching after refinement.
{
  public:
    typedef std::pair<EdgeCL*,EdgeCL*>  IdentifiedEdgesT;
    typedef std::list<IdentifiedEdgesT> PerEdgeContT;
    typedef PerEdgeContT::iterator      iterator;
    typedef MultiGridCL::EdgeIterator   EdgeIterator;

  private:
    PerEdgeContT      list_;
    MultiGridCL&      mg_;

    /// recompute data structure
    void Recompute( EdgeIterator begin, EdgeIterator end);
    /// accumulate local MFR counters of periodic edges and store the sum in the MFR counter
    void Accumulate();
    /// delete all data
    void Shrink();

  public:
    PeriodicEdgesCL( MultiGridCL& mg) : mg_(mg) {}
    // standard dtor

    BoundaryCL::BndType GetBndType( const EdgeCL& e) const;
    void AccumulateMFR( int lvl);
    /// print out list of identified edges for debugging
    void DebugInfo(std::ostream&);
};


template <class SimplexT>
struct TriangFillCL
{
  static void // not defined
  fill (MultiGridCL& mg, typename TriangCL<SimplexT>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<VertexCL>
{
    static void fill (MultiGridCL& mg, TriangCL<VertexCL>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<EdgeCL>
{
    static void fill (MultiGridCL& mg, TriangCL<EdgeCL>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<FaceCL>
{
    static void fill (MultiGridCL& mg, TriangCL<FaceCL>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<TetraCL>
{
    static void fill (MultiGridCL& mg, TriangCL<TetraCL>::LevelCont& c, int lvl);
};


#define DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) \
for (DROPS::TriangVertexCL::iterator it( mg.GetTriangVertexBegin( lvl)), end__( mg.GetTriangVertexEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) \
for (DROPS::TriangVertexCL::const_iterator it( mg.GetTriangVertexBegin( lvl)), end__( mg.GetTriangVertexEnd( lvl)); it != end__; ++it)


#define DROPS_FOR_TRIANG_EDGE( mg, lvl, it) \
for (DROPS::TriangEdgeCL::iterator it( mg.GetTriangEdgeBegin( lvl)), end__( mg.GetTriangEdgeEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) \
for (DROPS::TriangEdgeCL::const_iterator it( mg.GetTriangEdgeBegin( lvl)), end__( mg.GetTriangEdgeEnd( lvl)); it != end__; ++it)


#define DROPS_FOR_TRIANG_FACE( mg, lvl, it) \
for (DROPS::TriangFaceCL::iterator it( mg.GetTriangFaceBegin( lvl)), end__( mg.GetTriangFaceEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_ALL_TRIANG_CONST_FACE( mg, lvl, it) \
for (DROPS::TriangFaceCL::const_iterator FaceCL* it( mg.GetTriangFaceBegin( lvl)), end__( mg.GetTriangFaceEnd( lvl)); it != end__; ++it)


#define DROPS_FOR_TRIANG_TETRA( mg, lvl, it) \
for (DROPS::TriangTetraCL::iterator it( mg.GetTriangTetraBegin( lvl)), end__( mg.GetTriangTetraEnd( lvl)); it != end__; ++it)

#define DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) \
for (DROPS::TriangTetraCL::const_iterator it( mg.GetTriangTetraBegin( lvl)), end__( mg.GetTriangTetraEnd( lvl)); it != end__; ++it)


class MGBuilderCL
{
  protected:
    MultiGridCL::VertexCont& GetVertices(MultiGridCL* _MG) const { return _MG->_Vertices; }
    MultiGridCL::EdgeCont&   GetEdges   (MultiGridCL* _MG) const { return _MG->_Edges; }
    MultiGridCL::FaceCont&   GetFaces   (MultiGridCL* _MG) const { return _MG->_Faces; }
    MultiGridCL::TetraCont&  GetTetras  (MultiGridCL* _MG) const { return _MG->_Tetras; }
    BoundaryCL::SegPtrCont&  GetBnd     (MultiGridCL* _MG) const { return _MG->_Bnd.Bnd_; }
    void PrepareModify  (MultiGridCL* _MG) const { _MG->PrepareModify(); }
    void FinalizeModify (MultiGridCL* _MG) const { _MG->FinalizeModify(); }
    void AppendLevel    (MultiGridCL* _MG) const { _MG->AppendLevel(); }
    void RemoveLastLevel(MultiGridCL* _MG) const { _MG->RemoveLastLevel(); }

  public:
    // default ctor
    virtual ~MGBuilderCL() {}
    virtual void buildBoundary(MultiGridCL* _MG) const = 0;
    virtual void build(MultiGridCL*) const = 0;
};

class LocatorCL;

/// \brief Class for combining a tetrahedra wtih its bary center
class LocationCL
{
  private:
    const TetraCL* _Tetra;
    SVectorCL<4>   _Coord;

  public:
    LocationCL()
        : _Tetra(0), _Coord() {}
    LocationCL(const TetraCL* t, const SVectorCL<4>& p)
        : _Tetra(t), _Coord(p) {}
    LocationCL(const LocationCL& loc)
        : _Tetra(loc._Tetra), _Coord(loc._Coord) {}

#ifndef _PAR
    bool IsValid() const                        ///< Check if the tetrahedra is set
        { return _Tetra; }
#else
    bool IsValid(Uint lvl) const                ///< Check if the tetrahedra is set
        { return _Tetra && _Tetra->IsInTriang(lvl)/* && _Tetra->MayStoreUnk()*/; }
#endif
    const TetraCL& GetTetra() const             ///< Get a reference to the tetrahedra
        { return *_Tetra; }
    const SVectorCL<4>& GetBaryCoord() const    ///< Get bary center coordinates of the tetrahedra
        { return _Coord; }

    friend class LocatorCL;
};

/// \brief Find a tetrahedra that surrounds a given Point
class LocatorCL
{
  private:
    static bool InTetra(const SVectorCL<4>& b, double tol= 0)
        { return b[0] >= -tol && b[0] <= 1.+tol
              && b[1] >= -tol && b[1] <= 1.+tol
              && b[2] >= -tol && b[2] <= 1.+tol
              && b[3] >= -tol && b[3] <= 1.+tol; }

    static void MakeMatrix(const TetraCL& t, SMatrixCL<4,4>& M)
    {
        for (Uint j=0; j<4; ++j)
        {
            for (Uint i=0; i<3; ++i)
                M(i, j)= t.GetVertex(j)->GetCoord()[i];
            M(3, j)= 1.;
        }
    }

  public:
    // default ctor, copy-ctor, dtor, assignment-op

    /// \brief Locate a point with a given tetrahedra
    static void
    LocateInTetra(LocationCL&, Uint, const Point3DCL&, double tol= 0);
    /// \brief Find the tetrahedra that surounds a point
    static void
    Locate(LocationCL&, const MultiGridCL& MG, int, const Point3DCL&, double tol= 1e-14);
};
// inline functions

template <class SimplexT>
  TriangCL<SimplexT>::TriangCL (MultiGridCL& mg)
      : triang_( mg.GetNumLevel()), mg_( mg)
{}

template <class SimplexT>
  inline int
  TriangCL<SimplexT>::StdIndex(int lvl) const
{
    return lvl >= 0 ? lvl : lvl + mg_.GetNumLevel();
}

template <class SimplexT>
  inline void
  TriangCL<SimplexT>::MaybeCreate(int lvl) const
{
    const int level= StdIndex( lvl);
    Assert ( level >= 0 && level < static_cast<int>( mg_.GetNumLevel()),
        DROPSErrCL( "TriangCL::MaybeCreate: Wrong level."), DebugContainerC);
    if (triang_.size() != mg_.GetNumLevel()) {
        triang_.clear();
        triang_.resize( mg_.GetNumLevel());
    }
    if (triang_[level].empty()) {
        TriangFillCL<SimplexT>::fill( mg_, triang_[level], level);
        // Append a zero-pointer as explicit end-iterator of the sequence.
        triang_[level].push_back( 0);
    }
}

#ifdef _PAR
inline bool TetraCL::IsExclusive( Priority prio) const
/** A tetra is exclusive, if the tetra is local stored or the tetra is master and owned by
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
    Assert(minproc != ProcCL::Size(), DROPSErrCL("TetraCL::IsExclusive: no exclusive proc for vertex!"), DebugParallelC);
    if (me==minproc)
        return true;
    else
        return false;
}
#endif


void circumcircle(const TetraCL& t, Point3DCL& c, double& r);
void circumcircle(const TetraCL& t, Uint face, Point3DCL& c, double& r);

/// calculates the transpose of the transformation  Tetra -> RefTetra
inline void GetTrafoTr( SMatrixCL<3,3>& T, double& det, const Point3DCL pt[4])
{
    double M[3][3];
    const Point3DCL& pt0= pt[0];
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[j][i]= pt[i+1][j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    T(0,0)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    T(0,1)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    T(0,2)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    T(1,0)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    T(1,1)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    T(1,2)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    T(2,0)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    T(2,1)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    T(2,2)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
}

/// calculates the transpose of the transformation  Tetra -> RefTetra
inline void GetTrafoTr( SMatrixCL<3,3>& T, double& det, const TetraCL& t)
{
    double M[3][3];
    const Point3DCL& pt0= t.GetVertex(0)->GetCoord();
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[j][i]= t.GetVertex(i+1)->GetCoord()[j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    T(0,0)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    T(0,1)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    T(0,2)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    T(1,0)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    T(1,1)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    T(1,2)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    T(2,0)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    T(2,1)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    T(2,2)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
}


void MarkAll (MultiGridCL&);
void UnMarkAll (MultiGridCL&);

#ifdef _PAR
std::string PrioToString(Uint prio);
Ulint GetExclusiveVerts (const MultiGridCL&, Priority prio=PrioMaster, int lvl=-1);
Ulint GetExclusiveEdges (const MultiGridCL&, Priority prio=PrioMaster, int lvl=-1);
#endif

} // end of namespace DROPS

#endif
