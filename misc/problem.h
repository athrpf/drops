//**************************************************************************
// File:    problem.h                                                      *
// Content: classes that constitute a problem                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

#ifndef DROPS_PROBLEM_H
#define DROPS_PROBLEM_H


#include "geom/multigrid.h"
#include "geom/builder.h"
#include "num/spmat.h"


namespace DROPS
{

void BndCondInfo( BndCondT, std::ostream&);

class IdxDescCL
{
  private:
    static std::vector<bool> IdxFree;
    Uint                     Idx;
    
    Uint GetFreeIdx();
    
  public:
    Uint TriangLevel;
    Uint NumUnknownsVertex;
    Uint NumUnknownsEdge;
    Uint NumUnknownsFace;
    Uint NumUnknownsTetra;
    IdxT NumUnknowns;

    IdxDescCL( Uint unkVertex= 0, Uint unkEdge= 0, Uint unkFace= 0, Uint unkTetra= 0) 
      : Idx( GetFreeIdx()), NumUnknownsVertex( unkVertex), NumUnknownsEdge( unkEdge),
        NumUnknownsFace( unkFace), NumUnknownsTetra( unkTetra), NumUnknowns( 0) {}
    IdxDescCL( const IdxDescCL& orig);              
        // WARNING:  "orig" will be invalidated as the private "Idx" should not be the same for two different objects...
    ~IdxDescCL() { if (Idx!=NoIdx) IdxFree[Idx]= true; }
    
    IdxDescCL& operator= ( const IdxDescCL&);  // not implemented as the private "Idx" should not be the same for two different objects...
    void swap( IdxDescCL&);
    
    void Set( Uint unkVertex, Uint unkEdge= 0, Uint unkFace= 0, Uint unkTetra= 0);
    Uint GetIdx() const { if (Idx==NoIdx) throw DROPSErrCL("IdxDescCL::GetIdx: invalid index. Probably using copy instead of original IdxDescCL-object."); return Idx; }

    // Compare two IdxDescCL-objects. If a multigrid is given via mg, the
    // unknown-numbers on it are compared, too.
    static bool
    Equal(IdxDescCL& i, IdxDescCL& j, const MultiGridCL* mg= 0);
};


inline void
GetLocalNumbP1NoBnd(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx)
// Copies P1-unknown-indices from idx on s into Numb; assumes that all
// vertices have unknowns (NoBndDataCL and the like).
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i) {
        Numb[i]= s.GetVertex( i)->Unknowns( sys);
    }
}


//**************************************************************************
// Class:   LocalNumbP2CL                                                  *
// Purpose: Use this class to collect the indices of unknowns, boundary-   *
//          conditions and boundary-segment-numbers.                       *
// Description:                                                            *
//    num: field of unknown-indices; NoIdx, iff the degree of freedom lies *
//         on a boundary without unknowns.                                 *
//    bndnum: On boundaries, the number of the relevant BndSegDataCL-object*
//            in the corresponding BndDataCL-object, else NoBndC.          *
//    bc: The relevant BndCondT, NoBC in the interior dofs.                *
//**************************************************************************
class LocalNumbP2CL
{
  public:
    IdxT     num   [10]; // formerly Numb
    BndIdxT  bndnum[10];
    BndCondT bc    [10];

    LocalNumbP2CL() {}
    template<class BndDataT>
      LocalNumbP2CL(const TetraCL&, const IdxDescCL&, const BndDataT&);
    
    template<class BndDataT>
      void
      assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd);

    bool WithUnknowns(IdxT i) const { return num[i] != NoIdx; }
};


template<class T>
class VecDescBaseCL
{
  public:
    typedef T DataType;
    
    VecDescBaseCL()
        :RowIdx(0) {}
    VecDescBaseCL( IdxDescCL* idx) { SetIdx( idx); }
    IdxDescCL* RowIdx;
    DataType  Data;

    Uint GetLevel() const { return RowIdx->TriangLevel; }
    void SetIdx(IdxDescCL*);
    void Clear();
    void Reset();
};

typedef VecDescBaseCL<VectorCL> VecDescCL;


class MatDescCL
{
  public:
    typedef MatrixCL DataType;
    
    MatDescCL()
        :RowIdx(0), ColIdx(0) {}
    MatDescCL(IdxDescCL* r, IdxDescCL* c)
        { SetIdx( r, c); }
    IdxDescCL* RowIdx;
    IdxDescCL* ColIdx;
    DataType  Data;

    Uint GetRowLevel() const { return RowIdx->TriangLevel; }
    Uint GetColLevel() const { return ColIdx->TriangLevel; }
    void SetIdx(IdxDescCL*, IdxDescCL*);
    void Reset();
};


template <class Coeff, class BndData>
class ProblemCL
{
  public:
    typedef Coeff    CoeffCL;
    typedef BndData  BndDataCL;

  protected:
    bool         _myMG;
    MultiGridCL& _MG;
    CoeffCL      _Coeff;      // rechte Seite, Koeffizienten der PDE
    BndDataCL    _BndData;    // Randwerte

  public:
    ProblemCL(const MGBuilderCL& mgbuilder, const CoeffCL& coeff, const BndDataCL& bnddata)
    : _myMG( true), _MG( *new MultiGridCL( mgbuilder)), _Coeff( coeff), _BndData( bnddata) {}
    ProblemCL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bnddata)
    : _myMG( false), _MG( mg), _Coeff( coeff), _BndData( bnddata) {}
    ~ProblemCL() { if (_myMG) delete &_MG; }

    MultiGridCL&       GetMG()            { return _MG; }
    const MultiGridCL& GetMG()      const { return _MG; }
    const CoeffCL&     GetCoeff()   const { return _Coeff; }
    const BndDataCL&   GetBndData() const { return _BndData; }
};


template<class BndDataT>
  LocalNumbP2CL::LocalNumbP2CL(const TetraCL& s, const IdxDescCL& idx,
      const BndDataT& bnd)
{
    this->assign( s, idx, bnd);
}

template<class BndDataT>
  void
  LocalNumbP2CL::assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd)
{
    BndIdxT bidx= 0;
    const Uint sys= idx.GetIdx();
    
    for (Uint i= 0; i < NumVertsC; ++i)
        if (NoBC == (bc[i]= bnd.GetBC( *s.GetVertex( i), bidx))) {
            bndnum[i]= NoBndC;
            num[i]= s.GetVertex( i)->Unknowns( sys);
        }
        else {
            bndnum[i]= bidx;
            num[i]= (bnd.GetBndSeg( bidx).WithUnknowns())
                ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
        }
    for (Uint i= 0; i< NumEdgesC; ++i)
        if (NoBC == (bc[i+NumVertsC]= bnd.GetBC( *s.GetEdge( i), bidx))) {
            bndnum[i+NumVertsC]= NoBndC;
            num[i+NumVertsC]= s.GetEdge( i)->Unknowns( sys);
        }
        else {
            bndnum[i+NumVertsC]= bidx;
            num[i+NumVertsC]= (bnd.GetBndSeg( bidx).WithUnknowns())
                ? s.GetEdge( i)->Unknowns( sys) : NoIdx;
        }
}

template<class T>
void VecDescBaseCL<T>::SetIdx(IdxDescCL* idx)
{
    // create new vector of unknowns and set up the vector description
    RowIdx = idx;
    Data.resize(0);
    Data.resize(idx->NumUnknowns);
}

template<class T>
void VecDescBaseCL<T>::Clear()
{
    Data.resize(0);
    Data.resize(RowIdx->NumUnknowns);
}

template<class T>
void VecDescBaseCL<T>::Reset()
{
    RowIdx = 0;
    Data.resize(0);
}


// ----------------------------------------------------------------------------
//                      Routines for numbering of unknowns
// ----------------------------------------------------------------------------

template<class BndDataT>
void CreateNumbOnVertex( const Uint idx, IdxT& counter, Uint stride,
                         const MultiGridCL::TriangVertexIteratorCL& begin,
                         const MultiGridCL::TriangVertexIteratorCL& end,
                         const BndDataT& Bnd)
// Allocates memory for the Unknown-indices in system idx on all vertices
// between begin and end.
// The first number used is the initial value of counter, the next numbers are counter+stride,
// counter+2*stride, and so on.
// Upon return, counter contains the first number, that was not used,
// that is #Unknowns+stride. 
// Vertices on Dirichlet boundaries are skipped.
{
    if (stride == 0) return;
    for (MultiGridCL::TriangVertexIteratorCL it= begin; it != end; ++it)
    {
        if ( !Bnd.IsOnDirBnd( *it) )
        {        
            it->Unknowns.Prepare( idx);
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
}

template<class BndDataT>
void CreateNumbOnEdge( const Uint idx, IdxT& counter, Uint stride,
                       const MultiGridCL::TriangEdgeIteratorCL& begin,
                       const MultiGridCL::TriangEdgeIteratorCL& end,
                       const BndDataT& Bnd)
// Allocates memory for the Unknown-indices in system idx on all edges
// between begin and end.
// The first number used is the initial value of counter, the next numbers are counter+stride,
// counter+2*stride, and so on.
// Upon return, counter contains the first number, that was not used,
// that is #Unknowns+stride. 
// Edges on Dirichlet boundaries are skipped.
{
    if (stride == 0) return;
    for (MultiGridCL::TriangEdgeIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {        
            it->Unknowns.Prepare( idx);
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
}

template<class BndDataT>
void CreateNumbOnFace( const Uint idx, IdxT& counter, Uint stride,
                       const MultiGridCL::TriangFaceIteratorCL& begin,
                       const MultiGridCL::TriangFaceIteratorCL& end,
                       const BndDataT& Bnd)
// Allocates memory for the Unknown-indices in system idx on all faces
// between begin and end.
// The first number used is the initial value of counter, the next numbers are counter+stride,
// counter+2*stride, and so on.
// Upon return, counter contains the first number, that was not used,
// that is #Unknowns+stride. 
// Faces on Dirichlet boundaries are skipped.
{
    if (stride == 0) return;
    for (MultiGridCL::TriangFaceIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {        
            it->Unknowns.Prepare( idx);
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
}


// Allocates memory for the Unknown-indices in system idx on all tetras
// between begin and end.
// The first number used is the initial value of counter, the next numbers are counter+stride,
// counter+2*stride, and so on.
// Upon return, counter contains the first number, that was not used,
// that is #Unknowns+stride. 
void CreateNumbOnTetra( const Uint idx, IdxT& counter, Uint stride,
                        const MultiGridCL::TriangTetraIteratorCL& begin,
                        const MultiGridCL::TriangTetraIteratorCL& end);

template <class Iter>
inline void
DeleteNumbOnSimplex( Uint idx, const Iter& begin, const Iter& end)
// TODO: This does nothing anymore. It might be useful to have a marker for an invalid
// index, e. g. to be able to tell, which indices are new after a refinement step.
// deletes the memory for the Unknown-indices allocated by CreateNumbOn<Simplex>
{

    for (Iter it=begin; it!=end; ++it) 
        if (it->Unknowns.Exist() && it->Unknowns.Exist( idx) ) 
            it->Unknowns.Invalidate( idx);
}

typedef bool (*match_fun)(const Point3DCL&, const Point3DCL&);

template<class BndDataT>
void CreatePeriodicNumbOnVertex( Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const MultiGridCL::TriangVertexIteratorCL& begin,
                        const MultiGridCL::TriangVertexIteratorCL& end,
                        const BndDataT& Bnd)
{
    if (stride == 0) return;
    
    typedef std::list<VertexCL*> psetT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (MultiGridCL::TriangVertexIteratorCL it= begin; it!=end; ++it)
    {
        if ( Bnd.IsOnDirBnd( *it) ) continue;
        it->Unknowns.Prepare( idx);
        if (Bnd.IsOnPerBnd( *it))
        {
            if (Bnd.GetBC( *it)==Per1BC)
            {
                s1.push_back( &*it);
                it->Unknowns( idx)= counter;
                counter+= stride;
            }
            else
                s2.push_back( &*it);
        }
        else
        {
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
    if (s1.size() != s2.size())
        throw DROPSErrCL( "CreatePeriodicNumbOnVertex: Periodic boundaries do not match!");
    // match objects in s1 and s2
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; ++it2)
            if (match( (*it1)->GetCoord(), (*it2)->GetCoord()))
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2 and stop search
                s2.erase( it2);
                break;
            }
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnVertex: Periodic boundaries do not match!");
}

template<class BndDataT>
void CreatePeriodicNumbOnEdge( Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const MultiGridCL::TriangEdgeIteratorCL& begin,
                        const MultiGridCL::TriangEdgeIteratorCL& end,
                        const BndDataT& Bnd)
{
    if (stride == 0) return;
    
    typedef std::list<EdgeCL*> psetT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (MultiGridCL::TriangEdgeIteratorCL it= begin; it!=end; ++it)
    {
        if ( Bnd.IsOnDirBnd( *it) ) continue;
        it->Unknowns.Prepare( idx);
        if (Bnd.IsOnPerBnd( *it))
        {
            if (Bnd.GetBC( *it)==Per1BC)
            {
                s1.push_back( &*it);
                it->Unknowns( idx)= counter;
                counter+= stride;
            }
            else
                s2.push_back( &*it);
        }
        else
        {
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
    if (s1.size() != s2.size())
        throw DROPSErrCL( "CreatePeriodicNumbOnEdge: Periodic boundaries do not match!");
    // match objects in s1 and s2
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; ++it2)
            if (match( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2 and stop search
                s2.erase( it2);
                break;
            }
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnEdge: Periodic boundaries do not match!");
}

template<class BndDataT>
void CreatePeriodicNumbOnFace( Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const MultiGridCL::TriangFaceIteratorCL& begin,
                        const MultiGridCL::TriangFaceIteratorCL& end,
                        const BndDataT& Bnd)
{
    if (stride == 0) return;
    
    typedef std::list<FaceCL*> psetT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (MultiGridCL::TriangFaceIteratorCL it= begin; it!=end; ++it)
    {
        if ( Bnd.IsOnDirBnd( *it) ) continue;
        it->Unknowns.Prepare( idx);
        if (Bnd.IsOnPerBnd( *it))
        {
            if (Bnd.GetBC( *it)==Per1BC)
            {
                s1.push_back( &*it);
                it->Unknowns( idx)= counter;
                counter+= stride;
            }
            else
                s2.push_back( &*it);
        }
        else
        {
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
    if (s1.size() != s2.size())
        throw DROPSErrCL( "CreatePeriodicNumbOnFace: Periodic boundaries do not match!");
    // match objects in s1 and s2
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; ++it2)
            if (match( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2 and stop search
                s2.erase( it2);
                break;
            }
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnFace: Periodic boundaries do not match!");
}

template<class BndDataT>
void CreateNumb( Uint level, IdxDescCL& idx, MultiGridCL& mg, const BndDataT& Bnd, match_fun match= 0)
// used for numbering of the Unknowns depending on the index idx.
// if a matching function is given, numbering for periodic boundaries is considered.
// sets up the description of the index idx,
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, that IdxDesc[idxnr].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx.TriangLevel = level;
    idx.NumUnknowns = 0;

    const Uint idxnum= idx.GetIdx();

    // allocate space for indices; number unknowns in TriangLevel level
    if (match)
    {
        if (idx.NumUnknownsVertex)
            CreatePeriodicNumbOnVertex( idxnum, idx.NumUnknowns, idx.NumUnknownsVertex, match,
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd);
        if (idx.NumUnknownsEdge)
            CreatePeriodicNumbOnEdge( idxnum, idx.NumUnknowns, idx.NumUnknownsEdge, match,
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd);
        if (idx.NumUnknownsFace)
            CreatePeriodicNumbOnFace( idxnum, idx.NumUnknowns, idx.NumUnknownsFace, match,
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd);
        if (idx.NumUnknownsTetra)
            CreateNumbOnTetra( idxnum, idx.NumUnknowns, idx.NumUnknownsTetra,
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level));
    }
    else
    {
        if (idx.NumUnknownsVertex)
            CreateNumbOnVertex( idxnum, idx.NumUnknowns, idx.NumUnknownsVertex,
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd);
        if (idx.NumUnknownsEdge)
            CreateNumbOnEdge( idxnum, idx.NumUnknowns, idx.NumUnknownsEdge,
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd);
        if (idx.NumUnknownsFace)
            CreateNumbOnFace( idxnum, idx.NumUnknowns, idx.NumUnknownsFace,
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd);
        if (idx.NumUnknownsTetra)
            CreateNumbOnTetra( idxnum, idx.NumUnknowns, idx.NumUnknownsTetra,
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level));
    }
}

} // end of namespace DROPS


#endif
