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
    Uint GetIdx() const { return Idx; }
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
    IdxDescCL* RowIdx;
    IdxDescCL* ColIdx;
    DataType  Data;

    void SetIdx(IdxDescCL*, IdxDescCL*);
    void Reset();
};


template<class Data= double>
class NoBndDataCL
{
  public:
    // default ctor, dtor, whatever
    typedef Data bnd_type;

    static inline bool IsOnDirBnd (const VertexCL&) { return false; }
    static inline bool IsOnNeuBnd (const VertexCL&) { return false; }
    static inline bool IsOnDirBnd (const EdgeCL&)   { return false; }
    static inline bool IsOnNeuBnd (const EdgeCL&)   { return false; }
    
    static inline bnd_type GetDirBndValue (const VertexCL&)
        { throw DROPSErrCL("NoBndDataCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
    static inline bnd_type GetDirBndValue (const EdgeCL&)
        { throw DROPSErrCL("NoBndDataCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
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

} // end of namespace DROPS


#endif
