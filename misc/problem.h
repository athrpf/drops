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


typedef SparseMatBaseCL<double>    MatrixCL;
typedef SparseMatBuilderCL<double> MatrixBuilderCL;
typedef VectorBaseCL<double>       VectorCL;


class IdxDescCL
{
  public:
    Uint TriangLevel;
    Uint Idx;
    Uint NumUnknownsVertex;
    Uint NumUnknownsEdge;
    Uint NumUnknownsFace;
    Uint NumUnknownsTetra;
    IdxT NumUnknowns;
    // Uit _NumElimUnknowns;
    void Set(Uint idxnum, Uint unkVertex, Uint unkEdge, Uint unkFace, Uint unkTetra);
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


template <class MGB, class Coeff, class BndData>
class ProblemCL
{
  public:
    typedef MGB      MultiGridBuilderCL;
    typedef Coeff    CoeffCL;
    typedef BndData  BndDataCL;

  protected:
    bool         _myMG;
    MultiGridCL& _MG;
    CoeffCL      _Coeff;      // rechte Seite, Koeffizienten der PDE
    BndDataCL    _BndData;    // Randwerte

  public:
    ProblemCL(const MultiGridBuilderCL& mgbuilder, const CoeffCL& coeff, const BndDataCL& bnddata)
    : _myMG(true), _MG(*new MultiGridCL(mgbuilder)), _Coeff(coeff), _BndData(bnddata) {}
    ProblemCL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bnddata)
    : _myMG(false), _MG(mg), _Coeff(coeff), _BndData(bnddata) {}
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
void CreateNumbOnVertex( const Uint idx, IdxT& counter, Uint NumUnknown,
                         const MultiGridCL::TriangVertexIteratorCL& begin,
                         const MultiGridCL::TriangVertexIteratorCL& end,
                         const BndDataT& Bnd)
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence. Simplices on Dirichlet boundaries are skipped.
{
    if (NumUnknown == 0) return;
    for ( MultiGridCL::TriangVertexIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {        
            it->Unknowns.Prepare( idx, NumUnknown);
            for (Uint i=0; i<NumUnknown; ++i) it->Unknowns(idx)[i]= counter++;
        }
    }
}

template<class BndDataT>
void CreateNumbOnEdge( const Uint idx, IdxT& counter, Uint NumUnknown,
                       const MultiGridCL::TriangEdgeIteratorCL& begin,
                       const MultiGridCL::TriangEdgeIteratorCL& end,
                       const BndDataT& Bnd)
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence. Simplices on Dirichlet boundaries are skipped.
{
    if (NumUnknown == 0) return;
    for ( MultiGridCL::TriangEdgeIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {        
            it->Unknowns.Prepare( idx, NumUnknown);
            for (Uint i=0; i<NumUnknown; ++i) it->Unknowns(idx)[i]= counter++;
        }
    }
}

template<class BndDataT>
void CreateNumbOnFace( const Uint idx, IdxT& counter, Uint NumUnknown,
                       const MultiGridCL::TriangFaceIteratorCL& begin,
                       const MultiGridCL::TriangFaceIteratorCL& end,
                       const BndDataT& Bnd)
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence. Simplices on Dirichlet boundaries are skipped.
{
    if (NumUnknown == 0) return;
    for ( MultiGridCL::TriangFaceIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {        
            it->Unknowns.Prepare( idx, NumUnknown);
            for (Uint i=0; i<NumUnknown; ++i) it->Unknowns(idx)[i]= counter++;
        }
    }
}

void CreateNumbOnTetra( const Uint idx, IdxT& counter, Uint NumUnknown,
                        const MultiGridCL::TriangTetraIteratorCL& begin,
                        const MultiGridCL::TriangTetraIteratorCL& end);
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence.

template <class Iter>
inline void
DeleteNumbOnSimplex(Uint idx, const Iter& begin, const Iter& end)
    // deletes the memory for the Unknown-indices allocated by CreateNumbOn<Simplex>
{
    for (Iter it=begin; it!=end; ++it) 
        if (it->Unknowns.Exist() && it->Unknowns.Exist(idx) ) 
            it->Unknowns.Get()->Dealloc(idx);
}

} // end of namespace DROPS


#endif
