//**************************************************************************
// File:    problem.cpp                                                    *
// Content: classes that constitute a problem                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

#include "misc/problem.h"

namespace DROPS
{

std::vector<bool> IdxDescCL::IdxFree;

void BndCondInfo( BndCondT bc, std::ostream& os)
{
    switch(bc)
    {
      case DirBC:  os << "inhom. Dirichlet BC / inflow\n"; break; 
      case Dir0BC: os << "hom. Dirichlet BC / wall\n"; break; 
      case NatBC:  os << "inhom. Natural BC\n"; break; 
      case Nat0BC: os << "hom. Natural BC / outflow\n"; break; 
      default:     os << "WARNING! unknown BC\n";
    }
}


Uint IdxDescCL::GetFreeIdx()
{
    size_t sysnum= 0;
    for (; sysnum<IdxFree.size(); ++sysnum)
        if (IdxFree[sysnum]) break;
    if (sysnum>=IdxFree.size())
        IdxFree.push_back( false);
    else
        IdxFree[sysnum]= false;
    return sysnum;
}

IdxDescCL::IdxDescCL( const IdxDescCL& orig)
{
    Idx= orig.Idx;
    TriangLevel= orig.TriangLevel;
    NumUnknownsVertex= orig.NumUnknownsVertex;
    NumUnknownsEdge= orig.NumUnknownsEdge;
    NumUnknownsFace= orig.NumUnknownsFace;
    NumUnknownsTetra= orig.NumUnknownsTetra;
    NumUnknowns= orig.NumUnknowns;
    // invalidate orig
    const_cast<IdxDescCL&>(orig).Idx= NoIdx;
}

void IdxDescCL::swap( IdxDescCL& obj)
// Swaps the contents of obj and *this. Note that std::swap cannot be used
// for IdxDescCL-objects as the assignment operator is not implemented.
{
    std::swap( Idx, obj.Idx);
    std::swap( TriangLevel,       obj.TriangLevel);
    std::swap( NumUnknownsVertex, obj.NumUnknownsVertex);
    std::swap( NumUnknownsEdge,   obj.NumUnknownsEdge);
    std::swap( NumUnknownsFace,   obj.NumUnknownsFace);
    std::swap( NumUnknownsTetra,  obj.NumUnknownsTetra);
    std::swap( NumUnknowns,       obj.NumUnknowns);
}

void IdxDescCL::Set( Uint unkVertex, Uint unkEdge, Uint unkFace, Uint unkTetra)
// sets up the number of unknowns in every subsimplex for a certain index.
// Remark: expects _IdxDesc to be long enough
{
    NumUnknownsVertex = unkVertex;
    NumUnknownsEdge   = unkEdge;
    NumUnknownsFace   = unkFace;
    NumUnknownsTetra  = unkTetra;
}

void MatDescCL::SetIdx(IdxDescCL* row, IdxDescCL* col)
{
    // create new matrix and set up the matrix description
    RowIdx= row;
    ColIdx= col;
}

void MatDescCL::Reset()
{
    RowIdx = 0;
    ColIdx = 0;
    Data.clear();
}

void CreateNumbOnTetra( const Uint idx, IdxT& counter, Uint stride,
                        const MultiGridCL::TriangTetraIteratorCL& begin,
                        const MultiGridCL::TriangTetraIteratorCL& end)
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence.
{
    if (stride == 0) return;
    for (MultiGridCL::TriangTetraIteratorCL it=begin; it!=end; ++it)
    {
        it->Unknowns.Prepare( idx);
        it->Unknowns(idx)= counter;
        counter+= stride;
    }
}

} // end of namespace DROPS
