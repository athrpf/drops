/// \file
/// \brief Classes that constitute a problem.

#include "misc/problem.h"

namespace DROPS
{

const Uint        IdxDescCL::InvalidIdx = std::numeric_limits<Uint>::max();
std::vector<bool> IdxDescCL::IdxFree;

void BndCondInfo( BndCondT bc, std::ostream& os)
/// \param bc Value of type BndCondT, which shall be described.
/// \param os Stream, to which the description is written.
{
    switch(bc)
    {
      case DirBC:  os << "inhom. Dirichlet BC / inflow\n"; break;
      case Dir0BC: os << "hom. Dirichlet BC / wall\n"; break;
      case NatBC:  os << "inhom. Natural BC\n"; break;
      case Nat0BC: os << "hom. Natural BC / outflow\n"; break;
      case Per1BC: os << "periodic BC\n"; break;
      case Per2BC: os << "periodic BC, correspondent\n"; break;
      default:     os << "WARNING! unknown BC\n";
    }
}


//void FE_InfoCL::Set( Uint unkVertex, Uint unkEdge, Uint unkFace, Uint unkTetra)
//{
//    NumUnknownsVertex_ = unkVertex;
//    NumUnknownsEdge_   = unkEdge;
//    NumUnknownsFace_   = unkFace;
//    NumUnknownsTetra_  = unkTetra;
//}

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
 : FE_InfoCL(orig), Idx_(orig.Idx_), Bnd_(orig.Bnd_), match_(orig.match_),
   TriangLevel_(orig.TriangLevel_), NumUnknowns_(orig.NumUnknowns_)
{
    // invalidate orig
    const_cast<IdxDescCL&>(orig).Idx_= InvalidIdx;
}

void IdxDescCL::swap( IdxDescCL& obj)
/// Note, that std::swap cannot be used for IdxDescCL-objects as the
/// assignment operator is not implemented.
{
    Assert( GetFE()==obj.GetFE(), DROPSErrCL("IdxDescCL::swap: FE-types differ"), ~0);
	std::swap( Idx_, obj.Idx_);
    std::swap( TriangLevel_,        obj.TriangLevel_);
    std::swap( NumUnknowns_,        obj.NumUnknowns_);
}

bool IdxDescCL::Equal(IdxDescCL& i, IdxDescCL& j, const MultiGridCL* mg)
/// \param i The left IdxDescCL.
/// \param j The left IdxDescCL.
/// \param mg Optional pointer to a multigrid. If it is given the numbers
///     on the simplices are compared, too. This is rather expensive and
///     only needed for some correctness tests.
/// \return If mg==0: True, iff all members of i and j have the same value.
///     If mg!=0: True, iff all members of i and j have the same value and
///     all numbers on the simplices of the given trianbgulation are equal.
{
    const Uint lvl= i.TriangLevel();
    if (lvl != j.TriangLevel()) {
        std::cerr << "Compare_Indices: Indices on different levels.\n";
        return false;
    }
    if (i.NumUnknowns() != j.NumUnknowns()) {
        std::cerr << "Compare_Indices: NumUnknowns different.\n";
        return false;
    }
    if (i.GetFE() !=  j.GetFE()) {
        std::cerr << "Compare_Indices: FE types different.\n";
        return false;
    }
    if (i.NumUnknownsVertex_ != j.NumUnknownsVertex_) {
        std::cerr << "Compare_Indices: NumUnknownsVertex different.\n";
        return false;
    }
    if (i.NumUnknownsEdge_ != j.NumUnknownsEdge_) {
        std::cerr << "Compare_Indices: NumUnknownsEdge different.\n";
        return false;
    }
    if (i.NumUnknownsFace_ != j.NumUnknownsFace_) {
        std::cerr << "Compare_Indices: NumUnknownsFace different.\n";
        return false;
    }
    if (i.NumUnknownsTetra_ != j.NumUnknownsTetra_) {
        std::cerr << "Compare_Indices: NumUnknownsTetra different.\n";
        return false;
    }
    if (!mg)
        return true;

    const Uint iidx= i.GetIdx(),
               jidx= j.GetIdx();
    if (iidx == jidx)
        return true;
    if ( i.NumUnknownsVertex_ != 0)
        for (MultiGridCL::const_TriangVertexIteratorCL it= mg->GetTriangVertexBegin( lvl),
             theend= mg->GetTriangVertexEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cerr << "Compare_Indices: Vertex difference.\n";
                    return false;
                }
    if (i.NumUnknownsEdge_ != 0)
        for (MultiGridCL::const_TriangEdgeIteratorCL it= mg->GetTriangEdgeBegin( lvl),
             theend= mg->GetTriangEdgeEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cerr << "Compare_Indices: Edge difference.\n";
                    return false;
                }
    if (i.NumUnknownsFace_ != 0)
        for (MultiGridCL::const_TriangFaceIteratorCL it= mg->GetTriangFaceBegin( lvl),
             theend= mg->GetTriangFaceEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cerr << "Compare_Indices: Face difference.\n";
                    return false;
                }
    if ( i.NumUnknownsTetra_ != 0)
        for (MultiGridCL::const_TriangTetraIteratorCL it= mg->GetTriangTetraBegin( lvl),
             theend= mg->GetTriangTetraEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cerr << "Compare_Indices: Tetra difference.\n";
                    return false;
                }
    return true;
}

void CreateNumbOnTetra( const Uint idx, IdxT& counter, Uint stride,
                        const MultiGridCL::TriangTetraIteratorCL& begin,
                        const MultiGridCL::TriangTetraIteratorCL& end)
{
    if (stride == 0) return;
    for (MultiGridCL::TriangTetraIteratorCL it=begin; it!=end; ++it)
    {
        it->Unknowns.Prepare( idx);
        it->Unknowns(idx)= counter;
        counter+= stride;
    }
}

void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg)
/// Memory for the Unknown-Indices on TriangLevel level is allocated
/// and the unknowns are numbered.
/// If a matching function is specified, numbering on periodic boundaries
/// is performed, too.
{
    const Uint idxnum= GetIdx();
    TriangLevel_= level; 
    NumUnknowns_ = 0;

    // allocate space for indices; number unknowns in TriangLevel level
    if (match_)
    {
        if (NumUnknownsVertex())
            CreatePeriodicNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsVertex(), match_,
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd_);
        if (NumUnknownsEdge())
            CreatePeriodicNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsEdge(), match_,
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd_);
        if (NumUnknownsFace())
            CreatePeriodicNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsFace(), match_,
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd_);
        if (NumUnknownsTetra())
            CreateNumbOnTetra( idxnum, NumUnknowns_, NumUnknownsTetra(),
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level));
    }
    else
    {
        if (NumUnknownsVertex())
            CreateNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsVertex(),
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd_);
        if (NumUnknownsEdge())
            CreateNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsEdge(),
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd_);
        if (NumUnknownsFace())
            CreateNumbOnSimplex( idxnum, NumUnknowns_, NumUnknownsFace(),
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd_);
        if (NumUnknownsTetra())
            CreateNumbOnTetra( idxnum, NumUnknowns_, NumUnknownsTetra(),
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level));
    }
}

void IdxDescCL::DeleteNumbering(MultiGridCL& MG)
/// This routine writes NoIdx as unknown-index for all indices of the
/// given index-description. NumUnknowns will be set to zero.
{
    const Uint idxnum = GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = TriangLevel_;
    NumUnknowns_ = 0;

    // delete memory allocated for indices
    if (NumUnknownsVertex())
        DeleteNumbOnSimplex( idxnum, MG.GetAllVertexBegin(level), MG.GetAllVertexEnd(level) );
    if (NumUnknownsEdge())
        DeleteNumbOnSimplex( idxnum, MG.GetAllEdgeBegin(level), MG.GetAllEdgeEnd(level) );
    if (NumUnknownsFace())
        DeleteNumbOnSimplex( idxnum, MG.GetAllFaceBegin(level), MG.GetAllFaceEnd(level) );
    if (NumUnknownsTetra())
        DeleteNumbOnSimplex( idxnum, MG.GetAllTetraBegin(level), MG.GetAllTetraEnd(level) );
}

} // end of namespace DROPS
