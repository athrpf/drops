/// \file
/// \brief Classes that constitute a problem.

#include "misc/problem.h"
#include "num/interfacePatch.h"

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
 : FE_InfoCL(orig), Idx_(orig.Idx_), TriangLevel_(orig.TriangLevel_), NumUnknowns_(orig.NumUnknowns_),
   Bnd_(orig.Bnd_), match_(orig.match_), extIdx_(orig.extIdx_)

{
    // invalidate orig
    const_cast<IdxDescCL&>(orig).Idx_= InvalidIdx;
}

void IdxDescCL::swap( IdxDescCL& obj)
/// Note, that std::swap cannot be used for IdxDescCL-objects as the
/// assignment operator is not implemented.
{
    Assert( GetFE()==obj.GetFE(), DROPSErrCL("IdxDescCL::swap: FE-types differ"), ~0);
	std::swap( Idx_,         obj.Idx_);
    std::swap( TriangLevel_, obj.TriangLevel_);
    std::swap( NumUnknowns_, obj.NumUnknowns_);
    std::swap( Bnd_,         obj.Bnd_);
    std::swap( match_,       obj.match_);
    std::swap( extIdx_,      obj.extIdx_);
}

bool IdxDescCL::Equal(IdxDescCL& i, IdxDescCL& j, const MultiGridCL* mg)
/// \param i The left IdxDescCL.
/// \param j The right IdxDescCL.
/// \param mg Optional pointer to a multigrid. If it is given the numbers
///     on the simplices are compared, too. This is rather expensive and
///     only needed for some correctness tests.
/// \return If mg==0: True, iff all members of i and j have the same value.
///     If mg!=0: True, iff all members of i and j have the same value and
///     all numbers on the simplices of the given triangulation are equal.
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

void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const VecDescCL* lsetp)
/// Memory for the Unknown-Indices on TriangLevel level is allocated
/// and the unknowns are numbered.
/// If a matching function is specified, numbering on periodic boundaries
/// is performed, too.
/// After that the extended DoFs are numbered for extended FE.
{
	CreateNumbStdFE( level, mg);
	if (IsExtended()) {
		if (lsetp == 0) throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for XFEM numbering given");
		NumUnknowns_= extIdx_.UpdateXNumbering( this, mg, *lsetp, true);
	}
}

void IdxDescCL::CreateNumbStdFE( Uint level, MultiGridCL& mg)
// numbering of standard FE
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
    extIdx_.DeleteXNumbering();
}

IdxT ExtIdxDescCL::UpdateXNumbering( const IdxDescCL* Idx, const MultiGridCL& mg, const VecDescCL& lset, bool NumberingChanged)
{
    const Uint sysnum= Idx->GetIdx(),
        level= Idx->TriangLevel();
    IdxT extIdx= NumberingChanged ? Idx->NumUnknowns() : Xidx_.size();
    Xidx_old_.assign( extIdx, NoIdx);
    Xidx_.swap( Xidx_old_);
    InterfacePatchCL cut;

    LocalP2CL<> hat_sq[4]; // values of phi_i*phi_i

    for (int i=0; i<4; ++i) //initialize hat_ii
    {
        hat_sq[i][i]=1.;
        for (int j=0; j<4; ++j)
            if (i!=j)
                hat_sq[i][EdgeByVert(i,j)+4]=0.25;
    }
    LocalP2CL<> locPhi;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, level, it)
    {
        const double h3= it->GetVolume()*6,
            h= cbrt( h3), h5= h*h*h3, // h^5
            limit= h5*omit_bound_;
        locPhi.assign( *it, lset, NoBndDataCL<>());
        cut.Init( *it, locPhi);
        SVectorCL<4> loc_int; // stores integrals \int_T p^2 dx, where p = p_i^\Gamma. Extended DoF is omitted
                              // iff this integral falls below a certain bound (omit_bound*h^5) for all tetrahedra T.
        if (cut.Intersects() )
        {

            for (int ch=0; ch<8; ++ch)
            {
                cut.ComputeCutForChild(ch);
                for (Uint i=0; i<4; ++i)
                    loc_int[i]+= cut.quad( hat_sq[i], h3, cut.GetSign(i) != 1); // integrate on other part
            }
            if (cut.IntersectsInterior())
                for (Uint i=0; i<4; ++i)
                { // extend all DoFs
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance, which lead to unstable solutions
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx_[nr]==NoIdx) // not extended yet
                        Xidx_[nr]= extIdx++;
                }
            else
                for (Uint i=0; i<4; ++i)
                { // extend only DoFs on interface
                    if (cut.GetSign(i)!=0) continue;
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance, which lead to unstable solutions
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx_[nr]==NoIdx) // not extended yet
                        Xidx_[nr]= extIdx++;
                }
        }
    }
    return extIdx;
}

void ExtIdxDescCL::Old2New(VecDescCL* v)
{
    VectorCL tmp( v->Data);
    v->SetIdx( v->RowIdx);
    std::slice p1_part( 0, Xidx_.size(), 1);
    v->Data[p1_part]= tmp[p1_part];

    if (Xidx_.size() != Xidx_old_.size()) {
        std::cerr << "ExtIdxDescCL::Old2New: Xidx: " << Xidx_.size()
                  << "\tXidx_old: " << Xidx_old_.size()
                  << "Extended Unknowns set to 0.\n";
        return;
    }

    IdxT ni= 0, di=0, ri= 0, ci= 0;
    for (size_t i= 0; i < Xidx_.size(); ++i) {
        if ( Xidx_old_[i] == NoIdx && Xidx_[i] != NoIdx) ++ni;
        if ( Xidx_old_[i] != NoIdx && Xidx_[i] == NoIdx) ++di;
        if ( Xidx_old_[i] != NoIdx && Xidx_[i] != NoIdx) {
            if ( Xidx_old_[i] != Xidx_[i])
                ++ri;
            v->Data[Xidx_[i]]= tmp[Xidx_old_[i]];
            ++ci;
        }
    }
    std::cerr << "ExtIdxDescCL::Old2New: #P1-unknowns: " << Xidx_.size()
              << "\t#new dof: " << ni
              << "\t#deleted dof: " << di
              << "\t#renumbered dof: " << ri
              << "\t#copied extended-dof: " << ci
              << '\n';
}

} // end of namespace DROPS
