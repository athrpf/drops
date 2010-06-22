/// \file problem.cpp
/// \brief Classes that constitute a problem.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "misc/problem.h"
#include "num/interfacePatch.h"
#ifdef _PAR
#  include "parallel/interface.h"
#  include "parallel/exchange.h"
#endif

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

IdxDescCL::IdxDescCL( FiniteElementT fe, const BndCondCL& bnd, match_fun match, double omit_bound)
    : FE_InfoCL( fe), Idx_( GetFreeIdx()), TriangLevel_( 0), NumUnknowns_( 0), Bnd_(bnd), match_(match),
      extIdx_( omit_bound != -99 ? omit_bound : IsExtended() ? 1./32. : -1.) // default value is 1./32. for XFEM and -1 otherwise
{
#ifdef _PAR
    ex_= new ExchangeCL();
#endif
}

IdxDescCL::~IdxDescCL()
{
    if (Idx_!=InvalidIdx)
        IdxFree[Idx_]= true;
#ifdef _PAR
    delete ex_;
#endif
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
#ifdef _PAR
    ex_= new ExchangeCL(*orig.ex_);
#endif
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
#ifdef _PAR
    std::swap( ex_,          obj.ex_);
#endif
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
        std::cout << "Compare_Indices: Indices on different levels.\n";
        return false;
    }
    if (i.NumUnknowns() != j.NumUnknowns()) {
        std::cout << "Compare_Indices: NumUnknowns different.\n";
        return false;
    }
    if (i.GetFE() !=  j.GetFE()) {
        std::cout << "Compare_Indices: FE types different.\n";
        return false;
    }
    if (i.NumUnknownsVertex_ != j.NumUnknownsVertex_) {
        std::cout << "Compare_Indices: NumUnknownsVertex different.\n";
        return false;
    }
    if (i.NumUnknownsEdge_ != j.NumUnknownsEdge_) {
        std::cout << "Compare_Indices: NumUnknownsEdge different.\n";
        return false;
    }
    if (i.NumUnknownsFace_ != j.NumUnknownsFace_) {
        std::cout << "Compare_Indices: NumUnknownsFace different.\n";
        return false;
    }
    if (i.NumUnknownsTetra_ != j.NumUnknownsTetra_) {
        std::cout << "Compare_Indices: NumUnknownsTetra different.\n";
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
                    std::cout << "Compare_Indices: Vertex difference.\n";
                    return false;
                }
    if (i.NumUnknownsEdge_ != 0)
        for (MultiGridCL::const_TriangEdgeIteratorCL it= mg->GetTriangEdgeBegin( lvl),
             theend= mg->GetTriangEdgeEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Edge difference.\n";
                    return false;
                }
    if (i.NumUnknownsFace_ != 0)
        for (MultiGridCL::const_TriangFaceIteratorCL it= mg->GetTriangFaceBegin( lvl),
             theend= mg->GetTriangFaceEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Face difference.\n";
                    return false;
                }
    if ( i.NumUnknownsTetra_ != 0)
        for (MultiGridCL::const_TriangTetraIteratorCL it= mg->GetTriangTetraBegin( lvl),
             theend= mg->GetTriangTetraEnd( lvl); it != theend; ++it)
            if (it->Unknowns.Exist())
                if ( (it->Unknowns.Exist( iidx) != it->Unknowns.Exist( jidx)) ) {
                    std::cout << "Compare_Indices: Tetra difference.\n";
                    return false;
                }
    return true;
}


#ifdef _PAR
/// \brief Count number of exclusive unknowns on this proc
/** Get number of unknowns that are exclusive (master copy of simplex and proc with smalles proc-id)*/
IdxT IdxDescCL::GetExclusiveNumUnknowns(const MultiGridCL &mg, int lvl) const
{
    IdxT ret=0;
    // If there are unknowns on vertices, the proc with the smallest proc-id
    // who owns a master copy of the vertex counts the unknowns of the chosen
    // vertex.
    if (NumUnknownsVertex())
        for (MultiGridCL::const_TriangVertexIteratorCL it(mg.GetTriangVertexBegin(lvl)), end(mg.GetTriangVertexEnd(lvl)); it != end; ++it)
            if ( it->Unknowns.Exist() && it->Unknowns.Exist(GetIdx()))
                if (it->IsExclusive())
                    ret += NumUnknownsVertex();
    if (NumUnknownsEdge())
        for (MultiGridCL::const_TriangEdgeIteratorCL it(mg.GetTriangEdgeBegin(lvl)), end(mg.GetTriangEdgeEnd(lvl)); it != end; ++it)
            if ( it->Unknowns.Exist() && it->Unknowns.Exist(GetIdx()))
                if (it->IsExclusive())
                    ret += NumUnknownsEdge();
    if (NumUnknownsTetra())
        for (MultiGridCL::const_TriangTetraIteratorCL it(mg.GetTriangTetraBegin(lvl)), end(mg.GetTriangTetraEnd(lvl)); it != end; ++it)
            if ( it->Unknowns.Exist() && it->Unknowns.Exist(GetIdx()))
                if (it->IsExclusive())
                    ret += NumUnknownsTetra();
    return ret;
}

/// \brief Count global number of unknowns
/** Get global number over all procs of unknowns. Each unknown is just count one time*/
IdxT IdxDescCL::GetGlobalNumUnknowns(const MultiGridCL &mg, int lvl) const
{
    return ProcCL::GlobalSum(GetExclusiveNumUnknowns(mg,lvl));
}
#endif

void P1XtoP1 (const IdxDescCL& xidx, const VectorCL& p1x, const IdxDescCL& idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg)
{
    const Uint lvl= idx.TriangLevel(),
               p1idxnum= idx.GetIdx(),
                 idxnum= xidx.GetIdx(),
               lsidxnum= lset.RowIdx->GetIdx();
    const ExtIdxDescCL& extIdx= xidx.GetXidx();
    const size_t p1unknowns = extIdx.GetNumUnknownsStdFE();
    if (p1unknowns != idx.NumUnknowns())
        throw DROPSErrCL( "P1XtoP1: inconsistent indices\n");

    negPart.resize(p1unknowns);
    posPart.resize(p1unknowns);
    posPart = negPart = p1x[std::slice(0, p1unknowns, 1)];

    // add extended pressure
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        const IdxT   nr= it->Unknowns( idxnum);
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT p1nr= it->Unknowns( p1idxnum);
        if (extIdx[nr]==NoIdx) continue;

        if (InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)]) == 1)
            negPart[p1nr]= p1x[nr] - p1x[extIdx[nr]];
        else
            posPart[p1nr]= p1x[nr] + p1x[extIdx[nr]];
    }
}

void P1toP1X (const IdxDescCL& xidx, VectorCL& p1x, const IdxDescCL& idx, const VectorCL& posPart, const VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg)
{
    const Uint lvl= idx.TriangLevel(),
                idxnum= xidx.GetIdx(),
                p1idxnum= idx.GetIdx(),
                lsidxnum= lset.RowIdx->GetIdx();

    p1x.resize(xidx.NumUnknowns());
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT p1nr= it->Unknowns(p1idxnum);
        const bool is_pos= InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)])==1;
        if (is_pos)
            p1x[nr]= posPart[p1nr];
        else
            p1x[nr]= negPart[p1nr];
    }
    // add extended pressure
    const ExtIdxDescCL& extIdx= xidx.GetXidx();
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (!it->Unknowns.Exist( idxnum)) continue;
        const IdxT p1nr= it->Unknowns(p1idxnum);
        if (extIdx[nr]==NoIdx) continue;
        p1x[extIdx[nr]]= (posPart[p1nr] - negPart[p1nr]);
    }
}

void ExtractComponent( const VectorCL& vecFE, VectorCL& scalarFE, Uint comp, Uint stride)
{
    Assert( vecFE.size()==scalarFE.size()*stride, DROPSErrCL("ExtractComponent: vector sizes do not match"), DebugNumericC);
    for (size_t i=0, s=scalarFE.size(); i<s; ++i)
        scalarFE[i]= vecFE[i*stride+comp];
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

/// \brief Routine to number unknowns on the vertices surrounding an
/// interface.
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are cut by the zero level of lset. If InterfacePatchCL::IntersectsInterior() is true,
/// all vertices are numbered, if only InterfacePatchCL::Intersects() is true, only the
/// vertices with InterfacePatchCL::GetSign(vertex) == 0 are numbered. The other vertices in
/// such tetrahedra obtain NoIdx as number, but they are not counted as unknowns.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// A more user friendly interface is provided by IdxDescCL::CreateNumbOnInterface.
void CreateNumbOnInterfaceVertex (const Uint idx, IdxT& counter, Uint stride,
        const MultiGridCL::TriangVertexIteratorCL& vbegin,
        const MultiGridCL::TriangVertexIteratorCL& vend,
        const MultiGridCL::TriangTetraIteratorCL& begin,
        const MultiGridCL::TriangTetraIteratorCL& end,
    const VecDescCL& ls, double omit_bound= -1./*default to using all dof*/)
{
    if (stride == 0) return;

    LocalP2CL<> hat_sq[4]; // values of phi_i*phi_i
    for (int i= 0; i < 4; ++i)  {
        hat_sq[i][i]=1.;
        for (int j = 0; j < 4; ++j)
            if (i != j) hat_sq[i][EdgeByVert(i,j)+4]= 0.25;
    }
    // first set NoIdx in all vertices
    for (MultiGridCL::TriangVertexIteratorCL vit= vbegin; vit != vend; ++vit) {
        vit->Unknowns.Prepare(idx);
        vit->Unknowns.Invalidate(idx);
    }
    // then create numbering of vertices at the interface
    InterfaceTriangleCL p;
    for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
        p.Init( *it, ls);
        if (!p.Intersects()) continue;

        const double h3= it->GetVolume()*6, h= cbrt( h3), h4= h*h3, limit= h4*omit_bound;
        SVectorCL<4> loc_int; // stores integrals \int_{\Gamma_T} p^2 dx, with p1-dof p.
        for (int ch= 0; ch < 8; ++ch) {
            if (!p.ComputeForChild( ch)) continue;// no patch for this child
            for (int tri= 0; tri < p.GetNumTriangles(); ++tri) {
                for (int i= 0; i < 4; ++i) {
                    loc_int[i]+= p.quad2D( hat_sq[i], tri);
                }
            }
        }

        const bool innercut( p.IntersectsInterior());
        for (Uint i= 0; i < NumVertsC; ++i) {
            UnknownHandleCL& u= const_cast<VertexCL*>( it->GetVertex( i))->Unknowns;
            if (innercut || p.GetSign( i) == 0) {
                if ( u( idx) == NoIdx) {
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance
                    u( idx)= counter;
                    counter+= stride;
                }
            }
        }
    }
}

void IdxDescCL::CreateNumbOnInterface(Uint level, MultiGridCL& mg, const VecDescCL& ls, double omit_bound)
/// Uses CreateNumbOnInterfaceVertex on the triangulation with level \p level on the multigrid \p mg.
/// One can only create P1-elements.
{
    // set up the index description
    const Uint idxnum= GetIdx();
    TriangLevel_= level;
    NumUnknowns_= 0;

    // allocate space for indices; number unknowns in TriangLevel level
    if (NumUnknownsVertex() != 0)
        CreateNumbOnInterfaceVertex( idxnum, NumUnknowns_, NumUnknownsVertex(),
            mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level),
            mg.GetTriangTetraBegin( level), mg.GetTriangTetraEnd( level), ls, omit_bound);

    if (NumUnknownsEdge() != 0 || NumUnknownsFace() != 0 || NumUnknownsTetra() != 0)
        throw DROPSErrCL( "CreateNumbOnInterface: Only vertex unknowns are implemented\n" );
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

void IdxDescCL::CreateNumbering( Uint level, MultiGridCL& mg, const VecDescCL* lsetp)
/// Memory for the Unknown-Indices on TriangLevel level is allocated
/// and the unknowns are numbered.
/// If a matching function is specified, numbering on periodic boundaries
/// is performed, too.
/// After that the extended DoFs are numbered for extended FE.
{
    if (IsOnInterface())
    {
#ifdef _PAR
        throw DROPSErrCL("IdxDescCL::CreateNumbering: Check first, if numbering on interface works in parDROPS.");
#endif
        if (lsetp == 0) throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for interface numbering given");
        CreateNumbOnInterface( level, mg, *lsetp, GetXidx().GetBound());
    }
    else {
        CreateNumbStdFE( level, mg);
        if (IsExtended()) {
            if (lsetp == 0) throw DROPSErrCL("IdxDescCL::CreateNumbering: no level set function for XFEM numbering given");
            NumUnknowns_= extIdx_.UpdateXNumbering( this, mg, *lsetp, true);
        }
    }
#ifdef _PAR
    ex_->CreateList(mg, this, true, true);
#endif
}

void IdxDescCL::UpdateXNumbering( MultiGridCL& mg, const VecDescCL& lset)
{
    if (IsExtended()) {
        NumUnknowns_= extIdx_.UpdateXNumbering( this, mg, lset, false);
#ifdef _PAR
        ex_->CreateList(mg, this, true, true);
#endif
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
#ifdef _PAR
    ex_->clear();
#endif
}

IdxT ExtIdxDescCL::UpdateXNumbering( IdxDescCL* Idx, const MultiGridCL& mg, const VecDescCL& lset, bool NumberingChanged)
{
    const Uint sysnum= Idx->GetIdx(),
        level= Idx->TriangLevel(),
        stride= Idx->IsScalar() ? 1: 3; // scalar or vector-valued FE
    IdxT extIdx= NumberingChanged ? Idx->NumUnknowns() : Xidx_.size();
    Xidx_old_.assign( extIdx, NoIdx);
    Xidx_.swap( Xidx_old_);
    lset_= &lset;
    InterfaceTetraCL cut;
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
                    if (!it->GetVertex(i)->Unknowns.Exist( sysnum)) continue;
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx_[nr]==NoIdx) // not extended yet
                        for (Uint k=0; k<stride; ++k)
                            Xidx_[nr+k]= extIdx++;
                }
            else
                for (Uint i=0; i<4; ++i)
                { // extend only DoFs on interface
                    if (cut.GetSign(i)!=0) continue;
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance, which lead to unstable solutions
                    if (!it->GetVertex(i)->Unknowns.Exist( sysnum)) continue;
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx_[nr]==NoIdx) // not extended yet
                        for (Uint k=0; k<stride; ++k)
                            Xidx_[nr+k]= extIdx++;
                }
        }
    }
#ifdef _PAR
    // communicate extended dofs on vertices
    current_Idx_= Idx;
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<VertexCL>::GetIF(),  // exchange datas over distributed vertices
                   sizeof(bool),                    // number of datas to be exchanged
                   HandlerGatherUpdateXNumbC,       // how to gather datas
                   HandlerScatterUpdateXNumbC       // how to scatter datas
                  );
    current_Idx_= 0;

    // number all extended dofs from other procs (where extended dof is flagged by NoIdx-1)
    for (size_t i=0; i<Xidx_.size(); ++i)
        if (Xidx_[i] == NoIdx-1)
            Xidx_[i]= extIdx++;
#endif
    return extIdx;
}

#ifdef _PAR
IdxDescCL* ExtIdxDescCL::current_Idx_= 0;

int ExtIdxDescCL::HandlerGatherUpdateXNumb( OBJT objp, void* buf)
{
    VertexCL* const sp= ddd_cast<VertexCL*>(objp);
    bool* buffer= static_cast<bool*>(buf);
    if (sp->Unknowns.Exist(current_Idx_->GetIdx()))
        *buffer= current_Idx_->IsExtended( sp->Unknowns(current_Idx_->GetIdx()));
    else
        *buffer= false;
    return 0;
}

int ExtIdxDescCL::HandlerScatterUpdateXNumb( OBJT objp, void* buf)
{
    VertexCL* const sp= ddd_cast<VertexCL*>(objp);
    bool RemoteExtended= *static_cast<bool*>(buf);
    if (!sp->Unknowns.Exist(current_Idx_->GetIdx()))
        return 0;
    const IdxT dof= sp->Unknowns(current_Idx_->GetIdx());

    if (!current_Idx_->IsExtended( dof) && RemoteExtended)
        current_Idx_->GetXidx()[dof]= NoIdx-1;
    return 0;
}
#endif

void ExtIdxDescCL::Old2New(VecDescCL* v)
{
    VectorCL tmp( v->Data);
    // set v to zero vector with appropriate length
    v->SetIdx( v->RowIdx);
    // copy standard FE part
    const IdxT nStd= GetNumUnknownsStdFE();
#if defined(__SUNPRO_CC) || defined(DROPS_WIN)
    for (size_t i=0; i<nStd; ++i)
        v->Data[i] = tmp[i];
#else
    v->Data[std::slice( 0, nStd, 1)]= tmp[std::slice( 0, nStd, 1)];
#endif

    if (Xidx_.size() != Xidx_old_.size()) { // standard FE index changed (e.g., grid changed)
#ifndef _PAR
          std::cout << "ExtIdxDescCL::Old2New: Xidx: " << Xidx_.size()
                    << "\tXidx_old: " << Xidx_old_.size()
                    << "Extended Unknowns set to 0.\n";
#endif
        return;
    }

    // treat extended part
    IdxT ni= 0, di=0, ri= 0, ci= 0;
    for (size_t i= 0; i < Xidx_.size(); ++i) {
        if ( Xidx_old_[i] == NoIdx && Xidx_[i] != NoIdx) ++ni;
        if ( Xidx_old_[i] != NoIdx && Xidx_[i] == NoIdx) ++di;
        if ( Xidx_old_[i] != NoIdx && Xidx_[i] != NoIdx) { // extended dof was also extended before
            if ( Xidx_old_[i] != Xidx_[i])
                ++ri;
            v->Data[Xidx_[i]]= tmp[Xidx_old_[i]];
            ++ci;
        }
    }

#ifndef _PAR
    std::cout << "ExtIdxDescCL::Old2New: #P1-unknowns: " <<Xidx_.size()
              << "\t#new dof: " << ni
              << "\t#deleted dof: " << di
              << "\t#renumbered dof: " << ri
              << "\t#copied extended-dof: " << ci
              << '\n';
#endif
}

} // end of namespace DROPS
