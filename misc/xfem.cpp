/// \file
/// \brief
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser, IGPM

#include "misc/xfem.h"

namespace DROPS
{

void ExtIdxDescCL::UpdateXNumbering(IdxDescCL* idx, const LevelsetP2CL& lset, bool NumberingChanged)
/// Has to be called in two situations:
/// - when numbering of \p idx has changed, i.e. \p CreateNumbering was called before (set \p NumberingChanged=true)
/// - whenever level set function has changed to account for the moving interface (set \p NumberingChanged=false).
{
    if (idx != this->Idx) throw DROPSErrCL("ExtIdxDescCL::UpdateXNumbering: Wrong Index.\n");

    const Uint sysnum= idx->GetIdx(),
        level= idx->TriangLevel();
    IdxT extIdx= NumberingChanged ? idx->NumUnknowns() : Xidx.size();
    Xidx_old.assign( extIdx, NoIdx);
    Xidx.swap( Xidx_old);
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

    DROPS_FOR_TRIANG_CONST_TETRA( lset.GetMG(), level, it)
    {
        const double h3= it->GetVolume()*6,
            h= cbrt( h3), h5= h*h*h3, // h^5
            limit= h5*omit_bound_;
        locPhi.assign( *it, lset.Phi, NoBndDataCL<>());
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
                    if (Xidx[nr]==NoIdx) // not extended yet
                        Xidx[nr]= extIdx++;
                }
            else
                for (Uint i=0; i<4; ++i)
                { // xtend only DoFs on interface
                    if (cut.GetSign(i)!=0) continue;
                    if (loc_int[i] < limit) continue; // omit DoFs of minor importance, which lead to unstable solutions
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
                    if (Xidx[nr]==NoIdx) // not extended yet
                        Xidx[nr]= extIdx++;
                }
        }
    }
    idx->SetNumUnknowns( extIdx);
}

void ExtIdxDescCL::Old2New(VecDescCL* v)
{
    if (v->RowIdx != this->Idx) throw DROPSErrCL("ExtIdxDescCL::Old2New: Wrong Index.\n");

    VectorCL tmp( v->Data);
    v->SetIdx( v->RowIdx);
    std::slice p1_part( 0, Xidx.size(), 1);
    v->Data[p1_part]= tmp[p1_part];

    if (Xidx.size() != Xidx_old.size()) {
        std::cerr << "ExtIdxDescCL::Old2New: Xidx: " << Xidx.size()
                  << "\tXidx_old: " << Xidx_old.size()
                  << "Extended Unknowns set to 0.\n";
        return;
    }

    IdxT ni= 0, di=0, ri= 0, ci= 0;
    for (size_t i= 0; i < Xidx.size(); ++i) {
        if ( Xidx_old[i] == NoIdx && Xidx[i] != NoIdx) ++ni;
        if ( Xidx_old[i] != NoIdx && Xidx[i] == NoIdx) ++di;
        if ( Xidx_old[i] != NoIdx && Xidx[i] != NoIdx) {
            if (Xidx_old[i] != Xidx[i])
                ++ri;
            v->Data[Xidx[i]]= tmp[Xidx_old[i]];
            ++ci;
        }
    }
    std::cerr << "ExtIdxDescCL::Old2New: #P1-unknowns: " << Xidx.size()
              << "\t#new dof: " << ni
              << "\t#deleted dof: " << di
              << "\t#renumbered dof: " << ri
              << "\t#copied extended-dof: " << ci
              << '\n';
}

void P1XtoP1 (const ExtIdxDescCL& xidx, const VectorCL& p1x, const IdxDescCL& idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg)
{
    const Uint lvl= idx.TriangLevel(),
                idxnum= xidx.Idx->GetIdx(),
                lsidxnum= lset.RowIdx->GetIdx();
    const size_t p1unknowns = xidx.GetNumUnknownsP1();
    negPart.resize(p1unknowns);
    posPart.resize(p1unknowns);

    posPart = negPart = p1x[std::slice(0, p1unknowns, 1)];

    // add extended pressure
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (xidx[nr]==NoIdx) continue;

        const bool is_pos= InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)])==1;
        if (is_pos)
            negPart[nr]-= p1x[xidx[nr]];
        else
            posPart[nr]+= p1x[xidx[nr]];
    }
}

void P1toP1X (const ExtIdxDescCL& xidx, VectorCL& p1x, const IdxDescCL& idx, const VectorCL& posPart, const VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg)
{
    const Uint lvl= idx.TriangLevel(),
                idxnum= xidx.Idx->GetIdx(),
                p1idxnum= idx.GetIdx(),
                lsidxnum= lset.RowIdx->GetIdx();

    p1x.resize(xidx.Idx->NumUnknowns());
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        const IdxT p1nr= it->Unknowns(p1idxnum);
        const bool is_pos= InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)])==1;
        if (is_pos)
            p1x[nr]= posPart[p1nr];
        else
            p1x[nr]= negPart[p1nr];
    }
    // add extended pressure
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        const IdxT p1nr= it->Unknowns(p1idxnum);
        if (xidx[nr]==NoIdx) continue;
        p1x[xidx[nr]]+= (posPart[p1nr] - negPart[p1nr]);
    }
}

} // end of namespace DROPS

