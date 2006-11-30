/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#include "stokes/instatstokes2phase.h"

namespace DROPS
{

void ExtIdxDescCL::UpdateXNumbering(IdxDescCL* idx, const LevelsetP2CL& lset, bool NumberingChanged)
/// Has to be called in two situations:
/// - when numbering of \p idx has changed, i.e. \p CreateNumbering was called before (set \p NumberingChanged=true)
/// - whenever level set function has changed to account for the moving interface (set \p NumberingChanged=false).
{
    if (idx != this->Idx) throw DROPSErrCL("ExtIdxDescCL::UpdateXNumbering: Wrong Index.\n");

    const Uint sysnum= idx->GetIdx(),
        level= idx->TriangLevel;
    IdxT extIdx= NumberingChanged ? idx->NumUnknowns : Xidx.size();
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
        
    DROPS_FOR_TRIANG_CONST_TETRA( lset.GetMG(), level, it)
    {
        const double h3= it->GetVolume()*6,
            h= cbrt( h3), h5= h*h*h3, // h^5
            limit= h5*omit_bound_;
        cut.Init( *it, lset.Phi);
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
    idx->NumUnknowns= extIdx;
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


} // end of namespace DROPS
