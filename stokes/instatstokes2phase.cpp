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
    std::swap( Xidx, Xidx_old);
    InterfacePatchCL cut;

    DROPS_FOR_TRIANG_CONST_TETRA( lset.GetMG(), level, it)
    {
//const double absdet= it->GetVolume()*6.,
//min_ratio=1./64;
        cut.Init( *it, lset.Phi);
        if (cut.Intersects() )
        {
/*
// compute volume of part
double vol= 0;

for (int ch=0; ch<8; ++ch)
{
    cut.ComputeCutForChild(ch);
    vol+= cut.quad( ones, absdet, false); // integrate on other part
}
double ratio[2];
ratio[0]= std::abs( vol)/absdet*6.;
ratio[1]= 1-ratio[0];
*/
            if (cut.IntersectsInterior())
                for (Uint i=0; i<4; ++i)
                { // extend all DoFs
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
//if (ratio[cut.GetSign(i)==1] < min_ratio) { std::cerr << "    >>> Omitted DoF (sign="<<cut.GetSign(i)<<") with ratio = " << ratio[cut.GetSign(i)==1] << std::endl; continue; }
                    if (Xidx[nr]==NoIdx) // not extended yet
                        Xidx[nr]= extIdx++;
                }
            else
                for (Uint i=0; i<4; ++i)
                { // xtend only DoFs on interface
                    if (cut.GetSign(i)!=0) continue;
                    const IdxT nr= it->GetVertex(i)->Unknowns(sysnum);
//if (ratio[cut.GetSign(i)==1] < min_ratio) { std::cerr << "    >>> Omitted DoF (sign="<<cut.GetSign(i)<<") with ratio = " << ratio[cut.GetSign(i)==1] << std::endl; continue; }
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
