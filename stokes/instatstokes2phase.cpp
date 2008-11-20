/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#include "stokes/instatstokes2phase.h"

namespace DROPS
{

P1XRepairCL::P1XRepairCL (bool UsesXFEM, MultiGridCL& mg, VecDescCL& p,
    ExtIdxDescCL& extidx)
    : UsesXFEM_( UsesXFEM), mg_( mg), idx_( P1X_FE), extidx_( extidx), p_( p)
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    size_t extbegin( extidx.GetNumUnknownsP1());
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p.RowIdx->GetIdx());
    // Save the extended unknown-values.
    extData_.resize( p.RowIdx->NumUnknowns() - extbegin);
    extData_= p.Data[std::slice( extbegin, extData_.size(), 1)];

    // Attach the extended index to the vertex, so that it survives grid modifications.
    // We assume that all vertices have p-unknowns (like e. g. the pressure).
    DROPS_FOR_TRIANG_VERTEX( mg, p.RowIdx->TriangLevel(), it) {
        if ( (extunknown= extidx.Xidx[it->Unknowns( pidx)]) != NoIdx ) {
            it->Unknowns.Prepare( repairidx);
            it->Unknowns( repairidx)= extunknown - extbegin;
        }
    }
}

void P1XRepairCL::operator() (const LevelsetP2CL& lset)
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    // We assume that the caller has repaired p as a P1-FE-variable.
    // Thus we extend the index and save the P1-part first.
    VectorCL tmp( p_.Data);
    extidx_.UpdateXNumbering( p_.RowIdx, lset, /*NumberingChanged*/ true);
    p_.SetIdx( p_.RowIdx);
    p_.Data[std::slice( 0, tmp.size(), 1)]= tmp;

    size_t ci= 0;
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p_.RowIdx->GetIdx());
    // We assume that all vertices in p's level hold a p1-value. Thus, it->Unknowns.Exist()
    // can be spared.
    DROPS_FOR_TRIANG_VERTEX( mg_, p_.RowIdx->TriangLevel(), it) {
        if ( ((extunknown= extidx_.Xidx[it->Unknowns( pidx)]) != NoIdx)
            && it->Unknowns.Exist( repairidx) ) {
            p_.Data[extunknown]= extData_[it->Unknowns( repairidx)];
            ++ci;
        }
    }
    std::cerr << "P1XRepairCL::(): #P1-unknowns: " << extidx_.GetNumUnknownsP1()
              << "\t#copied extended-dof: " << ci << '\n';
}


} // end of namespace DROPS
