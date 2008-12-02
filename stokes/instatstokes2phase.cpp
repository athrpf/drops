/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#include "stokes/instatstokes2phase.h"

namespace DROPS
{

P1XRepairCL::P1XRepairCL (MultiGridCL& mg, VecDescCL& p)
    : UsesXFEM_( p.RowIdx->IsExtended()), mg_( mg), idx_( p.RowIdx->GetFE()), p_( p)
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    const ExtIdxDescCL& extidx= p_.RowIdx->GetXidx();
    size_t extbegin( extidx.GetNumUnknownsStdFE());
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p.RowIdx->GetIdx());
    // Save the extended unknown-values.
    extData_.resize( p.RowIdx->NumUnknowns() - extbegin);
    extData_= p.Data[std::slice( extbegin, extData_.size(), 1)];

    // Attach the extended index to the vertex, so that it survives grid modifications.
    // We assume that all vertices have p-unknowns (like e. g. the pressure).
    DROPS_FOR_TRIANG_VERTEX( mg, p.RowIdx->TriangLevel(), it) {
        if ( (extunknown= extidx[it->Unknowns( pidx)]) != NoIdx ) {
            it->Unknowns.Prepare( repairidx);
            it->Unknowns( repairidx)= extunknown - extbegin;
        }
    }
}

void P1XRepairCL::operator() ()
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    // We assume that the caller has repaired p as a P1-FE-variable.
    // Thus we only repair the extended part of p.
    size_t ci= 0;
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p_.RowIdx->GetIdx());
    const ExtIdxDescCL& extidx= p_.RowIdx->GetXidx();
    // We assume that all vertices in p's level hold a p1-value. Thus, it->Unknowns.Exist()
    // can be spared.
    DROPS_FOR_TRIANG_VERTEX( mg_, p_.RowIdx->TriangLevel(), it) {
        if ( ((extunknown= extidx[it->Unknowns( pidx)]) != NoIdx)
            && it->Unknowns.Exist( repairidx) ) {
            p_.Data[extunknown]= extData_[it->Unknowns( repairidx)];
            ++ci;
        }
    }
    std::cerr << "P1XRepairCL::(): #P1-unknowns: " << extidx.GetNumUnknownsStdFE()
              << "\t#copied extended-dof: " << ci << '\n';
}


} // end of namespace DROPS
