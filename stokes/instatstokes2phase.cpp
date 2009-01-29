/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#include "stokes/instatstokes2phase.h"

namespace DROPS
{

P1XRepairCL::P1XRepairCL (MultiGridCL& mg, VecDescCL& p)
    : UsesXFEM_( p.RowIdx->IsExtended()), mg_( mg), idx_( P1_FE), ext_( &idx_), p_( p)
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    const ExtIdxDescCL& extidx= p_.RowIdx->GetXidx();
//     size_t extbegin( extidx.GetNumUnknownsStdFE());
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p.RowIdx->GetIdx());

    idx_.CreateNumbering( p.RowIdx->TriangLevel(), mg);
    // Save the extended unknown-values.
    ext_.SetIdx( &idx_);

    // Attach the extended index to the vertex, so that it survives grid modifications.
    // We assume that all vertices have p-unknowns (like e. g. the pressure).
    DROPS_FOR_TRIANG_VERTEX( mg, p.RowIdx->TriangLevel(), it) {
        if ( (extunknown= extidx[it->Unknowns( pidx)]) != NoIdx ) {
            ext_.Data[ it->Unknowns( repairidx)]= p.Data[extunknown];
        }
        else{
            it->Unknowns( repairidx)= NoIdx;
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
        extunknown= extidx[it->Unknowns( pidx)];
        if ( extunknown!=NoIdx && it->Unknowns.Exist( repairidx) ) {
            p_.Data[extunknown]= ext_.Data[it->Unknowns( repairidx)];
            ++ci;
        }
    }
//     std::cerr << "P1XRepairCL::(): #P1-unknowns: " << extidx.GetNumUnknownsStdFE()
//               << "\t#copied extended-dof: " << ci << '\n';
}


void SetupMassDiag_P1(const MultiGridCL& MG, VectorCL& M, IdxDescCL& RowIdx, const BndDataCL<>& bnd)
{
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL Numb;

    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        Numb.assign( *sit, RowIdx, bnd);
        for(int i=0; i<4; ++i)
            if (Numb.WithUnknowns( i))
                M[Numb.num[i]]+= P1DiscCL::GetMass( i, i)*absdet;
    }
}

void SetupMassDiag_P1X (const MultiGridCL& MG, VectorCL& M, IdxDescCL& RowIdx, const LevelsetP2CL& lset,
                        const BndDataCL<>& bnd)
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL Numb;
    double coup[4], coupT2[4];

    double integralp;
    InterfacePatchCL cut;
    bool sign[4];

    // The 4 squares of the P1-shape-functions
    LocalP2CL<> pi2[4];
    for(int i= 0; i < 4; ++i) {
        pi2[i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
            pi2[i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }

    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        loc_phi.assign( *sit, lset.Phi, NoBndDataCL<> ());
        cut.Init( *sit, loc_phi);
        const bool nocut= !cut.Intersects();
        Numb.assign( *sit, RowIdx, bnd);
        if (nocut) {
            for(int i= 0; i < 4; ++i)
                if ( Numb.WithUnknowns( i))
                    M[Numb.num[i]]+= P1DiscCL::GetMass( i, i)*absdet;
        }
        else { // extended basis functions have only support on tetra intersecting Gamma!
            for(int i=0; i<4; ++i) {
                sign[i]= cut.GetSign(i) == 1;
                // compute the integrals
                // \int_{T_2} p_i^2 dx,    where T_2 = T \cap \Omega_2
                integralp= 0.;
                for (int ch= 0; ch < 8; ++ch) {
                    cut.ComputeCutForChild( ch);
                    integralp+= cut.quad( pi2[i], absdet, true);  // integrate on positive part
                }
                coup[i]= P1DiscCL::GetMass( i, i)*absdet;
                coupT2[i]= integralp;
            }

            // write values into matrix
            for(int i=0; i<4; ++i) {
                if (!Numb.WithUnknowns( i)) continue;
                M[Numb.num[i]]+= coup[i];
                const IdxT xidx_i= Xidx[Numb.num[i]];
                if (xidx_i!=NoIdx)
                    M[xidx_i]+= coupT2[i]*(1 - 2*sign[i]) + sign[i]*coup[i];
            }
        }
    }
}

} // end of namespace DROPS
