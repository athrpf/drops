/// \file
/// \brief
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser, IGPM

#include "misc/xfem.h"
#include "levelset/levelset.h"

namespace DROPS
{

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
        const IdxT p1nr= it->Unknowns( p1idxnum);
        if (extIdx[nr]==NoIdx) continue;

        if (InterfacePatchCL::Sign( lset.Data[it->Unknowns(lsidxnum)]) == 1)
            negPart[p1nr]-= p1x[extIdx[nr]];
        else
            posPart[p1nr]+= p1x[extIdx[nr]];
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
        const IdxT p1nr= it->Unknowns(p1idxnum);
        if (extIdx[nr]==NoIdx) continue;
        p1x[extIdx[nr]]+= (posPart[p1nr] - negPart[p1nr]);
    }
}

} // end of namespace DROPS

