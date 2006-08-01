//**************************************************************************
// File:    fastmarch.h                                                    *
// Content: fast marching method for reparametrization                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_FASTMARCH_H
#define DROPS_FASTMARCH_H

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/spmat.h"
#include <set>

namespace DROPS
{

class FastMarchCL
{
  public:
    enum { Finished= 1, Close= 2, Far= 0};

    typedef SArrayCL<IdxT,4>        ReprTetraT;
    typedef std::vector<ReprTetraT> VertexNeighT;

  private:
    MultiGridCL&               MG_;
    VecDescCL&                 v_;
    const IdxT                 size_;
    std::set<IdxT>             Close_;
    VectorBaseCL<byte>         Typ_;
    VectorBaseCL<VertexNeighT> neigh_;
    VectorBaseCL<Point3DCL>    Coord_;
    VectorCL                   Old_;
    VectorBaseCL<IdxT>         map_;

    void   InitClose();
    IdxT   FindTrial() const;
    void   Update( const IdxT);
    double CompValueProj( IdxT, int num, const IdxT upd[3]) const;

    // variants for periodic boundaries
    void   InitClosePer();
    IdxT   FindTrialPer() const;
    void   UpdatePer( const IdxT);
    double CompValueProjPer( IdxT, int num, const IdxT upd[3]) const;
    inline IdxT Map( IdxT i) const { return i<size_ ? i: map_[i-size_]; }

  public:
    FastMarchCL( MultiGridCL& mg, VecDescCL& v)
      : MG_(mg), v_(v), size_(v.RowIdx->NumUnknowns), Typ_(Far, size_) {}

    void InitZero( bool ModifyZero= true);
    void RestoreSigns();
    void Reparam( bool ModifyZero= true);

    // variants for periodic boundaries
    void InitZeroPer( const BndDataCL<>&, bool ModifyZero= true);
    void ReparamPer( const BndDataCL<>&, bool ModifyZero= true);
};

} // end of namespace DROPS

#endif
