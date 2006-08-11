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
/// \brief Class for performing the fast marching algorithm for reparametrization of the levelset function
class FastMarchCL
{
  public:
    enum { Finished= 1, Close= 2, Far= 0};              ///< Types for verts (handled, next to handled, far away)

    typedef SArrayCL<IdxT,4>        ReprTetraT;         ///< Tetras stored by sysnum
    typedef std::vector<ReprTetraT> VertexNeighT;       ///< Tetras that are owned by a vert

  private:
    MultiGridCL&               MG_;                     // Multigrid on which the fast marching algorithm should be applied
    VecDescCL&                 v_;                      // Reference on the levelset VecDescCL
    const IdxT                 size_;                   // Number of unknowns of the levelset function
    std::set<IdxT>             Close_;                  // Verts that have finished neighbors
    VectorBaseCL<byte>         Typ_;                    // Each vert is Finished, Close or Far
    VectorBaseCL<VertexNeighT> neigh_;                  // Neighbor tetras for each vert
    VectorBaseCL<Point3DCL>    Coord_;                  // Coordinates for each vert
    VectorCL                   Old_;                    // Old values of the levelset function before reparametrization
    VectorBaseCL<IdxT>         map_;

    void   InitClose();                                 // Mark verts next to finished verts as close update them. Initialize also neigh_
    IdxT   FindTrial() const;                           // Find vert with minimal value in Close_
    void   Update( const IdxT);                         // Update value on a vert and put this vert into Close_
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

    /// \brief Finds the zero level of the levelset function and handle verts aroud it.
    void InitZero( bool ModifyZero= true);
    /// \brief Restore signes of the levelset function according to old signes.
    void RestoreSigns();
    /// \brief Reparametrization of the levelset function with the fast marching algorithm.
    /// This function also calls InitZero and RestorSigns.
    void Reparam( bool ModifyZero= true);

    /// \name variants for periodic boundaries
    //@{
    void InitZeroPer( const BndDataCL<>&, bool ModifyZero= true);
    void ReparamPer( const BndDataCL<>&, bool ModifyZero= true);
    //@}
};

} // end of namespace DROPS

#endif
