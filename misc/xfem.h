/// \file
/// \brief
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser, IGPM

#ifndef DROPS_XFEM_H_
#define DROPS_XFEM_H_

#include "levelset/levelset.h"

namespace DROPS
{


/// \brief Extended index for P1X_FE-elements.
///
/// Depending on the position of the zero-level of the levelset-function
/// additional dof are needed in the vertices. Internally, a std::vector
/// with a component for each vertex is kept. If the vertex has an
/// extended dof, the index is stored in this component; otherwise,
/// NoIdx is stored.
class ExtIdxDescCL
{
  private:
    double omit_bound_; ///< constant for stabilization of XFEM, controls omission of extended DoFs
  public:
    typedef std::vector<IdxT> ExtendedIdxT;

    IdxDescCL* Idx; ///< Pointer to the index-description.

    ExtendedIdxT   Xidx;
    ExtendedIdxT   Xidx_old;

    ExtIdxDescCL(IdxDescCL* idx, double omit_bound= 1./32.) : omit_bound_(omit_bound), Idx( idx) {}
    ExtIdxDescCL() : omit_bound_(1./32.), Idx(0){}
    void SetIdx(IdxDescCL* idx) {Idx = idx;}
    void SetBound(double omit_bound) {omit_bound_= omit_bound;}

    IdxT operator[](const IdxT i) const { return Xidx[i]; }
    IdxT GetNumUnknownsP1() const { return Xidx.size(); }

    void UpdateXNumbering(IdxDescCL*, const LevelsetP2CL&, bool NumberingChanged= false);
    void Old2New(VecDescCL*);
};

} // end of namespace DROPS

#endif /* DROPS_XFEM_H_ */
