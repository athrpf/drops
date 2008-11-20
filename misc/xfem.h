/// \file
/// \brief
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser, IGPM

#ifndef DROPS_XFEM_H_
#define DROPS_XFEM_H_

#include "levelset/levelset.h"
#include "misc/container.h"

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

    ExtIdxDescCL( IdxDescCL* idx, double omit_bound= 1./32. ) : omit_bound_( omit_bound ), Idx( idx ) {}
    ExtIdxDescCL() : omit_bound_( 1./32. ), Idx( 0 ){}
    void   SetIdx( IdxDescCL* idx ) { Idx = idx; }
    void   SetBound( double omit_bound ) { omit_bound_= omit_bound; }
    double GetBound() { return omit_bound_; }

    IdxT operator[]( const IdxT i ) const { return Xidx[i]; }
    IdxT GetNumUnknownsP1() const { return Xidx.size(); }

    void UpdateXNumbering( IdxDescCL*, const LevelsetP2CL&, bool NumberingChanged= false );
    void Old2New( VecDescCL* );
};

class MLExtIdxDescCL : public MLDataCL<ExtIdxDescCL>
{
  public:
    MLExtIdxDescCL( MLIdxDescCL* idx, double omit_bound= 1./32. )
    {
        for (MLIdxDescCL::iterator it = idx->begin(); it != idx->end(); ++it)
            this->push_back(ExtIdxDescCL(&(*it), omit_bound));
    }
    
    void UpdateXNumbering( IdxDescCL* idx,   const LevelsetP2CL& lset, bool NumberingChanged= false )
    {
        this->GetFinest().UpdateXNumbering( idx, lset, NumberingChanged);
    }
    
    void UpdateXNumbering( MLIdxDescCL* idx, const LevelsetP2CL& lset, bool NumberingChanged= false )
    {
        MLIdxDescCL::iterator itIdx = idx->begin();
        for (MLExtIdxDescCL::iterator it= this->begin(); it != this->end(); ++it)
        {
            it->UpdateXNumbering( &(*itIdx), lset, NumberingChanged);
            itIdx++;
        }            
    }
    
    void resize( MLIdxDescCL* idx, double omit_bound= 1./32.)
    {
        this->clear();
        for (MLIdxDescCL::iterator it = idx->begin(); it != idx->end(); ++it)
            this->push_back(ExtIdxDescCL(&(*it), omit_bound));
    }
    
    void Old2New( VecDescCL* v) { this->GetFinest().Old2New(v); }
    IdxT operator[]( const IdxT i ) const { return this->GetFinest().Xidx[i]; }
    IdxT GetNumUnknownsP1() const { return this->GetFinest().GetNumUnknownsP1(); }
};

/// merges two p1-VectorCL into a p1x-VectorCL
void P1toP1X ( const ExtIdxDescCL& xidx, VectorCL& p1x, const IdxDescCL& idx, const VectorCL& posPart, const VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg );

/// splits a p1x-VectorCL into two p1-VectorCL
void P1XtoP1 ( const ExtIdxDescCL& xidx, const VectorCL& p1x, const IdxDescCL& idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg );

} // end of namespace DROPS

#endif /* DROPS_XFEM_H_ */
