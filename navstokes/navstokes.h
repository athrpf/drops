//**************************************************************************
// File:    navstokes.h                                                    *
// Content: classes that constitute the navier-stokes-problem              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 30 2001                                         *
//**************************************************************************

#ifndef DROPS_NAVSTOKES_H
#define DROPS_NAVSTOKES_H

#include "stokes/stokes.h"


namespace DROPS
{

template <class MGB, class Coeff>
class NavierStokesP2P1CL : public StokesP2P1CL<MGB, Coeff>
{
  private:
    typedef StokesP2P1CL<MGB, Coeff> _base;
    using                            _base::_MG;
    using                            _base::_BndData;
    using                            _base::b;
    using                            _base::c;
    using                            _base::A;
    using                            _base::B;

  public:
    typedef MGB                          MultiGridBuilderCL;
    typedef Coeff                        CoeffCL;
    typedef typename _base::BndDataCL    BndDataCL;
    typedef typename _base::VelVecDescCL VelVecDescCL;

    MatDescCL    N;
    VelVecDescCL cplN;
  
    NavierStokesP2P1CL(const MultiGridBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : StokesP2P1CL<MGB, Coeff>(mgb, coeff, bdata) {}  
    NavierStokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : StokesP2P1CL<MGB, Coeff>(mg, coeff, bdata) {}  

    // Set up matrix for nonlinearity
    void SetupNonlinear(MatDescCL*, const VelVecDescCL*, VelVecDescCL*) const;

    // Check system and computed solution
    void GetDiscError (vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr);
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, vector_fun_ptr, scalar_fun_ptr);
};


} // end of namespace DROPS

#include "navstokes/navstokes.tpp"

#endif
