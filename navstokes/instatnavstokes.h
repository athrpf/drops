//**************************************************************************
// File:    instatnavstokes.h                                              *
// Content: classes that constitute the instationary navier-stokes-problem *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 22 2001                                           *
//**************************************************************************

#ifndef DROPS_INSTATNAVSTOKES_H
#define DROPS_INSTATNAVSTOKES_H

#include "stokes/instatstokes.h"


namespace DROPS
{

template <class MGB, class Coeff>
class InstatNavierStokesP2P1CL : public InstatStokesP2P1CL<MGB, Coeff>
{
  private:
    typedef InstatStokesP2P1CL<MGB, Coeff> BaseCL;
    
  public:
    typedef MGB                           MultiGridBuilderCL;
    typedef Coeff                         CoeffCL;
    typedef typename BaseCL::BndDataCL    BndDataCL;
    typedef typename BaseCL::VelVecDescCL VelVecDescCL;
  
    MatDescCL    N;
    VelVecDescCL cplN;
  
    InstatNavierStokesP2P1CL(const MultiGridBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : InstatStokesP2P1CL<MGB, Coeff>(mgb, coeff, bdata) {}  

    // Set up matrix for nonlinearity
    void SetupNonlinear(MatDescCL*, const VelVecDescCL*, VelVecDescCL*, double) const;
    void SetupNonlinear(MatDescCL* N, const VelVecDescCL* v, VelVecDescCL* cplN) const
        {SetupNonlinear(

    // Check system and computed solution
    void GetDiscError (vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr);
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, vector_fun_ptr, scalar_fun_ptr);
};


//double ResidualErrEstimator(const TetraCL& t, const VecDescCL& sol);
//double Estimator (const TetraCL& t, const VecDescCL& lsg);

} // end of namespace DROPS

//======================================
//  definition of template functions 
//======================================
#include "navstokes/navstokes.tpp"


#endif

