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
  
    MatDescCL    N;
    VelVecDescCL cplN;
    VelVecDescCL cplM;
    double _t; // Hack to allow the direct use of stationary NavStokes-Solvers by stripping the time
               // argument from SetupNonlinear.
	       // The base class already contains t for the same purpose!
  
    InstatNavierStokesP2P1CL(const MultiGridBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : InstatStokesP2P1CL<MGB, Coeff>(mgb, coeff, bdata) {}  

    // Set up matrix and rhs for nonlinearity
    void SetupNonlinear(MatDescCL*, const VelVecDescCL*, VelVecDescCL*, double, double) const;
    // Set up matrix for nonlinearity, use time _t
    void SetupNonlinear(MatDescCL* matN, const VelVecDescCL* velvec, VelVecDescCL* vecb) const
    { this->SetupNonlinear(matN, velvec, vecb, _t, _t); }
    // Set up only rhs for nonlinearity: use time t1 for the velocity in N,
    // t2 for the boundary-data in the velocity unknowns
//    void SetupNonlinearRhs(const VelVecDescCL*, VelVecDescCL*, double t1, double t2) const;
//    void SetupNonlinearRhs(const VelVecDescCL* velvec, VelVecDescCL* vecb) const
//    { this->SetupNonlinearRhs(velvec, vecb, _t, _t); }

    // Set time for use with stationary NavStokes-Solvers. This shall be the new time t_old+dt!!!!!!!!!!!!!!!!!!
    void SetTime (double tt) { _t= tt; }
    
    // Check system and computed solution
    void GetDiscError (vector_instat_fun_ptr LsgVel, vector_instat_fun_ptr DtLsgVel,
                       scalar_instat_fun_ptr LsgPr, double t);
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, vector_instat_fun_ptr, vector_instat_fun_ptr, scalar_instat_fun_ptr, double);
};


} // end of namespace DROPS

#include "navstokes/instatnavstokes.tpp"

#endif
