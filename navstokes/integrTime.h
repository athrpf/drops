//**************************************************************************
// File:    integrTime.h                                                   *
// Content: classes that perform time-integration steps                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 19 2001                                           *
//**************************************************************************

#ifndef DROPS_NS_INTEGRTIME_H
#define DROPS_NS_INTEGRTIME_H

#include "navstokes/instatnavstokes.h"

// TODO: FracStepScheme fuer instat. NavierStokes

namespace DROPS
{

template <class NavStokesT, class SolverT>
class InstatNavStokesThetaSchemeCL
/*****************************************************************************
*   for solving the instat. Navier-Stokes equation of type NavStokesT with a 
*   1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
*                        theta=1/2 -> Crank-Nicholson (Trapezregel)
*
*   Inner stat. Navier-Stokes-type problems are solved with a SolverT-solver.
*   The matrices A, B, M, N and the rhs b, c, cplN of the NavStokes class have
*   to be set properly! After construction, SetTimeStep has to be called once.
*   Then every DoStep performs one step in time. Changing time steps require 
*   further calls to SetTimeStep.
******************************************************************************/
{
  private:
//    typedef typename NavStokesT::VelVecDescCL VelVecDescCL;
    
    NavStokesT& _NS;
    SolverT&    _solver;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VelVecDescCL *_cplN;              // couplings with nonlinearity N
    VectorCL      _rhs;
    MatrixCL      _L;                 // M + theta*dt*A  = linear part 
    
    double _theta, _dt;
    
  public:
    InstatNavStokesThetaSchemeCL( NavStokesT& NS, SolverT& solver, const VelVecDescCL* v, double theta= 0.5)
        : _NS( NS), _solver( solver), _b( &NS.b), _old_b( new VelVecDescCL),
          _cplM( &NS.cplM), _old_cplM( new VelVecDescCL),
	  _cplN( &NS.cplN), _rhs( NS.b.RowIdx->NumUnknowns), _theta( theta)
    { 
        _old_b->SetIdx( _b->RowIdx);
        _old_cplM->SetIdx( _b->RowIdx);
        // Redundant for _NS.c but does not change its value
        _NS.SetupInstatRhs( _old_b, &_NS.c, _old_cplM, 0., _old_b, 0.);
    }

    ~InstatNavStokesThetaSchemeCL()
    {
        if (_old_b == &_NS.b)
            delete _b;
        else
            delete _old_b; 
        delete _old_cplM;
    }
    
    double GetTheta()    const { return _theta; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
        _dt= dt;
        _L.LinComb( 1., _NS.M.Data, _theta*_dt, _NS.A.Data);
    }
       
    void DoStep( VecDescCL& v, VectorCL& p);
};


//=================================
//     template definitions
//=================================

template <class NavStokesT, class SolverT>
void InstatNavStokesThetaSchemeCL<NavStokesT,SolverT>::DoStep( VecDescCL& v, VectorCL& p)
{
    // _NS._t contains the new time!
    _NS.SetupInstatRhs( _b, &_NS.c, _cplM, _NS._t, _b, _NS._t);
    const double alpha= _theta*_dt;
    const double beta= (_theta - 1.)*_dt;
    _rhs=  alpha*( _b->Data  + _cplN->Data)
         - beta*_old_b->Data
         + beta*( _NS.A.Data*_NS.v.Data + _NS.N.Data*_NS.v.Data )
         +_cplM->Data - _old_cplM->Data + _NS.M.Data*_NS.v.Data;
    p*= _dt;
    _solver.Solve( _L, _NS.B.Data, v, p, _rhs, *_cplN, _NS.c.Data, alpha);
    p/= _dt;
    _old_b->Data= _b->Data;
    _old_cplM->Data= _cplM->Data;
}


} // end of namespace DROPS

#endif
