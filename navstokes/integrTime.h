//**************************************************************************
// File:    integrTime.h                                                   *
// Content: classes that perform time-integration steps                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 19 2001                                           *
//**************************************************************************

#ifndef DROPS_NS_INTEGRTIME_H
#define DROPS_NS_INTEGRTIME_H

#include "misc/problem.h"

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
    NavStokesT& _NS;
    SolverT&    _solver;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MatrixCL      _L;                 // M + theta*dt*A  = linear part 
    
    double _theta, _dt;
    
  public:
    InstatNavStokesThetaSchemeCL( NavStokesT& NS, SolverT& solver, double theta= 0.5)
        : _NS( NS), _solver( solver), _b( &NS.b), _old_b( new VelVecDescCL),
          _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
          _rhs( NS.b.RowIdx->NumUnknowns), _theta( theta)
    { 
        _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
        _NS.SetupInstatRhs( _old_b, &_NS.c, _old_cplM, _NS.t, _old_b, _NS.t);
    }

    ~InstatNavStokesThetaSchemeCL()
    {
        if (_old_b == &_NS.b)
            delete _b;
        else
            delete _old_b; 
        delete _cplM; delete _old_cplM;
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
void InstatNavStokesFracStepSchemeCL<NavStokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
        double macrostep= _theta*_dt;
        _NS.SetupInstatRhs( _cplA, &_NS.c, _cplM, _NS.t+macrostep, _NS.b, _NS.t);
        _rhs= _NS.A->Data * v;
        _rhs*= -(1-alpha)*macrostep;
        _rhs+= _NS.M->Data*v + macrostep*_NS.b.Data;
        _rhs+= _cplM->Data - _old_cplM->Data + macrostep * ( _alpha*_cplA->Data + (1.-_alpha)*_old_cplA->Data );
        p*= _dt;
        _solver.Solve( _L, _NS.B.Data, v, p, _rhs, _NS.c.Data, _alpha*_theta*_dt);
        p/= _dt;
        std::swap( _cplA, _old_cplA);
        std::swap( _cplM, _old_cplM);

        macrostep= (1.-2.*_theta)*_dt;
        _NS.SetupInstatRhs( _cplA, &_NS.c, _cplm, _NS.t+(1.-_theta)*_dt, &_NS.b, _NS.t+(1.-_theta)*_dt);
        _rhs= _NS.A.Data * v;
        _rhs*= -_alpha*macrostep;
        _rhs+= _NS.M.Data*v + macrostep*_NS.b.Data;
        _rhs+= _cplM->Data - _old_cplM->Data + macrostep * ( _alpha*_old_cplA->Data + (1.-_alpha)*_cplA->Data );
        p*= _dt;
        _solver.Solve( _L, _NS.B.Data, v, p, _rhs, _NS.c.Data, _alpha*_theta*_dt);
        p/= _dt;
        std::swap( _cplA, _old_cplA);
        std::swap( _cplM, _old_cplM);

        macrostep= _theta*_dt;
        _NS.SetupInstatRhs( _cplA, &_NS.c, _cplM, _NS.t+_dt, &_NS.b, _NS.t+(1.-_theta)*_dt);
        _rhs= _NS.A.Data * v;
        _rhs*= -(1.-_alpha)*macrostep;
        _rhs+= _NS.M.Data*v + macrostep*_NS.b.Data;
        _rhs+= _cplM->Data - _old_cplM->Data + macrostep * ( _alpha*_cplA->Data + (1.-_alpha)*_old_cplA->Data );
        p*= _dt;
        _solver.Solve( _L, _NS.B.Data, v, p, _rhs, _NS.c.Data, _alpha*_theta*_dt);
        p/= _dt;
        std::swap( _cplA, _old_cplA);
        std::swap( _cplM, _old_cplM);

        _NS.t+= _dt;
}

}    // end of namespace DROPS

#endif
