//**************************************************************************
// File:    integrTime.h                                                   *
// Content: classes that perform time-integration steps                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - Nov, 26 2002                                           *
//**************************************************************************

#ifndef DROPS_POI_INTEGRTIME_H
#define DROPS_POI_INTEGRTIME_H

#include "misc/problem.h"

// TODO: FracStepScheme fuer instat. Poisson


namespace DROPS
{


/*****************************************************************************
*   for solving the instat. Poisson equation of type PoissonT with a 
*   1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
*                        theta=1/2 -> Crank-Nicholson (Trapezregel)
*
*   Inner stat. Poisson-type problems are solved with a SolverT-solver.
*   The matrices A, M and the rhs b of the Poisson class have
*   to be set properly! After construction, SetTimeStep has to be called once.
*   Then every DoStep performs one step in time. Changing time steps require 
*   further calls to SetTimeStep.
******************************************************************************/

template <class PoissonT, class SolverT>
class InstatPoissonThetaSchemeCL
{
  private:
    PoissonT& 	_Poisson;
    SolverT& 	_solver;
	 
    VecDescCL *_b, *_old_b;        	// rhs 
    VecDescCL *_cplA, *_old_cplA;	  // couplings with poisson matrix A
    VecDescCL *_cplM, *_old_cplM;  	// couplings with mass matrix M
    VectorCL 	_rhs;
    MatrixCL  _Lmat;                // M + theta*dt*A  = linear part 
    
    double _theta, _dt, _nu;
	
	
    
  public:
    InstatPoissonThetaSchemeCL( PoissonT& Poisson, SolverT& solver, double theta= 0.5)
    : _Poisson( Poisson), _solver( solver),
      _b( &Poisson.b), _old_b( new VecDescCL),
      _cplA( new VecDescCL), _old_cplA( new VecDescCL),
      _cplM( new VecDescCL), _old_cplM( new VecDescCL),
      _rhs( Poisson.b.RowIdx->NumUnknowns), _theta( theta)
    { 
      _old_b->SetIdx( _b->RowIdx);
      _cplA->SetIdx( _b->RowIdx); _old_cplA->SetIdx( _b->RowIdx);
      _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
      _Poisson.SetupInstatRhs( *_old_cplA, *_old_cplM, _Poisson.t, *_old_b, _Poisson.t);
    }

    ~InstatPoissonThetaSchemeCL()
    {
      if (_old_b == &_Poisson.b)
        delete _b;
      else
        delete _old_b;
      delete _cplA; delete _old_cplA; 
      delete _cplM; delete _old_cplM;
    }
    
    double GetTheta()    const { return _theta; }
    double GetNu()       const { return _nu; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt, double nu= 1.0)
    {
      _dt= dt;
      _nu= nu;
      ConvexComb( _Lmat, 1., _Poisson.M.Data, _theta*_dt*_nu, _Poisson.A.Data);
    }
       
    void DoStep( VecDescCL& v);
};




//=================================
//     template definitions
//=================================


template <class PoissonT, class SolverT>
void InstatPoissonThetaSchemeCL<PoissonT,SolverT>::DoStep( VecDescCL& v)
{
  _Poisson.t+= _dt;
  _Poisson.SetupInstatRhs( *_cplA, *_cplM, _Poisson.t, *_b, _Poisson.t);

  _rhs = _Poisson.A.Data * v.Data;
  _rhs*= -_dt*(1.0-_theta)*_nu;
  _rhs+= _Poisson.M.Data*v.Data 
         + _dt*( _theta*_b->Data + (1.0-_theta)*(_old_b->Data))
         + _dt*(1.0-_theta)*_nu*(_old_cplA->Data) + _dt*_theta*_nu*(_cplA->Data)
         - _old_cplM->Data + _cplM->Data;

  _solver.Solve( _Lmat, v.Data, _rhs);

  std::swap( _b, _old_b);
  std::swap( _cplA, _old_cplA);
  std::swap( _cplM, _old_cplM);
}


}    // end of namespace DROPS

#endif
