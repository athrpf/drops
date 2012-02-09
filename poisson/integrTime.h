/// \file integrTime.h
/// \brief classes that perform time-integration steps
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_POI_INTEGRTIME_H
#define DROPS_POI_INTEGRTIME_H

#include "misc/problem.h"

/// \todo FracStepScheme fuer instat. Poisson


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
    PoissonT&   _Poisson;
    SolverT&    _solver;

    VecDescCL *_b, *_old_b;             // rhs
    VecDescCL *_cplA, *_old_cplA;       // couplings with poisson matrix A
    VecDescCL *_cplM, *_old_cplM;       // couplings with mass matrix M
    VecDescCL *_cplU;                   // couplings with convection matrix U
    VectorCL  _rhs;
    MLMatrixCL  _Lmat;                  // M + theta*dt*nu*A  = linear part

    double _theta, _dt;
    bool   _Convection;
    bool   _SUPG;
    instat_scalar_fun_ptr _presol;
    instat_scalar_fun_ptr _delta;
  public:
    InstatPoissonThetaSchemeCL( PoissonT& Poisson, SolverT& solver, double theta= 0.5, bool Convection= false, 
                                bool SUPG=false, instat_scalar_fun_ptr presol = NULL, instat_scalar_fun_ptr delta = NULL)
    : _Poisson( Poisson), _solver( solver),
      _b( &Poisson.b), _old_b( new VecDescCL),
      _cplA( new VecDescCL), _old_cplA( new VecDescCL),
      _cplM( new VecDescCL), _old_cplM( new VecDescCL),
      _cplU( new VecDescCL),
      _rhs( Poisson.b.RowIdx->NumUnknowns()), _theta( theta), _Convection( Convection), _SUPG(SUPG), 
      _presol(presol), _delta(delta)
    {
      _old_b->SetIdx( _b->RowIdx);
      _cplA->SetIdx( _b->RowIdx); _old_cplA->SetIdx( _b->RowIdx);
      _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
      _Poisson.SetupInstatRhs( *_old_cplA, *_old_cplM, _Poisson.x.t, *_old_b, _Poisson.x.t, _SUPG);
      
      if(_presol != NULL)
      {
        VecDescCL *tmp;
        tmp =new VecDescCL;
        tmp->SetIdx( _b->RowIdx);
        _Poisson.SetupGradSrc( *tmp, _presol, _delta, _Poisson.x.t);
        _old_b->Data += tmp->Data; 
      }
      if (Convection)
      {
        _cplU->SetIdx( _b->RowIdx);
        _Poisson.SetupConvection( _Poisson.U, *_cplU, _Poisson.x.t);
      }
      if(_SUPG)
      {
        std::cout << "----------------------------------------------------------------------------------\n"
                  <<"The SUPG stabilization has been added ...\n";  
      }
    }

    ~InstatPoissonThetaSchemeCL()
    {
      if (_old_b == &_Poisson.b)
        delete _b;
      else
        delete _old_b;
      delete _cplA; delete _old_cplA;
      delete _cplM; delete _old_cplM;
      delete _cplU;
    }

    double GetTheta()    const { return _theta; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
      _dt= dt;
      if (!_Convection)
        _Lmat.LinComb( 1., _Poisson.M.Data, _theta*_dt, _Poisson.A.Data);
    }

    void DoStep( VecDescCL& v);
};




//=================================
//     template definitions
//=================================


template <class PoissonT, class SolverT>
void InstatPoissonThetaSchemeCL<PoissonT,SolverT>::DoStep( VecDescCL& v)
{
  _Poisson.x.t+= _dt;
  
  if(_SUPG)
  _Poisson.SetupInstatSystem( _Poisson.A, _Poisson.M, _Poisson.x.t, _SUPG );
  
  _Poisson.SetupInstatRhs( *_cplA, *_cplM, _Poisson.x.t, *_b, _Poisson.x.t, _SUPG);
  if(_presol != NULL)
  {
      VecDescCL *tmp;
      tmp = new VecDescCL;
      tmp->SetIdx( _b->RowIdx);
      _Poisson.SetupGradSrc( *tmp, _presol, _delta, _Poisson.x.t);
      _b->Data += tmp->Data; 
  }

  _rhs = _Poisson.A.Data * v.Data;
  _rhs*= -_dt*(1.0-_theta);
  _rhs+= _Poisson.M.Data*v.Data
         + _dt*( _theta*_b->Data + (1.0-_theta)*(_old_b->Data))
         + (_dt*(1.0-_theta)) * _old_cplA->Data + (_dt*_theta) * _cplA->Data
         - _old_cplM->Data + _cplM->Data;

  if (_Convection)
  {
      _rhs+= (_dt*(1.0-_theta)) * (_cplU->Data - _Poisson.U.Data * v.Data );
      _Poisson.SetupConvection( _Poisson.U, *_cplU, _Poisson.x.t);
      _rhs+= (_dt*_theta) * _cplU->Data;
      MLMatrixCL AU;
      AU.LinComb( 1, _Poisson.A.Data, 1, _Poisson.U.Data);
      _Lmat.LinComb( 1, _Poisson.M.Data, _dt*_theta, AU);
  }
  _solver.Solve( _Lmat, v.Data, _rhs);

  std::swap( _b, _old_b);
  std::swap( _cplA, _old_cplA);
  std::swap( _cplM, _old_cplM);
}


}    // end of namespace DROPS

#endif
