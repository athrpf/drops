//**************************************************************************
// File:    integrTime.h                                                   *
// Content: classes that perform time-integration steps                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 19 2001                                           *
//**************************************************************************

#ifndef _INTEGRTIME_H_
#define _INTEGRTIME_H_

namespace DROPS
{

template <class StokesT, class SolverT>
class InstatStokesThetaSchemeCL
/*****************************************************************************
*   for solving the instationary Stokes equation of type StokesT with a 
*   1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
*                        theta=1/2 -> Crank-Nicholson (Trapezregel)
*
*   Inner stationary Stokes-type problems are solved with a SolverT-solver.
*   The matrices A, B, M and the rhs b, c of the Stokes class have to be set
*   properly! After construction, SetTimeStep has to be called once. Then
*   every DoStep performs one step in time. Changing time steps require 
*   further calls to SetTimeStep.
******************************************************************************/
{
  private:
    StokesT& _Stokes;
    SolverT& _solver;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    InstatStokesThetaSchemeCL( StokesT& Stokes, SolverT& solver, double theta= 0.5)
        : _Stokes( Stokes), _solver( solver), _b( &Stokes.b), _old_b( new VelVecDescCL),
          _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns), 
          _theta( theta)
    { 
        _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
        _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    }

    ~InstatStokesThetaSchemeCL()
    {
        if (_old_b == &_Stokes.b)
            delete _b;
        else
            delete _old_b; 
        delete _cplM; delete _old_cplM;
    }
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
        _dt= dt;
        ConvexComb( _mat, 1., _Stokes.M.Data, _theta*dt, _Stokes.A.Data);
    }
       
    void DoStep( VectorCL& v, VectorCL& p);
};


template <class StokesT, class SolverT>
class InstatStokesFracStepSchemeCL
/*****************************************************************************
*   for solving the instationary Stokes equation of type StokesT with a 
*   fractional-step-scheme.
*
*   Inner stationary Stokes-type problems are solved with a SolverT-solver.
*   The matrices A, B, M and the rhs b, c of the Stokes class have to be set
*   properly! After construction, SetTimeStep has to be called once. Then
*   every DoStep performs one step in time. Changing time steps require 
*   further calls to SetTimeStep.
******************************************************************************/
{
  private:
    StokesT& _Stokes;
    SolverT& _solver;
    
    VelVecDescCL *_cplA, *_old_cplA;  // couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MatrixCL      _L;                 // M + alpha*theta*dt*A
    
    const double _theta, _alpha;
    double       _dt;
    
  public:
    InstatStokesFracStepSchemeCL( StokesT& Stokes, SolverT& solver)
        : _Stokes( Stokes), _solver( solver), _cplA( new VelVecDescCL), _old_cplA( new VelVecDescCL), 
          _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns), 
          _theta( 1. - sqrt(2.)/2.), _alpha( (1. - 2.*theta)/(1. - theta))
    { 
        _cplA->SetIdx( _b->RowIdx); _old_cplA->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
        _Stokes.SetupInstatRhs( _old_cplA, &_Stokes.c, _old_cplM, _Stokes.t, &_Stokes.b, _Stokes.t);
    }

    ~InstatStokesFracStepSchemeCL()
    {
        delete _cplA; delete _old_cplA;
        delete _cplM; delete _old_cplM;
    }
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
        _dt= dt;
        ConvexComb( _L, 1., _Stokes.M.Data, _alpha*_theta*dt, _Stokes.A.Data);
    }
       
    void DoStep( VectorCL& v, VectorCL& p);
};


//=================================
//     template definitions
//=================================

template <class StokesT, class SolverT>
void InstatStokesThetaSchemeCL<StokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * v;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*v + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*_old_b->Data);

    p*= _dt;
    _solver.solve( _mat, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);
    p/= _dt;

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
}


template <class NavStokesT, class SolverT>
void InstatNavStokesThetaSchemeCL<NavStokesT,SolverT>::DoStep( VecDescCL& v, VectorCL& p)
{
//        _NS.SetupNonlinar( &_NS.N, v, &_NS.cplN, _t);
// sollte vorher schon geschehen sein: im 1. Zeitschritt durch User, danach durch _solver!

    _NS.t+= _dt;
    _NS.SetupInstatRhs( _b, &_NS.c, _cplM, _NS.t, _b, _NS.t);

    _rhs=  _NS.A.Data * v + _NS.N.Data * v;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _NS.M.Data*v + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*(_old_b->Data + _NS.cplN.Data);

    p*= _dt;
    _solver.solve( _L, _NS.B.Data, v, p, _rhs, _NS.c.Data, _dt*_theta);
    p/= _dt;

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
}

template <class StokesT, class SolverT>
void InstatStokesFracStepSchemeCL<StokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
        double macrostep= _theta*_dt;
        _Stokes.SetupInstatRhs( _cplA, &_Stokes.c, _cplM, _Stokes.t+macrostep, _Stokes.b, _Stokes.t);
        _rhs= _Stokes.A->Data * v;
        _rhs*= -(1-alpha)*macrostep;
        _rhs+= _Stokes.M->Data*v + macrostep*_Stokes.b.Data;
        _rhs+= _cplM->Data - _old_cplM->Data + macrostep * ( _alpha*_cplA->Data + (1.-_alpha)*(_old_cplA->Data + _NS.cplN.Data) );
        p*= _dt;
        _solver.Solve( _L, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);
        p/= _dt;
        std::swap( _cplA, _old_cplA);
        std::swap( _cplM, _old_cplM);

        macrostep= (1.-2.*_theta)*_dt;
        _Stokes.SetupInstatRhs( _cplA, &_Stokes.c, _cplm, _Stokes.t+(1.-_theta)*_dt, &_Stokes.b, _Stokes.t+(1.-_theta)*_dt);
        _rhs= _Stokes.A.Data * v;
        _rhs*= -_alpha*macrostep;
        _rhs+= _Stokes.M.Data*v + macrostep*_Stokes.b.Data;
        _rhs+= _cplM->Data - _old_cplM->Data + macrostep * ( _alpha*(_old_cplA->Data + _NS.cplN.Data) + (1.-_alpha)*_cplA->Data );
        p*= _dt;
        _solver.Solve( _L, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);
        p/= _dt;
        std::swap( _cplA, _old_cplA);
        std::swap( _cplM, _old_cplM);

        macrostep= _theta*_dt;
        _Stokes.SetupInstatRhs( _cplA, &_Stokes.c, _cplM, _Stokes.t+_dt, &_Stokes.b, _Stokes.t+(1.-_theta)*_dt);
        _rhs= _Stokes.A.Data * v;
        _rhs*= -(1.-_alpha)*macrostep;
        _rhs+= _Stokes.M.Data*v + macrostep*_Stokes.b.Data;
        _rhs+= _cplM->Data - _old_cplM->Data + macrostep * ( _alpha*_cplA->Data + (1.-_alpha)*(_old_cplA->Data + _NS.cplN.Data) );
        p*= _dt;
        _solver.Solve( _L, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);
        p/= _dt;
        std::swap( _cplA, _old_cplA);
        std::swap( _cplM, _old_cplM);

        _Stokes.t+= _dt;
}

}    // end of namespace DROPS

#endif
