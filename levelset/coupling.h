//**************************************************************************
// File:    coupling.h                                                     *
// Content: coupling of levelset and (Navier-)Stokes equations             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_COUPLING_H
#define DROPS_COUPLING_H

#include "levelset/levelset.h"

namespace DROPS
{

template <class StokesT, class SolverT>
class CouplStokesLevelsetCL
{
  private:
    StokesT&      _Stokes;
    SolverT&      _solver;
    LevelsetP2CL& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    CouplStokesLevelsetCL( StokesT& Stokes, LevelsetP2CL& ls, 
                           SolverT& solver, double theta= 0.5);
    ~CouplStokesLevelsetCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; _mat.LinComb( 1., _Stokes.M.Data, _theta*dt, _Stokes.A.Data); }
       
    void InitStep();
    // perform fixed point iteration
    void DoFPIter();
    void CommitStep();
    
    void DoStep( int maxFPiter= -1);  // combines the 3 former routines
};
    

template <class StokesT, class SolverT>
class CouplLevelsetStokesCL
{
  private:
    StokesT&      _Stokes;
    SolverT&      _solver;
    LevelsetP2CL& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    CouplLevelsetStokesCL( StokesT& Stokes, LevelsetP2CL& ls, 
                           SolverT& solver, double theta= 0.5);
    ~CouplLevelsetStokesCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; _mat.LinComb( 1., _Stokes.M.Data, _theta*dt, _Stokes.A.Data); }
       
    void InitStep();
    // perform fixed point iteration
    void DoFPIter();
    void CommitStep();
    
    void DoStep( int maxFPiter= -1);  // combines the 3 former routines
};


template <class StokesT, class SolverT>
class CouplLevelsetStokes2PhaseCL
{
  private:
    StokesT&      _Stokes;
    SolverT&      _solver;
    LevelsetP2CL& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    CouplLevelsetStokes2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, 
                           SolverT& solver, double theta= 0.5);
    ~CouplLevelsetStokes2PhaseCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; }
       
    void InitStep();
    // perform fixed point iteration
    void DoFPIter();
    void CommitStep();
    
    void DoStep( int maxFPiter= -1);  // combines the 3 former routines
    
    // update after grid has changed
    void Update();
};


template <class StokesT, class SolverT>
class CouplLevelsetNavStokes2PhaseCL
{
  private:
    StokesT&      _Stokes;
    SolverT&      _solver;
    LevelsetP2CL& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _AN,                // A + N
                  _mat;               // M + theta*dt*(A+N)
    
    double _theta, _dt;
    
  public:
    CouplLevelsetNavStokes2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, 
                           SolverT& solver, double theta= 0.5);
    ~CouplLevelsetNavStokes2PhaseCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; }
       
    void InitStep();
    // perform fixed point iteration
    void DoFPIter();
    void CommitStep();
    
    void DoStep( int maxFPiter= -1);  // combines the 3 former routines
    
    // update after grid has changed
    void Update();
};


} // end of namespace DROPS

#include "levelset/coupling.tpp"

#endif

