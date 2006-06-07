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
    MatrixCL      _mat;               // 1./dt*M + theta*A
    
    double _theta, _dt;
    
  public:
    CouplLevelsetStokesCL( StokesT& Stokes, LevelsetP2CL& ls, 
                           SolverT& solver, double theta= 0.5);
    ~CouplLevelsetStokesCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; _mat.LinComb( 1./dt, _Stokes.M.Data, _theta, _Stokes.A.Data); _LvlSet.SetTimeStep( dt); }
       
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
    MatrixCL*     _mat;               // 1./dt*M + theta*A

    bool          _usematMG;          // MG-hierachy for _mat
    MGDataCL*     _matMG;
        
    double _theta, _dt;
    
  public:
    CouplLevelsetStokes2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, 
        SolverT& solver, double theta= 0.5, bool usematMG= false, MGDataCL* matMG= 0);
    ~CouplLevelsetStokes2PhaseCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; _LvlSet.SetTimeStep( dt); }
       
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
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with stiff matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VelVecDescCL *_cplN, *_old_cplN;  // couplings with convection matrix N
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // 1./dt*M + theta*(A+N)
    
    double _theta, _dt;
    const double _nonlinear;
    
  public:
    CouplLevelsetNavStokes2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls, 
                           SolverT& solver, double theta= 0.5, double nonlinear= 1);
    ~CouplLevelsetNavStokes2PhaseCL();
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; _LvlSet.SetTimeStep( dt); }
       
    void InitStep();
    // perform fixed point iteration
    void DoFPIter();
    void CommitStep();
    
    void DoStep( int maxFPiter= -1);  // combines the 3 former routines
    
    // update after grid has changed
    void Update();
};


template <class StokesT, class SolverT>
class CouplLsNsBaenschCL
{
  private:
    StokesT&      _Stokes;
    SolverT&      _solver;
    SSORPcCL      _pc;
    GMResSolverCL<SSORPcCL> _gm;
    
    LevelsetP2CL& _LvlSet;
    
    VelVecDescCL *_b;                 // rhs
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VelVecDescCL *_cplA, *_old_cplA;  // couplings with stiff matrix A
    VelVecDescCL *_cplN, *_old_cplN;  // couplings with convection matrix N
    VecDescCL    *_curv;              // curvature term
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _AN,                // A + N
                  _mat;               // 1./dt*M + theta*(A+N)
    
    double        _theta, _alpha, _dt;
    const double  _nonlinear;
    Uint          _iter_nonlinear;
    
  public:
    CouplLsNsBaenschCL( StokesT& Stokes, LevelsetP2CL& ls, 
        SolverT& solver, int gm_iter, double gm_tol, double nonlinear= 1);
    ~CouplLsNsBaenschCL();
    
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt) { _dt= dt; }
       
    void InitStep( bool StokesStep= true);
    // perform fixed point iteration
    void DoStokesFPIter();
    void DoNonlinearFPIter();
    void CommitStep();
    
    void DoStep( int maxFPiter= -1);  // combines the 3 former routines
    
    // update after grid has changed
    void Update();
};


} // end of namespace DROPS

#include "levelset/coupling.tpp"

#endif

