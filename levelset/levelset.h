//**************************************************************************
// File:    levelset.h                                                     *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_LEVELSET_H
#define DROPS_LEVELSET_H

#include "stokes/stokes.h"
#include "num/spmat.h"

namespace DROPS
{

// This class is more or less a hack to describe, that no boundary conditions 
// are imposed. (Note, that P*EvalCL expects a certain BndDataCL.)
class DummyBndDataCL
{
  public:
    // default ctor, dtor, whatever
    typedef double bnd_type;

    static inline bool IsOnDirBnd (const VertexCL&) { return false; }
    static inline bool IsOnNeuBnd (const VertexCL&) { return false; }
    static inline bool IsOnDirBnd (const EdgeCL&)   { return false; }
    static inline bool IsOnNeuBnd (const EdgeCL&)   { return false; }
    
    static inline bnd_type GetDirBndValue (const VertexCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
    static inline bnd_type GetDirBndValue (const EdgeCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on edge."); }
};

template<class StokesProblemT>
class LevelsetP2CL
// P2-discretization and solution of the levelset equation for two phase
// flow problems.
{
  public:
    typedef P2EvalCL<double, const DummyBndDataCL, const VecDescCL> DiscSolCL;
    typedef typename StokesProblemT::DiscVelSolCL                   DiscVelSolCL;

    IdxDescCL           idx;
    VecDescCL           Phi;
    double              sigma;   // surface tension

  private:
    MultiGridCL&        _MG;
    double              _SD,     // streamline diffusion
                        _theta, _dt;  
    MatrixCL            _E, _H, _L;
    DummyBndDataCL      _dummyBnd;
    DummyPcCL           _pc;
    GMResSolverCL<DummyPcCL>  _gm;

    void SetupReparamSystem( MatrixCL&, const VectorCL&, VectorCL&);
    
  public:
    LevelsetP2CL( MultiGridCL& mg, double sig= 0, double theta= 0.5, double SD= 0., Uint iter=1000, double tol=1e-7)
      : idx( 1, 1), sigma( sig), _MG( mg), _SD( SD), _theta( theta), _dt( 0.),  
        _gm( _pc, 10, iter, tol)
    {}
    
    GMResSolverCL<DummyPcCL>& GetSolver() { return _gm; }
    
    const DummyBndDataCL& GetBndData() const { return _dummyBnd; }
    
    void CreateNumbering( Uint level, IdxDescCL*);
    void DeleteNumbering( IdxDescCL*);

    void Init( scalar_fun_ptr);

    void SetTimeStep( double dt) { _dt= dt; _L.LinComb( 1., _E, _theta*_dt, _H); }
    // call SetupSystem *before* calling SetTimeStep!
    void SetupSystem( const DiscVelSolCL&);
    void DoStep();
    void Reparam( Uint steps, double dt);
//    void Reparam2();
    
    bool   Intersects( const TetraCL&) const;
    double GetMass() const;
    void   AccumulateBndIntegral( VecDescCL& f) const;
    
    DiscSolCL GetSolution() const
        { return DiscSolCL( &Phi, &_dummyBnd, &_MG); }
        
    // the following member functions are added to enable an easier implementation
    // of the coupling navstokes-levelset. They should not be called by a common user.
    void ComputeRhs( VectorCL&) const;
    void DoStep    ( const VectorCL&);
};


template <class StokesT, class SolverT>
class CouplStokesLevelsetCL
{
  private:
    StokesT&               _Stokes;
    SolverT&               _solver;
    LevelsetP2CL<StokesT>& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    CouplStokesLevelsetCL( StokesT& Stokes, LevelsetP2CL<StokesT>& ls, 
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
    StokesT&               _Stokes;
    SolverT&               _solver;
    LevelsetP2CL<StokesT>& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    CouplLevelsetStokesCL( StokesT& Stokes, LevelsetP2CL<StokesT>& ls, 
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
    StokesT&               _Stokes;
    SolverT&               _solver;
    LevelsetP2CL<StokesT>& _LvlSet;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VecDescCL    *_curv, *_old_curv;  // curvature terms
    VectorCL      _rhs, _ls_rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    CouplLevelsetStokes2PhaseCL( StokesT& Stokes, LevelsetP2CL<StokesT>& ls, 
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


} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

