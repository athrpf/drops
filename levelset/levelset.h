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

class LevelsetP2CL
// P2-discretization and solution of the levelset equation for two phase
// flow problems.
{
  public:
    typedef P2EvalCL<double, const DummyBndDataCL, const VecDescCL> DiscSolCL;

    IdxDescCL           idx;
    VecDescCL           Phi;
    double              sigma;   // surface tension

  private:
    MultiGridCL&        _MG;
    double              _diff,   // amount of diffusion in reparametrization
                        _SD,     // streamline diffusion
                        _theta, _dt;  
    MatrixCL            _E, _H, _L;
    DummyBndDataCL      _dummyBnd;
    SSORPcCL            _pc;
    GMResSolverCL<SSORPcCL>  _gm;

    void SetupReparamSystem( MatrixCL&, const VectorCL&, VectorCL&);
    
  public:
    LevelsetP2CL( MultiGridCL& mg, double sig= 0, double theta= 0.5, double SD= 0., 
                  double diff= 0., Uint iter=1000, double tol=1e-7)
      : idx( 1, 1), sigma( sig), _MG( mg), _diff(diff), _SD( SD), _theta( theta), _dt( 0.),  
        _gm( _pc, 10, iter, tol)
    {}
    
    GMResSolverCL<SSORPcCL>& GetSolver() { return _gm; }
    
    const DummyBndDataCL& GetBndData() const { return _dummyBnd; }
    
    void CreateNumbering( Uint level, IdxDescCL*);
    void DeleteNumbering( IdxDescCL*);

    void Init( scalar_fun_ptr);

    void SetTimeStep( double dt) { _dt= dt; _L.LinComb( 1., _E, _theta*_dt, _H); }
    // call SetupSystem *before* calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&);
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


} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

