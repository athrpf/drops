//**************************************************************************
// File:    coupling.h                                                     *
// Content: coupling of levelset and (Navier-)Stokes equations             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_COUPLING_H
#define DROPS_COUPLING_H

#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "num/MGsolver.h"

namespace DROPS
{

template <class StokesT>
class CouplLevelsetBaseCL
{
  protected:
    StokesT&        _Stokes;
    LevelsetP2CL&   _LvlSet;

    VelVecDescCL *_b, *_old_b;       // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM; // couplings with mass matrix M
    VecDescCL    *_curv;             // curvature term, the old vector is not used in all derived classes.
    VectorCL      _rhs, _ls_rhs;
    MatrixCL*     _mat;              // 1./dt*M + theta*A

    double _theta, _dt;

  public:
    CouplLevelsetBaseCL( StokesT& Stokes, LevelsetP2CL& ls, double theta= 0.5);
    virtual ~CouplLevelsetBaseCL();

    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    virtual void SetTimeStep ( double dt) {_dt= dt;}

    virtual void DoStep( int maxFPiter= -1) = 0;

    // update after grid has changed
    virtual void Update() = 0;
};

template <class StokesT, class SolverT>
class CouplLevelsetStokesCL: public CouplLevelsetBaseCL<StokesT> 
{
  private:
    typedef CouplLevelsetBaseCL<StokesT> _base;
    using _base::_Stokes;
    using _base::_LvlSet;
    using _base::_b;       using _base::_old_b;
    using _base::_cplM;    using _base::_old_cplM;
    using _base::_curv;
    using _base::_rhs;
    using _base::_ls_rhs;
    using _base::_mat;
    using _base::_theta;
    using _base::_dt;

    SolverT&   _solver;
    VecDescCL* _old_curv;

  public:
    CouplLevelsetStokesCL( StokesT& Stokes, LevelsetP2CL& ls,
                           SolverT& solver, double theta= 0.5);
    ~CouplLevelsetStokesCL();

    void SetTimeStep( double dt) { 
        _base::SetTimeStep( dt); 
        _mat->LinComb( 1./dt, _Stokes.M.Data, _theta, _Stokes.A.Data); 
        _LvlSet.SetTimeStep( dt);
    }

    void InitStep();
    // perform fixed point iteration
    void DoFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1); // combines the 3 former routines

    void Update() {};
};

template <class StokesT, class SolverT>
class CouplLevelsetStokes2PhaseCL: public CouplLevelsetBaseCL<StokesT>
{
  private:
    typedef CouplLevelsetBaseCL<StokesT> _base;
    using _base::_Stokes;
    using _base::_LvlSet;
    using _base::_b;       using _base::_old_b;
    using _base::_cplM;    using _base::_old_cplM;
    using _base::_curv;
    using _base::_rhs;
    using _base::_ls_rhs;
    using _base::_mat;
    using _base::_theta;
    using _base::_dt;
    
    SolverT&   _solver;
    VecDescCL* _old_curv;
    VecDescCL  cplLB_;
    MatDescCL  LB_;
    bool       withProj_; // preceding projection step
    bool       _usematMG; // MG-hierachy for _mat
    MGDataCL*  _matMG;
    
  public:
    CouplLevelsetStokes2PhaseCL(StokesT& Stokes, LevelsetP2CL& ls,  
        SolverT& solver, double theta= 0.5, bool withProjection= false, bool usematMG= false, MGDataCL* matMG= 0);
    ~CouplLevelsetStokes2PhaseCL();

    void SetTimeStep( double dt) {
        _base::SetTimeStep(dt);
        _LvlSet.SetTimeStep( dt);
    }

    void InitStep();
    void DoProjectionStep();
    void DoFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

template <class StokesT, class SolverT>
class CouplLevelsetNavStokes2PhaseCL: public CouplLevelsetBaseCL<StokesT>
{
  private:
    typedef CouplLevelsetBaseCL<StokesT> _base;
    using _base::_Stokes;
    using _base::_LvlSet;
    using _base::_b;       using _base::_old_b;
    using _base::_cplM;    using _base::_old_cplM;
    using _base::_curv;
    using _base::_rhs;
    using _base::_ls_rhs;
    using _base::_mat; // 1./dt*M + theta*A + stab_*_theta*_dt*LB
    using _base::_theta;
    using _base::_dt;

    SolverT&     _solver;
    VelVecDescCL *_cplN, *_old_cplN;  // couplings with convection matrix N
    VecDescCL    *_old_curv;
    VecDescCL    cplLB_;
    MatDescCL    LB_;
    const double _nonlinear;
    bool         withProj_;
    const double stab_;

    void MaybeStabilize (VectorCL&);

  public:
    CouplLevelsetNavStokes2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                           SolverT& solver, double theta= 0.5, double nonlinear= 1, bool withProjection= false, double stab= 0.0);
    ~CouplLevelsetNavStokes2PhaseCL();

    void SetTimeStep( double dt, double theta= -1) {
        _base::SetTimeStep( dt);
        if (theta >=0 ) _theta= theta;
        _LvlSet.SetTimeStep( dt, theta);
    }

    void InitStep();
    void DoProjectionStep();
    void DoFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

template <class StokesT, class SolverT>
class CouplLevelsetNavStokesDefectCorr2PhaseCL: public CouplLevelsetBaseCL<StokesT>
{
  private:
    typedef CouplLevelsetBaseCL<StokesT> _base;
    using _base::_Stokes;
    using _base::_LvlSet;
    using _base::_b;       using _base::_old_b;
    using _base::_cplM;    using _base::_old_cplM;
    using _base::_curv;
    using _base::_rhs;
    using _base::_ls_rhs;
    using _base::_mat; // 1./dt*M + theta*A
    using _base::_theta;
    using _base::_dt;

    SolverT&     _solver;
    VelVecDescCL *_cplN, *_old_cplN;  // couplings with convection matrix N
    VecDescCL    *_old_curv;
    const double _nonlinear;
    const double stab_;

    void MaybeStabilize (VectorCL&);

  public:
    CouplLevelsetNavStokesDefectCorr2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                           SolverT& solver, double theta= 0.5, double nonlinear= 1, double stab= 0.0);
    ~CouplLevelsetNavStokesDefectCorr2PhaseCL();

    void SetTimeStep( double dt, double theta= -1) {
        _base::SetTimeStep( dt);
        if (theta >=0 ) _theta= theta;
        _LvlSet.SetTimeStep( dt, theta);
    }

    void InitStep();
    void DoFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};

template <class StokesT, class SolverT>
class CouplLsNsFracStep2PhaseCL : public CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>
{
  private:
    static const double facdt_[3];
    static const double theta_[3];

    typedef CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT> _base;

    double dt3_;
    int step_;

  public:
    CouplLsNsFracStep2PhaseCL( StokesT& Stokes, LevelsetP2CL& ls,
                               SolverT& solver, double nonlinear= 1, bool withProjection= false, double stab= 0.0, int step = -1)
        : _base( Stokes, ls, solver, 0.5, nonlinear, withProjection, stab), step_((step >= 0) ? step%3 : 0) {}

    double GetSubTimeStep() const { return facdt_[step_]*dt3_; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep( double dt, int step = -1) {
        dt3_= dt;
        if (step>=0) step_= step%3;
    }

    void DoSubStep( int maxFPiter= -1) {
        std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fractional Step Method: Substep " << step_ << '\n';
        _base::SetTimeStep( GetSubTimeStep(), GetSubTheta());
        _base::DoStep( maxFPiter);
        step_= (step_ + 1)%3;
    }

    void DoStep( int maxFPiter= -1) {
        DoSubStep( maxFPiter);
        DoSubStep( maxFPiter);
        DoSubStep( maxFPiter);
    }

    void Update() { _base::Update(); }
};

template <class NavStokesT, class SolverT>
const double CouplLsNsFracStep2PhaseCL<NavStokesT,SolverT>::facdt_[3]
//  = { 1./3, 1./3, 1./3 };
//  = { 1./3, 1./3, 1./3 };
  = { 1.0 - std::sqrt( 0.5), std::sqrt( 2.0) - 1.0, 1.0 - std::sqrt( 0.5) };

template <class NavStokesT, class SolverT>
const double CouplLsNsFracStep2PhaseCL<NavStokesT,SolverT>::theta_[3]
//  = { 1.0, 1.0, 1.0 };
//  = { 1./3, 5./6, 1./3 };
  = { 2.0 - std::sqrt( 2.0), std::sqrt( 2.0) - 1.0, 2.0 - std::sqrt( 2.0) };


template <class StokesT, class SolverT>
class CouplLsNsBaenschCL : public CouplLevelsetBaseCL<StokesT> 
{
  private:
    typedef CouplLevelsetBaseCL<StokesT> _base;
    using _base::_Stokes;
    using _base::_LvlSet;
    using _base::_b;    using _base::_old_b;
    using _base::_cplM;    using _base::_old_cplM;
    using _base::_curv;
    using _base::_rhs;
    using _base::_ls_rhs;
    using _base::_mat; // 1./dt*M + theta*A
    using _base::_theta;
    using _base::_dt;

    SolverT&                _solver;
    SSORPcCL                _pc;
    GMResSolverCL<SSORPcCL> _gm;

    MatrixCL      _AN;                // A + N
    VelVecDescCL *_cplA, *_old_cplA;  // couplings with stiff matrix A
    VelVecDescCL *_cplN, *_old_cplN;  // couplings with convection matrix N

    double       _alpha;
    const double _nonlinear;
    Uint         _iter_nonlinear;


  public:
    CouplLsNsBaenschCL( StokesT& Stokes, LevelsetP2CL& ls,
        SolverT& solver, int gm_iter, double gm_tol, double nonlinear= 1);
    ~CouplLsNsBaenschCL();

    void SetTimeStep( double dt) { _base::SetTimeStep(dt); }

    void InitStep( bool StokesStep= true);
    // perform fixed point iteration
    void DoStokesFPIter();
    void DoNonlinearFPIter();
    void CommitStep();

    void DoStep( int maxFPiter= -1);

    void Update();
};


} // end of namespace DROPS

#include "levelset/coupling.tpp"

#endif

