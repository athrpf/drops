//**************************************************************************
// File:    coupling.tpp                                                   *
// Content: coupling of levelset and (Navier-)Stokes equations             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

namespace DROPS
{

// ==============================================
//              CouplStokesLevelsetCL
// ==============================================

template <class StokesT, class SolverT>
CouplStokesLevelsetCL<StokesT,SolverT>::CouplStokesLevelsetCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _rhs( Stokes.b.RowIdx->NumUnknowns), _ls_rhs( ls.Phi.RowIdx->NumUnknowns),
    _theta( theta)
{ 
    _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
    _old_curv->SetIdx( _b->RowIdx); _curv->SetIdx( _b->RowIdx);
    _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    ls.AccumulateBndIntegral( *_old_curv);
}

template <class StokesT, class SolverT>
CouplStokesLevelsetCL<StokesT,SolverT>::~CouplStokesLevelsetCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*(_old_b->Data + _old_curv->Data));
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _dt*_theta*_curv->Data, _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving Stokes took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    // setup system for levelset eq. as velocity has changed
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);
    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::CommitStep()
{
    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_LvlSet.GetSolver().GetIter()==0) // no change of Phi -> no change of vel
            break;
    }
    CommitStep();
}


// ==============================================
//              CouplLevelsetStokesCL
// ==============================================

template <class StokesT, class SolverT>
CouplLevelsetStokesCL<StokesT,SolverT>::CouplLevelsetStokesCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _rhs( Stokes.b.RowIdx->NumUnknowns), _ls_rhs( ls.Phi.RowIdx->NumUnknowns),
    _theta( theta)
{ 
    _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
    _old_curv->SetIdx( _b->RowIdx); _curv->SetIdx( _b->RowIdx);
    _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    ls.AccumulateBndIntegral( *_old_curv);
}

template <class StokesT, class SolverT>
CouplLevelsetStokesCL<StokesT,SolverT>::~CouplLevelsetStokesCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*(_old_b->Data + _old_curv->Data));
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    // setup system for levelset eq.
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);
    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _dt*_theta*_curv->Data, _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving Stokes took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::CommitStep()
{
    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_solver.GetIter()==0) // no change of vel -> no change of Phi
            break;
    }
    CommitStep();
}


// ==============================================
//              CouplLevelsetStokes2PhaseCL
// ==============================================

template <class StokesT, class SolverT>
CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::CouplLevelsetStokes2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _theta( theta)
{ 
    Update(); 
}

template <class StokesT, class SolverT>
CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::~CouplLevelsetStokes2PhaseCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupRhs2( &_Stokes.c, _Stokes.t);

    _rhs=  _Stokes.A.Data * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data - _old_cplM->Data
         + (_dt*(1.-_theta))*(_old_b->Data + _old_curv->Data);
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    // setup system for levelset eq.
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);

    time.Stop();
    std::cerr << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();

    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _b, _cplM, _LvlSet, _Stokes.t);
    _mat.LinComb( 1., _Stokes.M.Data, _theta*_dt, _Stokes.A.Data);

    time.Stop();
    std::cerr << "Discretizing Stokes/Curv took "<<time.GetTime()<<" sec.\n";
    time.Reset();

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _cplM->Data + _dt*_theta*(_curv->Data + _b->Data), _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving Stokes took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::CommitStep()
{
    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_solver.GetIter()==0 && _LvlSet.GetSolver().GetResid()<_LvlSet.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cerr << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::Update()
{
    IdxDescCL* const vidx= &_Stokes.vel_idx;
    IdxDescCL* const pidx= &_Stokes.pr_idx;
    TimerCL time;
    time.Reset();
    time.Start();

    std::cerr << "Updating discretization...\n";
    // IndexDesc setzen
    _b->SetIdx( vidx);       _old_b->SetIdx( vidx); 
    _cplM->SetIdx( vidx);    _old_cplM->SetIdx( vidx);
    _curv->SetIdx( vidx);    _old_curv->SetIdx( vidx);
    _rhs.resize( vidx->NumUnknowns);
    _ls_rhs.resize( _LvlSet.idx.NumUnknowns);
    _Stokes.c.SetIdx( pidx);
    _Stokes.A.SetIdx( vidx, vidx);
    _Stokes.B.SetIdx( pidx, vidx);
    _Stokes.M.SetIdx( vidx, vidx);

    // Diskretisierung
    _LvlSet.AccumulateBndIntegral( *_old_curv);
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _old_b, _old_cplM, _LvlSet, _Stokes.t);
    _Stokes.SetupSystem2( &_Stokes.B, &_Stokes.c, _Stokes.t);
    
    time.Stop();
    std::cerr << "Discretizing took " << time.GetTime() << " sec.\n";
}


// ==============================================
//              CouplLevelsetNavStokes2PhaseCL
// ==============================================

template <class StokesT, class SolverT>
CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::CouplLevelsetNavStokes2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _rhs( Stokes.v.RowIdx->NumUnknowns), _ls_rhs( ls.Phi.RowIdx->NumUnknowns),
    _theta( theta)
{ 
    Update(); 
}

template <class StokesT, class SolverT>
CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::~CouplLevelsetNavStokes2PhaseCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupRhs2( &_Stokes.c, _Stokes.t);

    _rhs=  _AN * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data - _old_cplM->Data
         + (_dt*(1.-_theta))*(_old_b->Data + _old_curv->Data);
}

template <class StokesT, class SolverT>
void CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    // setup system for levelset eq.
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);

    time.Stop();
    std::cerr << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();

    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _b, _cplM, _LvlSet, _Stokes.t);
    _Stokes.SetupNonlinear( &_Stokes.N, &_Stokes.v, _b, _LvlSet, _Stokes.t);
    _AN.LinComb( 1., _Stokes.A.Data, 1., _Stokes.N.Data);
    _mat.LinComb( 1., _Stokes.M.Data, _theta*_dt, _AN);

    time.Stop();
    std::cerr << "Discretizing NavierStokes/Curv took "<<time.GetTime()<<" sec.\n";
    time.Reset();

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _cplM->Data + _dt*_theta*(_curv->Data + _b->Data), _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving NavierStokes took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::CommitStep()
{
    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_solver.GetIter()==0 && _LvlSet.GetSolver().GetResid()<_LvlSet.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cerr << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();
}

template <class StokesT, class SolverT>
void CouplLevelsetNavStokes2PhaseCL<StokesT,SolverT>::Update()
{
    IdxDescCL* const vidx= &_Stokes.vel_idx;
    IdxDescCL* const pidx= &_Stokes.pr_idx;
    TimerCL time;
    time.Reset();
    time.Start();

    std::cerr << "Updating discretization...\n";
    // IndexDesc setzen
    _b->SetIdx( vidx);       _old_b->SetIdx( vidx); 
    _cplM->SetIdx( vidx);    _old_cplM->SetIdx( vidx);
    _curv->SetIdx( vidx);    _old_curv->SetIdx( vidx);
    _rhs.resize( vidx->NumUnknowns);
    _ls_rhs.resize( _LvlSet.idx.NumUnknowns);
    _Stokes.c.SetIdx( pidx);
    _Stokes.A.SetIdx( vidx, vidx);
    _Stokes.B.SetIdx( pidx, vidx);
    _Stokes.M.SetIdx( vidx, vidx);
    _Stokes.N.SetIdx( vidx, vidx);

    // Diskretisierung
    _LvlSet.AccumulateBndIntegral( *_old_curv);
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _old_b, _old_cplM, _LvlSet, _Stokes.t);
    _Stokes.SetupSystem2( &_Stokes.B, &_Stokes.c, _Stokes.t);
    _Stokes.SetupNonlinear( &_Stokes.N, &_Stokes.v, _b, _LvlSet, _Stokes.t);
    _AN.LinComb( 1., _Stokes.A.Data, 1., _Stokes.N.Data);
    
    time.Stop();
    std::cerr << "Discretizing took " << time.GetTime() << " sec.\n";
}

} // end of namespace DROPS
