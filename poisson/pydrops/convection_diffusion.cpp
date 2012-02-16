/// \file convection_diffusion.cpp
/// \brief Solver for Poisson problem with P1 functions for python interface
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

// include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "geom/geomselect.h"
#include "misc/bndmap.h"                //include function container

// include numeric computing!
#include "num/fe.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include "poisson/integrTime.h"
#include "num/poissonsolverfactory.h"

// include problem class
#include "misc/params.h"
//#include "poisson/poissonCoeff.h"      // Coefficient-Function-Container poissonCoeffCL
#include "poisson/poisson.h"      // setting up the Poisson problem
#include "num/bndData.h"

// include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

// include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "num/parsolver.h"              // various parallel solvers
#include "num/parprecond.h"             // various parallel preconditioners
#endif

#include "pyconnect.h"
#include "pdefunction.h"
#include <math.h>

using namespace std;

const char line[] ="--------------------------------------------------------------------\n";

namespace DROPS
{
  template<class CoeffCL, class SolverT>
  void SolveStatProblem( PoissonP1CL<CoeffCL>& Poisson, SolverT& solver, ParamCL& P, PythonConnectCL::Ptr PyC)
  {
    // time measurements
#ifndef _PAR
    TimerCL timer;
    const bool doErrorEstimate= P.get<int>("Err.DoErrorEstimate");
#else
    const bool doErrorEstimate= false;
    if (P.get<int>("Err.DoErrorEstimate"))
      *(PyC->outfile) << "Skipping Error-Estimation ..." << std::endl;
    ParTimerCL timer;
#endif

    if ( !doErrorEstimate) {
      std::string Gradstr ("IA2Gradient");
      std::string Sensstr ("IA2Sensitivity");
      std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
      bool GradProb = (Gradstr.compare(IAProbstr) == 0);
      bool SensProb = (Sensstr.compare(IAProbstr) == 0);
      Poisson.SetupSystem( Poisson.A, Poisson.b, P.get<int>("Stabilization.SUPG"), GradProb);
      if(P.get<int>("PoissonCoeff.Convection"))
        {
	  Poisson.vU.SetIdx( &Poisson.idx);
	  Poisson.SetupConvection(Poisson.U, Poisson.vU, 0.0);
	  Poisson.A.Data.LinComb(1., Poisson.A.Data, 1., Poisson.U.Data);
	  Poisson.b.Data+=Poisson.vU.Data;
        }
      if(GradProb)
        {
	  VecDescCL b1;
	  b1.SetIdx(&Poisson.idx);
	  Poisson.SetupL2ProjGrad(b1,*PyC->GetPresol, *PyC->GetDelPsi, NULL);
	  Poisson.b.Data+=b1.Data;
	  *(PyC->outfile) << line << "We are solving Gradient problem in IA2 ...\n"<<line;
        }
      if(SensProb)
        {
	  VecDescCL b1;
	  b1.SetIdx(&Poisson.idx);
	  Poisson.SetupGradSrc(b1,*PyC->GetPresol, *PyC->GetDelPsi);
	  Poisson.b.Data+=b1.Data;
	  *(PyC->outfile) << line << "We are solving sensetivity problem in IA2 ...\n"<<line;
        }
      if(P.get<int>("Stabilization.SUPG"))
        {
	  //CoeffCL::Show_Pec();
	  *(PyC->outfile) << line << "The SUPG stabilization has been added ...\n"<<line;
        }
      timer.Reset();
      solver.Solve( Poisson.A.Data, Poisson.x.Data, Poisson.b.Data);
      timer.Stop();
#ifndef _PAR
      double realresid = norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data));
#else
      double realresid = Poisson.idx.GetEx().Norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data), false);
#endif
      //Setup solution for python interface
      PyC->SetSol3D(Poisson.GetSolution());
      *(PyC->outfile) << " o Solved system with:\n"
		<< "   - time          " << timer.GetTime()   << " s\n"
		<< "   - iterations    " << solver.GetIter()  << '\n'
		<< "   - residuum      " << solver.GetResid() << '\n'
		<< "   - real residuum " << realresid         << std::endl;
      if(solver.GetResid() > P.get<double>("Poisson.Tol"))
        {
	  throw (PyDropsErr(PyC, P, 0, "The residual is bigger than tolerance"));
        }
      if(solver.GetResid()!=solver.GetResid())
        {
	  std::stringstream ss;
	  ss << "The linear solver in DROPS did not converge.\nIterations: " << solver.GetIter() << "\nMaxIter: " << solver.GetMaxIter() << std::endl;
	  throw (PyDropsErr(PyC, P, 0, ss.str()));
        }
      if (P.get<int>("Poisson.SolutionIsKnown")) {
	*(PyC->outfile) << line << "Check result against known solution ...\n";
	exit(1);
	//Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
      }
    }
    /*else{
      MultiGridCL& MG= Poisson.GetMG();
      const typename PoissonP1CL<CoeffCL>::BndDataCL& BndData= Poisson.GetBndData();

      MLIdxDescCL  loc_idx;
      VecDescCL  loc_x;
      MLIdxDescCL* new_idx= &Poisson.idx;
      MLIdxDescCL* old_idx= &loc_idx;
      VecDescCL* new_x= &Poisson.x;
      VecDescCL* old_x= &loc_x;

      //DoerflerMarkCL<typename PoissonP1CL<CoeffCL>::est_fun, typename PoissonP1CL<CoeffCL>::base_>
      //Estimator( P.get<double>("Err.RelReduction"), P.get<double>("Err.MinRatio"), P.get<double>("Err.Threshold"), P.get<double>("Err.Meas"), P.get<int>("Err.DoMark"),
      //&PoissonP1CL<CoeffCL>::ResidualErrEstimator, *static_cast<typename PoissonP1CL<CoeffCL>::base_*>(&Poisson) );

      int step= 0;
      bool new_marks;

      new_idx->SetFE( P1_FE);
      old_idx->SetFE( P1_FE);

      do{
	timer.Reset();
	//*(PyC->outfile) << DROPS::SanityMGOutCL(MG) << std::endl;
	MG.Refine();

	Poisson.CreateNumbering( MG.GetLastLevel(), new_idx);    // create numbering for this idx
	*(PyC->outfile) << "new triangLevel: " << Poisson.idx.TriangLevel() << std::endl;
	Poisson.b.SetIdx( new_idx);                              // tell b about numbering
	new_x->SetIdx( new_idx);                    			 // second vector with the same idx

	*(PyC->outfile) << line << "Problem size\no number of unknowns             " << new_x->Data.size() << std::endl;

	MG.SizeInfo(*(PyC->outfile));
	if ( step == 0)
	  //Estimator.Init( typename PoissonP1CL<CoeffCL>::DiscSolCL( new_x, &BndData, &MG));

	if ( old_x->RowIdx)
	  {
	    P1EvalCL<double, const PoissonBndDataCL, const VecDescCL>  oldx( old_x, &BndData, &MG);
	    P1EvalCL<double, const PoissonBndDataCL, VecDescCL>        newx( new_x, &BndData, &MG);
	    Interpolate( newx, oldx);
	    //            CheckSolution(*new_x, &::Lsg);
	    old_x->Reset();
	  }

	Poisson.A.SetIdx( new_idx, new_idx);             // tell A about numbering
	Poisson.SetupSystem( Poisson.A, Poisson.b);
	timer.Stop();
	timer.Reset();
	solver.Solve( Poisson.A.Data, new_x->Data, Poisson.b.Data);
	timer.Stop();
	//Setup solution for python interface
	PyC->SetSol3D(Poisson.GetSolution());
	double realresid = norm( VectorCL(Poisson.A.Data*new_x->Data-Poisson.b.Data));
	*(PyC->outfile) << " o Solved system with:\n"
		  << "   - time          " << timer.GetTime()   << " s\n"
		  << "   - iterations    " << solver.GetIter()  << '\n'
		  << "   - residuum      " << solver.GetResid() << '\n'
		  << "   - real residuum " << realresid         << std::endl;
	if(solver.GetResid() > P.get<double>("Poisson.Tol"))
	  {
	    throw (PyDropsErr(PyC, P, 0, "The residual is bigger than tolerence"));
	  }
	if(solver.GetResid()!=solver.GetResid())
	  {
	    std::stringstream ss;
	    ss << "The linear solver in DROPS did not converge.\nIterations: " << solver.GetIter() << "\nMaxIter: " << solver.GetMaxIter() << std::endl;
	    throw (PyDropsErr(PyC, P, 0, ss.str()));
	  }
	Poisson.A.Reset();
	Poisson.b.Reset();
	if (P.get<int>("Poisson.SolutionIsKnown")) {
	  *(PyC->outfile) << line << "Check result against known solution ...\n";
	  exit(1);
	  //Poisson.CheckSolution( *new_x, CoeffCL::Solution);
	}
	//new_marks = Estimator.Estimate( typename PoissonP1CL<CoeffCL>::const_DiscSolCL( new_x, &BndData, &MG) );

	std::swap( old_x, new_x);
	std::swap( old_idx, new_idx);
      } while ( new_marks && step++ < P.get<int>("Err.NumRef"));
      // I want the solution to be in Poisson.x
      if ( old_x == &loc_x)
        {
	  Poisson.idx.swap( loc_idx);
	  Poisson.x.SetIdx( &Poisson.idx);

	  Poisson.x.Data.resize( loc_x.Data.size());
	  Poisson.x.Data = loc_x.Data;
        }
	}*/
  }

  /// \brief Strategy to solve the Poisson problem on a given triangulation
  template<class CoeffCL>
  void Strategy( PoissonP1CL<CoeffCL>& Poisson, ParamCL& P, PythonConnectCL::Ptr PyC)
  {
    // time measurements
#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif

    // the triangulation
    MultiGridCL& mg= Poisson.GetMG();

    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    *(PyC->outfile) << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    Poisson.idx.SetFE( P1_FE);                                  // set quadratic finite elements
    //see class for explanation: template didnt work
    if ( PoissonSolverFactoryHelperCL().MGUsed(P))
      Poisson.SetNumLvl ( mg.GetNumLevel());
    Poisson.CreateNumbering( mg.GetLastLevel(), &Poisson.idx);  // number vertices and edges
    Poisson.b.SetIdx( &Poisson.idx);                            // tell b about numbering
    Poisson.x.SetIdx( &Poisson.idx);                            // tell x about numbering
    Poisson.A.SetIdx( &Poisson.idx, &Poisson.idx);              // tell A about numbering
    Poisson.M.SetIdx( &Poisson.idx, &Poisson.idx);              // tell M about numbering
    Poisson.U.SetIdx( &Poisson.idx, &Poisson.idx);

    timer.Stop();
    *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;


    // display problem size
    // -------------------------------------------------------------------------
    *(PyC->outfile) << line << "Problem size\n";
#ifdef _PAR
    std::vector<size_t> UnkOnProc= ProcCL::Gather( Poisson.x.Data.size(), 0);
    const IdxT numUnk   = Poisson.idx.GetGlobalNumUnknowns( mg),
      numAccUnk= std::accumulate(UnkOnProc.begin(), UnkOnProc.end(), 0);
#else
    std::vector<size_t> UnkOnProc( 1);
    UnkOnProc[0]  = Poisson.x.Data.size();
    IdxT numUnk   = Poisson.x.Data.size(),
      numAccUnk= Poisson.x.Data.size();
#endif
    *(PyC->outfile) << " o number of unknowns             " << numUnk    << '\n'
              << " o number of accumulated unknowns " << numAccUnk << '\n'
              << " o number of unknowns on proc\n";
    for (size_t i=0; i<UnkOnProc.size(); ++i)
      *(PyC->outfile) << " - Proc " << i << ": " << UnkOnProc[i]<< '\n';

    // discretize (setup linear equation system)
    // -------------------------------------------------------------------------
    *(PyC->outfile) << line << "Discretize (setup linear equation system) ...\n";

    timer.Reset();
    if (P.get<int>("Time.NumSteps") != 0)
      Poisson.SetupInstatSystem( Poisson.A, Poisson.M, Poisson.x.t, P.get<int>("Stabilization.SUPG"));
    timer.Stop();
    *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;


    // solve the linear equation system
    // -------------------------------------------------------------------------
    *(PyC->outfile) << line << "Solve the linear equation system ...\n";

    // type of preconditioner and solver
    PoissonSolverFactoryCL<> factory( P, Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    if ( factory.GetProlongation() != 0)
      SetupP1ProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);

    // Solve the linear equation system
    if(P.get<int>("Time.NumSteps") !=0)
      Poisson.Init( Poisson.x, *PyC->GetInitial, 0.0);
    else
      SolveStatProblem( Poisson, *solver, P, PyC);

    if(P.get<int>("Time.NumSteps")!=0)
      {


      std::string IA12Sensstr ("IA12Sensitivity");
      std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
      bool IA12SensProb = (IA12Sensstr.compare(IAProbstr) == 0);
      typedef boost::shared_ptr<InstatPoissonThetaSchemeCL<PoissonP1CL<CoeffCL>, PoissonSolverBaseCL> > ThetaPtr;
      ThetaPtr ThetaScheme;
      if(IA12SensProb)
      {
        ThetaScheme = ThetaPtr(new InstatPoissonThetaSchemeCL<PoissonP1CL<CoeffCL>, PoissonSolverBaseCL>( Poisson, *solver, P.get<double>("Time.Theta") , P.get<int>("PoissonCoeff.Convection"), P.get<int>("Stabilization.SUPG"), *PyC->GetPresol, *PyC->GetDelPsi));
      }
      else
      {
        ThetaScheme = ThetaPtr(new InstatPoissonThetaSchemeCL<PoissonP1CL<CoeffCL>, PoissonSolverBaseCL>( Poisson, *solver, P.get<double>("Time.Theta") , P.get<int>("PoissonCoeff.Convection"), P.get<int>("Stabilization.SUPG")));
      }
      ThetaScheme->SetTimeStep(P.get<double>("Time.StepSize") );
      for ( int step = 1; step <= P.get<int>("Time.NumSteps") ; ++step) {
	  timer.Reset();
	  *(PyC->outfile) << line << "Step: " << step << std::endl;
	  ThetaScheme->DoStep( Poisson.x);
	  //Setup solutions for python interface
	  PyC->SetSol3D(Poisson.GetSolution(), Poisson.x.t);

	  timer.Stop();
	  *(PyC->outfile) << " o Solved system with:\n"
		    << "   - time          " << timer.GetTime()    << " s\n"
		    << "   - iterations    " << solver->GetIter()  << '\n'
		    << "   - residuum      " << solver->GetResid() << '\n';
	  if (isnan(solver->GetResid())) {
	    std::stringstream ss;
	    ss << "The residual of the solver is NAN!";
	    throw(PyDropsErr(PyC, P, step-1, ss.str()));
	  }
	  if (isinf(solver->GetResid())) {
	    std::stringstream ss;
	    ss << "The residual of the solver is infinite!";
	    throw(PyDropsErr(PyC, P, step-1, ss.str()));
	  }
	  if(solver->GetResid() > P.get<double>("Poisson.Tol")) {
	    std::stringstream ss;
	    ss << "After " << solver->GetIter() << " iterations, the residual " << solver->GetResid() << " is still greater than the tolerance " << P.get<double>("Poisson.Tol");
	    throw (PyDropsErr(PyC, P, step-1, "The residual is greater than tolerence"));
	  }
	  if (P.get("Poisson.SolutionIsKnown", 0)) {
	    *(PyC->outfile) << line << "Check result against known solution ...\n";
	    exit(1);
	    //Poisson.CheckSolution( Poisson.x, CoeffCL::Solution, Poisson.x.t);
	  }
        }
      }

    delete solver;
  }

} // end of namespace DROPS

void SetMissingParameters(DROPS::ParamCL& P){
  P.put_if_unset<int>("Stabilization.SUPG",0);
  P.put_if_unset<std::string>("PoissonCoeff.IAProb", "IA1Direct");
}
//mainly used to solve a direct, sensetivity or adjoint problem in IA1
void convection_diffusion(std::ofstream& outfile, DROPS::ParamCL& P, PdeFunction::ConstPtr C0, PdeFunction::ConstPtr b_in, PdeFunction::ConstPtr b_interface, PdeFunction::ConstPtr source, PdeFunction::ConstPtr Dw, double* C_sol)
{
  PythonConnectCL::Ptr PyC(new PythonConnectCL());
  PyC->Init(&outfile, P, C0, b_in, source, Dw, b_interface, C_sol);
  SetMissingParameters(P);
#ifdef _PAR
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
#endif
  //try
  //{
  // time measurements
#ifndef _PAR
  DROPS::TimerCL timer;
#else
  DROPS::ParTimerCL timer;
#endif

  // set up data structure to represent a poisson problem
  // ---------------------------------------------------------------------
  *(PyC->outfile) << line << "Set up data structure to represent a Poisson problem ...\n";
  timer.Reset();

  //create geometry
  DROPS::MultiGridCL* mg= 0;

  const bool isneumann[6]=
    { false, true,              // inlet, outlet
      true,  false,             // wall, interface
      true,  true };            // in Z direction

  //DROPS::PoissonCoeffCL<DROPS::ParamCL> PoissonCoeff(P, *PyC->GetDiffusion, *PyC->GetSource, *PyC->GetInitial);

  const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
    { *PyC->GetInflow, &Zero, &Zero, *PyC->GetInterfaceValue, &Zero, &Zero};

  DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);

  //only for measuring cell, not used here
  double r = 1;
  std::string serfile = "none";

  DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);

  // Setup the problem
  //DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), bdata);
  std::string adstr ("IA1Adjoint");
  std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
  //DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, PoissonCoeff, bdata, adstr.compare(IAProbstr) == 0);
  DROPS::PoissonP1CL< PythonConnectCL > prob( *mg, *PyC, bdata, adstr.compare(IAProbstr) == 0);

#ifdef _PAR
  // Set parallel data structures
  DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
  pmg.AttachTo( *mg);                                  // handling of parallel multigrid
  DROPS::LoadBalHandlerCL lb( *mg, DROPS::metis);     // loadbalancing
  lb.DoInitDistribution( DROPS::ProcCL::Master());    // distribute initial grid
  lb.SetStrategy( DROPS::Recursive);                  // best distribution of data
#endif
  timer.Stop();
  *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;

  // Refine the grid
  // ---------------------------------------------------------------------
  *(PyC->outfile) << "Refine the grid " << P.get<int>("DomainCond.RefineSteps") << " times regulary ...\n";
  timer.Reset();
  // Create new tetrahedra
  for ( int ref=1; ref <= P.get<int>("DomainCond.RefineSteps"); ++ref){
    *(PyC->outfile) << " refine (" << ref << ")\n";
    DROPS::MarkAll( *mg);
    mg->Refine();
  }
  // do loadbalancing
#ifdef _PAR
  lb.DoMigration();
#endif

  timer.Stop();
  *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;
  mg->SizeInfo(*(PyC->outfile));

  //prepare MC
  PyC->SetMG(mg);
  PyC->ClearMaps();
  PyC->setFaceMap();
  PyC->setTetraMap();

  // Solve the problem
  DROPS::Strategy( prob, P, PyC);

  delete mg;
}

//Used to solve a direct, sensetivity, adjoint and gradient problem in IA2; source for direct and adjoint; presol for sensetivity and gradient; DelPsi for sensetivity and gradient
void CoefEstimation(std::ofstream& outfile, DROPS::ParamCL& P, PdeFunction::ConstPtr b_in, PdeFunction::ConstPtr b_interface, PdeFunction::ConstPtr source, PdeFunction::ConstPtr presol, PdeFunction::ConstPtr DelPsi, PdeFunction::ConstPtr Dw, double* C_sol)
{
  PythonConnectCL::Ptr PyC(new PythonConnectCL());
  PyC->Init(&outfile, P, b_in, b_interface, source, presol, DelPsi,  Dw, C_sol);
  SetMissingParameters(P);
#ifdef _PAR
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
#endif
  //try
  //{
  // time measurements
#ifndef _PAR
  DROPS::TimerCL timer;
#else
  DROPS::ParTimerCL timer;
#endif

  // set up data structure to represent a poisson problem
  // ---------------------------------------------------------------------
  *(PyC->outfile) << line << "Set up data structure to represent a Poisson problem ...\n";
  timer.Reset();

  //create geometry
  DROPS::MultiGridCL* mg= 0;
  //DROPS::PoissonBndDataCL* bdata = 0;

  const bool isneumann[6]=
    { false, true,              // inlet, outlet
      true,  false,             // wall, interface
      true,  true };            // in Z direction
  //  DROPS::PoissonCoeffCL<DROPS::ParamCL> PoissonCoeff(P, *PyC->GetDiffusion, *PyC->GetSource, &Zero);
  const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
    { *PyC->GetInflow, &Zero, &Zero, *PyC->GetInterfaceValue, &Zero, &Zero};

  DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);
  //only for measuring cell, not used here
  double r = 1;
  std::string serfile = "none";

  DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);
  // Setup the problem
  //DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), bdata);
  DROPS::PoissonP1CL< PythonConnectCL > prob( *mg, *PyC, bdata);
#ifdef _PAR
  // Set parallel data structures
  DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
  pmg.AttachTo( *mg);                                  // handling of parallel multigrid
  DROPS::LoadBalHandlerCL lb( *mg, DROPS::metis);     // loadbalancing
  lb.DoInitDistribution( DROPS::ProcCL::Master());    // distribute initial grid
  lb.SetStrategy( DROPS::Recursive);                  // best distribution of data
#endif
  timer.Stop();
  *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;
  // Refine the grid
  // ---------------------------------------------------------------------
  *(PyC->outfile) << "Refine the grid " << P.get<int>("DomainCond.RefineSteps") << " times regulary ...\n";
  timer.Reset();
  // Create new tetrahedra
  for ( int ref=1; ref <= P.get<int>("DomainCond.RefineSteps"); ++ref){
    *(PyC->outfile) << " refine (" << ref << ")\n";
    DROPS::MarkAll( *mg);
    mg->Refine();
  }
  // do loadbalancing
#ifdef _PAR
  lb.DoMigration();
#endif
  timer.Stop();
  *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;
  mg->SizeInfo(*(PyC->outfile));
  //prepare MC
  PyC->SetMG(mg);
  PyC->ClearMaps();
  PyC->setFaceMap();
  PyC->setTetraMap();
  // Solve the problem
  DROPS::Strategy( prob, P, PyC);
  delete mg;
}

//used to solve the sensetivity problem in IA12
void CoefCorrection(std::ofstream& outfile, DROPS::ParamCL& P, PdeFunction::ConstPtr C0, PdeFunction::ConstPtr b_in, PdeFunction::ConstPtr b_interface,
                    PdeFunction::ConstPtr source, PdeFunction::ConstPtr Dw, PdeFunction::ConstPtr presol, PdeFunction::ConstPtr Del,double* C_sol)
{
  PythonConnectCL::Ptr PyC(new PythonConnectCL());
  PyC->Init(&outfile, P, C0, b_in, source, Dw, b_interface, C_sol);
  PyC->SetupIA12SensiRhs(presol, Del);
  SetMissingParameters(P);
#ifdef _PAR
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
#endif
  //try
  //{
  // time measurements
#ifndef _PAR
  DROPS::TimerCL timer;
#else
  DROPS::ParTimerCL timer;
#endif

  // set up data structure to represent a poisson problem
  // ---------------------------------------------------------------------
  *(PyC->outfile) << line << "Set up data structure to represent a Poisson problem ...\n";
  timer.Reset();

  //create geometry
  DROPS::MultiGridCL* mg= 0;

  const bool isneumann[6]=
    { false, true,              // inlet, outlet
      true,  false,             // wall, interface
      true,  true };            // in Z direction

  //DROPS::PoissonCoeffCL<DROPS::ParamCL> PoissonCoeff(P, *PyC->GetDiffusion, *PyC->GetSource, *PyC->GetInitial);

  const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
    { *PyC->GetInflow, &Zero, &Zero, *PyC->GetInterfaceValue, &Zero, &Zero};

  DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);

  //only for measuring cell, not used here
  double r = 1;
  std::string serfile = "none";

  DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);

  // Setup the problem
  //DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), bdata);
  std::string adstr ("IA1Adjoint");
  std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
  //DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, PoissonCoeff, bdata, adstr.compare(IAProbstr) == 0);
  DROPS::PoissonP1CL< PythonConnectCL > prob( *mg, *PyC, bdata, adstr.compare(IAProbstr) == 0);

#ifdef _PAR
  // Set parallel data structures
  DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
  pmg.AttachTo( *mg);                                  // handling of parallel multigrid
  DROPS::LoadBalHandlerCL lb( *mg, DROPS::metis);     // loadbalancing
  lb.DoInitDistribution( DROPS::ProcCL::Master());    // distribute initial grid
  lb.SetStrategy( DROPS::Recursive);                  // best distribution of data
#endif
  timer.Stop();
  *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;

  // Refine the grid
  // ---------------------------------------------------------------------
  *(PyC->outfile) << "Refine the grid " << P.get<int>("DomainCond.RefineSteps") << " times regulary ...\n";
  timer.Reset();
  // Create new tetrahedra
  for ( int ref=1; ref <= P.get<int>("DomainCond.RefineSteps"); ++ref){
    *(PyC->outfile) << " refine (" << ref << ")\n";
    DROPS::MarkAll( *mg);
    mg->Refine();
  }
  // do loadbalancing
#ifdef _PAR
  lb.DoMigration();
#endif

  timer.Stop();
  *(PyC->outfile) << " o time " << timer.GetTime() << " s" << std::endl;
  mg->SizeInfo(*(PyC->outfile));

  //prepare MC
  PyC->SetMG(mg);
  PyC->ClearMaps();
  PyC->setFaceMap();
  PyC->setTetraMap();

  // Solve the problem
  DROPS::Strategy( prob, P, PyC);

  delete mg;
}
