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
#include "poisson/poissonCoeff.h"      // Coefficient-Function-Container poissonCoeffCL
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

const char line[] ="----------------------------------------------------------------------------------\n";

PythonConnectCL PyC;

DROPS::ParamCL P;

double GetInitial(const DROPS::Point3DCL& p, double t)
{
  return PyC.GetInitial(p,t);
}
double GetDiffusion(const DROPS::Point3DCL& p, double t){return PyC.GetDiffusion(p,t);}
double GetSource(const DROPS::Point3DCL& p, double t){return PyC.GetSource(p,t);}
double GetInflow(const DROPS::Point3DCL& p, double t){return PyC.GetInflow(p,t);}
double GetInterfaceFlux( const DROPS::Point3DCL& p, double t){return PyC.GetInterfaceFlux(p,t);}
double GetInterfaceValue( const DROPS::Point3DCL& p, double t){return PyC.GetInterfaceValue(p,t);}
//for IA2
double GetPresol( const DROPS::Point3DCL& p, double t){return PyC.GetPresol(p,t);}
double GetDelPsi( const DROPS::Point3DCL& p, double t){return PyC.GetDelPsi(p,t);}



double Zero(const DROPS::Point3DCL&, double) { return 0.0; }
double One(const DROPS::Point3DCL&, double) { return 1.0; }
// -> Here: f = 128*(xy(1-x)(1-y)+xz(1-x)(1-z)+yz(1-y)(1-z))
double Source_test(const DROPS::Point3DCL& p, double){
        return 128*(p[1]*p[2]*(1.-p[1])*(1.-p[2])
                + p[0]*p[2]*(1.-p[0])*(1.-p[2])
                + p[0]*p[1]*(1.-p[0])*(1.-p[1]));
       }
double Inflow(const DROPS::Point3DCL& p, double t) { return PyC.GetInflow(p,t); }
//double Interface(const DROPS::Point3DCL& p, double t) { return PyC.GetInterfaceFlux(p,t); }

namespace DROPS
{


template<class CoeffCL, class SolverT>
void SolveStatProblem( PoissonP1CL<CoeffCL>& Poisson, SolverT& solver, ParamCL& P)
{
    // time measurements
#ifndef _PAR
    TimerCL timer;
    const bool doErrorEstimate= P.get<int>("Err.DoErrorEstimate");
#else
    const bool doErrorEstimate= false;
    if (P.get<int>("Err.DoErrorEstimate"))
        std::cout << "Skipping Error-Estimation ..." << std::endl;
    ParTimerCL timer;
#endif

    if ( !doErrorEstimate) {
        std::string Gradstr ("IA2Gradient");
        std::string Sensstr ("IA2Sense");
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
            Poisson.SetupL2ProjGrad(b1,GetPresol, GetDelPsi, NULL);
            Poisson.b.Data+=b1.Data;
            std::cout << line << "We are solving Gradient problem in IA2 ...\n"<<line;  
        }
        if(SensProb)
        {
            VecDescCL b1;
            b1.SetIdx(&Poisson.idx);
            Poisson.SetupGradSrc(b1,GetPresol, GetDelPsi);
            Poisson.b.Data+=b1.Data;
            std::cout << line << "We are solving sensetivity problem in IA2 ...\n"<<line;  
        }
        if(P.get<int>("Stabilization.SUPG"))
        {
            //CoeffCL::Show_Pec();
            std::cout << line << "The SUPG stabilization has been added ...\n"<<line;           
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
        PyC.SetSol3D(Poisson.GetSolution());
        std::cout << " o Solved system with:\n"
                  << "   - time          " << timer.GetTime()   << " s\n"
                  << "   - iterations    " << solver.GetIter()  << '\n'
                  << "   - residuum      " << solver.GetResid() << '\n'
                  << "   - real residuum " << realresid         << std::endl;
        if(solver.GetResid() > P.get<double>("Poisson.Tol"))
        {  
            std::cout <<"The residual is bigger than tolerence ...\n";
            abort();
        }
        if (P.get<int>("Poisson.SolutionIsKnown")) {
            std::cout << line << "Check result against known solution ...\n";
            Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
        }
    }
    else{
        MultiGridCL& MG= Poisson.GetMG();
        const typename PoissonP1CL<CoeffCL>::BndDataCL& BndData= Poisson.GetBndData();

        MLIdxDescCL  loc_idx;
        VecDescCL  loc_x;
        MLIdxDescCL* new_idx= &Poisson.idx;
        MLIdxDescCL* old_idx= &loc_idx;
        VecDescCL* new_x= &Poisson.x;
        VecDescCL* old_x= &loc_x;

        DoerflerMarkCL<typename PoissonP1CL<CoeffCL>::est_fun, typename PoissonP1CL<CoeffCL>::base_>
            Estimator( P.get<double>("Err.RelReduction"), P.get<double>("Err.MinRatio"), P.get<double>("Err.Threshold"), P.get<double>("Err.Meas"), P.get<int>("Err.DoMark"),
                       &PoissonP1CL<CoeffCL>::ResidualErrEstimator, *static_cast<typename PoissonP1CL<CoeffCL>::base_*>(&Poisson) );

        int step= 0;
        bool new_marks;

        new_idx->SetFE( P1_FE);
        old_idx->SetFE( P1_FE);

        do{
            timer.Reset();
            //std::cout << DROPS::SanityMGOutCL(MG) << std::endl;
            MG.Refine();

            Poisson.CreateNumbering( MG.GetLastLevel(), new_idx);    // create numbering for this idx
            std::cout << "new triangLevel: " << Poisson.idx.TriangLevel() << std::endl;
            Poisson.b.SetIdx( new_idx);                              // tell b about numbering
            new_x->SetIdx( new_idx);                    			 // second vector with the same idx

            std::cout << line << "Problem size\no number of unknowns             " << new_x->Data.size() << std::endl;

            MG.SizeInfo(std::cout);
            if ( step == 0)
                Estimator.Init( typename PoissonP1CL<CoeffCL>::DiscSolCL( new_x, &BndData, &MG));

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
            PyC.SetSol3D(Poisson.GetSolution());
            double realresid = norm( VectorCL(Poisson.A.Data*new_x->Data-Poisson.b.Data));
            std::cout << " o Solved system with:\n"
                      << "   - time          " << timer.GetTime()   << " s\n"
                      << "   - iterations    " << solver.GetIter()  << '\n'
                      << "   - residuum      " << solver.GetResid() << '\n'
                      << "   - real residuum " << realresid         << std::endl;
            if(solver.GetResid() > P.get<double>("Poisson.Tol"))
            {  
                std::cout <<"The residual is bigger than tolerence ...\n";
                abort();
            }
            Poisson.A.Reset();
            Poisson.b.Reset();
            if (P.get<int>("Poisson.SolutionIsKnown")) {
                std::cout << line << "Check result against known solution ...\n";
                Poisson.CheckSolution( *new_x, CoeffCL::Solution);
            }
            new_marks = Estimator.Estimate( typename PoissonP1CL<CoeffCL>::const_DiscSolCL( new_x, &BndData, &MG) );

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
    }
}

/// \brief Strategy to solve the Poisson problem on a given triangulation
  template<class CoeffCL>
void Strategy( PoissonP1CL<CoeffCL>& Poisson, ParamCL& P)
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
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
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
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;


    // display problem size
    // -------------------------------------------------------------------------
    std::cout << line << "Problem size\n";
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
    std::cout << " o number of unknowns             " << numUnk    << '\n'
              << " o number of accumulated unknowns " << numAccUnk << '\n'
              << " o number of unknowns on proc\n";
    for (size_t i=0; i<UnkOnProc.size(); ++i)
        std::cout << " - Proc " << i << ": " << UnkOnProc[i]<< '\n';

    // discretize (setup linear equation system)
    // -------------------------------------------------------------------------
    std::cout << line << "Discretize (setup linear equation system) ...\n";

    timer.Reset();
    if (P.get<int>("Time.NumSteps") != 0)
        Poisson.SetupInstatSystem( Poisson.A, Poisson.M, Poisson.x.t, P.get<int>("Stabilization.SUPG") );
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;


    // solve the linear equation system
    // -------------------------------------------------------------------------
    std::cout << line << "Solve the linear equation system ...\n";

    // type of preconditioner and solver
    PoissonSolverFactoryCL<> factory( P, Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    if ( factory.GetProlongation() != 0)
        SetupP1ProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);

    // Solve the linear equation system
    if(P.get<int>("Time.NumSteps") !=0)
        Poisson.Init( Poisson.x, CoeffCL::InitialCondition, 0.0);
    else
        SolveStatProblem( Poisson, *solver, P);

    if(P.get<int>("Time.NumSteps")!=0)
    {
        //CoeffCL::Show_Pec();
        InstatPoissonThetaSchemeCL<PoissonP1CL<CoeffCL>, PoissonSolverBaseCL>
        ThetaScheme( Poisson, *solver, P.get<double>("Time.Theta") , P.get<int>("PoissonCoeff.Convection"), P.get<int>("Stabilization.SUPG"));
        ThetaScheme.SetTimeStep(P.get<double>("Time.StepSize") );
        for ( int step = 1; step <= P.get<int>("Time.NumSteps") ; ++step) {
            timer.Reset();
            std::cout << line << "Step: " << step << std::endl;
            ThetaScheme.DoStep( Poisson.x);
            //Setup solutions for python interface
            PyC.SetSol3D(Poisson.GetSolution(), Poisson.x.t);

            timer.Stop();
            std::cout << " o Solved system with:\n"
                      << "   - time          " << timer.GetTime()    << " s\n"
                      << "   - iterations    " << solver->GetIter()  << '\n'
                      << "   - residuum      " << solver->GetResid() << '\n';
            if(solver->GetResid() > P.get<double>("Poisson.Tol"))
            {  
                std::cout <<"The residual is bigger than tolerence ...\n";
                abort();
            }

            if (P.get("Poisson.SolutionIsKnown", 0)) {
                std::cout << line << "Check result against known solution ...\n";
                Poisson.CheckSolution( Poisson.x, CoeffCL::Solution, Poisson.x.t);
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
void convection_diffusion(DROPS::ParamCL& P, const PdeFunction* C0, const PdeFunction* b_in, const PdeFunction* b_interface, const PdeFunction* source, const PdeFunction* Dw, double* C_sol)
{
        PyC.Init(P, C0, b_in, source, Dw, b_interface, C_sol);
        SetMissingParameters(P);
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
#endif
    try
    {
    // time measurements
#ifndef _PAR
        DROPS::TimerCL timer;
#else
        DROPS::ParTimerCL timer;
#endif

        // set up data structure to represent a poisson problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Poisson problem ...\n";
        timer.Reset();

        //create geometry
        DROPS::MultiGridCL* mg= 0;
        //DROPS::PoissonBndDataCL* bdata = 0;

/*        const bool isneumann[6]=
        { false, true,              // inlet, outlet
          true,  false,             // wall, interface
          true,  true };            // in Z direction*/

        //for testing
        const bool isneumann[6]=
        { false, false,              // inlet, outlet
          false,  false,             // wall, interface
          false,  false };            // in Z direction

        DROPS::PoissonCoeffCL<DROPS::ParamCL> PoissonCoeff(P, GetDiffusion, GetSource, GetInitial);
        //DROPS::PoissonCoeffCL<DROPS::ParamCL> PoissonCoeff(P);

        //for testing
/*        const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { &Inflow, &Zero, &Zero, &GetInterfaceValue, &Zero, &Zero};*/

        const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { &Inflow, &One, &One, &GetInterfaceValue, &One, &One};

        DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);

        //only for measuring cell, not used here
        double r = 1;
        std::string serfile = "none";

        DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);
        
        // Setup the problem
        //DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), bdata);
        std::string adstr ("IA1Adjoint");
        std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
        DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, PoissonCoeff, bdata, adstr.compare(IAProbstr) == 0);

#ifdef _PAR
        // Set parallel data structures
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
        pmg.AttachTo( *mg);                                  // handling of parallel multigrid
        DROPS::LoadBalHandlerCL lb( *mg, DROPS::metis);     // loadbalancing
        lb.DoInitDistribution( DROPS::ProcCL::Master());    // distribute initial grid
        lb.SetStrategy( DROPS::Recursive);                  // best distribution of data
#endif
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        // Refine the grid
        // ---------------------------------------------------------------------
        std::cout << "Refine the grid " << P.get<int>("DomainCond.RefineSteps") << " times regulary ...\n";
        timer.Reset();
        // Create new tetrahedra
        for ( int ref=1; ref <= P.get<int>("DomainCond.RefineSteps"); ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( *mg);
            mg->Refine();
        }
        // do loadbalancing
#ifdef _PAR
        lb.DoMigration();
#endif

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        mg->SizeInfo(cout);

        //prepare MC
        PyC.SetMG(mg);
        PythonConnectCL::ClearMaps();
        PythonConnectCL::setFaceMap();
        PythonConnectCL::setTetraMap();

        // Solve the problem
	    DROPS::Strategy( prob, P);
        //std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;

        delete mg;
        //delete bdata;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}

int main(int argc, char** argv)
{
  try {
    using namespace DROPS;

    std::ifstream param;
    if (argc!=2){
      std::cout << "Using default parameter file: poissonex1.json\n";
      param.open( "poissonex1.json");
    }
    else
      param.open( argv[1]);
    if (!param){
      std::cerr << "error while opening parameter file\n";
      return 1;
    }
    param >> P;
    param.close();
    std::cout << P << std::endl;
    
    
    
    
    instat_scalar_fun_ptr initial;
    instat_scalar_fun_ptr rhs;
    instat_scalar_fun_ptr bnd;
    instat_scalar_fun_ptr diffusion;
    initial = One;
    bnd     = One;
    rhs     = Source_test;
    diffusion = Zero;
    // set up data structure to represent a poisson problem
    // ---------------------------------------------------------------------
    typedef const PdeFunction* PdeFunPtr;
    PdeFunPtr C0(new TestPdeFunction(P,initial));
    PdeFunPtr b_in(new TestPdeFunction(P,bnd));
    PdeFunPtr b_interface(new TestPdeFunction(P, bnd));
    PdeFunPtr source(new TestPdeFunction(P, rhs));
    PdeFunPtr Dw(new TestPdeFunction(P, diffusion));
    
    double lx, ly, lz;
    int nx, ny, nz;
    int refinesteps= P.get<int>("DomainCond.RefineSteps");
    
    std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
    size_t idx_;
    while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx_]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> lx >> ly>> lz >> nx >> ny >> nz;
    int Nx = nx * pow (2, refinesteps)+1;
    int Ny = ny * pow (2, refinesteps)+1;
    int Nz = nz * pow (2, refinesteps)+1;
    int Nxyz= Nx*Ny*Nz;
    double* C_sol;
    if(P.get<double>("Time.NumSteps")==0)
        C_sol=new double[Nxyz];
    else
    {
      int nt = P.get<double>("Time.NumSteps");
      int Ns= Nxyz*nt;
        C_sol=new double[Ns];
    }
    

    convection_diffusion(P, C0, b_in, b_interface, source, Dw, C_sol);

    delete[] C0;
    delete[] b_in;
    delete[] b_interface;
    delete[] source;
    delete[] C_sol;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
