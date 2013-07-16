/// \file poissonP1.cpp
/// \brief Solver for Poisson problem with P1 functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Yuanjun Zhang, Thorolf Schulte, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

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
#include "out/output.h"
#include "geom/geomselect.h"
#include "misc/bndmap.h"

// include numeric computing!
#include "num/fe.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include "poisson/integrTime.h"

 // include problem class
#include "misc/params.h"
#include "poisson/poissonCoeff.h"      // Coefficient-Function-Container poissonCoeffCL
#include "poisson/poisson.h"           // setting up the Poisson problem
#include "num/bndData.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

//include output
#include "out/ensightOut.h"
#include "out/vtkOut.h"

// include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "num/parsolver.h"              // various parallel solvers
#include "num/parprecond.h"             // various parallel preconditioners
#endif

// include function container
#include "misc/bndmap.h"

#include "num/poissonsolverfactory.h"

#include "poisson/ale.h"

using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

DROPS::ParamCL P;   //Parameter class, read in json file in main function

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
        // discretize (setup linear equation system)
        std::cout << line << "Discretize (setup linear equation system) in stationary problem...\n";
        timer.Reset();
        Poisson.SetupSystem( Poisson.A, Poisson.b, P.get<int>("Stabilization.SUPG"));
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        //If we need to add convection
        if(P.get<int>("PoissonCoeff.Convection"))
        {
            std::cout << line << "Setup convection...\n";
            timer.Reset();
            Poisson.vU.SetIdx( &Poisson.idx);
            Poisson.SetupConvection(Poisson.U, Poisson.vU, 0.0);
            Poisson.A.Data.LinComb(1., Poisson.A.Data, 1., Poisson.U.Data);
            Poisson.b.Data+=Poisson.vU.Data;
            timer.Stop();
            std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        }
        timer.Reset();
        solver.Solve( Poisson.A.Data, Poisson.x.Data, Poisson.b.Data);
        timer.Stop();

#ifndef _PAR
        double realresid = norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data));
#else
        double realresid = Poisson.idx.GetEx().Norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data), false);
#endif
        std::cout << " o Solved system with:\n"
                  << "   - time          " << timer.GetTime()   << " s\n"
                  << "   - iterations    " << solver.GetIter()  << '\n'
                  << "   - residuum      " << solver.GetResid() << '\n'
                  << "   - real residuum " << realresid         << std::endl;

        if (P.get<int>("Poisson.SolutionIsKnown")) {
            std::cout << line << "Check result against known solution ...\n";
            timer.Reset();
            Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
            timer.Stop();
            std::cout << " o time " << timer.GetTime() << " s" << std::endl;
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

        DoerflerMarkCL<typename PoissonP1CL<CoeffCL>::est_fun, typename PoissonP1CL<CoeffCL>::base_>
            Estimator( P.get<double>("Err.RelReduction"), P.get<double>("Err.MinRatio"), P.get<double>("Err.Threshold"), P.get<double>("Err.Meas"), P.get<int>("Err.DoMark"),
                       &PoissonP1CL<CoeffCL>::ResidualErrEstimator, *static_cast<typename PoissonP1CL<CoeffCL>::base_*>(&Poisson) );

        int step= 0;
        bool new_marks;

        new_idx->SetFE( P1_FE);
        old_idx->SetFE( P1_FE);

        do{
            timer.Reset();
            std::cout << DROPS::SanityMGOutCL(MG) << std::endl;
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
            double realresid = norm( VectorCL(Poisson.A.Data*new_x->Data-Poisson.b.Data));
            std::cout << " o Solved system with:\n"
                      << "   - time          " << timer.GetTime()   << " s\n"
                      << "   - iterations    " << solver.GetIter()  << '\n'
                      << "   - residuum      " << solver.GetResid() << '\n'
                      << "   - real residuum " << realresid         << std::endl;
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
    }*/
}

/// \brief Strategy to solve the Poisson problem on a given triangulation
template<class CoeffCL>
void Strategy( PoissonP1CL<CoeffCL>& Poisson)
{
    // time measurements
#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif

    // the triangulation
    MultiGridCL& mg= Poisson.GetMG();
    ALECL ALE(P, mg);
    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    Poisson.idx.SetFE( P1_FE);                                  // set quadratic finite elements
    Poisson.vel_idx.SetFE( vecP1_FE);
    //see class for explanation: template didnt work
    if ( PoissonSolverFactoryHelperCL().MGUsed(P))
        Poisson.SetNumLvl ( mg.GetNumLevel());
    Poisson.CreateNumbering( mg.GetLastLevel(), &Poisson.idx);      // number vertices and edges;
    Poisson.CreateVelNumbering( mg.GetLastLevel(), &Poisson.vel_idx);      // number vertices and edges;
    Poisson.b.SetIdx( &Poisson.idx);                            // tell b about numbering
    Poisson.x.SetIdx( &Poisson.idx);                            // tell x about numbering
    Poisson.velocity.SetIdx( &Poisson.vel_idx);                            // tell x about numbering
    Poisson.A.SetIdx( &Poisson.idx, &Poisson.idx);              // tell A about numbering
    Poisson.M.SetIdx( &Poisson.idx, &Poisson.idx);              // tell M about numbering
    Poisson.U.SetIdx( &Poisson.idx, &Poisson.idx);              // tell U about numbering

    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

    // display problem size
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

    std::cout << line << "Choose the poisson solver...\n";
    timer.Reset();
    // type of preconditioner and solver
    PoissonSolverFactoryCL<> factory( P, Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    if ( factory.GetProlongation() != 0)
        SetupP1ProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);

    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

    //If it is NOT a stationary problem, set up the system and the initial condition
    if (P.get<int>("Time.NumSteps") != 0)
    {
        // discretize (setup linear equation system)
        std::cout << line << "Discretize (setup linear equation system) for instationary problem...\n";
        timer.Reset();
        if(Poisson.ALE_)
            ALE.InitGrid();
        Poisson.SetupInstatSystem( Poisson.A, Poisson.M, Poisson.x.t);
        Poisson.Init( Poisson.x, CoeffCL::InitialCondition, 0.0);
        Poisson.SetupVel(Poisson.velocity, CoeffCL::Vel, 0.0);
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
    }
    else
    {//if it is a stationary problem, call SolveStatProblem
        std::cout << line << "Solve the linear equation system ...\n";
        SolveStatProblem( Poisson, *solver, P);
    }

    // Output-Registrations:
    //Ensight format
    Ensight6OutCL* ensight = NULL;
    if (P.get<int>("Ensight.EnsightOut",0)){
        // Initialize Ensight6 output
        const std::string filename= P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase");
        ensight = new Ensight6OutCL(P.get<std::string>("Ensight.EnsCase")+".case", P.get<int>("Time.NumSteps")+1,
                                    P.get<int>("Ensight.Binary"), P.get<int>("Ensight.MasterOut"));
        ensight->Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(),
                                                P.get<std::string>("Ensight.GeomName"), filename + ".geo"));
        ensight->Register( make_Ensight6Scalar( Poisson.GetSolution(), "Temperatur", filename + ".tp", true));
        ensight->Write();
    }
    //VTK format
    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.VTKOut",0)){
        vtkwriter = new VTKOutCL(mg, "DROPS data",
                                 P.get<int>("Time.NumSteps")+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"),
                                 P.get<int>("VTK.Binary"), true );
        vtkwriter->Register( make_VTKVector( Poisson.GetVelocity(), "Velocity"));
        vtkwriter->Register( make_VTKScalar( Poisson.GetSolution(), "ConcenT"));
        vtkwriter->Write( Poisson.x.t);
    }
    //Do we have an instationary problem?
    if(P.get<int>("Time.NumSteps")!=0)
    {
        //Creat instationary ThetaschemeCL to handle time integration for instationary problem and set time steps
        InstatPoissonThetaSchemeCL<PoissonP1CL<CoeffCL>, PoissonSolverBaseCL>
        ThetaScheme( Poisson, *solver, P);
        ThetaScheme.SetTimeStep(P.get<double>("Time.StepSize") );
        //Solve linear systerm in each time step
        for ( int step = 1; step <= P.get<int>("Time.NumSteps") ; ++step) {
            timer.Reset();
            std::cout << line << "Step: " << step << std::endl;
            if(Poisson.ALE_)
                ALE.MovGrid(Poisson.x.t);
            ThetaScheme.DoStep( Poisson.x);

            timer.Stop();
            std::cout << " o Solved system with:\n"
                      << "   - time          " << timer.GetTime()    << " s\n"
                      << "   - iterations    " << solver->GetIter()  << '\n'
                      << "   - residuum      " << solver->GetResid() << '\n';


        Poisson.SetupVel( Poisson.velocity, CoeffCL::Vel, Poisson.x.t);
            if (P.get("Poisson.SolutionIsKnown", 0)) {
                std::cout << " o Check result against known solution ...\n";
                timer.Reset();
                Poisson.CheckSolution( Poisson.x, CoeffCL::Solution, Poisson.x.t);
                timer.Stop();
                std::cout << " o -time " << timer.GetTime() << " s" << std::endl;
            }

            if (ensight && step%P.get<int>("Ensight.EnsightOut", 0)==0){
                std::cout << " o Ensight output ...\n";
                timer.Reset();
                ensight->Write( Poisson.x.t);
                timer.Stop();
                std::cout << " o -time " << timer.GetTime() << " s" << std::endl;
            }
            if (vtkwriter && step%P.get<int>("VTK.VTKOut", 0)==0){
                std::cout << " o VTK output ...\n";
                timer.Reset();
                vtkwriter->Write( Poisson.x.t);
                timer.Stop();
                std::cout << " o -time " << timer.GetTime() << " s" << std::endl;
            }
        }
    }

    if (vtkwriter) delete vtkwriter;
    if (ensight) delete ensight;
    delete solver;
}

} // end of namespace DROPS

/// \brief Set Default parameters here s.t. they are initialized.
/// The result can be checked when Param-list is written to the output.
void SetMissingParameters(DROPS::ParamCL& P){
    P.put_if_unset<int>("Stabilization.SUPG",0);
    P.put_if_unset<double>("Stabilization.Magnitude",1.0);
    P.put_if_unset<int>("Stabilization.Grids",1);
    P.put_if_unset<int>("ALE.wavy",0);
    P.put_if_unset<std::string>("ALE.Interface","Zero");
}

int main (int argc, char** argv)
{
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

        std::ifstream param;
        if (argc!=2){
            std::cout << "Using default parameter file: statpoissonEx.json\n";
            param.open( "statpoissonEx.json");
        }
        else
            param.open( argv[1]);
        if (!param){
            std::cerr << "error while opening parameter file\n";
            return 1;
        }
        //output all the parameters
        param >> P;
        param.close();
        //Setup missing parameters
        SetMissingParameters(P);
        std::cout << P << std::endl;

        // set up data structure to represent a poisson problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Poisson problem ...\n";
        timer.Reset();

        //create geometry
        DROPS::MultiGridCL* mg= 0;
        DROPS::PoissonBndDataCL* bdata = 0;

        //only for measuring cell, not used here
        double r = 1;
        std::string serfile = "none";
        //build computational domain
        DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);
        //Setup boundary conditions
        DROPS::BuildBoundaryData( mg, bdata, P.get<std::string>("DomainCond.BoundaryType"), P.get<std::string>("DomainCond.BoundaryFncs"));
        //Initialize SUPGCL class
        DROPS::SUPGCL supg;
        if(P.get<int>("Stabilization.SUPG"))
        {
            supg.init(P);
            std::cout << line << "The SUPG stabilization will be added ...\n"<<line;
        }
        // Setup the problem
        DROPS::PoissonCoeffCL tmp = DROPS::PoissonCoeffCL( P);
        DROPS::PoissonP1CL<DROPS::PoissonCoeffCL> prob( *mg, tmp, *bdata, supg, P.get<int>("ALE.wavy"));

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

        // Solve the problem
        DROPS::Strategy(prob);
        //Check if Multigrid is sane
        std::cout << line << "Check if multigrid works properly...\n";
        timer.Reset();
        if(P.get<int>("ALE.wavy"))
            std::cout << "Because of ALE method, we don't check the sanity of multigrid here!" << std::endl;
        else
            std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        // delete dynamically allocated objects
        delete mg;
        delete bdata;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
