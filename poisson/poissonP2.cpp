/// \file poissonP2.cpp
/// \brief Solver for Poisson problem with P2 functions
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

/** We solve \f$ -\Delta u = f\;\mbox{in}\; \Omega:=[0,1]^3 \f$ for the given
    solution \f$ u(x,y,z):= 64 \cdot xyz (1-x) (1-y) (1-z) \f$, i.e. homogeneous
    Dirichlet conditions are used. A uniform tetrahedral grid is applied as
    a triangulation of \f$ \Omega \f$. GMRES is used as a linear solver for the
    discretized linear equation system. Note, that CG-type methods can be used
    as well because the resulting linear equation system is s.p.d. However,
    since this program acts as a base performance test, GMRES is used here.
*/

 // include parallel computing!
#ifdef _PAR
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-measurement
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "num/parsolver.h"              // various parallel solvers
#include "num/parprecond.h"             // various parallel preconditioners
#endif

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "geom/geomselect.h"
#include "misc/bndmap.h"

 // include outputs
#include "out/output.h"
#include "out/vtkOut.h"

 // include numeric computing!
#include "num/fe.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include "poisson/integrTime.h"
#include "num/poissonsolverfactory.h"

 // include problem class
#include "misc/params.h"
#include "poisson/poisson.h"      // setting up the Poisson problem
#include "poisson/poissonCoeff.h"      // Coefficient-Function-Container poissonCoeffCL
#include "num/bndData.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

DROPS::ParamCL P;

namespace DROPS
{


/// \brief Strategy to solve the Poisson problem on a given triangulation
template<class CoeffCL>
void Strategy( PoissonP2CL<CoeffCL>& Poisson, SUPGCL& supg)
{
    // time measurement
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

    Poisson.idx.SetFE( P2_FE);                                  // set quadratic finite elements
    if ( PoissonSolverFactoryHelperCL().MGUsed(P))
        Poisson.SetNumLvl ( mg.GetNumLevel());
    Poisson.CreateNumbering( mg.GetLastLevel(), &Poisson.idx);  // number vertices and edges
    Poisson.b.SetIdx( &Poisson.idx);                            // tell b about numbering
    Poisson.x.SetIdx( &Poisson.idx);                            // tell x about numbering
    Poisson.A.SetIdx( &Poisson.idx, &Poisson.idx);              // tell A about numbering
    Poisson.M.SetIdx( &Poisson.idx, &Poisson.idx);              // tell M about numbering
    Poisson.U.SetIdx( &Poisson.idx, &Poisson.idx);              // tell U about numbering

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
    
    if(P.get<int>("Time.NumSteps") !=0)
        Poisson.SetupInstatSystem(Poisson.A, Poisson.M, 0., supg);    //IntationarySystem
    else
    {
        Poisson.SetupSystem( Poisson.A, Poisson.b, supg);         //StationarySystem
        if(P.get<int>("PoissonCoeff.Convection"))
          {
            Poisson.vU.SetIdx( &Poisson.idx); 
            Poisson.SetupConvection(Poisson.U, Poisson.vU, 0.0, supg);                 //Setupconvection
            Poisson.A.Data.LinComb(1., Poisson.A.Data, 1., Poisson.U.Data); //Combination with convection
            Poisson.b.Data+=Poisson.vU.Data;
          }
    }
        
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;


    // solve the linear equation system
    // -------------------------------------------------------------------------
    std::cout << line << "Solve the linear equation system ...\n";

    // type of preconditioner and solver
    PoissonSolverFactoryCL<> factory( P, Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    timer.Reset();
    if ( factory.GetProlongation() != 0)
        SetupP2ProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);
    
    if(P.get<int>("Time.NumSteps") !=0)
         Poisson.Init(Poisson.x, CoeffCL::InitialCondition, 0.0);
         
   
 

   //Solve stationary problem
   if(P.get<int>("Time.NumSteps") ==0)
    {          
        timer.Reset();
        solver->Solve( Poisson.A.Data, Poisson.x.Data, Poisson.b.Data);
        timer.Stop();
    double realresid;
#ifndef _PAR
    realresid= norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data));
#else
    realresid= Poisson.idx.GetEx().Norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data),false);
#endif

    std::cout << " o Solved system with:\n"
              << "   - time          " << timer.GetTime()    << " s\n"
              << "   - iterations    " << solver->GetIter()  << '\n'
              << "   - residuum      " << solver->GetResid() << '\n'
              << "   - real residuum " << realresid          << std::endl;
        if (P.get<int>("Poisson.SolutionIsKnown")) {
        std::cout << line << "Check result against known solution ...\n";
        Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
        }
    }
    


    // Output-Registrations:
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


    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.VTKOut",0)){
        vtkwriter = new VTKOutCL(mg, "DROPS data", 
                                 P.get<int>("Time.NumSteps")+1, 
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"), 
                                 P.get<int>("VTK.Binary") );
        vtkwriter->Register( make_VTKScalar( Poisson.GetSolution(), "ConcenT"));
        vtkwriter->Write( Poisson.x.t);
    }

    if (P.get<int>("Time.NumSteps") != 0){
        InstatPoissonThetaSchemeCL<PoissonP2CL<CoeffCL>, PoissonSolverBaseCL>
           ThetaScheme(Poisson, *solver, supg, P.get<double>("Time.Theta"), P.get<double>("PoissonCoeff.Convection"));
        ThetaScheme.SetTimeStep(P.get<double>("Time.StepSize"));
   
        for ( int step = 1; step <= P.get<int>("Time.NumSteps"); ++step)
        {
            timer.Reset();

            std::cout << line << "Step: " << step << std::endl;
            ThetaScheme.DoStep( Poisson.x);

            timer.Stop();
            std::cout << " o Solved system with:\n"
                    << "   - time          " << timer.GetTime()    << " s\n"
                    << "   - iterations    " << solver->GetIter()  << '\n'
                    << "   - residuum      " << solver->GetResid() << '\n';

            // check the result
            // -------------------------------------------------------------------------
            if (P.get("Poisson.SolutionIsKnown", 0)) {
                std::cout << line << "Check result against known solution ...\n";
                Poisson.CheckSolution( Poisson.x, CoeffCL::Solution, Poisson.x.t);
            }

            if (ensight && step%P.get<int>("Ensight.EnsightOut", 0)==0)
                ensight->Write( Poisson.x.t);
            if (vtkwriter && step%P.get<int>("VTK.VTKOut", 0)==0)
                vtkwriter->Write( Poisson.x.t);
        }
    }

    if (vtkwriter) delete vtkwriter;
    if (ensight) delete ensight;
    delete solver;
}

} // end of namespace DROPS

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
#endif
    try
    {

        std::ifstream param;
        if (argc!=2){
            std::cout << "Using default parameter file: PoissonEx.json\n";
            param.open( "PoissonEx.json");
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

        // time measurement
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
        DROPS::PoissonBndDataCL* bdata = 0;

        //only for measuring cell, not used here
        double r = 1;
        std::string serfile = "none";

        DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);
        DROPS::BuildBoundaryData( mg, bdata, P.get<std::string>("DomainCond.BoundaryType"), P.get<std::string>("DomainCond.BoundaryFncs"));
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        mg->SizeInfo(cout);
        std::cout << line << "Set up load balancing ...\n";
        // Setup the problem
        DROPS::PoissonP2CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), *bdata);
        timer.Reset();
#ifdef _PAR
        // Set parallel data structures
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
        pmg.AttachTo( *mg);                                  // handling of parallel multigrid
        DROPS::LoadBalHandlerCL lb( *mg, DROPS::metis);                    // loadbalancing
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
        for ( int ref=1; ref<=P.get<int>("DomainCond.RefineSteps"); ++ref){
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
        DROPS::SUPGCL supg;
        // Solve the problem
        DROPS::Strategy( prob, supg);
        std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;

        // maple/geomview-output
//        DROPS::RBColorMapperCL colormap;
//        std::ofstream maple("maple.txt");
//        DROPS::Point3DCL e3(0.0); e3[2]= 1.0;
//        maple << DROPS::MapleMGOutCL(*mg, -1, false, true, DROPS::PlaneCL(e3, 0.6)) << std::endl;
//        std::ofstream fil("geom.off");
//        fil << DROPS::GeomSolOutCL<DROPS::PoissonP1CL<PoissonCoeffCL<DROPS::Params> >::DiscSolCL>( *mg, prob.GetSolution(), &colormap, -1, false, 0.0, prob.x.Data.min(), prob.x.Data.max()) << std::endl;
//        std::cout << DROPS::GeomMGOutCL(*mg, -1, true) << std::endl;
        delete mg;
        delete bdata;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
