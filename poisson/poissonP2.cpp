/// \file poissonP2.cpp
/// \brief Solver for Poisson problem with P2 functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier
/// \todo: rename this to poissonP2.cpp as soon as repository has moved to git
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
#include "poisson/params.h"
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

DROPS::ParamPoissonProblemCL C;

namespace DROPS
{

typedef ParamPoissonProblemCL Params;

/// \brief Strategy to solve the Poisson problem on a given triangulation
template<class CoeffCL>
void Strategy( PoissonP2CL<CoeffCL>& Poisson)
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
    if ( PoissonSolverFactoryHelperCL<Params>().MGUsed(C))
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
    
    if(C.tm_NumSteps !=0)
        Poisson.SetupInstatSystem(Poisson.A, Poisson.M);    //IntationarySystem
    else
    {
        Poisson.SetupSystem( Poisson.A, Poisson.b);         //StationarySystem
        if(C.tm_Convection)
          {
            Poisson.vU.SetIdx( &Poisson.idx); 
            Poisson.SetupConvection(Poisson.U, Poisson.vU, 0.0);                 //Setupconvection
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
    PoissonSolverFactoryCL< Params> factory( C, Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    timer.Reset();
    if ( factory.GetProlongation() != 0)
        SetupP2ProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);
    
    if(C.tm_NumSteps !=0)
         Poisson.Init(Poisson.x, CoeffCL::InitialCondition, 0.0);
         
   
 

   //Solve stationary problem
   if(C.tm_NumSteps ==0)
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
        if (C.pos_SolutionIsKnown) {
        std::cout << line << "Check result against known solution ...\n";
        Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
        }
    }
    
    //write for ensight-format
/*    Ensight6OutCL  ens(C.ens_EnsCase+".case", C.tm_NumSteps+1, C.ens_Binary, C.ens_MasterOut);
    if ( C.ens_EnsightOut){
        const std::string filename= C.ens_EnsDir + "/" + C.ens_EnsCase;
        ens.Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(), C.ens_GeomName,       filename + ".geo"));
        ens.Register( make_Ensight6Scalar( Poisson.GetSolution(), "Temperatur", filename + ".tp", true));
        ens.Write();
    }*/
    
        Ensight6OutCL  ens(C.ens_EnsCase+".case", C.tm_NumSteps+1, C.ens_Binary, C.ens_MasterOut);
    if ( C.ens_EnsightOut){
        const std::string filename= C.ens_EnsDir + "/" + C.ens_EnsCase;
        ens.Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(), C.ens_GeomName,       filename + ".geo"));
        ens.Register( make_Ensight6Scalar( Poisson.GetSolution(), "Temperatur", filename + ".tp", true));
        ens.Write();
    }
    
              
    //write for vtk-format
    VTKOutCL vtkwriter(mg, "DROPS data", C.tm_NumSteps+1, std:: string(C.vtk_VTKDir+"/"+C.vtk_VTKName), C.vtk_Binary );//??
    if (C.vtk_VTKOut){
        vtkwriter.Register( make_VTKScalar( Poisson.GetSolution(), "ConcenT"));
        vtkwriter.Write( Poisson.x.t);
    }
   
    if (C.tm_NumSteps != 0){
        InstatPoissonThetaSchemeCL<PoissonP2CL<CoeffCL>, PoissonSolverBaseCL>
           ThetaScheme(Poisson, *solver, C.tm_Theta, C.tm_Convection);
        ThetaScheme.SetTimeStep(C.tm_StepSize);
   
        for ( int step = 1; step <= C.tm_NumSteps; ++step) 
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
            if (C.pos_SolutionIsKnown) {
                std::cout << line << "Check result against known solution ...\n";
                Poisson.CheckSolution( Poisson.x, CoeffCL::Solution, Poisson.x.t);
            }

            if ( C.vtk_VTKOut && step%C.vtk_VTKOut==0)
                vtkwriter.Write( Poisson.x.t);
        }
    }
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
        if (argc!=2)
        {
           std::cout << "Using default parameter file: poissonex1.param\n";
           param.open( "poissonex1.param");
        }
        else
           param.open( argv[1]);
        if (!param)
        {
           std::cerr << "error while opening parameter file\n";
           return 1;
        }
        param >> C;
        param.close();
        std::cout << C << std::endl;

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

        DROPS::BuildDomain( mg, C.dmc_MeshFile, C.dmc_GeomType, serfile, r);
        DROPS::BuildBoundaryData( mg, bdata, C.dmc_BoundaryType, C.dmc_BoundaryFncs);
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        mg->SizeInfo(cout);
        std::cout << line << "Set up load balancing ...\n";
        // Setup the problem
        DROPS::PoissonP2CL<DROPS::PoissonCoeffCL<DROPS::Params> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::Params>(C), *bdata);
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
        std::cout << "Refine the grid " << C.dmc_InitialCond << " times regulary ...\n";
        timer.Reset();

        // Create new tetrahedra
        for ( int ref=1; ref<=C.dmc_InitialCond; ++ref){
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
        DROPS::Strategy( prob);
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
