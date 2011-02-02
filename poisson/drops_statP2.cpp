/// \file drops_statP2.cpp
/// \brief Solver for Poisson problem with P2 functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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
#include "parallel/partime.h"           // parallel time-messurement
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

/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/

template<class ParamsT>
class PoissonCoeffCL
{
  public:
    PoissonCoeffCL(ParamsT&) {}
    /// \brief Reaction
    static double q(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL&, double)
         { return DROPS::Point3DCL(0.); } // no convection
    /// \brief Right-hand side
    static double f(const DROPS::Point3DCL& p, double= 0.0){
        return 128*(p[1]*p[2]*(1.-p[1])*(1.-p[2])
                + p[0]*p[2]*(1.-p[0])*(1.-p[2])
                + p[0]*p[1]*(1.-p[0])*(1.-p[1]));
    }
    /// \brief Diffusion
    static double alpha(const DROPS::Point3DCL&, double=0.0){
        return 1.;
    }
    /// \brief Solution
    static double Solution( const DROPS::Point3DCL& p, double=0.0){
        return 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
    }
};

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

    // writer for vtk-format
    VTKOutCL vtkwriter(mg, "DROPS data", 1, std::string(C.vtk_VTKDir + "/" + C.vtk_VTKName), C.vtk_Binary);
    vtkwriter.Register( make_VTKScalar( Poisson.GetSolution(), "velocity") );

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
    Poisson.SetupSystem( Poisson.A, Poisson.b);
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

    // Solve the linear equation system
    solver->Solve( Poisson.A.Data, Poisson.x.Data, Poisson.b.Data);
    timer.Stop();
    double realresid;
#ifndef _PAR
    realresid= norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data));
#else
    realresid= Poisson.idx.GetEx().Norm( VectorCL(Poisson.A.Data*Poisson.x.Data-Poisson.b.Data),false);
#endif

    std::cout << " o Solved system with:\n"
              << "   - time          " << timer.GetTime()   << " s\n"
              << "   - iterations    " << solver->GetIter()  << '\n'
              << "   - residuum      " << solver->GetResid() << '\n'
              << "   - real residuum " << realresid         << std::endl;
    if (C.pos_SolutionIsKnown) {
        std::cout << line << "Check result against known solution ...\n";
        Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
    }

    if ( C.vtk_VTKOut){
        std::cout << line << "Write solution as VTK file" << std::endl;
        vtkwriter.Write(0.0, false);
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
           std::cout << "Using default parameter file: drops.param\n";
           param.open( "drops.param");
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
        DROPS::PoissonP2CL<PoissonCoeffCL<DROPS::Params> > prob( *mg, PoissonCoeffCL<DROPS::Params>(C), *bdata);
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
