/// \file drops_statP2.cpp
/// \brief Solver for Poisson problem with P2 functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
#include "geom/builder.h"               // construuct the initial multigrid

 // include numeric computing!
#include "num/fe.h"
#include "num/solver.h"

 // include problem class
#include "poisson/poisson.h"      // setting up the Poisson problem
#include "num/bndData.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/
class PoissonCoeffCL
{
  public:
    /// \brief Reaction
    static double q(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    /// \brief Convection
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL, double){
        DROPS::Point3DCL ret; //ret[0]=1.; ret[1]=1.; ret[2]=1.;
        return ret;
    }
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
};

/// \brief Solution of the Poisson problem on the unit cube
inline double Solution( const DROPS::Point3DCL& p, double=0.0){
    return 64*p[0]*p[1]*p[2]*(1.-p[0])*(1.-p[1])*(1.-p[2]);
}

namespace DROPS
{
/// \brief Solve a linear equation system by GMRES and Jacobi as preconditioner
template <typename Mat, typename Vec>
void Solve(const Mat &A, Vec &x, const Vec &b, __UNUSED__ const IdxDescCL& idx)
{
    // parameter for solver:
    const int    restart =  100;
    const int    maxiter = 1000;
    const double tol     = 1e-10;
    const bool   relative= false;

#ifndef _PAR
    // time measurement
    TimerCL timer;

    // type of preconditioner and solver
    typedef JACPcCL                 PrecondT;
    typedef GMResSolverCL<PrecondT> SolverT;

    // preconditioner and solver
    PrecondT pc;
    SolverT solver( pc, restart, maxiter, tol, relative);
#else
    // time measurement
    ParTimerCL timer;

    // type of preconditioner and solver
    typedef ParJac0CL                     PrecondT;
    typedef ParPreGMResSolverCL<PrecondT> SolverT;

    // preconditioner and solver
    PrecondT pc( idx);
    SolverT solver( restart, maxiter, tol, idx, pc, relative);
#endif

    std::cout << " o Solving system with Jacobi-GMRes:\n"
              << "   - maxiter   " << maxiter << '\n'
              << "   - tolerance " << tol     << '\n'
              << "   - restart   " << restart << std::endl;

    // Solve the linear equation system
    timer.Reset();
    solver.Solve( A, x, b);
    timer.Stop();

    double realresid;
#ifndef _PAR
    realresid= norm( VectorCL(A*x-b));
#else
    realresid= idx.GetEx().Norm( VectorCL(A*x-b),false);
#endif

    std::cout << " o Solved system with:\n"
              << "   - time          " << timer.GetTime()   << " s\n"
              << "   - iterations    " << solver.GetIter()  << '\n'
              << "   - residuum      " << solver.GetResid() << '\n'
              << "   - real residuum " << realresid         << std::endl;
}

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
    Solve( Poisson.A.Data.GetFinest(), Poisson.x.Data, Poisson.b.Data, Poisson.idx.GetFinest());


    // check the result
    // -------------------------------------------------------------------------
    std::cout << line << "Check result against known solution ...\n";

    timer.Reset();
    Poisson.CheckSolution( Poisson.x, Solution);
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;
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
        if (argc!=5){
            std::cout << "Usage " << argv[0] << " <ref x> <ref y> <ref z> <refinement>\n"
                      << " with\n"
                      << " o ref x: spatial resulotion in x direction\n"
                      << " o ref y: spatial resulotion in y direction\n"
                      << " o ref z: spatial resulotion in z direction\n"
                      << " o refinement: number of regular refinements of each tetrahedron\n\n"
                      << "Ending program" << std::endl;
            return 0;
        }

        // time measurement
#ifndef _PAR
        DROPS::TimerCL timer;
#else
        DROPS::ParTimerCL timer;
#endif

        // parameter for geometry
        const int refX  = atoi( argv[1]),
                  refY  = atoi( argv[2]),
                  refZ  = atoi( argv[3]),
                  refAll= atoi( argv[4]);
        DROPS::Point3DCL orig, e1, e2, e3;  // origin and orientation of unit cube
        e1[0]= e2[1]= e3[2]= 1.0;

        // set up data structure to represent a poisson problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Poisson problem ...\n";
        timer.Reset();

        // create builder for geometry
        DROPS::MGBuilderCL * mgb;
#ifdef _PAR
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
        if ( !DROPS::ProcCL::IamMaster())
            mgb = new DROPS::EmptyBrickBuilderCL( orig, e1, e2, e3);
#endif
        IF_MASTER
            mgb = new DROPS::BrickBuilderCL( orig, e1, e2, e3, refX, refY, refZ);

        // boundary conditions
        DROPS::BndCondT bndcond[6] = { DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC,
                                       DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC };
        // boundary function
        DROPS::BndDataCL<>::bnd_val_fun bndfunc[6] = { DROPS::Zero, DROPS::Zero, DROPS::Zero,
                                               DROPS::Zero, DROPS::Zero, DROPS::Zero};
        // boundary data ( = condition & function)
        DROPS::PoissonBndDataCL bdata( 6, bndcond, bndfunc);

        // Setup the problem
        DROPS::PoissonP2CL<PoissonCoeffCL> prob( *mgb, PoissonCoeffCL(), bdata);
        DROPS::MultiGridCL& mg= prob.GetMG();

#ifdef _PAR
        // Set parallel data structures
        pmg.AttachTo( mg);                                  // handling of parallel multigrid
        DROPS::LoadBalHandlerCL lb( mg);                    // loadbalancing
        lb.DoInitDistribution( DROPS::ProcCL::Master());    // distribute initial grid
        lb.SetStrategy( DROPS::Recursive);                  // best distribution of data
#endif

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        // Refine the grid
        // ---------------------------------------------------------------------
        std::cout << "Refine the grid " << refAll << " times regulary ...\n";
        timer.Reset();

        // Create new tetrahedra
        for ( int ref=1; ref<=refAll; ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( mg);
            mg.Refine();
        }

        // do loadbalancing
#ifdef _PAR
        lb.DoMigration();
#endif
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s\n"
                  << " o distribution of elements" << '\n';
        mg.SizeInfo(cout);

        // Solve the problem
        // ---------------------------------------------------------------------
        DROPS::Strategy( prob);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
