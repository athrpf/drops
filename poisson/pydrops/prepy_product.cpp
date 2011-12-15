/// \file  prepy_product.cpp
/// \brief Setup problem for py_product function
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

double GetProductF1(const DROPS::Point3DCL& p, double t){return PyC.GetProductF1(p,t);}
double GetProductF2(const DROPS::Point3DCL& p, double t){return PyC.GetProductF2(p,t);}
//for testing
double One(const DROPS::Point3DCL&, double) { return 1.0; }

namespace DROPS
{

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
        Poisson.SetupInstatSystem( Poisson.A, Poisson.M, 0.0, 0 );
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;

}

} // end of namespace DROPS

template<class CoeffCL>
void prepy_product(DROPS::ParamCL& P, DROPS::PoissonP1CL<CoeffCL>& Poisson)
{
  Strategy(Poisson, P);      
}

int main(int argc, char** argv)
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
            std::cout << "Using default parameter file: preproduct.json\n";
            param.open( "preproduct.json");
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

        // Setup the problem
        DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, DROPS::PoissonCoeffCL<DROPS::ParamCL>(P), *bdata);    
      
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
        
        prepy_product(P, prob);                  //Setup System 

        DROPS::instat_scalar_fun_ptr f1;
        DROPS::instat_scalar_fun_ptr f2;
        DROPS::instat_scalar_fun_ptr test;
        test = One;
        
        typedef PdeFunction* PdeFunPtr;
        PdeFunPtr PDEf1(new TestPdeFunction(P,test));
        PdeFunPtr PDEf2(new TestPdeFunction(P,test));
        PyC.Init(P);                             //Inital Pyc
        PyC.SetProductFun(PDEf1, PDEf2);         //initalize two functions
        f1 = GetProductF1;                       //get function 1
        f2 = GetProductF2;                       //get function 2
    
        double result = DROPS::Py_product( prob.GetMG(), prob.idx, prob.A, prob.M, f1, f2, 0., false);
        std::cout<<"The result of py_product is"<<result<<std::endl;

        delete[] PDEf1;
        delete[] PDEf2;
        return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


