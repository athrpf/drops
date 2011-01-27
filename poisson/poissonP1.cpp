/// \file poissonP1.cpp
/// \brief Solver for Poisson problem with P1 functions
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


using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

/// \brief Coefficients of the Poisson problem
/** The coefficients of the Poisson problem are:
    \f$ - \alpha \cdot \Delta u + Vel.(\nabla u) +q \cdot u = f \f$
*/

/// \brief stationary case + error estimator
/*
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

    /// \brief Solution
    static double Solution( const DROPS::Point3DCL& p, double=0.0){
        return 64.*p[0]*p[1]*p[2]*(1.-p[0])*(1.-p[1])*(1.-p[2]);
    }
};
*/

// du/dt - nu*laplace u + Vel grad u + q*u = f

/// \ instationary case
/*
template<class ParamsT>
class PoissonCoeffCL
{
  public:
    PoissonCoeffCL(ParamsT&) {}
    static double q(const DROPS::Point3DCL&, double) { return 0.0; }
    static double alpha(const DROPS::Point3DCL&, double)
      { return 1; }
    static double f(const DROPS::Point3DCL& p, double t)
      { return (-2.0*std::exp(t)*std::exp(p[0]+p[1]+p[2])); }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL&, double)
      { return DROPS::Point3DCL(0.); } // no convection

    static double Solution(const DROPS::Point3DCL& p, double t)
    { return (std::exp(t)*std::exp(p[0]+p[1]+p[2])); }
};
*/
/// \brief instationary film heat transfer case
// du/dt - nu*laplace u + Vel grad u + q*u = f

template<class ParamsT>
class PoissonCoeffCL
{
  private:
    static ParamsT C_;
    static double dx_, dy_;

  public:
	PoissonCoeffCL( ParamsT& params) {
		C_ = params;
		int nx_,ny_,nz_;
		double dz_;
		std::string mesh( C_.dmc_MeshFile), delim("x@");
		size_t idx_;
		while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
		mesh[idx_]= ' ';
		std::istringstream brick_info( mesh);
		brick_info >> dx_ >> dy_ >> dz_ >> nx_ >> ny_ >> nz_;
	}
    static double q(const DROPS::Point3DCL&, double) { return 0.0; }
    static double alpha(const DROPS::Point3DCL&, double)
      { return 1; }
    static double f(const DROPS::Point3DCL& p, double t)
    {
        const double u= C_.exp_Rho*9.81*dy_*dy_/2/C_.exp_Mu*1e-3;
        return std::cos((p[0] + t*u)/dx_*2*M_PI);
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL& p, double)
    {
        DROPS::Point3DCL ret;
        const double d= p[1]/dy_,
            u= C_.exp_Rho*9.81*dy_*dy_/2/C_.exp_Mu*1e-3;
        ret[0]= u*(2-d)*d; // Nusselt
        return ret;
    }
    static double Solution(const DROPS::Point3DCL& , double )
    { return (0.0); }
};

template<class ParamsT>
ParamsT PoissonCoeffCL<ParamsT>::C_;

template<class ParamsT>
double PoissonCoeffCL<ParamsT>::dx_;

template<class ParamsT>
double PoissonCoeffCL<ParamsT>::dy_;


/*
/// \brief drops.cpp
// -laplace u + q*u = f
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
    static double f(const DROPS::Point3DCL& p, double= 0.0)
        { return 128.0*( p[0]*p[1]*(1-p[0])*(1-p[1]) + p[0]*p[2]*(1-p[0])*(1-p[2])
                + p[1]*p[2]*(1-p[1])*(1-p[2]) ); }
    /// \brief Diffusion
    static double alpha(const DROPS::Point3DCL&, double=0.0){
        return 1.;
    }
    /// \brief Solution
    static double Solution( const DROPS::Point3DCL& p, double=0.0){
        return 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
    }
};
*/

namespace DROPS
{

typedef ParamPoissonProblemCL Params;
Params C;

template<class CoeffCL, class SolverT, class ParamsT>
void SolveStatProblem( PoissonP1CL<CoeffCL>& Poisson, SolverT& solver, ParamsT& param)
{
    // time measurements
#ifndef _PAR
    TimerCL timer;
    const bool doErrorEstimate= param.err_DoErrorEstimate;
#else
    const bool doErrorEstimate= false;
    if (param.err_DoErrorEstimate)
        std::cout << "Skipping Error-Estimation ..." << std::endl;
    ParTimerCL timer;
#endif

	if ( !doErrorEstimate) {
        Poisson.SetupSystem( Poisson.A, Poisson.b);
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
        if (C.pos_SolutionIsKnown) {
        	std::cout << line << "Check result against known solution ...\n";
        	Poisson.CheckSolution( Poisson.x, CoeffCL::Solution);
        }
	}
	else {
	    MultiGridCL& MG= Poisson.GetMG();
	    const typename PoissonP1CL<CoeffCL>::BndDataCL& BndData= Poisson.GetBndData();

	    MLIdxDescCL  loc_idx;
	    VecDescCL  loc_x;

	    MLIdxDescCL* new_idx= &Poisson.idx;
	    MLIdxDescCL* old_idx= &loc_idx;

	    VecDescCL* new_x= &Poisson.x;
	    VecDescCL* old_x= &loc_x;

        DoerflerMarkCL<typename PoissonP1CL<CoeffCL>::est_fun, typename PoissonP1CL<CoeffCL>::base_>
            Estimator( param.err_RelReduction, param.err_MinRatio, param.err_Threshold, param.err_Meas, param.err_DoMark,
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
            if (C.pos_SolutionIsKnown) {
            	std::cout << line << "Check result against known solution ...\n";
            	Poisson.CheckSolution( *new_x, CoeffCL::Solution);
            }
            new_marks = Estimator.Estimate( typename PoissonP1CL<CoeffCL>::const_DiscSolCL( new_x, &BndData, &MG) );

            std::swap( old_x, new_x);
            std::swap( old_idx, new_idx);
        } while ( new_marks && step++ < param.err_NumRef);
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

    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    Poisson.idx.SetFE( P1_FE);                                  // set quadratic finite elements
    if ( PoissonSolverFactoryHelperCL<Params>().MGUsed(C))
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
    if (C.tm_NumSteps != 0)
    	Poisson.SetupInstatSystem( Poisson.A, Poisson.M, 0.0);
    timer.Stop();
    std::cout << " o time " << timer.GetTime() << " s" << std::endl;


    // solve the linear equation system
    // -------------------------------------------------------------------------
    std::cout << line << "Solve the linear equation system ...\n";
    
    // type of preconditioner and solver
    PoissonSolverFactoryCL< Params> factory( C, Poisson.idx);
    PoissonSolverBaseCL* solver = factory.CreatePoissonSolver();

    if ( factory.GetProlongation() != 0)
        SetupP1ProlongationMatrix( mg, *(factory.GetProlongation()), &Poisson.idx, &Poisson.idx);

    // Solve the linear equation system
    Poisson.Init( Poisson.x, CoeffCL::Solution, 0.0);
    InstatPoissonThetaSchemeCL<PoissonP1CL<CoeffCL>, PoissonSolverBaseCL>
        ThetaScheme( Poisson, *solver, C.tm_Theta, C.tm_Convection);
    ThetaScheme.SetTimeStep(C.tm_StepSize, C.tm_Nu);

    if (C.tm_NumSteps == 0) {
    	SolveStatProblem( Poisson, *solver, C);
    }

    Ensight6OutCL  ens(C.ens_EnsCase+".case", C.tm_NumSteps+1, C.ens_Binary, C.ens_MasterOut);
    if ( C.ens_EnsightOut){
        const std::string filename= C.ens_EnsDir + "/" + C.ens_EnsCase;
        ens.Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(), C.ens_GeomName,       filename + ".geo"));
        ens.Register( make_Ensight6Scalar( Poisson.GetSolution(), "Temperatur", filename + ".tp", true));
        ens.Write();
    }

    // writer for vtk-format
    VTKOutCL vtkwriter(mg, "DROPS data", C.tm_NumSteps/C.vtk_VTKOut+1 , std::string(C.vtk_VTKDir + "/" + C.vtk_VTKName), C.vtk_Binary);    

    if (C.vtk_VTKOut){
        vtkwriter.Register( make_VTKScalar( Poisson.GetSolution(), "velocity") );
        vtkwriter.Write( Poisson.x.t);
    }

    for ( int step = 1; step <= C.tm_NumSteps; ++step) {
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
        if ( C.ens_EnsightOut && step%C.ens_EnsightOut==0)
            ens.Write( step*C.tm_StepSize);
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
    // time measurements
#ifndef _PAR
        DROPS::TimerCL timer;
#else
        DROPS::ParTimerCL timer;
#endif

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
    	param >> DROPS::C;
    	param.close();
    	std::cout << DROPS::C << std::endl;

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

        // DROPS::C.dmc_BoundaryType = 4: solving head transport problem (ipfilm.cpp)
        DROPS::CreateGeomPoisson (mg, bdata, DROPS::C.dmc_MeshFile, DROPS::C.dmc_GeomType, DROPS::C.dmc_BoundaryType, DROPS::C.dmc_BoundaryFncs, serfile, r);

        // Setup the problem
        DROPS::PoissonP1CL<PoissonCoeffCL<DROPS::Params> > prob( *mg, PoissonCoeffCL<DROPS::Params>(DROPS::C), *bdata);
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
        std::cout << "Refine the grid " << DROPS::C.dmc_InitialCond << " times regulary ...\n";
        timer.Reset();

        // Create new tetrahedra
        for ( int ref=1; ref<=DROPS::C.dmc_InitialCond; ++ref){
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
