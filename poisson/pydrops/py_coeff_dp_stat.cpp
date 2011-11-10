// STL
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>

// DROPS
#include "poisson/poisson.h"
#include "num/solver.h"

#include "functors.hpp"
#include "drops_utils.hpp"

#include "py_journalist.hpp"

#include "py_coeff_dp_stat.hpp"

namespace DROPS {
  void coeff_stat(Journalist& jnlst, MassTransferBrick& brick, PoissonCoeffCL& pcl, SolutionContainer* sol_container, double tol, int maxiter)
  {
    PoissonP1CL<PoissonCoeffCL>
      Poisson(brick.get_brick(), pcl, brick.get_bdata(), 0);

    MultiGridCL& MG= Poisson.GetMG();
    IdxDescCL& idx= Poisson.idx;
    VecDescCL& x= Poisson.x;
    VecDescCL& b= Poisson.b;
    MatDescCL& A= Poisson.A;
    MatDescCL& M= Poisson.M;

    VecDescCL cplA, dummy;

    idx.Set(1);
    // erzeuge Nummerierung zu diesem Index
    Poisson.CreateNumbering(MG.GetLastLevel(), &idx);

    // Vektoren mit Index idx
    x.SetIdx(&idx);
    b.SetIdx(&idx);
    cplA.SetIdx( &idx);
    dummy.SetIdx( &idx);
    // Matrizen mit Index idx (Zeilen und Spalten)
    A.SetIdx(&idx, &idx);
    M.SetIdx(&idx, &idx);

    jnlst << "Number of Unknowns: " << (int)x.Data.size() << "\n";
    jnlst << "Tolerance GMRES: " << tol <<"\t";
    jnlst << "max. Num. GMRES-Iterations: " << maxiter << "\n";

    // stationaerer Anteil
    Poisson.SetupInstatSystem(A, M, Poisson.t);
    if (false) {
      Poisson.SetupInstatRhs( cplA, dummy, 0, dummy, 0);
      //Poisson.SetupGradSrc(b, directproblem_solution, perturbation, NULL);
    } else {
      Poisson.SetupInstatRhs( cplA, dummy, 0, b, 0);
    }
    b.Data+= cplA.Data;

    //initialize solver
    SSORPcCL pc(1.0);
    int numiter;
    double resid;
    PCG_SsorCL solver(pc, maxiter, tol);
    //solve
    solver.Solve(A.Data, x.Data, b.Data);
    if (solver.GetResid()>tol) {
      std::stringstream ss;
      ss << "The solver tolerance could not be met.\n"
	 << "residual: " << solver.GetResid() << ", tol = " << tol
	 << "\nnumiter: " << solver.GetIter() << ", maxiter = " << maxiter;
      throw(DROPS::DROPSErrCL(ss.str().c_str()));
    }

    sol_container->set_solution(Poisson, 0.0);

    A.Reset();
    // M.Reset(); we still need this matrix...
    b.Reset();
  }
}
