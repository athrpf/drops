#include "num/solver.h"
#include "num/spmat.h"
#include "misc/utils.h"

#include <iostream>
#include <fstream>

using namespace DROPS;

int
main()
{
    std::ifstream s0( "Mvel.txt");
    MatrixCL M;
    in( s0, M);
    std::cerr << "supnorm( M): " << supnorm( M) << '\n';

    std::ifstream s1( "Mvel_old.txt");
    MatrixCL M1;
    in( s1, M1);
    M1*= 1107.;
    std::cerr << "supnorm( M1): " << supnorm( M1) << '\n';

    MatrixCL D;
    D.LinComb( 1., M, -1., M1);
    std::cerr << "supnorm( D): " << supnorm( D) << '\n';

    MatrixCL I( std::valarray<double>( 1.0, M.num_rows()));

    MatrixCL N;
    N.LinComb(1e-2*M.val().apply( std::fabs).max(), I, 1.0, M);
//    MatrixCL N( M1);

    // GMResSolverCL<SSORPcCL> ExactAsolver( SSORPcCL( 1.0), 500, 1000, 1e-6, /*relative=*/ true);
    PCGSolverCL<DummyPcCL> ExactAsolver( DummyPcCL(), 1000, 1e-6, /*relative=*/ true);

    VectorCL v( M.num_cols());
    VectorCL r( 1.0, M.num_rows());
    VectorCL r2( r/norm( r));
    ExactAsolver.Solve( N, v, r2);
    std::cerr << "\niterations: " << ExactAsolver.GetIter()
              << "\tresidual: " << ExactAsolver.GetResid()
              << '\n' << norm( v)
              << '\n' << norm( N*v - r2)
              << '\n';
    return 0;
}
