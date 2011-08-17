/// \file mass.cpp
/// \brief check whether given mass matrices are pos. def.
/// \author LNM RWTH Aachen: Joerg Grande, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
    std::cout << "supnorm( M): " << supnorm( M) << '\n';

    std::ifstream s1( "Mvel_old.txt");
    MatrixCL M1;
    in( s1, M1);
    M1*= 1107.;
    std::cout << "supnorm( M1): " << supnorm( M1) << '\n';

    MatrixCL D;
    D.LinComb( 1., M, -1., M1);
    std::cout << "supnorm( D): " << supnorm( D) << '\n';

    MatrixCL I( std::valarray<double>( 1.0, M.num_rows()));

    std::valarray<double> Mval( M.raw_val(), M.num_nonzeros());
    MatrixCL N;
    N.LinComb(1e-2*Mval.apply( std::abs).max(), I, 1.0, M);
//    MatrixCL N( M1);

    // SSORPcCL ssorpc;
    // GMResSolverCL<SSORPcCL> ExactAsolver( ssorpc, 500, 1000, 1e-6, /*relative=*/ true);
    DummyPcCL dummypc;
    PCGSolverCL<DummyPcCL> ExactAsolver( dummypc, 1000, 1e-6, /*relative=*/ true);

    VectorCL v( M.num_cols());
    VectorCL r( 1.0, M.num_rows());
    VectorCL r2( r/norm( r));
    ExactAsolver.Solve( N, v, r2);
    std::cout << "\niterations: " << ExactAsolver.GetIter()
              << "\tresidual: " << ExactAsolver.GetResid()
              << '\n' << norm( v)
              << '\n' << norm( N*v - r2)
              << '\n';
    return 0;
}
