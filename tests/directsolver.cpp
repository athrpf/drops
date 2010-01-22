/// \file directsolver.cpp
/// \brief tests implementation of direct solvers
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande; SC RWTH Aachen:

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

#include "num/directsolver.h"
#include <iostream>

int Test()
{
    std::cout << "DirectSolver 4x4:\n" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 4, 4);
    AB( 0, 0)= -249.;
    AB( 0, 1)= -453.;
    AB( 0, 2)= -397.;
    AB( 0, 3)= -52.;

    AB( 1, 0)= -453.;
    AB( 1, 1)= -731.;
    AB( 1, 2)= -601.;
    AB( 1, 3)= -78.;

    AB( 2, 0)= -397.;
    AB( 2, 1)= -601.;
    AB( 2, 2)= -648.;
    AB( 2, 3)= -91.;

    AB( 3, 0)= -52.;
    AB( 3, 1)= -78.;
    AB( 3, 2)= -91.;
    AB( 3, 3)= -13.;
    AB.Build();
    DROPS::VectorCL b( 0., 4);
    b[0]= -1277./2/1.;
    b[1]= -2015./2.;
    b[2]= -3907./4.;
    b[3]= -533./4.;
    DROPS::VectorCL x(0., 4);

    std::cout << "A\n" << A << "b\n" << b << std::endl;
    DROPS::DirectSymmSolverCL dsolver(A);
    dsolver.Solve(A,x,b);
    dsolver.Update(A);
    dsolver.Solve(A,x,b);
    std::cout << "x\n" << x << std::endl;
    DROPS::DirectNonSymmSolverCL dnsolver(A);
    dnsolver.Solve(A,x,b);
    std::cout << "x\n" << x << std::endl;
    return 0;
}


int main (int, char**)
{
  try {
    return Test();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
