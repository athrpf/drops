/// \file mattest.cpp
/// \brief tests Gauss solver and QR decomposition
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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

#include "misc/container.h"
#include "num/solver.h"

const double dat[]= { 0., 1., 1., 2. };
const double dat2[]= { 1., 0. };

const double dat3[]= { -9./2., -4., 3., -3., 1., 2. };
const double dat4[]= { -7., -2., -9. };


using namespace DROPS;

int Test1()
{
    SMatrixCL<2,2> M(dat+0);
    SVectorCL<2> v(dat2+0);
    std::cout << M << std::endl << v << std::endl;
#if 0
    std::cout << v-2.0*v << std::endl;
    std::cout << (M*= 0.5) << std::endl;
    std::cout << M+M/4. << std::endl;
    std::cout << M*v << std::endl;
    std::cout << M*M << std::endl;
#endif
    gauss_pivot(M, v);
    std::cout << v << std::endl;
    return 1;
}

int Test2()
{
    QRDecompCL<2> qr ( dat + 0);
    SVectorCL<2> r( dat2 + 0);
    qr.Solve( r);
    std::cout << r << std::endl;

    std::cout << "Least squares: solution should be (0.4 0.4) with residual 11. Drops says:" << std::endl;
    QRDecompCL<3,2> qrls ( dat3 + 0);
    SVectorCL<3> rls( dat4 + 0);
    qrls.Solve( rls);
    std::cout << rls << std::endl;

    return 0;
}

int main()
{
    return Test1() + Test2();
}
