#include "num/solver.h"
#include <iostream>


int Test()
{
    std::cout << "BiCGStab 2x2:\n" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    DROPS::VectorCL b( 1., 2);
    b[0]= 1.; b[1]= 2.;
    DROPS::VectorCL x( 0., 2);

    std::cout << "A\n" << A << "b\n" << b << std::endl;
    int mi= 10;
    double tol= 1e-10;
    DROPS::BICGSTAB( A, x, b, DROPS::DummyPcCL(), mi, tol);
    std::cout << x << DROPS::VectorCL( A*x - b) << '\n' << mi << '\n' << tol << std::endl;
    return 0;
}

int Test2()
{
    std::cout << "BiCGStab 4x4:\n" << std::endl;
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
    DROPS::VectorCL x( 0., 4);

    std::cout << "A\n" << A << "b\n" << b << std::endl;
    int mi= 10;
    double tol= 1e-10;
    DROPS::BICGSTAB( A, x, b, DROPS::DummyPcCL(), mi, tol);
    std::cout << x << DROPS::VectorCL( A*x - b) << '\n' << mi << '\n' << tol << std::endl;
    return 0;
}

int main (int, char**)
{
  try {
    return Test() + Test2();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
