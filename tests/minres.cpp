#include "num/solver.h"
#include <iostream>

class TrivialPreCL
{
  public:
    TrivialPreCL() {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& K, Vec& x, const Vec& b) const {
        x= b; }
};

int TestLanczos()
{
    std::cout << "Lanczos:\n" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    DROPS::VectorCL r0( 1., 2);
    r0[0]= 1.; r0[1]= 2.;

    std::cout << "A\n" << A << "r0\n" << r0 << std::endl;
    DROPS::LanczosONBCL<DROPS::MatrixCL, DROPS::VectorCL> onb;
    onb.new_basis(A, r0);
    std::vector<DROPS::VectorCL> basis;
    do {
        basis.push_back( onb.q[0]);
        std::cout << onb.q[0] << std::endl;
        std::cout << "a " << onb.a0 << " b " << onb.b[0] << std::endl;
    } while (onb.next());
    basis.push_back( onb.q[0]);
    std::cout << onb.q[0] << std::endl;
    std::cout << "a " << onb.a0 << " b " << onb.b[0] << std::endl;
    std::cout << "lucky breakdown after " << basis.size() << " vectors."
              <<std::endl;
    for (unsigned int i= 0; i<basis.size(); ++i) {
        for (unsigned int j= 0; j<basis.size(); ++j) {
            std::cout << basis[i]*basis[j] << '\t';
        }
        std::cout << std::endl;
    }
    return 0;
}

int TestMinres()
{
    std::cout << "Minres 2x2:\n" << std::endl;
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
    MINRES( A, x, b, mi, tol);
    std::cout << x << A*x - b << '\n' << mi << '\n' << tol << std::endl;
    return 0;
}

int TestMinres2()
{
    std::cout << "Minres 4x4:\n" << std::endl;
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
    MINRES( A, x, b, mi, tol);
    std::cout << x << A*x - b << '\n' << mi << '\n' << tol << std::endl;
    return 0;
}

int TestPMinres()
{
    std::cout << "PMinres 4x4:\n" << std::endl;
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
    DROPS::LanczosONBCL<DROPS::MatrixCL, DROPS::VectorCL> q;
    q.new_basis( A, b);
    DROPS::PMResSolverCL<DROPS::LanczosONBCL<DROPS::MatrixCL, DROPS::VectorCL> > pmr( q, mi, tol);
    pmr.Solve( A, x, b);
    std::cout << x << A*x - b << '\n' << pmr.GetIter() << '\n' << pmr.GetResid() << std::endl;
    return 0;
}

int TestPMinres2()
{
    std::cout << "PMinres2 4x4:\n" << std::endl;
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
    DROPS::DummyPcCL pc;
    DROPS::PLanczosONBCL<DROPS::MatrixCL, DROPS::VectorCL, DROPS::DummyPcCL> q( pc);
    DROPS::PMResSolverCL<DROPS::PLanczosONBCL<DROPS::MatrixCL, DROPS::VectorCL, DROPS::DummyPcCL> > pmr( q, mi, tol);
    pmr.Solve( A, x, b);
    std::cout << x << A*x - b << '\n' << pmr.GetIter() << '\n' << pmr.GetResid() << std::endl;
    return 0;
}

int main (int, char**)
{
  try {
    return TestLanczos() + TestMinres() + TestMinres2()
           + TestPMinres() + TestPMinres2();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
