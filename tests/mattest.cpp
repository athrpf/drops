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
