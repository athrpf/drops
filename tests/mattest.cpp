#include "misc/container.h"
#include "num/solver.h"

const double dat[]= { 0., 1., 1., 2. };
const double dat2[]= { 1., 0. };

using namespace DROPS;


int main()
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
