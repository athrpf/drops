#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "geom/builder.h"
#include "num/discretize.h"
#include "num/fe.h"

using namespace DROPS;

int degx, degy;

double f(const Point3DCL& p, double)
{
    return  (degx==0 ? 1. : std::pow( p[0], degx))
           *(degy==0 ? 1. : std::pow( p[1], degy));
}

double exactint[21] = { // with maple
  0.50000000000000000000, 0.16666666666666666667, 0.083333333333333333333,

  0.050000000000000000000, 0.033333333333333333333, 0.023809523809523809524,

  0.16666666666666666667, 0.041666666666666666667, 0.016666666666666666667,

  0.0083333333333333333333, 0.0047619047619047619048, 0.083333333333333333333,

  0.016666666666666666667, 0.0055555555555555555556, 0.0023809523809523809524,

  0.050000000000000000000, 0.0083333333333333333333, 0.0023809523809523809524,

  0.033333333333333333333, 0.0047619047619047619048, 0.023809523809523809524
};

void TestExactness()
{
    DROPS::TetraBuilderCL tet( 0);
    DROPS::MultiGridCL mg( tet);
    BaryCoordCL b[3];
    b[0]= DROPS::std_basis<4>(1);
    b[1]= DROPS::std_basis<4>(2);
    b[2]= DROPS::std_basis<4>(3);
    TetraCL& s= *mg.GetAllTetraBegin();
//    s.DebugInfo( std::cerr);
    std::cout.precision(25);

    Quad5_2DCL<> q;
    size_t c= 0;
         for (degy= 0; degy <= 5; ++degy) {
            for (degx= 0; degx + degy <= 5; ++degx) {
                q.assign( s, b, f);
                std::cout <<  "\tdegy: " << degy << "\tdegx: " << degx
                          << "\t\tI-Q_h: " << exactint[c++] - q.quad( 1.)
                          << "\tIntegral: " << q.quad( 1.) <<  "             ";
//                 for (size_t i= 0; i < q.size(); ++i)
//                     std::cout << '\t' << q[i];
               std::cout << std::endl;
            }
        }
}


int main ()
{
  try {
    TestExactness();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
