#include "geom/builder.h"

DROPS::Uint Rule(DROPS::Uint r)
{
    return r < 64 ? r : 127;
}

int TestReMark()
{
    int ret= 0;
    for (DROPS::Uint i= 0; i<=64; ++i) {
        for (DROPS::Uint j= 0; j<=64; ++j) {
            DROPS::TetraBuilderCL tet( Rule( i));
            DROPS::MultiGridCL mg( tet);
            tet.BogoReMark( mg, Rule( j));
            if ( mg.GetTetrasBegin( 0)->GetRefRule() != Rule( j)) {
                std::cout << "Aerger mit Regel " << Rule( i) << std::endl;
                ++ret;
            }
        }
    }
    return ret;
}

int main ()
{
  try {
    int ret= 0;
    std::cout << "Testing building single rules." << std::endl;
    for (DROPS::Uint i= 0; i<64; ++i) {
        DROPS::TetraBuilderCL brick( i);
        DROPS::MultiGridCL mg( brick);
        if ( mg.GetTetrasBegin( 0)->GetRefRule() != i) {
            std::cout << "Aerger mit Regel " << i << std::endl;
            ++ret;
        }
    }

    std::cout << "Testing remarking between different rules." << std::endl;
    ret+= TestReMark();
    return ret;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
