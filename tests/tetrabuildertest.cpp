#include "geom/builder.h"

int main ()
{
  try {
    int ret= 0;
    for (DROPS::Uint i= 0; i<64; ++i) {
        DROPS::TetraBuilderCL brick( i);
        DROPS::MultiGridCL mg( brick);
        if ( mg.GetTetrasBegin( 0)->GetRefRule() != i) {
            std::cerr << "Aerger mit Regel " << i << std::endl;
            ++ret;
        }
    }

    return ret;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
