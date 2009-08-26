#include "geom/multigrid.h"
#include "geom/builder.h"
#include "out/output.h"
#include <fstream>
#include <cmath>
#include <cstdio>

int main()
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    DROPS::MultiGridCL mg(brick);

    mg.SplitMultiBoundaryTetras();
    mg.MakeConsistentNumbering();
    std::cout << "Checking Sanity...\n" << std::flush;
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    std::cout << "...ok" << std::endl;

    std::ofstream geomview("splitboundary.geo");
    DROPS::GeomMGOutCL ggg( mg, -1, false, 3.0);
    geomview << ggg << std::flush;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
