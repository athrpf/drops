#include "misc/container.h"
#include "out/output.h"
#include "geom/builder.h"
#include "geom/multigrid.h"
#include "geom/boundary.h"
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>



const char meshfile[]= {
    "(0 \"GAMBIT to Rampant Filtered File\")\n"
    "(1 \"Joerg: DROPS-Checking.\")\n"
    "(0 \"Dimension:\")\n"
    "(2 3)\n"
    "\n"
    "(0 \"Nodes:\")\n"
    "(10 (0 1 9 1))\n"
    "(10 (a 1 9 1) (\n"
    " 1.000000000e+00  0.000000000e+00  1.000000000e+00\n"
    " 1.000000000e+00  0.000000000e+00  0.000000000e+00\n"
    " 1.000000000e+00  1.000000000e+00  0.000000000e+00\n"
    " 1.000000000e+00  1.000000000e+00  1.000000000e+00\n"
    " 0.000000000e+00  1.000000000e+00  0.000000000e+00\n"
    " 0.000000000e+00  1.000000000e+00  1.000000000e+00\n"
    " 0.000000000e+00  0.000000000e+00  0.000000000e+00\n"
    " 0.000000000e+00  0.000000000e+00  1.000000000e+00\n"
    " 9.999599457e-01  5.000000000e-01  5.000000000e-01\n"
    "))\n"
    "\n"
    "(0 \"Faces:\")\n"
    "(13 (0 1 1e 0))\n"
    "(13 (2 1 4 3 3) (\n"
    "6 5 3 4 0\n"
    "6 3 4 8 0\n"
    "4 2 1 9 0\n"
    "4 3 2 a 0\n"
    "))\n"
    "(13 (3 5 8 24 3) (\n"
    "1 7 8 1 0\n"
    "8 7 5 3 0\n"
    "8 5 6 5 0\n"
    "1 2 7 b 0\n"
    "))\n"
    "(13 (4 9 a a 3) (\n"
    "3 5 7 2 0\n"
    "3 7 2 c 0\n"
    "))\n"
    "(13 (5 b c e 3) (\n"
    "1 8 6 6 0\n"
    "1 6 4 7 0\n"
    "))\n"
    "(13 (7 d 1e 2 3) (\n"
    "5 7 9 3 2\n"
    "7 8 9 3 1\n"
    "3 5 9 4 2\n"
    "6 5 9 5 4\n"
    "5 8 9 5 3\n"
    "6 8 9 6 5\n"
    "8 1 9 6 1\n"
    "6 1 9 7 6\n"
    "3 6 9 8 4\n"
    "6 4 9 8 7\n"
    "4 1 9 9 7\n"
    "3 4 9 a 8\n"
    "4 2 9 a 9\n"
    "2 1 9 b 9\n"
    "1 7 9 b 1\n"
    "2 7 9 c b\n"
    "7 3 9 c 2\n"
    "3 2 9 c a\n"
    "))\n"
    "\n"
    "(0 \"Cells:\")\n"
    "(12 (0 1 c 0))\n"
    "(12 (1 1 c 1 2))\n"
    "\n"
    "(0 \"Zones:\")\n"
    "(45 (1 fluid fluid)())\n"
    "(45 (2 wall wall.1)())\n"
    "(45 (3 outflow outflow.2)())\n"
    "(45 (4 velocity-inlet velocity_inlet.3)())\n"
    "(45 (5 fan fan.4)())\n"
    "(45 (7 interior default-interior)())\n"
};

/* Make member functions of ReadMeshBuilderCL public to test them.
int TestPrimitives()
{
    const char* s[]= {
        "\"Bla1\"", "\"\"", "\" \"", "\"\n\"", "  \n  \"blub\"", "\"ha\"llo\"",
        "\"bla bla\" zweiter Teil"
    };
    const char* s2[]= {
        "blub(12345", "(1\"(falle)\"(2(3(4", "(2(3 )4)5"
    };
    const char* s3[]= {
        " (1 1a )x", "()"
    };
    DROPS::ReadMeshBuilderCL b( std::cin, false);
    std::cout << "0: " << b.Symbolic( 0) << "\n-1u: " << b.Symbolic( -1)
              << std::endl;
    for (unsigned i=0; i<7; ++i) {
        std::cout << s[i] << std::endl;
        std::istringstream is( s[i]);
        std::string st;
        DROPS::MeshStringCL mst;
        is >> mst >> st;
        std::cout << st << std::endl;
        if (i==6) { st.clear(); is >> st; std::cout << st; }
        std::cout << std::endl;
    }
    std::cout << s[6] << std::endl;
    std::istringstream is( s[6]);
    std::string st;
    DROPS::MeshStringCL mst;
    is >> mst >> st;
    std::cout << st << std::endl;
    std::cout << "--------------------------" << std::endl;
    {
        std::cout << s2[0] << std::endl;
        std::istringstream is( s2[0]);
        DROPS::ReadMeshBuilderCL builder( is, false);
        builder.NextSection();
        std::istream_iterator<char> i1( is), i2;
        std::copy( i1, i2,
                   std::ostream_iterator<char>( std::cout));
        std::cout << std::endl;
    }
    {
        std::cout << s2[1] << std::endl;
        std::istringstream is( s2[1]);
        DROPS::ReadMeshBuilderCL builder( is, false);
        builder.NextSection();

        std::istream_iterator<char> i2;
//        std::copy( std::istream_iterator<char>( is), i2,
//                   std::ostream_iterator<char>( std::cout));
        is.get();
        builder.NextSection();
//        std::copy( std::istream_iterator<char>( is), i2,
//                   std::ostream_iterator<char>( std::cout));
        is.get();
        builder.NextSection();
        std::copy( std::istream_iterator<char>( is), i2,
                   std::ostream_iterator<char>( std::cout));
        std::cout << std::endl;
    }
    {
        std::cout << s2[2] << std::endl;
        std::istringstream is( s2[2]);
        DROPS::ReadMeshBuilderCL builder( is, false);
        builder.SkipSection();
        std::istream_iterator<char> i2;
        std::copy( std::istream_iterator<char>( is), i2,
                   std::ostream_iterator<char>( std::cout));
        std::cout << std::endl;
    }
    {
        std::cout << s3[0] << std::endl;
        std::istringstream is( s3[0]);
        DROPS::ReadMeshBuilderCL builder( is, false);
        builder.NextSection();
        std::vector<unsigned> hgh= builder.ReadHeaderInfoHex();
        std::istream_iterator<char> i1( is), i2;
        std::copy( i1, i2,
                   std::ostream_iterator<char>( std::cout));
        std::cout << std::endl;
        std::copy( hgh.begin(), hgh.end(),
                   std::ostream_iterator<unsigned>( std::cout));
        std::cout << std::endl;
    }
    return 0;
}
*/

int TestReading()
{
    std::cout << "---------------------------------------------------\n";
    std::istringstream is( meshfile);
//    std::copy( std::istream_iterator<char>( is), std::istream_iterator<char>(),
//               std::ostream_iterator<char>( std::cerr));
    DROPS::ReadMeshBuilderCL builder( is);
    DROPS::MultiGridCL mg( builder);
    std::cerr << DROPS::DumpMGCL( mg) << std::endl;
    std::cerr << DROPS::SanityMGOutCL( mg) << std::endl;
    
    return 0;
}

int Test()
{
//    std::cout << "---------------------------------------------------\n";
    DROPS::ReadMeshBuilderCL builder( std::cin);
    DROPS::MultiGridCL mg( builder);
//    std::cerr << DROPS::DumpMGCL( mg) << std::endl;
    std::cerr << DROPS::SanityMGOutCL( mg) << std::endl;
    std::cout << DROPS:: GeomMGOutCL( mg) << std::flush;   
    return 0;
}

int main (int, char**)
{
/*
    int i;
    std::cin >> std::hex >> i;
    std::cout << i << std::endl;
    std::cin >>  i;
    std::cout << i << std::endl;
*/    

  try {
    return Test();
//    return TestPrimitives() + TestReading();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
