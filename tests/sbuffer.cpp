#include "misc/container.h"
#include <iostream>


int TestPushBack()
{
    DROPS::SBufferCL<int, 3> b;
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    b[0]= 0; b[1]= 1; b[2]= 2;
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    b.push_back( 3);
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    std::cout << b[-1] << '\t' << b[-2] << '\t' << b[-3] << '\t' << std::endl;
    std::cout << b[3] << '\t' << b[4] << '\t' << b[5] << '\t' << std::endl;
    std::cout << std::endl;
    return 0;
}

int TestRotate()
{
    DROPS::SBufferCL<int, 3> b;
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    b[0]= 0; b[1]= 1; b[2]= 2;
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    b.rotate( 2);
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    b.rotate( -1);
    b.rotate( -1);
    std::cout << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << std::endl;
    return 0;
}

int main (int, char**)
{
  try {
    return TestPushBack() + TestRotate();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
