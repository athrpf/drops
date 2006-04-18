#include "num/spmat.h"
#include <iostream>


int Test()
{
    std::cout << "Blockmatrix 2x2 & 2x1 \\ 1x2 & 0s:\n" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    std::cerr << "A:\n" << A << '\n';
    DROPS::MatrixCL B;
    DROPS::MatrixBuilderCL BB(&B, 1, 2);
    BB( 0, 0)= 1.; BB( 0, 1)= 1.;
    BB.Build();
    std::cerr << "B:\n" << B << '\n';
    DROPS::BlockMatrixCL S( &A, DROPS::MUL, &B, DROPS::TRANSP_MUL, &B, DROPS::MUL, 0, DROPS::MUL);
    DROPS::VectorCL b( 3);
    b[0]= 1.; b[1]= 2.; b[2]= 3.;

    std::cerr << "nr0: "   << S.num_rows( 0) << "\tnr1: " << S.num_rows( 1) << "\tnr: " << S.num_rows()
	      << "\nnc0: " << S.num_cols( 0) << "\tnc1: " << S.num_cols( 1) << "\tnc: " << S.num_cols()
	      << "\nop0: " << S.GetOperation( 0) << '\t' << S.GetBlock( 0) << '\n'
	      << "op1: " << S.GetOperation( 1) << '\t' << S.GetBlock( 1) << '\n'
	      << "op2: " << S.GetOperation( 2) << '\t' << S.GetBlock( 2) << '\n'
	      << "op3: " << S.GetOperation( 3) << '\t' << S.GetBlock( 3) << '\n';
    std::cerr << "A:" << &A << '\t' << A.num_rows() << "x" << A.num_cols() << '\n'
	      << "B:" << &B << '\t' << B.num_rows() << "x" << B.num_cols() << '\n';
    DROPS::VectorCL x( S*b);
    std::cerr << "x0: " << x[0] << "\tx1: " << x[1] << "\tx2: " << x[2] << '\n';
    DROPS::VectorCL y( 3);
    y= transp_mul( S, b);
    std::cerr << "y0: " << y[0] << "\ty1: " << y[1] << "\ty2: " << y[2] << '\n';
    return 0;
}


int
main(int, char**)
{
  try {
    return Test();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
