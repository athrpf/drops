/// \file blockmat.cpp
/// \brief tests implementation of BlockMatrixCL
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "num/spmat.h"
#include "num/spblockmat.h"
#include <iostream>


int Test()
{
    std::cout << "Blockmatrix 2x2 & 2x1 \\ 1x2 & 0s:\n" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    std::cout << "A:\n" << A << '\n';
    DROPS::MatrixCL B;
    DROPS::MatrixBuilderCL BB(&B, 1, 2);
    BB( 0, 0)= 1.; BB( 0, 1)= 1.;
    BB.Build();
    std::cout << "B:\n" << B << '\n';
    DROPS::BlockMatrixCL S( &A, DROPS::MUL, &B, DROPS::TRANSP_MUL, &B, DROPS::MUL, 0, DROPS::MUL);
    DROPS::VectorCL b( 3);
    b[0]= 1.; b[1]= 2.; b[2]= 3.;

    std::cout << "nr0: "   << S.num_rows( 0) << "\tnr1: " << S.num_rows( 1) << "\tnr: " << S.num_rows()
              << "\nnc0: " << S.num_cols( 0) << "\tnc1: " << S.num_cols( 1) << "\tnc: " << S.num_cols()
              << "\nop0: " << S.GetOperation( 0) << '\t' << S.GetBlock( 0) << '\n'
              << "op1: " << S.GetOperation( 1) << '\t' << S.GetBlock( 1) << '\n'
              << "op2: " << S.GetOperation( 2) << '\t' << S.GetBlock( 2) << '\n'
              << "op3: " << S.GetOperation( 3) << '\t' << S.GetBlock( 3) << '\n';
    std::cout << "A:" << &A << '\t' << A.num_rows() << "x" << A.num_cols() << '\n'
              << "B:" << &B << '\t' << B.num_rows() << "x" << B.num_cols() << '\n';
    DROPS::VectorCL x( S*b);
    std::cout << "x0: " << x[0] << "\tx1: " << x[1] << "\tx2: " << x[2] << '\n';
    DROPS::VectorCL y( 3);
    y= transp_mul( S, b);
    std::cout << "y0: " << y[0] << "\ty1: " << y[1] << "\ty2: " << y[2] << '\n';
    return 0;
}

int TestComposite()
{
    std::cout << "\n\nCompositeMatrixCL:\n" << std::endl;
    DROPS::MatrixCL A;
    DROPS::MatrixBuilderCL AB(&A, 2, 2);
    AB( 0, 1)= 1.; AB( 1, 0)= 1.;
    AB.Build();
    std::cout << "A:\n" << A << '\n';
    DROPS::MatrixCL B;
    DROPS::MatrixBuilderCL BB(&B, 1, 2);
    BB( 0, 0)= 1.; BB( 0, 1)= 1.;
    BB.Build();
    std::cout << "B:\n" << B << '\n';
    DROPS::CompositeMatrixCL S( &A, DROPS::TRANSP_MUL, &B, DROPS::MUL);
    DROPS::VectorCL b( 2);
    b[0]= 1.; b[1]= 2.;

    std::cout << "nr: "   << S.num_rows()
              << "\nnc: " << S.num_cols()
              << "\nni: " << S.intermediate_dim()
              << "\nop0: " << S.GetOperation( 0) << '\t' << S.GetBlock0() << '\n'
              << "op1: " << S.GetOperation( 1) << '\t' << S.GetBlock1() << '\n';
    std::cout << "A:" << A << '\n' << A.num_rows() << "x" << A.num_cols() << '\n'
              << "B:" << B << '\n' << B.num_rows() << "x" << B.num_cols() << '\n';
    DROPS::VectorCL x( S*b);
    std::cout << "x.size(): " << x.size() << " x0: " << x[0] << "\n";
    DROPS::VectorCL b2( 1);
    b2[0]= 3.;
    std::cout << "b2: " << b2 << '\n' << "transp_mul( B, b2): " << transp_mul( B, b2) << '\n';

    DROPS::VectorCL y( transp_mul( S, b2));
    std::cout << "y.size(): " << y.size()  << " y0: " << y[0] << "\ty1: " << y[1] << "\n";

    DROPS::CompositeMatrixCL St( S.GetTranspose());
    std::cout << "nr: "   << St.num_rows()
              << "\nnc: " << St.num_cols()
              << "\nni: " << St.intermediate_dim()
              << "\nop0: " << St.GetOperation( 0) << '\t' << St.GetBlock0() << '\n'
              << "op1: " << St.GetOperation( 1) << '\t' << St.GetBlock1() << '\n';
    std::cout << "A:" << &A << '\n' << A.num_rows() << "x" << A.num_cols() << '\n'
              << "B:" << &B << '\n' << B.num_rows() << "x" << B.num_cols() << '\n';

    typedef DROPS::CompositeMatrixBaseCL<DROPS::CompositeMatrixCL, DROPS::MatrixCL> CompositeMatrix3CL;
    CompositeMatrix3CL S2( &S, DROPS::MUL, &B, DROPS::TRANSP_MUL);
    std::cout << "nr: "   << S2.num_rows()
              << "\nnc: " << S2.num_cols()
              << "\nni: " << S2.intermediate_dim()
              << "\nop0: " << S2.GetOperation( 0) << '\t' << S2.GetBlock0() << '\n'
              << "op1: " << S2.GetOperation( 1) << '\t' << S2.GetBlock1() << '\n';
    std::cout << "S2*b: " << S2*b << '\n';

    typedef DROPS::CompositeMatrixBaseCL<DROPS::CompositeMatrixCL, DROPS::VectorAsDiagMatrixCL> CompositeMatrix4CL;
    DROPS::VectorCL v( 1);
    DROPS::VectorAsDiagMatrixCL Mv( &v);
    CompositeMatrix4CL S3( &S, DROPS::MUL, &Mv, DROPS::MUL);
    std::cout << "nr: "   << S3.num_rows()
              << "\nnc: " << S3.num_cols()
              << "\nni: " << S3.intermediate_dim()
              << "\nop0: " << S3.GetOperation( 0) << '\t' << S3.GetBlock0() << '\n'
              << "op1: " << S3.GetOperation( 1) << '\t' << S3.GetBlock1() << '\n';
    std::cout << "S2*b: " << S3*b << '\n';

    return 0;
}

int
main(int, char**)
{
  try {
    return Test() + TestComposite();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
