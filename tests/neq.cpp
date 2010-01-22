/// \file neq.cpps
/// \brief  Used to test/implement CG for the normal equations and the NEGSPcCL.
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/utils.h"
#include "num/spmat.h"
#include "num/solver.h"

using namespace DROPS;

///\brief Used to test/implement CG for the normal equations and the NEGSPcCL.
int main ()
{
  try {

    MatrixCL m, m2, m3;
    std::fstream in( "mker.txt"); // discretization of divergence-operator
    in >> m;
    std::fstream in2( "m2.txt");  // pressure mass-matrix
    in2 >> m2;
    std::fstream in3( "m3.txt");  // velocity mass-matrix
    in3 >> m3;
    ScaleRows( m, VectorCL( 1./std::sqrt( m2.GetDiag())));
    //ScaleCols( m, VectorCL( 1./std::sqrt( m3.GetDiag())));

    CompositeMatrixCL M( &m, TRANSP_MUL, &m, MUL);
    VectorCL b( M.num_rows());
    //add_col_to_vec( m, 1., b, 10); // Test with a consistent rhs.
    //b/= norm( b);
    b[0]= 1.; // This gives an inconsistent system of linear equations.

    VectorCL ker0( std::sqrt( m2.GetDiag()));  // Kernel of the scaled matrix.
    ker0/= std::sqrt( dot( ker0, ker0));
    typedef NEGSPcCL SPcT;
    SPcT spc( true);
    VectorCL ker( spc.mul( m, ker0));
    ker/= std::sqrt( dot( ker, ker0)); // ker is now an spc-unit-norm vector.
    std::cout << "Kernel: rel. residual norm: " << dot( M*ker, ker)/norm_sq( ker) << " norm ker: " << norm(ker) << '\n';
    double alpha= std::sqrt( 1.); // std::sqrt( BBTDiag( m).min());
    VectorCL scaledKer( alpha * ker); // insertion into A leads to AA^T + alpha^2*ee^T.
    m.insert_col( m.num_cols(), scaledKer);
    std::cout << "Kernel: after stab: rel. residual norm: " << dot( M*ker, ker)/norm_sq( ker) << " alpha: " << alpha << '\n';

    //VectorCL D( 1.0/BBTDiag( m));
    //typedef DiagPcCL SPcT;
    //SPcT spc( D);
    //typedef DummyPcCL SPcT;
    //SPcT spc;
    typedef SSORPcCL SPc2T;
    typedef GSPcCL SPc3T;
    SPc2T spc2;
    SPc2T spc3;
    //PCGSolverCL<SPcT> solver( spc, 200, 1e-6, true);
    PCGNESolverCL<SPcT> solver( spc, 400, 1e-6, true);
    //GCRSolverCL<SPcT> solver( spc, 500, 500, 1e-6, true, &std::cout);

    VectorCL x( M.num_cols());
    TimerCL t;
    solver.Solve( m, x, b);
    //solver.Solve( M, x, b);
    t.Stop();
    std::cout << "Time: " << t.GetTime() << "Iter: " << solver.GetIter()
              << " Norm from solver: " << solver.GetResid()
              << " Norm x: " << norm( x) << " Residual: " << norm( M*x - b)
              << " Norm b: " << norm( b)
              << " ker^T*x " << dot( ker0, x)/norm( ker0) //<< dot( ker, x)/norm( ker)
              << " orth-part of x: " << norm( x - dot( ker0, x)*ker0) << '\n';
    std::cout << "dim: " << m.num_rows() << " x " << m.num_cols() << '\n';

  } catch (DROPSErrCL d) {
    d.what( std::cout);
  }

    return 0;
}
