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
    //add_col_to_vec( m, 1., b, 10);
    //b/= norm( b);
    b[0]= 1.;

    VectorCL ker( std::sqrt( m2.GetDiag()));  //( 1./std::sqrt(m.num_rows()), m.num_rows());
    ker/= norm( ker);
    std::cerr << "Kernel: rel. residual norm: " << dot( M*ker, ker)/norm_sq( ker) << " norm ker: " << norm(ker) << '\n';
    double alpha= std::sqrt( 10.*BBTDiag( m).min());
    VectorCL scaledKer( alpha * ker); // insertion leads to AA^T + alpha^2*ee^T
    m.insert_col( m.num_cols(), scaledKer);
    std::cerr << "Kernel: after stab: rel. residual norm: " << dot( M*ker, ker)/norm_sq( ker) << " alpha: " << alpha << '\n';

    //VectorCL D( 1.0/BBTDiag( m));
    //typedef DiagPcCL SPcT;
    //SPcT spc( D);
    typedef NEGSPcCL SPcT;
    SPcT spc( true);
    //typedef DummyPcCL SPcT;
    //SPcT spc;
    typedef SSORPcCL SPc2T;
    typedef GSPcCL SPc3T;
    SPc2T spc2;
    SPc2T spc3;
    //PCGSolverCL<SPcT> solver( spc, 200, 1e-6, true);
    PCGNESolverCL<SPcT> solver( spc, 200, 1e-6, true);
    //GCRSolverCL<SPcT> solver( spc, 500, 500, 1e-6, true, &std::cerr);
  
    VectorCL x( M.num_cols());
    TimerCL t;
    solver.Solve( m, x, b);
    //solver.Solve( M, x, b);
    t.Stop();
     std::cout << "Time: " << t.GetTime() << "Iter: " << solver.GetIter()
               << " Norm from solver: " << solver.GetResid()
               << " Norm x: " << norm( x) << " Residual: " << norm( M*x - b)
               << " Norm b: " << norm( b)
               << " ker^T*x " << dot( ker, x)/norm( ker)
               << " orth-part of x: " << norm( x - dot( ker, x)*ker) << '\n';

  } catch (DROPSErrCL d) {
    d.what( std::cerr);
  }

    return 0;
}
