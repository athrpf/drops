#include "geom/multigrid.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "geom/builder.h"
#include "num/stokessolver.h"
#include "stokes/instatstokes.h"
#include "stokes/integrTime.h"
#include <fstream>
#include <sstream>


struct StokesCL
{
    static double g_;
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]=    sin(p[0])*sin(p[1])*sin(p[2]);
        ret[1]=  - cos(p[0])*cos(p[1])*sin(p[2]);
        ret[2]= 2.*cos(p[0])*sin(p[1])*cos(p[2]);
        return ret/3.;
    }

    // Jacobi-matrix od exact solution
    static inline DROPS::SMatrixCL<3, 3> DLsgVel(const DROPS::Point3DCL& p)
    {
        DROPS::SMatrixCL<3, 3> ret;
        ret(0,0)= cos(p[0])*sin(p[1])*sin(p[2])/3.;
        ret(0,1)= sin(p[0])*cos(p[1])*sin(p[2])/3.;
        ret(0,2)= sin(p[0])*sin(p[1])*cos(p[2])/3.;

        ret(1,0)=   sin(p[0])*cos(p[1])*sin(p[2])/3.;
        ret(1,1)=   cos(p[0])*sin(p[1])*sin(p[2])/3.;
        ret(1,2)= - cos(p[0])*cos(p[1])*cos(p[2])/3.;

        ret(2,0)= -2.*sin(p[0])*sin(p[1])*cos(p[2])/3.;
        ret(2,1)=  2.*cos(p[0])*cos(p[1])*cos(p[2])/3.;
        ret(2,2)= -2.*cos(p[0])*sin(p[1])*sin(p[2])/3.;
        return ret;
    }

    static double LsgPr(const DROPS::Point3DCL& p)
    {
        return cos(p[0])*sin(p[1])*sin(p[2])
               -(sin( 1.) -2.*sin( 1.)*cos( 1.) + sin( 1.)*pow( cos( 1.), 2)); // (...)==0.1778213062
    }

    // q*u - nu*laplace u + Dp = f
    //                  -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return StokesCL::g_; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p)
        { 
            const double g= StokesCL::g_;
            DROPS::SVectorCL<3> ret;
            ret[0]= g/3.*sin(p[0])*sin(p[1])*sin(p[2]);
            ret[1]= -g/3.*cos(p[0])*cos(p[1])*sin(p[2]);
            ret[2]= cos(p[0])*sin(p[1])*cos(p[2])*(2./3.*g + 3.);
            return ret;
        }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };
    static StokesCoeffCL Coeff;
};

double StokesCL::g_;
StokesCL::StokesCoeffCL StokesCL::Coeff;


namespace DROPS
{

class ScaledMGPreCL
{
  private:
    MGDataCL& A_;
    Uint iter_;
    double s_;
    mutable SSORPcCL directpc_;
    mutable PCG_SsorCL solver_;
    Uint sm_; // how many smoothing steps?
    int lvl_; // how many levels? (-1=all)
    double omega_; // relaxation parameter for smoother
    mutable SSORsmoothCL smoother_;  // Symmetric-Gauss-Seidel with over-relaxation

  public:
    ScaledMGPreCL( MGDataCL& A, Uint iter, double s= 1.0)
        :A_( A), iter_( iter), s_( s), solver_( directpc_, 200, 1e-12), sm_( 1),
         lvl_( -1), omega_( 1.0), smoother_( omega_)
    {}

    template <class Mat, class Vec>
    inline void
    Apply( const Mat&, Vec& x, const Vec& r) const;
};


template <class Mat, class Vec>
inline void
ScaledMGPreCL::Apply( const Mat&, Vec& x, const Vec& r) const
{
    x= 0.0;
    const double oldres= norm(r - A_.back().A.Data*x);
    for (Uint i= 0; i < iter_; ++i) {
        MGM( A_.begin(), --A_.end(), x, r, smoother_, sm_, solver_, lvl_, -1);
        x*= s_;
    }
    const double res= norm(r - A_.back().A.Data*x);
    std::cerr << "ScaledMGPreCL: it: " << iter_ << "\treduction: " << (oldres==0.0 ? res : res/oldres) << '\n';
}

class DiagMatrixPCCL
{
  private:
    MatrixCL& M_;

  public:
    DiagMatrixPCCL( MatrixCL& M)
        :M_( M) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& , Vec& x, const Vec& b) const
    {
        for (Uint i= 0; i<M_.num_rows(); ++i) {
            x[i]= b[i]/M_.val( i); // M_ is a diagonal-matrix: exact inversion
        }
    }
};


template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
SchurAR(const Mat& A, const Mat& B, Vec& xu, Vec& xp, const Vec& f, const Vec& g,
    PC1& Apc, PC2& Spc,
    int& max_iter, double& tol)
{
    VectorCL ru= f - A*xu - transp_mul( B, xp);
    VectorCL w( f.size());
    VectorCL z( g.size());
    VectorCL a( f.size());
    VectorCL xuneu( f.size());
    VectorCL b( f.size());
    ApproximateSchurComplMatrixCL<PC1> asc( A, Apc, B);
    PCGSolverCL<PC2> pcgsolver( Spc, 100, 0.6);
    const double resid0= std::sqrt( norm_sq( ru) + norm_sq( g - B*xu));
    double resid= 0.0;
    std::cerr << "residual (2-norm): " << resid0 << '\n';

    for (int k= 1; k<=max_iter; ++k) {
        w= 0.0;
        Apc.Apply( A, w, ru);
        w+= xu;
        z= 0.0;
//        std::cerr << "B*w : " << norm( B*w) << "\tprojektion auf 1: "
//                  << (B*w)*VectorCL( 1.0/std::sqrt( (double)g.size()), g.size())
//                  << "\n";
        pcgsolver.Solve( asc, z, VectorCL( B*w - g));
        std::cerr << "pcgsolver: iterations: " << pcgsolver.GetIter() 
                  << "\tresid: " << pcgsolver.GetResid() << '\n';

        b= transp_mul( B, z);
        a= 0.0;
        Apc.Apply( A, a, b);
        z_xpay( xuneu, w, -1.0, a); // xuneu= w - a;
        xp+= z;
        z_xpaypby2(ru, ru, -1.0, A*VectorCL(xuneu - xu), -1.0, transp_mul( B, z)); // ru-= A*(xuneu - xu) + transp_mul( B, z);
        xu= xuneu;
        resid= std::sqrt( norm_sq( f - A*xu - transp_mul( B, xp)) + norm_sq( g - B*xu));
        std::cerr << "relative residual (2-norm): " << resid/resid0 
                  << "\tv: " << norm( f - A*xu - transp_mul( B, xp))
                  << "\tp: " << norm( g - B*xu)
                  << '\n';
        if (resid<=tol*resid0)
        {
            tol= resid0==0.0 ? 0.0 : resid/resid0;
            max_iter= k;
            return true;
        }
    }
    tol= resid0==0.0 ? 0.0 : resid/resid0;
    return false;
}

class PSchur_AR_CL: public SolverBaseCL
{
  private:
    MGPreCL Apc_;
    ISMGPreCL& Spc_;

  public:
    PSchur_AR_CL(MGDataCL& A_MG, ISMGPreCL& Spc, int outer_iter, double outer_tol)
        :SolverBaseCL( outer_iter, outer_tol), Apc_( A_MG, 1), Spc_( Spc)
        {}
    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        SchurAR( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
};

class PSchur_Diag_AR_CL: public SolverBaseCL
{
  private:
    MGPreCL Apc_;
    DiagMatrixPCCL& Spc_;

  public:
    PSchur_Diag_AR_CL(MGDataCL& A_MG, DiagMatrixPCCL& Spc, int outer_iter, double outer_tol)
        :SolverBaseCL( outer_iter, outer_tol), Apc_( A_MG, 1), Spc_( Spc)
        {}
    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        SchurAR( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
};

double
EigenValueMaxMG(const MGDataCL& A, VectorCL& x, int iter, double  tol= 1e-3)
{
    const_MGDataIterCL finest= --A.end();
    Uint   sm   =  1; // how many smoothing steps?
    int    lvl  = -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
    double l= -1.0;
    double l_old= -1.0;
    VectorCL tmp( x.size());
    VectorCL z( x.size());
//    JORsmoothCL smoother( omega); // Jacobi
//    GSsmoothCL smoother( omega); // Gauss-Seidel
//    SGSsmoothCL smoother( omega); // symmetric Gauss-Seidel
//    SORsmoothCL smoother( omega); // Gauss-Seidel with over-relaxation
    SSORsmoothCL smoother( omega); // symmetric Gauss-Seidel with over-relaxation
//    CGSolverCL  solver( 200, tol); //CG-Verfahren
    SSORPcCL directpc; PCG_SsorCL solver( directpc, 200, 1e-15);
    x/= norm( x);
    std::cerr << "EigenValueMaxMG:\n";
    for (int i= 0; i<iter; ++i) {
        tmp= 0.0;
        MGM( A.begin(), finest, tmp, A.back().A.Data*x, smoother, sm, solver, lvl, -1);
        z= x - tmp;
        l= dot( x, z);
        std::cerr << "iteration: " << i  << "\tlambda: " << l << "\trelative_change= : " << (i==0 ? -1 : std::fabs( (l-l_old)/l_old)) << '\n';
        if (i > 0 && std::fabs( (l-l_old)/l_old) < tol) break;
        l_old= l;
        x= z/norm( z);
    }
    std::cerr << "maximal value for lambda: " << l << '\n';
    return l;
}

//-----------------------------------------------------------------------------
// CG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// max_iter - number of iterations performed before tolerance was reached
//      tol - residual after the final iteration
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
UzawaCGEff(const Mat& A, const Mat& B, Vec& xu, Vec& xp, const Vec& f, const Vec& g,
    PC1& Apc, PC2& Spc,
    int& max_iter, double& tol)
{
    double err= std::sqrt( norm_sq( f - ( A*xu + transp_mul( B, xp))) + norm_sq( g - B*xu));
    const double err0= err;
    Vec rbaru= f - (A*xu + transp_mul(  B, xp));
    Vec rbarp= g - B*xu;
    Vec ru( f.size());
    Apc.Apply( A, ru, rbaru);
    Vec rp= B*ru - rbarp;
    Vec a( f.size()), b( f.size()), s( f.size()), pu( f.size()), qu( f.size());
    Vec z( g.size()), pp( g.size()), qp( g.size()), t( g.size());
    double alpha= 0.0, initialbeta=0.0, beta= 0.0, beta0= 0.0, beta1= 0.0;

    for (int i= 0; i < max_iter; ++i) {
        z= 0.0;
        Spc.Apply( B, z, rp);
        a= A*ru;
        beta1= dot(a, ru) - dot(rbaru,ru) + dot(z,rp);
        if (i==0) initialbeta= beta1;
//        std::cerr << "UzawaCGEff: beta1: " << beta1 << '\n';
        if (beta1 <= 0.0) throw DROPSErrCL( "UzawaCGEff: Matrix is not spd.\n");
        // This is for fair comparisons of different solvers:
        err= std::sqrt( norm_sq( f - (A*xu + transp_mul( B, xp))) + norm_sq( g - B*xu));
        std::cerr << "relative residual (2-norm): " << err/err0 
                  << "\t(problem norm): " << std::sqrt( beta1/initialbeta) << '\n';
//        if (beta1/initialbeta <= tol*tol) {
//            tol= std::sqrt( beta1/initialbeta);
        if (err/err0 <= tol) {
            tol= err/err0;
            max_iter= i;
            return true;
        }
        if (i != 0) {
            beta= beta1/beta0;
            z_xpay( pu, ru, beta, pu); // pu= ru + beta*pu;
            z_xpay( pp, z, beta, pp);  // pp= z  + beta*pp;
            z_xpay( b, a, beta, b);    // b= a + beta*b;
        }
        else {
            pu= ru;
            pp= z;
            b= a; // A*pu;
        }
        qu= b + transp_mul( B, pp);
        qp= B*pu;
        s= 0.0;
        Apc.Apply( A, s, qu);
        z_xpay( t, B*s, -1.0, qp); // t= B*s - qp;
        alpha= beta1/(dot( s, b) - dot(qu, pu) + dot(t, pp));
        axpy( alpha, pu, xu); // xu+= alpha*pu;
        axpy( alpha, pp, xp); // xp+= alpha*pp;
        axpy( -alpha, s, ru); // ru-= alpha*s;
        axpy( -alpha, t, rp); // rp-= alpha*t;
        axpy( -alpha, qu, rbaru);// rbaru-= alpha*qu;
        // rbarp-= alpha*qp; rbarp is not needed in this algo.
        beta0= beta1;
    }
    tol= std::sqrt( beta1);
    return false;
}

//-----------------------------------------------------------------------------
// CG: The return value indicates convergence within max_iter (input)
// iterations (true), or no convergence within max_iter iterations (false).
// max_iter - number of iterations performed before tolerance was reached
//      tol - residual after the final iteration
//-----------------------------------------------------------------------------

template <typename Mat, typename Vec, typename PC1, typename PC2>
bool
UzawaCG(const Mat& A, const Mat& B, Vec& u, Vec& p, const Vec& b, const Vec& c,
        PC1& Apc, PC2& Spc,
        int& max_iter, double& tol)
{
    double err= std::sqrt( norm_sq( b - (A*u + transp_mul( B, p))) + norm_sq( c - B*u));
    const double err0= err;
    Vec ru= b - ( A*u + transp_mul(  B, p));
    Vec rp= c - B*u;
    Vec s1( b.size()); // This is r2u...
    Apc.Apply( A, s1, ru);
    Vec s2= B*s1 - rp;
    Vec r2p( c.size());
    Spc.Apply( B, r2p, s2);
    double rho0= dot( s1, ( A*s1 - ru)) + dot( r2p, s2);
    const double initialrho= rho0;
//    std::cerr << "UzawaCG: rho: " << rho0 << '\n';
    if (rho0<=0.0) throw DROPSErrCL("UzawaCG: Matrix is not spd.\n");
//    tol*= tol*rho0*rho0; // For now, measure the relative error.
//    if (rho0<=tol) {
//        tol= std::sqrt( rho0);
//        max_iter= 0;
//        return true;
//    }
    Vec pu= s1; // s1 is r2u.
    Vec pp= r2p;
    Vec qu= A*pu + transp_mul( B, pp);
    Vec qp= B*pu;
    double rho1= 0.0;
    Vec t1( b.size());
    Vec t2( c.size());
    for (int i= 1; i<=max_iter; ++i) {
        Apc.Apply( A, t1, qu);
        z_xpay( t2, B*t1, -1.0, qp); // t2= B*t1 - qp;
        const double alpha= rho0/( dot(pu, ( A*t1 - qu)) + dot( pp, t2));
        axpy(alpha, pu, u);  // u+= alpha*pu;
        axpy(alpha, pp, p);  // p+= alpha*pp;
        axpy( -alpha, qu, ru);
        axpy( -alpha, qp, rp);
        s1= 0.0;
        Apc.Apply( A, s1, ru);
        z_xpay( s2, B*s1, -1.0, rp); // s2= B*s1 - rp;
        r2p= 0.0;
        Spc.Apply( B, r2p, s2);
        rho1= dot( s1, ( A*s1 - ru)) + dot( r2p, s2);
//        std::cerr << "UzawaCG: rho: " << rho1 << '\n';
        if (rho1<=0.0) throw DROPSErrCL("UzawaCG: Matrix is not spd.\n");
        // This is for fair comparisons of different solvers:
        err= std::sqrt( norm_sq( b - (A*u + transp_mul( B, p))) + norm_sq( c - B*u));
        std::cerr << "relative residual (2-norm): " << err/err0
                  << "\t(problem norm): " << std::sqrt( rho1/initialrho)<< '\n';
//        if (rho1 <= tol) {
//            tol= std::sqrt( rho1);
        if (err <= tol*err0) {
            tol= err/err0;
            max_iter= i;
            return true;
        }
        const double beta= rho1/rho0;
        z_xpay( pu, s1, beta, pu); // pu= s1 + beta*pu; // s1 is r2u.
        z_xpay( pp, r2p, beta, pp); // pp= r2p + beta*pp;
        qu= (A*pu) + transp_mul( B, pp);
        qp= (B*pu);
        rho0= rho1;
        t1= 0.0;
    }
    tol= std::sqrt( rho1);
    return false;
}

template <typename PC1, typename PC2>
class UzawaCGSolverCL : public SolverBaseCL
{
  private:
    PC1& Apc_;
    PC2& Spc_;

  public:
    UzawaCGSolverCL (PC1& Apc, PC2& Spc, int maxiter, double tol)
        : SolverBaseCL(maxiter, tol), Apc_( Apc), Spc_( Spc) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        UzawaCG( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
};

template <typename PC1, typename PC2>
class UzawaCGSolverEffCL : public SolverBaseCL
{
  private:
    PC1& Apc_;
    PC2& Spc_;

  public:
    UzawaCGSolverEffCL (PC1& Apc, PC2& Spc, int maxiter, double tol)
        : SolverBaseCL(maxiter, tol), Apc_( Apc), Spc_( Spc) {}

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        UzawaCGEff( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
};


template <typename PressurePreT, typename PoissonSolver2T>
class MyUzawaSolver2CL : public SolverBaseCL
{
  private:
    PressurePreT&    pr_pre_;
    PoissonSolver2T& poissonSolver2_;
    MatrixCL& M_;
    double    tau_;

  public:
    MyUzawaSolver2CL (PressurePreT& pre, PoissonSolver2T& solver2,
                      MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : SolverBaseCL( maxiter, tol), pr_pre_( pre), poissonSolver2_( solver2),
          M_( M), tau_( tau) {}

    double GetTau()            const { return tau_; }
    void   SetTau( double tau)       { tau_= tau; }

    void Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <class PoissonSolverT, class PoissonSolver2T>
void MyUzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::Solve(
    const MatrixCL& A, const MatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    VectorCL v_corr( v.size()),
             p_corr( p.size()),
             res1( v.size()),
             res2( p.size());
    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= norm_sq( res2);
        pr_pre_.Apply( M_, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(tau_, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2( res1, A*v, 1.0, transp_mul( B, p), -1.0, b);
        res1_norm= norm_sq( res1);
        if (res1_norm + res2_norm < tol) {
            _res= ::sqrt( res1_norm + res2_norm);
            return;
        }
        if( (_iter%output)==0)
            std::cerr << "step " << _iter << ": norm of 1st eq= " << ::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << ::sqrt( res2_norm) << std::endl;

        poissonSolver2_.Apply( A, v_corr, res1);
//        poissonSolver2_.SetTol( std::sqrt( res1_norm)/20.0);
//        poissonSolver2_.Solve( A, v_corr, res1);
//        std::cerr << "velocity: iterations: " << poissonSolver2_.GetIter()
//                  << "\tresidual: " << poissonSolver2_.GetResid() << std::endl;
        v-= v_corr;
    }
    _res= ::sqrt( res1_norm + res2_norm );
}

void
ZeroMean(DROPS::P1EvalCL< double,
                          const DROPS::StokesBndDataCL::PrBndDataCL,
                          DROPS::VecDescCL>& f)
{
    const DROPS::Uint lvl= f.GetSolution()->RowIdx->TriangLevel;
    DROPS::MultiGridCL& mg= const_cast<DROPS::MultiGridCL&>( f.GetMG());
    double MV= 0., vol= 0., sum;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin( lvl),
         send= mg.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        sum= 0.;
        for(int i=0; i<4; ++i)
            sum+= f.val( *sit->GetVertex( i));
        sum/= 120;
        sum+= 2./15.*f.val( *sit, .25, .25, .25);
        MV+= sum * sit->GetVolume()*6.;
        vol+= sit->GetVolume();
    }
    const double c= MV/vol;
    std::cerr << "\nconstant pressure offset: " << c << ", volume of domain: " << vol
              << std::endl;
    for (DROPS::MultiGridCL::TriangVertexIteratorCL sit= mg.GetTriangVertexBegin( lvl),
         send= mg.GetTriangVertexEnd( lvl); sit != send; ++sit) {
        f.SetDoF( *sit, f.val( *sit) - c);
    }
}


class PMinresSP_FullMG_CL : public PMResSPCL<PLanczosONB_SPCL<DROPS::MatrixCL, DROPS::VectorCL, ISMinresMGPreCL> >
{
  private:
    ISMinresMGPreCL pre_;
    PLanczosONB_SPCL<DROPS::MatrixCL, DROPS::VectorCL, ISMinresMGPreCL> q_;

  public:
    PMinresSP_FullMG_CL( DROPS::MGDataCL& MGAvel, DROPS::MGDataCL& MGApr,
                         DROPS::MGDataCL& Mpr, double k_pc,
                         int iter_vel, int iter_prA, int iter_prM, double /*s*/, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<DROPS::MatrixCL, DROPS::VectorCL, ISMinresMGPreCL> >( q_, maxiter, tol),
         pre_( MGAvel, MGApr, Mpr, k_pc, iter_vel, iter_prA, iter_prM, /*s,*/ tol), q_( pre_)
    {}
};


class MyDiagMGPreCL
{
  private:
    const MGDataCL& A_; // Preconditioner for A.
    const MatrixCL& M_; // Preconditioner for S.
    Uint iter_vel_;

  public:
    MyDiagMGPreCL(const MGDataCL& A, const MatrixCL& M, Uint iter_vel)
      :A_( A), M_( M), iter_vel_( iter_vel) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& /*A*/, const Mat& /*B*/, Vec& v, Vec& p, const Vec& b, const Vec& c) const {
//        PA_.SetMaxIter( 1); PA_.SetTol( (bb - K.A_*u).norm()*1e-4);
        Uint   sm   =  1; // how many smoothing steps?
        int    lvl  = -1; // how many levels? (-1=all)
        double omega= 1.; // relaxation parameter for smoother
        SSORsmoothCL smoother( omega);  // Gauss-Seidel with over-relaxation
        SSORPcCL P1;
        SSORsmoothCL P2;
        PCG_SsorCL solver( P1, 200, 1e-12);
        v= 0.0;
        for (DROPS::Uint i= 0; i<iter_vel_; ++i)
            MGM( A_.begin(), --A_.end(), v, b, smoother, sm, solver, lvl, -1);
        for (Uint i= 0; i<M_.num_rows(); ++i) {
            p[i]= c[i]/M_.val( i); // M_ is a diagonal-matrix: exact inversion
        }
    }
};

class MyPMinresSP_DiagMG_CL : public PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, MyDiagMGPreCL> >
{
  private:
    MyDiagMGPreCL pre_;
    PLanczosONB_SPCL<MatrixCL, VectorCL, MyDiagMGPreCL> q_;

  public:
    MyPMinresSP_DiagMG_CL(const MGDataCL& A, const MatrixCL& M, int iter_vel, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, MyDiagMGPreCL> >( q_, maxiter, tol),
         pre_( A, M, iter_vel), q_( pre_)
    {}
};

} // end of namespace DROPS

template<class Coeff>
void
SetupPoissonVelocityMG(
    DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MGDataCL& MGData)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::IdxDescCL* c_idx= 0;
    for(DROPS::Uint lvl= 0; lvl<=mg.GetLastLevel(); ++lvl) {
        MGData.push_back( DROPS::MGLevelDataCL());
        DROPS::MGLevelDataCL& tmp= MGData.back();
        std::cerr << "                        Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 3, 3);
        stokes.CreateNumberingVel( lvl, &tmp.Idx);
//        DROPS::MatDescCL A, M;
//        A.SetIdx( &tmp.Idx, &tmp.Idx);
//        M.SetIdx( &tmp.Idx, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns << std::endl;
        stokes.SetupStiffnessMatrix( &tmp.A);
//        stokes.SetupMassMatrix( &M);
//        tmp.A.Data.LinComb( gamma, M.Data, 1.0, A.Data);
        if(lvl!=0) {
            std::cerr << "                        Create Prolongation on Level " << lvl << std::endl;
            SetupP2ProlongationMatrix( mg, tmp.P, c_idx, &tmp.Idx);
//           std::cout << "    Matrix P " << tmp.P.Data << std::endl;
        }
        c_idx= &tmp.Idx;
    }
    CheckMGData( MGData.begin(), MGData.end());
}


// Assumes, that indices for A_pr are set up. We know, there are only natural
// boundary conditions.
void
SetupPoissonPressure( DROPS::MultiGridCL& mg, DROPS::MatDescCL& A_pr)
{
    DROPS::MatrixBuilderCL A( &A_pr.Data, A_pr.RowIdx->NumUnknowns, A_pr.ColIdx->NumUnknowns);
    const DROPS::Uint lvl= A_pr.RowIdx->TriangLevel;
    const DROPS::Uint idx= A_pr.RowIdx->GetIdx();
    DROPS::SMatrixCL<3,4> G;
    double coup[4][4];
    double det;
    double absdet;
    DROPS::IdxT UnknownIdx[4];

    for (DROPS::MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const DROPS::MultiGridCL&>( mg).GetTriangTetraBegin( lvl),
         send= const_cast<const DROPS::MultiGridCL&>( mg).GetTriangTetraEnd( lvl);
         sit != send; ++sit) {
        DROPS::P1DiscCL::GetGradients( G,det,*sit);
        absdet= fabs( det);
        for(int i=0; i<4; ++i) {
            for(int j=0; j<=i; ++j) {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )/6.0*absdet;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetVertex( i)->Unknowns( idx);
        }
        for(int i=0; i<4; ++i)    // assemble row i
            for(int j=0; j<4;++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i]; 
    }
    A.Build();
    std::cerr << A_pr.Data.num_nonzeros() << " nonzeros in A_pr.\n";
}

// We know, there are only natural boundary conditions.
template<class Coeff>
void
SetupPoissonPressureMG(DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MGDataCL& MGData)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::IdxDescCL* c_idx= 0;
    for(DROPS::Uint lvl= 0; lvl<=mg.GetLastLevel(); ++lvl) {
        MGData.push_back( DROPS::MGLevelDataCL());
        DROPS::MGLevelDataCL& tmp= MGData.back();
        std::cerr << "Pressure-MG:            Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 1);
        stokes.CreateNumberingPr( lvl, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
        SetupPoissonPressure( mg, tmp.A);
        if(lvl!=0) {
            std::cerr << "                        Create Prolongation on Level " << lvl << std::endl;
            DROPS::SetupP1ProlongationMatrix( mg, tmp.P, c_idx, &tmp.Idx);
        }
        c_idx= &tmp.Idx;
    }
    std::cerr << "Check MG-Data..." << std::endl;
    std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
    std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
    CheckMGData( MGData.begin(), MGData.end());
}


template<class Coeff>
void
SetupPressureMassMG(DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MGDataCL& MGData)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::IdxDescCL* c_idx= 0;
    for(DROPS::Uint lvl= 0; lvl<=mg.GetLastLevel(); ++lvl) {
        MGData.push_back( DROPS::MGLevelDataCL());
        DROPS::MGLevelDataCL& tmp= MGData.back();
        std::cerr << "Mass-Pressure-MG:       Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 1);
        stokes.CreateNumberingPr( lvl, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
//        stokes.SetupPrMass( &tmp.A);
        stokes.SetupMass( &tmp.A);
//        tmp.A.Data*= 1e4;
        std::cerr << tmp.A.Data.num_nonzeros() << " nonzeros in M_pr.\n";
        if(lvl!=0) {
            std::cerr << "                        Create Prolongation on Level " << lvl << std::endl;
            DROPS::SetupP1ProlongationMatrix( mg, tmp.P, c_idx, &tmp.Idx);
        }
        c_idx= &tmp.Idx;
    }
    std::cerr << "Check MG-Data..." << std::endl;
    std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
    std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
    CheckMGData( MGData.begin(), MGData.end());
}

template<class Coeff>
void
SetupLumpedPrMass(DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MatDescCL& matM)
{  
    const DROPS::IdxT num_unks_pr=  matM.RowIdx->NumUnknowns;

    DROPS::MatrixBuilderCL M(&matM.Data, num_unks_pr,  num_unks_pr);

    const DROPS::Uint lvl    = matM.RowIdx->TriangLevel;
    const DROPS::Uint pidx   = matM.RowIdx->GetIdx();

    DROPS::IdxT prNumb[4];

    // compute all couplings between HatFunctions on verts:
    // I( i, j) = int ( psi_i*psi_j, T_ref) * absdet
    const double coupl_ii= 1./120. + 2./15./16.,
                 coupl_ij=           2./15./16.;

    for (DROPS::MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const DROPS::MultiGridCL&>(stokes.GetMG()).GetTriangTetraBegin(lvl), send=const_cast<const DROPS::MultiGridCL&>(stokes.GetMG()).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6;
        
        for(int i=0; i<4; ++i)
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                    M( prNumb[i], prNumb[i])+= (i==j ? coupl_ii : coupl_ij) * absdet;
    }
    M.Build();
}

namespace DROPS
{

class PSchur2_PCG_Pr_CL: public PSchurSolver2CL<PCG_SsorCL,
                                                PCGSolverCL<ISPreCL> >
{
  private:
    PCG_SsorCL           PCGsolver_;
    PCGSolverCL<ISPreCL> PCGsolver2_;

  public:
    PSchur2_PCG_Pr_CL(ISPreCL& Spc, int outer_iter, double outer_tol,
                                    int inner_iter, double inner_tol)
        : PSchurSolver2CL<PCG_SsorCL, PCGSolverCL<ISPreCL> >(
              PCGsolver_, PCGsolver2_, outer_iter, outer_tol
          ),
          PCGsolver_( SSORPcCL( 1.), inner_iter, inner_tol),
          PCGsolver2_( Spc, outer_iter, outer_tol)
        {}
};

class PSchur2_PCG_Pr_MG_CL: public PSchurSolver2CL<PCG_SsorCL,
                                                   PCGSolverCL<ISMGPreCL> >
{
  private:
    PCG_SsorCL             PCGsolver_;
    PCGSolverCL<ISMGPreCL> PCGsolver2_;

  public:
    PSchur2_PCG_Pr_MG_CL(ISMGPreCL& Spc, int outer_iter, double outer_tol,
                                         int inner_iter, double inner_tol)
        : PSchurSolver2CL<PCG_SsorCL, PCGSolverCL<ISMGPreCL> >(
              PCGsolver_, PCGsolver2_, outer_iter, outer_tol
          ),
          PCGsolver_( SSORPcCL( 1.), inner_iter, inner_tol),
          PCGsolver2_( Spc, outer_iter, outer_tol)
        {}
};

class PSchur2_Full_MG_CL: public PSchurSolver2CL<MGSolverCL,
                                                 PCGSolverCL<ISMGPreCL> >
{
  private:
    MGSolverCL             solver_;
    PCGSolverCL<ISMGPreCL> solver2_;

  public:
    PSchur2_Full_MG_CL(MGDataCL& A_MG, ISMGPreCL& Spc,
                       int outer_iter, double outer_tol,
                       int inner_iter, double inner_tol)
        : PSchurSolver2CL<MGSolverCL, PCGSolverCL<ISMGPreCL> >(
              solver_, solver2_, outer_iter, outer_tol),
          solver_( A_MG, inner_iter, inner_tol),
          solver2_( Spc, outer_iter, outer_tol)
        {}
};

} // end of namespace DROPS


typedef StokesCL MyPdeCL;

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&);

int
CheckVel(DROPS::P2EvalCL< DROPS::SVectorCL<3>,
                          const DROPS::StokesVelBndDataCL,
                          DROPS::VelVecDescCL>& fun,
         fun_ptr f)
{
    using namespace DROPS;
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const DROPS::Uint trilevel= fun.GetSolution()->RowIdx->TriangLevel;
    std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0.;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= (fun.val( *sit) - f( sit->GetCoord())).norm();
        if ( std::abs( diff) > vmaxdiff) { ++ret; vmaxdiff= std::abs( diff); v= &*sit; }
    }
    std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin( trilevel),
         theend= mg.GetTriangEdgeEnd( trilevel); sit!=theend; ++sit) {
        diff = (fun.val( *sit, .5) - f( (sit->GetVertex( 0)->GetCoord()
                                        +sit->GetVertex( 1)->GetCoord())*0.5)).norm();
        if (std::abs( diff) > emaxdiff) { ++ret; emaxdiff= std::abs( diff); e= &*sit; }
    }
    {
        std::cout << "maximale Differenz Vertices: " << vmaxdiff << " auf\n";
        if (v) v->DebugInfo( std::cout);
        std::cout << "maximale Differenz Edges: " << emaxdiff << " auf\n";
        if (e) e->DebugInfo( std::cout);
        std::cout << std::endl;
    }
    return ret;
}


template<class Coeff>
void
SetMatVecIndices(DROPS::StokesP2P1CL<Coeff>& NS,
                 DROPS::IdxDescCL* const vidx,
                 DROPS::IdxDescCL* const pidx)
{
    std::cout << "#Druck-Unbekannte: " << pidx->NumUnknowns << std::endl;
    std::cout << "#Geschwindigkeitsunbekannte: " << vidx->NumUnknowns << std::endl;
    NS.b.SetIdx( vidx);
    NS.c.SetIdx( pidx);
    NS.A.SetIdx( vidx, vidx);
    NS.B.SetIdx( pidx, vidx);
}

template<class Coeff>
void
ResetSystem(DROPS::StokesP2P1CL<Coeff>& NS)
{
    NS.A.Reset(); NS.B.Reset();
    NS.b.Reset(); NS.c.Reset();
}

void
PrepareStart( DROPS::VelVecDescCL* v, DROPS::VecDescCL*p, DROPS::MatDescCL* M)
{
    for (DROPS::Uint i= 0; i< v->Data.size(); ++i) {
        v->Data[i]= std::cos( double( i));
    }
    v->Data/= norm( v->Data);
    for (DROPS::Uint i= 0; i< p->Data.size(); ++i) {
        p->Data[i]= std::sin( double( i));
    }
    // p must be in L_2^0(\Omega).
    DROPS::VectorCL ones( 1.0, p->Data.size());
    double c= std::sqrt( dot(M->Data*ones, ones));
    DROPS::VectorCL one= ones/c;
    DROPS::VectorCL oneM= M->Data*one;
    p->Data-= dot( oneM*p->Data, one);
    p->Data/= norm( p->Data);
    std::cout << "SP: " << dot( p->Data, M->Data*one) << ", " << dot( p->Data, one) << std::endl;
}


template<class Coeff>
void
StrategyMRes(DROPS::StokesP2P1CL<Coeff>& NS,
             int stokes_maxiter, double stokes_tol,
             double /*a*/, double /*b*/)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MatDescCL  M_pr;
    MatDescCL  ML_pr;
    MGDataCL MG_Mpr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;

    typedef MyPMinresSP_DiagMG_CL StatsolverCL;
//    typedef PMinresSP_FullMG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;

    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    SetupPoissonVelocityMG( NS, MG_vel);
    SetupPoissonPressureMG( NS, MG_pr);
    ML_pr.SetIdx( pidx1, pidx1);
    M_pr.SetIdx( pidx1, pidx1);
    SetupLumpedPrMass( NS, ML_pr);  
    NS.SetupMass( &M_pr);  
//    SetupPressureMassMG( NS, MG_Mpr);
    PrepareStart( v1, p1, &M_pr);
    statsolver= new StatsolverCL( MG_vel, ML_pr.Data, 1, stokes_maxiter, stokes_tol);
//    VectorCL xx( 1.0, vidx1->NumUnknowns);
//    const double rhoinv=  1.0 - EigenValueMaxMG( MG_vel, xx, 10);
//    statsolver= new PMinresSP_FullMG_CL( MG_vel, MG_pr, MG_Mpr, a,
//                     1, 1, 1, b, stokes_maxiter, stokes_tol);
    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    ML_pr.Reset();
    M_pr.Reset();
    MG_Mpr.clear();
}

template<class Coeff>
void
StrategyUzawaCGEff(DROPS::StokesP2P1CL<Coeff>& NS,
    int stokes_maxiter, double stokes_tol,
    double a, double /*b*/)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MatDescCL  ML_pr;
    MatDescCL  M_pr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    MGDataCL MG_Mpr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;
//    typedef UzawaCGSolverEffCL<ScaledMGPreCL, DiagMatrixPCCL> StatsolverCL;
    typedef UzawaCGSolverEffCL<ScaledMGPreCL, ISMGPreCL> StatsolverCL;
    StatsolverCL* statsolver= 0;
    ISMGPreCL* ispcp= 0;
//    DiagMatrixPCCL* ispcp= 0;
    ScaledMGPreCL* velprep= 0;
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    ML_pr.SetIdx( pidx1, pidx1);
//    SetupLumpedPrMass( NS, ML_pr);
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupMass( &M_pr);
    SetupPoissonVelocityMG( NS, MG_vel);
    SetupPoissonPressureMG( NS, MG_pr);
    SetupPressureMassMG( NS, MG_Mpr);
    ispcp= new ISMGPreCL( MG_pr, MG_Mpr, a, /*b,*/ 1);
//    ispcp= new DiagMatrixPCCL( ML_pr.Data);
    PrepareStart( v1, p1, &M_pr);
    VectorCL xx( 1.0, vidx1->NumUnknowns);
//    const double rhoinv= 0.99*( 1.0 - 1.1*EigenValueMaxMG( MG_vel, xx, 1000, 1e-6));
    const double rhoinv= 0.99*( 1.0 - 1.1*0.242869);
    velprep= new ScaledMGPreCL( MG_vel, 1, 1.0/rhoinv);
//    statsolver= new UzawaCGSolverEffCL<ScaledMGPreCL, DiagMatrixPCCL>(
//                        *velprep,
//                        *ispcp,
//                        stokes_maxiter, stokes_tol);
    statsolver= new UzawaCGSolverEffCL<ScaledMGPreCL, ISMGPreCL>(
                        *velprep,
                        *ispcp,
                        stokes_maxiter, stokes_tol);
    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
    DROPS::P1EvalCL<double, const DROPS::StokesBndDataCL::PrBndDataCL,
         DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
    ZeroMean( pr);
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver;
    delete ispcp;
    delete velprep;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
    ML_pr.Reset();
    M_pr.Reset();
}


template<class Coeff>
void
StrategyUzawaCG(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         double a, double /*b*/)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
//    MatDescCL  M_pr;
//    MatDescCL  A_pr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    MGDataCL MG_Mpr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;
    typedef UzawaCGSolverCL<ScaledMGPreCL, ISMGPreCL> StatsolverCL;
    StatsolverCL* statsolver= 0;
    ISMGPreCL* ispcp= 0;
    ScaledMGPreCL* velprep= 0;
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
//    M_pr.SetIdx( pidx1, pidx1);
//    NS.SetupMass( &M_pr);
//    M_pr.Data*= 1e2;
//    A_pr.SetIdx( pidx1, pidx1);
//    SetupPoissonPressure( mg, A_pr);
//    ispcp= new ISPreCL( A_pr.Data, M_pr.Data, k_pc, 1.0);
    SetupPoissonVelocityMG( NS, MG_vel);
    SetupPoissonPressureMG( NS, MG_pr);
    SetupPressureMassMG( NS, MG_Mpr);
    ispcp= new ISMGPreCL( MG_pr, MG_Mpr, a, /*b,*/ 1);
//    PCGSolverCL<ISPreCL> sol1( ispc, stokes_maxiter, stokes_tol);
//    PCG_SsorCL sol2( SSORPcCL( 1.0), stokes_maxiter, stokes_tol);
    VectorCL xx( 1.0, vidx1->NumUnknowns);
    const double rhoinv= 0.95*( 1.0 - 1.1*EigenValueMaxMG( MG_vel, xx, 50, 2e-2));
    velprep= new ScaledMGPreCL( MG_vel, 1, 1.0/rhoinv);
    statsolver= new UzawaCGSolverCL<ScaledMGPreCL, ISMGPreCL>(
                        *velprep,
                        *ispcp,
                        stokes_maxiter, stokes_tol);
    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
    DROPS::P1EvalCL<double, const DROPS::StokesBndDataCL::PrBndDataCL,
         DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
    ZeroMean( pr);
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver;
    delete ispcp;
    delete velprep;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
//    M_pr.Reset();
//    A_pr.Reset();
}


template<class Coeff>
void
StrategyUzawa(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         double a, double /*b*/)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MatDescCL  M_pr;
//    MatDescCL  A_pr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    MGDataCL MG_Mpr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;
//    typedef MyUzawaSolver2CL<ISPreCL, PCG_SsorCL> StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPreCL> > StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, MGPreCL> StatsolverCL;
    typedef MyUzawaSolver2CL<ISMGPreCL, MGPreCL> StatsolverCL;
    StatsolverCL* statsolver= 0;
    ISMGPreCL* ispcp= 0;
//    ISPreCL* ispcp= 0;
    MGPreCL* velprep= 0;
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupMass( &M_pr);  
//    A_pr.SetIdx( pidx1, pidx1);
//    SetupPoissonPressure( mg, A_pr);
//    ispcp= new ISPreCL( A_pr.Data, M_pr.Data, k_pc, 1.0);
    SetupPoissonVelocityMG( NS, MG_vel);
    SetupPoissonPressureMG( NS, MG_pr);
    SetupPressureMassMG( NS, MG_Mpr);
    ispcp= new ISMGPreCL( MG_pr, MG_Mpr, a, /*b,*/ 1);
//    PCGSolverCL<ISPreCL> sol1( ispc, stokes_maxiter, stokes_tol);
//    PCG_SsorCL sol2( SSORPcCL( 1.0), stokes_maxiter, stokes_tol);
    velprep= new MGPreCL( MG_vel, 1);
//    PCGSolverCL<MGPreCL> sol2( velpre, stokes_maxiter, stokes_tol);
//    statsolver= new MyUzawaSolver2CL<ISPreCL, PCG_SsorCL>(
//                       ispc,
//                        sol2,
//                        M_pr.Data, stokes_maxiter, stokes_tol);
//    statsolver= new MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPreCL> >(
//                        ispc,
//                        sol2,
//                        M_pr.Data, stokes_maxiter, stokes_tol);
//    statsolver= new MyUzawaSolver2CL<ISPreCL, MGPreCL>(
//                        *ispcp,
//                        *velprep,
//                        M_pr.Data, stokes_maxiter, stokes_tol);
    statsolver= new MyUzawaSolver2CL<ISMGPreCL, MGPreCL>(
                        *ispcp,
                        *velprep,
                        M_pr.Data, stokes_maxiter, stokes_tol, 1.0/a);
    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
    DROPS::P1EvalCL<double, const DROPS::StokesBndDataCL::PrBndDataCL,
         DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
    ZeroMean( pr);
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver;
    delete ispcp;
    delete velprep;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
    M_pr.Reset();
//    A_pr.Reset();
}

template<class Coeff>
void
StrategyAR(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         int /*poi_maxiter*/, double /*poi_tol*/,
         double a, double /*b*/)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MatDescCL  ML_pr;
    MatDescCL  M_pr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    MGDataCL MG_Mpr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;

    typedef PSchur_AR_CL StatsolverCL;
//    typedef PSchur_Diag_AR_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    ML_pr.SetIdx( pidx1, pidx1);
//    SetupLumpedPrMass( NS, ML_pr);  
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupMass( &M_pr);  
    SetupPoissonVelocityMG( NS, MG_vel);
    SetupPoissonPressureMG( NS, MG_pr);
    SetupPressureMassMG( NS, MG_Mpr);
    PrepareStart( v1, p1, &M_pr);
    ISMGPreCL ispc( MG_pr, MG_Mpr, a, /*b,*/ 1);
//    DiagMatrixPCCL ispc( ML_pr.Data);
    statsolver= new PSchur_AR_CL( MG_vel, ispc, stokes_maxiter, stokes_tol);
//    statsolver= new PSchur_Diag_AR_CL( MG_vel, ispc, stokes_maxiter, stokes_tol);
    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
    DROPS::P1EvalCL<double, const DROPS::StokesBndDataCL::PrBndDataCL,
         DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
    ZeroMean( pr);
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
    ML_pr.Reset();
    M_pr.Reset();
}


template<class Coeff>
void
Strategy(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         int poi_maxiter, double poi_tol,
         double k, double /*gamma*/)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MatDescCL  M_pr;
//    MatDescCL  A_pr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    MGDataCL MG_Mpr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;

//    typedef PSchur_PCG_CL StatsolverCL; 
//    typedef PSchur2_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_MG_CL StatsolverCL;
    typedef PSchur2_Full_MG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupMass( &M_pr);  
//    A_pr.SetIdx( pidx1, pidx1);
//    SetupPoissonPressure( mg, A_pr);
//    ISPreCL ispc( A_pr.Data, M_pr.Data, k_pc, 10);
    SetupPoissonVelocityMG( NS, MG_vel);
    SetupPoissonPressureMG( NS, MG_pr);
    SetupPressureMassMG( NS, MG_Mpr);
    ISMGPreCL ispc( MG_pr, MG_Mpr, k, /*gamma,*/ 1);
//    statsolver= new PSchur_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                   poi_maxiter, poi_tol);
//    statsolver= new PSchur2_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                    poi_maxiter, poi_tol);
//    statsolver= new PSchur2_PCG_Pr_CL( ispc, stokes_maxiter, stokes_tol,
//                                       poi_maxiter, poi_tol);
//    statsolver= new PSchur2_PCG_Pr_MG_CL( ispc, stokes_maxiter, stokes_tol,
//                                          poi_maxiter, poi_tol);
    statsolver= new PSchur2_Full_MG_CL( MG_vel, ispc, stokes_maxiter, stokes_tol,
                                        poi_maxiter, poi_tol);
    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    DROPS::P1EvalCL<double, const DROPS::StokesBndDataCL::PrBndDataCL,
         DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
    ZeroMean( pr);
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
    M_pr.Reset();
//    A_pr.Reset();
}


int main (int argc, char** argv)
{
  try
  {
    // No C-IO here
    std::ios::sync_with_stdio(false);
    if (argc!=10) {
        std::cerr <<
"Usage (stsdrops): <stokes_maxiter> <stokes_tol> <poi_maxiter> <poi_tol>\n"
"    <a> <b> <gamma> <level> <method>" << std::endl;
        return 1;
    }

    int stokes_maxiter= atoi( argv[1]);
    double stokes_tol= atof( argv[2]);
    int poi_maxiter= atoi( argv[3]);
    double poi_tol= atof( argv[4]);
    double a= atof( argv[5]);
    double b= atof( argv[6]);
    double gamma= atof( argv[7]);
    int level= atoi( argv[8]);
    int method= atoi( argv[9]);
    std::cerr << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cerr << "stokes_tol: " << stokes_tol << ", ";
    std::cerr << "poi_maxiter: " << poi_maxiter << ", ";
    std::cerr << "poi_tol: " << poi_tol << ", ";
    std::cerr << "a: " << a << ", ";
    std::cerr << "b: " << b << ", ";
    std::cerr << "gamma: " << gamma << ", ";
    std::cerr << "level: " << level << ", ";
    std::cerr << "method: " << method << std::endl;

    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
				DROPS::std_basis<3>(2),
				DROPS::std_basis<3>(3),
				2, 2, 2);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
	  &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };
    StokesCL::g_= gamma;
    typedef DROPS::StokesP2P1CL<MyPdeCL::StokesCoeffCL> 
    	    NSOnBrickCL;
    typedef NSOnBrickCL MyStokesCL;
    MyStokesCL prob( brick, MyPdeCL::StokesCoeffCL(),
                     DROPS::StokesBndDataCL( 6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    for (int i= 0; i<level; ++i) {
        DROPS::MarkAll( mg);
        mg.Refine();
    }

    switch (method) {
      case 0:
        StrategyUzawa( prob, stokes_maxiter, stokes_tol,
                       a, b);
        break;
      case 1:
        Strategy( prob, stokes_maxiter, stokes_tol, poi_maxiter, poi_tol,
                  a, b);
        break;
      case 2:
        StrategyMRes( prob, stokes_maxiter, stokes_tol,
                      a, b);
        break;
      case 3:
        StrategyUzawaCG( prob, stokes_maxiter, stokes_tol,
                         a, b);
        break;
      case 4:
        StrategyUzawaCGEff( prob, stokes_maxiter, stokes_tol,
                            a, b);
        break;
      case 5:
        StrategyAR( prob, stokes_maxiter, stokes_tol, poi_maxiter, poi_tol,
                    a, b);
        break;
      default:
        std::cerr << "Unknown method.\n";
        break;
    }
    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL( mg) << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
