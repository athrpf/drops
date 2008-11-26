#include "geom/multigrid.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "geom/builder.h"
#include "num/stokessolver.h"
#include "stokes/stokes.h"
#include "stokes/integrTime.h"
#include <fstream>
#include <sstream>


struct StokesCL
{
    static double g_;
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
        ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
        ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
        return ret/3.;
    }

    // Jacobi-matrix od exact solution
    static inline DROPS::SMatrixCL<3, 3> DLsgVel(const DROPS::Point3DCL& p)
    {
        DROPS::SMatrixCL<3, 3> ret;
        ret(0,0)= std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
        ret(0,1)= std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
        ret(0,2)= std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;

        ret(1,0)=   std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
        ret(1,1)=   std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
        ret(1,2)= - std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;

        ret(2,0)= -2.*std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
        ret(2,1)=  2.*std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;
        ret(2,2)= -2.*std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
        return ret;
    }

    static double LsgPr(const DROPS::Point3DCL& p)
    {
        return std::cos(p[0])*std::sin(p[1])*std::sin(p[2])
               -(std::sin( 1.) -2.*std::sin( 1.)*std::cos( 1.) + std::sin( 1.)*std::pow( std::cos( 1.), 2)); // (...)==0.1778213062
    }

    // q*u - nu*laplace u + Dp = f
    //                  -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return StokesCL::g_; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double)
        {
            const double g= StokesCL::g_;
            DROPS::SVectorCL<3> ret;
            ret[0]= g/3.*std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
            ret[1]= -g/3.*std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
            ret[2]= std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*(2./3.*g + 3.);
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
    VectorCL ru( f - A*xu - transp_mul( B, xp));
    VectorCL w( f.size());
    VectorCL z( g.size());
    VectorCL a( f.size());
    VectorCL xuneu( f.size());
    VectorCL b( f.size());
    ApproximateSchurComplMatrixCL<PC1, Mat> asc( A, Apc, B);
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
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc;
    SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > Apc_;
    ISMGPreCL& Spc_;

  public:
    PSchur_AR_CL( ISMGPreCL& Spc, int outer_iter, double outer_tol)
        :SolverBaseCL( outer_iter, outer_tol), smoother_(1.0), coarsesolver_ (SSORPcCL(1.0), 500, 1e-14),
         mgc( smoother_, coarsesolver_, 1, -1., false), Apc_( mgc), Spc_( Spc)
        {}
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        SchurAR( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
    MLMatrixCL* GetPVel() { return mgc.GetProlongation(); }
};

class PSchur_Diag_AR_CL: public SolverBaseCL
{
  private:
      SSORsmoothCL smoother_;
      PCG_SsorCL   coarsesolver_;
      MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc_;
      SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > Apc_;
      DiagMatrixPCCL& Spc_;

  public:
    PSchur_Diag_AR_CL( DiagMatrixPCCL& Spc, int outer_iter, double outer_tol)
        :SolverBaseCL( outer_iter, outer_tol), smoother_(1.0), coarsesolver_ (SSORPcCL(1.0), 500, 1e-14),
         mgc_( smoother_, coarsesolver_, 1, -1., false), Apc_( mgc_), Spc_( Spc)
        {}
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c) {
        _res=  _tol;
        _iter= _maxiter;
        SchurAR( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
    }
    MLMatrixCL* GetPVel() { return mgc_.GetProlongation(); }
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

    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <class PoissonSolverT, class PoissonSolver2T>
void MyUzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B,
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
            _res= std::sqrt( res1_norm + res2_norm);
            return;
        }
        if( (_iter%output)==0)
            std::cerr << "step " << _iter << ": norm of 1st eq= " << std::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << std::sqrt( res2_norm) << std::endl;

        poissonSolver2_.Apply( A, v_corr, res1);
//        poissonSolver2_.SetTol( std::sqrt( res1_norm)/20.0);
//        poissonSolver2_.Solve( A, v_corr, res1);
//        std::cerr << "velocity: iterations: " << poissonSolver2_.GetIter()
//                  << "\tresidual: " << poissonSolver2_.GetResid() << std::endl;
        v-= v_corr;
    }
    _res= std::sqrt( res1_norm + res2_norm );
}

void
ZeroMean(DROPS::P1EvalCL< double,
                          const DROPS::StokesBndDataCL::PrBndDataCL,
                          DROPS::VecDescCL>& f)
{
    const DROPS::Uint lvl= f.GetSolution()->RowIdx->TriangLevel();
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


typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > APcT;
typedef BlockPreCL<APcT, ISMGPreCL> PcT;
typedef PMResSolverCL<PLanczosONBCL<DROPS::VectorCL,PcT> > SolverT;
class PMinresSP_FullMG_CL : public BlockMatrixSolverCL<SolverT>
{
  private:
    SolverT solver_;
    PcT pre_;
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc;
    APcT Apc_;
    ISMGPreCL Spc_;
    PLanczosONBCL<DROPS::VectorCL, PcT> q_;

  public:
    PMinresSP_FullMG_CL( MLMatrixCL& MGApr,
                         MLMatrixCL& Mpr, MLMatrixCL& PPr, double kA, double kM,
                         int iter_vel, int iter_prA, int iter_prM, int maxiter, double tol)
        :BlockMatrixSolverCL<SolverT> ( solver_), solver_( q_, maxiter, tol),
         pre_( Apc_, Spc_), smoother_(1.0), coarsesolver_ ( SSORPcCL(1.0), 500, 1e-14),
         mgc( smoother_, coarsesolver_, iter_vel, 1e-14, false), Apc_( mgc),
         Spc_( MGApr, Mpr, PPr, kA, kM, iter_prA, iter_prM), q_( pre_)
    {}
    MLMatrixCL* GetPVel() { return mgc.GetProlongation(); }
};


class MyDiagPreCL
{
  private:
    const MatrixCL& M_; // Preconditioner for S.

  public:
    MyDiagPreCL(const MatrixCL& M)
      :M_( M) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& /*A*/, Vec& p, const Vec& c) const {
        for (Uint i= 0; i<M_.num_rows(); ++i) {
            p[i]= c[i]/M_.val( i); // M_ is a diagonal-matrix: exact inversion
        }
    }
};

typedef BlockPreCL<APcT, MyDiagPreCL> BlockMGDiagT;
typedef PMResSolverCL<PLanczosONBCL<DROPS::VectorCL, BlockMGDiagT> > MGDiagSolverT;
class MyPMinresSP_DiagMG_CL : public BlockMatrixSolverCL<MGDiagSolverT>
{
  private:
    MGDiagSolverT solver_;
    BlockMGDiagT pre_;
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc;
    APcT Apc_;
    MyDiagPreCL Spc_;
    PLanczosONBCL<DROPS::VectorCL, BlockMGDiagT> q_;

  public:
    MyPMinresSP_DiagMG_CL( const MatrixCL& M, int iter_vel, int maxiter, double tol)
        :BlockMatrixSolverCL<MGDiagSolverT> (solver_), solver_(q_, maxiter, tol),
         pre_( Apc_, Spc_), smoother_(1.0), coarsesolver_ (SSORPcCL(1.0), 500, 1e-14),
         mgc( smoother_, coarsesolver_, iter_vel, 1e-14, false), Apc_( mgc),
         Spc_(M), q_( pre_)
    {}
    MLMatrixCL* GetPVel() { return mgc.GetProlongation(); }
};

} // end of namespace DROPS

// Assumes, that indices for A_pr are set up. We know, there are only natural
// boundary conditions.
void
SetupPoissonPressure( DROPS::MultiGridCL& mg, DROPS::MatrixCL& A_pr, DROPS::IdxDescCL& RowIdx, DROPS::IdxDescCL& ColIdx)
{
    DROPS::MatrixBuilderCL A( &A_pr, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    const DROPS::Uint lvl= RowIdx.TriangLevel();
    const DROPS::Uint idx= RowIdx.GetIdx();
    DROPS::SMatrixCL<3,4> G;
    double coup[4][4];
    double det;
    double absdet;
    DROPS::IdxT UnknownIdx[4];

    for (DROPS::MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const DROPS::MultiGridCL&>( mg).GetTriangTetraBegin( lvl),
         send= const_cast<const DROPS::MultiGridCL&>( mg).GetTriangTetraEnd( lvl);
         sit != send; ++sit) {
        DROPS::P1DiscCL::GetGradients( G,det,*sit);
        absdet= std::fabs( det);
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
    std::cerr << A_pr.num_nonzeros() << " nonzeros in A_pr.\n";
}

// We know, there are only natural boundary conditions.
template<class Coeff>
void
SetupPoissonPressureMG(DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MLMatDescCL& MGData)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::MLIdxDescCL::iterator itRow= MGData.RowIdx->begin();
    DROPS::MLIdxDescCL::iterator itCol= MGData.ColIdx->begin();
    for (DROPS::MLMatrixCL::iterator itA = MGData.Data.begin(); itA != MGData.Data.end(); ++itA)
    {
        SetupPoissonPressure( mg, *itA, *itRow, *itCol);
        ++itRow; ++itCol;
    }
    std::cerr << "Check MG-Data..." << std::endl;
    std::cerr << "                begin     " << MGData.RowIdx->GetCoarsest().NumUnknowns() << std::endl;
    std::cerr << "                end       " << MGData.RowIdx->GetFinest().NumUnknowns() << std::endl;
    //CheckMGData( MGData.begin(), MGData.end());
}

template<class Coeff>
void
SetupLumpedPrMass(DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MLMatDescCL& matM)
{
    const DROPS::IdxT num_unks_pr=  matM.RowIdx->NumUnknowns();

    DROPS::MatrixBuilderCL M(&matM.Data.GetFinest(), num_unks_pr,  num_unks_pr);

    const DROPS::Uint lvl    = matM.RowIdx->TriangLevel();
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

class PSchur2_Full_MG_CL: public PSchurSolver2CL<MGSolverCL<SSORsmoothCL, PCG_SsorCL>,
                                                 PCGSolverCL<ISMGPreCL> >
{
  private:
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> solver_;
    PCGSolverCL<ISMGPreCL> solver2_;

  public:
    PSchur2_Full_MG_CL( ISMGPreCL& Spc,
                       int outer_iter, double outer_tol,
                       int inner_iter, double inner_tol)
        : PSchurSolver2CL<MGSolverCL<SSORsmoothCL, PCG_SsorCL>, PCGSolverCL<ISMGPreCL> >(
          solver_, solver2_, outer_iter, outer_tol),
          smoother_(1.0), coarsesolver_(SSORPcCL(1.0), 500, inner_tol),
          solver_( smoother_, coarsesolver_, inner_iter, inner_tol),
          solver2_( Spc, outer_iter, outer_tol)
        {}
    MLMatrixCL* GetPVel() { return solver_.GetProlongation(); }
};

} // end of namespace DROPS


typedef StokesCL MyPdeCL;

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&);

int
CheckVel(DROPS::P2EvalCL< DROPS::SVectorCL<3>,
                          const DROPS::StokesBndDataCL::VelBndDataCL,
                          DROPS::VelVecDescCL>& fun,
         fun_ptr f)
{
    using namespace DROPS;
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const DROPS::Uint trilevel= fun.GetSolution()->RowIdx->TriangLevel();
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
                 DROPS::MLIdxDescCL* const vidx,
                 DROPS::MLIdxDescCL* const pidx)
{
    std::cout << "#Druck-Unbekannte: " << pidx->NumUnknowns() << std::endl;
    std::cout << "#Geschwindigkeitsunbekannte: " << vidx->NumUnknowns() << std::endl;
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
PrepareStart( DROPS::VelVecDescCL* v, DROPS::VecDescCL*p, DROPS::MLMatDescCL* M)
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
    DROPS::VectorCL one( ones/c);
    DROPS::VectorCL oneM= M->Data*one;
    // XXX ??? p->Data-= dot( oneM*p->Data, one);
    p->Data-= dot( oneM, p->Data);
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
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
    MLMatDescCL  ML_pr, MG_pr;
    MLMatDescCL  PPr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;

    typedef MyPMinresSP_DiagMG_CL StatsolverCL;
//    typedef PMinresSP_FullMG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    NS.SetNumVelLvl( mg.GetNumLevel());
    NS.SetNumPrLvl ( mg.GetNumLevel());
    M_pr.Data.resize( mg.GetNumLevel());
    ML_pr.Data.resize( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    PPr.Data.resize( mg.GetNumLevel());

    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr ( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    PPr.SetIdx( pidx1, pidx1);
    MG_pr.SetIdx( pidx1, pidx1);

    SetupPoissonPressureMG( NS, MG_pr);
    ML_pr.SetIdx( pidx1, pidx1);
    M_pr.SetIdx( pidx1, pidx1);
    SetupLumpedPrMass( NS, ML_pr);
    NS.SetupPrMass( &M_pr);
//    SetupPressureMassMG( NS, MG_Mpr);
    PrepareStart( v1, p1, &M_pr);
    statsolver= new StatsolverCL( ML_pr.Data.GetFinest(), 1, stokes_maxiter, stokes_tol);
    SetupP2ProlongationMatrix( mg, *(statsolver->GetPVel()), vidx1, vidx1);
    SetupP1ProlongationMatrix( mg, PPr);

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
    //ML_pr.Reset();
    M_pr.Reset();
}

template<class Coeff>
void
StrategyUzawaCGEff(DROPS::StokesP2P1CL<Coeff>& NS,
    int stokes_maxiter, double stokes_tol,
    double kA, double kM)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  ML_pr;
    MLMatDescCL  M_pr;
    MLMatDescCL  MG_pr;
    MLMatDescCL  PVel, PPr;

    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
//    typedef UzawaCGSolverEffCL<ScaledMGPreCL<>, DiagMatrixPCCL> StatsolverCL;
    typedef UzawaCGSolverEffCL<ScaledMGPreCL<>, ISMGPreCL> StatsolverCL;
    StatsolverCL* statsolver= 0;
    ISMGPreCL* ispcp= 0;
//    DiagMatrixPCCL* ispcp= 0;
    ScaledMGPreCL<>* velprep= 0;
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    PVel.Data.resize ( mg.GetNumLevel());
    PPr.Data.resize  ( mg.GetNumLevel());
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr ( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    ML_pr.SetIdx( pidx1, pidx1);
//    SetupLumpedPrMass( NS, ML_pr);
    PPr.SetIdx( pidx1, pidx1);
    PVel.SetIdx( vidx1, vidx1);
    SetupP1ProlongationMatrix( mg, PPr);
    SetupP2ProlongationMatrix( mg, PVel);
    M_pr.SetIdx( pidx1, pidx1);
    MG_pr.SetIdx( pidx1, pidx1);
    NS.SetupPrMass( &M_pr);

    SetupPoissonPressureMG( NS, MG_pr);

    ispcp= new ISMGPreCL( MG_pr.Data, M_pr.Data, PPr.Data, kA, kM, 1);
//    ispcp= new DiagMatrixPCCL( ML_pr.Data);
    PrepareStart( v1, p1, &M_pr);
    VectorCL xx( 1.0, vidx1->NumUnknowns());
//    const double rhoinv= 0.99*( 1.0 - 1.1*EigenValueMaxMG( NS.A.Data, PVel.Data, xx, 1000, 1e-6));
    const double rhoinv= 0.99*( 1.0 - 1.1*0.242869);
    velprep= new ScaledMGPreCL<>( PVel.Data, 1, 1.0/rhoinv);
//    statsolver= new UzawaCGSolverEffCL<ScaledMGPreCL<>, DiagMatrixPCCL>(
//                        *velprep,
//                        *ispcp,
//                        stokes_maxiter, stokes_tol);
    statsolver= new UzawaCGSolverEffCL<ScaledMGPreCL<>, ISMGPreCL>(
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
}


template<class Coeff>
void
StrategyUzawaCG(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         double kA, double kM)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
    MLMatDescCL  MG_pr;
    MLMatDescCL  PVel, PPr;
//    MLMatDescCL  M_pr;
//    MLMatDescCL  A_pr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
    typedef UzawaCGSolverCL<ScaledMGPreCL<>, ISMGPreCL> StatsolverCL;
    StatsolverCL* statsolver= 0;
    ISMGPreCL* ispcp= 0;
    ScaledMGPreCL<>* velprep= 0;
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    PVel.Data.resize ( mg.GetNumLevel());
    PPr.Data.resize  ( mg.GetNumLevel());
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr ( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    PPr.SetIdx( pidx1, pidx1);
    PVel.SetIdx( vidx1, vidx1);
    SetupP1ProlongationMatrix( mg, PPr);
    SetupP2ProlongationMatrix( mg, PVel);
    M_pr.SetIdx( pidx1, pidx1);
    MG_pr.SetIdx( pidx1, pidx1);
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupPrMass( &M_pr);
//    M_pr.Data*= 1e2;
//    A_pr.SetIdx( pidx1, pidx1);
//    SetupPoissonPressure( mg, A_pr);
//    ispcp= new ISPreCL( A_pr.Data, M_pr.Data, kA, kM, 1.0);
    SetupPoissonPressureMG( NS, MG_pr);
    ispcp= new ISMGPreCL( MG_pr.Data, M_pr.Data, PPr.Data, kA, kM, 1);
//    PCGSolverCL<ISPreCL> sol1( ispc, stokes_maxiter, stokes_tol);
//    PCG_SsorCL sol2( SSORPcCL( 1.0), stokes_maxiter, stokes_tol);
    VectorCL xx( 1.0, vidx1->NumUnknowns());
    const double rhoinv= 0.95*( 1.0 - 1.1*EigenValueMaxMG( NS.A.Data, PVel.Data, xx, 50, 2e-2));
    velprep= new ScaledMGPreCL<>( PVel.Data, 1, 1.0/rhoinv);
    statsolver= new UzawaCGSolverCL<ScaledMGPreCL<>, ISMGPreCL>(
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
}


template<class Coeff>
void
StrategyUzawa(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         double kA, double kM)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  ML_pr;
    MLMatDescCL  M_pr;
    MLMatDescCL  MG_pr;
    MLMatDescCL  PPr;
//    MLMatDescCL  A_pr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
    SSORsmoothCL smoother(1.0);
    PCG_SsorCL   coarsesolver(SSORPcCL(1.0), 500, 1e-14);
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc (smoother, coarsesolver, 1, -1.0, false);
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > MGPCT;
    MGPCT MGPC (mgc);
//    typedef MyUzawaSolver2CL<ISPreCL, PCG_SsorCL> StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPCT> > StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, MGPCT> StatsolverCL;
    typedef MyUzawaSolver2CL<ISMGPreCL, MGPCT> StatsolverCL;
    StatsolverCL* statsolver= 0;
    ISMGPreCL* ispcp= 0;
//    ISPreCL* ispcp= 0;
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    PPr.Data.resize  ( mg.GetNumLevel());
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr ( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    PPr.SetIdx( pidx1, pidx1);
    SetupP1ProlongationMatrix( mg, PPr);
    SetupP2ProlongationMatrix( mg, *mgc.GetProlongation(), vidx1, vidx1);
    M_pr.SetIdx( pidx1, pidx1);
    MG_pr.SetIdx( pidx1, pidx1);
    ML_pr.SetIdx( pidx1, pidx1);
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupPrMass( &M_pr);
//    A_pr.SetIdx( pidx1, pidx1);
//    SetupPoissonPressure( mg, A_pr);
//    ispcp= new ISPreCL( A_pr.Data, M_pr.Data, kA, kM, 1.0);

    SetupPoissonPressureMG( NS, MG_pr);
    ispcp= new ISMGPreCL( MG_pr.Data, M_pr.Data, PPr.Data, kA, kM, 1);
//    PCGSolverCL<ISPreCL> sol1( ispc, stokes_maxiter, stokes_tol);
//    PCG_SsorCL sol2( SSORPcCL( 1.0), stokes_maxiter, stokes_tol);
//    PCGSolverCL<MGPCT> sol2( MGPC, stokes_maxiter, stokes_tol);
//    statsolver= new MyUzawaSolver2CL<ISPreCL, PCG_SsorCL>(
//                       ispc,
//                        sol2,
//                        M_pr.Data, stokes_maxiter, stokes_tol);
//    statsolver= new MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPCT> >(
//                        ispc,
//                        sol2,
//                        M_pr.Data, stokes_maxiter, stokes_tol);
//    statsolver= new MyUzawaSolver2CL<ISPreCL, MGPCT>(
//                        *ispcp,
//                        MGPC,
//                        M_pr.Data, stokes_maxiter, stokes_tol);
    statsolver= new MyUzawaSolver2CL<ISMGPreCL, MGPCT>(
                        *ispcp,
                        MGPC,
                        M_pr.Data.GetFinest(), stokes_maxiter, stokes_tol, 1.0/kA);
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
    ResetSystem( NS);
//    A_pr.Reset();
}

template<class Coeff>
void
StrategyAR(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         int /*poi_maxiter*/, double /*poi_tol*/,
         double kA, double kM)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  ML_pr;
    MLMatDescCL  M_pr;
    MLMatDescCL MG_pr;
    MLMatDescCL PPr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;

    typedef PSchur_AR_CL StatsolverCL;
//    typedef PSchur_Diag_AR_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    PPr.Data.resize  ( mg.GetNumLevel());

    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr ( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    PPr.SetIdx( pidx1, pidx1);
    M_pr.SetIdx( pidx1, pidx1);
    MG_pr.SetIdx( pidx1, pidx1);
    ML_pr.SetIdx( pidx1, pidx1);
//    SetupLumpedPrMass( NS, ML_pr);
    M_pr.SetIdx( pidx1, pidx1);
    NS.SetupPrMass( &M_pr);
    SetupPoissonPressureMG( NS, MG_pr);
    PrepareStart( v1, p1, &M_pr);
    ISMGPreCL ispc( MG_pr.Data, M_pr.Data, PPr.Data, kA, kM, 1);
//    DiagMatrixPCCL ispc( ML_pr.Data);
    statsolver= new PSchur_AR_CL( ispc, stokes_maxiter, stokes_tol);
    SetupP2ProlongationMatrix( mg, *(statsolver->GetPVel()), vidx1, vidx1);
    SetupP1ProlongationMatrix( mg, PPr);

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
}


template<class Coeff>
void
Strategy(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         int poi_maxiter, double poi_tol,
         double kA, double kM)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
//    MLMatDescCL  A_pr;
    MLMatDescCL MG_pr;
    MLMatDescCL PPr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;

//    typedef PSchur_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_MG_CL StatsolverCL;
    typedef PSchur2_Full_MG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    PPr.Data.resize  ( mg.GetNumLevel());
    
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.CreateNumberingPr ( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    SetMatVecIndices( NS, vidx1, pidx1);
    time.Reset(); time.Start();
    NS.SetupSystem( &NS.A, &NS.b, &NS.B, &NS.c);
    time.Stop();
    std::cerr << "SetupSystem: " << time.GetTime() << " seconds" << std::endl;
    time.Reset();
    PPr.SetIdx( pidx1, pidx1);
    M_pr.SetIdx( pidx1, pidx1);
    MG_pr.SetIdx( pidx1, pidx1);
    NS.SetupPrMass( &M_pr);
//    A_pr.SetIdx( pidx1, pidx1);
//    SetupPoissonPressure( mg, A_pr);
//    ISPreCL ispc( A_pr.Data, M_pr.Data, kA, kM, 1.0);
    SetupPoissonPressureMG( NS, MG_pr);
    
    ISMGPreCL ispc( MG_pr.Data, M_pr.Data, PPr.Data, kA, kM, 1);
//    statsolver= new PSchur_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                   poi_maxiter, poi_tol);
//    statsolver= new PSchur2_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                    poi_maxiter, poi_tol);
//    statsolver= new PSchur2_PCG_Pr_CL( ispc, stokes_maxiter, stokes_tol,
//                                       poi_maxiter, poi_tol);
//    statsolver= new PSchur2_PCG_Pr_MG_CL( ispc, stokes_maxiter, stokes_tol,
//                                          poi_maxiter, poi_tol);
    statsolver= new PSchur2_Full_MG_CL( ispc, stokes_maxiter, stokes_tol,
                                        poi_maxiter, poi_tol);
    SetupP2ProlongationMatrix( mg, *(statsolver->GetPVel()), vidx1, vidx1);
    SetupP1ProlongationMatrix( mg, PPr);

    std::cerr << "Before solve." << std::endl;
    statsolver->Solve( NS.A.Data, NS.B.Data, v1->Data, p1->Data, NS.b.Data, NS.c.Data);
    std::cerr << "After solve." << std::endl;
    DROPS::P1EvalCL<double, const DROPS::StokesBndDataCL::PrBndDataCL,
         DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
    ZeroMean( pr);
    NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DLsgVel, &MyPdeCL::LsgPr);
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
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
"    <kA> <kM> <gamma> <level> <method>" << std::endl;
        return 1;
    }
    int stokes_maxiter= std::atoi( argv[1]);
    double stokes_tol= std::atof( argv[2]);
    int poi_maxiter= std::atoi( argv[3]);
    double poi_tol= std::atof( argv[4]);
    double kA= std::atof( argv[5]);
    double kM= std::atof( argv[6]);
    double gamma= std::atof( argv[7]);
    int level= std::atoi( argv[8]);
    int method= std::atoi( argv[9]);
    std::cerr << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cerr << "stokes_tol: " << stokes_tol << ", ";
    std::cerr << "poi_maxiter: " << poi_maxiter << ", ";
    std::cerr << "poi_tol: " << poi_tol << ", ";
    std::cerr << "kA: " << kA << ", ";
    std::cerr << "kM: " << kM << ", ";
    std::cerr << "gamma: " << gamma << ", ";
    std::cerr << "level: " << level << ", ";
    std::cerr << "method: " << method << std::endl;

    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                2, 2, 2);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
          &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };
    StokesCL::g_= gamma;
    typedef DROPS::StokesP2P1CL<MyPdeCL::StokesCoeffCL> NSOnBrickCL;
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
                       kA, kM);
        break;
      case 1:
        Strategy( prob, stokes_maxiter, stokes_tol, poi_maxiter, poi_tol,
                  kA, kM);
        break;
      case 2:
        StrategyMRes( prob, stokes_maxiter, stokes_tol,
                      kA, kM);
        break;
      case 3:
        StrategyUzawaCG( prob, stokes_maxiter, stokes_tol,
                         kA, kM);
        break;
      case 4:
        StrategyUzawaCGEff( prob, stokes_maxiter, stokes_tol,
                            kA, kM);
        break;
      case 5:
        StrategyAR( prob, stokes_maxiter, stokes_tol, poi_maxiter, poi_tol,
                    kA, kM);
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
