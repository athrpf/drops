#include "geom/multigrid.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "geom/builder.h"
#include "num/stokessolver.h"
#include "stokes/stokes.h"
#include "stokes/integrTime.h"
#include <fstream>
#include <sstream>


struct InstatStokesCL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])*t*t;
        ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])*t*t;
        ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*t;
        return ret/3.;
    }

    static double LsgPr(const DROPS::Point3DCL& p, double t)
    {
        return std::cos(p[0])*std::sin(p[1])*std::sin(p[2])*t*t
               -(std::sin( 1.) -2.*std::sin( 1.)*std::cos( 1.) + std::sin( 1.)*std::pow( std::cos( 1.), 2))*t*t; // (...)==0.1778213062
    }

    // du/dt + q*u - nu*laplace u + Dp = f
    //                          -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double t)
        {
            DROPS::SVectorCL<3> ret;
            ret[0]= 2./3.*t*std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
            ret[1]= -2./3.*t*std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
            ret[2]= std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*(4./3. + 3.*t);
            return ret;
        }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };
    static StokesCoeffCL Coeff;
};

InstatStokesCL::StokesCoeffCL InstatStokesCL::Coeff;


namespace DROPS
{
template <class StokesT, class SolverT>
class MyInstatStokesThetaSchemeCL
//**************************************************************************
//  for solving the instationary Stokes equation of type StokesT with a
//  1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
//                       theta=1/2 -> Crank-Nicholson (Trapezregel)
//
//  Inner stationary Stokes-type problems are solved with a SolverT-solver.
//  The matrices A, B, M and the rhs b, c of the Stokes class have to be set
//  properly! After construction, SetTimeStep has to be called once. Then
//  every DoStep performs one step in time. Changing time steps require
//  further calls to SetTimeStep.
//**************************************************************************
{
  private:
    StokesT& _Stokes;
    SolverT& _solver;

    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MatrixCL      _mat;               // (1./dt)*M + theta*A

    double _theta, _dt;

  public:
    MyInstatStokesThetaSchemeCL( StokesT& Stokes, SolverT& solver, double theta= 0.5)
        : _Stokes( Stokes), _solver( solver), _b( &Stokes.b), _old_b( new VelVecDescCL),
          _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns),
          _theta( theta)
    {
        _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
        _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    }

    ~MyInstatStokesThetaSchemeCL()
    {
        if (_old_b == &_Stokes.b)
            delete _b;
        else
            delete _old_b;
        delete _cplM; delete _old_cplM;
    }

    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
        _dt= dt;
        _mat.LinComb( 1./dt, _Stokes.M.Data, _theta, _Stokes.A.Data);
    }

    void DoStep( VectorCL& v, VectorCL& p);
};

template <class StokesT, class SolverT>
void MyInstatStokesThetaSchemeCL<StokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * v;
    _rhs*= (_theta-1.);
    _rhs+= (1./_dt)*(_Stokes.M.Data*v + _cplM->Data - _old_cplM->Data)
         +  _theta*_b->Data + (1.-_theta)*_old_b->Data;

    _solver.Solve( _mat, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
}


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
                          const DROPS::StokesPrBndDataCL,
                          DROPS::VecDescCL>& f)
{
    const DROPS::Uint lvl= f.GetLevel();
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
                         DROPS::MGDataCL& Mpr, double kA, double kM,
                         int iter_vel, int iter_prA, int iter_prM, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<DROPS::MatrixCL, DROPS::VectorCL, ISMinresMGPreCL> >( q_, maxiter, tol),
         pre_( MGAvel, MGApr, Mpr, kA, kM, iter_vel, iter_prA, iter_prM, tol), q_( pre_)
    {}
};


} // end of namespace DROPS


template<class Coeff>
void
SetupPoissonVelocityMG(
    DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MGDataCL& MGData,
    const double theta, const double dt)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::IdxDescCL* c_idx= 0;
    for(DROPS::Uint lvl= 0; lvl<=mg.GetLastLevel(); ++lvl) {
        MGData.push_back( DROPS::MGLevelDataCL());
        DROPS::MGLevelDataCL& tmp= MGData.back();
        std::cerr << "                        Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 3, 3);
        stokes.CreateNumberingVel( lvl, &tmp.Idx);
        DROPS::MatDescCL A, M;
        A.SetIdx( &tmp.Idx, &tmp.Idx);
        M.SetIdx( &tmp.Idx, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns << std::endl;
        stokes.SetupStiffnessMatrix( &A);
        stokes.SetupMassMatrix( &M);
        tmp.A.Data.LinComb( 1./dt, M.Data, theta, A.Data);
        if(lvl!=0) {
            std::cerr << "                        Create Prolongation on Level " << lvl << std::endl;
            SetupP2ProlongationMatrix( mg, tmp.P, c_idx, &tmp.Idx);
//           std::cout << "    Matrix P " << tmp.P.Data << std::endl;
        }
        c_idx= &tmp.Idx;
    }
}


// Assumes, that indices for A_pr are set up. We know, there are only natural
// boundary conditions.
void
SetupPoissonPressure( DROPS::MultiGridCL& mg, DROPS::MatDescCL& A_pr)
{
    DROPS::MatrixBuilderCL A( &A_pr.Data, A_pr.RowIdx->NumUnknowns, A_pr.ColIdx->NumUnknowns);
    const DROPS::Uint lvl= A_pr.GetRowLevel();
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
    DROPS::CheckMGData( MGData.begin(), MGData.end());
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
        stokes.SetupPrMass( &tmp.A);
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
    DROPS::CheckMGData( MGData.begin(), MGData.end());
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


// Interface to EnsightP2SolOutCL for time-dependent geometry,
// velocity and pressure
class EnsightWriterCL
{
  private:
    std::string casefile_;
    std::string geomfile_;
    std::string prfile_;
    std::string velfile_;
    DROPS::MultiGridCL& MG_;
    DROPS::IdxDescCL ensightidx_;
    DROPS::EnsightP2SolOutCL ensight_;
    bool have_idx_;

  public:
    EnsightWriterCL(DROPS::MultiGridCL& MG, DROPS::Uint num_timestep,
                    std::string casefile, std::string geomfile,
                    std::string prfile, std::string velfile);
    ~EnsightWriterCL();

    // Call after a grid-change before writing data. The constructor calls this, too.
    void CreateNumbering(int level= -1);
    // Destroy the index before modifying the grid.
    void DeleteNumbering();

    // To write geometry, pressure and velocity at time t.
    template<class InstatNSCL>
    void
    WriteAtTime(const InstatNSCL& NS, const double t);
};

EnsightWriterCL::EnsightWriterCL(DROPS::MultiGridCL& MG, DROPS::Uint num_timestep,
                                 std::string casefile, std::string geomfile,
                                 std::string prfile, std::string velfile)
    :casefile_( casefile), geomfile_( geomfile), prfile_( prfile), velfile_( velfile),
     MG_( MG), ensight_( MG, &ensightidx_), have_idx_( false)
{
    ensightidx_.Set( 1,1,0,0);
    this->CreateNumbering();
    have_idx_= true;
    ensight_.CaseBegin( casefile_.c_str(), num_timestep);
    ensight_.DescribeGeom( "insa_geometry", geomfile_, true);
    ensight_.DescribeScalar( "p", prfile_, true);
    ensight_.DescribeVector( "v", velfile_, true);
}

EnsightWriterCL::~EnsightWriterCL()
{
    if (!have_idx_)
        std::cerr << "EnsightWriter::~EnsightWriterCL: Error; no index found.\n";
    ensight_.CaseEnd();
    this->DeleteNumbering();
}

template<class InstatNSCL>
void
EnsightWriterCL::WriteAtTime(const InstatNSCL& NS, const double t)
{
    if (!have_idx_)
        throw DROPS::DROPSErrCL( "EnsightWriter::WriteAtTime: Call CreateNumbering first.");
    ensight_.putGeom( geomfile_, t);
    typename InstatNSCL::const_DiscPrSolCL ensightp( &NS.p, &NS.GetBndData().Pr, &MG_);
    ensight_.putScalar( prfile_, ensightp, t);
    typename InstatNSCL::const_DiscVelSolCL ensightv( &NS.v, &NS.GetBndData().Vel, &MG_, t);
    ensight_.putVector( velfile_, ensightv, t);
}

void
EnsightWriterCL::CreateNumbering( int level)
{
    if (have_idx_)
        throw DROPS::DROPSErrCL( "EnsightWriter::CreateIndex: Already done.");
    DROPS::NoBndDataCL<> ensightbnd;
    ensightidx_.TriangLevel= level < 0 ? MG_.GetLastLevel(): level;
    DROPS::CreateNumbOnVertex( ensightidx_.GetIdx(), ensightidx_.NumUnknowns, 1,
                               MG_.GetTriangVertexBegin( ensightidx_.TriangLevel),
                               MG_.GetTriangVertexEnd( ensightidx_.TriangLevel),
                               ensightbnd);
    DROPS::CreateNumbOnEdge( ensightidx_.GetIdx(), ensightidx_.NumUnknowns, 1,
                             MG_.GetTriangEdgeBegin( ensightidx_.TriangLevel),
                             MG_.GetTriangEdgeEnd( ensightidx_.TriangLevel),
                             ensightbnd);
    have_idx_= true;
}

void
EnsightWriterCL::DeleteNumbering()
{
    if (!have_idx_)
        throw DROPS::DROPSErrCL( "EnsightWriter::WriteAtTime: Call CreateNumbering first.");
    DROPS::DeleteNumbOnSimplex( ensightidx_.GetIdx(),
                                MG_.GetAllVertexBegin( ensightidx_.TriangLevel),
                                MG_.GetAllVertexEnd( ensightidx_.TriangLevel));
    DROPS::DeleteNumbOnSimplex( ensightidx_.GetIdx(),
                                MG_.GetAllEdgeBegin( ensightidx_.TriangLevel),
                                MG_.GetAllEdgeEnd( ensightidx_.TriangLevel));
    ensightidx_.NumUnknowns= 0;
    have_idx_= false;
}


typedef InstatStokesCL MyPdeCL;

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&, double);

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
    const double t= fun.GetTime();
    const DROPS::Uint trilevel= fun.GetLevel();
    std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0.;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= (fun.val( *sit) - f( sit->GetCoord(), t)).norm();
        if ( std::abs( diff) > vmaxdiff) { ++ret; vmaxdiff= std::abs( diff); v= &*sit; }
    }
    std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin( trilevel),
         theend= mg.GetTriangEdgeEnd( trilevel); sit!=theend; ++sit) {
        diff = (fun.val( *sit, .5) - f( (sit->GetVertex( 0)->GetCoord()
                                        +sit->GetVertex( 1)->GetCoord())*0.5, t)).norm();
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


const double radiusorbit= 0.3; // Radius of the drops' orbit.
const double radiusdrop= 0.15; // Initial radius of the drop.

// positive outside the drop, negative inside the drop.
double
SignedDistToInterface(const DROPS::Point3DCL& p, double t)
{
   DROPS::Point3DCL c;
   c[0]= 0.5 + radiusorbit*std::cos( 2.*M_PI*t);
   c[1]= 0.5 + radiusorbit*std::sin( 2.*M_PI*t);
   c[2]= 0.5;
   return (p-c).norm() - radiusdrop;
}

typedef double (*signed_dist_fun)(const DROPS::Point3DCL& p, double t);


bool
ModifyGridStep(DROPS::MultiGridCL& mg,
               const signed_dist_fun Dist,
               const double width,         // Thickness of refined shell on each side of the interface
               const DROPS::Uint c_level,  // Outside the shell, use this level
               const DROPS::Uint f_level,  // Inside the shell, use this level
               const double t)             // Time of evaluation
// One step of grid change; returns true if modifications were necessary,
// false, if nothing changed.
{
    using namespace DROPS;
    bool shell_not_ready= false;
        for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
             theend= mg.GetTriangTetraEnd(); it!=theend; ++it) {
            double d= 1.;
            for (Uint j=0; j<4; ++j) {
                d= std::min( d, std::abs( Dist( it->GetVertex( j)->GetCoord(), t)));
            }
            const Uint l= it->GetLevel();
            if (d<=width) { // In the shell; level should be f_level.
                if (l < f_level) {
                    shell_not_ready= true;
                    it->SetRegRefMark();
                }
                else
                    if (l > f_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
            else { // Outside the shell; level should be c_level;
                if (l < c_level) { it->SetRegRefMark(); }
                else
                    if (l> c_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
        }
        mg.Refine();
    return shell_not_ready;
}

template<class Coeff>
void
UpdateTriangulation(DROPS::StokesP2P1CL<Coeff>& NS,
                    const signed_dist_fun Dist,
                    const double t,
                    const double width,         // Thickness of refined shell on eache side of the interface
                    const DROPS::Uint c_level,  // Outside the shell, use this level
                    const DROPS::Uint f_level,  // Inside the shell, use this level
                    DROPS::VelVecDescCL* v1,
                    DROPS::VecDescCL* p1)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    Assert( c_level<=f_level, "UpdateTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= v1->RowIdx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx1= p1->RowIdx;
    IdxDescCL* pidx2= &loc_pidx;
    VelVecDescCL  loc_v;
    VecDescCL     loc_p;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p2= &loc_p;
    vidx2->Set( 3, 3, 0, 0);
    pidx2->Set( 1, 0, 0, 0);
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    const StokesBndDataCL& BndData= NS.GetBndData();
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i) {
        shell_not_ready= ModifyGridStep( mg, Dist, width, c_level, f_level, t);
        // Repair velocity
        std::swap( v2, v1);
        std::swap( vidx2, vidx1);
        NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
        if ( mg.GetLastLevel() != vidx2->TriangLevel) {
            std::cout << "LastLevel: " << mg.GetLastLevel()
                      << " vidx2->TriangLevel: " << vidx2->TriangLevel << std::endl;
            throw DROPSErrCL( "Strategy: Sorry, not yet implemented.");
        }
        v1->SetIdx( vidx1);
        P2EvalCL< SVectorCL<3>, const StokesVelBndDataCL,
                  const VelVecDescCL> funv2( v2, &BndData.Vel, &mg, t);
        RepairAfterRefineP2( funv2, *v1);
        v2->Clear();
        NS.DeleteNumberingVel( vidx2);
//P2EvalCL< SVectorCL<3>, const StokesVelBndDataCL,
//          VelVecDescCL> funv1( v1, &BndData.Vel, &mg, t);
//CheckVel( funv1, &MyPdeCL::LsgVel);
        // Repair pressure
        std::swap( p2, p1);
        std::swap( pidx2, pidx1);
        NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
        p1->SetIdx( pidx1);
        typename StokesCL::const_DiscPrSolCL oldfunpr( p2, &BndData.Pr, &mg);
        RepairAfterRefineP1( oldfunpr, *p1);
        p2->Clear();
        NS.DeleteNumberingPr( pidx2);
    }
    // We want the solution to be where v1, p1 point to.
    if (v1 == &loc_v) {
        NS.vel_idx.swap( loc_vidx);
        NS.pr_idx.swap( loc_pidx);
        NS.v.SetIdx( &NS.vel_idx);
        NS.p.SetIdx( &NS.pr_idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
    }
    time.Stop();
    std::cout << "UpdateTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

void
MakeInitialTriangulation(DROPS::MultiGridCL& mg,
                         const signed_dist_fun Dist,
                         const double width,         // Thickness of refined shell on eache side of the interface
                         const DROPS::Uint c_level,  // Outside the shell, use this level
                         const DROPS::Uint f_level)  // Inside the shell, use this level
{
    using namespace DROPS;
    Assert( c_level<=f_level, "MakeInitialTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i)
        shell_not_ready=  ModifyGridStep( mg, Dist, width, c_level, f_level, 0.);
    time.Stop();
    std::cout << "MakeTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
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
    NS.M.SetIdx( vidx, vidx);
}

template<class Coeff>
void
ResetSystem(DROPS::StokesP2P1CL<Coeff>& NS)
{
    NS.A.Reset(); NS.B.Reset();
    NS.M.Reset();
    NS.b.Reset(); NS.c.Reset();
}


template<class Coeff>
void
StrategyMRes(DROPS::StokesP2P1CL<Coeff>& NS,
             int stokes_maxiter, double stokes_tol,
             double theta,
             DROPS::Uint num_timestep,
             double kA, double kM,
             double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
//    MatDescCL  M_pr;
    MGDataCL MG_pr;
    MGDataCL MG_vel;
    MGDataCL MG_Mpr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;

//    typedef PMinresSP_Id_CL StatsolverCL;
    typedef PMinresSP_FullMG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    typedef InstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
//    typedef MyInstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
    InstatsolverCL* instatsolver= 0;

    MakeInitialTriangulation( mg, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);
    EnsightWriterCL ensightout( mg, num_timestep, "insa.case", "insa_geo", "insa_pr", "insa_vel");
    ensightout.WriteAtTime( NS, 0.);
    for (; timestep<num_timestep; ++timestep, t+= dt) {
        std::cerr << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
//        if (timestep%(num_timestep/10) == 0) { // modify the grid
        if (timestep == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                MG_pr.clear();
                MG_vel.clear();
                MG_Mpr.clear();
//                M_pr.Reset();
                ResetSystem( NS);
                ensightout.DeleteNumbering();
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
                ensightout.CreateNumbering();
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cerr << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
//            M_pr.SetIdx( pidx1, pidx1);
//            NS.SetupPrMass( &M_pr);
            SetupPoissonVelocityMG( NS, MG_vel, theta, dt);
            SetupPoissonPressureMG( NS, MG_pr);
            SetupPressureMassMG( NS, MG_Mpr);
//            statsolver= new StatsolverCL( stokes_maxiter, stokes_tol);
            statsolver= new PMinresSP_FullMG_CL( MG_vel, MG_pr, MG_Mpr, kA, kM,
                                                 1, 1, 1, stokes_maxiter, stokes_tol);
            instatsolver= new InstatsolverCL( NS, *statsolver, theta);
        }
        instatsolver->SetTimeStep( dt);
        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,
             DROPS::VecDescCL> pr( p1, &NS.GetBndData().Pr, &mg);
        ZeroMean( pr);
        std::cerr << "Before timestep." << std::endl;
        instatsolver->DoStep( v1->Data, p1->Data);
//        instatsolver->DoStep( v1->Data, p1->Data, timestep==0);
        std::cerr << "After timestep." << std::endl;
        std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, t+dt);
        ensightout.WriteAtTime( NS, t+dt);
    }
    delete instatsolver; instatsolver= 0;
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
//    M_pr.Reset();
}

template<class Coeff>
void
StrategyUzawa(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         double theta,
         DROPS::Uint num_timestep,
         double kA, double kM,
         double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
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
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;
    MGSolverCL mgc (MG_vel, 1, -1.);
    typedef SolverAsPreCL<MGSolverCL> MGPCT;
    MGPCT MGPC (mgc);
//    typedef MyUzawaSolver2CL<ISPreCL, PCG_SsorCL> StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPCT> > StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, MGPCT> StatsolverCL;
    typedef MyUzawaSolver2CL<ISMGPreCL, MGPCT> StatsolverCL;
    StatsolverCL* statsolver= 0;
    typedef InstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
    InstatsolverCL* instatsolver= 0;
    ISMGPreCL* ispcp= 0;
//    ISPreCL* ispcp= 0;
    MakeInitialTriangulation( mg, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    EnsightWriterCL ensightout( mg, num_timestep, "insa.case", "insa_geo", "insa_pr", "insa_vel");
    ensightout.WriteAtTime( NS, 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt) {
        std::cerr << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
//        if (timestep%(num_timestep/10) == 0) { // modify the grid
        if (timestep == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                delete ispcp; ispcp= 0;
                MG_pr.clear();
                MG_vel.clear();
                MG_Mpr.clear();
                M_pr.Reset();
//                A_pr.Reset();
                ResetSystem( NS);
                ensightout.DeleteNumbering();
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
                ensightout.CreateNumbering();
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cerr << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);
//            A_pr.SetIdx( pidx1, pidx1);
//            SetupPoissonPressure( mg, A_pr);
//            ispcp= new ISPreCL( A_pr.Data, M_pr.Data, kA, kM, 1.0);
            SetupPoissonVelocityMG( NS, MG_vel, theta, dt);
            SetupPoissonPressureMG( NS, MG_pr);
            SetupPressureMassMG( NS, MG_Mpr);
            ispcp= new ISMGPreCL( MG_pr, MG_Mpr, kA, kM, 1);
//            PCGSolverCL<ISPreCL> sol1( ispc, stokes_maxiter, stokes_tol);
//            PCG_SsorCL sol2( SSORPcCL( 1.0), stokes_maxiter, stokes_tol);
//            PCGSolverCL<MGPCT> sol2( MGPC, stokes_maxiter, stokes_tol);
//            statsolver= new MyUzawaSolver2CL<ISPreCL, PCG_SsorCL>(
//                               ispc,
//                                sol2,
//                                M_pr.Data, stokes_maxiter, stokes_tol);
//            statsolver= new MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPCT> >(
//                                ispc,
//                                sol2,
//                                M_pr.Data, stokes_maxiter, stokes_tol);
//            statsolver= new MyUzawaSolver2CL<ISPreCL,MGPCT>(
//                                *ispcp,
//                                MGPC,
//                                M_pr.Data, stokes_maxiter, stokes_tol);
            statsolver= new MyUzawaSolver2CL<ISMGPreCL, MGPCT>(
                                *ispcp,
                                MGPC,
                                M_pr.Data, stokes_maxiter, stokes_tol, 1.0);
            instatsolver= new InstatsolverCL( NS, *statsolver, theta);
        }
        instatsolver->SetTimeStep( dt);
        std::cerr << "Before timestep." << std::endl;
        instatsolver->DoStep( v1->Data, p1->Data);
        std::cerr << "After timestep." << std::endl;
        std::cerr << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,
             DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
        ZeroMean( pr);
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, t+dt);
        ensightout.WriteAtTime( NS, t+dt);
    }
    delete instatsolver;
    delete statsolver;
    delete ispcp;
    ResetSystem( NS);
    MG_pr.clear();
    MG_vel.clear();
    MG_Mpr.clear();
    M_pr.Reset();
//    A_pr.Reset();
}

template<class Coeff>
void
Strategy(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         int poi_maxiter, double poi_tol,
         double theta,
         DROPS::Uint num_timestep,
         double kA, double kM,
         double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
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
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;

//    typedef PSchur_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_MG_CL StatsolverCL;
    typedef PSchur2_Full_MG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    typedef InstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
//    typedef MyInstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
    InstatsolverCL* instatsolver= 0;
    MakeInitialTriangulation( mg, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    EnsightWriterCL ensightout( mg, num_timestep, "insa.case", "insa_geo", "insa_pr", "insa_vel");
    ensightout.WriteAtTime( NS, 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt) {
        std::cerr << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
//        if (timestep%(num_timestep/10) == 0) { // modify the grid
        if (timestep == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                MG_pr.clear();
                MG_vel.clear();
                MG_Mpr.clear();
                M_pr.Reset();
//                A_pr.Reset();
                ResetSystem( NS);
                ensightout.DeleteNumbering();
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
                ensightout.CreateNumbering();
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cerr << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);
//            A_pr.SetIdx( pidx1, pidx1);
//            SetupPoissonPressure( mg, A_pr);
//            ISPreCL ispc( A_pr.Data, M_pr.Data, kA, kM, 1.0);
            SetupPoissonVelocityMG( NS, MG_vel, theta, dt);
            SetupPoissonPressureMG( NS, MG_pr);
            SetupPressureMassMG( NS, MG_Mpr);
            ISMGPreCL ispc( MG_pr, MG_Mpr, kA, kM, 1);
//            statsolver= new PSchur_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                           poi_maxiter, poi_tol);
//            statsolver= new PSchur2_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                            poi_maxiter, poi_tol);
//            statsolver= new PSchur2_PCG_Pr_CL( ispc, stokes_maxiter, stokes_tol,
//                                               poi_maxiter, poi_tol);
//            statsolver= new PSchur2_PCG_Pr_MG_CL( ispc, stokes_maxiter, stokes_tol,
//                                                  poi_maxiter, poi_tol);
            statsolver= new PSchur2_Full_MG_CL( MG_vel, ispc, stokes_maxiter, stokes_tol,
                                                poi_maxiter, poi_tol);
            instatsolver= new InstatsolverCL( NS, *statsolver, theta);
        }
        instatsolver->SetTimeStep( dt);
        std::cerr << "Before timestep." << std::endl;
        instatsolver->DoStep( v1->Data, p1->Data);
//        instatsolver->DoStep( v1->Data, p1->Data, timestep==0);
        std::cerr << "After timestep." << std::endl;
        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,
             DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
        ZeroMean( pr);
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, t+dt);
        ensightout.WriteAtTime( NS, t+dt);
    }
    delete instatsolver; instatsolver= 0;
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
    if (argc!=13) {
        std::cerr <<
"Usage (isadrops): <stokes_maxiter> <stokes_tol> <poi_maxiter> <poi_tol>\n"
"    <theta> <num_timestep> <kA> <kM> <shell_width> <c_level> <f_level>\n"
"    <method>" << std::endl;
        return 1;
    }
    // No C-IO here
    std::ios::sync_with_stdio(false);
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                2,2,2);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
          &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };

    int stokes_maxiter= std::atoi( argv[1]);
    double stokes_tol= std::atof( argv[2]);
    int poi_maxiter= std::atoi( argv[3]);
    double poi_tol= std::atof( argv[4]);
    double theta= std::atof( argv[5]);
    int num_timestep= std::atoi( argv[6]);
    double kA= std::atof( argv[7]);
    double kM= std::atof( argv[8]);
    double shell_width= std::atof( argv[9]);
    int c_level= std::atoi( argv[10]);
    int f_level= std::atoi( argv[11]);
    int method= std::atoi( argv[12]);
    std::cerr << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cerr << "stokes_tol: " << stokes_tol << ", ";
    std::cerr << "poi_maxiter: " << poi_maxiter << ", ";
    std::cerr << "poi_tol: " << poi_tol << ", ";
    std::cerr << "theta: " << theta << ", ";
    std::cerr << "num_timestep: " << num_timestep <<  ", ";
    std::cerr << "kA: " << kA <<  ", ";
    std::cerr << "kM: " << kM <<  ", ";
    std::cerr << "shell_width: " << shell_width <<  ", ";
    std::cerr << "c_level: " << c_level << ", ";
    std::cerr << "f_level: " << f_level << ", ";
    std::cerr << "method: " << method << std::endl;

    typedef DROPS::StokesP2P1CL<MyPdeCL::StokesCoeffCL> NSOnBrickCL;
    typedef NSOnBrickCL MyStokesCL;
    MyStokesCL prob( brick, MyPdeCL::StokesCoeffCL(),
                     DROPS::StokesBndDataCL( 6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();

    switch (method) {
      case 0:
        StrategyUzawa( prob, stokes_maxiter, stokes_tol,
                       theta, num_timestep, kA, kM, shell_width, c_level, f_level);
        break;
      case 1:
        Strategy( prob, stokes_maxiter, stokes_tol, poi_maxiter, poi_tol,
                  theta, num_timestep, kA, kM, shell_width, c_level, f_level);
        break;
      case 2:
        StrategyMRes( prob, stokes_maxiter, stokes_tol,
                      theta, num_timestep, kA, kM, shell_width, c_level, f_level);
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
