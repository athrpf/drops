#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include "num/MGsolver.h"
#include <fstream>


inline DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    sin(p[0])*sin(p[1])*sin(p[2])/3.;
    ret[1]=  - cos(p[0])*cos(p[1])*sin(p[2])/3.;
    ret[2]= 2.*cos(p[0])*sin(p[1])*cos(p[2])/3.;
    return ret;
}

// Jacobi-matrix od exact solution
inline DROPS::SMatrixCL<3, 3> DLsgVel(const DROPS::Point3DCL& p)
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


// Volume of the box: 0.484473073129685
// int(p)/vol = -0.125208551608365
inline double LsgPr(const DROPS::Point3DCL& p)
{
    return cos(p[0])*sin(p[1])*sin(p[2]) - 0.125208551608365;
}


// q*u - nu*laplace u + Dp = f
//                  -div u = 0
class StokesCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p)
        { DROPS::SVectorCL<3> ret(0.0); ret[2]= 3.*cos(p[0])*sin(p[1])*cos(p[2]); return ret; }
    const double nu;
    
    StokesCoeffCL() : nu(1.0) {}
};

typedef DROPS::StokesP2P1CL<StokesCoeffCL> 
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

namespace DROPS // for Strategy
{

class StokesVectorCL;

class StokesMatrixCL
{
  public:
    const MatrixCL& A_;
    const MatrixCL& B_;

    StokesMatrixCL(const MatrixCL& A, const MatrixCL& B)
      :A_( A), B_( B) {}

    friend StokesVectorCL operator*(const StokesMatrixCL&, const StokesVectorCL&);
};

class StokesVectorCL :public VectorCL
{
  private:
    unsigned char u_[sizeof( VectorCL)];
    unsigned char p_[sizeof( VectorCL)];

  public:
    StokesVectorCL(const VectorCL& u, const VectorCL& p)
      :VectorCL( u.size() + p.size()) {
        std::copy( Addr( u.raw()), Addr( u.raw()) + u.size(), &this->raw()[0]);
        // for (Uint i= 0; i<u.size(); ++i) (*this)[i]= u[i];
        std::copy( Addr( p.raw()), Addr( p.raw()) + p.size(), &this->raw()[0] + u.size());
        // for (Uint i= 0; i<p.size(); ++i) (*this)[i+u.size()]= p[i];
    }
    StokesVectorCL(const VectorCL& v) :VectorCL( v) {}
    StokesVectorCL() :VectorCL() {}
    StokesVectorCL(Uint s) :VectorCL( s) {}

    VectorCL& u(Uint sizeu);
    VectorCL& p(Uint sizep);
    const VectorCL& u(Uint sizeu) const;
    const VectorCL& p(Uint sizep) const;
    friend StokesVectorCL operator*(const StokesMatrixCL&, const StokesVectorCL&);
};

VectorCL& StokesVectorCL::u( Uint sizeu)
{
    *reinterpret_cast<size_t*>( u_)= sizeu;
    *reinterpret_cast<double**>( u_+ sizeof( size_t))= &this->raw()[0];
    return *reinterpret_cast<VectorCL*>( u_);
}

VectorCL& StokesVectorCL::p( Uint sizep)
{
    *reinterpret_cast<size_t*>( p_)= sizep;
    *reinterpret_cast<double**>( p_+ sizeof( size_t))= &this->raw()[0] + this->size() - sizep;
    return *reinterpret_cast<VectorCL*>( p_);
}

const VectorCL& StokesVectorCL::u( Uint sizeu) const
{
    *reinterpret_cast<size_t*>( const_cast<unsigned char*>( u_))= sizeu;
    *reinterpret_cast<double**>( const_cast<unsigned char*>( u_) + sizeof( size_t))= const_cast<double*>( &this->raw()[0]);
    return *reinterpret_cast<VectorCL*>( const_cast<unsigned char*>( u_));
}

const VectorCL& StokesVectorCL::p( Uint sizep) const
{
    *reinterpret_cast<size_t*>( const_cast<unsigned char*>( p_))= sizep;
    *reinterpret_cast<double**>( const_cast<unsigned char*>( p_) + sizeof( size_t))= const_cast<double*>( &this->raw()[0]) + this->size() - sizep;
    return *reinterpret_cast<VectorCL*>( const_cast<unsigned char*>( p_));
}


StokesVectorCL operator*(const StokesMatrixCL& K, const StokesVectorCL& x)
{
//    VectorCL u( x.raw()[std::slice( 0, K.B_.num_cols(), 1)]);
//    VectorCL p( x.raw()[std::slice( K.B_.num_cols(), K.B_.num_rows(), 1)]);
//    return StokesVectorCL( K.A_*u + transp_mul( K.B_, p), K.B_*u);
    StokesVectorCL ret( x.size());
    y_Ax( &ret.raw()[0],
          K.A_.num_rows(),
          K.A_.raw_val(),
          K.A_.raw_row(),
          K.A_.raw_col(),
          Addr( x.raw()));
    // y_ATx is +=, not = ...
    y_ATx( &ret.raw()[0],
           K.B_.num_rows(),
           K.B_.raw_val(),
           K.B_.raw_row(),
           K.B_.raw_col(),
           Addr( x.raw()) + K.A_.num_rows());
    y_Ax( &ret.raw()[0] + K.A_.num_rows(),
          K.B_.num_rows(),
          K.B_.raw_val(),
          K.B_.raw_row(),
          K.B_.raw_col(),
          Addr( x.raw()));
    return ret;
}


void
CheckMatVec( const StokesMatrixCL& K,
             const MatrixCL& A, const MatrixCL& B,
             const StokesVectorCL& x,
             const VectorCL& u, const VectorCL& p)
{
    StokesVectorCL tmp( A*u + transp_mul( B, p), B*u);
    std::cout << "CheckMatVec: " << (tmp - K*x).norm() << std::endl;
}


class FullPreCL
{
  private:
    mutable PCG_SsorCL PA_; // Preconditioner for A.
    const MatrixCL&   PS_; // Preconditioner for S. The pressure-mass matrix, system solved with above PCG.

  public:
    FullPreCL( const MatrixCL& ps)
      :PA_( SSORPcCL(), 5, 1e-20), PS_( ps) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& K, Vec& x, const Vec& b) const {
//        VectorCL u( x.raw()[std::slice( 0, K.B_.num_cols(), 1)]);
//        VectorCL p( x.raw()[std::slice( K.B_.num_cols(), K.B_.num_rows(), 1)]);
//        const VectorCL bb( b.raw()[std::slice( 0, K.B_.num_cols(), 1)]);
//        const VectorCL bc( b.raw()[std::slice( K.B_.num_cols(), K.B_.num_rows(), 1)]);
        SSORPcCL P1;
        VectorCL& u= x.u( K.B_.num_cols());
        VectorCL& p= x.p( K.B_.num_rows());
        const VectorCL& bb= b.u( K.B_.num_cols());
        const VectorCL& bc= b.p( K.B_.num_rows());

//        SchurComplMatrixCL S( K.A_, K.B_, 1e-10, 1.);
//        std::cerr << (bb - K.A_*u).norm() << '\t';
        PA_.SetMaxIter( 500);PA_.SetTol( (bb - K.A_*u).norm()*1e-4);
        PA_.Solve( K.A_, u, bb);
//        std::cerr << PA_.GetIter() << '\t' << PA_.GetResid() << '\n';
//        std::cerr << (bc - PS_*p).norm() << '\t';
//        PA_.SetTol( (bc - PS_*p).norm()*1e-1);
//        PA_.SetMaxIter( 3); PA_.SetTol( 1e-20);
//        PA_.Solve( S, p, bc);
//        CGSolverCL cgs( 50, 1e-20);
//        cgs.Solve( S, p, bc);
        P1.Apply( PS_, p, bc);
//        std::cerr << PA_.GetIter() << '\t' << PA_.GetResid() << "\n\n";
//        std::cerr << cgs.GetIter() << '\t' << cgs.GetResid() << "\n\n";
//        std::copy( &u[0], &u[0] + u.size(), &x[0]);
//        std::copy( &p[0], &p[0] + p.size(), &x[0] + u.size());
    }
};

class MGPreCL
{
  private:
//    mutable MGSolverCL PA_; // Preconditioner for A.
    const MGDataCL& mgd_; // Preconditioner for A.
    const MatrixCL& PS_; // Preconditioner for S.

  public:
//    MGPreCL( MGSolverCL& pa, const MatrixCL& ps)
//      :PA_( pa), PS_( ps) {}
    MGPreCL( const MGDataCL& mgd, const MatrixCL& ps)
      :mgd_( mgd), PS_( ps) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& K, Vec& x, const Vec& b) const {
        SSORPcCL P1;
        VectorCL& u= x.u( K.B_.num_cols());
        VectorCL& p= x.p( K.B_.num_rows());
        const VectorCL& bb= b.u( K.B_.num_cols());
        const VectorCL& bc= b.p( K.B_.num_rows());
//        PA_.SetMaxIter( 1); PA_.SetTol( (bb - K.A_*u).norm()*1e-4);
//        PA_.Solve( K.A_, u, bb);
        Uint   sm   =  2; // how many smoothing steps?
        int    lvl  = -1; // how many levels? (-1=all)
        double omega= 1.; // relaxation parameter for smoother
        SORsmoothCL smoother( omega);  // Gauss-Seidel with over-relaxation
        PCG_SsorCL solver( P1, 200, 1e-12);
        MGM( mgd_.begin(), --mgd_.end(), u, bb, smoother, sm, solver, lvl, -1);
        MGM( mgd_.begin(), --mgd_.end(), u, bb, smoother, sm, solver, lvl, -1);
        P1.Apply( PS_, p, bc);
    }
};


class SSORPreCL
{
  private:
    const MatrixCL&  PS_; // Preconditioner for S. The pressure-mass matrix, system solved with above PCG.

  public:
    SSORPreCL( const MatrixCL& ps)
      : PS_( ps) {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& K, Vec& x, const Vec& b) const {
        SSORPcCL P1;
        SSORsmoothCL P2;
//        VectorCL u( x.raw()[std::slice( 0, K.B_.num_cols(), 1)]);
//        VectorCL p( x.raw()[std::slice( K.B_.num_cols(), K.B_.num_rows(), 1)]);
//        const VectorCL bb( b.raw()[std::slice( 0, K.B_.num_cols(), 1)]);
//        const VectorCL bc( b.raw()[std::slice( K.B_.num_cols(), K.B_.num_rows(), 1)]);
        VectorCL& u= x.u( K.B_.num_cols());
        VectorCL& p= x.p( K.B_.num_rows());
        const VectorCL& bb= b.u( K.B_.num_cols());
        const VectorCL& bc= b.p( K.B_.num_rows());
        P1.Apply( K.A_, u, bb);
        for (int i=0; i<0; ++i) {
            P2.Apply( K.A_, u, bb);
        }
        P1.Apply( PS_, p, bc);
        for (int i=0; i<0; ++i) {
            P2.Apply( PS_, p, bc);
        }
//        std::copy( &u[0], &u[0] + u.size(), &x[0]);
//        std::copy( &p[0], &p[0] + p.size(), &x[0] + u.size());
    }
};


using ::MyStokesCL;

template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, double omega, double inner_iter_tol, double tol, int meth,
                                           Uint maxStep, double rel_red, double markratio,
                                           double tau, Uint uzawa_inner_iter)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();
    const typename MyStokesCL::BndDataCL::PrBndDataCL& PrBndData= Stokes.GetBndData().Pr;
    const typename MyStokesCL::BndDataCL::VelBndDataCL& VelBndData= Stokes.GetBndData().Vel;

    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= &Stokes.vel_idx;
    IdxDescCL* pidx1= &Stokes.pr_idx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &Stokes.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &Stokes.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;

    MatDescCL* A= &Stokes.A;
    MatDescCL* B= &Stokes.B;

    Uint step= 0;
    StokesDoerflerMarkCL<typename MyStokesCL::est_fun, MyStokesCL>
        Estimator(rel_red, markratio, .484473073129685, true, &MyStokesCL::ResidualErrEstimator, Stokes);
    bool new_marks= false;

    vidx1->Set(3, 3, 0, 0);
    vidx2->Set(3, 3, 0, 0);
    pidx1->Set(1, 0, 0, 0);
    pidx2->Set(1, 0, 0, 0);
    TimerCL time;

    do
    {
        MG.Refine();
        Stokes.CreateNumberingVel(MG.GetLastLevel(), vidx1);    
        Stokes.CreateNumberingPr(MG.GetLastLevel(), pidx1);    
        std::cerr << "altes und neues TriangLevel: " << vidx2->TriangLevel << ", "
                  << vidx1->TriangLevel << std::endl;
        MG.SizeInfo(std::cerr);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        std::cerr << "Anzahl der Druck-Unbekannten: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cerr << "Anzahl der Geschwindigkeitsunbekannten: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;
        if (p2->RowIdx)
        {
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &PrBndData, &MG);
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &PrBndData, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &VelBndData, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &VelBndData, &MG);
//            Interpolate(pr1, pr2);
//            Interpolate(vel1, vel2);
//            Stokes.CheckSolution(v1,p1,&LsgVel, &DLsgVel, &LsgPr);
            v1->Clear();
            p1->Clear();
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        time.Reset();
        time.Start();
        Stokes.SetupSystem(A, b, B, c);
        time.Stop();
        std::cerr << "SetupSystem: " << time.GetTime() << " seconds." << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cerr << " A*x: " << time.GetTime() << " seconds." << std::endl;
        time.Reset();
        time.Start();
        transp_mul( A->Data, v1->Data);
        time.Stop();
        std::cerr << "AT*x: " << time.GetTime() << " seconds." << std::endl;
        
//        { // write system in files for MatLab
//            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
//            Adat << A->Data;   Bdat << B->Data;    bdat << b->Data;    cdat << c->Data;
//        }
        Stokes.GetDiscError(&LsgVel, &LsgPr);
//std::cout << A->Data << std::endl << b->Data << std::endl
//          << B->Data << std::endl << c->Data << std::endl
//          << v1->Data << std::endl << p1->Data << std::endl;
        time.Reset();

        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        Stokes.SetupMass( &M);
        double err0= (A->Data*v1->Data + transp_mul( B->Data, p1->Data) - b->Data).norm2()
                    +(B->Data*v1->Data - c->Data).norm2();
        std::cerr << "000 residual: " << std::sqrt( err0) << std::endl;

        double outer_tol= tol;
        switch (meth) {
          case 1: { // Schur
            PSchur_PCG_CL schurSolver( M.Data, 200, outer_tol*std::sqrt( err0), 200, inner_iter_tol);
            time.Start();
//std::cout << M.Data << std::endl;
//std::cout << A->Data << std::endl << b->Data << std::endl
//          << B->Data << std::endl << c->Data << std::endl
//          << v1->Data << std::endl << p1->Data << std::endl;
            schurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            double err= (A->Data*v1->Data + transp_mul( B->Data, p1->Data) - b->Data).norm2()
                        +(B->Data*v1->Data - c->Data).norm2();
            std::cerr << "000 residual: " << std::sqrt( err)/std::sqrt( err0) << std::endl;
            break;
          }
          case 0: { // Uzawa
//            double tau;
//            Uint inner_iter;
//            std::cerr << "tau = "; std::cin >> tau;
//            std::cerr << "#PCG steps = "; std::cin >> inner_iter;
            Uzawa_PCG_CL uzawaSolver( M.Data, 5000, outer_tol, uzawa_inner_iter, inner_iter_tol, tau);
            time.Start();
            uzawaSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "iterations: " << uzawaSolver.GetIter()
                      << "\tresidual: " << uzawaSolver.GetResid() << std::endl;
            break;
          }
          case 2: { // Minres
            std::cerr << "MINRES!\n";
            StokesVectorCL x( v1->Data, p1->Data);
            StokesVectorCL rhs( b->Data, c->Data);
            StokesMatrixCL K( A->Data, B->Data);
//            CheckMatVec( K, A->Data, B->Data, x, v1->Data, p1->Data);
            MResSolverCL mressolver( uzawa_inner_iter, outer_tol*std::sqrt( err0));
            time.Start();
            mressolver.Solve( K, x, rhs);
            time.Stop();
            std::copy( &x.raw()[0], &x.raw()[0] + v1->Data.size(), Addr( v1->Data.raw()));
            std::copy( &x.raw()[0] + v1->Data.size(), &x.raw()[0] + x.size(), Addr( p1->Data.raw()));
            std::cerr << "iterations: " << mressolver.GetIter()
                      << "\terror: " << mressolver.GetResid()/std::sqrt( err0) << std::endl;
            break;
          }
          case 3: { // PMinres
            std::cerr << "PMINRES!\n";
            StokesVectorCL x( v1->Data, p1->Data);
            StokesVectorCL rhs( b->Data, c->Data);
            StokesMatrixCL K( A->Data, B->Data);
//            CheckMatVec( K, A->Data, B->Data, x, v1->Data, p1->Data);
            FullPreCL pc( M.Data);
            PLanczosONBCL<StokesMatrixCL, StokesVectorCL, FullPreCL> q( K, pc, rhs - K*x);
            PMResSolverCL<PLanczosONBCL<StokesMatrixCL, StokesVectorCL, FullPreCL> > pmr( q, uzawa_inner_iter, outer_tol*std::sqrt( err0));
//            SSORPreCL pc( M.Data);
//            PLanczosONBCL<StokesMatrixCL, StokesVectorCL, SSORPreCL> q( K, pc, rhs - K*x);
//            PMResSolverCL<PLanczosONBCL<StokesMatrixCL, StokesVectorCL, SSORPreCL> > pmr( q, uzawa_inner_iter, outer_tol*std::sqrt( err0));
//            DummyPcCL pc;
//            PLanczosONBCL<StokesMatrixCL, StokesVectorCL, DummyPcCL> q( K, pc, rhs - K*x);
//            PMResSolverCL<PLanczosONBCL<StokesMatrixCL, StokesVectorCL, DummyPcCL> > pmr( q, uzawa_inner_iter, outer_tol*std::sqrt( err0));
            time.Start();
            pmr.Solve( K, x, rhs);
            time.Stop();
            std::copy( &x.raw()[0], &x.raw()[0] + v1->Data.size(), Addr( v1->Data.raw()));
            std::copy( &x.raw()[0] + v1->Data.size(), &x.raw()[0] + x.size(), Addr( p1->Data.raw()));
            std::cerr << "iterations: " << pmr.GetIter()
                      << "\terror: " << pmr.GetResid()/std::sqrt( err0) << std::endl;
            break;
          }
          case 4: {	
            std::cerr << "MG_Schur!\n";
            MGDataCL MGData;
	    IdxDescCL* c_idx;
            time.Reset();
            time.Start();
            for(Uint lvl= 0; lvl<=MG.GetLastLevel(); ++lvl) {
                MGData.push_back( MGLevelDataCL());
                MGLevelDataCL& tmp= MGData.back();
                std::cerr << "                        Create MGData on Level " << lvl << std::endl;
                tmp.Idx.Set( 3, 3);
                Stokes.CreateNumberingVel(lvl, &tmp.Idx);
                tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
                std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
                Stokes.SetupStiffnessMatrix( &tmp.A );
                if(lvl!=0) {
                    std::cerr << "                       Create Prolongation on Level " << lvl << std::endl;
                    SetupP2ProlongationMatrix( MG, tmp.P, c_idx, &tmp.Idx);
//                   std::cout << "    Matrix P " << tmp.P.Data << std::endl;
	        }
                c_idx= &tmp.Idx;
            }
            time.Stop();
            std::cerr <<"MG-Setup: " << time.GetTime() << " seconds." << std::endl;
	    MGDataCL& MGD = MGData;
            std::cerr << "Check MG-Data..." << std::endl;
            std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
            std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
//            CheckMGData( MGData.begin(), MGData.end());
            PSchur_MG_CL MGschurSolver( M.Data, 200, outer_tol*std::sqrt( err0), MGD, 200, inner_iter_tol);
            time.Start();
            MGschurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            double err= (A->Data*v1->Data + transp_mul( B->Data, p1->Data) - b->Data).norm2()
                        +(B->Data*v1->Data - c->Data).norm2();
            std::cerr << "000 residual: " << std::sqrt( err)/std::sqrt( err0) << std::endl;
            break;
          }
          case 5: { // MG-PMinres
            std::cerr << "MG-PMINRES!\n";
            StokesVectorCL x( v1->Data, p1->Data);
            StokesVectorCL rhs( b->Data, c->Data);
            StokesMatrixCL K( A->Data, B->Data);
//            CheckMatVec( K, A->Data, B->Data, x, v1->Data, p1->Data);
            MGDataCL MGData;
	    IdxDescCL* c_idx;
            time.Reset();
            time.Start();
            for(Uint lvl= 0; lvl<=MG.GetLastLevel(); ++lvl) {
                MGData.push_back( MGLevelDataCL());
                MGLevelDataCL& tmp= MGData.back();
                std::cerr << "                        Create MGData on Level " << lvl << std::endl;
                tmp.Idx.Set( 3, 3);
                Stokes.CreateNumberingVel(lvl, &tmp.Idx);
                tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
                std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
                Stokes.SetupStiffnessMatrix( &tmp.A );
                if(lvl!=0) {
                    std::cerr << "                       Create Prolongation on Level " << lvl << std::endl;
                    SetupP2ProlongationMatrix( MG, tmp.P, c_idx, &tmp.Idx);
//                   std::cout << "    Matrix P " << tmp.P.Data << std::endl;
	        }
                c_idx= &tmp.Idx;
            }
            time.Stop();
            std::cerr <<"MG-Setup: " << time.GetTime() << " seconds." << std::endl;
            std::cerr << "Check MG-Data..." << std::endl;
            std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
            std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
//            CheckMGData( MGData.begin(), MGData.end());
            MGSolverCL MGAsolver( MGData, 2, 1e-10);
//            MGPreCL pc( MGAsolver, M.Data);
            MGPreCL pc( MGData, M.Data);
            PLanczosONBCL<StokesMatrixCL, StokesVectorCL, MGPreCL> q( K, pc, rhs - K*x);
            PMResSolverCL<PLanczosONBCL<StokesMatrixCL, StokesVectorCL, MGPreCL> > pmr( q, uzawa_inner_iter, outer_tol*std::sqrt( err0));
            time.Start();
            pmr.Solve( K, x, rhs);
            time.Stop();
            std::copy( &x.raw()[0], &x.raw()[0] + v1->Data.size(), Addr( v1->Data.raw()));
            std::copy( &x.raw()[0] + v1->Data.size(), &x.raw()[0] + x.size(), Addr( p1->Data.raw()));
            std::cerr << "iterations: " << pmr.GetIter()
                      << "\terror: " << pmr.GetResid()/std::sqrt( err0) << std::endl;
            break;
          }
          case 6: { // MG-Uzawa
            std::cerr << "MG_Uzawa!\n";
            MGDataCL MGData;
	    IdxDescCL* c_idx;
            time.Reset();
            time.Start();
            for(Uint lvl= 0; lvl<=MG.GetLastLevel(); ++lvl) {
                MGData.push_back( MGLevelDataCL());
                MGLevelDataCL& tmp= MGData.back();
                std::cerr << "                        Create MGData on Level " << lvl << std::endl;
                tmp.Idx.Set( 3, 3);
                Stokes.CreateNumberingVel(lvl, &tmp.Idx);
                tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
                std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
                Stokes.SetupStiffnessMatrix( &tmp.A );
                if(lvl!=0) {
                    std::cerr << "                       Create Prolongation on Level " << lvl << std::endl;
                    SetupP2ProlongationMatrix( MG, tmp.P, c_idx, &tmp.Idx);
//                   std::cout << "    Matrix P " << tmp.P.Data << std::endl;
	        }
                c_idx= &tmp.Idx;
            }
            time.Stop();
            std::cerr <<"MG-Setup: " << time.GetTime() << " seconds." << std::endl;
	    MGDataCL& MGD = MGData;
            std::cerr << "Check MG-Data..." << std::endl;
            std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
            std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
//            CheckMGData( MGData.begin(), MGData.end());
            Uzawa_MG_CL uzawaMG( M.Data, 5000, outer_tol*std::sqrt( err0), MGD, 1, inner_iter_tol, tau);
            time.Start();
            uzawaMG.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "iterations: " << uzawaMG.GetIter()
                      << "\tresidual: " << uzawaMG.GetResid() << std::endl;
            break;
          }
        }
        std::cerr << "Solver: "<<time.GetTime()<<" seconds.\n";
        Stokes.CheckSolution(v1, p1, &LsgVel, &DLsgVel, &LsgPr);
        if (step==0) {
            Estimator.Init(typename MyStokesCL::DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::DiscVelSolCL(v1, &VelBndData, &MG));
        }
        time.Reset();
        time.Start();
        char dummy;
        std::cin >> dummy;
        new_marks= Estimator.Estimate(typename MyStokesCL::DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::DiscVelSolCL(v1, &VelBndData, &MG) );
        time.Stop();
        std::cerr << "Estimation: " << time.GetTime() << " seconds.\n";
        A->Reset();
        B->Reset();
        b->Reset();
        c->Reset();
//        std::cerr << "Loesung Druck: " << p1->Data << std::endl;
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cerr << std::endl;
    }
    while (++step<maxStep);
    // we want the solution to be in Stokes.v, Stokes.pr
    if (v2 == &loc_v)
    {
        Stokes.vel_idx.swap( loc_vidx);
        Stokes.pr_idx.swap( loc_pidx);
        Stokes.v.SetIdx(&Stokes.vel_idx);
        Stokes.p.SetIdx(&Stokes.pr_idx);
        
        Stokes.v.Data= loc_v.Data;
        Stokes.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    std::cout << sizeof( DROPS::VectorCL) << std::endl;
    std::cout << sizeof( DROPS::StokesVectorCL) << std::endl;
    if (argc!=10)
    {
        std::cerr << "Usage: sdropsP2 <omega> <inner_iter_tol> <tol> <meth> <num_refinement> <rel_red> <markratio> <tau> <uz_inner_iter>"
                  << std::endl;
        return 1;
    }
    double omega= atof(argv[1]);
    double inner_iter_tol= atof(argv[2]);
    double tol= atof(argv[3]);
    int meth= atoi(argv[4]);
    int num_ref= atoi(argv[5]);
    double rel_red= atof(argv[6]);
    double markratio= atof(argv[7]);
    double tau= atof(argv[8]);
    unsigned int uz_inner_iter= atoi(argv[9]);
    std::cerr << "Omega: " << omega << ", "
              << "inner iter tol: " << inner_iter_tol << ", "
              << "tol: " << tol << ", "
              << "meth: " << meth << ", "
              << "refinements: " << num_ref << ", "
              << "relative error reduction: " << rel_red << ", "
              << "markratio: " << markratio << ", "
              << "tau: " << tau << ", "
              << "uzawa inner iter: " << uz_inner_iter
              << std::endl;

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= M_PI/4.;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 3, 3, 3);
    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel};
        
    StokesOnBrickCL prob(brick, StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;

    Strategy(prob, omega, inner_iter_tol, tol, meth, num_ref, rel_red, markratio, tau, uz_inner_iter);
    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    fil << DROPS::GeomSolOutCL<MyStokesCL::DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;

//    std::cout << DROPS::GeomMGOutCL(mg, -1, true) << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
