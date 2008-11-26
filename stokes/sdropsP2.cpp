#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include "num/MGsolver.h"
#include <fstream>


inline DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    return ret;
}

// Jacobi-matrix od exact solution
inline DROPS::SMatrixCL<3, 3> DLsgVel(const DROPS::Point3DCL& p)
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


// Volume of the box: 0.484473073129685
// int(p)/vol = -0.125208551608365
inline double LsgPr(const DROPS::Point3DCL& p, double)
{
    return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]) - 0.125208551608365;
}
inline double LsgPr(const DROPS::Point3DCL& p)
{
    return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]) - 0.125208551608365;
}


// q*u - nu*laplace u + Dp = f
//                  -div u = 0
class StokesCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double)
        { DROPS::SVectorCL<3> ret(0.0); ret[2]= 3.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]); return ret; }
    const double nu;

    StokesCoeffCL() : nu(1.0) {}
};

typedef DROPS::StokesP2P1CL<StokesCoeffCL>
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

namespace DROPS // for Strategy
{


using ::MyStokesCL;

class PSchur_PCG_CL: public PSchurSolverCL<PCG_SsorCL>
{
  private:
    PCG_SsorCL _PCGsolver;
  public:
    PSchur_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double omega= 1.)
        : PSchurSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol),
          _PCGsolver(SSORPcCL(omega), inner_iter, inner_tol)
        {}
};

class PSchur_MG_CL: public PSchurSolverCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> >
{
  private:
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> _MGsolver;
    SSORsmoothCL smoother_;
    PCG_SsorCL   solver_;
  public:
    PSchur_MG_CL( MatrixCL& M,      int outer_iter, double outer_tol,
                  int inner_iter, double inner_tol )
        : PSchurSolverCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> >( _MGsolver, M, outer_iter, outer_tol ),
          _MGsolver( smoother_, solver_, inner_iter, inner_tol ),
          smoother_(1.0), solver_(SSORPcCL(1.0), 500, inner_tol)
        {}
     MLMatrixCL* GetPVel() { return _MGsolver.GetProlongation(); }
};

class Uzawa_PCG_CL : public UzawaSolverCL<PCG_SsorCL>
{
  private:
    PCG_SsorCL _PCGsolver;
  public:
    Uzawa_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1., double omega=1.)
        : UzawaSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol, tau),
          _PCGsolver(SSORPcCL(omega), inner_iter, inner_tol)
        {}
};

class Uzawa_MG_CL : public UzawaSolver2CL<PCG_SsorCL, MGSolverCL<SSORsmoothCL, PCG_SsorCL> >
{
  private:
    PCG_SsorCL PCGsolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> MGsolver_;
    SSORsmoothCL smoother_;
    PCG_SsorCL   solver_;

  public:
    Uzawa_MG_CL(MatrixCL& M,      int outer_iter, double outer_tol,
                int inner_iter, double inner_tol, double tau= 1., double omega= 1.)
        : UzawaSolver2CL<PCG_SsorCL, MGSolverCL<SSORsmoothCL, PCG_SsorCL> >( PCGsolver_, MGsolver_, M,
                                                  outer_iter, outer_tol, tau),
          PCGsolver_( SSORPcCL(omega), inner_iter, inner_tol),
          MGsolver_( smoother_, solver_, inner_iter, inner_tol ),
          smoother_(1.), solver_(SSORPcCL(1.0), 500, inner_tol)
        {}
    MLMatrixCL* GetPVel() { return MGsolver_.GetProlongation(); }
};

typedef PMResSolverCL<LanczosONBCL<VectorCL> > MinresSPT;
class MinresSPCL : public BlockMatrixSolverCL<MinresSPT>
{
  private:
    MinresSPT solver_;
    LanczosONBCL<VectorCL> q_;

  public:
    MinresSPCL(int maxiter, double tol)
        :BlockMatrixSolverCL<MinresSPT>( solver_), solver_(q_, maxiter, tol),
         q_()
    {}
};

typedef SolverAsPreCL<PCG_SsorCL> PPcT;
typedef BlockPreCL<PPcT, ISPreCL> BlockDiagPCGPreCL;
typedef PMResSolverCL<PLanczosONBCL<VectorCL, BlockDiagPCGPreCL> > PMinresSP_DiagPCGT;
class PMinresSP_DiagPCG_CL : public BlockMatrixSolverCL<PMinresSP_DiagPCGT>
{
  private:
    PMinresSP_DiagPCGT solver_;
    PCG_SsorCL PPA_;
    MatrixCL& M_;
    PPcT PA_;
    ISPreCL PS_;
    BlockDiagPCGPreCL pre_;
    PLanczosONBCL<VectorCL, BlockDiagPCGPreCL> q_;

  public:
    PMinresSP_DiagPCG_CL(MatrixCL& M, int maxiter, double tol, double omega= 1.)
        :BlockMatrixSolverCL<PMinresSP_DiagPCGT>( solver_), solver_( q_, maxiter, tol),
         PPA_(SSORPcCL(omega), 8, 1e-20), M_(M),
         PA_(PPA_), PS_(M_, M_, 0.0, 1.0, omega),
         pre_( PA_, PS_), q_( pre_)
    {}
};

template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, double omega, double inner_iter_tol, double tol, int meth,
                                           Uint maxStep, double rel_red, double markratio,
                                           double tau, Uint uzawa_inner_iter)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();
    const typename MyStokesCL::BndDataCL::PrBndDataCL& PrBndData= Stokes.GetBndData().Pr;
    const typename MyStokesCL::BndDataCL::VelBndDataCL& VelBndData= Stokes.GetBndData().Vel;

    MLIdxDescCL  loc_vidx, loc_pidx;
    MLIdxDescCL* vidx1= &Stokes.vel_idx;
    MLIdxDescCL* pidx1= &Stokes.pr_idx;
    MLIdxDescCL* vidx2= &loc_vidx;
    MLIdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &Stokes.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &Stokes.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;

    MLMatDescCL* A= &Stokes.A;
    MLMatDescCL* B= &Stokes.B;

    Uint step= 0;
    StokesDoerflerMarkCL<typename MyStokesCL::est_fun, MyStokesCL>
        Estimator(rel_red, markratio, .484473073129685, true, &MyStokesCL::ResidualErrEstimator, Stokes);
    bool new_marks= false;

    vidx1->SetFE( vecP2_FE);
    vidx2->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    pidx2->SetFE( P1_FE);
    TimerCL time;

    do
    {
        MG.Refine();
        if (meth >=4) // mgused
        {
            match_fun match= MG.GetBnd().GetMatchFun();
            vidx1->resize( MG.GetNumLevel(), vecP2_FE, Stokes.GetBndData().Vel, match);
            pidx1->resize( MG.GetNumLevel(), P1_FE, Stokes.GetBndData().Pr, match);
            A->Data.resize( MG.GetNumLevel());
            B->Data.resize( MG.GetNumLevel());
        }
            
        Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx1);
        Stokes.CreateNumberingPr ( MG.GetLastLevel(), pidx1);
        std::cerr << "altes und neues TriangLevel: " << vidx2->TriangLevel() << ", "
                  << vidx1->TriangLevel() << std::endl;
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
//            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &PrBndData, &MG);
//            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &PrBndData, &MG);
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

        MLMatDescCL M;
        if ( meth>=4 ) M.Data.resize(pidx1->size());
        M.SetIdx( pidx1, pidx1);
        Stokes.SetupPrMass( &M);
        double err0= norm_sq( A->Data*v1->Data + transp_mul( B->Data, p1->Data) - b->Data)
                    +norm_sq( B->Data*v1->Data - c->Data);
        std::cerr << "000 residual: " << std::sqrt( err0) << std::endl;

        double outer_tol= tol;
        switch (meth) {
          case 1: { // Schur
            PSchur_PCG_CL schurSolver( M.Data.GetFinest(), 200, outer_tol*std::sqrt( err0), 200, inner_iter_tol, omega);
            time.Start();
//std::cout << M.Data << std::endl;
//std::cout << A->Data << std::endl << b->Data << std::endl
//          << B->Data << std::endl << c->Data << std::endl
//          << v1->Data << std::endl << p1->Data << std::endl;
            schurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            double err= norm_sq( A->Data*v1->Data + transp_mul( B->Data, p1->Data) - b->Data)
                        +norm_sq( B->Data*v1->Data - c->Data);
            std::cerr << "000 residual: " << std::sqrt( err)/std::sqrt( err0) << std::endl;
            break;
          }
          case 0: { // Uzawa
//            double tau;
//            Uint inner_iter;
//            std::cerr << "tau = "; std::cin >> tau;
//            std::cerr << "#PCG steps = "; std::cin >> inner_iter;
            Uzawa_PCG_CL uzawaSolver( M.Data.GetFinest(), 5000, outer_tol, uzawa_inner_iter, inner_iter_tol, tau, omega);
            time.Start();
            uzawaSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "iterations: " << uzawaSolver.GetIter()
                      << "\tresidual: " << uzawaSolver.GetResid() << std::endl;
            break;
          }
          case 2: { // Stokes-Minres
            std::cerr << "Stokes-Minres!\n";
            MinresSPCL solver( uzawa_inner_iter, outer_tol/**std::sqrt( err0)*/);
            time.Start();
            solver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "iterations: " <<solver.GetIter()
                      << "\tresidual: " << solver.GetResid()/**std::sqrt( err0)*/ << std::endl;
            break;
          }
          case 3: { // Stokes-PMinres
            std::cerr << "Stokes-PMinres!\n";
            PMinresSP_DiagPCG_CL solver( M.Data.GetFinest(), uzawa_inner_iter, outer_tol/**std::sqrt( err0)*/);
            time.Start();
            solver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "iterations: " <<solver.GetIter()
                      << "\tresidual: " << solver.GetResid()/**std::sqrt( err0)*/ << std::endl;
            break;
          }
          case 4: {
            std::cerr << "MG_Schur!\n";
            PSchur_MG_CL MGschurSolver( M.Data.GetFinest(), 200, outer_tol*std::sqrt( err0), 200, inner_iter_tol);
            MLMatrixCL* PVel = MGschurSolver.GetPVel();
            SetupP2ProlongationMatrix( MG, *PVel, vidx1, vidx1);
            std::cerr << "Check MG-Data..." << std::endl;
            std::cerr << "                begin     " << vidx1->GetCoarsest().NumUnknowns() << std::endl;
            std::cerr << "                end       " << vidx1->GetFinest().NumUnknowns() << std::endl;
            CheckMGData( Stokes.A.Data, *PVel);

            time.Start();
            MGschurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            double err= norm_sq( A->Data*v1->Data + transp_mul( B->Data, p1->Data) - b->Data)
                        +norm_sq( B->Data*v1->Data - c->Data);
            std::cerr << "000 residual: " << std::sqrt( err)/std::sqrt( err0) << std::endl;
            break;
          }
          case 6: { // MG-Uzawa
            std::cerr << "MG_Uzawa!\n";
            Uzawa_MG_CL uzawaMG( M.Data.GetFinest(), 5000, outer_tol*std::sqrt( err0), 1, inner_iter_tol, tau, omega);
            MLMatrixCL* PVel = uzawaMG.GetPVel();
            SetupP2ProlongationMatrix( MG, *PVel, vidx1, vidx1);
            std::cerr << "Check MG-Data..." << std::endl;
            std::cerr << "                begin     " << vidx1->GetCoarsest().NumUnknowns() << std::endl;
            std::cerr << "                end       " << vidx1->GetFinest().NumUnknowns() << std::endl;
            CheckMGData( Stokes.A.Data, *PVel);

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
            Estimator.Init(typename MyStokesCL::const_DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::const_DiscVelSolCL(v1, &VelBndData, &MG));
        }
        time.Reset();
        time.Start();
        char dummy;
        std::cin >> dummy;
        new_marks= Estimator.Estimate(typename MyStokesCL::const_DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::const_DiscVelSolCL(v1, &VelBndData, &MG) );
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
    if (argc!=10)
    {
        std::cerr << "Usage: sdropsP2 <omega> <inner_iter_tol> <tol> <meth> <num_refinement> <rel_red> <markratio> <tau> <uz_inner_iter>"
                  << std::endl;
        return 1;
    }
    double omega= std::atof(argv[1]);
    double inner_iter_tol= std::atof(argv[2]);
    double tol= std::atof(argv[3]);
    int meth= std::atoi(argv[4]);
    int num_ref= std::atoi(argv[5]);
    double rel_red= std::atof(argv[6]);
    double markratio= std::atof(argv[7]);
    double tau= std::atof(argv[8]);
    unsigned int uz_inner_iter= std::atoi(argv[9]);
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
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
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
    fil << DROPS::GeomSolOutCL<MyStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;

//    std::cout << DROPS::GeomMGOutCL(mg, -1, true) << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
