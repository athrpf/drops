#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
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
            P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>  pr2(p2, &PrBndData, &MG);
            P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>        pr1(p1, &PrBndData, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &VelBndData, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &VelBndData, &MG);
            Interpolate(pr1, pr2);
//            Interpolate(vel1, vel2);
//            Stokes.CheckSolution(v1,p1,&LsgVel, &DLsgVel, &LsgPr);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        time.Reset();
        time.Start();
        Stokes.SetupSystem(A, b, B, c);
        time.Stop();
        std::cerr << time.GetTime() << " seconds for setting up all systems!" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cerr << " A*x took " << time.GetTime() << " seconds!" << std::endl;
        time.Reset();
        time.Start();
        transp_mul( A->Data, v1->Data);
        time.Stop();
        std::cerr << "AT*x took " << time.GetTime() << " seconds!" << std::endl;
        
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

        double outer_tol= tol;
        if (meth)
        {
            PSchur_PCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
            time.Start();
//std::cout << M.Data << std::endl;
//std::cout << A->Data << std::endl << b->Data << std::endl
//          << B->Data << std::endl << c->Data << std::endl
//          << v1->Data << std::endl << p1->Data << std::endl;
            schurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
        }
        else // Uzawa
        {
//            double tau;
//            Uint inner_iter;
//            std::cerr << "tau = "; std::cin >> tau;
//            std::cerr << "#PCG steps = "; std::cin >> inner_iter;
            Uzawa_PCG_CL uzawaSolver( M.Data, 5000, outer_tol, uzawa_inner_iter, inner_iter_tol, tau);
            time.Start();
            uzawaSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "Iterationen: " << uzawaSolver.GetIter()
                      << "Norm des Res.: " << uzawaSolver.GetResid() << std::endl;
        }
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
        Stokes.CheckSolution(v1, p1, &LsgVel, &DLsgVel, &LsgPr);
        if (step==0)
        {
            Estimator.Init(typename MyStokesCL::DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::DiscVelSolCL(v1, &VelBndData, &MG));
        }
        time.Reset();
        time.Start();
    char dummy;
    std::cin >> dummy;
        new_marks= Estimator.Estimate(typename MyStokesCL::DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::DiscVelSolCL(v1, &VelBndData, &MG) );
        time.Stop();
        std::cerr << "Estimation took " << time.GetTime() << " seconds\n";
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
