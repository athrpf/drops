#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include <fstream>

// q*u - nu*laplace u + Dp = f
//                  -div u = 0
class DrivenCavityCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL&)
        { DROPS::SVectorCL<3> ret(0.0); return ret; }
    const double nu;
    
    DrivenCavityCL() : nu(1.0) {}
};

typedef DROPS::StokesP2P1CL<DrivenCavityCL> 
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

inline DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)   { return DROPS::SVectorCL<3>(0.); }
//inline DROPS::SMatrixCL<3, 3> NullMat( const DROPS::Point3DCL&)   { return DROPS::SMatrixCL<3, 3>(0.); }
inline double Nullsc( const DROPS::Point3DCL&)   { return 0.; }

const double st= 0.1;

inline DROPS::SVectorCL<3> Stroem( const DROPS::Point3DCL& p, double)
{
    const DROPS::SVectorCL<3> ret= DROPS::std_basis<3>(1);
    const double d0= fabs(p[0]-.5);
    const double d1= fabs(p[1]-.5);
    const double m= std::max(d0, d1);
    return (.5-st<m) ? ((.5-m)/st)*ret : ret;
}


void MarkLower( DROPS::MultiGridCL& mg, double tresh)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin()),
             ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
    {
        if (GetBaryCenter(*It)[2]<=tresh )
            It->SetRegRefMark();
    }
}

namespace DROPS // for Strategy
{

using ::MyStokesCL;


template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, double inner_iter_tol, double tol,
                                   int meth, Uint maxStep, double rel_red, double markratio,
                                   double tau, Uint uzawa_inner_iter)
// flow control
{
    typedef StokesP2P1CL<Coeff> StokesCL;
    
    MultiGridCL& MG= Stokes.GetMG();
    const typename MyStokesCL::BndDataCL::PrBndDataCL& PrBndData= Stokes.GetBndData().Pr;
    const typename MyStokesCL::BndDataCL::VelBndDataCL& VelBndData= Stokes.GetBndData().Vel;

    IdxDescCL  loc_vidx, loc_pidx, corr_vidx;
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
        Estimator(rel_red, markratio, 1., true, &MyStokesCL::ResidualErrEstimator, Stokes);
    bool new_marks= false;

    vidx1->Set( 3, 3, 0, 0);
    vidx2->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    pidx2->Set( 1, 0, 0, 0);
    
    corr_vidx.Set( 3, 3, 0, 0);

//    MarkLower(MG,0.25); 
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
/*
            P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>  pr2(p2, &PrBndData, &MG);
            P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>        pr1(p1, &PrBndData, &MG);
            Interpolate(pr1, pr2);

            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &VelBndData, &MG);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &VelBndData, &MG);
            Stokes.CreateNumberingVel( v2->GetLevel(), &corr_vidx);
            VelVecDescCL corr_v;
            corr_v.SetIdx( &corr_vidx);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>  corr_vel( &corr_v, &VelBndData, &MG);
            Adapt( corr_vel, vel2);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL>  corr_vel_read( &corr_v, &VelBndData, &MG);
            Interpolate( vel1, corr_vel_read );   */
//            CheckSolution(v1,p1,&LsgVel,&LsgPr);
//            Stokes.CheckSolution(v1,p1,&Null,&Nullsc);
            v2->Reset();
            p2->Reset();
        }
        A->Reset();
        B->Reset();
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
/*        
        { // write system in files for MatLab
            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
            Adat << A->Data;   Bdat << B->Data;    bdat << b->Data;    cdat << c->Data;
        }
*/ //        Stokes.GetDiscError(&LsgVel, &LsgPr);
//std::cout << A->Data << std::endl << b->Data << std::endl;
/*        double half= M_PI/8;
        MultiGridCL::TriangVertexIteratorCL vert= MG.GetTriangVertexBegin(A->GetRowLevel());
        while (vert->GetCoord()[0]!=half || vert->GetCoord()[1]!=half || vert->GetCoord()[2]!=half) ++vert;
        IdxT unk= vert->Unknowns(A->RowIdx->Idx);
        std::cerr << vert->GetCoord() << " has index " << unk << std::endl;
        std::cerr << "A(i,i) = " << A->Data(unk,unk) <<std::endl;    
        std::cerr << "B(i,j) = " << B->Data(vert->Unknowns(B->RowIdx->Idx),unk) << std::endl;
*/
//        Uint meth;
//        std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
        time.Reset();

        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        Stokes.SetupMass( &M);

        double outer_tol= tol;
//        std::cerr << "tol = "; std::cin >> outer_tol;

        if (meth)
        {
//            PSchur_PCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
            PSchur_GSPCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
            time.Start();
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
                      << "\tNorm des Res.: " << uzawaSolver.GetResid() << std::endl;
        }
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
//        Stokes.CheckSolution(v1, p1, &LsgVel, &LsgPr);
        if (step==0)
        {
            Estimator.Init(typename MyStokesCL::DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::DiscVelSolCL(v1, &VelBndData, &MG));
        }
        time.Reset();
        time.Start();
//    char dummy;
//    std::cin >> dummy;
        new_marks= Estimator.Estimate(typename MyStokesCL::DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::DiscVelSolCL(v1, &VelBndData, &MG) );
        time.Stop();
        std::cerr << "Estimation took " << time.GetTime() << " seconds\n";
        A->Reset();
        B->Reset();
        b->Reset();
        c->Reset();
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cerr << std::endl;
    }
    while (new_marks && ++step<maxStep);
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
    if (argc!=9)
    {
        std::cerr << "Usage: drvcavad <inner_iter_tol> <tol> <meth> <num_refinement> <rel_red> <markratio> <tau> <uz_inner_iter>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Stroem};
        
    StokesOnBrickCL prob(brick, DrivenCavityCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;
    double inner_iter_tol= atof(argv[1]);
    double tol= atof(argv[2]);
    int meth= atoi(argv[3]);
    int num_ref= atoi(argv[4]);
    double rel_red= atof(argv[5]);
    double markratio= atof(argv[6]);
    double tau= atof(argv[7]);
    unsigned int uz_inner_iter= atoi(argv[8]);
    std::cerr << "inner iter tol: " << inner_iter_tol << ", ";
    std::cerr << "tol: " << tol << ", ";
    std::cerr << "meth: " << meth << ", ";
    std::cerr << "refinements: " << num_ref << ", ";
    std::cerr << "relative error reduction: " << rel_red << ", ";
    std::cerr << "markratio: " << markratio << ", ";
    std::cerr << "tau: " << tau << ", ";
    std::cerr << "uzawa inner iter: " << uz_inner_iter << std::endl;
    Strategy(prob, inner_iter_tol, tol, meth, num_ref, rel_red, markratio, tau, uz_inner_iter);
//    std::cerr << "hallo" << std::endl;
//    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;
    fil << DROPS::GeomSolOutCL<MyStokesCL::DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, -10, 10) << std::endl;

    DROPS::IdxDescCL tecIdx;
    tecIdx.Set( 1, 0, 0, 0);
    prob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);    
    
    std::ofstream v2d("data2D.dat");
    DROPS::TecPlot2DSolOutCL< MyStokesCL::DiscVelSolCL, MyStokesCL::DiscPrSolCL>
        tecplot2d( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx, -1, 1, 0.5); // cutplane is y=0.5
    v2d << tecplot2d;
    v2d.close();

    std::ofstream v3d("data3D.dat");
    DROPS::MapleSolOutCL<MyStokesCL::DiscVelSolCL, MyStokesCL::DiscPrSolCL>
        mapleplot( mg, prob.GetVelSolution(), prob.GetPrSolution(), -1, DROPS::PlaneCL(DROPS::std_basis<3>(2), .5));
    v3d << mapleplot;
    v3d.close();

    v3d.open("data3D_tec.dat");
    DROPS::TecPlotSolOutCL< MyStokesCL::DiscVelSolCL, MyStokesCL::DiscPrSolCL>
        tecplot( mg, prob.GetVelSolution(), prob.GetPrSolution() );
    v3d << tecplot;
    v3d.close();

//    std::ofstream g3d("gdata3D.dat");
//    DROPS::MapleMGOutCL mgplot(mg);
//    g3d << mgplot;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
