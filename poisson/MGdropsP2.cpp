#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "poisson/poisson.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include "num/fe.h"
#include "misc/utils.h"
#include <fstream>

// laplace u + q*u = f
class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&, double= 0.0) { return 0.0; }
    static double f(const DROPS::Point3DCL& p, double= 0.0)
        { return 128.0*( p[0]*p[1]*(1-p[0])*(1-p[1]) + p[0]*p[2]*(1-p[0])*(1-p[2])
                                                     + p[1]*p[2]*(1-p[1])*(1-p[2]) ); }
//    static double f(const Point3DCL& p, double= 0.0) { return p[2]>0.49?-15.:0; }
//    static double f(const DROPS::Point3DCL& p, double= 0.0) { return 0.0; }
};


inline double Lsg( const DROPS::Point3DCL& p, double= 0.0)
{
    return 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
}

namespace DROPS // for strategy
{

template<class Coeff>
void Strategy(PoissonP2CL<Coeff>& Poisson, double omega)
{
    typedef PoissonP2CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    TimerCL time;

    MG.SizeInfo(std::cout);

    // Initialize the idx, A, P
    time.Reset();
    Poisson.SetNumLvl( MG.GetNumLevel());
    Poisson.CreateNumbering( MG.GetLastLevel(), &Poisson.idx);
    MLIdxDescCL* idx = &Poisson.idx;
    Poisson.A.SetIdx( idx, idx);
    Poisson.b.SetIdx( idx);
    Poisson.x.SetIdx( idx);
    Poisson.SetupSystem( Poisson.A, Poisson.b);
    MLMatDescCL P;
    P.SetIdx( idx, idx);
    P.Data.resize( MG.GetNumLevel());
    SetupP2ProlongationMatrix( MG, P);

    time.Stop();
    std::cout << "Setting up all stuff took " << time.GetTime()
              << " seconds. " << std::endl;
//    std::cout << "Check Data...\n";
//    CheckMGData( MGData.begin(), MGData.end() );
    MLMatrixCL::const_iterator finest = --Poisson.A.Data.end();
    MLMatrixCL::const_iterator finestP= --P.Data.end();

    Uint sm= 0;
    int lvl= 0;
    Uint nit;
    double tol= 0.0;
    std::cout << "tolerance: "; std::cin >> tol;
    std::cout << "how many levels? (-1=all) > "; std::cin >> lvl;

    double resid, old_resid;
//    JORsmoothCL smoother(omega);  // Jacobi
//    GSsmoothCL smoother(omega);  // Gauss-Seidel
//    SGSsmoothCL smoother(omega);  // symmetric Gauss-Seidel
     SORsmoothCL smoother(omega);  // Gauss-Seidel with over-relaxation
//    SSORsmoothCL smoother(omega);  // symmetric Gauss-Seidel with over-relaxation
    CGSolverCL  solver(200, tol); //CG-Verfahren
    do
    {
        std::cout << "Smoothing steps (0=Quit): "; std::cin >> sm;
        if (sm<=0) return;
        // delete former solution
        Poisson.x.Data.resize(0);
        Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns());
        std::cout << "initial error:" << std::endl;
        Poisson.CheckSolution(&::Lsg);
        resid= norm( Poisson.b.Data - Poisson.A.Data * Poisson.x.Data);
        std::cout << "initial residuum: " << resid <<std::endl;
        nit = 0;
        time.Reset();
        do
        {
            MGM( Poisson.A.Data.begin(), finest, finestP, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, lvl, -1);
            Poisson.CheckSolution(&::Lsg);
            old_resid= resid;
            resid= norm( Poisson.b.Data - Poisson.A.Data * Poisson.x.Data);
            nit = nit+1;
            std::cout << "iteration: " << nit
                      << "\tresiduum: " << resid
                      << "\tred. " << resid/old_resid << std::endl;
        } while ( resid > tol);
        time.Stop();
        std::cout << " time = "<< time.GetTime() <<std::endl;
    } while (sm>0);
}


} //end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

// boundary functions (neumann, dirichlet type)
int main (int argc, char** argv)
{
  try
  {
    if (argc<2)
    {
        std::cout << "missing argument! Usage: MGdrops <omega>" << std::endl;
        return 1;
    }
    double omega;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    typedef DROPS::PoissonP2CL<PoissonCoeffCL> PoissonOnBCL;
    typedef PoissonOnBCL                       MyPoissonCL;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);

    bool isneumann[6]= {false, false, false, false, false, false};
    DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg };
    DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);

    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;

#if 0
    mg.Refine();
    DROPS::MarkAll(mg);
    mg.Refine();

    for (DROPS::Uint i=0; i<1; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
        mg.Refine();
    }
#endif
    std::cout << "Creating Grid..." << std::endl;
    for (DROPS::Uint i=0; i<4; ++i)
    {
        DROPS::MarkAll(mg);
        mg.Refine();
    }
//    omega= 1.0;
    omega= std::atof(argv[1]);
    std::cout << omega << std::endl;
    DROPS::Strategy(prob, omega);
    std::cout << "hallo" << std::endl;
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
//    fil << DROPS::GeomSolOutCL<MyPoissonCL::DiscSolCL>(mg, prob.GetSolution(), &colormap, -1, false, 0.0, prob.x.Data.min(), prob.x.Data.max()) << std::endl;
    mg.SizeInfo(std::cout);
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
