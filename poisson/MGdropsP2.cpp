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
    IdxDescCL* c_idx;
    TimerCL time;
    MGDataCL MGData;

    MG.SizeInfo(std::cerr);
    
    // Initialize the MGData: Idx, A, P, R
    time.Reset();    
    for(Uint i=0, lvl= 0; lvl<=MG.GetLastLevel(); ++lvl, ++i)
    {
        MGData.push_back(MGLevelDataCL());
        MGLevelDataCL& tmp= MGData.back();
     
        std::cerr << "Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 1, 1);
        Poisson.CreateNumbering(lvl, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        if(lvl==MG.GetLastLevel())
        {
            Poisson.x.SetIdx( &tmp.Idx);
            Poisson.b.SetIdx( &tmp.Idx);
            std::cerr << "Create System --->>> " << std::endl;
            Poisson.SetupSystem( tmp.A, Poisson.b);
//	    std::cout << "P2: SetupSystem: A.Data =" << tmp.A.Data << std::endl;
            std::cerr << "Create System <<<--- " << std::endl;
        }
        else        
        {
            std::cerr << "Create StiffMatrix" << std::endl;
            Poisson.SetupStiffnessMatrix( tmp.A);
        }
            
        if(i!=0)
        {
            std::cerr << "Create Prolongation on Level " << lvl << std::endl;
//	    Poisson.SetupProlongation( tmp.P, c_idx, &tmp.Idx);
	    SetupP2ProlongationMatrix( MG, tmp.P, c_idx, &tmp.Idx);
       }
       c_idx= &tmp.Idx;
    }
    time.Stop();
    std::cerr << "Setting up all stuff took " << time.GetTime() 
	      << " seconds. " << std::endl;
//    std::cerr << "Check Data...\n";
//    CheckMGData( MGData.begin(), MGData.end() );
    const_MGDataIterCL finest= --MGData.end();

    Uint sm;
    int lvl;
    Uint nit;
    double tol;
    std::cerr << "tolerance: "; std::cin >> tol;
    std::cerr << "how many levels? (-1=all) > "; std::cin >> lvl;

    double resid, old_resid;
//    JORsmoothCL smoother(omega);  // Jacobi
//    GSsmoothCL smoother(omega);  // Gauss-Seidel
//    SGSsmoothCL smoother(omega);  // symmetric Gauss-Seidel
     SORsmoothCL smoother(omega);  // Gauss-Seidel with over-relaxation
//    SSORsmoothCL smoother(omega);  // symmetric Gauss-Seidel with over-relaxation
    CGSolverCL  solver(200, tol); //CG-Verfahren
    do
    {
        std::cerr << "Smoothing steps (0=Quit): "; std::cin >> sm;
        if (sm<=0) return;
        // delete former solution
        Poisson.x.Data.resize(0);
        Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns);
        std::cerr << "initial error:" << std::endl;
        Poisson.CheckSolution(&::Lsg);
        resid= norm( Poisson.b.Data - finest->A.Data * Poisson.x.Data);
        std::cerr << "initial residuum: " << resid <<std::endl;
	nit = 0;
        time.Reset();    
        do
        {
            MGM( MGData.begin(), finest, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, lvl, -1);
            Poisson.CheckSolution(&::Lsg);
            old_resid= resid;
            resid= norm( Poisson.b.Data - finest->A.Data * Poisson.x.Data);
	    nit = nit+1;
            std::cerr << "iteration: " << nit 
	              << "\tresiduum: " << resid 
	              << "\tred. " << resid/old_resid << std::endl;
        } while ( resid > tol);
        time.Stop();
        std::cerr << " time = "<< time.GetTime() <<std::endl;
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
inline double Null(const DROPS::Point3DCL&, double= 0) { return 0.0; }

int main (int argc, char** argv)
{
  try
  {
    if (argc<2) 
    {
        std::cerr << "missing argument! Usage: MGdrops <omega>" << std::endl;
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
    std::cerr << "Creating Grid..." << std::endl;
    for (DROPS::Uint i=0; i<4; ++i)
    {
        DROPS::MarkAll(mg);  
        mg.Refine();
    }
//    omega= 1.0;
    omega= atof(argv[1]);
    std::cerr << omega << std::endl;
    DROPS::Strategy(prob, omega);
    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
//    fil << DROPS::GeomSolOutCL<MyPoissonCL::DiscSolCL>(mg, prob.GetSolution(), &colormap, -1, false, 0.0, prob.x.Data.min(), prob.x.Data.max()) << std::endl;
    mg.SizeInfo(std::cerr);
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
