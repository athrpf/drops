#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "poisson/poisson.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include <fstream>

// laplace u + q*u = f
class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double f(const DROPS::Point3DCL& p)
        { return 128.0*( p[0]*p[1]*(1-p[0])*(1-p[1]) + p[0]*p[2]*(1-p[0])*(1-p[2])
	                                             + p[1]*p[2]*(1-p[1])*(1-p[2]) ); }
//    static double f(const Point3DCL& p) { return p[2]>0.49?-15.:0; }
};


inline double Lsg( const DROPS::Point3DCL& p)
{
    return 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
}

namespace DROPS // for strategy
{

template<class Coeff>
void Strategy(PoissonP1CL<Coeff>& Poisson, double omega)
{
    typedef PoissonP1CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    IdxDescCL* c_idx=0;
    MGDataCL MGData;

    
    // Initialize the MGData: Idx, A, P, R
    for(Uint i=0, lvl= 1; lvl<=MG.GetLastLevel(); ++lvl, ++i)
    {
        MGData.push_back(MGLevelDataCL());
        MGLevelDataCL& tmp= MGData.back();
        
        std::cerr << "Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 1);
        Poisson.CreateNumbering(lvl, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        if(lvl==MG.GetLastLevel())
        {
            Poisson.x.SetIdx( &tmp.Idx);
            Poisson.b.SetIdx( &tmp.Idx);
            std::cerr << "Create System " << std::endl;
            Poisson.SetupSystem( tmp.A, Poisson.b);
        }
        else        
        {
            std::cerr << "Create StiffMatrix" << std::endl;
            Poisson.SetupStiffnessMatrix( tmp.A);
        }
            
        if(i!=0)
        {
            std::cerr << "Create Prolongation on Level " << lvl << std::endl;
            Poisson.SetupProlongation( tmp.P, c_idx, &tmp.Idx);
        }
        c_idx= &tmp.Idx;
    }
    std::cerr << "Check Data...\n";
    CheckMGData( MGData.begin(), MGData.end() );
    const_MGDataIterCL finest= --MGData.end();

    Uint sm;
    int lvl;
    double tol;
    std::cerr << "tolerance: "; std::cin >> tol;
    std::cerr << "how many levels? (-1=all) > "; std::cin >> lvl;
    
    double resid, old_resid;
    SORsmoothCL smoother(omega);  //gewichtetes Gauss-Seidel
    CGSolverCL  solver(200, tol); //CG-Verfahren
    do
    {
        std::cerr << "Smoothing steps (0=Quit): "; std::cin >> sm;
        if (sm<=0) return;
        // delete former solution
        Poisson.x.Data.resize(0);
        Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns);
//        std::cerr << "initial error:" << std::endl;
//        Poisson.CheckSolution();
        resid= norm( Poisson.b.Data - finest->A.Data * Poisson.x.Data);
        std::cerr << "initial residuum: " << resid <<std::endl;
        do
        {
            MGM( MGData.begin(), finest, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, lvl, -1);
            Poisson.CheckSolution(&::Lsg);
            old_resid= resid;
                resid= norm( Poisson.b.Data - finest->A.Data * Poisson.x.Data);
            std::cerr << "residuum: " << resid << "\tred. " << resid/old_resid << std::endl;
        } while ( resid > tol);
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
// used for BndSegCL-object of a UnitCube
inline double neu_val(const DROPS::Point2DCL& p) { return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]); }
inline double dir_val(const DROPS::Point2DCL&) { return 0.0; }

// dirichlet value for planes of cube, that has been cut out
inline double dir_val0(const DROPS::Point2DCL& p) { return (1. - p[0]*p[0])*(1. - p[1]*p[1]); }


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

    typedef DROPS::PoissonP1CL<PoissonCoeffCL> PoissonOnBCL;
    typedef PoissonOnBCL                       MyPoissonCL;
//    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);

    bool isneumann[24]= {false, false, false, false, false, false, false, false,
                         false, false, false, false, false, false, false, false,
                         false, false, false, false, false, false, false, false};
    DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[24]=
        { &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val0,
          &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val0,
          &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val0 };
    DROPS::PoissonBndDataCL bdata(24, isneumann, bnd_fun);

    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;

#if 0
    mg.Refine();
    DROPS::MarkAll(mg);
    mg.Refine();

    for (DROPS::Uint i=0; i<4; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        MarkAll(mg);
//        MarkSome(mg);
//        std::cerr << DROPS::SanityMGOutCL(mg);
//        std::cerr << DROPS::DumpMGCL(mg);
//        std::cerr << i << std::endl;
        mg.Refine();
    }
#endif
    std::cerr << "Creating Grid..." << std::endl;
    for (DROPS::Uint i=0; i<3; ++i)
    {
//        MarkDrop(mg,mg.GetLastLevel());
        DROPS::MarkAll(mg);  
        mg.Refine();
    }
    omega= atof(argv[1]);
    std::cerr << omega << std::endl;
    DROPS::Strategy(prob, omega);
    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    fil << DROPS::GeomSolOutCL<MyPoissonCL::DiscSolCL>(mg, prob.GetSolution(), &colormap, -1, false, 0.0, prob.x.Data.min(), prob.x.Data.max()) << std::endl;
//    mg.SizeInfo(std::cerr);
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
