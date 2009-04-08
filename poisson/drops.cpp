#include "geom/multigrid.h"
#include "geom/builder.h"
#include "poisson/poisson.h"
#include "out/output.h"
#include "num/solver.h"
#include <fstream>

// -laplace u + q*u = f
class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&, double= 0.0) { return 0.0; }
    static double f(const DROPS::Point3DCL& p, double= 0.0)
        { return 128.0*( p[0]*p[1]*(1-p[0])*(1-p[1]) + p[0]*p[2]*(1-p[0])*(1-p[2])
                                                     + p[1]*p[2]*(1-p[1])*(1-p[2]) ); }
//    static double f(const Point3DCL& p, double= 0.0) { return p[2]>0.49?-15.:0; }
};


inline double Lsg( const DROPS::Point3DCL& p, double= 0.0)
{
    return 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
}


namespace DROPS // for Strategy
{

template<class Coeff>
void Strategy(PoissonP1CL<Coeff>& Poisson, double omega, double rel_red)
{
    typedef PoissonP1CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    const typename MyPoissonCL::BndDataCL& BndData= Poisson.GetBndData();

    MLIdxDescCL  loc_idx;
    VecDescCL  loc_x;

    MLIdxDescCL* new_idx= &Poisson.idx;
    MLIdxDescCL* old_idx= &loc_idx;
//    IdxDescCL* err_idx= &_err_idx;
    VecDescCL* new_x= &Poisson.x;
    VecDescCL* old_x= &loc_x;
    VecDescCL* b= &Poisson.b;
//    VecDescCL* err= &_err;
    MLMatDescCL* A= &Poisson.A;
    SSORPcCL   pc(omega);
    DoerflerMarkCL<typename MyPoissonCL::est_fun, typename MyPoissonCL::_base>
        Estimator(rel_red, 0.0, 0.6, 0.875, true, &MyPoissonCL::ResidualErrEstimator, *static_cast<typename MyPoissonCL::_base*>(&Poisson) );
    Uint step= 0;
    bool new_marks;

    new_idx->SetFE( P1_FE);
    old_idx->SetFE( P1_FE);
//    err_idx->Set( 2, 0, 0, 0, 1);
    do
    {
        std::cout << DROPS::SanityMGOutCL(MG) << std::endl;
        MG.Refine();
        Poisson.CreateNumbering( MG.GetLastLevel(), new_idx);    // erzeuge Nummerierung zu diesem Index
        std::cout << "altes und neues TriangLevel: " << old_idx->TriangLevel() << ", "
                  << new_idx->TriangLevel() << std::endl;
        b->SetIdx( new_idx);                        // Erster Vektor aus new_idx
        new_x->SetIdx( new_idx);                    // Zweiter Vektor ebenfalls aus new_idx
        std::cout << "Anzahl der Unbekannten: " << old_x->Data.size() << ", "
                  << new_x->Data.size() << std::endl;
        MG.SizeInfo(std::cout);
        if (step==0)
        {
            Estimator.Init(typename MyPoissonCL::DiscSolCL(new_x, &BndData, &MG) );
        }
        if (old_x->RowIdx)
        {
            P1EvalCL<double, const PoissonBndDataCL, const VecDescCL>  oldx(old_x, &BndData, &MG);
            P1EvalCL<double, const PoissonBndDataCL, VecDescCL>        newx(new_x, &BndData, &MG);
// TODO: what happens, if refinement creates new vertices in old_level? FIXED in fe.h!
            Interpolate(newx, oldx);
//            CheckSolution(*new_x, &::Lsg);
            old_x->Reset();
        }
        A->SetIdx(new_idx, new_idx);               // Erste Matrix aus Index 0 (Zeilen und Spalten)
        Poisson.SetupSystem(*A, *b);
//        Poisson.GetDiscError(&::Lsg);
//std::cout << A->Data << std::endl << b->Data << std::endl;
        double tol= 1.0e-7;
        int max_iter= 100;
        PCG(A->Data, new_x->Data, b->Data, pc, max_iter, tol);
//        CG(A->Data, new_x->Data, b->Data, max_iter, tol);
        A->Reset();
        b->Reset();
        std::cout << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;
//        std::cout << "Loesung: " << new_x->Data << std::endl;
        Poisson.CheckSolution(*new_x, &::Lsg);
//        CreateNumbering(new_idx->TriangLevel, err_idx);
//        err->SetIdx( err_idx);
//        NumMarked= EstimateError(*new_x, 0.2, &Estimator);
        new_marks= Estimator.Estimate( typename MyPoissonCL::const_DiscSolCL(new_x, &BndData, &MG) );
//        err->Reset();
//        DeleteNumbering(err_idx);
        std::swap(old_x, new_x);
        std::swap(old_idx, new_idx);
        std::cout << std::endl;
    }
    while (new_marks && step++<6);
    // I want the solution to be in Poisson.x
    if (old_x == &loc_x)
    {
        Poisson.idx.swap( loc_idx);
        Poisson.x.SetIdx(&Poisson.idx);

        Poisson.x.Data.resize(loc_x.Data.size());
        Poisson.x.Data= loc_x.Data;
    }
}

} // end of namespace DROPS

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
    if (argc!=3)
    {
        std::cout << "You have to specify two parameters:\n\tdrops <omega> <err_red>" << std::endl;
        return 1;
    }
    double omega;
    double rel_red;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    typedef DROPS::PoissonP1CL<PoissonCoeffCL> PoissonOnBCL;
    typedef PoissonOnBCL                       MyPoissonCL;

//    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    const bool isneumann[24]=
        {false, false, false, false, false, false, false, false,
         false, false, false, false, false, false, false, false,
         false, false, false, false, false, false, false, false};
    const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[24]=
        { &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg,
          &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg,
          &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg, &Lsg };

    DROPS::PoissonBndDataCL bdata(24, isneumann, bnd_fun);
    PoissonOnBCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;

    std::ofstream maple("ttt.txt");
    maple << DROPS::MapleMGOutCL(mg, -1, false, true, DROPS::PlaneCL(e3, 0.6)) << std::endl;

//    MarkAll(mg);


//    mg.Refine();
//   MarkAll(mg);
//    mg.Refine();
/*
    for (DROPS::Uint i=0; i<4; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        MarkAll(mg);
//        MarkSome(mg);
        std::cout << DROPS::SanityMGOutCL(mg);
//        std::cout << DROPS::DumpMGCL(mg);
//        std::cout << i << std::endl;
        mg.Refine();
    }
*/
//    MarkAll(mg);
    omega= std::atof(argv[1]);
    rel_red= std::atof(argv[2]);
    std::cout << "Omega: " << omega << " rel_red: " << rel_red << std::endl;
    DROPS::Strategy(prob, omega, rel_red);
    std::cout << "hallo" << std::endl;
    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    fil << DROPS::GeomSolOutCL<MyPoissonCL::DiscSolCL>(mg, prob.GetSolution(), &colormap, -1, false, 0.0, prob.x.Data.min(), prob.x.Data.max()) << std::endl;
//    std::cout.flags(std::ios::fixed|std::ios::showpoint);
//    for ( DROPS::MultiGridCL::const_TriangVertexIteratorCL tit=const_cast<const DROPS::MultiGridCL&>(mg).GetTriangVertexBegin(prob.GetSolution().GetTriangLevel()),
//          theend=const_cast<const DROPS::MultiGridCL&>(mg).GetTriangVertexEnd(prob.GetSolution().GetTriangLevel()); tit!=theend; ++tit )
//        std::cout << prob.GetSolution().val(*tit) << '\n';
//    std::cout << DROPS::GeomSolOutCL<DROPS::PoissonLinearDirectCL::PoissonDiscSolCL>(mg, prob.GetSolution(), &colormap, -1, false, 0.0, 0.0, 1.1) << std::endl;
//    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
//    std::cout << "Verts: "   << mg.GetVertices().size()
//              << " Edges: "  << mg.GetEdges().size()
//              << " Tetras: " << mg.GetTetras().size()
//              << std::endl;
/*
//int wait;
//cin>>wait;
    for (DROPS::Uint i=0; i<6; ++i)
    {
        UnMarkDrop(mg, min(mg.GetLastLevel(),10u));
//        std::cout << DROPS::SanityMGOutCL(mg);
//        std::cout << i << std::endl;
        mg.Refine();
    }
    UnMarkAll(mg);
    mg.Refine();
    UnMarkAll(mg);
    mg.Refine();
    UnMarkAll(mg);
    mg.Refine();
    UnMarkAll(mg);
    mg.Refine();
*/
//    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
/*
    std::cout << "Verts: "   << mg.GetVertices().size()
         << " Edges: "  << mg.GetEdges().size()
         << " Tetras: " << mg.GetTetras().size()
         << std::endl;
*/
//    std::cout << DROPS::GeomMGOutCL(mg, -1, true) << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
