#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "poisson/poisson.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include <fstream>




// -laplace u + q*u = f
class PoissonCoeffCL
{
  public:
    static const double a, b;
    static double q(const DROPS::Point3DCL&, double= 0.0) { return 0.0; }
    static double f(const DROPS::Point3DCL& p, double= 0.0)
    {
        const double t0= p.norm();
        const double t1= std::exp(a*(t0-b));
        if (t0<1.e-6 || t1 > 1.e10)
            return 0.;
        else
        {
            const double t2= 1.0+t1;
            return -a*t1*(-2.0*t2 -a*t0+a*t1*t0)/(std::pow(t2,3)*t0);
        }
    }
};

const double PoissonCoeffCL::a= -60.;
const double PoissonCoeffCL::b= .3;


inline double Lsg( const DROPS::Point3DCL& p, double= 0.0)
{
    return 1/(1.0+std::exp(PoissonCoeffCL::a*(p.norm()-PoissonCoeffCL::b)));
}


double root( double base, double deg)
// computes the <deg>-th root of <base>
{
    return std::pow(base, 1./deg);
}

namespace DROPS // for Strategies
{

template <class DiscSol>
class GeomSolOutReport1CL : public MGOutCL
// output of solution in GeomView format
{
  private:
    Uint   _level;
    bool   _onlyBnd;
    double _explode;
    double _min;
    double _max;
    const ColorMapperCL* _color;
    DiscSol _discsol;

  public:
    GeomSolOutReport1CL (const MultiGridCL& MG, const DiscSol& discsol, const ColorMapperCL* colmap, int TriLevel=-1, bool onlyBnd=false,
                         double explode=0.5, double min=0., double max=1.)
        : MGOutCL(&MG), _level( TriLevel<0 ? MG.GetLastLevel() : TriLevel ),
          _onlyBnd(onlyBnd), _explode(explode), _min(min), _max(max), _color(colmap), _discsol(discsol) {}

    void   SetExplode (double explode) { _explode = explode; }
    double GetExplode () const         { return _explode; }
    void   SetMinMax  (double min, double max) { _min= min; _max= max; }
    void   SetColorMap (const ColorMapperCL* colmap) { _color= colmap; }
    virtual std::ostream& put (std::ostream&) const;
};

template <class DiscSol>
std::ostream&
GeomSolOutReport1CL<DiscSol>::put(std::ostream &os) const
{
//    const double val_diff= _max-_min;
//    ColorMapperCL::RGBAType rgba;
//    std::ios_base::fmtflags my_format= std::ios_base::fixed|std::ios_base::showpoint;
//    std::ios_base::fmtflags old_format= os.flags(my_format);
    std::ios::fmtflags my_format= std::ios::fixed|std::ios::showpoint;
    std::ios::fmtflags old_format= os.flags(my_format);
    Assert(_level==_discsol.GetLevel(), DROPSErrCL("GeomSolOutCL::put: wrong level"), ~0);
    os << "appearance {\n-concave\nshading smooth\n}\n";
    os << "LIST {\n";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level);
          tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)) )
            continue;
//        if (GetBaryCenter(*tit)[2]>0.55) continue;
        std::vector<Uint> verts;
        verts.reserve(3);
        for (Uint i=0; i<4; ++i)
        {
            if ( std::fabs(tit->GetVertex(i)->GetCoord()[2] ) < 1.e-10 )
                verts.push_back(i);
        }
        if (verts.size() != 3) continue;

        std::vector<double> val( NumVertsC);
        _discsol.GetDoF(*tit, val);

        os << "geom { OFF 3 1 3\n";
        for ( int i=0; i<3; i++ )
        {
            os << tit->GetVertex(verts[i])->GetCoord()[0] << ' '
               << tit->GetVertex(verts[i])->GetCoord()[1] << ' '
               << val[verts[i]]
               << '\n';
        }
        os <<   "3 0 1 2"
           << "\n}" << std::endl;
    }
    os.flags(old_format);
    return os << '}' << std::endl;
}


template<class Coeff>
void Strategy(PoissonP1CL<Coeff>& Poisson, double omega, double tol, int meth, int sm)
{
    typedef PoissonP1CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    TimerCL time, Ltime;

    for(MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(MG.GetLastLevel()), send=MG.GetTriangVertexEnd(MG.GetLastLevel());
        sit != send; ++sit)
    {
        sit->Unknowns.Destroy();
    }

    // Initialize the MGData: Idx, A, P, R
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
    SetupP1ProlongationMatrix( MG, P);

    time.Stop();
    std::cerr << "Setting up all stuff took " << time.GetTime()
              << " seconds including " << Ltime.GetTime() << " seconds for the largest system." << std::endl;
//    std::cerr << "Check Data...\n";
//    CheckMGData( MGData.begin(), MGData.end() );
    MLMatrixCL::const_iterator finest = --Poisson.A.Data.end();
    MLMatrixCL::const_iterator finestP= --P.Data.end();

    std::cerr << idx->NumUnknowns() << " unknowns on finest grid!" << std::endl;
    MG.SizeInfo(std::cerr);
    std::cerr << "Checking the discretization error: ";
    Poisson.GetDiscError( Poisson.A, &::Lsg);

    for(Uint level=1; level<=MG.GetLastLevel(); ++level)
    {
        Uint num= std::distance(MG.GetTriangTetraBegin(level), MG.GetTriangTetraEnd(level));
        std::cerr << num << " Tetras on Level " << level << std::endl;
    }
//    int meth;
//    for(;;) // ever
    {
//        std::cerr << "Which method? 0=MG, 1=PCG, -1=Quit > ";  std::cin>>meth;
//        if(meth==-1) break;
//        double tol;
//        std::cerr <<"tol = "; std::cin >> tol;
        time.Reset();
        if (meth)
        {
            int max_iter= 200;
            SSORPcCL pc(omega);

            // delete former solution
            Poisson.x.Data.resize(0);
            Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns());
//            std::cerr << "initial error:" << std::endl;
//            Poisson.CheckSolution(&::Lsg);
            double resid= norm( Poisson.b.Data - Poisson.A.Data*Poisson.x.Data);
            std::cerr << "initial residuum = " << resid << std::endl;
            time.Start();
            PCG(Poisson.A.Data, Poisson.x.Data, Poisson.b.Data, pc, max_iter, tol);
            time.Stop();
//            std::cerr << "residuum = "<<tol<<" (av.red. "<<root(tol/resid,max_iter)
//                      <<"), #iterations = "<<max_iter<< ", time = "<<time.GetTime()<<std::endl;
            std::cerr << "residuum = "<<tol<<", av.red. "<<root(tol/resid,max_iter)
                      <<", #iterations = "<<max_iter<< ", time = "<<time.GetTime()<<std::endl;
            Poisson.CheckSolution(&::Lsg);
        }
        else
        {
            double resid, old_resid;
            SORsmoothCL smoother(omega);  //gewichtetes Gauss-Seidel
            CGSolverCL  solver(200, tol); //CG-Verfahren
//            Uint sm;
//            do
//            {
//                std::cerr << "Smoothing steps (0=Quit MG): "; std::cin >> sm;
//                if (sm<=0) continue;
                // delete former solution
                Poisson.x.Data.resize(0);
                Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns());
                old_resid= resid= norm( Poisson.b.Data - Poisson.A.Data*Poisson.x.Data);
//                std::cerr << "initial error:" << std::endl;
//                Poisson.CheckSolution(&::Lsg);
                std::cerr << "initial residuum = " << old_resid << std::endl;
                time.Reset();
                Uint step= 0;
                time.Start();
                do
                {
//                    time.Start();
                    MGM( Poisson.A.Data.begin(), finest, finestP, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, -1, -1);
//                    time.Stop();
//                    Poisson.CheckSolution(&::Lsg);
//                    old_resid= resid;
                    resid= norm( Poisson.b.Data - Poisson.A.Data*Poisson.x.Data);
                    ++step;
//                    std::cerr << "Step "<< step <<": residuum = " << resid << ", red. " << resid/old_resid << ", time = "<<time.GetTime()<<std::endl;
                } while ( resid > tol);
                time.Stop();
                std::cerr << "Step "<< step <<": residuum = " << resid << ", av.red. " << root(resid/old_resid,step) << ", time = "<<time.GetTime()<<std::endl;
                Poisson.CheckSolution(&::Lsg);
//            } while (sm>0);
        }
    }
}


template<class Coeff>
void StrategyAdaptive(PoissonP1CL<Coeff>& Poisson, double omega,
                                                   double tol, int meth, int sm,
                                                   double stoperr, double markratio,
                                                   double minratio)
{
    typedef PoissonP1CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    const typename MyPoissonCL::BndDataCL& BndData= Poisson.GetBndData();
    Uint step= 0;
    double true_err= 0.0;
    TimerCL time, time2;

    MLIdxDescCL* idx = &Poisson.idx;
    MLMatDescCL P;
    MLMatrixCL::const_iterator finest;
    MLMatrixCL::const_iterator finestP;

    bool new_marks;
    DoerflerMarkCL<typename MyPoissonCL::est_fun, typename MyPoissonCL::_base>
        Estimator(1e-10, minratio, markratio, 8, true, &MyPoissonCL::ResidualErrEstimator, *static_cast<typename MyPoissonCL::_base*>(&Poisson) );
    DoerflerMarkCL<typename MyPoissonCL::est_fun, typename MyPoissonCL::_base>
        EstimatorL2(1e-10, 0., markratio, 8, false, &MyPoissonCL::ResidualErrEstimatorL2, *static_cast<typename MyPoissonCL::_base*>(&Poisson) );

//    std::cerr << "Which method? 0=MG, 1=PCG > ";  std::cin>>meth;
//    double tol;
//    std::cerr <<"tol = "; std::cin >> tol;
    do{
        std::cout << "Global Step " << step << ". ";
        time.Reset();
        time2.Start();
        time.Start();
        MG.Refine();
        time.Stop();
        time2.Stop();
        std::cerr << time.GetTime() <<" seconds for refining!" << std::endl;

        time.Reset();
        time2.Start();
        time.Start();
        if (!meth) // !PCG
        {
            for(MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(MG.GetLastLevel()), send=MG.GetTriangVertexEnd(MG.GetLastLevel());
                sit != send; ++sit)
            {
                // kill all numberings
                if (sit->Unknowns.Exist()) sit->Unknowns.Destroy();
            }
            std::cerr << "Create System " << std::endl;
            Poisson.SetNumLvl( MG.GetNumLevel());
        }
        Poisson.DeleteNumbering( &Poisson.idx);
        Poisson.CreateNumbering( MG.GetLastLevel(), &Poisson.idx);
        Poisson.A.SetIdx( idx, idx);
        Poisson.b.SetIdx( idx);
        Poisson.x.SetIdx( idx);
        Poisson.SetupSystem( Poisson.A, Poisson.b);
        MLMatrixCL P;
        if (!meth) SetupP1ProlongationMatrix( MG, P, idx, idx);
        finest = --Poisson.A.Data.end();
        finestP= --P.end();

        time.Stop();
        time2.Stop();
        std::cerr << "Setting up all stuff took " << time.GetTime() << " seconds!" << std::endl;
        std::cerr << idx->NumUnknowns() << " unknowns on finest grid!" << std::endl;
        MG.SizeInfo(std::cerr);

        time.Reset();
        if (meth) // PCG
        {
            int max_iter= 200;
            double mytol= tol;
            SSORPcCL pc(omega);
            // delete former solution
            Poisson.x.Data.resize(0);
            Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns());
            double resid= norm( Poisson.b.Data - Poisson.A.Data * Poisson.x.Data);
//            std::cerr << "initial error:" << std::endl;
//            Poisson.CheckSolution(&::Lsg);
            std::cerr << "initial residuum = " << resid << std::endl;
            time2.Start();
            time.Start();
            PCG(Poisson.A.Data, Poisson.x.Data, Poisson.b.Data, pc, max_iter, mytol);
            time.Stop();
            time2.Stop();
            std::cerr << "residuum = "<<mytol<<", av.red. "<<root(tol/resid,max_iter)
                      <<", #iterations = "<<max_iter<< ", time = "<<time.GetTime()<<std::endl;
        }
        else // MGM
        {
            double resid, old_resid;
            SORsmoothCL smoother(omega);  //gewichtetes Gauss-Seidel
            CGSolverCL  solver(200, tol); //CG-Verfahren
//            Uint sm;
//            do
//            {
//                std::cerr << "Smoothing steps (0=Quit MG): "; std::cin >> sm;
//                if (sm<=0) continue;
                // delete former solution
                Poisson.x.Data.resize(0);
                Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns());
                old_resid= resid= norm( Poisson.b.Data - Poisson.A.Data * Poisson.x.Data);
//                std::cerr << "initial error:" << std::endl;
//                Poisson.CheckSolution(&::Lsg);
                std::cerr << "initial residuum = " << resid << std::endl;
                time.Reset();
                Uint step2= 0;
                time2.Start();
                time.Start();
                do
                {
//                        time.Start();
                    MGM( Poisson.A.Data.begin(), finest, finestP, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, -1, -1);
//                        time.Stop();
//                        Poisson.CheckSolution(&::Lsg);
//                        old_resid= resid;
                    resid= norm( Poisson.b.Data - Poisson.A.Data * Poisson.x.Data);
                    ++step2;
//                        std::cerr << "Step "<< step2 <<": residuum = " << resid << " (red. " << resid/old_resid << "), time = "<<time.GetTime()<<std::endl;
                } while ( resid > tol);
                time.Stop();
                time2.Stop();
                std::cerr << "Step "<< step2 <<": residuum = " << resid << ", av.red. " << root(resid/old_resid,step2) << ", time = "<<time.GetTime()<<std::endl;
//            } while (sm>0);

        }
        if (step==0)
        {
            Estimator.Init(typename MyPoissonCL::DiscSolCL(&Poisson.x, &BndData, &MG));
            EstimatorL2.Init(typename MyPoissonCL::DiscSolCL(&Poisson.x, &BndData, &MG));
        }
//        std::cerr << "maximum is " << x.Data.max() << std::endl;
        // Fehler schaetzen, Verfeinerung
        true_err= Poisson.CheckSolution(&::Lsg);
        time.Reset();
        time2.Start();
        time.Start();
        new_marks= Estimator.Estimate(typename MyPoissonCL::const_DiscSolCL(&Poisson.x, &BndData, &MG) );
        time.Stop();
        time2.Stop();
        std::cerr << time.GetTime() <<" seconds for estimating" << std::endl;
        EstimatorL2.Estimate(typename MyPoissonCL::const_DiscSolCL(&Poisson.x, &BndData, &MG) );

//        std::cerr << "Checking the discretization error: ";
//        Poisson.GetDiscError(finest->A, &::Lsg);

        ++step;
    } while ( true_err > stoperr /*&& step<14*/ );
    std::cout << "cumulative solving time: " << time2.GetTime() << std::endl;
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

void MarkPart( DROPS::MultiGridCL& mg, DROPS::Uint maxLevel, double ratio)
// Marks ratio * 100% of the Tetras for refinement ( 0<=ratio<=1 )
{
     if (ratio<=0) return;

     DROPS::Uint numMark= static_cast<DROPS::Uint>(ratio*mg.GetTetras().size()),
          num= 1;
     for (DROPS::MultiGridCL::TriangTetraIteratorCL It= mg.GetTriangTetraBegin(maxLevel),
             ItEnd= mg.GetTriangTetraEnd(maxLevel); It!=ItEnd; ++It, ++num)
     {
         if ( num > numMark ) return;
// TODO: Shouldn't we mark only leaves?
         It->SetRegRefMark();
     }
}

// boundary functions (neumann, dirichlet type)
inline double dir_val(const DROPS::Point3DCL&, double= 0.0) { return 1.; }


int main (int argc, char** argv)
{
  try
  {
    if (argc!=11)
    {
        std::cerr << "missing argument! Usage: testdrops <omega> <file> <tol> <meth> <sm> <adap> <grids> <stoperr> <markratio> <minratio>" << std::endl;
        return 1;
    }
    double omega= std::atof(argv[1]);
    const char* filename= argv[2];
    double tol= std::atof(argv[3]);
    int meth= std::atoi(argv[4]);
    int sm= std::atoi(argv[5]);
    int adaptiv= std::atoi(argv[6]);
    int grids= std::atoi(argv[7]);
    double stoperr= std::atof(argv[8]);
    double markratio= std::atof(argv[9]);
    double minratio= std::atof(argv[10]);
    std::cout << "omega " << omega << ", ";
    std::cout << "filename " << filename << ", ";
    std::cout << "tol " << tol << ", ";
    std::cout << "meth " << meth << ", ";
    std::cout << "sm " << sm << ", ";
    std::cout << "adap " << adaptiv << ", ";
    std::cout << "grids " << grids << ", ";
    std::cout << "stoperr " << stoperr << ", ";
    std::cout << "markratio " << markratio << ", ";
    std::cout << "minratio " << minratio << std::endl;
    DROPS::Point3DCL orig(-1.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 2.0;

    typedef DROPS::PoissonP1CL<PoissonCoeffCL> PoissonOnBrickCL;
    typedef PoissonOnBrickCL                   MyPoissonCL;

    DROPS::BrickBuilderCL domain(orig, e1, e2, e3, 4, 4, 4);
//    DROPS::LBuilderCL domain(null, e1, e2, e3, 2, 2, 2, 2, 2);
//    DROPS::BBuilderCL domain(null, e1, e2, e3, 2, 2, 2, 1, 1, 1);

    bool isneumann[6]= { false, false, false, false, false, false };
    DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val };
    DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);

    MyPoissonCL MGprob(domain, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = MGprob.GetMG();
    DROPS::RBColorMapperCL colormap;

//    std::cerr << "Refinement: 0=regular, 1=adaptive > ";   std::cin >> adaptiv;
    if(adaptiv)
    {
//        double rel_red;
//        std::cerr << "error reduction: ";   std::cin >> rel_red;
        std::cerr << "Creating Grid..." << std::endl;
        for (int i=0; i<1; ++i)
        {
    //        MarkDrop(mg,mg.GetLastLevel());
            DROPS::MarkAll(mg);
            mg.Refine();
        }
        StrategyAdaptive(MGprob, omega, tol, meth, sm, stoperr, markratio, minratio);
    }
    else
    {
//        std::cerr << "Creating Grid..." << std::endl;
        for (int i=0; i<grids; ++i)
        {
    //        MarkDrop(mg,mg.GetLastLevel());
            DROPS::MarkAll(mg);
            mg.Refine();
        }
        Strategy(MGprob, omega, tol, meth, sm);
    }
//    std::cerr << "hallo" << std::endl;
//    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil(filename);
    fil << DROPS::GeomSolOutReport1CL<MyPoissonCL::DiscSolCL>(mg, MGprob.GetSolution(), &colormap, -1, false, 0.0, MGprob.x.Data.min(), MGprob.x.Data.max()) << std::endl;
//    std::ofstream fil2("ttt2.off");
//    fil2 << DROPS::GeomSolOutCL<MyPoissonCL::DiscSolCL>(mg, MGprob.GetSolution(), &colormap, -1, false, 0.0, MGprob.x.Data.min(), MGprob.x.Data.max()) << std::endl;
//    mg.SizeInfo(std::cerr);
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
