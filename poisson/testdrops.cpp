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
    static const double a, b;
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double f(const DROPS::Point3DCL& p)
    {
        const double t0= p.norm();
        const double t1= exp(a*(t0-b));
        if (t0<1.e-5 || t1 > 1.e10)
            return 0.;
        else
        {
            const double t2= 1.0+t1;
            return a*t1*(-2.0*t2 -a*t0+a*t1*t0)/(pow(t2,3)*t0);
        }
    }

};

const double PoissonCoeffCL::a= -60.;
const double PoissonCoeffCL::b= .5;

inline double Lsg( const DROPS::Point3DCL& p)
{
    return 1/(1.0+exp(PoissonCoeffCL::a*(p.norm()-PoissonCoeffCL::b)));
}



/*
class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double f(const DROPS::Point3DCL& p) { return -128.0*( p[0]*p[1]*(1-p[0])*(1-p[1]) + p[0]*p[2]*(1-p[0])*(1-p[2])
                                                        + p[1]*p[2]*(1-p[1])*(1-p[2]) ); }
//    static double f(const Point3DCL& p) { return p[2]>0.49?-15.:0; }
};


inline double Lsg( const DROPS::Point3DCL& p)
{
    return 64.*p[0]*p[1]*p[2]*(1-p[0])*(1-p[1])*(1-p[2]);
}
*/

double root( double base, double deg)
// computes the <deg>-th root of <base>
{
    return exp(log(base)/deg);
}

namespace DROPS // for Strategies
{

template <class DiscSol>
class GeomSolOutReportCL : public MGOutCL
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
    GeomSolOutReportCL (const MultiGridCL& MG, const DiscSol& discsol, const ColorMapperCL* colmap, int TriLevel=-1, bool onlyBnd=false,
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
GeomSolOutReportCL<DiscSol>::put(std::ostream &os) const
{
    const double val_diff= _max-_min;
    ColorMapperCL::RGBAType rgba;
//    std::ios_base::fmtflags my_format= std::ios_base::fixed|std::ios_base::showpoint;
//    std::ios_base::fmtflags old_format= os.flags(my_format);
    std::ios::fmtflags my_format= std::ios::fixed|std::ios::showpoint;
    std::ios::fmtflags old_format= os.flags(my_format);
    Assert(_level==_discsol.GetSolution()->RowIdx->TriangLevel, DROPSErrCL("GeomSolOutCL::put: wrong level"), -1);
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
            if ( fabs(tit->GetVertex(i)->GetCoord()[2] ) < 1.e-10 )
                verts.push_back(i);
        }
        if (verts.size() != 3) continue;
         
        std::vector<double> val;
        val.reserve(4);
        _discsol.GetDoF(*tit, val);

        os << "geom { COFF 3 1 3\n";
//        os << "geom { OFF 3 1 3\n";
        for ( int i=0; i<3; i++ )
        {
            os << tit->GetVertex(verts[i])->GetCoord()[0] << ' '
               << tit->GetVertex(verts[i])->GetCoord()[1] << ' '
               << val[verts[i]]
               << ' '; //'\n';
            rgba= _color->map( (val[i]-_min)/val_diff );
            os << rgba[0] << ' ' << rgba[1] << ' ' << rgba[2] << ' ' << rgba[3] << '\n';
        }
        os <<   "3 0 1 2"
           << "\n}\n";
    }
    os.flags(old_format);
    return os << '}' << std::endl;
}


template<class MGB, class Coeff>
void Strategy(PoissonP1CL<MGB,Coeff>& Poisson, double omega)
{
    typedef PoissonP1CL<MGB,Coeff> MyPoissonCL;
    
    MultiGridCL& MG= Poisson.GetMG();
    IdxDescCL* c_idx;
    TimerCL time;
    MGDataCL MGData;

    for(MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(MG.GetLastLevel()), send=MG.GetTriangVertexEnd(MG.GetLastLevel());
        sit != send; ++sit)
    {
        sit->Unknowns.Destroy();
    }

    // Initialize the MGData: Idx, A, P, R
    time.Reset();
    for(Uint i=0, lvl= 1; lvl<=MG.GetLastLevel(); ++lvl, ++i)
    {
        MGData.push_back(MGLevelDataCL());
        MGLevelDataCL& tmp= MGData.back();
        
        std::cerr << "Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set(i, 1, 0, 0, 0);
        Poisson.CreateNumbering(lvl, &tmp.Idx);
        tmp.A.SetIdx(&tmp.Idx, &tmp.Idx);
        if(lvl==MG.GetLastLevel())
        {
	    time.GReset();
	    time.GStart();
            Poisson.x.SetIdx(&tmp.Idx);
            Poisson.b.SetIdx(&tmp.Idx);
            std::cerr << "Create System " << std::endl;
            Poisson.SetupSystem( tmp.A, Poisson.b);
//            std::cerr << "A(0,0)= " << tmp.A.Data(unk,unk) << std::endl;
	    time.GStop();
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
    time.Stop();
    std::cerr << "Setting up all stuff took " << time.GetTime() 
	      << " seconds including " << time.GetGTime() << " seconds for the largest system." << std::endl;
    std::cerr << "Check Data...\n";
    CheckMGData( MGData.begin(), MGData.end() );
    const_MGDataIterCL finest= --MGData.end();
    std::cerr << finest->Idx.NumUnknowns << " unknowns on finest grid!" << std::endl;
    MG.SizeInfo(std::cerr);
    
    for(Uint level=1; level<=MG.GetLastLevel(); ++level)
    {
        Uint num= std::distance(MG.GetTriangTetraBegin(level), MG.GetTriangTetraEnd(level));
        std::cerr << num << " Tetras on Level " << level << std::endl;
    }
    int meth;
    for(;;) // ever
    {
        std::cerr << "Which method? 0=MG, 1=PCG, -1=Quit > ";  std::cin>>meth;
        if(meth==-1) break;
        double tol;
        std::cerr <<"tol = "; std::cin >> tol;
        time.Reset();
        if (meth)
        {
            int max_iter= 200;
            SsorPcCL<VectorCL, double> pc(omega);
            // delete former solution
            Poisson.x.Data.resize(0);
            Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns);
            std::cerr << "initial error:" << std::endl;
//            Poisson.CheckSolution(&::Lsg);
            double resid= (Poisson.b.Data - finest->A.Data*Poisson.x.Data).norm();
            std::cerr << "initial residuum = " << resid << std::endl;
            time.Start();
            PCG(finest->A.Data, Poisson.x.Data, Poisson.b.Data, pc, max_iter, tol);
            time.Stop();
            std::cerr << "residuum = "<<tol<<" (av.red. "<<root(tol/resid,max_iter)
                      <<"), #iterations = "<<max_iter<< ", time = "<<time.GetTime()<<std::endl;
//            Poisson.CheckSolution(&::Lsg);
        }
        else
        {
            double resid, old_resid;
            WGSSmootherCL<VectorCL, double> smoother(omega);  //gewichtetes Gauss-Seidel
            CGSolverCL   <VectorCL, double> solver(tol, 200); //CG-Verfahren
            Uint sm;
            do
            {
                std::cerr << "Smoothing steps (0=Quit MG): "; std::cin >> sm;
                if (sm<=0) continue;
                // delete former solution
                Poisson.x.Data.resize(0);
                Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns);
                old_resid= resid= (Poisson.b.Data - finest->A.Data*Poisson.x.Data).norm();
                std::cerr << "initial error:" << std::endl;
//                Poisson.CheckSolution(&::Lsg);
                std::cerr << "initial residuum = " << old_resid << std::endl;
                time.Reset();
                Uint step= 0;
                time.Start();
                do
                {
//                    time.Start();
                    MGM( MGData.begin(), finest, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, -1, -1);
//                    time.Stop();
//                    Poisson.CheckSolution(&::Lsg);
//                    old_resid= resid;
                    resid= (Poisson.b.Data - finest->A.Data*Poisson.x.Data).norm();
                    ++step;
//                    std::cerr << "Step "<< step <<": residuum = " << resid << " (red. " << resid/old_resid << "), time = "<<time.GetTime()<<std::endl;
                } while ( resid > tol);
                time.Stop();
                std::cerr << "Step "<< step <<": residuum = " << resid << " (av.red. " << root(resid/old_resid,step) << "), time = "<<time.GetTime()<<std::endl;
//                Poisson.CheckSolution(&::Lsg);
            } while (sm>0);
        }
    }
    std::cerr << "Checking the discretization error: ";
    Poisson.GetDiscError(finest->A, &::Lsg);
    MGData.resize(0);
}


template<class MGB, class Coeff>
void StrategyAdaptive(PoissonP1CL<MGB,Coeff>& Poisson, double omega, double rel_red)
{
    typedef PoissonP1CL<MGB,Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    const typename MyPoissonCL::BndDataCL& BndData= Poisson.GetBndData();
    Uint step= 0, meth;
    IdxDescCL *c_idx;
    TimerCL time;
    MGDataCL MGData;
    
    const_MGDataIterCL finest;
    Uint NumRefined;
    bool new_marks;
    PoissonErrEstCL<typename MyPoissonCL::est_fun, typename MyPoissonCL::_base>
        Estimator(rel_red, 1, 1, true, &Poisson.ResidualErrEstimator, *static_cast<typename MyPoissonCL::_base*>(&Poisson) );

    std::cerr << "Which method? 0=MG, 1=PCG > ";  std::cin>>meth;
    double tol;
    std::cerr <<"tol = "; std::cin >> tol;
    do{
        time.Reset();
        time.Start();
        if (meth) // PCG
        {
            MGData.push_back(MGLevelDataCL());
            MGLevelDataCL& tmp= MGData.back();
            tmp.Idx.Set(1+step%2, 1, 0, 0, 0);
            Poisson.CreateNumbering(MG.GetLastLevel(), &tmp.Idx);
            tmp.A.SetIdx(&tmp.Idx, &tmp.Idx);
            Poisson.b.SetIdx(&tmp.Idx);
            Poisson.x.SetIdx(&tmp.Idx);
            std::cerr << "Create System " << std::endl;
            Poisson.SetupSystem( tmp.A, Poisson.b);
            finest= MGData.begin();
        }
        else
        {
            for(MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(MG.GetLastLevel()), send=MG.GetTriangVertexEnd(MG.GetLastLevel());
                sit != send; ++sit)
            {
                // kill all numberings
                if (sit->Unknowns.Exist()) sit->Unknowns.Destroy();
            }
            // Initialize the MGData: Idx, A, P, R
            for(Uint lvl= 1; lvl<=MG.GetLastLevel(); ++lvl)
            {
                MGData.push_back(MGLevelDataCL());
                MGLevelDataCL& tmp= MGData.back();

                std::cerr << "Create MGData on Level " << lvl << std::endl;
                tmp.Idx.Set(lvl, 1, 0, 0, 0);
                Poisson.CreateNumbering(lvl, &tmp.Idx);
                tmp.A.SetIdx(&tmp.Idx, &tmp.Idx);
                if(lvl!=1)
                {
                    std::cerr << "Create Prolongation on Level " << lvl << std::endl;
                    Poisson.SetupProlongation( tmp.P, c_idx, &tmp.Idx);
                }
                if(lvl==MG.GetLastLevel())
                {
                    Poisson.b.SetIdx(&tmp.Idx);
                    Poisson.x.SetIdx(&tmp.Idx);
                    std::cerr << "Create System " << std::endl;
                    Poisson.SetupSystem( tmp.A, Poisson.b);
                }
                else        
                {
                    std::cerr << "Create StiffMatrix" << std::endl;
                    Poisson.SetupStiffnessMatrix( tmp.A);
                }

                c_idx= &tmp.Idx;
            }
            finest= --MGData.end();
        }
        time.Stop();
        std::cerr << "Setting up all stuff took " << time.GetTime() << " seconds!" << std::endl;
        std::cerr << finest->Idx.NumUnknowns << " unknowns on finest grid!" << std::endl;
        MG.SizeInfo(std::cerr);

        time.Reset();
        if (meth) // PCG
        {
            int max_iter= 200;
            SsorPcCL<VectorCL, double> pc(omega);
            // delete former solution
            Poisson.x.Data.resize(0);
            Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns);
            double resid= (Poisson.b.Data - finest->A.Data * Poisson.x.Data).norm();
            std::cerr << "initial error:" << std::endl;
//            Poisson.CheckSolution(&::Lsg);
            std::cerr << "initial residuum = " << resid << std::endl;
            time.Start();
            PCG(finest->A.Data, Poisson.x.Data, Poisson.b.Data, pc, max_iter, tol);
            time.Stop();
            std::cerr << "residuum = "<<tol<<" (av.red. "<<root(tol/resid,max_iter)
                      <<"), #iterations = "<<max_iter<< ", time = "<<time.GetTime()<<std::endl;
//            Poisson.CheckSolution(&::Lsg);
        }
        else // MGM
        {
            double resid, old_resid;
            WGSSmootherCL<VectorCL, double> smoother(omega);  //gewichtetes Gauss-Seidel
            CGSolverCL   <VectorCL, double> solver(tol, 200); //CG-Verfahren
            Uint sm;
            do
            {
                std::cerr << "Smoothing steps (0=Quit MG): "; std::cin >> sm;
                if (sm<=0) continue;
                // delete former solution
                Poisson.x.Data.resize(0);
                Poisson.x.Data.resize(Poisson.x.RowIdx->NumUnknowns);
                old_resid= resid= (Poisson.b.Data - finest->A.Data * Poisson.x.Data).norm();
                std::cerr << "initial error:" << std::endl;
//                Poisson.CheckSolution(&::Lsg);
                std::cerr << "initial residuum = " << resid << std::endl;
                time.Reset();
                Uint step= 0;
                time.Start();
                do
                {
//                        time.Start();
                    MGM( MGData.begin(), finest, Poisson.x.Data, Poisson.b.Data, smoother, sm, solver, -1, -1);
//                        time.Stop();
//                        Poisson.CheckSolution(&::Lsg);
//                        old_resid= resid;
                    resid= (Poisson.b.Data - finest->A.Data * Poisson.x.Data).norm();
                    ++step;
//                        std::cerr << "Step "<< step <<": residuum = " << resid << " (red. " << resid/old_resid << "), time = "<<time.GetTime()<<std::endl;
                } while ( resid > tol);
                time.Stop();
                std::cerr << "Step "<< step <<": residuum = " << resid << " (av.red. " << root(resid/old_resid,step) << "), time = "<<time.GetTime()<<std::endl;
            } while (sm>0);

        }
        time.Reset();
        time.Start();
        if (step==0)
        {
            Estimator.Init(typename MyPoissonCL::DiscSolCL(&Poisson.x, &BndData, &MG));
        }
//        std::cerr << "maximum is " << x.Data.max() << std::endl;
        // Fehler schaetzen, Verfeinerung
        new_marks= Estimator.Estimate(typename MyPoissonCL::DiscSolCL(&Poisson.x, &BndData, &MG) );
//        NumRefined= EstimateError( x, rel_red, SimpleGradEstimator);
//        std::cerr << NumRefined << " Tetras marked for refinement!" << std::endl;
        time.Stop();
        std::cerr << time.GetTime() <<" seconds for estimating" << std::endl;
        time.Reset();
        time.Start();
        MG.Refine();

//        std::cerr << "Checking the discretization error: ";
//        Poisson.GetDiscError(finest->A, &::Lsg);
        time.Stop();
        std::cerr << time.GetTime() <<" seconds for refining!" << std::endl;
        if(!new_marks)
        {
            // save numbering of solution
            Poisson.idx.Set(0,1,0,0,0);
            Poisson.CreateNumbering( Poisson.x.RowIdx->TriangLevel, &Poisson.idx);
            Poisson.x.RowIdx= &Poisson.idx;
        }
        MGData.resize(0);
        
        ++step;
    } while (new_marks && step<9);
//    } while (NumRefined>0 && step<9);

}


} // end of namespace DROPS

void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

void MarkPart( DROPS::MultiGridCL& mg, DROPS::Uint maxLevel, double percent)
// Marks percent * 100% of the Tetras for refinement ( 0<=percent<=1 )
{
     if (percent<=0) return;

     DROPS::Uint numMark= static_cast<DROPS::Uint>(percent*mg.GetTetras().GetFullSize()),
          num= 1;
     for (DROPS::MultiGridCL::TriangTetraIteratorCL It= mg.GetTriangTetraBegin(maxLevel),
             ItEnd= mg.GetTriangTetraEnd(maxLevel); It!=ItEnd; ++It, ++num)
     {
         if ( num > numMark ) return;
         It->SetRegRefMark();
     }
}

// boundary functions (neumann, dirichlet type)
// used for BndSegCL-object of a UnitCube
inline double neu_val(const DROPS::Point2DCL& p) { return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]); }
inline double dir_val(const DROPS::Point2DCL&) { return 0.0; }


int main (int argc, char** argv)
{
  try
  {
    if (argc<2) 
    {
        std::cerr << "missing argument! Usage: testdrops <omega>" << std::endl;
        return 1;
    }
    double omega= atof(argv[1]);
    DROPS::Point3DCL null(-1.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 2.0;

    typedef DROPS::PoissonP1CL<DROPS::BrickBuilderCL, PoissonCoeffCL> PoissonOnBrickCL;
    typedef PoissonOnBrickCL                                          MyPoissonCL;

    DROPS::BrickBuilderCL domain(null, e1, e2, e3, 4, 4, 4);

    bool isneumann[6]= { false, false, false, false, false, false };
    DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { &dir_val, &dir_val, &dir_val, &dir_val, &dir_val, &dir_val };
    DROPS::PoissonBndDataCL bdata(24, isneumann, bnd_fun);

    MyPoissonCL MGprob(domain, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = MGprob.GetMG();
    DROPS::RBColorMapperCL colormap;

    int adaptiv;
    std::cerr << "Refinement: 0=regular, 1=adaptive > ";   std::cin >> adaptiv;
    if(adaptiv)
    {
        double rel_red;
        std::cerr << "error reduction: ";   std::cin >> rel_red;
        std::cerr << "Creating Grid..." << std::endl;
        for (DROPS::Uint i=0; i<1; ++i)
        {
    //        MarkDrop(mg,mg.GetLastLevel());
            DROPS::MarkAll(mg);  
            mg.Refine();
        }
        StrategyAdaptive(MGprob, omega, rel_red);
    }
    else
    {
        std::cerr << "Creating Grid..." << std::endl;
        for (DROPS::Uint i=0; i<4; ++i)
        {
    //        MarkDrop(mg,mg.GetLastLevel());
            DROPS::MarkAll(mg);  
            mg.Refine();
        }
        Strategy(MGprob, omega);
        double ref;
        std::cerr << "\namount of refinement in % > "; std::cin >> ref;
        if (ref<=0) return 0;
        MarkPart( mg, mg.GetLastLevel(), ref/100);
        mg.Refine();
        Strategy(MGprob, omega);
    }
//    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("grid.off"), fil1("sol2d.off");
    fil << DROPS::GeomSolOutCL<MyPoissonCL::DiscSolCL>(mg, MGprob.GetSolution(), &colormap, -1, false, 0.0, MGprob.x.Data.min(), MGprob.x.Data.max()) << std::endl;
    fil1 << DROPS::GeomSolOutReportCL<MyPoissonCL::DiscSolCL>(mg, MGprob.GetSolution(), &colormap, -1, false, 0.0, MGprob.x.Data.min(), MGprob.x.Data.max()) << std::endl;
//    mg.SizeInfo(std::cerr);
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
