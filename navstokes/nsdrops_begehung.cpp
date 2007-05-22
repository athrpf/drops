#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include "navstokes/navstokes.h"
#include <fstream>


struct NSDrCavCL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL&, double)
    {
        return DROPS::SVectorCL<3>(0.);
    }

    static double LsgPr(const DROPS::Point3DCL&, double)
    {
        return 0;
    }
    static const double st;
    static inline DROPS::SVectorCL<3> Stroem( const DROPS::Point3DCL& p, double)
    {
        const DROPS::SVectorCL<3> ret= 5.*DROPS::std_basis<3>(1);
        const double d0= std::fabs(p[0]-.5);
        const double d1= std::fabs(p[1]-.5);
        const double m= std::max(d0, d1);
        return (.5-st<m) ? ((.5-m)/st)*ret : ret;
    }


    // q*u - nu*laplace u + (u*D)u + Dp = f
    //                           -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL&, double)
            { DROPS::SVectorCL<3> ret(0.0); return ret; }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };

};

const double NSDrCavCL::st= .1;

typedef DROPS::StokesP2P1CL<NSDrCavCL::StokesCoeffCL>
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;
typedef NSDrCavCL MyPdeCL;

namespace DROPS // for Strategy
{

using ::MyStokesCL;

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
//    Assert(_level==_discsol.GetLevel(), DROPSErrCL("GeomSolOutCL::put: wrong level"), ~0);
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
            if ( std::fabs(tit->GetVertex(i)->GetCoord()[1] -0.5 ) < 1.e-10 )
                verts.push_back(i);
        }
        if (verts.size() != 3) continue;

        std::vector<double> val( NumVertsC);
        _discsol.GetDoF(*tit, val);

        os << "geom { OFF 3 1 3\n";
        for ( int i=0; i<3; i++ )
        {
            os << tit->GetVertex(verts[i])->GetCoord()[0] << ' '
               << tit->GetVertex(verts[i])->GetCoord()[2] << ' '
               << std::log(std::fabs(val[verts[i]])+0.5)
               << '\n';
        }
        os <<   "3 0 1 2"
           << "\n}" << std::endl;
    }
    os.flags(old_format);
    return os << '}' << std::endl;
}

template <class DiscVel>
class PlotMTVSolOutCL : public MGOutCL
// output of solution in plotmtv format
{
  private:
    Uint   _level;
    bool   _onlyBnd;
    double _explode;
    const ColorMapperCL* _color;
    DiscVel _discsol;

  public:
    PlotMTVSolOutCL (const MultiGridCL& MG, const DiscVel& discsol, const ColorMapperCL* colmap, int TriLevel=-1, bool onlyBnd=false,
                     double vscale=1.)
        : MGOutCL(&MG), _level( TriLevel<0 ? MG.GetLastLevel() : TriLevel ),
          _onlyBnd(onlyBnd), _explode(vscale), _color(colmap), _discsol(discsol) {}

    void   SetVScale (double explode) { _explode = explode; }
    double GetVScale () const         { return _explode; }
    void   SetColorMap (const ColorMapperCL* colmap) { _color= colmap; }
    virtual std::ostream& put (std::ostream&) const;
};

template <class DiscVel>
std::ostream&
PlotMTVSolOutCL<DiscVel>::put(std::ostream &os) const
{
//    std::ios_base::fmtflags my_format= std::ios_base::fixed|std::ios_base::showpoint;
//    std::ios_base::fmtflags old_format= os.flags(my_format);
    std::ios::fmtflags my_format= std::ios::fixed|std::ios::showpoint;
    std::ios::fmtflags old_format= os.flags(my_format);
//    Assert(_level==_discsol.GetLevel(), DROPSErrCL("GeomSolOutCL::put: wrong level"), ~0);
    os << "# Drops Vector Field\n$ data=vector\n";
    os << "% linewidth = 1 vscale = " << _explode << " toplabel = \"Velocity Field\"\n"
       << "% xlabel = \"X\" ylabel = \"Z\"\n"
       << "% xmin = -0.05 xmax = 2.0 ymin = -0.05 ymax = 1.4\n";
    for ( MultiGridCL::const_TriangVertexIteratorCL tit=_MG->GetTriangVertexBegin(_level);
          tit!=_MG->GetTriangVertexEnd(_level); ++tit )
    {
        if ( std::fabs(tit->GetCoord()[1] -0.5 ) > 1.e-10 )
            continue;

            os << tit->GetCoord()[0] << ' ' << tit->GetCoord()[2] << " 0.0\t"
               << _discsol.val(*tit)[0] << ' ' << _discsol.val(*tit)[2] << _discsol.val(*tit)[1]
               << '\n';
    }
    os.flags(old_format);
    return os  << std::endl;
}

template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, double inner_iter_tol, double tol,
                                   int meth, Uint maxStep, double rel_red, double markratio,
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
        Estimator(rel_red, markratio, 1., true, &MyStokesCL::ResidualErrEstimator, Stokes );
    bool new_marks= false;

    vidx1->Set( 3, 3, 0, 0);
    vidx2->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    pidx2->Set( 1, 0, 0, 0);

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
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &PrBndData, &MG);
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &PrBndData, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &VelBndData, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &VelBndData, &MG);
            Interpolate(pr1, pr2);
//            Interpolate(vel1, vel2);
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
*/
        time.Reset();

        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        Stokes.SetupPrMass( &M);

        double outer_tol= tol;

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
//            Uzawa_PCG_CL uzawaSolver(M.Data, 5000, outer_tol, uzawa_inner_iter, inner_iter_tol, tau);
            Uzawa_IPCG_CL uzawaSolver(M.Data, 5000, outer_tol, uzawa_inner_iter, inner_iter_tol, tau);
            uzawaSolver.Init_A_Pc(A->Data); // only for Uzawa_IPCG_CL.
            time.Start();
            uzawaSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "Iterationen: " << uzawaSolver.GetIter()
                      << "\tNorm des Res.: " << uzawaSolver.GetResid() << std::endl;
        }
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
        if (step==0)
        {
            Estimator.Init(typename MyStokesCL::const_DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::const_DiscVelSolCL(v1, &VelBndData, &MG));
        }
        time.Reset();
        time.Start();
        if (step+1 == maxStep) Estimator.SwitchMark();
        new_marks= Estimator.Estimate(typename MyStokesCL::const_DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::const_DiscVelSolCL(v1, &VelBndData, &MG) );
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


template<class Coeff>
void StrategyNavSt(NavierStokesP2P1CL<Coeff>& NS, int maxStep, double fp_tol, int fp_maxiter,
                                              double uzawa_red, double poi_tol, int poi_maxiter)
// flow control
{
    MultiGridCL& MG= NS.GetMG();

    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &NS.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &NS.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &NS.b;
    VelVecDescCL* c= &NS.c;

    MatDescCL* A= &NS.A;
    MatDescCL* B= &NS.B;
    MatDescCL* N= &NS.N;
    int step= 0;

    vidx1->Set( 3, 3, 0, 0);
    vidx2->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    pidx2->Set( 1, 0, 0, 0);

    TimerCL time;
    do
    {
        MG.Refine();
        NS.CreateNumberingVel(MG.GetLastLevel(), vidx1);
        NS.CreateNumberingPr(MG.GetLastLevel(), pidx1);
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
        if (v2->RowIdx)
        {
            const StokesBndDataCL& BndData= NS.GetBndData();
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
            Interpolate(pr1, pr2);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        N->SetIdx(vidx1, vidx1);
        time.Reset();
        time.Start();
        NS.SetupSystem(A, b, B, c);
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
//            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), Ndat("Nmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
//            Adat << A->Data;   Bdat << B->Data; Ndat << N->Data;   bdat << b->Data;    cdat << c->Data;
//        }

        // adaptive fixedpoint defect correction
        //---------------------------------------
        time.Reset();
        VectorCL d( vidx1->NumUnknowns), e( pidx1->NumUnknowns),
                 w( vidx1->NumUnknowns), q( pidx1->NumUnknowns);
        VelVecDescCL rhsN( vidx1), v_omw( vidx1);
        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        NS.SetupPrMass( &M);
        double omega= 1, res; // initial value (no damping)
        Uzawa_IPCG_CL uzawaSolver(M.Data, 500, -1., poi_maxiter, poi_tol, 1.);
        for(int fp_step=0; fp_step<fp_maxiter; ++fp_step)
        {
            NS.SetupNonlinear( N, v1, &rhsN);
            MatrixCL AN;
            AN.LinComb( 1., A->Data, 1., N->Data);

            // calculate defect:
            d= AN*v1->Data + transp_mul( B->Data, p1->Data) - b->Data - rhsN.Data;
//            e= B->Data*v1->Data                             - c->Data;
            z_xpay(e, B->Data*v1->Data, -1.0, c->Data);

            std::cerr << "fp_step: " << fp_step << ", res = " << (res= std::sqrt( norm_sq( d) + norm_sq( e))) << std::endl;
            if (res < fp_tol )
                break;

            // solve correction:
            double uzawa_tol= res/uzawa_red;
            if (uzawa_tol < fp_tol) uzawa_tol= fp_tol;
            uzawaSolver.SetTol(uzawa_tol);
            uzawaSolver.Init_A_Pc(AN); // only for Uzawa_IPCG_CL.
            uzawaSolver.Solve(AN, B->Data, w, q, d, e);
            std::cerr << "iteration stopped after step " << uzawaSolver.GetIter()
                      << " with res = " << uzawaSolver.GetResid() << std::endl;

            // calculate adaption:
            N->Data.clear();
//            v_omw.Data= v1->Data - omega*w;
            z_xpay(v_omw.Data, v1->Data, -omega, w);
            NS.SetupNonlinear( N, &v_omw, &rhsN);

//            d= A->Data*w + N->Data*w + transp_mul( B->Data, q);
            z_xpaypby2(d, A->Data*w, 1.0, N->Data*w, 1.0, transp_mul(B->Data, q) );
            e= B->Data*w;
            omega= dot( d, VectorCL( A->Data*v1->Data + N->Data*v1->Data
                + transp_mul( B->Data, p1->Data) - b->Data - rhsN.Data))
                + dot( e, VectorCL( B->Data*v1->Data - c->Data));
            omega/= norm_sq( d) + norm_sq( e);
            std::cerr << "omega = " << omega << std::endl;

            // update solution:
//            v1->Data-= omega*w;
            axpy(-omega, w, v1->Data);
//            p1->Data-= omega*q;
            axpy(-omega, q, p1->Data);
        }
        time.Stop();
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";

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
    while (++step<maxStep);
    // we want the solution to be in NS.v, NS.pr
    if (v2 == &loc_v)
    {
        NS.vel_idx.swap( loc_vidx);
        NS.pr_idx.swap( loc_pidx);
        NS.v.SetIdx(&NS.vel_idx);
        NS.p.SetIdx(&NS.pr_idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=9)
    {
        std::cerr << "Usage (stokes): nsdrops <inner_iter_tol> <tol> <meth> <num_refinement> <rel_red> <markratio> <tau> <uz_inner_iter>" << std::endl;
        std::cerr << "Usage (navstokes): <fp_tol> <poi_tol> <fp_maxiter> <poi_maxiter> <uzawa_red> <num_refinement>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    const bool IsNeumann[6]=
        {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &NSDrCavCL::Stroem};

    StokesOnBrickCL stokesprob(brick, NSDrCavCL::StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = stokesprob.GetMG();
    DROPS::RBColorMapperCL colormap;

    {
        double inner_iter_tol= std::atof(argv[1]);
        double tol= std::atof(argv[2]);
        int meth= std::atoi(argv[3]);
        int num_ref= std::atoi(argv[4]);
        double rel_red= std::atof(argv[5]);
        double markratio= std::atof(argv[6]);
        double tau= std::atof(argv[7]);
        unsigned int uz_inner_iter= std::atoi(argv[8]);
        std::cerr << "inner iter tol: " << inner_iter_tol << ", ";
        std::cerr << "tol: " << tol << ", ";
        std::cerr << "meth: " << meth << ", ";
        std::cerr << "refinements: " << num_ref << ", ";
        std::cerr << "relative error reduction: " << rel_red << ", ";
        std::cerr << "markratio: " << markratio << ", ";
        std::cerr << "tau: " << tau << ", ";
        std::cerr << "uzawa inner iter: " << uz_inner_iter << std::endl;
        Strategy(stokesprob, inner_iter_tol, tol, meth, num_ref, rel_red, markratio, tau, uz_inner_iter);

        double min= stokesprob.p.Data.min(),
               max= stokesprob.p.Data.max();
        std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;
        std::ofstream fil("stokespr.off");
        fil << DROPS::GeomSolOutCL<MyStokesCL::const_DiscPrSolCL>(mg, stokesprob.GetPrSolution(), &colormap, -1, false, 0.0, -10, 10) << std::endl;
        std::ofstream fil2("stokespr_cut.off");
        fil2 << DROPS::GeomSolOutReport1CL<MyStokesCL::const_DiscPrSolCL>(mg, stokesprob.GetPrSolution(), &colormap, -1, false, 0.0, 1., 2.) << std::endl;

        DROPS::IdxDescCL tecIdx;
        tecIdx.Set( 1, 0, 0, 0);
        stokesprob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);

        std::ofstream v2d("stokestec2D.dat");
        DROPS::TecPlot2DSolOutCL< MyStokesCL::const_DiscVelSolCL, MyStokesCL::const_DiscPrSolCL>
            tecplot2d( mg, stokesprob.GetVelSolution(), stokesprob.GetPrSolution(), tecIdx, -1, 1, 0.5); // cutplane is y=0.5
        v2d << tecplot2d;
        v2d.close();

        v2d.open("stokes_vel.mtv");
        v2d << DROPS::PlotMTVSolOutCL<MyStokesCL::const_DiscVelSolCL>(mg, stokesprob.GetVelSolution(), &colormap, -1, false, 0.19) << std::endl;
        v2d.close();

        std::ofstream v3d("stokesmaple3D.dat");
        DROPS::MapleSolOutCL<MyStokesCL::const_DiscVelSolCL, MyStokesCL::const_DiscPrSolCL>
            mapleplot( mg, stokesprob.GetVelSolution(), stokesprob.GetPrSolution(), -1, DROPS::PlaneCL(DROPS::std_basis<3>(2), .5));
        v3d << mapleplot;
        v3d.close();

        v3d.open("stokestec3D.dat");
        DROPS::TecPlotSolOutCL< MyStokesCL::const_DiscVelSolCL, MyStokesCL::const_DiscPrSolCL>
            tecplot( mg, stokesprob.GetVelSolution(), stokesprob.GetPrSolution() );
        v3d << tecplot;
        v3d.close();
    }
{
        std::cerr << "Enter nsdrops parameters: <fp_tol> <poi_tol> <fp_maxiter> <poi_maxiter> <uzawa_red> <num_refinement>" << std::endl;
        double fp_tol= 0.0; std::cin >> fp_tol;
        double poi_tol= 0.0; std::cin >> poi_tol;
        int fp_maxiter= 0; std::cin >> fp_maxiter;
        int poi_maxiter= 0; std::cin >> poi_maxiter;
        double uzawa_red= 0.0; std::cin >> uzawa_red;
        int num_ref= 0; std::cin >> num_ref;
        std::cerr << "fp_tol: " << fp_tol<< ", ";
        std::cerr << "poi_tol: " << poi_tol << ", ";
        std::cerr << "fp_maxiter: " << fp_maxiter << ", ";
        std::cerr << "poi_maxiter: " << poi_maxiter << ", ";
        std::cerr << "uzawa_red: " << uzawa_red << ", ";
        std::cerr << "num_ref: " << num_ref << std::endl;

        typedef DROPS::NavierStokesP2P1CL<MyPdeCL::StokesCoeffCL>
                NSOnBrickCL;
        typedef NSOnBrickCL MyNavierStokesCL;

        MyNavierStokesCL prob(stokesprob.GetMG(), MyPdeCL::StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
        DROPS::MultiGridCL& mg = prob.GetMG();

        StrategyNavSt(prob, num_ref, fp_tol, fp_maxiter, uzawa_red, poi_tol, poi_maxiter);

        std::cerr << "hallo" << std::endl;
        std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
        std::ofstream fil("navstokespr.off");
        double min= prob.p.Data.min(),
               max= prob.p.Data.max();
        fil << DROPS::GeomSolOutCL<MyNavierStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;
        fil.close();
        std::ofstream fil2("navstokespr_cut.off");
        fil2 << DROPS::GeomSolOutReport1CL<MyNavierStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, 1., 2.) << std::endl;

        DROPS::IdxDescCL tecIdx;
        tecIdx.Set( 1, 0, 0, 0);
        prob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);

        std::ofstream v2d("navstokestec2D.dat");
        DROPS::TecPlot2DSolOutCL< MyNavierStokesCL::const_DiscVelSolCL, MyNavierStokesCL::const_DiscPrSolCL>
            tecplot2d( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx, -1, 1, 0.5); // cutplane is y=0.5
        v2d << tecplot2d;
        v2d.close();

        v2d.open("navstokes_vel.mtv");
        v2d << DROPS::PlotMTVSolOutCL<MyStokesCL::const_DiscVelSolCL>(mg, stokesprob.GetVelSolution(), &colormap, -1, false, 0.19) << std::endl;
        v2d.close();
    }
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
