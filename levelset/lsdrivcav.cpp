//**************************************************************************
// File:    lsdrivcav.cpp                                                  *
// Content: drop in driven cavity flow + surface tension                   *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************


#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/levelset.h"
#include <fstream>

// q*u - nu*laplace u + Dp = f - okn
//                  -div u = 0
class ZeroFlowCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL&, double)
        { DROPS::SVectorCL<3> ret(0.0); return ret; }
    const double nu;
    
    ZeroFlowCL() : nu(1.0) {}
};

// Tropfendaten:
DROPS::Point3DCL Mitte(0.5);
double           Radius= 0.3;

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)   { return DROPS::SVectorCL<3>(0.); }
DROPS::SVectorCL<3> Stroem( const DROPS::Point3DCL&, double) { DROPS::SVectorCL<3> ret(0.); ret[0]= 1.; return ret; }

double DistanceFct( const DROPS::Point3DCL& p)
{
    return (Mitte-p).norm()-Radius;
}


namespace DROPS // for Strategy
{

template<class Coeff>
void Strategy( StokesP2P1CL<Coeff>& Stokes, double inner_iter_tol, double sigma)
// flow control
{
    typedef StokesP2P1CL<Coeff> StokesProblemT;
    
    MultiGridCL& MG= Stokes.GetMG();
    LevelsetP2CL lset( MG, sigma);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    VelVecDescCL* v= &Stokes.v;
    VecDescCL*    p= &Stokes.p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;
    MatDescCL* A= &Stokes.A;
    MatDescCL* B= &Stokes.B;
    
    TimerCL time;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);    
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);    
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    MG.SizeInfo(std::cerr);
    b->SetIdx(vidx);
    c->SetIdx(pidx);
    p->SetIdx(pidx);
    v->SetIdx(vidx);
    std::cerr << "Anzahl der Druck-Unbekannten: " << p->Data.size() << std::endl;
    std::cerr << "Anzahl der Geschwindigkeitsunbekannten: " << v->Data.size() << std::endl;
    A->Reset();
    B->Reset();
    A->SetIdx(vidx, vidx);
    B->SetIdx(pidx, vidx);
    time.Reset();
    time.Start();
    Stokes.SetupSystem(A, b, B, c);
    time.Stop();
    std::cerr << time.GetTime() << " seconds for setting up all systems!" << std::endl;
    time.Reset();
    time.Start();
    lset.AccumulateBndIntegral( *b);
    time.Stop();
    std::cerr << time.GetTime() << " seconds for coupling term!" << std::endl;

    MatDescCL M;
    M.SetIdx( pidx, pidx);
    Stokes.SetupPrMass( &M);

    Uint meth;
    std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
    time.Reset();
    double outer_tol;
    std::cerr << "tol = "; std::cin >> outer_tol;

    if (meth)
    {
//            PSchur_PCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
        PSchur_GSPCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
        time.Start();
        schurSolver.Solve( A->Data, B->Data, v->Data, p->Data, b->Data, c->Data);
        time.Stop();
    }
    else // Uzawa
    {
        double tau;
        Uint inner_iter;
        std::cerr << "tau = "; std::cin >> tau;
        std::cerr << "#PCG steps = "; std::cin >> inner_iter;
        Uzawa_PCG_CL uzawaSolver( M.Data, 5000, outer_tol, inner_iter, inner_iter_tol, tau);
        time.Start();
        uzawaSolver.Solve( A->Data, B->Data, v->Data, p->Data, b->Data, c->Data);
        time.Stop();
        std::cerr << "Iterationen: " << uzawaSolver.GetIter()
                  << "\tNorm des Res.: " << uzawaSolver.GetResid() << std::endl;
    }
    std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
    std::cerr << std::endl;
    
    EnsightP2SolOutCL ensight( MG, lidx);
    
    const char datgeo[]= "ensight/curv.geo", 
               datpr[] = "ensight/curv.pr",
               datvec[]= "ensight/curv.vel",
               datscl[]= "ensight/curv.scl";
    ensight.CaseBegin( "curv.case");
    ensight.DescribeGeom( "curvature driven flow", datgeo);
    ensight.DescribeScalar( "Levelset", datscl); 
    ensight.DescribeScalar( "Pressure", datpr ); 
    ensight.DescribeVector( "Velocity", datvec); 
    ensight.putGeom( datgeo);
    ensight.putVector( datvec, Stokes.GetVelSolution());
    ensight.putScalar( datpr,  Stokes.GetPrSolution());
    ensight.putScalar( datscl, lset.GetSolution());
    ensight.CaseEnd();
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=4)
    {
        std::cerr << "You have to specify three parameters:\n\tlscoupl <inner_iter_tol> <num_subdiv> <surf.tension>" << std::endl;
        return 1;
    }

    double inner_iter_tol= atof(argv[1]);
    int sub_div= atoi(argv[2]);
    double sigma= atof(argv[3]);
    std::cerr << "inner iter tol:  " << inner_iter_tol << std::endl;
    std::cerr << "sub divisions:   " << sub_div << std::endl;
    std::cerr << "surface tension: " << sigma << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    typedef DROPS::StokesP2P1CL<ZeroFlowCL> MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, sub_div, sub_div, sub_div);
    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Null }; 
        
    MyStokesCL prob(brick, ZeroFlowCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();

    Strategy(prob, inner_iter_tol, sigma);

    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
