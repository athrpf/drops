//**************************************************************************
// File:    lsshear.cpp                                                    *
// Content: drop in shear flow                                             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include <fstream>

const double      delta_t= 0.01;
const DROPS::Uint num_steps= 50;
const int         FPsteps= -1;

// du/dt - q*u - nu*laplace u + Dp = f - okn
//                          -div u = 0
//                               u = u0, t=t0


class ShearFlowCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL&, double)
        { DROPS::SVectorCL<3> ret(0.0); return ret; }
    const double nu;
    
    ShearFlowCL() : nu(1.0) {}
};

// Randdaten: x=0, x=1, y=0, y=1:  Dirichlet 0
//            z=0 und x<0.5        Neumann   0   (aus Impl.gruenden: Dir.)
//            z=0 und x>0.5        Inflow Dirichlet  parabol.
//            z=1 und x<0.5        Inflow Dirichlet  parabol.
//            z=1 und x>0.5        Neumann   0   (aus Impl.gruenden: Dir.)


// Tropfendaten:
DROPS::Point3DCL Mitte(0.5);
double           Radius= 0.2;

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)   { return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Parabol( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.); 
    if (p[0]<0.5)
        ret[2]= 4*p[0]*(p[0]-0.5);
    else
        ret[2]= 4*(1-p[0])*(p[0]-0.5);
    return ret;
}

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
    LevelsetP2CL lset( MG, sigma, 0.5, 0.1);
    
    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    VelVecDescCL* v= &Stokes.v;
    VecDescCL*    p= &Stokes.p;
    VelVecDescCL* b= &Stokes.b;
    VecDescCL* c= &Stokes.c;
    VelVecDescCL cpl_M;
    MatDescCL* A= &Stokes.A;
    MatDescCL* B= &Stokes.B;
    MatDescCL* M= &Stokes.M;
    MatDescCL prM;

    TimerCL time;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);    
    Stokes.CreateNumberingPr( MG.GetLastLevel(), pidx);    
    lset.CreateNumbering(     MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    MG.SizeInfo( std::cerr);
    b->SetIdx( vidx);
    c->SetIdx( pidx);
    cpl_M.SetIdx( vidx);
    p->SetIdx( pidx);
    v->SetIdx( vidx);
    std::cerr << "Anzahl der Druck-Unbekannten: " << p->Data.size() << std::endl;
    std::cerr << "Anzahl der Geschwindigkeitsunbekannten: " << v->Data.size() << std::endl;
    A->Reset();
    B->Reset();
    M->Reset();
    A->SetIdx(vidx, vidx);
    B->SetIdx(pidx, vidx);
    M->SetIdx(vidx, vidx);
    prM.SetIdx( pidx, pidx);
    time.Reset();
    time.Start();
    Stokes.SetupInstatSystem( A, B, M);
    Stokes.SetupInstatRhs( b, c, &cpl_M, 0, b, 0);
    Stokes.SetupPrMass( &prM);
    time.Stop();
    std::cerr << time.GetTime() << " seconds for setting up all systems!" << std::endl;
    
    Stokes.InitVel( v, Null);
    lset.SetupSystem( Stokes.GetVelSolution() );

    Uint meth;
    std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
    time.Reset();

    double outer_tol;
    std::cerr << "tol = "; std::cin >> outer_tol;

    {
        std::cerr << "Computing initial velocity..." << std::endl;
        PSchur_GSPCG_CL schurSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);
        schurSolver.Solve( A->Data, B->Data, v->Data, p->Data, b->Data, c->Data);
    }    

    EnsightP2SolOutCL ensight( MG, lidx);
    
    const char datgeo[]= "ensight/shear.geo", 
               datpr[] = "ensight/shear.pr",
               datvec[]= "ensight/shear.vec",
               datscl[]= "ensight/shear.scl";
    ensight.CaseBegin( "shear.case", num_steps+1);
    ensight.DescribeGeom( "shear flow field", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true); 
    ensight.DescribeScalar( "Pressure", datpr,  true); 
    ensight.DescribeVector( "Velocity", datvec, true); 
    ensight.putGeom( datgeo);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);

    if (meth)
    {
//            PSchur_PCG_CL schurSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);
        PSchur_GSPCG_CL schurSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);
        CouplLevelsetStokesCL<StokesProblemT, PSchur_GSPCG_CL> 
            cpl( Stokes, lset, schurSolver);
        cpl.SetTimeStep( delta_t);

        for (Uint step= 1; step<=num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( FPsteps);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*delta_t);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*delta_t);
            ensight.putScalar( datscl, lset.GetSolution(), step*delta_t);
        }
    }
    else // Uzawa
    {
        double tau;
        Uint inner_iter;
        tau=  0.5*delta_t;
        std::cerr << "#PCG steps = "; std::cin >> inner_iter;
        Uzawa_PCG_CL uzawaSolver( prM.Data, 5000, outer_tol, inner_iter, inner_iter_tol, tau);
        CouplLevelsetStokesCL<StokesProblemT, Uzawa_PCG_CL> 
            cpl( Stokes, lset, uzawaSolver);
        cpl.SetTimeStep( delta_t);

        for (Uint step= 1; step<=num_steps; ++step)
        {
            std::cerr << "============= Schritt " << step << ":\n";
            cpl.DoStep( FPsteps);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*delta_t);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*delta_t);
            ensight.putScalar( datscl, lset.GetSolution(), step*delta_t);
        }
        std::cerr << "Iterationen: " << uzawaSolver.GetIter()
                  << "\tNorm des Res.: " << uzawaSolver.GetResid() << std::endl;
    }

    ensight.CaseEnd();
    
    std::cerr << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=4)
    {
        std::cerr << "You have to specify three parameters:\n\tlsshear <inner_iter_tol> <num_subdiv> <surf.tension>" << std::endl;
        return 1;
    }

    double inner_iter_tol= std::atof(argv[1]);
    int sub_div= std::atoi(argv[2]);
    double sigma= std::atof(argv[3]);
    std::cerr << "inner iter tol:  " << inner_iter_tol << std::endl;
    std::cerr << "sub divisions:   " << sub_div << std::endl;
    std::cerr << "surface tension: " << sigma << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    typedef DROPS::StokesP2P1CL<ShearFlowCL> MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, sub_div, sub_div, sub_div);

    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null,  &Parabol, &Parabol }; 
    // parabol. Einstroembedingungen bei z=0 und z=1 
        
    MyStokesCL prob(brick, ShearFlowCL(), DROPS::StokesBndDataCL(24, IsNeumann, bnd_fun));
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
