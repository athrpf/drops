//**************************************************************************
// File:    risingBubble.cpp                                               *
// Content: gravity driven flow of a rising bubble                         *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include <fstream>

double      delta_t= 0.1;
DROPS::Uint num_steps= 5;
const int   FPsteps= -1;

// rho*du/dt - mu/Re*laplace u + Dp = f + rho*g - okn
//                          -div u = 0
//                               u = u0, t=t0

// Tropfendaten:
DROPS::Point3DCL Mitte(0.5);
double           Radius= 0.25;

// Glaettungszone fuer Dichte-/Viskositaetssprung
const double sm_eps= 0.05;


class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    DROPS::SmoothedJumpCL rho, mu;
    const double Re, We;
    DROPS::Point3DCL g;

    ZeroFlowCL() 
      : rho( DROPS::JumpCL( 1, 10), DROPS::H_sm, sm_eps),
         mu( DROPS::JumpCL( 1, 2), DROPS::H_sm, sm_eps),
        Re(1), We(1) 
    { g[2]= -9.81; }
};


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)   { return DROPS::SVectorCL<3>(0.); }

double DistanceFct( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= Mitte-p;
    return d.norm()-Radius;
}


namespace DROPS // for Strategy
{


template<class Coeff>
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes, double inner_iter_tol, double sigma)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    LevelsetP2CL lset( MG, sigma, 0.5, 0.1);
    
    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    MatDescCL prM;

    TimerCL time;
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);    
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);    
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    MG.SizeInfo( std::cerr);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";
    Stokes.A.SetIdx(vidx, vidx);
    Stokes.B.SetIdx(pidx, vidx);
    Stokes.M.SetIdx(vidx, vidx);
    prM.SetIdx( pidx, pidx);
    
    Stokes.InitVel( &Stokes.v, Null);

    Stokes.SetupPrMass( &prM);

    Uint meth;
    std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
    time.Reset();

    double outer_tol;
    std::cerr << "tol = "; std::cin >> outer_tol;

    EnsightP2SolOutCL ensight( MG, lidx);
    
    const char datgeo[]= "ensight/risebubble.geo", 
               datpr[] = "ensight/risebubble.pr",
               datvec[]= "ensight/risebubble.vec",
               datscl[]= "ensight/risebubble.scl";
    ensight.CaseBegin( "risebubble.case", num_steps+1);
    ensight.DescribeGeom( "zero flow", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true); 
    ensight.DescribeScalar( "Pressure", datpr,  true); 
    ensight.DescribeVector( "Velocity", datvec, true); 
    ensight.putGeom( datgeo);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    if (meth)
    {
        PSchur_GSPCG_CL schurSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);
        CouplLevelsetStokes2PhaseCL<StokesProblemT, PSchur_GSPCG_CL> 
            cpl( Stokes, lset, schurSolver);
        
        cpl.SetTimeStep( delta_t);

        for (Uint step= 1; step<=num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( FPsteps);
//            if ((step%10)==0) lset.Reparam( 5, 0.01);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*delta_t);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*delta_t);
            ensight.putScalar( datscl, lset.GetSolution(), step*delta_t);

//            Stokes.SetupPrMass( &prM, lset);
        }
    }
    else // Uzawa
    {
        double tau;
        Uint inner_iter;
        tau=  0.5*delta_t;
        std::cerr << "#PCG steps = "; std::cin >> inner_iter;
        Uzawa_PCG_CL uzawaSolver( prM.Data, 5000, outer_tol, inner_iter, inner_iter_tol, tau);
        CouplLevelsetStokes2PhaseCL<StokesProblemT, Uzawa_PCG_CL> 
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


void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(Radius,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(Radius,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

int main (int argc, char** argv)
{
  try
  {
    if (argc<4)
    {
        std::cerr << "You have to specify at least three parameters:\n\t" 
                  << argv[0] << " <inner_iter_tol> <num_subdiv> <surf.tension> [<dt> <num_steps>]" << std::endl;
        return 1;
    }
    double inner_iter_tol= atof(argv[1]);
    int sub_div= atoi(argv[2]);
    double sigma= atof(argv[3]);
    if (argc>4) delta_t= atof(argv[4]);
    if (argc>5) num_steps= atoi(argv[5]);

    std::cerr << "inner iter tol:  " << inner_iter_tol << std::endl;
    std::cerr << "sub divisions:   " << sub_div << std::endl;
    std::cerr << "surface tension: " << sigma << std::endl;
    std::cerr << num_steps << " time steps of size " << delta_t << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= 1.; e3[2]= 2.;

    typedef DROPS::InstatStokes2PhaseP2P1CL<ZeroFlowCL> 
            StokesOnBrickCL;
    typedef StokesOnBrickCL MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, sub_div, sub_div, DROPS::Uint(sub_div*e3[2]));

    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::InstatStokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Null }; 
        
    StokesOnBrickCL prob(brick, ZeroFlowCL(), DROPS::InstatStokesBndDataCL(6, IsNeumann, bnd_fun));
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
