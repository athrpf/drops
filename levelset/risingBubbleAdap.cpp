//**************************************************************************
// File:    risingBubbleAdap.cpp                                           *
// Content: gravity driven flow of a rising bubble, grid adaptivity        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/adaptriang.h"
#include <fstream>

double      delta_t= 0.05;
DROPS::Uint num_steps= 5;
const int   FPsteps= -1;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

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
    DROPS::Point3DCL g;
    const double SurfTens;

    ZeroFlowCL()
      : rho( DROPS::JumpCL( 1, 10), DROPS::H_sm, sm_eps),
         mu( DROPS::JumpCL( 2, 1), DROPS::H_sm, sm_eps), SurfTens(0.)
    { g[2]= -9.81; }
};

double DistanceFct( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= Mitte-p;
    return d.norm()-Radius;
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } 


namespace DROPS // for Strategy
{

template<class Coeff>
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adap, double inner_iter_tol)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, 1, 0.1); // impl. Euler, SD=0.1

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    TimerCL time;

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";
    Stokes.prM.SetIdx( pidx, pidx);
    Stokes.SetupPrMass( &Stokes.prM, lset);

    Stokes.InitVel( &Stokes.v, ZeroVel);
    lset.Init( DistanceFct);

    time.Reset();

    double outer_tol;
    std::cerr << "tol = "; std::cin >> outer_tol;

    lset.GetSolver().SetTol( 1e-14);
    lset.GetSolver().SetMaxIter( 50000);

    IdxDescCL ens_idx( 1, 1);
    lset.CreateNumbering( MG.GetLastLevel(), &ens_idx);
    EnsightP2SolOutCL ensight( MG, &ens_idx);

    const char datgeo[]= "ensight/risebubbleadap.geo",
               datpr[] = "ensight/risebubbleadap.pr",
               datvec[]= "ensight/risebubbleadap.vec",
               datscl[]= "ensight/risebubbleadap.scl";
    ensight.CaseBegin( "risebubbleadap.case", num_steps+1);
    ensight.DescribeGeom( "zero flow", datgeo,  true);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr,  true);
    ensight.DescribeVector( "Velocity", datvec, true);
    ensight.putGeom( datgeo, 0);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);

    PSchur_GSPCG_CL schurSolver( Stokes.prM.Data, 200, outer_tol, 200, inner_iter_tol);

    CouplLevelsetStokes2PhaseCL<StokesProblemT, PSchur_GSPCG_CL>
        cpl( Stokes, lset, schurSolver, 1.0); // impl. Euler
    cpl.SetTimeStep( delta_t);

    for (Uint step= 1; step<=num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( FPsteps);
        if (step%10 == 0)
            lset.ReparamFastMarching();
        ensight.putGeom( datgeo, step*delta_t);
        ensight.putScalar( datpr, Stokes.GetPrSolution(), step*delta_t);
        ensight.putVector( datvec, Stokes.GetVelSolution(), step*delta_t);
        ensight.putScalar( datscl, lset.GetSolution(), step*delta_t);
        ensight.Commit();
        if (step<num_steps) // omit in last step
        {
//            LevelsetP2CL::DiscSolCL sol= lset.GetSolution();
            lset.DeleteNumbering( &ens_idx);
            adap.UpdateTriang( Stokes, lset);
            lset.CreateNumbering( MG.GetLastLevel(), &ens_idx);
/*
        ensight.putGeom( datgeo, (step+0.01)*delta_t);
        ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.01)*delta_t);
        ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.01)*delta_t);
        ensight.putScalar( datscl, lset.GetSolution(), (step+0.01)*delta_t);
*/
            if (adap.WasModified() )
            {
                cpl.Update();
                // don't forget to update the pr mass matrix for the schur compl. preconditioner!!
                Stokes.prM.SetIdx( pidx, pidx);
                Stokes.SetupPrMass( &Stokes.prM, lset);
            }
        }
    }
    ensight.CaseEnd();

    std::cerr << std::endl;
}

} // end of namespace DROPS


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
    double inner_iter_tol= std::atof(argv[1]);
    int sub_div= std::atoi(argv[2]);
    sigma= std::atof(argv[3]);
    if (argc>4) delta_t= std::atof(argv[4]);
    if (argc>5) num_steps= std::atoi(argv[5]);

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
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel };

    StokesOnBrickCL prob(brick, ZeroFlowCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::AdapTriangCL adap( mg, 0.1, 0, 3);

    adap.MakeInitialTriang( DistanceFct);
    Strategy( prob, adap, inner_iter_tol);
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
