//**************************************************************************
// File:    surfTens.cpp                                                   *
// Content: effect of surface tension                                      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "num/nssolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/mzelle_hdr.h"
#include <fstream>


DROPS::ParamMesszelleCL C;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

namespace DROPS // for Strategy
{

class PSchur_PCG_CL: public PSchurSolverCL<PCG_SsorCL>
{
  private:
    PCG_SsorCL _PCGsolver;
  public:
    PSchur_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double omega= 1.)
        : PSchurSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol),
          _PCGsolver(SSORPcCL(omega), inner_iter, inner_tol)
        {}
};

class ISPSchur_PCG_CL: public PSchurSolver2CL<PCGSolverCL<SSORPcCL>, PCGSolverCL<ISPreCL> >
{
  public:
    typedef PCGSolverCL<SSORPcCL> innerSolverT;
    typedef PCGSolverCL<ISPreCL>  outerSolverT;

  private:
    innerSolverT innerSolver_;
    outerSolverT outerSolver_;

  public:
    ISPSchur_PCG_CL(ISPreCL& Spc, int outer_iter, double outer_tol,
                                  int inner_iter, double inner_tol)
        : PSchurSolver2CL<innerSolverT, outerSolverT>(
              innerSolver_, outerSolver_, outer_iter, outer_tol
          ),
          innerSolver_( SSORPcCL( 1.), inner_iter, inner_tol),
          outerSolver_( Spc, outer_iter, outer_tol)
         {}
};

template<class Coeff>
void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes)
// flow control
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    // Levelset-Disc.: Crank-Nicholson
    sigma= C.sigma;
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);

    EnsightP2SolOutCL ensight( MG, lidx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datcrv= filename+".crv",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr,  true);
    ensight.DescribeVector( "Velocity", datvec, true);
    ensight.DescribeVector( "Curvature", datcrv, true);
    ensight.putGeom( datgeo);

    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);
    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

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
    Stokes.N.SetIdx(vidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx);
    Stokes.prA.SetIdx( pidx, pidx);

    Stokes.InitVel( &Stokes.v, ZeroVel);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    PSchur_PCG_CL   schurSolver( Stokes.prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    VelVecDescCL curv( vidx);

    if (C.IniCond != 0)
    {
        // solve stationary problem for initial velocities
        TimerCL time;
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        time.Reset();
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
            Stokes.v.Data, Stokes.p.Data, VectorCL( Stokes.b.Data + curv.Data), Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
    }
    curv.Clear();
    lset.AccumulateBndIntegral( curv);

    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putVector( datcrv, Stokes.GetVelSolution( curv), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    ISPreCL ispc( Stokes.prA.Data, Stokes.prM.Data, 1./C.dt, C.theta);
    ISPSchur_PCG_CL ISPschurSolver( ispc,  C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    ISPschurSolver.SetTol( C.outer_tol);

    typedef NSSolverBaseCL<StokesProblemT> SolverT;
    SolverT navstokessolver(Stokes, ISPschurSolver);

    LinThetaScheme2PhaseCL<StokesProblemT, SolverT>
        cpl( Stokes, lset, navstokessolver, C.theta, 0.);
    cpl.SetTimeStep( C.dt);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( C.cpl_iter);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);

        ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
        ensight.putVector( datcrv, Stokes.GetVelSolution( curv), step*C.dt);
        ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
        ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
        ensight.Commit();
        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        if (C.VolCorr)
        {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        }

        if (C.RepFreq && step%C.RepFreq==0)
        {
            lset.ReparamFastMarching( C.RepMethod);
            curv.Clear();
            lset.AccumulateBndIntegral( curv);

            ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.1)*C.dt);
            ensight.putVector( datcrv, Stokes.GetVelSolution( curv), (step+0.1)*C.dt);
            ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.1)*C.dt);
            ensight.putScalar( datscl, lset.GetSolution(), (step+0.1)*C.dt);
            ensight.Commit();
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr)
            {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
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
    if (argc>2)
    {
        std::cerr << "You have to specify at most one parameter:\n\t"
                  << argv[0] << " [<param_file>]" << std::endl;
        return 1;
    }
    std::ifstream param;
    if (argc>1)
        param.open( argv[1]);
    else
        param.open( "surfTens.param");
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<ZeroFlowCL>    MyStokesCL;

    const double L= 3e-3;
    DROPS::Point3DCL orig(-L), e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 2*L;

    const int n= std::atoi( C.meshfile.c_str());
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, n, n, n);

    const DROPS::BndCondT bc[6]=
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel};

    MyStokesCL prob(builder, ZeroFlowCL(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun));

    DROPS::MultiGridCL& mg = prob.GetMG();

    EllipsoidCL::Init( C.Mitte, C.Radius);
    for (int i=0; i<C.ref_flevel; ++i)
    {
        DROPS::MarkInterface( EllipsoidCL::DistanceFct, C.ref_width, mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    Strategy( prob);  // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
