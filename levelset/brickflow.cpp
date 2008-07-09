//****************************************************************************
// File:    brickflow.cpp                                                    *
// Content: test case: one-phase flow in square pipe with oscillating inflow *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen      *
//****************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/mzelle_hdr.h"
#include <fstream>

DROPS::ParamMesszelleNsCL C;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = p[0]*(2*C.r_inlet-p[0]) / (C.r_inlet*C.r_inlet),
                 z = p[2]*(2*C.r_inlet-p[2]) / (C.r_inlet*C.r_inlet),
                 freq= 50, ampl= 0.1;

    ret[1]= x * z * C.Anstroem * (1-ampl*std::cos(2*M_PI*freq*t));  // Rohr
    return ret;
}

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

class Schur_GMRes_CL: public PSchurSolver2CL<GMResSolverCL<SSORPcCL>, GMResSolverCL<DummyPcCL> >
{
  public:
    typedef GMResSolverCL<SSORPcCL>   innerSolverT;
    typedef GMResSolverCL<DummyPcCL>  outerSolverT;

  private:
    innerSolverT innerSolver_;
    outerSolverT outerSolver_;

  public:
    Schur_GMRes_CL( int outer_iter, double outer_tol,
                    int inner_iter, double inner_tol)
        : PSchurSolver2CL<innerSolverT, outerSolverT>(
              innerSolver_, outerSolver_, outer_iter, outer_tol
          ),
          innerSolver_( SSORPcCL( 1.), 150, inner_iter, inner_tol, false),  // absolute tolerances
          outerSolver_( DummyPcCL(), 150, outer_iter, outer_tol, false)
        {}
};

class ISPSchur_GMRes_CL: public PSchurSolver2CL<GMResSolverCL<SSORPcCL>, GMResSolverCL<ISPreCL> >
{
  public:
    typedef GMResSolverCL<SSORPcCL> innerSolverT;
    typedef GMResSolverCL<ISPreCL>  outerSolverT;

  private:
    innerSolverT innerSolver_;
    outerSolverT outerSolver_;

  public:
    ISPSchur_GMRes_CL(ISPreCL& Spc, int outer_iter, double outer_tol,
                                    int inner_iter, double inner_tol)
        : PSchurSolver2CL<innerSolverT, outerSolverT>(
              innerSolver_, outerSolver_, outer_iter, outer_tol
          ),
          innerSolver_( SSORPcCL( 1.), 150, inner_iter, inner_tol, false),  // absolute tolerances
          outerSolver_( Spc, 150, outer_iter, outer_tol, false)
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

    sigma=Stokes.GetCoeff().SurfTens;
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    VecDescCL cplN;

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);

    EnsightP2SolOutCL ensight( MG, lidx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr, true);
    ensight.DescribeVector( "Velocity", datvec, true);
    ensight.putGeom( datgeo);

    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);

    MG.SizeInfo( std::cerr);
    Stokes.b.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    cplN.SetIdx( vidx);
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

    switch (C.IniCond)
    {
      case 1: case 2: // stationary flow with/without drop
      {
        PSchur_PCG_CL schurSolver( Stokes.prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    //    Schur_GMRes_CL schurSolver( C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

        if (C.IniCond==2) // stationary flow without drop
            lset.Init( &One);
        // solve stationary problem for initial velocities
        TimerCL time;
        VelVecDescCL curv( vidx);
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";
        time.Reset();
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
                Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";

        if (C.IniCond==2)
            lset.Init( EllipsoidCL::DistanceFct);
      } break;

      case 3: // read from file
      {
        ReadEnsightP2SolCL reader( MG);
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
      } break;

      default:
        lset.Init( EllipsoidCL::DistanceFct);
    }

    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    if (C.scheme)
    {
        // PC for instat. Schur complement
        typedef ISPreCL SPcT;
        SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, 1./C.dt, C.theta);
//        typedef ISNonlinearPreCL SPcT;
//        SPcT ispc( prA.Data, prM.Data, 1./C.dt, C.theta);

        // PC for A-Block-PC
//      typedef  DummyPcCL APcPcT;
        typedef JACPcCL  APcPcT;
//        typedef SSORPcCL APcPcT;
//        typedef GSPcCL   APcPcT;
        APcPcT Apcpc;

        // PC for A-block
        typedef BiCGStabSolverCL<APcPcT> ASolverT;        // BiCGStab-based APcT
        ASolverT Asolver( Apcpc, 500, 0.02, /*relative=*/ true);
//        typedef GMResSolverCL<APcPcT>    ASolverT;        // GMRes-based APcT
//        ASolverT Asolver( Apcpc, 500, /*restart*/ 20, 0.02, /*relative=*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver);

        // PC for Oseen solver
//        typedef DummyPcCL OseenPcT;
//        OseenPcT oseenpc;
        typedef BlockPreCL<APcT, SPcT> OseenPcT;
        OseenPcT oseenpc( Apc, ispc);

        // Oseen solver
//        typedef InexactUzawaCL<APcT, SPcT, APC_OTHER> OseenSolverT;
//        OseenSolverT oseensolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.1);
        typedef GCRSolverCL<OseenPcT> OseenBaseSolverT;
        OseenBaseSolverT oseensolver0( oseenpc, /*truncate*/ 50, C.outer_iter, C.outer_tol, /*relative*/ false);
        typedef BlockMatrixSolverCL<OseenBaseSolverT> OseenSolverT;
        OseenSolverT oseensolver( oseensolver0);

        // Navier-Stokes solver
        typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
        NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red);

        FracStepScheme2PhaseCL<StokesProblemT, NSSolverT>
            cpl( Stokes, lset, nssolver, C.nonlinear);

        cpl.SetTimeStep( C.dt);

        for (int step= 1; step<=C.num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( C.cpl_iter);
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr)
            {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
            ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
            ensight.Commit();

            if (C.RepFreq && step%C.RepFreq==0) // reparam levelset function
            {
                lset.ReparamFastMarching( C.RepMethod);
                std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                if (C.VolCorr)
                {
                    double dphi= lset.AdjustVolume( Vol, 1e-9);
                    std::cerr << "volume correction is " << dphi << std::endl;
                    lset.Phi.Data+= dphi;
                    std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                }
                ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.1)*C.dt);
                ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.1)*C.dt);
                ensight.putScalar( datscl, lset.GetSolution(), (step+0.1)*C.dt);
                ensight.Commit();
            }
        }
    }
    else // Baensch scheme
    {
        // PC for instat. Schur complement
        const double theta= 1-std::sqrt(2.)/2, frac_dt= C.dt*theta, alpha= (1-2*theta)/(1-theta);
        ISPreCL ispc( Stokes.prA.Data, Stokes.prM.Data, 1./frac_dt, alpha);

        // PC for A-block
        typedef PCG_SsorCL ASolverT;
        ASolverT Asolver( SSORPcCL( 1.0), C.inner_iter, 0.01, /*relative*/true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT pcg( Asolver);

        // Stokes solver
//        ISPSchur_PCG_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
        typedef InexactUzawaCL<APcT, ISPreCL, APC_SYM> InexactUzawaT;
        InexactUzawaT inexactUzawaSolver( pcg, ispc, C.outer_iter, C.outer_tol);

//        CouplLsNsBaenschCL<StokesProblemT, ISPSchur_PCG_CL>
//            cpl( Stokes, lset, ISPschurSolver, C.inner_iter, C.inner_tol, C.nonlinear);
        OperatorSplitting2PhaseCL<StokesProblemT, InexactUzawaT>
            cpl( Stokes, lset, inexactUzawaSolver, C.inner_iter, C.inner_tol, C.nonlinear);

        cpl.SetTimeStep( C.dt);

        for (int step= 1; step<=C.num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( C.cpl_iter);
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr)
            {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
            ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
            ensight.Commit();

            if (C.RepFreq && step%C.RepFreq==0)
            {
                lset.ReparamFastMarching( C.RepMethod);
                std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                if (C.VolCorr)
                {
                    double dphi= lset.AdjustVolume( Vol, 1e-9);
                    std::cerr << "volume correction is " << dphi << std::endl;
                    lset.Phi.Data+= dphi;
                    std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                }
                ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.1)*C.dt);
                ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.1)*C.dt);
                ensight.putScalar( datscl, lset.GetSolution(), (step+0.1)*C.dt);
                ensight.Commit();
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
    if (argc!=2)
    {
        std::cerr << "You have to specify one parameter:\n\t"
                  << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::Point3DCL e1(0.), e2(0.), e3(0.), orig;
    e1[0]=2*C.r_inlet; e2[1]=0.03; e3[2]= 2*C.r_inlet;
    DROPS::BrickBuilderCL builder(orig,e1,e2,e3,5,30,5);

    const bool bc[6]=
        {false, false, true, false, false, false};  // Rohr
//        {false, false, true, false, true, true};      // Kanal
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &Inflow, &DROPS::ZeroVel, &DROPS::ZeroVel};

    DROPS::MultiGridCL mg( builder);
    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

    EllipsoidCL::Init( C.Mitte, C.Radius );
    for (int i=0; i<C.ref_flevel; ++i)
    {
        MarkAll( mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun), DROPS::P1_FE, C.XFEMStab);

    Strategy( prob);    // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
