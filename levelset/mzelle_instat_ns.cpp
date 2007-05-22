//**************************************************************************
// File:    mzelle_instat.cpp                                              *
// Content: flow in drop cell                                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include <fstream>

DROPS::ParamMesszelleNsCL C;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( C.rhoD, C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muD,  C.muF),   DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( 1., C.rhoF/C.rhoD ), DROPS::H_sm, C.sm_eps),
        mu ( DROPS::JumpCL( 1., C.muF/C.muD),    DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma/C.rhoD), g( C.g)    {}
};


DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

double DistanceFct( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= C.Mitte-p;
    return d.norm()-C.Radius;
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } 


namespace DROPS // for Strategy
{

template<class Coeff>
void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes)
// flow control
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();

    sigma= Stokes.GetCoeff().SurfTens;
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.lset_theta, C.lset_SD, C.RepDiff, C.lset_iter, C.lset_tol, C.CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    VecDescCL cplN;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);

    EnsightP2SolOutCL ensight( MG, lidx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr,  true);
    ensight.DescribeVector( "Velocity", datvec, true);
    ensight.putGeom( datgeo);


    MG.SizeInfo( std::cerr);
    Stokes.b.SetIdx( vidx);
    Stokes.v.SetIdx( vidx);
    cplN.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";
    Stokes.A.SetIdx(vidx, vidx);
    Stokes.B.SetIdx(pidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx);
    Stokes.prA.SetIdx( pidx, pidx);
    Stokes.M.SetIdx(vidx, vidx);
    Stokes.N.SetIdx(vidx, vidx);

    Stokes.InitVel( &Stokes.v, &ZeroVel);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    switch (C.IniCond)
    {
      case 1: case 2: // stationary flow with/without drop
      {
        PSchur_PCG_CL schurSolver( Stokes.prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    //    Schur_GMRes_CL schurSolver( C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

        const double old_Radius= C.Radius;
        if (C.IniCond==2) // stationary flow without drop
            C.Radius= -10;
        lset.Init( DistanceFct);

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

        double theta= C.stat_theta, nl= C.stat_nonlinear;
        time.Reset();
        do
        {
            Stokes.SetupNonlinear( &Stokes.N, &Stokes.v, &cplN, lset, Stokes.t);
            MatrixCL mat;
            mat.LinComb( 1, Stokes.A.Data, nl*theta, Stokes.N.Data);
            cplN.Data-= (1-theta) * (Stokes.N.Data * Stokes.v.Data);
            schurSolver.Solve( mat, Stokes.B.Data,
                Stokes.v.Data, Stokes.p.Data, VectorCL( Stokes.b.Data + nl*cplN.Data), Stokes.c.Data);
        } while (schurSolver.GetIter() > 0);
        time.Stop();
        std::cerr << "Solving NavStokes for initial velocities took "<< time.GetTime()<<" sec.\n";

        if (C.IniCond==2)
        {
            C.Radius= old_Radius;
            lset.Init( DistanceFct);
        }
      } break;

      case 3: // read from file
      { // FIXME: only the P1-part of P1X-elements can be read
        ReadEnsightP2SolCL reader( MG);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset, /*NumberingChanged*/ false);
        Stokes.p.SetIdx( pidx); // Zero-vector for now.
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr); // reads the P1-part of the pressure
      } break;

      default:
        lset.Init( DistanceFct);
    }
    const double Vol= 4./3.*M_PI*std::pow(C.Radius,3);
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    if (C.scheme)
    {
        // PC for instat. Schur complement
        typedef ISPreCL SPcT;
        SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM*/ C.theta);
//        typedef ISNonlinearPreCL SPcT;
//        SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM*/ C.theta);

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
//        typedef DiagBlockPreCL<APcT, SPcT> OseenPcT;
//        OseenPcT oseenpc( Apc, ispc);

        // Oseen solver
        typedef InexactUzawaCL<APcT, SPcT, APC_OTHER> OseenSolverT;
        OseenSolverT oseensolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.1);
//        typedef GCRSolverCL<OseenPcT> OseenBaseSolverT;
//        OseenBaseSolverT oseensolver0( oseenpc, /*truncate*/ 50, C.outer_iter, C.outer_tol, /*relative*/ false);
//        typedef BlockMatrixSolverCL<OseenBaseSolverT> OseenSolverT;
//        OseenSolverT oseensolver( oseensolver0);

        // Navier-Stokes solver
        typedef AdaptFixedPtDefectCorrCL<StokesProblemT, OseenSolverT> NSSolverT;
        NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red);

        // Time-integration and coupling
//        typedef CouplLevelsetNavStokes2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//        CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear);
        typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
        CouplingT cpl( Stokes, lset, nssolver, C.nonlinear);

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
                if (C.RepMethod>1)
                    lset.Reparam( C.RepSteps, C.RepTau);
                else
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
        ISPreCL ispc( Stokes.prA.Data, Stokes.prM.Data, C.dt*(3 - 2*std::sqrt(2.)));

        // PC for A-block
        typedef PCG_SsorCL ASolverT;
        ASolverT Asolver( SSORPcCL( 1.0), C.inner_iter, 0.01, /*relative*/true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver);

        // Stokes solver
//        ISPSchur_PCG_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
        typedef InexactUzawaCL<APcT, ISPreCL, APC_SYM> InexactUzawaT;
        InexactUzawaT inexactUzawaSolver( Apc, ispc, C.outer_iter, C.outer_tol);

//        CouplLsNsBaenschCL<StokesProblemT, ISPSchur_PCG_CL>
//            cpl( Stokes, lset, ISPschurSolver, C.inner_iter, C.inner_tol, C.nonlinear);
        CouplLsNsBaenschCL<StokesProblemT, InexactUzawaT>
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
                if (C.RepMethod>1)
                    lset.Reparam( C.RepSteps, C.RepTau);
                else
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

    std::ifstream meshfile( C.meshfile.c_str());
    if (!meshfile)
    {
        std::cerr << "error while opening mesh file " << C.meshfile << "\n";
        return 1;
    }
    DROPS::ReadMeshBuilderCL builder( meshfile);
    DROPS::MultiGridCL mg( builder);
    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

    if (num_bnd>10) { std::cerr << "Increase size of BndSegs in main() for proper use!\n"; return 1; }
    DROPS::BndCondT bc[10];
    DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[10];
    for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
    {
        bnd_fun[i]= (bc[i]= builder.GetBC( i))==DROPS::DirBC ? &Inflow : &DROPS::ZeroVel;
        std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }

    for (int i=0; i<C.num_dropref; ++i)
    {
        DROPS::MarkInterface( DistanceFct, 0.5*C.Radius, mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL(
        num_bnd, bc, bnd_fun), DROPS::P1_FE);

    Strategy( prob);    // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
