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
#include "levelset/adaptriang.h"
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


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = p[0]*(2*C.r_inlet-p[0]) / (C.r_inlet*C.r_inlet),
                 z = p[2]*(2*C.r_inlet-p[2]) / (C.r_inlet*C.r_inlet);
    ret[1]= x * z * C.Anstroem;
    return ret;
}

double DistanceFct( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.Mitte-p;
//    return d.norm()-C.Radius;
    const double avgRad= cbrt(C.Radius[0]*C.Radius[1]* C.Radius[2]);
    d/= C.Radius/avgRad;
    return d.norm()-avgRad;
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
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    VecDescCL cplN;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);
    lset.SetSurfaceForce( SF_ImprovedLB);

    double maxGradPhi, Volume;
    Point3DCL bary_drop, min_drop, max_drop;
    std::ofstream infofile( string(C.EnsCase + ".info").c_str());
    lset.GetInfo( maxGradPhi, Volume, bary_drop, min_drop, max_drop);
    infofile << Stokes.t << '\t' << maxGradPhi << '\t' << Volume << '\t' << bary_drop << '\t' << min_drop << '\t' << max_drop << std::endl; 

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);

    EnsightP2SolOutCL ensight( MG, lidx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo, true);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr,  true);
    ensight.DescribeVector( "Velocity", datvec, true);

    AdapTriangCL adap( MG, C.ref_width, 0, C.ref_flevel);

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

    Stokes.InitVel( &Stokes.v, Null);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    switch (C.IniCond)
    {
      case 1: case 2: // stationary flow with/without drop
      {
          typedef CGSolverCL SPcSolverT;
          SPcSolverT CGsolver( 50, 0.02, /*relative*/ true);
          typedef ISNonlinearPreCL<SPcSolverT> SPcT;
          SPcT ispc( CGsolver, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 0., /*kM*/ 1.);
//            typedef ISPreCL SPcT;
//            SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 0., /*kM*/ 1.);

        // PC for A-Block-PC
//          typedef  DummyPcCL APcPcT;
          typedef JACPcCL  APcPcT;
//          typedef SSORPcCL APcPcT;
//          typedef GSPcCL   APcPcT;
          APcPcT Apcpc;

        // PC for A-block
          typedef PCGSolverCL<APcPcT> ASolverT;        // CG-based APcT
          ASolverT Asolver( Apcpc, 500, 0.02, /*relative=*/ true);
//          typedef GMResSolverCL<APcPcT>    ASolverT;        // GMRes-based APcT
//          ASolverT Asolver( Apcpc, 500, /*restart*/ 100, 0.02, /*relative=*/ true);
          typedef SolverAsPreCL<ASolverT> APcT;
          APcT Apc( Asolver);

        // PC for Oseen solver
//        typedef DummyPcCL OseenPcT;
//        OseenPcT oseenpc;
//        typedef BlockPreCL<APcT, SPcT> OseenPcT;
//        OseenPcT oseenpc( Apc, ispc);

        // Oseen solver
          typedef InexactUzawaCL<APcT, SPcT, APC_SYM> OseenSolverT;
          OseenSolverT schurSolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.1);

//        PSchur_PCG_CL schurSolver( Stokes.prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    //    Schur_GMRes_CL schurSolver( C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

        const Point3DCL old_Radius= C.Radius;
        if (C.IniCond==2) // stationary flow without drop
            C.Radius= Point3DCL(-10.);
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
        time.Reset();
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
            Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<< time.GetTime()<<" sec.\n";

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
    const double Vol= 4./3.*M_PI*C.Radius[0]*C.Radius[1]*C.Radius[2];
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    ensight.putGeom( datgeo, 0);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    if (C.scheme)
    {
        // PC for instat. Schur complement
        typedef DummyPcCL SPcT;
        SPcT ispc;
//         typedef ISPreCL SPcT;
//         SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));
//         typedef ISPreJacCL SPcT;
//         SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));
//       typedef PCGSolverCL<JACPcCL> SPcSolverT;
//       SPcSolverT JacCG( JACPcCL(), 50, 0.02, /*relative*/ true);
//       typedef ISNonlinearPreCL<SPcSolverT> SPcT;
//       SPcT ispc( JacCG, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));
//         typedef CGSolverCL SPcSolverT;
//         SPcSolverT CGsolver( 50, 0.02, /*relative*/ true);
//         typedef ISNonlinearPreCL<SPcSolverT> SPcT;
//         SPcT ispc( CGsolver, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));

        // PC for A-Block-PC
//      typedef  DummyPcCL APcPcT;
        typedef JACPcCL  APcPcT;
//        typedef SSORPcCL APcPcT;
//        typedef GSPcCL   APcPcT;
        APcPcT Apcpc;

        // PC for A-block
//        typedef BiCGStabSolverCL<APcPcT> ASolverT;        // BiCGStab-based APcT
//        ASolverT Asolver( Apcpc, 500, 0.02, /*relative=*/ true);
        typedef GMResSolverCL<APcPcT>    ASolverT;        // GMRes-based APcT
        ASolverT Asolver( Apcpc, 500, /*restart*/ 100, 2e-5, /*relative=*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver);

        // PC for Oseen solver
//        typedef DummyPcCL OseenPcT;
//        OseenPcT oseenpc;
        typedef BlockPreCL<APcT, SPcT> OseenPcT;
        OseenPcT oseenpc( Apc, ispc);

        // Oseen solver
//        typedef InexactUzawaCL<APcT, SPcT, APC_OTHER> OseenSolverT;
//        OseenSolverT oseensolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.2);
       typedef GCRSolverCL<OseenPcT> OseenBaseSolverT;
       OseenBaseSolverT oseensolver0( oseenpc, /*truncate*/ C.outer_iter, C.outer_iter, C.outer_tol, /*relative*/ false);
       typedef BlockMatrixSolverCL<OseenBaseSolverT> OseenSolverT;
       OseenSolverT oseensolver( oseensolver0);

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
            lset.GetInfo( maxGradPhi, Volume, bary_drop, min_drop, max_drop);
            infofile << Stokes.t << '\t' << maxGradPhi << '\t' << Volume << '\t' << bary_drop << '\t' << min_drop << '\t' << max_drop << std::endl; 
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr)
            {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }

            if (C.RepFreq && step%C.RepFreq==0) // reparam levelset function
            {
                lset.ReparamFastMarching( C.RepMethod);
                
                adap.UpdateTriang( Stokes, lset);
                if (adap.WasModified() )
                {
                    cpl.Update();
                    // don't forget to update the pr mass/stiff matrix for the schur compl. preconditioner!!
                    Stokes.prM.SetIdx( pidx, pidx);
                    Stokes.SetupPrMass( &Stokes.prM, lset);
                    Stokes.prA.SetIdx( pidx, pidx);
                    Stokes.SetupPrStiff( &Stokes.prA, lset);
                }
                
                std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                
                if (C.VolCorr)
                {
                    double dphi= lset.AdjustVolume( Vol, 1e-9);
                    std::cerr << "volume correction is " << dphi << std::endl;
                    lset.Phi.Data+= dphi;
                    std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                }
            }
            ensight.putGeom( datgeo, step*C.dt);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
            ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
            ensight.Commit();
        }
    } 
    // no Baensch scheme (FracStep + OpSplitting) 
    
    ensight.CaseEnd();
    std::cerr << std::endl;
}

} // end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, int maxLevel= -1)
{
    const double avgRadius= cbrt( C.Radius[0]*C.Radius[1]*C.Radius[2]);
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (DistanceFct(GetBaryCenter(*It)) + avgRadius)<=std::max(1.5*avgRadius,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}


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

    int nx, ny, nz;
    double dx, dy, dz;
    std::string mesh( C.meshfile), delim("x@");
    size_t idx;
    while ((idx= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
    if (!brick_info || dx!=dz)
    {
        std::cerr << "error while reading geometry information: " << mesh << "\n";
        return 1;
    }
    C.r_inlet= dx/2;
    DROPS::Point3DCL orig, px, py, pz;
    px[0]= dx; py[1]= dy; pz[2]= dz;
    DROPS::BrickBuilderCL builder( orig, px, py, pz, nx, ny, nz);
    DROPS::MultiGridCL mg( builder);
    const DROPS::BndIdxT num_bnd= 6;

    DROPS::BndCondT bc[6];
    DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6];
    // y=0,dy: inflow, other: wall
    for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
    {
        const bool wall= i!=2 && i!=3;
        bc[i]= wall ? DROPS::WallBC : DROPS::DirBC;
        bnd_fun[i]= wall ? &Null : &Inflow;
    }

    for (int i=0; i<C.ref_flevel; ++i)
    {
        MarkDrop( mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL(
        num_bnd, bc, bnd_fun), DROPS::P1X_FE, C.XFEMStab);

    Strategy( prob);    // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}