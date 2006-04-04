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


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

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


namespace DROPS // for Strategy
{

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
    LevelsetP2CL lset( MG, Stokes.GetCoeff().SurfTens, C.lset_theta, C.lset_SD, C.RepDiff, C.lset_iter, C.lset_tol, C.CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    MatDescCL prM, prA;
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
    ensight.DescribeScalar( "Pressure", datpr,  true); 
    ensight.DescribeVector( "Velocity", datvec, true); 
    ensight.putGeom( datgeo);

    lset.Phi.SetIdx( lidx);

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
    prM.SetIdx( pidx, pidx);
    prA.SetIdx( pidx, pidx);
    
    Stokes.InitVel( &Stokes.v, Null);
    Stokes.SetupPrMass(  &prM, lset);
    Stokes.SetupPrStiff( &prA, lset);
   
    switch (C.IniCond)
    {
      case 1: case 2: // stationary flow with/without drop
      {
        PSchur_PCG_CL schurSolver( prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
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
        std::cerr << "Solving NavStokes for initial velocities took "<<time.GetTime()<<" sec.\n";

        if (C.IniCond==2)
        { 
            C.Radius= old_Radius;
            lset.Init( DistanceFct);
        }
      } break;
      
      case 3: // read from file
      {
        ReadEnsightP2SolCL reader( MG);
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
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
        ISPreCL ispc( prA.Data, prM.Data, C.theta*C.dt);

        // PC for A-Block-PC
//      typedef  DummyPcCL APcPcT;
        typedef JACPcCL APcPcT;
//        typedef SSORPcCL APcPcT;
//        typedef GSPcCL APcPcT;
        APcPcT Apcpc;

        // PC for A-block
        typedef BiCGStab_PreCL<APcPcT> APcT;        // BiCGStab-based APcT
        APcT Asolver( Apcpc, 500, 0.02);
//        typedef GMRes_PreCL<APcPcT> APcT;        // GMRes-based APcT
//        APcT Asolver( Apcpc, 500, /*restart*/ 20, 0.02);

        // Oseen solver
        typedef InexactUzawaCL<APcT, ISPreCL, APC_OTHER> InexactUzawaSolverT;
        InexactUzawaSolverT inexactUzawaSolver( Asolver, ispc, C.outer_iter, C.outer_tol, 0.1);

        // Navier-Stokes solver
        typedef AdaptFixedPtDefectCorrCL<StokesProblemT, InexactUzawaSolverT> NSSolverT;
        NSSolverT nssolver( Stokes, inexactUzawaSolver, C.ns_iter, C.ns_tol, C.ns_red);

        CouplLevelsetNavStokes2PhaseCL<StokesProblemT, NSSolverT> 
            cpl( Stokes, lset, nssolver, C.theta, C.nonlinear);

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
        ISPreCL ispc( prA.Data, prM.Data, C.dt*(3 - 2*std::sqrt(2.)));

        // PC for A-block
        SSORPCG_PreCL pcg( C.inner_iter, 0.01);
        
        // Stokes solver
//        ISPSchur_PCG_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
        typedef InexactUzawaCL<SSORPCG_PreCL, ISPreCL, APC_SYM> InexactUzawaT;
        InexactUzawaT inexactUzawaSolver( pcg, ispc, C.outer_iter, C.outer_tol);

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


void MarkDrop (DROPS::MultiGridCL& mg, int maxLevel= -1)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-C.Mitte).norm()<=std::max(1.5*C.Radius,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
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
        bnd_fun[i]= (bc[i]= builder.GetBC( i))==DROPS::DirBC ? &Inflow : &Null;
        std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }
    
    for (int i=0; i<C.num_dropref; ++i)
    {
        MarkDrop( mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun));

    Strategy( prob);    // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
