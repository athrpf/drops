//**************************************************************************
// File:    mzelle_instat.cpp                                              *
// Content: flow in drop cell                                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "out/output.h"
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

// rho*du/dt - mu/Re*laplace u + Dp = f + rho*g - okn
//                          -div u = 0
//                               u = u0, t=t0

class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double Re, We;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleCL& c) 
      : rho( DROPS::JumpCL( c.rhoD, c.rhoF ), DROPS::H_sm, c.sm_eps),
         mu( DROPS::JumpCL( c.muD,  c.muF),   DROPS::H_sm, c.sm_eps),
        Re(1.), We(1.), g( c.g)    {}

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
          innerSolver_( SSORPcCL( 1.), 10, inner_iter, inner_tol),
          outerSolver_( DummyPcCL(), 10, outer_iter, outer_tol)
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
          innerSolver_( SSORPcCL( 1.), 10, inner_iter, inner_tol),
          outerSolver_( Spc, 10, outer_iter, outer_tol)
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
    LevelsetP2CL lset( MG, C.sigma, C.theta, C.lset_SD); 

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
    lset.Init( DistanceFct);
    const double Vol= 4./3.*M_PI*std::pow(C.Radius,3);
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
    

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

    Stokes.SetupPrMass(  &prM);
    Stokes.SetupPrStiff( &prA);
    MatrixCL prM_A;
    ISPreCL ispc( prA.Data, prM.Data, C.theta*C.dt*C.muF/C.rhoF);
   
    lset.GetSolver().SetTol( C.lset_tol);
    lset.GetSolver().SetMaxIter( C.lset_iter);

    PSchur_PCG_CL schurSolver( prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
//    Schur_GMRes_CL schurSolver( C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

    // solve stationary problem for initial velocities    
    TimerCL time;
    VelVecDescCL curv( vidx);
    time.Reset();
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
    Stokes.SetupSystem2( &Stokes.B, &Stokes.c, Stokes.t);
    curv.Clear();
    lset.AccumulateBndIntegral( curv);
    time.Stop();
    std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

    double theta= C.stat_theta, nl= C.stat_nonlinear;
//    std::cerr << "theta = "; std::cin >> theta;
//    std::cerr << "nonlinear = "; std::cin >> nl;
    time.Reset();
    do
    {
        Stokes.SetupNonlinear( &Stokes.N, &Stokes.v, &cplN, lset, Stokes.t);
        MatrixCL mat;
        mat.LinComb( 1, Stokes.A.Data, nl*theta, Stokes.N.Data);
        cplN.Data-= (1-theta)*nl * (Stokes.N.Data * Stokes.v.Data);
        schurSolver.Solve( mat, Stokes.B.Data, 
            Stokes.v.Data, Stokes.p.Data, Stokes.b.Data + curv.Data + nl*cplN.Data, Stokes.c.Data);
    } while (schurSolver.GetIter() > 0);
    time.Stop();
    std::cerr << "Solving NavStokes for initial velocities took "<<time.GetTime()<<" sec.\n";

    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    if (C.scheme)
    {
//        Schur_GMRes_CL ISPschurSolver( C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
        ISPSchur_GMRes_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
//        ISPSchur_PCG_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

//        CouplLevelsetNavStokes2PhaseCL<StokesProblemT, Schur_GMRes_CL> 
        CouplLevelsetNavStokes2PhaseCL<StokesProblemT, ISPSchur_GMRes_CL> 
//        CouplLevelsetNavStokes2PhaseCL<StokesProblemT, ISPSchur_PCG_CL> 
            cpl( Stokes, lset, ISPschurSolver, C.theta, C.nonlinear);

        cpl.SetTimeStep( C.dt);

        for (int step= 1; step<=C.num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( C.FPsteps);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
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
                lset.Reparam( C.RepSteps, C.RepTau);
                ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.1)*C.dt);
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
    }
    else // Baensch scheme
    {
        ISPSchur_PCG_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

        CouplLsNsBaenschCL<StokesProblemT, ISPSchur_PCG_CL> 
            cpl( Stokes, lset, ISPschurSolver, C.nonlinear);

        cpl.SetTimeStep( C.dt);

        for (int step= 1; step<=C.num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
            cpl.DoStep( C.FPsteps);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
            ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
            ensight.Commit();
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

            if (C.RepFreq && step%C.RepFreq==0)
            {
                lset.Reparam( C.RepSteps, C.RepTau);
                ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.1)*C.dt);
                ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.1)*C.dt);
                ensight.putScalar( datscl, lset.GetSolution(), (step+0.1)*C.dt);
                ensight.Commit();
                std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
        }
    }
    ensight.CaseEnd();
    std::cerr << std::endl;
}

} // end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel= ~0)
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
        param.open( "NMRmzi.param");
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<ZeroFlowCL>    MyStokesCL;

    std::ifstream meshfile( C.meshfile.c_str());
    if (!meshfile)
    {
        std::cerr << "error while opening mesh file " << C.meshfile << "\n";
        return 1;
    }
    DROPS::ReadMeshBuilderCL builder( meshfile);
    
    
    const DROPS::BndCondT bc[3]= 
        { DROPS::OutflowBC, DROPS::WallBC, DROPS::DirBC};
    //    bottom,           side,          top
    const DROPS::InstatStokesVelBndDataCL::bnd_val_fun bnd_fun[3]= 
        { &Null, &Null, &Inflow}; 
        
    MyStokesCL prob(builder, ZeroFlowCL(C), DROPS::InstatStokesBndDataCL( 3, bc, bnd_fun));

    DROPS::MultiGridCL& mg = prob.GetMG();
    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    
    for (DROPS::BndIdxT i=0, num= bnd.GetNumBndSeg(); i<num; ++i)
    {
        std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }
    
    for (int i=0; i<C.num_dropref; ++i)
    {
        MarkDrop( mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    Strategy( prob);    // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
