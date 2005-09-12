//**************************************************************************
// File:    surfTens.cpp                                                   *
// Content: effect of surface tension                                      *
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
#include "levelset/params.h"
#include <fstream>


DROPS::ParamMesszelleCL C;

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

class DimLessCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const double L, U;
    const DROPS::SmoothedJumpCL rho, mu;
    const double Re, We;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( double l, double u, const DROPS::ParamMesszelleCL& c) 
      : L( l), U( u),
        rho( DROPS::JumpCL( 1, c.rhoF/c.rhoD ), DROPS::H_sm, c.sm_eps/L),
         mu( DROPS::JumpCL( 1, c.muF/c.muD),   DROPS::H_sm, c.sm_eps/L),
        Re(c.rhoD*L*U/c.muD), We(c.rhoD*L*U*U/c.sigma), g( (L/U/U)*c.g)    {}
};


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

double DistanceFct( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= C.Mitte-p;
    return d.norm()-C.Radius;
}

namespace DROPS // for Strategy
{

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
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    // Levelset-Disc.: Crank-Nicholson
    LevelsetP2CL lset( MG, C.sigma, C.theta, C.lset_SD, C.RepDiff, C.lset_iter, C.lset_tol, C.CurvDiff); 

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    MatDescCL prM, prA;

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
    lset.Init( DistanceFct);
    const double Vol= 4./3.*M_PI*std::pow(C.Radius,3);
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
    prM.SetIdx( pidx, pidx);
    prA.SetIdx( pidx, pidx);
    
    Stokes.InitVel( &Stokes.v, Null);
    Stokes.SetupPrMass(  &prM, lset);
    Stokes.SetupPrStiff( &prA, lset);

    PSchur_PCG_CL   schurSolver( prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
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

    ISPreCL ispc( prA.Data, prM.Data, C.theta*C.dt);
    ISPSchur_PCG_CL ISPschurSolver( ispc,  C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    ISPschurSolver.SetTol( C.outer_tol);
    
    CouplLevelsetStokes2PhaseCL<StokesProblemT, ISPSchur_PCG_CL> 
        cpl( Stokes, lset, ISPschurSolver, C.theta);
    cpl.SetTimeStep( C.dt);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( C.FPsteps);
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
            lset.Reparam( C.RepSteps, C.RepTau);
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


void MarkDrop (DROPS::MultiGridCL& mg, int maxLevel= -1)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
            int neg= 0;
            for (int i=0; i<4; ++i)
            {
                const double val= DistanceFct( It->GetVertex(i)->GetCoord() );
                neg+= val<0; 
            }   
//            const double val= DistanceFct( GetBaryCenter(*It));
//            neg+= val<0; 
            
            if (neg!=0 && neg!=4)
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
        param.open( "surfTens.param");
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    typedef DROPS::InstatStokes2PhaseP2P1CL<ZeroFlowCL>    MyStokesCL;
    
    const double L= 3e-3;
    DROPS::Point3DCL orig(-L), e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 2*L;
    
    const int n= atoi( C.meshfile.c_str());
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, n, n, n);
    
    const DROPS::BndCondT bc[6]= 
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Null}; 
        
    MyStokesCL prob(builder, ZeroFlowCL(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun));

    DROPS::MultiGridCL& mg = prob.GetMG();
    
    for (int i=0; i<C.num_dropref; ++i)
    {
        MarkDrop( mg);
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
