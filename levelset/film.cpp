//**************************************************************************
// File:    mzelle_instat.cpp                                              *
// Content: flow in drop cell                                              *
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
#include "levelset/params.h"
#include <fstream>


DROPS::ParamFilmCL C;

// rho*du/dt - mu/Re*laplace u + Dp = f + rho*g - okn
//                          -div u = 0
//                               u = u0, t=t0


class ZeroFlowCL
{
// \Omega_1 = Film,    \Omega_2 = Gasphase
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double Re, SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamFilmCL& C) 
      : rho( DROPS::JumpCL( C.rhoF, C.rhoG ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muF,  C.muG),   DROPS::H_sm, C.sm_eps),
        Re( 1.), SurfTens( C.sigma), g( C.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Film,    \Omega_2 = Gasphase
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double Re, SurfTens;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( const DROPS::ParamFilmCL& C) 
      : rho( DROPS::JumpCL( 1., C.rhoG/C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu ( DROPS::JumpCL( 1., C.muG/C.muF),    DROPS::H_sm, C.sm_eps),
        Re( C.rhoF/C.muF), SurfTens( C.sigma/C.rhoF), g( C.g)    {}
};


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double t)
{ 
    DROPS::SVectorCL<3> ret(0.); 
    const double d= p[1]/C.Filmdicke;
    static const double u= C.rhoF*C.g[0]*C.Filmdicke*C.Filmdicke/C.muF/2;
    ret[0]= d<=1 ? (2*d-d*d)*u * (1 + C.PumpAmpl*sin(2*M_PI*t*C.PumpFreq))
                 : (C.mesh_size[1]-p[1])/(C.mesh_size[1]-C.Filmdicke)*u;
    return ret; 
}

double DistanceFct( const DROPS::Point3DCL& p)
{
    // wave length = 100 x film width
    const double wave= C.PumpAmpl*sin(2*M_PI*p[0]/C.mesh_size[0]),
        z= p[2]/C.mesh_size[2]*2; // z \in [-1,1]
    return p[1] - C.Filmdicke * (1 + wave*cos(z*M_PI/2.));
}


namespace DROPS // for Strategy
{

bool periodic_x( const Point3DCL& p, const Point3DCL& q)
{ // matching y-z- or x-y-coords, resp.
    return std::abs( p[1]-q[1]) + std::abs( p[2]-q[2]) < 1e-12
      ||   std::abs( p[0]-q[0]) + std::abs( p[1]-q[1]) < 1e-12;
}

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
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    IdxDescCL ens_idx( 1, 1);
    MatDescCL prM, prA;

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx, periodic_x);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, periodic_x);    
    lset.CreateNumbering(      MG.GetLastLevel(), lidx, periodic_x);
    CreateNumb( MG.GetLastLevel(), ens_idx, MG, NoBndDataCL<>());
    
    lset.Phi.SetIdx( lidx);
    
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
    switch (C.IniCond)
    {
      case 1: // stationary flow
      {
        lset.Init( DistanceFct);
        TimerCL time;
        VelVecDescCL curv( vidx);
        time.Reset();
        Stokes.SetupPrMass(  &prM, lset/*, C.muF, C.muG*/);
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        time.Reset();
        PSchur_PCG_CL   schurSolver( prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data, 
            Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
      } break;
      
      case 2: // Nusselt solution
      {
        lset.Init( DistanceFct);
        Stokes.InitVel( &Stokes.v, Inflow);
      } break;
      
      case -1: // read from file
      {
        ReadEnsightP2SolCL reader( MG);
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
      } break;
      
      default:  
        lset.Init( DistanceFct);
        Stokes.InitVel( &Stokes.v, Null);
    }
    
    const double Vol= C.Filmdicke * C.mesh_size[0] * C.mesh_size[2];
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    EnsightP2SolOutCL ensight( MG, &ens_idx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo", 
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "falling film", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true); 
    ensight.DescribeScalar( "Pressure", datpr,  true); 
    ensight.DescribeVector( "Velocity", datvec, true); 

    ensight.putGeom( datgeo);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    Stokes.SetupPrMass(  &prM, lset/*, C.muF, C.muG*/);
    Stokes.SetupPrStiff( &prA);
    ISPreCL ispc( prA.Data, prM.Data, C.theta*C.dt*C.muF/C.rhoF);
   
    ISPSchur_PCG_CL ISPschurSolver( ispc,  C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    
    CouplLevelsetStokes2PhaseCL<StokesProblemT, ISPSchur_PCG_CL> 
        cpl( Stokes, lset, ISPschurSolver, C.theta);

    cpl.SetTimeStep( C.dt);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( C.FPsteps);
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
            lset.ReparamFastMarching( C.RepMethod, true);
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

    ensight.CaseEnd();
    std::cerr << std::endl;
}

} // end of namespace DROPS


void MarkFilm (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel= ~0)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        bool ref= false;
        int num_pos= 0;
        for (int i=0; i<4; ++i)
        {
            const double d= DistanceFct( It->GetVertex(i)->GetCoord());
            if (d<1e-4)
                ref= true;
            num_pos+= d>0;
        }
        if ( DistanceFct( GetBaryCenter(*It))<1e-4 )
            ref= true;
        if (num_pos!=4 && num_pos!=0)
            ref= true;
        if (ref)
            It->SetRegRefMark();
    }
}


void MarkHalf (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel= ~0)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( GetBaryCenter(*It)[1] < C.mesh_size[1]/2. )
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

    typedef ZeroFlowCL                              CoeffT;
    typedef DROPS::InstatStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::Point3DCL orig, e1, e2, e3;
    orig[2]= -C.mesh_size[2]/2;
    e1[0]= C.mesh_size[0];
    e2[1]= C.mesh_size[1];
    e3[2]= C.mesh_size[2];
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, int( C.mesh_res[0]), int( C.mesh_res[1]), int( C.mesh_res[2]) );
    
    if (C.BndCond.size()!=6)
    {
        std::cerr << "too many/few bnd conditions!\n"; return 1;
    }
    DROPS::BndCondT bc[6], bc_ls[6];
    DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6];
    
    for (int i=0; i<6; ++i)
    {
        bc_ls[i]= DROPS::Nat0BC;
        switch(C.BndCond[i])
        {
            case 'w': case 'W': 
                bc[i]= DROPS::WallBC;    bnd_fun[i]= &Null;   break;
            case 'i': case 'I': 
                bc[i]= DROPS::DirBC;     bnd_fun[i]= &Inflow; break;
            case 'o': case 'O': 
                bc[i]= DROPS::OutflowBC; bnd_fun[i]= &Null;   break;
            case '1': 
                bc_ls[i]= bc[i]= DROPS::Per1BC;    bnd_fun[i]= &Null;   break;
            case '2': 
                bc_ls[i]= bc[i]= DROPS::Per2BC;    bnd_fun[i]= &Null;   break;
            default: 
                std::cerr << "Unknown bnd condition \"" << C.BndCond[i] << "\"\n";
                return 1;
        }
    }

/*
    const DROPS::BndCondT bc[6]= 
//        { DROPS::WallBC, DROPS::DirBC, DROPS::OutflowBC, DROPS::OutflowBC, DROPS::DirBC, DROPS::OutflowBC};
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    //    foil, air_infty, side, side, top, bottom
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
//        { &Null, &Inflow, &Null, &Null, &Inflow, &Null}; 
        { &Null, &Null, &Null, &Null, &Null, &Null}; 
    */
        
    MyStokesCL prob(builder, CoeffT(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun, bc_ls));

    DROPS::MultiGridCL& mg = prob.GetMG();
    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    
    DROPS::LevelsetP2CL lset( mg, DROPS::LevelsetP2CL::BndDataT( 6, bc_ls), 
        prob.GetCoeff().SurfTens, C.lset_theta, C.lset_SD, 0, C.lset_iter, C.lset_tol, C.CurvDiff);

    for (DROPS::BndIdxT i=0, num= bnd.GetNumBndSeg(); i<num; ++i)
    {
        std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }
    
    for (int i=0; i<C.num_ref; ++i)
    {
//        MarkFilm( mg);
//        MarkHalf( mg);
        MarkAll( mg);
        mg.Refine();
    } 

    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    mg.SizeInfo( std::cerr);
    std::cerr << "Film Reynolds number Re_f = " 
              << C.rhoF*C.rhoF*C.g[0]*std::pow(C.Filmdicke,3)/C.muF/C.muF/3 << std::endl;
    std::cerr << "max. inflow velocity at film surface = "
              << C.rhoF*C.g[0]*C.Filmdicke*C.Filmdicke/C.muF/2 << std::endl;
    Strategy( prob, lset);  // do all the stuff
    
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
