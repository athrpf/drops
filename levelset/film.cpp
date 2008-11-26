//**************************************************************************
// File:    film.cpp                                                       *
// Content: flow of falling film                                           *
// Author:  Sven Gross, Joerg Grande, Patrick Esser, IGPM RWTH Aachen      *
//**************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolverfactory.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/adaptriang.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/mgobserve.h"
#include <fstream>


DROPS::ParamFilmCL C;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0


class ZeroFlowCL
{
// \Omega_1 = Film,    \Omega_2 = Gasphase
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamFilmCL& C)
      : rho( DROPS::JumpCL( C.rhoF, C.rhoG ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muF,  C.muG),   DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Film,    \Omega_2 = Gasphase
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( const DROPS::ParamFilmCL& C)
      : rho( DROPS::JumpCL( 1., C.rhoG/C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu ( DROPS::JumpCL( 1., C.muG/C.muF),    DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma/C.rhoF), g( C.g)    {}
};


DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret(0.);
    const double d= p[1]/C.Filmdicke;
    static const double u= C.rhoF*C.g[0]*C.Filmdicke*C.Filmdicke/C.muF/2;
    ret[0]= d<=1 ? (2*d-d*d)*u * (1 + C.PumpAmpl*std::sin(2*M_PI*t*C.PumpFreq))
                 : (C.mesh_size[1]-p[1])/(C.mesh_size[1]-C.Filmdicke)*u;
    return ret;
}

double DistanceFct( const DROPS::Point3DCL& p)
{
    // wave length = 100 x film width
    const double wave= std::sin(2*M_PI*p[0]/C.mesh_size[0]),
        z= p[2]/C.mesh_size[2]*2; // z \in [-1,1]
//    return p[1] - C.Filmdicke * (1 + C.PumpAmpl*wave);
//    return p[1] - C.Filmdicke * (1 + C.PumpAmpl*(wave + C.AmplZ*std::cos(z*M_PI)));
    const double z_fac=  (1 + C.AmplZ/2*std::cos(z*M_PI));  // (z=+-1) 1-C.AmplZ <= z_fac <= 1+C.AmplZ (z=0)
    return p[1] - C.Filmdicke * (1 + C.PumpAmpl*wave) * z_fac;
}

bool periodic_xz( const DROPS::Point3DCL& p, const DROPS::Point3DCL& q)
{ // matching y-z- or x-y-coords, resp.
    const DROPS::Point3DCL d= fabs(p-q),
                           L= fabs(C.mesh_size);
    return (d[1] + d[2] < 1e-12 && std::abs( d[0] - L[0]) < 1e-12)  // dy=dz=0 and dx=Lx
      ||   (d[0] + d[1] < 1e-12 && std::abs( d[2] - L[2]) < 1e-12)  // dx=dy=0 and dz=Lz
      ||   (d[1] < 1e-12 && std::abs( d[0] - L[0]) < 1e-12 && std::abs( d[2] - L[2]) < 1e-12);  // dy=0 and dx=Lx and dz=Lz
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }

namespace DROPS // for Strategy
{

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Ensight index.
///
/// The actual work is done in post_refine().
class EnsightIdxRepairCL: public MGObserverCL
{
  private:
	MultiGridCL& mg_;
	IdxDescCL&   idx_;

  public:
    EnsightIdxRepairCL( MultiGridCL& mg, IdxDescCL& idx)
      : mg_(mg), idx_(idx) {}

    void pre_refine  () {}
    void post_refine () { idx_.DeleteNumbering( mg_); idx_.CreateNumbering( mg_.GetLastLevel(), mg_); }

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
};


template<class StokesProblemT>
void Strategy( StokesProblemT& Stokes, LevelsetP2CL& lset, AdapTriangCL& adap)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    IdxDescCL ens_idx( P2_FE, NoBndDataCL<>());

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset);
    adap.push_back( &prrepair);
    EnsightIdxRepairCL ensrepair( MG, ens_idx);
    adap.push_back( &ensrepair);

    lset.CreateNumbering(      MG.GetLastLevel(), lidx, periodic_xz);
    lset.Phi.SetIdx( lidx);
    lset.Init( DistanceFct);
    if ( StokesSolverFactoryHelperCL<ParamFilmCL>().VelMGUsed(C))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL<ParamFilmCL>().PrMGUsed(C))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx, periodic_xz);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, periodic_xz, &lset);
    ens_idx.CreateNumbering( MG.GetLastLevel(), MG);

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
    switch (C.IniCond)
    {
      case 1: // stationary flow
      {
        TimerCL time;
        VelVecDescCL curv( vidx);
        time.Reset();
        Stokes.SetupPrMass(  &Stokes.prM, lset/*, C.muF, C.muG*/);
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        time.Reset();
        PCG_SsorCL PCGsolver( SSORPcCL(), C.inner_iter, C.inner_tol);
        PSchurSolverCL<PCG_SsorCL> schurSolver( PCGsolver, Stokes.prM.Data, C.outer_iter, C.outer_tol);

        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
            Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
      } break;

      case 2: // Nusselt solution
      {
        Stokes.InitVel( &Stokes.v, Inflow);
      } break;

      case -1: // read from file
      {
        ReadEnsightP2SolCL reader( MG);
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
        Stokes.UpdateXNumbering( pidx, lset, /*NumberingChanged*/ false);
        Stokes.p.SetIdx( pidx); // Zero-vector for now.
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr); // reads the P1-part of the pressure
      } break;

      default:
        Stokes.InitVel( &Stokes.v, ZeroVel);
    }

    const double Vol= lset.GetVolume(); // approx. C.Filmdicke * C.mesh_size[0] * C.mesh_size[2];
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    // Initialize Ensight6 output
    std::string ensf( C.EnsDir + "/" + C.EnsCase);
    Ensight6OutCL ensight( C.EnsCase + ".case", C.num_steps + 1);
    ensight.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(),   "falling film", ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),      "Levelset",     ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar( Stokes.GetPrSolution(),  "Pressure",     ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector( Stokes.GetVelSolution(), "Velocity",     ensf + ".vel", true));

    ensight.Write();

    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    // Stokes-Solver
    StokesSolverFactoryCL<StokesProblemT, ParamFilmCL> stokessolverfactory(Stokes, C);
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver();

    // Navier-Stokes-Solver
    typedef NSSolverBaseCL<StokesProblemT> SolverT;
    SolverT * navstokessolver = 0;
    if (C.nonlinear==0.0)
        navstokessolver = new NSSolverBaseCL<StokesProblemT>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<StokesProblemT>(Stokes, *stokessolver, C.ns_iter, C.ns_tol, C.ns_red);

    LinThetaScheme2PhaseCL<StokesProblemT, SolverT>
        cpl( Stokes, lset, *navstokessolver, C.theta, /*nonlinear*/ 0., /*implicitCurv*/ true);

    cpl.SetTimeStep( C.dt);
    stokessolverfactory.SetMatrixA( &navstokessolver->GetAN()->GetFinest());

    //for Stokes-MGM
    stokessolverfactory.SetMatrices( &navstokessolver->GetAN()->GetCoarsest(), &Stokes.B.Data.GetCoarsest(),
                                     &Stokes.M.Data.GetCoarsest(), &Stokes.prM.Data.GetCoarsest());

    UpdateProlongationCL PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);

    bool secondSerial= false;
    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( C.cpl_iter);
        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

        bool forceVolCorr= false, forceUpdate= false,
             doReparam= C.RepFreq && step%C.RepFreq == 0,
             doGridMod= C.ref_freq && step%C.ref_freq == 0;

        // volume correction before reparam/grid modification
        if (C.VolCorr && (doReparam || doGridMod)) {
                double dphi= lset.AdjustVolume( Vol, 1e-9, C.mesh_size[0] * C.mesh_size[2]);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                forceUpdate= true; // volume correction modifies the level set
        }

        // reparam levelset function
        if (doReparam) {
        	double lsetmaxGradPhi, lsetminGradPhi;
            lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
            std::cerr << "checking level set func: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
            if (lsetmaxGradPhi > 10 || lsetminGradPhi < 0.1) {
                lset.ReparamFastMarching( C.RepMethod, /*periodic*/ true);
                lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
                std::cerr << "after reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
                forceVolCorr= forceUpdate= true; // volume correction and update after reparam
            }
            else
            	std::cerr << "Gradient does not exceed bounds, reparametrization skipped.\n";
        }

        // grid modification
        if (doGridMod) {
            adap.UpdateTriang( lset);
            forceUpdate  |= adap.WasModified();
            forceVolCorr |= adap.WasModified();
            if (C.serialization_file != "none") {
                std::stringstream filename;
                filename << C.serialization_file;
                if (secondSerial) filename << "0";
                secondSerial = !secondSerial;
                MGSerializationCL ser( MG, filename.str().c_str());
                ser.WriteMG();
                filename << ".time";
                std::ofstream serTime( filename.str().c_str());
                serTime << "Serialization info:\ntime step = " << step << "\t\tt = " << step*C.dt << "\n";
                serTime.close();
            }
        }

        // volume correction
        if (C.VolCorr && (step%C.VolCorr==0 || forceVolCorr)) {
            double dphi= lset.AdjustVolume( Vol, 1e-9, C.mesh_size[0] * C.mesh_size[2]);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            forceUpdate= true;
        }

        // update
        if (forceUpdate)
            cpl.Update();

if (step%10==0)
        ensight.Write( step*C.dt);
    }

    std::cerr << std::endl;
    delete stokessolver;
    delete navstokessolver;
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


void MarkLower (DROPS::MultiGridCL& mg, double y_max, DROPS::Uint maxLevel= ~0)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( GetBaryCenter(*It)[1] < y_max )
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

    DROPS::Point3DCL orig, e1, e2, e3;
    orig[2]= -C.mesh_size[2]/2;
    e1[0]= C.mesh_size[0];
    e2[1]= C.mesh_size[1];
    e3[2]= C.mesh_size[2];
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, int( C.mesh_res[0]), int( C.mesh_res[1]), int( C.mesh_res[2]) );
    DROPS::MultiGridCL* mgp;
    if (C.deserialization_file == "none")
        mgp= new DROPS::MultiGridCL( builder);
    else {
        DROPS::FileBuilderCL filebuilder( C.deserialization_file, &builder);
        mgp= new DROPS::MultiGridCL( filebuilder);
    }

    if (C.BndCond.size()!=6)
    {
        std::cerr << "too many/few bnd conditions!\n"; return 1;
    }
    DROPS::BndCondT bc[6], bc_ls[6];
    DROPS::BoundaryCL::BndTypeCont bndType;
    DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6];

    for (int i=0; i<6; ++i)
    {
        bc_ls[i]= DROPS::Nat0BC;
        switch(C.BndCond[i])
        {
            case 'w': case 'W':
                bc[i]= DROPS::WallBC;    bnd_fun[i]= &DROPS::ZeroVel; bndType.push_back( DROPS::BoundaryCL::OtherBnd); break;
            case 'i': case 'I':
                bc[i]= DROPS::DirBC;     bnd_fun[i]= &Inflow;         bndType.push_back( DROPS::BoundaryCL::OtherBnd); break;
            case 'o': case 'O':
                bc[i]= DROPS::OutflowBC; bnd_fun[i]= &DROPS::ZeroVel; bndType.push_back( DROPS::BoundaryCL::OtherBnd); break;
            case '1':
                bc_ls[i]= bc[i]= DROPS::Per1BC;    bnd_fun[i]= &DROPS::ZeroVel; bndType.push_back( DROPS::BoundaryCL::Per1Bnd); break;
            case '2':
                bc_ls[i]= bc[i]= DROPS::Per2BC;    bnd_fun[i]= &DROPS::ZeroVel; bndType.push_back( DROPS::BoundaryCL::Per2Bnd); break;
            default:
                std::cerr << "Unknown bnd condition \"" << C.BndCond[i] << "\"\n";
                return 1;
        }
    }

    MyStokesCL prob( *mgp, CoeffT(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun, bc_ls), DROPS::P1X_FE, C.XFEMStab);

    const DROPS::BoundaryCL& bnd= mgp->GetBnd();
    bnd.SetPeriodicBnd( bndType, periodic_xz);

    sigma= prob.GetCoeff().SurfTens;
    DROPS::LevelsetP2CL lset( *mgp, DROPS::LevelsetP2CL::BndDataT( 6, bc_ls),
        &sigmaf, /*grad sigma*/ 0, C.lset_theta, C.lset_SD, 0, C.lset_iter, C.lset_tol, C.CurvDiff);

    for (DROPS::BndIdxT i=0, num= bnd.GetNumBndSeg(); i<num; ++i)
    {
        std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }

    DROPS::AdapTriangCL adap( *mgp, C.ref_width, 0, C.ref_flevel);
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (C.deserialization_file == "none")
        adap.MakeInitialTriang( DistanceFct);

    std::cerr << DROPS::SanityMGOutCL(*mgp) << std::endl;
    mgp->SizeInfo( std::cerr);
    std::cerr << "Film Reynolds number Re_f = "
              << C.rhoF*C.rhoF*C.g[0]*std::pow(C.Filmdicke,3)/C.muF/C.muF/3 << std::endl;
    std::cerr << "max. inflow velocity at film surface = "
              << C.rhoF*C.g[0]*C.Filmdicke*C.Filmdicke/C.muF/2 << std::endl;
    Strategy( prob, lset, adap);  // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
