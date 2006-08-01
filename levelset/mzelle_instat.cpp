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


DROPS::ParamMesszelleCL C;
enum StokesMethod { schur= 0,       inexactuzawa= 1,       minres= 2,
                    schurMG= 3,     inexactuzawaMG= 4,     minresMG= 5,
                    schurfullMG= 6, inexactuzawafullMG= 7, minresfullMG= 8};

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

class ISPSchur2_MG_CL: public PSchurSolver2CL<PCGSolverCL<MGPreCL>,
                                              PCGSolverCL<ISPreCL> >
{
  public:
    typedef PCGSolverCL<MGPreCL> innerSolverT;
    typedef PCGSolverCL<ISPreCL> outerSolverT;

  private:
    innerSolverT innerSolver_;
    outerSolverT outerSolver_;

  public:
    ISPSchur2_MG_CL(MGPreCL& Apc, ISPreCL& Spc,
        int outer_iter, double outer_tol,
        int inner_iter, double inner_tol)
        : PSchurSolver2CL<innerSolverT, outerSolverT>(
            innerSolver_, outerSolver_, outer_iter, outer_tol
          ),
          innerSolver_( Apc, inner_iter, inner_tol),
          outerSolver_( Spc, outer_iter, outer_tol)
        {}
};

class ISPSchur2_fullMG_CL: public PSchurSolver2CL<PCGSolverCL<MGPreCL>,
                                                  PCGSolverCL<ISMGPreCL> >
{
  public:
    typedef PCGSolverCL<MGPreCL> innerSolverT;
    typedef PCGSolverCL<ISMGPreCL> outerSolverT;

  private:
    innerSolverT innerSolver_;
    outerSolverT outerSolver_;

  public:
    ISPSchur2_fullMG_CL(MGPreCL& Apc, ISMGPreCL& Spc,
        int outer_iter, double outer_tol,
        int inner_iter, double inner_tol)
        : PSchurSolver2CL<innerSolverT, outerSolverT>(
            innerSolver_, outerSolver_, outer_iter, outer_tol
          ),
          innerSolver_( Apc, inner_iter, inner_tol),
          outerSolver_( Spc, outer_iter, outer_tol)
        {}
};


typedef DiagPreCL<MGPreCL,ISPreCL>  MyISMinresMGPreCL;
typedef DiagPreCL<MGPreCL,ISMGPreCL>  MyISMinresfullMGPreCL;

//=============================================================================
// PMinres solver for the instationary Stokes-Equations with MG for velocities
//=============================================================================
class MyPMinresSP_Diag_CL: public PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL,
    MyISMinresMGPreCL> >
{
  private:
    MyISMinresMGPreCL pre_;
    PLanczosONB_SPCL<MatrixCL, VectorCL, MyISMinresMGPreCL> q_;

  public:
    MyPMinresSP_Diag_CL(MGPreCL& Apc, ISPreCL& Spc, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, MyISMinresMGPreCL> >( q_, maxiter, tol),
         pre_( Apc, Spc), q_( pre_)
    {}
};

//=============================================================================
// PMinres solver for the instationary Stokes-Equations with full MG
//=============================================================================
class MyPMinresSP_fullMG_CL: public PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL,
    MyISMinresfullMGPreCL> >
{
  private:
    MyISMinresfullMGPreCL pre_;
    PLanczosONB_SPCL<MatrixCL, VectorCL, MyISMinresfullMGPreCL> q_;

  public:
    MyPMinresSP_fullMG_CL(MGPreCL& Apc, ISMGPreCL& Spc, int maxiter, double tol)
        :PMResSPCL<PLanczosONB_SPCL<MatrixCL, VectorCL, MyISMinresfullMGPreCL> >( q_, maxiter, tol),
         pre_( Apc, Spc), q_( pre_)
    {}
};

// copied from isdrops.cpp
// We know, there are only natural boundary conditions.
template<class Coeff>
void
SetupPrStiffMG(DROPS::InstatStokes2PhaseP2P1CL<Coeff>& stokes,
    DROPS::MGDataCL& MGData, const LevelsetP2CL& lset)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::IdxDescCL* c_idx= 0;
    for(DROPS::Uint lvl= 0; lvl<=mg.GetLastLevel(); ++lvl) {
        MGData.push_back( DROPS::MGLevelDataCL());
        DROPS::MGLevelDataCL& tmp= MGData.back();
        std::cerr << "Pressure-MG:            Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 1);
        stokes.CreateNumberingPr( lvl, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
        stokes.SetupPrStiff( &tmp.A, lset);
        if(lvl!=0) {
            std::cerr << "                        Create Prolongation on Level " << lvl << std::endl;
            DROPS::SetupP1ProlongationMatrix( mg, tmp.P, c_idx, &tmp.Idx);
        }
        c_idx= &tmp.Idx;
    }
    std::cerr << "Check MG-Data..." << std::endl;
    std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
    std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
    DROPS::CheckMGData( MGData.begin(), MGData.end());
}

template<class Coeff>
void
SetupPrMassMG(DROPS::InstatStokes2PhaseP2P1CL<Coeff>& stokes,
    DROPS::MGDataCL& MGData, const LevelsetP2CL& lset)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::IdxDescCL* c_idx= 0;
    for(DROPS::Uint lvl= 0; lvl<=mg.GetLastLevel(); ++lvl) {
        MGData.push_back( DROPS::MGLevelDataCL());
        DROPS::MGLevelDataCL& tmp= MGData.back();
        std::cerr << "Mass-Pressure-MG:       Create MGData on Level " << lvl << std::endl;
        tmp.Idx.Set( 1);
        stokes.CreateNumberingPr( lvl, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "                        Create StiffMatrix     " << (&tmp.Idx)->NumUnknowns <<std::endl;
        stokes.SetupPrMass( &tmp.A, lset);
        std::cerr << tmp.A.Data.num_nonzeros() << " nonzeros in M_pr.\n";
        if(lvl!=0) {
            std::cerr << "                        Create Prolongation on Level " << lvl << std::endl;
            DROPS::SetupP1ProlongationMatrix( mg, tmp.P, c_idx, &tmp.Idx);
        }
        c_idx= &tmp.Idx;
    }
    std::cerr << "Check MG-Data..." << std::endl;
    std::cerr << "                begin     " << MGData.begin()->Idx.NumUnknowns << std::endl;
    std::cerr << "                end       " << (--MGData.end())->Idx.NumUnknowns << std::endl;
    DROPS::CheckMGData( MGData.begin(), MGData.end());
}

template<class Coeff>
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    LevelsetP2CL lset( MG, Stokes.GetCoeff().SurfTens, C.lset_theta, C.lset_SD, C.RepDiff, C.lset_iter, C.lset_tol, C.CurvDiff);

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
    ISPreCL ispc( prA.Data, prM.Data, 1./C.dt, C.theta);
    typedef PCG_SsorCL SPcSolverT;
    SPcSolverT spcsolver( SSORPcCL( 1.0), 100, 0.02, /*relative*/ true);
    ISNonlinearPreCL<SPcSolverT> isnonlinpc( spcsolver, prA.Data, prM.Data, 1./C.dt, C.theta); // May be used for inexact Uzawa.
    MGDataCL prA_MG;
    SetupPrStiffMG( Stokes, prA_MG, lset);
    MGDataCL prM_MG;
    SetupPrMassMG( Stokes, prM_MG, lset);
    ISMGPreCL ispcMG( prA_MG, prM_MG, 1./C.dt, C.theta, 1);
    typedef PCG_SsorCL ASolverT;
    ASolverT Asolver( SSORPcCL( 1.0), 500, 0.02, /*relative*/true);
    typedef SolverAsPreCL<ASolverT> APcT;
    APcT velpc( Asolver);

    // Available Stokes-solver
    ISPSchur_PCG_CL ISPschurSolver( ispc, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
//    typedef InexactUzawaCL<APcT, ISPreCL, APC_SYM> InexactUzawa_CL;
//    InexactUzawa_CL inexactUzawaSolver( velpc, ispc, C.outer_iter, C.outer_tol);
    typedef InexactUzawaCL<APcT, ISNonlinearPreCL<SPcSolverT>, APC_SYM> InexactUzawaNonlinear_CL;
    InexactUzawaNonlinear_CL inexactUzawaSolver( velpc, isnonlinpc, C.outer_iter, C.outer_tol);
    PMinresSP_Diag_CL stokessolver( ispc, 1, C.outer_iter, C.outer_tol);

    // Available Stokes-solver: MG for velocities
    MGDataCL VelMGPreData;
    MGPreCL velmgpc( VelMGPreData);
    ISPSchur2_MG_CL ISPschur2SolverMG( velmgpc, ispc,
        C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    typedef InexactUzawaCL<MGPreCL, ISPreCL, APC_SYM_LINEAR> InexactUzawaMG_CL;
    InexactUzawaMG_CL inexactUzawaSolverMG( velmgpc, ispc, C.outer_iter, C.outer_tol);
    MyPMinresSP_Diag_CL stokessolverMG( velmgpc, ispc, C.outer_iter, C.outer_tol);

    // Available Stokes-solver: full MG
    ISPSchur2_fullMG_CL ISPschur2SolverfullMG( velmgpc, ispcMG,
        C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    typedef InexactUzawaCL<MGPreCL, ISMGPreCL, APC_SYM_LINEAR> InexactUzawaFullMG_CL;
    InexactUzawaFullMG_CL inexactUzawaSolverFullMG( velmgpc, ispcMG, C.outer_iter, C.outer_tol);
    MyPMinresSP_fullMG_CL stokessolverfullMG( velmgpc, ispcMG, C.outer_iter, C.outer_tol);

    PSchur_PCG_CL   schurSolver( prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

    switch (C.IniCond)
    {
      case 1: case 2: // stationary flow with/without drop
      {
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

        time.Reset();
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
            Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";

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

    // XXX: Wozu? Das passiert doch schon bei der Konstruktion.
    ISPschurSolver.SetTol( C.outer_tol);
    inexactUzawaSolver.SetTol( C.outer_tol);
    stokessolver.SetTol( C.outer_tol);

    CouplLevelsetStokes2PhaseCL<StokesProblemT, ISPSchur_PCG_CL>
        cpl1( Stokes, lset, ISPschurSolver, C.theta);
//    CouplLevelsetStokes2PhaseCL<StokesProblemT, InexactUzawa_CL>
//        cpl2( Stokes, lset, inexactUzawaSolver, C.theta);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, InexactUzawaNonlinear_CL>
        cpl2( Stokes, lset, inexactUzawaSolver, C.theta);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, PMinresSP_Diag_CL>
        cpl3( Stokes, lset, stokessolver, C.theta);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, ISPSchur2_MG_CL>
        cpl4( Stokes, lset, ISPschur2SolverMG, C.theta,
        C.StokesMethod==schurMG, &VelMGPreData);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, InexactUzawaMG_CL>
        cpl5( Stokes, lset, inexactUzawaSolverMG, C.theta,
        C.StokesMethod==inexactuzawaMG, &VelMGPreData);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, MyPMinresSP_Diag_CL>
        cpl6( Stokes, lset, stokessolverMG, C.theta,
        C.StokesMethod==minresMG, &VelMGPreData);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, ISPSchur2_fullMG_CL>
        cpl7( Stokes, lset, ISPschur2SolverfullMG, C.theta,
        C.StokesMethod==schurfullMG, &VelMGPreData);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, InexactUzawaFullMG_CL>
        cpl8( Stokes, lset, inexactUzawaSolverFullMG, C.theta,
        C.StokesMethod==inexactuzawafullMG, &VelMGPreData);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, MyPMinresSP_fullMG_CL>
        cpl9( Stokes, lset, stokessolverfullMG, C.theta,
        C.StokesMethod==minresfullMG, &VelMGPreData);

    switch (C.StokesMethod) {
      case schur:              cpl1.SetTimeStep( C.dt); break;
      case inexactuzawa:       cpl2.SetTimeStep( C.dt); break;
      case minres:             cpl3.SetTimeStep( C.dt); break;
      case schurMG:            cpl4.SetTimeStep( C.dt); break;
      case inexactuzawaMG:     cpl5.SetTimeStep( C.dt); break;
      case minresMG:           cpl6.SetTimeStep( C.dt); break;
      case schurfullMG:        cpl7.SetTimeStep( C.dt); break;
      case inexactuzawafullMG: cpl8.SetTimeStep( C.dt); break;
      case minresfullMG:       cpl9.SetTimeStep( C.dt); break;
      default: std::cerr << "Strategy: Please Choose a Stokes-solver." << std::endl;
    }
    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        switch (C.StokesMethod) {
          case schur:              cpl1.DoStep( C.cpl_iter); break;
          case inexactuzawa:       cpl2.DoStep( C.cpl_iter); break;
          case minres:             cpl3.DoStep( C.cpl_iter); break;
          case schurMG:            cpl4.DoStep( C.cpl_iter); break;
          case inexactuzawaMG:     cpl5.DoStep( C.cpl_iter); break;
          case minresMG:           cpl6.DoStep( C.cpl_iter); break;
          case schurfullMG:        cpl7.DoStep( C.cpl_iter); break;
          case inexactuzawafullMG: cpl8.DoStep( C.cpl_iter); break;
          case minresfullMG:       cpl9.DoStep( C.cpl_iter); break;
	}
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

    typedef ZeroFlowCL                              CoeffT;
    typedef DROPS::InstatStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

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

    Strategy( prob);  // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
