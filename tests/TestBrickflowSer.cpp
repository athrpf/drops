/// \file TestBrickflowPar.cpp
/// \brief Testing parallel solver for a simple 2-phase problem
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross, Joerg Peters, Trung Hieu Nguyen, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "num/solver.h"
#include "num/stokessolver.h"
#include "stokes/integrTime.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"

 // include problem class
#include "navstokes/instatnavstokes2phase.h"
#include "levelset/coupling.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#endif

DROPS::ParParamMesszelleNsCL C;

const char line[] ="------------------------------------------------------------";

/// \brief Display number of unknowns
void DisplayNumUnknowns(const DROPS::MultiGridCL&, const DROPS::VecDescCL& x)
/// accumulated unknowns: accumulation of the unknowns of all processors. Some unknowns
///   are counted multiple due to the overlapping of edges and vertices
/// global unknowns: all unknowns are counted just once
{
    const DROPS::Ulint glo_num_unk = x.Data.size();
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    std::cout << "  + Number of DOF of index "<<idx<<":  " <<glo_num_unk<< std::endl;
}

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

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = (p[0]/C.r_inlet)*(p[0]/C.r_inlet)-1,
                 z = (p[2]/C.r_inlet)*(p[2]/C.r_inlet)-1;

    ret[1]= x * z * C.Anstroem;// * (1-ampl*std::cos(2*M_PI*freq*t));  // Rohr
    return ret;
}

// droplet
double DistanceFct1( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= C.Mitte-p;
    return d.norm()-C.Radius;
}

// middle
double DistanceFct( const DROPS::Point3DCL& p)
{
    return p[1]-0.015;
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }

namespace DROPS{

class EnsightOutCL
{
  private:
    double             time_;
    double             timeInc_;
    int                counter_;
    EnsightP2SolOutCL* ensight_;
    MultiGridCL*       mg_;
    std::string        datgeo,
                       datscl,
                       datpr,
                       datvel;

  public:
    EnsightOutCL() : time_(0.), ensight_(0), mg_(0) {}
    ~EnsightOutCL() { ensight_->CaseEnd(); }

    void init(MultiGridCL& mg, const IdxDescCL* idx, int timesteps, double dt)
    {
        counter_=0;
        timeInc_= dt;
        mg_     = &mg;
        ensight_= new EnsightP2SolOutCL( *mg_, idx, /*binary=*/true);
        const string filename= C.EnsDir + "/" + C.EnsCase;
        datgeo= filename+".geo",
        datpr = filename+".pr" ,
        datvel= filename+".vel",
        datscl= filename+".scl";

        ensight_->CaseBegin( string(C.EnsCase+".case").c_str(), timesteps);
        ensight_->DescribeGeom( "Messzelle",  datgeo,  true);
        ensight_->DescribeScalar( "Levelset", datscl, true);
        ensight_->DescribeScalar( "Pressure", datpr,  true);
        ensight_->DescribeVector( "Velocity", datvel, true);
    }

    template<class Coeff>
    void write(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
    {
        std::cout <<"- Writing ensight file no "<<counter_<<" for time "<<time_<<std::endl;
        ensight_->putGeom(   datgeo, time_);
        ensight_->putScalar( datpr, Stokes.GetPrSolution(), time_);
        ensight_->putVector( datvel, Stokes.GetVelSolution(), time_);
        ensight_->putScalar( datscl, lset.GetSolution(), time_);
        ensight_->Commit();
        time_+=timeInc_;
        ++counter_;
    }

} ensight;

template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset,
                          std::ofstream *infofile)
{
    TimerCL time;
    double duration;
    // droplet-information
    double maxGradPhi, Volume;
    Point3DCL bary_drop, min_drop, max_drop;

    // Initial velocity is zero
    Stokes.InitVel( &Stokes.v, Null);
    lset.Init( DistanceFct1);

    lset.GetInfo( maxGradPhi, Volume, bary_drop, min_drop, max_drop);
    (*infofile) << Stokes.t << '\t' << maxGradPhi << '\t' << Volume << '\t' << bary_drop << '\t' << min_drop << '\t' << max_drop << std::endl;

    switch (C.IniCond)
    {
      // stationary flow with/without drop
      case 1: case 2:
      {
        // Setting up solver
        typedef DummyPcCL SPcT;
        SPcT ispc;
        typedef JACPcCL  APcPcT;
        APcPcT Apcpc;
        typedef PCGSolverCL<APcPcT> ASolverT;
        ASolverT Asolver( Apcpc, 500, 0.02, /*relative=*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver/*, &std::cout*/);
        typedef InexactUzawaCL<APcT, SPcT, APC_SYM> OseenSolverT;
        OseenSolverT schurSolver( Apc, ispc, C.outer_iter, C.outer_tol, C.stokes_inner_red, 500);

        VelVecDescCL curv(&Stokes.vel_idx),
                    cplN(&Stokes.vel_idx);

        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        curv.Data=0.;
        lset.AccumulateBndIntegral( curv);
        time.Stop(); duration=time.GetTime();
        std::cout << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";

        //Solve initial problem
        double theta= C.stat_theta, nl= C.stat_nonlinear;
        time.Reset();
        int step=0;
        int iters=0;
        do
        {
            Stokes.SetupNonlinear( &Stokes.N, &Stokes.v, &cplN, lset, Stokes.t);
            MatrixCL mat;
            mat.LinComb( 1, Stokes.A.Data, nl*theta, Stokes.N.Data);
            cplN.Data-= (1-theta) * (Stokes.N.Data * Stokes.v.Data);
            schurSolver.Solve( mat, Stokes.B.Data,
                                    Stokes.v.Data, Stokes.p.Data,
                                    VectorCL( Stokes.b.Data + nl*cplN.Data), Stokes.c.Data);
            std::cout << "- Solving lin. Stokes ("<<step<<"): iter "<<schurSolver.GetIter()
                      <<", resid "<<schurSolver.GetResid()<<std::endl;
            ++step; iters+= schurSolver.GetIter();
        } while (schurSolver.GetIter() > 0);
        time.Stop(); duration=time.GetTime();
        std::cout << "- Solving Stokes for initialization took "<<duration<<" sec, "
                  << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<schurSolver.GetResid()<<'\n';
      }break;
      case 3:
        {
            ReadEnsightP2SolCL reader( Stokes.GetMG(), false);
            reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
            reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
            reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
            std::cout << "- Initial Conditions successfull read\n";
        } break;
    }
}

template<typename Coeff>
  void SolveCoupledNS(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset,
                      AdapTriangCL& adap, std::ofstream* infofile)
{
    TimerCL time;
    double duration;

    // droplet information
    double maxGradPhi, Volume;
    Point3DCL bary_drop, min_drop, max_drop;

    // type of the problem
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    const double Vol= 4./3.*M_PI*std::pow(C.Radius,3);
    double relVol = lset.GetVolume()/Vol;

    // linear solvers
    typedef DummyPcCL  SPcT;
    SPcT ispc;
    typedef JACPcCL  APcPcT;
    APcPcT Apcpc;
    typedef GMResSolverCL<APcPcT>    ASolverT;        // GMRes-based APcT
    ASolverT Asolver( Apcpc, C.navstokes_pc_restart, C.navstokes_pc_iter, C.navstokes_pc_reltol,
                       /*relative=*/ true, RightPreconditioning);

    typedef SolverAsPreCL<ASolverT> APcT;
    APcT Apc( Asolver/*,&std::cout*/);

    // stokes solver
    typedef InexactUzawaCL<APcT, SPcT> OseenSolverT;
    OseenSolverT oseensolver( Apc, ispc, C.navstokes_outer_iter, C.navstokes_outer_tol, C.navstokes_inner_red);

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT, OseenSolverT> NSSolverT;
    NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red, &std::cout);

    // coupling levelset NavStokes
    time.Reset();
    // Coupling Navier-Stokes with Levelset
    typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
    CouplingT cpl( Stokes, lset, nssolver, C.nonlinear);

    time.Stop();
    duration=time.GetTime();
    std::cout << "- Updating discretization took "<<duration<<" sec.\n";

    // Set time step and create matrices
    cpl.SetTimeStep( C.dt);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        TimerCL step_time;
        step_time.Reset();
        std::cout << "=================================================================================== Schritt " << step << ":\n"
                  << " Idx for vel  "<<Stokes.v.RowIdx->GetIdx()
                  << "\n Idx for pr   "<<Stokes.p.RowIdx->GetIdx()
                  << "\n Idx for lset "<<lset.Phi.RowIdx->GetIdx()<<std::endl;
        time.Reset();

        if (C.ref_freq && step%C.ref_freq==0)
        {
            std::cout << "==> Adaptive Refinement of MultiGrid"<<std::endl;

            adap.UpdateTriang( Stokes, lset);
            if (adap.WasModified() )
            {
                cpl.Update();
                // don't forget to update the pr mass/stiff matrix for the schur compl. preconditioner!!
                Stokes.prM.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);
                Stokes.SetupPrMass( &Stokes.prM, lset);
                Stokes.prA.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);
                Stokes.SetupPrStiff( &Stokes.prA, lset);
            }
        }

        DisplayNumUnknowns(adap.GetMG(), Stokes.v);
        DisplayNumUnknowns(adap.GetMG(), Stokes.p);
        DisplayNumUnknowns(adap.GetMG(), lset.Phi);

        std::cout << "==> Solving coupled Levelset-Navier-Stokes problem ....\n";


        cpl.DoStep( C.cpl_iter);
        if (C.ensight && step%C.ensight==0)
            ensight.write(Stokes, lset);

        // Write out droplet information
        lset.GetInfo( maxGradPhi, Volume, bary_drop, min_drop, max_drop);
        (*infofile) << Stokes.t << '\t' << maxGradPhi << '\t' << Volume << '\t' << bary_drop << '\t' << min_drop << '\t' << max_drop << std::endl;

        time.Stop(); duration=time.GetTime();
        std::cout << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec."<<std::endl;
        relVol = lset.GetVolume()/Vol;
        std::cout << "- rel. Volume: " << relVol << std::endl;
        if (C.VolCorr)
        {
            std::cout << "\n==> Adjust volume ...\n";
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            lset.Phi.Data+= dphi;
            relVol = lset.GetVolume()/Vol;
            std::cout << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< std::endl;
        }

        // Reparametrization of levelset function
        if (C.RepFreq && step%C.RepFreq==0)
        {
            std::cout << "\n==> Reparametrization with FastMarching algorithm"<<std::endl;
            time.Reset();
            lset.ReparamFastMarching();
            time.Stop(); duration=time.GetTime();
            relVol = lset.GetVolume()/Vol;
            std::cout << "- FastMarching took "<<duration<<" sec."<<std::endl;
            std::cout << "- rel. Volume: " << relVol << std::endl;
            if (C.VolCorr)
            {
                std::cout << "\n==> Adjust volume ...\n";

                double dphi= lset.AdjustVolume( Vol, 1e-9);
                lset.Phi.Data+= dphi;
                relVol = lset.GetVolume()/Vol;
                std::cout << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< std::endl;
            }
        }

        step_time.Stop();
        duration=step_time.GetTime();
        std::cout <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
    }

    infofile->close();
    delete infofile;
}

template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, MultiGridCL& mg)
{
    sigma= Stokes.GetCoeff().SurfTens;
    TimerCL time;
    SurfaceTensionCL sf( sigmaf, 0);
    LevelsetP2CL lset( mg, sf, C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff);
    lset.SetSurfaceForce( SF_ImprovedLB);

    AdapTriangCL adapt( mg, C.ref_width, 0, C.ref_flevel);
    adapt.MakeInitialTriang(::DistanceFct1);

    std::ofstream *infofile=0;
    infofile = new std::ofstream( string(C.EnsCase + ".info").c_str());
    IdxDescCL *vidx= &Stokes.vel_idx,
              *pidx= &Stokes.pr_idx,
              *lidx= &lset.idx;
    VecDescCL *v   = &Stokes.v,
              *p   = &Stokes.p,
              *l   = &lset.Phi;

    Stokes.CreateNumberingVel( mg.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), pidx);
    lset.CreateNumbering     ( mg.GetLastLevel(), lidx);

    // Init ensight
    ensight.init(mg, lidx, 2*C.num_steps+1, C.dt);

    // Tell matrices and vectors about the numbering
    v->SetIdx( vidx);               p->SetIdx( pidx);
    l->SetIdx( lidx);
    Stokes.b.SetIdx( vidx);         Stokes.c.SetIdx( pidx);
    Stokes.A.SetIdx(vidx, vidx);    Stokes.B.SetIdx(pidx, vidx);
    Stokes.M.SetIdx(vidx, vidx);    Stokes.N.SetIdx(vidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx); Stokes.prA.SetIdx( pidx, pidx);


    //Setup initial problem
    std::cout << "=================================================================================== Init:\n"
              << "==> Initialize Problem\n";
    InitProblemWithDrop(Stokes, lset, infofile);

    ensight.write(Stokes, lset);

    SolveCoupledNS(Stokes, lset, adapt, infofile);
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=2)
    {
        std::cout << "You have to specify one parameter:\n\t"
                  << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cout << C << std::endl;

    DROPS::TimerCL alltime;

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
        std::cout << "error while reading geometry information: " << mesh << "\n";
        return 1;
    }

    C.r_inlet= dx/2;
    DROPS::Point3DCL orig, px, py, pz;
    px[0]= dx; py[1]= dy; pz[2]= dz;


    DROPS::TimerCL time;
    DROPS::MGBuilderCL * mgb;
    mgb = new DROPS::BrickBuilderCL( orig, px, py, pz, nx, ny, nz);

    const bool bc[6]=
      {false, false, true, false, false, false};    // Rohr
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
      { &Null, &Null, &Null, &Inflow, &Null, &Null};

    DROPS::MultiGridCL mg(*mgb);

    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();


    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun));

    Strategy( prob, mg);    // do all the stuff

    alltime.Stop();
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

