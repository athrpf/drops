//**************************************************************************
// File:    TestSedPar.cpp                                                 *
// Content: Test case for simple sedimentation                             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestSedPar.cpp
/// \brief Test case for simple sedimentation

// include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include "parallel/exchange.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "num/parsolver.h"
#include "num/parprecond.h"
#include "num/stokessolver.h"
#include "num/parstokessolver.h"
#include "stokes/integrTime.h"
#include "levelset/adaptriang.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"

 // include problem class
#include "navstokes/instatnavstokes2phase.h"
#include "levelset/coupling.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

//#define SAMPLE
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#  ifdef SAMPLE
#    include <collectorAPI.h>
#  endif
#endif


enum TimePart{
    T_create_pmg,
    T_disc_init,
    T_create_ex,
    T_solve_NS_init
};

DROPS::TimeStoreCL Times(4);
DROPS::Point3DCL   dim_of_brick;

void SetDescriber()
{
    Times.SetDescriber(T_create_pmg,    "Initialization of parallel multigrid");
    Times.SetDescriber(T_disc_init,     "Discretizing Stokes/Curv for initial velocities");
    Times.SetDescriber(T_create_ex,     "Create ExchangeCL");
    Times.SetDescriber(T_solve_NS_init, "Solving NavStokes for initial velocities");
    Times.SetCounterDescriber("Moved MultiNodes");
}

DROPS::ParParamMesszelleNsCL C;

const char line[] ="------------------------------------------------------------";

/****************************************************************************
* P R I N T  M U L T I G R I D  I N  A S C I I  F I L E                     *
****************************************************************************/


bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream check(dat);
    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);
    check.close();
    return DROPS::Check(pmg_sane && mg_sane);
}

/// \brief Display number of unknowns
void DisplayNumUnknowns(const DROPS::MultiGridCL& MG, const DROPS::VecDescCL& x)
/// accumulated unknowns: accumulation of the unknowns of all processors. Some unknowns
///   are counted multiple due to the overlapping of edges and vertices
/// global unknowns: all unknowns are counted just once
{
    const DROPS::Ulint acc_num_unk = DROPS::GlobalSum(x.Data.size()),
                       glo_num_unk = x.RowIdx->GetGlobalNumUnknowns(MG);
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    if (DROPS::ProcCL::IamMaster())
        std::cout << "  + Number of DOF of index "<<idx<<" (accumulated/global):  "
                  <<acc_num_unk<< "/" <<glo_num_unk<< std::endl;
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

DROPS::SVectorCL<3> One( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(1.); }


DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = (p[0])*((p[0])-dim_of_brick[0]),
                 z = (p[2])*((p[2])-dim_of_brick[2]);

    ret[1]= x * z + 1.;
    ret[0]=1.; ret[2]=1.;
    return ret;
}

// droplet
double DistanceFct1( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.Mitte-p;
    const double avgRad = cbrt(C.Radius[0]*C.Radius[1]* C.Radius[2]);
    d/= C.Radius/avgRad;
    return d.norm()-avgRad;
}

double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }

namespace DROPS{

double offset=0;
double MovingDroplet(const Point3DCL& p)
{
    Point3DCL Mitte=C.Mitte; Mitte[1]+=offset;
    Point3DCL d=Mitte-p;
    const double avgRad = cbrt(C.Radius[0]*C.Radius[1]* C.Radius[2]);
    d/= C.Radius/avgRad;
    return d.norm()-avgRad;
}

double Pressure1(const Point3DCL& p, double=0.)
{
    return p[1]/0.03;
}

double Pressure2(const Point3DCL& p, double=0.)
{
    return Pressure1( p) - 100;
}

Point3DCL Velocity(const Point3DCL& p, double=0.)
{
    return Inflow(p,0);
}

void InitPr(VecDescCL& pr, const MultiGridCL& mg, double (*f)(const Point3DCL&, double))
{
    VectorCL& lsg= pr.Data;
    Uint lvl     = pr.GetLevel(),
         idx     = pr.RowIdx->GetIdx();

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(mg).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(mg).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
        if (sit->Unknowns.Exist(idx))
            lsg[sit->Unknowns(idx)]=f(sit->GetCoord(),0);
}

/// \brief Check interolated values and give (more or less usefull) information about the interpolation
template<class Coeff>
  void CheckInterPol(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
{
    std::valarray<double> dist_vals(-1., 10);
    std::valarray<double> global_dist_vals(-1., 10);
    double &vel_dist_vert = dist_vals[0], &abs_vel_dist_vert = dist_vals[1],
           &vel_dist_edge = dist_vals[2], &abs_vel_dist_edge = dist_vals[3],
           &pr_dist       = dist_vals[4], &abs_pr_dist       = dist_vals[5],
           &lset_dist_vert= dist_vals[6], &abs_lset_dist_vert= dist_vals[7],
           &lset_dist_edge= dist_vals[8], &abs_lset_dist_edge= dist_vals[9];
    double orig_pr_vert=-1., comp_pr_vert=-2.,
           orig_lset_vert=-1., comp_lset_vert=-2., orig_lset_edge=-1., comp_lset_edge=-2.,
           dist, min;
    Point3DCL orig_vel_vert, orig_vel_edge, comp_vel_vert, comp_vel_edge;

    char dat[30];
    std::sprintf(dat,"output/%i_difference.diff",DROPS::ProcCL::MyRank());
    static std::ofstream check(dat);
    static int time=0;

    typename InstatNavierStokes2PhaseP2P1CL<Coeff>::const_DiscPrSolCL  funP=Stokes.GetPrSolution();
    typename InstatNavierStokes2PhaseP2P1CL<Coeff>::const_DiscVelSolCL funV=Stokes.GetVelSolution();
    LevelsetP2CL::const_DiscSolCL                                      funL=lset.GetSolution();

    const Uint pr_idx  = funP.GetSolution()->RowIdx->GetIdx();
    const Uint vel_idx = funV.GetSolution()->RowIdx->GetIdx();
    const Uint lset_idx= funL.GetSolution()->RowIdx->GetIdx();

    const VertexCL *lset_vert=0, *pr_vert=0, *vel_vert=0;
    const EdgeCL   *lset_edge=0, *vel_edge=0;

    const MultiGridCL& mg=funP.GetMG();

    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()); sit!=mg.GetTriangVertexEnd(); ++sit){
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(pr_idx)){
            min  = std::min(funP.val(*sit), Pressure(sit->GetCoord()));
            dist =    std::fabs(funP.val(*sit)-Pressure(sit->GetCoord()))
                    / (min==0 ? 1. : min);
            if (dist>pr_dist){
                orig_pr_vert= Pressure(sit->GetCoord());
                comp_pr_vert= funP.val(*sit);
                pr_dist     = dist;
                abs_pr_dist = std::fabs(comp_pr_vert-Pressure(sit->GetCoord()));
                pr_vert     = &*sit;
            }
        }
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(vel_idx)){
            min  = std::min(funV.val(*sit).norm(),Velocity(sit->GetCoord()).norm());
            dist =    (funV.val(*sit)-Velocity(sit->GetCoord())).norm()
                    / (min==0 ? 1. : min);
            if (dist>abs_vel_dist_vert){
                orig_vel_vert    = Velocity(sit->GetCoord());
                comp_vel_vert    = funV.val(*sit);
                vel_dist_vert    = dist;
                abs_vel_dist_vert= (comp_vel_vert-orig_vel_vert).norm();
                vel_vert         = &*sit;
            }
        }
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(lset_idx)){
            min  =    std::min(funL.val(*sit),MovingDroplet(sit->GetCoord()));
            dist =    std::fabs(funL.val(*sit)-MovingDroplet(sit->GetCoord()))
                    / (min==0 ? 1. : min);
            if (dist>lset_dist_vert) {
                orig_lset_vert    = MovingDroplet(sit->GetCoord());
                comp_lset_vert    = funL.val(*sit);
                lset_dist_vert    = dist;
                abs_lset_dist_vert= std::fabs(comp_lset_vert-orig_lset_vert);
                lset_vert         = &*sit;
            }
        }
    }


    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()); sit!=mg.GetTriangEdgeEnd(); ++sit){

        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(vel_idx)){
            min  =   std::min(funV.val(*sit).norm(),Velocity(GetBaryCenter(*sit)).norm());
            dist =    (funV.val(*sit)-Velocity(GetBaryCenter(*sit))).norm()
                    / (min==0 ? 1. : min);
            if (dist>abs_vel_dist_edge){
                orig_vel_edge    = Velocity(GetBaryCenter(*sit));
                comp_vel_edge    = funV.val(*sit);
                vel_dist_edge    = dist;
                abs_vel_dist_edge= (comp_vel_edge-orig_vel_edge).norm();
                vel_edge         = &*sit;
            }
        }
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(lset_idx)){
            min  =    std::min(funL.val(*sit),MovingDroplet(GetBaryCenter(*sit)));
            dist =    std::fabs(funL.val(*sit)-MovingDroplet(GetBaryCenter(*sit)))
                    / (min==0 ? 1. : min);
            if (dist>lset_dist_edge) {
                orig_lset_edge    = MovingDroplet(GetBaryCenter(*sit));
                comp_lset_edge    = funL.val(*sit);
                lset_dist_edge    = dist;
                abs_lset_dist_edge= std::fabs(comp_lset_edge-orig_lset_edge);
                lset_edge         = &*sit;
            }
        }
    }

     GlobalMax(Addr(dist_vals), Addr(global_dist_vals), dist_vals.size(), ProcCL::Master());

     if (ProcCL::IamMaster())
         std::cout  << "Diferences (at time "<<time<<"):"
                    << "\npr            : "<< std::setw(10) << pr_dist        << ", absolut: "<<abs_pr_dist
                    << "\nvel (Vertex)  : "<< std::setw(10) << vel_dist_vert  << ", absolut: "<<abs_vel_dist_vert
                    << "\nvel (Edge)    : "<< std::setw(10) << vel_dist_edge  << ", absolut: "<<abs_vel_dist_edge
                    << "\nlset (Vertex) : "<< std::setw(10) << lset_dist_vert << ", absolut: "<<abs_lset_dist_vert
                    << "\nlset (Edge)   : "<< std::setw(10) << lset_dist_edge << ", absolut: "<<abs_lset_dist_edge
                    << std::endl<<std::endl;

     check  << "["<<ProcCL::MyRank()<<"] Diferences (at time "<<time++<<"):"
            << "\npr            : "<< std::setw(10) << pr_dist << " (" << std::setw(10) << std::fabs(comp_pr_vert-orig_pr_vert) << ')'
            << " on vertex "<< std::setw(8) <<pr_vert->GetGID()<<", "<<pr_vert->GetCoord()
            << ": "<<comp_pr_vert<<" instead of "<<orig_pr_vert

            << "\nvel (Vertex)  : "<< std::setw(10) << vel_dist_vert << " (" << std::setw(10) << (comp_vel_vert-orig_vel_vert).norm() << ')'
            << " on vertex "<< std::setw(8) <<vel_vert->GetGID()<<", "<<vel_vert->GetCoord()
            << ": "<<comp_vel_vert<<" instead of "<<orig_vel_vert

            << "\nvel (Edge)    : "<< std::setw(10) << vel_dist_edge << " (" << std::setw(10) << (comp_vel_edge-orig_vel_edge).norm() << ')'
            << " on Edge   "<< std::setw(8) <<vel_edge->GetGID()<<", "<<GetBaryCenter(*vel_edge)
            << ": "<<comp_vel_edge<<" instead of "<<orig_vel_edge

            << "\nlset (Vertex) : "<< std::setw(10)<< lset_dist_vert << " (" << std::setw(10) << std::fabs(comp_lset_vert-orig_lset_vert) << ')'
            << " on vertex "<< std::setw(8) <<lset_vert->GetGID()<<", "<<lset_vert->GetCoord()
            << ": "<<comp_lset_vert<<" instead of "<<orig_lset_vert

            << "\nlset (Edge)   : "<< std::setw(10)<< lset_dist_edge << " (" << std::setw(10) << std::fabs(comp_lset_edge-orig_lset_edge) << ')'
            << " on Edge   "<< std::setw(8) <<lset_edge->GetGID()<<", "<<GetBaryCenter(*lset_edge)
            << ": "<<comp_lset_edge<<" instead of "<<orig_lset_edge
            << std::endl<<std::endl;

}


template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
{
    // Init velocity by a known function: Velocity
    Stokes.InitVel( &Stokes.v, Velocity);
    offset=0;
    MultiGridCL& mg= Stokes.GetMG();
    lset.Init(MovingDroplet);
    IdxDescCL p1Idx( P1_FE);
    p1Idx.CreateNumbering( mg.GetLastLevel(), mg, 0, &lset);
    VecDescCL p1(p1Idx), p2(p1Idx);
    InitPr( p1, mg, Pressure1);
    InitPr( p2, mg, Pressure2);
    P1toP1X( *Stokes.p.RowIdx, Stokes.p.Data, p1Idx, p1, p2, lset, mg);
}

template<typename Coeff>
  void SolveCoupledNS(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset,
                      AdapTriangCL& adap)
{
    // type of the problem
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    const MultiGridCL& mg=Stokes.GetMG();

    // output
    const bool adaptive=true;
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightT;
    typedef TwoPhaseVTKCL<StokesProblemT, LevelsetP2CL> VTKT;
    EnsightT ensight( mg, lset.Phi.RowIdx, Stokes, lset, C.ensDir, C.ensCase, C.geomName, adaptive,
                      (C.ensight? C.num_steps/C.ensight+1 : 0), C.binary, C.masterOut);
    VTKT     vtk    ( mg, Stokes, lset, (C.vtk ? C.num_steps/C.vtk+1 : 0),
                      std::string(C.vtkDir + "/" + C.vtkName), C.vtkBinary, true);

    if (C.ensight)
        ensight.write();
    if (C.vtk)
        vtk.write();

    for (int step= 1; step<=C.num_steps; ++step)
    {
#ifdef __SUNPRO_CC
# ifdef SAMPLE
        char sample[25];
        std::sprintf(sample,"sample_of_step_%i",step);
        collector_sample(sample);
# endif
#endif

        ParTimerCL time;
        if (ProcCL::IamMaster())
            std::cout << "=================================================================================== Schritt " << step << ":\n"
                      << "==> Solving coupled Levelset-Navier-Stokes problem ....\n"
                      << " Idx for vel  "<<Stokes.v.RowIdx->GetIdx()
                      << "\n Idx for pr   "<<Stokes.p.RowIdx->GetIdx()
                      << "\n Idx for lset "<<lset.Phi.RowIdx->GetIdx()<<std::endl;
        DisplayNumUnknowns(adap.GetMG(), Stokes.v);
        DisplayNumUnknowns(adap.GetMG(), Stokes.p);
        DisplayNumUnknowns(adap.GetMG(), lset.Phi);

        offset+=C.Anstroem;
        lset.Init(MovingDroplet);
        Stokes.UpdateXNumbering( &Stokes.pr_idx, lset);
        Stokes.UpdatePressure  ( &Stokes.p);
        double dummy1, dummy2;
        Point3DCL dummy3, dummy4, dummy5;
        Point3DCL bary;

        lset.GetInfo( dummy1, dummy2, bary, dummy3, Stokes.GetVelSolution(), dummy4, dummy5);
        IF_MASTER
                std::cout << "Position of the drop "<<bary[1]<<std::endl;


        if (ProcCL::IamMaster())
            std::cout << "==> Adaptive Refinement of MultiGrid"<<std::endl;

        adap.UpdateTriang( lset);
        if (C.ensight && step%C.ensight==0)
            ensight.write();
        if (C.vtk && step%C.vtk==0)
            vtk.write();

        time.Stop();
        double duration=time.GetMaxTime();
        CheckInterPol(Stokes, lset);
        int cplmem=GlobalSum(DDD_InfoCplMemory());
        Uint numdistobj=GlobalSum(Stokes.GetMG().GetNumDistributedObjects());
        if (ProcCL::IamMaster()){
             std::cout << "Memory for couplings: "<<cplmem<<"\nnumber of distributed objects "<<numdistobj<<'\n'
                       << "--> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
        }
    }
}

template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    MultiGridCL& mg= pmg.GetMG();
    sigma= Stokes.GetCoeff().SurfTens;

    ParTimerCL time;

    LevelsetP2CL lset( mg, &sigmaf, /*grad sigma*/ 0,
                       C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff, C.NarrowBand);
    lset.SetSurfaceForce( SF_ImprovedLB);

    AdapTriangCL adapt( pmg, lb, C.ref_width, 0, C.ref_flevel);

    MLIdxDescCL *vidx= &Stokes.vel_idx,
                *pidx= &Stokes.pr_idx;
    IdxDescCL   *lidx= &lset.idx;
    VecDescCL *v   = &Stokes.v,
              *p   = &Stokes.p,
              *l   = &lset.Phi;

    int cplmem;
    Uint numdistobj;
    double dur;

    offset=0;
    time.Reset();
    adapt.MakeInitialTriang(::DistanceFct1);
    time.Stop();
    dur=time.GetTime();
    cplmem=GlobalSum(DDD_InfoCplMemory());
    numdistobj=GlobalSum(mg.GetNumDistributedObjects());
    if (ProcCL::IamMaster()){
            std::cout << "Memory for couplings: "<<cplmem<<"\nnumber of distributed objects "<<numdistobj<<'\n'
                    << "--> Step 0 took "<<dur<<" sec."<<std::endl;
    }

    LevelsetRepairCL lsetrepair( lset, pmg);
    adapt.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes, pmg);
    adapt.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset, pmg);
    adapt.push_back( &prrepair);

    // Init levelset
    lset.CreateNumbering( mg.GetLastLevel(), lidx);
    l->SetIdx( lidx);
    lset.Init(MovingDroplet);

    Stokes.CreateNumberingVel( mg.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), pidx, 0, &lset);

    if (DROPS::ProcCL::IamMaster())
        std::cout << " - Distribution of elements of the multigrid of level "<<mg.GetLastLevel()<<"\n";
    mg.SizeInfo(std::cout);

    // Tell matrices and vectors about the numbering
    v->SetIdx( vidx);               p->SetIdx( pidx);
    Stokes.b.SetIdx( vidx);         Stokes.c.SetIdx( pidx);
    Stokes.A.SetIdx(vidx, vidx);    Stokes.B.SetIdx(pidx, vidx);
    Stokes.M.SetIdx(vidx, vidx);    Stokes.N.SetIdx(vidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx); Stokes.prA.SetIdx( pidx, pidx);

    //Setup initial problem
    if (ProcCL::IamMaster())
        std::cout << "=================================================================================== Init:\n"
                    << "==> Initialize Problem\n";

    InitProblemWithDrop(Stokes, lset);
    SolveCoupledNS(Stokes, lset, adapt);
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
  try
  {
    if (argc!=2)
    {
        IF_MASTER
          std::cout << "You have to specify one parameter:\n\t"
                    << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        IF_MASTER
          std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    IF_MASTER
      std::cout << C << std::endl;

    DROPS::ParTimerCL alltime;
    SetDescriber();
    DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
//     DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MSGSALL);

    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    int nx, ny, nz;
    std::string mesh( C.meshfile), delim("x@");
    size_t idx;
    while ((idx= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> dim_of_brick[0] >> dim_of_brick[1] >> dim_of_brick[2] >> nx >> ny >> nz;
    if (!brick_info)
    {
        std::cout << "error while reading geometry information: " << mesh << "\n";
        return 1;
    }

    C.r_inlet= dim_of_brick[0]/2;
    DROPS::Point3DCL orig, px, py, pz;
    px[0]= dim_of_brick[0]; py[1]= dim_of_brick[1]; pz[2]= dim_of_brick[2];

    DROPS::ParTimerCL time;
    DROPS::MGBuilderCL * mgb;
    if (DROPS::ProcCL::IamMaster())
        mgb = new DROPS::BrickBuilderCL(orig, px, py, pz, nx, ny, nz);
    else
        mgb = new DROPS::EmptyBrickBuilderCL(orig, px, py, pz);

    const bool bc[6]=
      {false, false, true, false, false, false};    // Rohr
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
      { &One, &One, &One, &Inflow, &One, &One};
//      { &Null, &Null, &Null, &Inflow, &Null, &Null};

    DROPS::MultiGridCL mg(*mgb);
    pmg.AttachTo(mg);

    DROPS::LoadBalHandlerCL lb(mg);
    lb.DoInitDistribution(DROPS::ProcCL::Master());
    switch (C.refineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
        case 3 : lb.SetStrategy(DROPS::Identity);  break;
    }

    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun), DROPS::P1X_FE /*default omit_bound*/);
    Strategy( prob, pmg, lb);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cout);
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
  return 0;
}

