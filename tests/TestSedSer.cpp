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

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"
#include "partests/params.h"

 // include numeric computing!
#include "levelset/adaptriang.h"

 // include in- and output
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
    const DROPS::Ulint acc_num_unk = x.Data.size();
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    std::cout << "  + Number of DOF of index "<<idx<<" (accumulated/global):  "
              <<acc_num_unk<< std::endl;
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
    const double x = (p[0])*(p[0])-1,
                 z = (p[2])*(p[2])-1;

//     ret[1]= x * z * C.Anstroem;// * (1-ampl*std::cos(2*M_PI*freq*t));  // Rohr
    ret[1]= x * z;
    ret[0]=1.; ret[2]=1.;
    return ret;
}

// droplet
double DistanceFct1( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= C.Mitte-p;
    return d.norm()-C.Radius;
}

namespace DROPS{

double offset=0;
double MovingDroplet(const Point3DCL& p)
{
    Point3DCL Mitte=C.Mitte; Mitte[1]+=offset;
    const Point3DCL d=Mitte-p;
    return d.norm()-C.Radius;
}

double Pressure(const Point3DCL p, double=0.)
{
    return p[1]/0.03;
}

Point3DCL Velocity(const Point3DCL& p, double=0.)
{
    return Inflow(p,0);
}

void InitPr(VecDescCL& pr, const MultiGridCL& mg)
{
    VectorCL& lsg= pr.Data;
    Uint lvl     = pr.GetLevel(),
         idx     = pr.RowIdx->GetIdx();

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(mg).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(mg).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
        if (sit->Unknowns.Exist(idx))
            lsg[sit->Unknowns(idx)]=Pressure(sit->GetCoord());
}

template<class Coeff>
  void CheckInterPol(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
{
    double vel_dist=-1, pr_dist=-1, lset_dist_vert=-1, lset_dist_edge=-1,
           orig_val_vert, orig_val_edge, comp_val_vert, comp_val_edge,
           dist;

    char dat[30];
    std::sprintf(dat,"output/difference.diff");
    static std::ofstream check(dat);
    static int time=0;


    typename InstatNavierStokes2PhaseP2P1CL<Coeff>::const_DiscPrSolCL  funP=Stokes.GetPrSolution();
    typename InstatNavierStokes2PhaseP2P1CL<Coeff>::const_DiscVelSolCL funV=Stokes.GetVelSolution();
    LevelsetP2CL::const_DiscSolCL                                      funL=lset.GetSolution();

    const Uint pr_idx  = funP.GetSolution()->RowIdx->GetIdx();
    const Uint vel_idx = funV.GetSolution()->RowIdx->GetIdx();
    const Uint lset_idx= funL.GetSolution()->RowIdx->GetIdx();

    const VertexCL* max_vert;
    const EdgeCL*   max_edge;

    const MultiGridCL& mg=funP.GetMG();

    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()); sit!=mg.GetTriangVertexEnd(); ++sit){
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(pr_idx)){
            dist =    std::fabs(funP.val(*sit)-Pressure(sit->GetCoord()))
                    / std::min(funP.val(*sit),Pressure(sit->GetCoord()));
            if (dist>pr_dist) pr_dist=dist;
        }
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(vel_idx)){
            dist =    (funV.val(*sit)-Velocity(sit->GetCoord())).norm()
                    / std::min(funV.val(*sit).norm(),Velocity(sit->GetCoord()).norm());
            if (dist>vel_dist) vel_dist=dist;
        }
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(lset_idx)){
            dist =    std::fabs(funL.val(*sit)-MovingDroplet(sit->GetCoord()))
                    / std::min(funL.val(*sit),MovingDroplet(sit->GetCoord()));
            if (dist>lset_dist_vert) {
                lset_dist_vert=dist;
                orig_val_vert=MovingDroplet(sit->GetCoord());
                comp_val_vert=funL.val(*sit);
                max_vert=&(*sit);
            }
        }
    }

    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()); sit!=mg.GetTriangEdgeEnd(); ++sit){
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(vel_idx)){
            dist =    (funV.val(*sit)-Velocity(GetBaryCenter(*sit))).norm()
                    / std::min(funV.val(*sit).norm(),Velocity(GetBaryCenter(*sit)).norm());
            if (dist>vel_dist) vel_dist=dist;
        }
        if (sit->Unknowns.Exist() && sit->Unknowns.Exist(lset_idx)){
            dist =    std::fabs(funL.val(*sit)-MovingDroplet(GetBaryCenter(*sit)))
                    / std::min(funL.val(*sit),MovingDroplet(GetBaryCenter(*sit)));
            if (dist>lset_dist_edge) {
                lset_dist_edge=dist;
                orig_val_edge=MovingDroplet(GetBaryCenter(*sit));
                comp_val_edge=funL.val(*sit);
                max_edge=&(*sit);
            }
        }
    }
    check  << "Diferences (at time "<<(time++)<<"):"
           << "\npr            : "<<pr_dist
           << "\nvel           : "<<vel_dist
           << "\nlset (Vertex) : "<<lset_dist_vert<<" on vertex "<<max_vert->GetCoord()
           << ": "<<comp_val_vert<<" instead of "<<orig_val_vert
           << "\nlset (Edge)   : "<<lset_dist_edge<<" on Edge   "<<GetBaryCenter(*max_edge)
           << ": "<<comp_val_edge<<" instead of "<<orig_val_edge
           <<std::endl<<std::endl;
}


template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset)
{
    // Initial velocity is zero
    Stokes.InitVel( &Stokes.v, Velocity);
    offset=0;
    lset.Init(MovingDroplet);
    InitPr(Stokes.p, Stokes.GetMG());
}

template<typename Coeff>
  void SolveCoupledNS(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset,
                      AdapTriangCL& adap)
{
    // type of the problem
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    for (int step= 1; step<=C.num_steps; ++step)
    {
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

        std::cout << "==> Adaptive Refinement of MultiGrid"<<std::endl;

        adap.UpdateTriang( Stokes, lset);
        CheckInterPol(Stokes, lset);
    }
}

template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, MultiGridCL& mg)
{
    LevelsetP2CL lset( mg, Stokes.GetCoeff().SurfTens,
                       C.lset_theta, C.lset_SD, C.RepDiff, C.lset_iter, C.lset_tol, C.CurvDiff);

    AdapTriangCL adapt( mg, C.ref_width, 0, C.ref_flevel);
    adapt.MakeInitialTriang(::DistanceFct1);

    std::cout << " - Distribution of elements of the multigrid of level "<<mg.GetLastLevel()<<"\n";
    mg.SizeInfo(std::cout);

    IdxDescCL *vidx= &Stokes.vel_idx,
              *pidx= &Stokes.pr_idx,
              *lidx= &lset.idx;
    VecDescCL *v   = &Stokes.v,
              *p   = &Stokes.p,
              *l   = &lset.Phi;

    Stokes.CreateNumberingVel( mg.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), pidx);
    lset.CreateNumbering     ( mg.GetLastLevel(), lidx);

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
    InitProblemWithDrop(Stokes, lset);

    SolveCoupledNS(Stokes, lset, adapt);

}
} // end of namespace DROPS


int main (int argc, char** argv)
{
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

    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::Point3DCL e1(0.), e2(0.), e3(0.), orig(0.);
//     e1[0]=e3[2]=2*C.r_inlet; e2[1]=10;
//     orig[0]=orig[2]=-C.r_inlet;
    e1[0]=e3[2]=1.; e2[1]=5.;

    DROPS::MGBuilderCL * mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, 4,20,4);

    const bool bc[6]=
      {false, false, true, false, false, false};    // Rohr
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
      { &Null, &Null, &Null, &Inflow, &Null, &Null};

    DROPS::MultiGridCL mg(*mgb);

    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun));

    Strategy( prob, mg);    // do all the stuff

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

