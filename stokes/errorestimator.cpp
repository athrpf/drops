/// \file errorestimator.cpp
/// \brief stationary stokes problem with error estimator
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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
 * Copyright 2010 LNM/SC RWTH Aachen, Germany
*/

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construct the initial multigrid
#include "geom/geomselect.h"
#include "out/output.h"

 // include numeric computing!
#include "num/fe.h"
#include "num/solver.h"
#include "num/MGsolver.h"
#include "stokes/integrTime.h"
#include "num/stokessolverfactory.h"

 // include problem class
#include "stokes/stokes.h"
#include "misc/params.h"
#include "num/bndData.h"

//include coefficient class
#include "stokes/stokesCoeff.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "misc/bndmap.h"

using namespace std;

const char line[] ="----------------------------------------------------------------------------------\n";

void MarkLower( DROPS::MultiGridCL& mg, double tresh)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin()),
             ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
    {
        if (GetBaryCenter(*It)[2]<=tresh )
            It->SetRegRefMark();
    }
}

void
ZeroMean(DROPS::P1EvalCL< double,
                          const DROPS::StokesPrBndDataCL,
                          DROPS::VecDescCL>& f)
{
    const DROPS::Uint lvl= f.GetLevel();
    DROPS::MultiGridCL& mg= const_cast<DROPS::MultiGridCL&>( f.GetMG());
    double MV= 0., vol= 0., sum;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin( lvl),
         send= mg.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        sum= 0.;
        for(int i=0; i<4; ++i)
            sum+= f.val( *sit->GetVertex( i));
        sum/= 120;
        sum+= 2./15.*f.val( *sit, .25, .25, .25);
        MV+= sum * sit->GetVolume()*6.;
        vol+= sit->GetVolume();
    }
    const double c= MV/vol;
    std::cout << "\nconstant pressure offset: " << c << ", volume of domain: " << vol
              << std::endl;
    for (DROPS::MultiGridCL::TriangVertexIteratorCL sit= mg.GetTriangVertexBegin( lvl),
         send= mg.GetTriangVertexEnd( lvl); sit != send; ++sit) {
        f.SetDoF( *sit, f.val( *sit) - c);
    }
}

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&, double);

int
CheckVel(DROPS::P2EvalCL< DROPS::SVectorCL<3>,
         const DROPS::StokesVelBndDataCL,
         DROPS::VelVecDescCL>& fun,
         fun_ptr f)
{
    using namespace DROPS;
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const double t= fun.GetTime();
    const DROPS::Uint trilevel= fun.GetLevel();
    std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0.;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= (fun.val( *sit) - f( sit->GetCoord(), t)).norm();
        if ( std::abs( diff) > vmaxdiff) { ++ret; vmaxdiff= std::abs( diff); v= &*sit; }
    }
    std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin( trilevel),
         theend= mg.GetTriangEdgeEnd( trilevel); sit!=theend; ++sit) {
        diff = (fun.val( *sit, .5) - f( (sit->GetVertex( 0)->GetCoord()
                                        +sit->GetVertex( 1)->GetCoord())*0.5, t)).norm();
        if (std::abs( diff) > emaxdiff) { ++ret; emaxdiff= std::abs( diff); e= &*sit; }
    }
    {
        std::cout << "maximale Differenz Vertices: " << vmaxdiff << " auf\n";
        if (v) v->DebugInfo( std::cout);
        std::cout << "maximale Differenz Edges: " << emaxdiff << " auf\n";
        if (e) e->DebugInfo( std::cout);
        std::cout << std::endl;
    }
    return ret;
}

const double radiusorbit= 0.3; // Radius of the drops' orbit.
const double radiusdrop= 0.15; // Initial radius of the drop.

// positive outside the drop, negative inside the drop.
double
SignedDistToInterface(const DROPS::Point3DCL& p, double t)
{
   DROPS::Point3DCL c;
   c[0]= 0.5 + radiusorbit*std::cos( 2.*M_PI*t);
   c[1]= 0.5 + radiusorbit*std::sin( 2.*M_PI*t);
   c[2]= 0.5;
   return (p-c).norm() - radiusdrop;
}

typedef double (*signed_dist_fun)(const DROPS::Point3DCL& p, double t);

bool
ModifyGridStep(DROPS::MultiGridCL& mg,
               const signed_dist_fun Dist,
               const double width,         // Thickness of refined shell on each side of the interface
               const DROPS::Uint c_level,  // Outside the shell, use this level
               const DROPS::Uint f_level,  // Inside the shell, use this level
               const double t)             // Time of evaluation
// One step of grid change; returns true if modifications were necessary,
// false, if nothing changed.
{
    using namespace DROPS;
    bool shell_not_ready= false;
        for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
             theend= mg.GetTriangTetraEnd(); it!=theend; ++it) {
            double d= 1.;
            for (Uint j=0; j<4; ++j) {
                d= std::min( d, std::abs( Dist( it->GetVertex( j)->GetCoord(), t)));
            }
            const Uint l= it->GetLevel();
            if (d<=width) { // In the shell; level should be f_level.
                if (l < f_level) {
                    shell_not_ready= true;
                    it->SetRegRefMark();
                }
                else
                    if (l > f_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
            else { // Outside the shell; level should be c_level;
                if (l < c_level) { it->SetRegRefMark(); }
                else
                    if (l> c_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
        }
        mg.Refine();
    return shell_not_ready;
}

template<class Coeff>
void
UpdateTriangulation(DROPS::StokesP2P1CL<Coeff>& NS,
                    const signed_dist_fun Dist,
                    const double t,
                    const double width,         // Thickness of refined shell on each side of the interface
                    const DROPS::Uint c_level,  // Outside the shell, use this level
                    const DROPS::Uint f_level,  // Inside the shell, use this level
                    DROPS::VelVecDescCL* v1,
                    DROPS::VecDescCL* p1)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    Assert( c_level<=f_level, "UpdateTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    double time_v = v1->t, time_p = p1->t;
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= v1->RowIdx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx1= p1->RowIdx;
    IdxDescCL* pidx2= &loc_pidx;
    VelVecDescCL  loc_v;
    VecDescCL     loc_p;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p2= &loc_p;
    vidx2->SetFE( vecP2_FE);
    pidx2->SetFE( P1_FE);
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    const StokesBndDataCL& BndData= NS.GetBndData();
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i) {
        shell_not_ready= ModifyGridStep( mg, Dist, width, c_level, f_level, t);
        // Repair velocity
        std::swap( v2, v1);
        std::swap( vidx2, vidx1);
        match_fun match= mg.GetBnd().GetMatchFun();
        vidx1->CreateNumbering( mg.GetLastLevel(), mg, NS.GetBndData().Vel, match);
        if ( mg.GetLastLevel() != vidx2->TriangLevel()) {
            std::cout << "LastLevel: " << mg.GetLastLevel()
                      << " vidx2->TriangLevel: " << vidx2->TriangLevel() << std::endl;
            throw DROPSErrCL( "Strategy: Sorry, not yet implemented.");
        }
        v1->SetIdx( vidx1);
        P2EvalCL< SVectorCL<3>, const StokesVelBndDataCL,
                  const VelVecDescCL> funv2( v2, &BndData.Vel, &mg);
        RepairAfterRefineP2( funv2, *v1);
        v2->Clear( time_v);
        vidx2->DeleteNumbering( mg);
//P2EvalCL< SVectorCL<3>, const StokesVelBndDataCL,
//          VelVecDescCL> funv1( v1, &BndData.Vel, &mg, t);
//CheckVel( funv1, &MyPdeCL::LsgVel);
        // Repair pressure
        std::swap( p2, p1);
        std::swap( pidx2, pidx1);
        pidx1->CreateNumbering( mg.GetLastLevel(), mg, NS.GetBndData().Pr, match);
        p1->SetIdx( pidx1);
        typename StokesCL::const_DiscPrSolCL oldfunpr( p2, &BndData.Pr, &mg);
        RepairAfterRefineP1( oldfunpr, *p1);
        p2->Clear( time_p);
        pidx2->DeleteNumbering( mg);
    }
    // We want the solution to be where v1, p1 point to.
    if (v1 == &loc_v) {
        NS.vel_idx.GetFinest().swap( loc_vidx);
        NS.pr_idx.GetFinest().swap( loc_pidx);
        NS.v.SetIdx( &NS.vel_idx);
        NS.p.SetIdx( &NS.pr_idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
    }
    time.Stop();
    std::cout << "UpdateTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

void
MakeInitialTriangulation(DROPS::MultiGridCL& mg,
                         const signed_dist_fun Dist,
                         const double width,         // Thickness of refined shell on eache side of the interface
                         const DROPS::Uint c_level,  // Outside the shell, use this level
                         const DROPS::Uint f_level)  // Inside the shell, use this level
{
    using namespace DROPS;
    Assert( c_level<=f_level, "MakeInitialTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i)
        shell_not_ready=  ModifyGridStep( mg, Dist, width, c_level, f_level, 0.);
    time.Stop();
    std::cout << "MakeTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

typedef DROPS::StokesP2P1CL<DROPS::StokesFlowCoeffCL>
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;

DROPS::ParamCL P;

namespace DROPS // for Strategy
{

using ::MyStokesCL;

template< class StokesProblemT>
void Strategy( StokesProblemT& Stokes)
// flow control
{
    //Timer function
    TimerCL timer;

    //the triangulation
    MultiGridCL& MG= Stokes.GetMG();
    const typename MyStokesCL::BndDataCL::PrBndDataCL& PrBndData= Stokes.GetBndData().Pr;
    const typename MyStokesCL::BndDataCL::VelBndDataCL& VelBndData= Stokes.GetBndData().Vel;

    // connection triangulation and vectors
    // -------------------------------------------------------------------------
    std::cout << line << "Connecting triangulation and matrices/vectors ...\n";
    timer.Reset();

    Stokes.vel_idx.SetFE( vecP2_FE);
    Stokes.pr_idx.SetFE( P1_FE);

    //Modify Triangulation
    if( P.get("AdaptRef.ModifyGrid", 0) == 1)
        MakeInitialTriangulation( MG, &SignedDistToInterface, P.get<double>("AdaptRef.Width"), P.get<int>("AdaptRef.CoarsestLevel"), P.get<int>("AdaptRef.FinestLevel"));

    // solve the linear equation system
    // -------------------------------------------------------------------------
    std::cout << line << "Solve the linear equation system ...\n";

    // type of preconditioner and solver
    StokesSolverFactoryCL<StokesProblemT>         factory( Stokes, P);
    StokesSolverFactoryObsoleteCL<StokesProblemT> obsoletefactory( Stokes, P);
    StokesSolverBaseCL* solver = (P.get<int>("Stokes.StokesMethod")< 500000) ? factory.CreateStokesSolver() : obsoletefactory.CreateStokesSolver();

    MLIdxDescCL  loc_vidx, loc_pidx;
    MLIdxDescCL* vidx1= &Stokes.vel_idx;
    MLIdxDescCL* pidx1= &Stokes.pr_idx;
    MLIdxDescCL* vidx2= &loc_vidx;
    MLIdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &Stokes.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &Stokes.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL  cplM(vidx1);

    vidx2->SetFE( vecP2_FE);
    pidx2->SetFE( P1_FE);

    int step= 0;
    StokesDoerflerMarkCL<typename MyStokesCL::est_fun, MyStokesCL>
        Estimator( P.get<double>("Error.RelReduction"), P.get<double>("Error.Threshold"), P.get<double>("Error.Meas"), true, &MyStokesCL::ResidualErrEstimator, Stokes);
    bool new_marks= false;

    do
    {
        if( P.get<double>("Error.MarkLower") != 0)
        	MarkLower( MG, P.get<double>("Error.MarkLower"));

        MG.Refine();

        if( StokesSolverFactoryHelperCL().VelMGUsed(P) || StokesSolverFactoryObsoleteHelperCL().VelMGUsed(P)){
            Stokes.SetNumVelLvl( MG.GetNumLevel());
            vidx1->resize( MG.GetNumLevel(), vecP2_FE);
        }
        if( StokesSolverFactoryHelperCL().PrMGUsed(P) || StokesSolverFactoryObsoleteHelperCL().PrMGUsed(P)){
            Stokes.SetNumPrLvl( MG.GetNumLevel());
            pidx1->resize( MG.GetNumLevel(), P1_FE);
        }

        Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx1);
        Stokes.CreateNumberingPr ( MG.GetLastLevel(), pidx1);
        std::cout << "old and new TriangLevel: " << vidx2->TriangLevel() << ", "
                  << vidx1->TriangLevel() << std::endl;
        MG.SizeInfo(std::cout);

        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        std::cout << "Number of pressure unknowns: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cout << "Number of velocity unknowns " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;

        if( StokesSolverFactoryHelperCL().VelMGUsed(P) || StokesSolverFactoryObsoleteHelperCL().VelMGUsed(P))
        {
            MLMatrixCL* PVel = (P.get<int>("Stokes.StokesMethod")< 500000) ? factory.GetPVel() : obsoletefactory.GetPVel();
            SetupP2ProlongationMatrix( MG, *PVel, vidx1, vidx1);
        }

        if( StokesSolverFactoryHelperCL().PrMGUsed(P) || StokesSolverFactoryObsoleteHelperCL().PrMGUsed(P))
        {
            MLMatrixCL* PPr = (P.get<int>("Stokes.StokesMethod") < 500000) ? factory.GetPPr() : obsoletefactory.GetPPr();
            SetupP1ProlongationMatrix( MG, *PPr, pidx1, pidx1);
        }

        if (P.get<int>("Stokes.StokesMethod") < 500000) {
            factory.SetMatrixA( &Stokes.A.Data.GetFinest());
            factory.SetMatrices( &Stokes.A.Data, &Stokes.B.Data, &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
        }

        if ( v2->RowIdx)
        {
            P2EvalCL<Point3DCL, const MyStokesCL::BndDataCL::VelBndDataCL, const VecDescCL>  oldx( v2, &VelBndData, &MG);
            VelVecDescCL loc_v;
            IdxDescCL loc_vidx( vecP2_FE);

            loc_vidx.CreateNumbering( vidx2->TriangLevel(), MG, VelBndData);
            std::cout << "vidx2->TriangLevel(): " << vidx2->TriangLevel() << std::endl;
            loc_v.SetIdx( &loc_vidx);
            RepairAfterRefineP2( oldx, loc_v);

            v1->Clear( 0.0);
            vidx2->DeleteNumbering( MG);

            vidx2->GetFinest().swap( loc_vidx);
            v2->SetIdx( vidx2);
            v2->Data= loc_v.Data;
            P2EvalCL<Point3DCL, const MyStokesCL::BndDataCL::VelBndDataCL, VecDescCL>        newx( v1, &VelBndData, &MG);

//            Interpolate( newx, oldx);

            v2->Reset();
         }

        if ( p2->RowIdx)
        {
            P1EvalCL<double, const MyStokesCL::BndDataCL::PrBndDataCL, const VecDescCL>  oldx( p2, &PrBndData, &MG);
            P1EvalCL<double, const MyStokesCL::BndDataCL::PrBndDataCL, VecDescCL>        newx( p1, &PrBndData, &MG);
  //          RepairAfterRefineP1( oldx, *p1);
  //          Interpolate( newx, oldx);
            p1->Clear(0.0);
            p2->Reset();
         }

        Stokes.A.SetIdx  ( vidx1, vidx1);
        Stokes.M.SetIdx  ( vidx1, vidx1);
        Stokes.B.SetIdx  ( pidx1, vidx1);
        Stokes.prM.SetIdx( pidx1, pidx1);
        Stokes.prA.SetIdx( pidx1, pidx1);
        Stokes.b.SetIdx  ( vidx1);
        Stokes.c.SetIdx  ( pidx1);
        cplM.SetIdx( vidx1);
        timer.Reset();
        timer.Start();

        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &cplM, 0.0);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, 0.0);
        Stokes.SetupPrMass( &Stokes.prM);
        Stokes.SetupPrStiff( &Stokes.prA);


        timer.Stop();

        std::cout << "SetupSystem: " << timer.GetTime() << " seconds." << std::endl;
        timer.Reset();

        if( P.get<std::string>("StokesCoeff.Solution_Vel").compare("None")!=0)  // check whether solution is given
            Stokes.GetDiscError( StokesFlowCoeffCL::LsgVel, StokesFlowCoeffCL::LsgPr);
        timer.Reset();

        double err0= norm_sq( Stokes.A.Data*v1->Data + transp_mul( Stokes.B.Data, p1->Data) - Stokes.b.Data)
                    +norm_sq( Stokes.B.Data*v1->Data - Stokes.c.Data);
        std::cout << "000 residual: " << std::sqrt( err0) << std::endl;

        timer.Start();
        solver->Solve( Stokes.A.Data, Stokes.B.Data, v1->Data, p1->Data, Stokes.b.Data, Stokes.c.Data);
        timer.Stop();
        double err= norm_sq( Stokes.A.Data*v1->Data + transp_mul( Stokes.B.Data, p1->Data) - Stokes.b.Data)
                    +norm_sq( Stokes.B.Data*v1->Data - Stokes.c.Data);
        std::cout << "Solver: "<<timer.GetTime()<<" seconds.\n";
        std::cout << "iter: " << solver->GetIter() << "\tresid: " << solver->GetResid() << std::endl;
        std::cout << "000 residual: " << std::sqrt( err)/std::sqrt( err0) << std::endl;

        if( P.get<std::string>("StokesCoeff.Solution_Vel").compare("None")!=0)  // check whether solution is given
          Stokes.CheckSolution( v1, p1, StokesFlowCoeffCL::LsgVel, StokesFlowCoeffCL::DLsgVel, StokesFlowCoeffCL::LsgPr, true);

        if( step==0 && P.get<int>("Error.DoErrorEstimate"))
            Estimator.Init(typename MyStokesCL::const_DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::const_DiscVelSolCL(v1, &VelBndData, &MG));

        timer.Reset();
        timer.Start();

        if ( P.get<int>("Error.DoErrorEstimate"))
            new_marks= Estimator.Estimate(typename MyStokesCL::const_DiscPrSolCL(p1, &PrBndData, &MG), typename MyStokesCL::const_DiscVelSolCL(v1, &VelBndData, &MG) );
        timer.Stop();
        std::cout << "Estimation: " << timer.GetTime() << " seconds.\n";
        Stokes.A.Reset();
        Stokes.M.Reset();
        Stokes.prA.Reset();
        Stokes.prM.Reset();
        Stokes.B.Reset();
        Stokes.b.Reset();
        Stokes.c.Reset();

        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
    } while (++step< P.get<int>("Error.NumRef") && new_marks);

    // we want the solution to be in Stokes.v, Stokes.pr
    if (v2 == &loc_v)
    {
        Stokes.vel_idx.swap( loc_vidx);
        Stokes.pr_idx.swap( loc_pidx);
        Stokes.v.SetIdx( &Stokes.vel_idx);
        Stokes.p.SetIdx( &Stokes.pr_idx);
        Stokes.v.Data= loc_v.Data;
        Stokes.p.Data= loc_p.Data;
    }

    delete solver;
}

} // end of namespace DROPS

int main ( int argc, char** argv)
{
    try
    {
        std::ifstream param;
        if (argc!=2)
        {
            std::cout << "Using default parameter file: drivcav.json\n";
            param.open( "MGsdropsP2.json");
        }
        else
            param.open( argv[1]);
        if (!param)
        {
            std::cerr << "error while opening parameter file\n";
            return 1;
        }
        param >> P;
        param.close();
        std::cout << P << std::endl;

        // Check MarkLower value
        if( P.get<int>("DomainCond.GeomType") == 0) P.put("Error.MarkLower", 0);
        else {
          int nx, ny, nz;
          double dx, dy, dz;
          std::string mesh( P.get<string>("DomainCond.MeshFile")), delim("x@");
          size_t idx;
          while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
          std::istringstream brick_info( mesh);
          brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
          if (P.get("Error.MarkLower", 0)<0 || P.get("Error.MarkLower", 0) > dy)
          {
        	  std::cerr << "Wrong value of MarkLower\n";
        	  return 1;
          }
        }

        // time measurement
        DROPS::TimerCL timer;

        // set up data structure to represent a Stokes problem
        // ---------------------------------------------------------------------
        std::cout << line << "Set up data structure to represent a Stokes problem ...\n";
        timer.Reset();

        //create geometry
        DROPS::MultiGridCL* mg= 0;
        DROPS::StokesBndDataCL* bdata = 0;

        //only for measuring cell, not used here
        double r = 1;
        std::string serfile = "none";

        DROPS::BuildDomain( mg, P.get<string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), serfile, r);
        DROPS::BuildBoundaryData( mg, bdata, P.get<string>("DomainCond.BoundaryType"), P.get<string>("DomainCond.BoundaryFncs"));

        // Setup the problem
        DROPS::StokesFlowCoeffCL tmp = DROPS::StokesFlowCoeffCL( P);
        StokesOnBrickCL prob(*mg, tmp, *bdata);
        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;

        // Refine the grid
        // ---------------------------------------------------------------------
        std::cout << "Refine the grid " << P.get<int>("DomainCond.InitialCond") << " times regulary ...\n";
        timer.Reset();

        // Create new tetrahedra
        for ( int ref=1; ref<=P.get<int>("DomainCond.InitialCond"); ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( *mg);
            mg->Refine();
        }

        timer.Stop();
        std::cout << " o time " << timer.GetTime() << " s" << std::endl;
        mg->SizeInfo(std::cout);

        // Solve the problem
        DROPS::Strategy( prob); //do all the stuff
        std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;

        double min= prob.p.Data.min(),
               max= prob.p.Data.max();
        std::cout << "pressure min/max: "<<min<<", "<<max<<std::endl;

        delete mg;
        delete bdata;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
