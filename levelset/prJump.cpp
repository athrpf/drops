//**************************************************************************
// File:    prJump.cpp    (cf. surfTens.cpp)                               *
// Content: test FE spaces for pressure                                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "num/MGsolver.h"
#include "num/solver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/adaptriang.h"
#include "levelset/params.h"
#include "levelset/mzelle_hdr.h"
#include "num/stokessolverfactory.h"
#include <fstream>
#include <iomanip>
#include <vector>


DROPS::ParamMesszelleNsCL C;

// program for testing various FE pressure spaces
// using a special constant surface force:
//     \sigma \int_\Gamma v n ds
// => implies a constant pressure jump [p] = \sigma across \Gamma.

// rho*du/dt - mu*laplace u + Dp = f + rho*g - ovn
//                        -div u = 0
//                             u = u0, t=t0

//   - different interfaces by choosing DistanceFct (planar/sperical)
//     -> set avg_ex: average of exact pressure solution
//   - different surface forces by applying lset.SetSurfaceForce(...)
//     -> set prJump: height of pressure jump

/*
double DistanceFct( const DROPS::Point3DCL& p)
{ // plane perpendicular to n=PosDrop with distance Radius[0] from origin.
    return inner_prod( C.Mitte/norm(C.Mitte), p) - C.Radius[0];
}
*/

namespace DROPS // for Strategy
{

/*=====================================================================================================
init testcase
=====================================================================================================*/


void InitPr( VecDescCL& p, double delta_p, const MultiGridCL& mg, const FiniteElementT prFE, const ExtIdxDescCL& Xidx)
{
    const Uint lvl= p.RowIdx->TriangLevel(),
        idxnum= p.RowIdx->GetIdx();

    delta_p/= 2;
    switch (prFE)
    {
      case P0_FE:
        for( MultiGridCL::const_TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl),
            end= mg.GetTriangTetraEnd(lvl); it != end; ++it)
        {
            const double dist= EllipsoidCL::DistanceFct( GetBaryCenter( *it));
            p.Data[it->Unknowns(idxnum)]= dist > 0 ? -delta_p : delta_p;
        }
        break;
      case P1X_FE:
        for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
            end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
        {
            const IdxT idx= it->Unknowns(idxnum);
            if (Xidx[idx]==NoIdx) continue;
            p.Data[Xidx[idx]]= -2*delta_p; // jump height
        }
      case P1_FE: // and P1X_FE
        for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
            end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
        {
            const double dist= EllipsoidCL::DistanceFct( it->GetCoord());
            p.Data[it->Unknowns(idxnum)]= InterfacePatchCL::Sign(dist)==1 ? -delta_p : delta_p;
        }
        break;
      default:
        std::cout << "InitPr not implemented yet for this FE type!\n";
    }
}

double my_abs( double x) { return std::abs(x); }

void L2ErrorPr( const VecDescCL& p, const LevelsetP2CL& lset, const MatrixCL& prM, double delta_p, const MultiGridCL& mg, const FiniteElementT prFE, const ExtIdxDescCL& Xidx, double p_ex_avg)
{
    const double min= p.Data.min(), max= p.Data.max();
    std::cout << "pressure min/max/diff:\t" << min << "\t" << max << "\t" << (max-min-delta_p) << "\n";

    VectorCL ones( 1.0, p.Data.size());
    if (prFE==P1X_FE)
        for (int i=Xidx.GetNumUnknownsStdFE(), n=ones.size(); i<n; ++i)
            ones[i]= 0;
    const double Vol= dot( prM*ones, ones)*C.mat_ViscDrop; // note that prM is scaled by 1/mu !!
// std::cout << "Vol = " << Vol << '\n';
    const double p_avg= dot( prM*p.Data, ones)*C.mat_ViscDrop/Vol; // note that prM is scaled by 1/mu !!
    VectorCL diff( p.Data - p_avg*ones);
    const double p0_avg= dot( prM*diff, ones)*C.mat_ViscDrop/Vol;
    std::cout << "average of pressure:\t" << p_avg << std::endl;
    std::cout << "avg. of scaled pr:\t" << p0_avg << std::endl;

    if (prFE==P1X_FE)
    {
        VecDescCL p_exakt( p.RowIdx);
        InitPr( p_exakt, delta_p, mg, prFE, Xidx);
        const double p_ex_avg2= dot( prM*p_exakt.Data, ones)*C.mat_ViscDrop/Vol;
        std::cout << "avg. of exact pr:\t" << p_ex_avg2 << std::endl;
        diff-= VectorCL( p_exakt.Data - p_ex_avg*ones);
        const double L2= std::sqrt( C.mat_ViscDrop*dot( prM*diff, diff));
        std::cout << "*************\n"
                  << "assuming avg(p*)==" << p_ex_avg
                  << "  ===>  \t||e_p||_L2 = " << L2 << std::endl
                  << "*************\n";
        return;
    }
    IdxT prNumb[4];
    InterfacePatchCL cut;
    double L2= 0, L1= 0;
    const Uint lvl= p.RowIdx->TriangLevel();
    for (MultiGridCL::const_TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(lvl),
         send= mg.GetTriangTetraEnd(lvl); sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6.;
        cut.Init( *sit, lset.Phi);
        const bool nocut= !cut.Intersects();

        LocalP2CL<> p0;
        if (prFE==P0_FE)
            p0= p.Data[sit->Unknowns(p.RowIdx->GetIdx())] - p_avg;
        else
        {
            GetLocalNumbP1NoBnd( prNumb, *sit, *p.RowIdx);
            for (int i=0; i<4; ++i)
                p0[i]= p.Data[prNumb[i]] - p_avg;
            for (int j=0; j<6; ++j)
                p0[j+4]= 0.5*(p0[VertOfEdge(j,0)] + p0[VertOfEdge(j,1)]);
        }

        if (nocut)
        {
            Quad2CL<> diff( p0);
            diff-= (EllipsoidCL::DistanceFct( GetBaryCenter(*sit)) >0 ? -delta_p/2 : delta_p/2) - p_ex_avg;
            L2+= Quad2CL<>( diff*diff).quad( absdet);
            diff.apply( my_abs);
            L1+= diff.quad( absdet);
        }
        else
        {
            LocalP2CL<> diffpos( p0 + delta_p/2 + p_ex_avg), diffneg( p0 - delta_p/2 + p_ex_avg);
            for (int ch=0; ch<8; ++ch)
            {
                cut.ComputeCutForChild(ch);
                L2+= cut.quad( LocalP2CL<>(diffpos*diffpos), absdet, true); // integrate on positive part
                L2+= cut.quad( LocalP2CL<>(diffneg*diffneg), absdet, false); // integrate on negative part
                diffpos.apply( my_abs);
                diffneg.apply( my_abs);
                L1+= cut.quad( diffpos, absdet, true); // integrate on positive part
                L1+= cut.quad( diffneg, absdet, false); // integrate on negative part
            }
        }
    }
    L2= std::sqrt(L2);
    std::cout << "*************\n"
              << "assuming avg(p*)==" << p_ex_avg
              << "  ===>  \t||e_p||_L2 = " << L2 << std::endl
              << "                           \t||e_p||_L1 = " << L1 << std::endl
              << "*************\n";
}

void PostProcessPr( const VecDescCL& p, VecDescCL& new_p, const MultiGridCL& mg)
// as Ensight output is for cont. data only, put values of P0 DoFs in tetras
// into P1 DoFs in vertices (s.t. the average over each tetra gives
// the former P0 value)
{
    VectorCL num( new_p.Data.size()), sum( new_p.Data.size());

    const Uint lvl= p.RowIdx->TriangLevel(),
        idxnum= p.RowIdx->GetIdx(),
        idxnum2= new_p.RowIdx->GetIdx();

    if (p.RowIdx->NumUnknownsTetra())
        for( MultiGridCL::const_TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl),
            end= mg.GetTriangTetraEnd(lvl); it != end; ++it)
        {
            const double val= p.Data[it->Unknowns(idxnum)];
            for (int i=0; i<4; ++i)
            {
                const IdxT nr2= it->GetVertex(i)->Unknowns(idxnum2);
                sum[nr2]+= val;
                num[nr2]+= 1;
            }
        }

    for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
        end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
    {
        const IdxT nr= it->Unknowns(idxnum),
            nr2= it->Unknowns(idxnum2);
        double val= p.Data[nr];
        new_p.Data[nr2]= p.RowIdx->NumUnknownsTetra() ? val + sum[nr2]/num[nr2]
                                                      : val;
    }
}

void PrintNorm( string name, const VectorCL& v)
{
    std::cout << name << ":\t2-norm: "
        << norm( v) << "\tmax: " << supnorm( v) << std::endl;
}

template<class Coeff>
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adap)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    sigma= C.sft_SurfTension;
    // Levelset-Disc.: Crank-Nicholson
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.lvs_SD, C.lvs_CurvDiff);

//    lset.SetSurfaceForce( SF_LB);
    lset.SetSurfaceForce( SF_ImprovedLB);
//    lset.SetSurfaceForce( SF_Const);
    const double Vol= 8.,
//        prJump= C.sigma, // for SF_Const force
        prJump= C.sft_SurfTension*2/C.exp_RadDrop[0], // for SF_*LB force
        avg_ex= prJump/2.*(8./3.*M_PI*C.exp_RadDrop[0]*C.exp_RadDrop[0]*C.exp_RadDrop[0] - Vol)/Vol; // for spherical interface
//        avg_ex= 0; // for planar interface

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    if ( StokesSolverFactoryHelperCL<ParamMesszelleNsCL>().VelMGUsed(C))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL<ParamMesszelleNsCL>().PrMGUsed(C))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());

    VecDescCL new_pr;  // for pressure output in Ensight

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( MG.GetLastLevel(), pidx, NULL, &lset);
    MG.SizeInfo( std::cout);
    Stokes.SetIdx();
    Stokes.v.SetIdx(vidx);
    Stokes.p.SetIdx(pidx);
    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << lset.Phi.Data.size() << " levelset unknowns.\n";

    new_pr.SetIdx( lidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset); // makes no sense for P0

    SSORPcCL ssor;
    PCG_SsorCL PCG( ssor, C.stk_InnerIter, C.stk_InnerTol); // for computing curvature force

    VelVecDescCL curv( vidx);
    VelVecDescCL curvForce( vidx);

    switch (C.dmc_InitialCond)
    {
      case 2:
      {
        // solve initial velocities for given pressure field
        TimerCL time;
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cout << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        InitPr( Stokes.p, prJump, MG, Stokes.GetPrFE(), Stokes.GetXidx());
        VectorCL surf( Stokes.b.Data + curv.Data), BTp( transp_mul( Stokes.B.Data, Stokes.p.Data));
        PrintNorm( "surf. force", curv.Data);
        PrintNorm( "BT p", BTp);
        PrintNorm( "Diff.", VectorCL(curv.Data - BTp));
        std::cout << "Solving velocity for exact pressure given...\n";
        PCG.Solve( Stokes.A.Data, Stokes.v.Data, VectorCL( curv.Data - transp_mul( Stokes.B.Data, Stokes.p.Data)) );
      } break;

    case 1:
      {
        // solve stationary problem for initial velocities
        TimerCL time;
        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        time.Stop();
        std::cout << "Discretizing Stokes/Surf.Force for initial velocities took "<<time.GetTime()<<" sec.\n";

        //InitPr( Stokes.p, prJump, MG, Stokes.GetPrFE(), Stokes.GetXidx().GetFinest());
        time.Reset();

        StokesSolverFactoryCL<StokesProblemT, ParamMesszelleNsCL> stokessolverfactory(Stokes, C);
        StokesSolverBaseCL* solver = stokessolverfactory.CreateStokesSolver();

        // initializes prolongation matrices
        UpdateProlongationCL PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
        adap.push_back( &PVel);
        UpdateProlongationCL PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
        adap.push_back( &PPr);

        // for MinComm
        stokessolverfactory.SetMatrixA( &Stokes.A.Data.GetFinest());
        //for Stokes-MGM
        stokessolverfactory.SetMatrices( &Stokes.A.Data.GetCoarsest(), &Stokes.B.Data.GetCoarsest(),
                                         &Stokes.M.Data.GetCoarsest(), &Stokes.prM.Data.GetCoarsest(), &Stokes.pr_idx.GetCoarsest());

        solver->Solve( Stokes.A.Data, Stokes.B.Data, Stokes.v.Data, Stokes.p.Data,
                       curv.Data, Stokes.c.Data);
        std::cout << "iter: " << solver->GetIter()
                  << "\tresid: " << solver->GetResid() << std::endl;
        time.Stop();
        std::cout << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
        delete solver;
      }
    }

    const VectorCL& u= Stokes.v.Data;
    std::cout << "\n----------------\n || u ||_oo = " << supnorm(u)
              << "\n || u ||_M  = " << std::sqrt( dot( Stokes.M.Data*u, u))
              << "\n || u ||_A  = " << std::sqrt( dot( Stokes.A.Data*u, u))
              << "\n----------------\n";
    if (Stokes.UsesXFEM())
    {
        const ExtIdxDescCL& Xidx= Stokes.GetXidx();
        const size_t n= Stokes.p.Data.size();

        const double limtol= 10,
            lim_min= -prJump - limtol*prJump,
            lim_max= -prJump + limtol*prJump;
        double xmin= 1e99, xmax= -1e99, sum= 0, sum_lim= 0;
        IdxT num= 0;
        for (size_t i=Xidx.GetNumUnknownsStdFE(); i<n; ++i)
        {
            const double pr= Stokes.p.Data[i];
            sum+= pr;
            ++num;
            if (pr>xmax) xmax= pr;
            if (pr<xmin) xmin= pr;
            if (pr>lim_max) Stokes.p.Data[i]= lim_max;
            if (pr<lim_min) Stokes.p.Data[i]= lim_min;
            sum_lim+= Stokes.p.Data[i];
        }
        std::cout << "extended pr: min/max/avg = " << xmin << ", " << xmax << ", " << sum/num << std::endl;
        std::cout << "limited pr:  min/max/avg = " << lim_min << ", " << lim_max << ", " << sum_lim/num << std::endl;
    }

    L2ErrorPr( Stokes.p, lset, Stokes.prM.Data.GetFinest(), prJump, MG, Stokes.GetPrFE(), Stokes.GetXidx(), avg_ex);

    PostProcessPr( Stokes.p, new_pr, MG);

    // Initialize Ensight6 output
    std::string ensf( C.ens_EnsDir + "/" + C.ens_EnsCase);
    Ensight6OutCL ensight( C.ens_EnsCase + ".case", C.tm_NumSteps + 1);
    ensight.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(), C.ens_GeomName,           ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),             "Levelset",  ensf + ".scl"));
    ensight.Register( make_Ensight6Scalar( Stokes.GetPrSolution( new_pr),  "Pressure",  ensf + ".pr"));
    ensight.Register( make_Ensight6Vector( Stokes.GetVelSolution(),        "Velocity",  ensf + ".vel"));
    ensight.Register( make_Ensight6Vector(  Stokes.GetVelSolution( curvForce),
                                                                           "Curvature", ensf + ".crv"));
    if (Stokes.UsesXFEM())
        ensight.Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p,  "XPressure", ensf + ".pr"));

    if (C.ens_EnsightOut) ensight.Write();

    std::cout << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc>2)
    {
        std::cout << "You have to specify at most one parameter:\n\t"
                  << argv[0] << " [<param_file>]" << std::endl;
        return 1;
    }
    std::ifstream param;
    if (argc>1)
        param.open( argv[1]);
    else
        param.open( "prJump.param");
    if (!param)
    {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cout << C << std::endl;

    typedef DROPS::InstatStokes2PhaseP2P1CL<DROPS::ZeroFlowCL>    MyStokesCL;

    const double L= 1; // Vol= 8*L*L*L;
    DROPS::Point3DCL orig(-L), e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 2*L;

    const int n= std::atoi( C.dmc_MeshFile.c_str());
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, n, n, n);

    const DROPS::BndCondT bc[6]=
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel};

    MyStokesCL prob(builder, DROPS::ZeroFlowCL(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun),
                    C.stk_XFEMStab<0 ? DROPS::P1_FE : DROPS::P1X_FE, C.stk_XFEMStab);

    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::EllipsoidCL::Init( C.exp_PosDrop, C.exp_RadDrop );
    DROPS::AdapTriangCL adap( mg, C.ref_Width, C.ref_CoarsestLevel, C.ref_FinestLevel);
    adap.MakeInitialTriang( DROPS::EllipsoidCL::DistanceFct);

    std::cout << DROPS::SanityMGOutCL(mg) << std::endl;

    Strategy( prob, adap);  // do all the stuff

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

