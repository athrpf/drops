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
#include <fstream>
#include <iomanip>
#include <vector>


DROPS::ParamMesszelleCL C;

enum StokesMethod { 
	minres                   =  1, // Minres without PC

        pminresmgssor            =  2, // MG-PC for A, SSOR for S
                                       // <SolverAsPreCL<MGSolverCL>, ISPreCL>

        pminresmgpcg             =  3, // MG-PC for A, PCG for S
                                       // <SolverAsPreCL<MGSolverCL>, ISBBT>

        pminrespcgssor           =  5, // PCG for A, SSOR for S
                                       // <SolverAsPreCL<PCGSolverCL<SSORPcCL>>, ISPreCL>

        pminrespcgpcg            =  6, // PCG for A, PCG for S
                                       // <SolverAsPreCL<PCGSolverCL<SSORPcCL>>, ISBBT>

        pminresmglumped          =  8, // MG-PC for A,    DiagPrMassMatrix for S

        pminrespcglumped         =  9, // PCG for A,      DiagPrMassMatrix for S

      //---------------------------------------------------------------------------------------------------

        inexactuzawamgssor       = 10, // MG-PC for A, SSOR for S
                                       // <SolverAsPreCL<MGSolverCL>, ISPreCL>

        inexactuzawamgpcg        = 11, // MG-PC for A, PCG for S
                                       // <SolverAsPreCL<MGSolverCL>, ISBBT>

        inexactuzawapcgssor      = 12, // PCG for A, SSOR for S
                                       // <SolverAsPreCL<PCGSolverCL<SSORPcCL>>, ISPreCL>

        inexactuzawapcgpcg       = 13, // PCG for A, PCG for S
                                       // <SolverAsPreCL<PCGSolverCL<SSORPcCL>>, ISBBT>

        inexactuzawamglumped     = 14, // MG-PC for A,    DiagPrMassMatrix for S

        inexactuzawapcglumped    = 15, // PCG for A,      DiagPrMassMatrix for S

      //---------------------------------------------------------------------------------------------------

        pcgmgssor                = 16, // MG-PC for A, SSOR for S
                                       // <SolverAsPreCL<MGSolverCL>, ISPreCL>

        pcgmgpcg                 = 17, // MG-PC for A, PCG for S
                                       // <SolverAsPreCL<MGSolverCL>, ISBBT>

        pcgmglumped              = 18  // MG-PC for A,    DiagPrMassMatrix for B
};

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

// Diagonal of a matrix as preconditioner
class DiagMatrixPCCL
{
  private:
    VectorCL& M_;
    bool sq_;

  public:
    DiagMatrixPCCL( VectorCL& M, bool sq = false)
        :M_( M), sq_(sq) {}

    template <typename Mat, typename Vec>
    void Apply(const Mat& , Vec& x, const Vec& b) const
    {
        for (Uint i= 0; i<M_.size(); ++i) {
            if (sq_) x[i]= b[i]/std::sqrt(M_[i]);
            else x[i]= b[i]/M_[i]; // M_ is a diagonal-matrix: exact inversion
        }
    }
};

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
        std::cerr << "InitPr not implemented yet for this FE type!\n";
    }
}

double my_abs( double x) { return std::abs(x); }

void L2ErrorPr( const VecDescCL& p, const LevelsetP2CL& lset, const MatrixCL& prM, double delta_p, const MultiGridCL& mg, const FiniteElementT prFE, const ExtIdxDescCL& Xidx, double p_ex_avg)
{
    const double min= p.Data.min(), max= p.Data.max();
    std::cerr << "pressure min/max/diff:\t" << min << "\t" << max << "\t" << (max-min-delta_p) << "\n";

    VectorCL ones( 1.0, p.Data.size());
    if (prFE==P1X_FE)
        for (int i=Xidx.GetNumUnknownsP1(), n=ones.size(); i<n; ++i)
            ones[i]= 0;
    const double Vol= dot( prM*ones, ones)*C.muD; // note that prM is scaled by 1/mu !!
// std::cerr << "Vol = " << Vol << '\n';
    const double p_avg= dot( prM*p.Data, ones)*C.muD/Vol; // note that prM is scaled by 1/mu !!
    VectorCL diff( p.Data - p_avg*ones);
    const double p0_avg= dot( prM*diff, ones)*C.muD/Vol;
    std::cerr << "average of pressure:\t" << p_avg << std::endl;
    std::cerr << "avg. of scaled pr:\t" << p0_avg << std::endl;

    if (prFE==P1X_FE)
    {
        VecDescCL p_exakt( p.RowIdx);
        InitPr( p_exakt, delta_p, mg, prFE, Xidx);
        const double p_ex_avg2= dot( prM*p_exakt.Data, ones)*C.muD/Vol;
        std::cerr << "avg. of exact pr:\t" << p_ex_avg2 << std::endl;
        diff-= VectorCL( p_exakt.Data - p_ex_avg*ones);
        const double L2= std::sqrt( C.muD*dot( prM*diff, diff));
        std::cerr << "*************\n"
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
    std::cerr << "*************\n"
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
    std::cerr << name << ":\t2-norm: "
        << norm( v) << "\tmax: " << supnorm( v) << std::endl;
}

template<class Coeff>
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    sigma= C.sigma;
    // Levelset-Disc.: Crank-Nicholson
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff);

//    lset.SetSurfaceForce( SF_LB);
    lset.SetSurfaceForce( SF_ImprovedLB);
//    lset.SetSurfaceForce( SF_Const);
    const double Vol= 8.,
//        prJump= C.sigma, // for SF_Const force
        prJump= C.sigma*2/C.Radius[0], // for SF_*LB force
        avg_ex= prJump/2.*(8./3.*M_PI*C.Radius[0]*C.Radius[0]*C.Radius[0] - Vol)/Vol; // for spherical interface
//        avg_ex= 0; // for planar interface

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    Stokes.SetNumVelLvl( MG.GetNumLevel());
    //Stokes.SetNumPrLvl ( MG.GetNumLevel());

    VecDescCL new_pr;  // for pressure output in Ensight

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr ( MG.GetLastLevel(), pidx, NULL, &lset);
    MG.SizeInfo( std::cerr);
    Stokes.SetIdx();
    Stokes.v.SetIdx(vidx);
    Stokes.p.SetIdx(pidx);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";

    new_pr.SetIdx( lidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset); // makes no sense for P0

    //LumpedMassMatrix von abs(M)
/*    VectorCL prMLDiag (0.0, Stokes.p.Data.size());
    for (Uint i=0; i< prMLDiag.size(); i++) {
        for (Uint j=0; j< prMLDiag.size(); j++) {
            prMLDiag[i]+=my_abs(prM.Data(i,j));
        }
    }
*/

    //LumpedMassMatrix = diag(M)
    VectorCL prMLDiag (0.0, Stokes.p.Data.size());
    prMLDiag=Stokes.prM.Data.GetDiag();
//    for (Uint i=0; i< prMLDiag.size(); i++) {
//	prMLDiag[i]=Stokes.prM.Data(i,i);
//    }

    //LumpedMassMatrix
/*    VectorCL prMLDiag(0.0, prM.Data.num_rows());
    VectorCL e( 1.0, prM.Data.num_rows());
    e[std::slice( Stokes.GetXidx().GetNumUnknownsP1(),
            prM.Data.num_rows() - Stokes.GetXidx().GetNumUnknownsP1(), 1)]= 0.0;
  //  prMLDiag = prM.Data * e;
    }*/

    SSORPcCL ssor;
    PCG_SsorCL PCG( ssor, C.inner_iter, C.inner_tol); // for computing curvature force

    VelVecDescCL curv( vidx);
    VelVecDescCL curvForce( vidx);

    switch (C.IniCond)
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
        std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        InitPr( Stokes.p, prJump, MG, Stokes.GetPrFE(), Stokes.GetXidx().GetFinest());
        VectorCL surf( Stokes.b.Data + curv.Data), BTp( transp_mul( Stokes.B.Data, Stokes.p.Data));
        PrintNorm( "surf. force", curv.Data);
        PrintNorm( "BT p", BTp);
        PrintNorm( "Diff.", VectorCL(curv.Data - BTp));
        std::cerr << "Solving velocity for exact pressure given...\n";
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
        std::cerr << "Discretizing Stokes/Surf.Force for initial velocities took "<<time.GetTime()<<" sec.\n";

        //InitPr( Stokes.p, prJump, MG, Stokes.GetPrFE(), Stokes.GetXidx().GetFinest());
        time.Reset();

        const double kA=0;
        const double kM=1;
        // Preconditioner for A
            //Multigrid
        Stokes.SetupProlongations();
        CheckMGData( Stokes.A.Data, Stokes.PPr.Data);
        SSORsmoothCL smoother(1.0);
        PCG_SsorCL   coarsesolver( SSORPcCL(1.0), 500, C.inner_tol);
        MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc ( Stokes.PVel.Data, smoother, coarsesolver, 1, -1.0, false);
        typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > MGPCT;
        MGPCT MGPC (mgc);
        VectorCL xx( 1.0, vidx->NumUnknowns());
        double rhoinv = 0.99*( 1.0-1.1*0.363294);
/*        if ( C.StokesMethod == 16 || C.StokesMethod == 17 || C.StokesMethod == 18)
            rhoinv= 0.99*( 1.0 - 1.1*EigenValueMaxMG( Stokes.A.Data, Stokes.PVel.Data, xx, 1000, 1e-4));*/
        ScaledMGPreCL velprep( Stokes.PVel.Data, 1, 1.0/rhoinv);

            //PCG
        typedef SSORPcCL APcPcT;
        APcPcT Apcpc;
        typedef PCGSolverCL<APcPcT> ASolverT;        // CG-based APcT
        ASolverT Asolver( Apcpc, 500, 0.02, true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver);

        // Preconditioner for instat. Schur complement
        typedef ISBBTPreCL ISBBT;
        ISBBT isbbt (&Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), kA, kM);
        ISPreCL ispc( Stokes.prA.Data, Stokes.prM.Data, kA, kM);
        DiagMatrixPCCL lumped( prMLDiag);

        typedef BlockPreCL<APcT, DiagMatrixPCCL> LanczosPcT;

        // Preconditioner for PMINRES
        typedef BlockPreCL<MGPCT, ISPreCL> Lanczos2PCT;
        typedef PLanczosONBCL<VectorCL, Lanczos2PCT> Lanczos2T;
        typedef BlockPreCL<MGPCT, ISBBT> Lanczos3PCT;
        typedef PLanczosONBCL<VectorCL, Lanczos3PCT> Lanczos3T;
        typedef BlockPreCL<APcT, ISPreCL> Lanczos5PCT;
        typedef PLanczosONBCL< VectorCL, Lanczos5PCT> Lanczos5T;
        typedef BlockPreCL<APcT, ISBBT> Lanczos6PCT;
        typedef PLanczosONBCL<VectorCL, Lanczos6PCT> Lanczos6T;
        typedef BlockPreCL<MGPCT,DiagMatrixPCCL> Lanczos8PCT;
        typedef PLanczosONBCL<VectorCL, Lanczos8PCT> Lanczos8T;
        typedef BlockPreCL<APcT,DiagMatrixPCCL> Lanczos9PCT;
        typedef PLanczosONBCL<VectorCL, Lanczos9PCT> Lanczos9T;

        Lanczos2PCT lanczos2pc (MGPC, ispc);
        Lanczos2T lanczos2 (lanczos2pc);

        Lanczos3PCT lanczos3pc (MGPC, isbbt);
        Lanczos3T lanczos3 (lanczos3pc);

        Lanczos5PCT lanczos5pc (Apc, ispc);
        Lanczos5T lanczos5 (lanczos5pc);

        Lanczos6PCT lanczos6pc (Apc, isbbt);
        Lanczos6T lanczos6 (lanczos6pc);

        Lanczos8PCT lanczos8pc (MGPC, lumped);
        Lanczos8T lanczos8 (lanczos8pc);

        Lanczos9PCT lanczos9pc (Apc, lumped);
        Lanczos9T lanczos9 (lanczos9pc);

        // available Stokes Solver
        typedef MResSolverCL PMinres1T; // Minres
        PMinres1T minressolver        (C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres1T> blockminressolver(minressolver);

        typedef PMResSolverCL<Lanczos2T> PMinres2T; // PMinRes - MG-ISPreCL
        PMinres2T pminresmgssorsolver (lanczos2, C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres2T> blockpminresmgssorsolver(pminresmgssorsolver);

        typedef PMResSolverCL<Lanczos3T> PMinres3T; // PMinRes - MG-ISBBTCL
        PMinres3T pminresmgpcgsolver  (lanczos3, C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres3T> blockpminresmgpcgsolver(pminresmgpcgsolver);

        typedef PMResSolverCL<Lanczos5T> PMinres5T; // PMinRes - PCG-ISPreCL
        PMinres5T pminrespcgssorsolver (lanczos5, C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres5T> blockpminrespcgssorsolver(pminrespcgssorsolver);

        typedef PMResSolverCL<Lanczos6T> PMinres6T; // PMinRes - PCG-ISBBT
        PMinres6T pminrespcgpcgsolver (lanczos6, C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres6T> blockpminrespcgpcgsolver(pminrespcgpcgsolver);

        typedef PMResSolverCL<Lanczos8T> PMinres8T; 
        PMinres8T pminresmgdiagsolver (lanczos8, C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres8T> blockpminresmgdiagsolver(pminresmgdiagsolver);

        typedef PMResSolverCL<Lanczos9T> PMinres9T; 
        PMinres9T pminrespcgdiagsolver (lanczos9, C.outer_iter, C.outer_tol);
        BlockMatrixSolverCL<PMinres9T> blockpminrespcgdiagsolver(pminrespcgdiagsolver);

        typedef InexactUzawaCL<MGPCT, ISPreCL, APC_SYM> InexactUzawa10T;
        InexactUzawa10T inexactuzawamgssorsolver( MGPC, ispc, C.outer_iter, C.outer_tol, 0.5);

        typedef InexactUzawaCL<MGPCT, ISBBT, APC_SYM> InexactUzawa11T;
        InexactUzawa11T inexactuzawamgpcgsolver( MGPC, isbbt, C.outer_iter, C.outer_tol, 0.5);

        typedef InexactUzawaCL<APcT, ISPreCL, APC_SYM> InexactUzawa12T;
        InexactUzawa12T inexactuzawapcgssorsolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.5);

        typedef InexactUzawaCL<APcT, ISBBT, APC_SYM> InexactUzawa13T;
        InexactUzawa13T inexactuzawapcgpcgsolver( Apc, isbbt, C.outer_iter, C.outer_tol, 0.5);

        typedef InexactUzawaCL<MGPCT, DiagMatrixPCCL, APC_SYM> InexactUzawa14T;
        InexactUzawa14T inexactuzawamgdiagsolver( MGPC, lumped, C.outer_iter, C.outer_tol, 0.5);

        typedef InexactUzawaCL<APcT, DiagMatrixPCCL, APC_SYM> InexactUzawa15T;
        InexactUzawa15T inexactuzawapcgdiagsolver( Apc, lumped, C.outer_iter, C.outer_tol, 0.5);

        typedef UzawaCGSolverEffCL<ScaledMGPreCL, ISPreCL> UzawaCGEff16T;
        UzawaCGEff16T pcgmgssorsolver (velprep, ispc, C.outer_iter, C.outer_tol);

        typedef UzawaCGSolverEffCL<ScaledMGPreCL, ISBBT> UzawaCGEff17T;
        UzawaCGEff17T pcgmgpcgsolver (velprep, isbbt, C.outer_iter, C.outer_tol);

        typedef UzawaCGSolverEffCL<ScaledMGPreCL, DiagMatrixPCCL> UzawaCGEff18T;
        UzawaCGEff18T pcgmgdiagsolver (velprep, lumped, C.outer_iter, C.outer_tol);

        StokesSolverBaseCL* solver=0;
        switch (C.StokesMethod) {
            case minres:                solver = &blockminressolver;         break;
            case pminresmgssor:         solver = &blockpminresmgssorsolver;  break;
            case pminresmgpcg:          solver = &blockpminresmgpcgsolver;   break;
            case pminrespcgssor:        solver = &blockpminrespcgssorsolver; break;
            case pminrespcgpcg:         solver = &blockpminrespcgpcgsolver;  break;
            case pminresmglumped:       solver = &blockpminresmgdiagsolver;  break;
            case pminrespcglumped:      solver = &blockpminrespcgdiagsolver; break;
            case inexactuzawamgssor:    solver = &inexactuzawamgssorsolver;  break;
            case inexactuzawamgpcg:     solver = &inexactuzawamgpcgsolver;   break;
            case inexactuzawapcgssor:   solver = &inexactuzawapcgssorsolver; break;
            case inexactuzawapcgpcg:    solver = &inexactuzawapcgpcgsolver;  break;
            case inexactuzawamglumped:  solver = &inexactuzawamgdiagsolver;  break;
            case inexactuzawapcglumped: solver = &inexactuzawapcgdiagsolver; break;
            case pcgmgssor:             solver = &pcgmgssorsolver;           break;
            case pcgmgpcg:              solver = &pcgmgpcgsolver;            break;
            case pcgmglumped:           solver = &pcgmgdiagsolver;           break;

            default: throw DROPSErrCL( "unknown method\n");
        }
        if (solver != 0) {
            solver->Solve( Stokes.A.Data, Stokes.B.Data, Stokes.v.Data, Stokes.p.Data,
                           curv.Data, Stokes.c.Data);
            std::cerr << "iter: " << solver->GetIter()
                      << "\tresid: " << solver->GetResid() << std::endl;
        }
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
      }
    }

    const VectorCL& u= Stokes.v.Data;
    std::cerr << "\n----------------\n || u ||_oo = " << supnorm(u)
              << "\n || u ||_M  = " << std::sqrt( dot( Stokes.M.Data*u, u))
              << "\n || u ||_A  = " << std::sqrt( dot( Stokes.A.Data*u, u))
              << "\n----------------\n";
    if (Stokes.GetPrFE()==P1X_FE)
    {
        const MLExtIdxDescCL& Xidx= Stokes.GetXidx();
        const size_t n= Stokes.p.Data.size();

        const double limtol= 10,
            lim_min= -prJump - limtol*prJump,
            lim_max= -prJump + limtol*prJump; 
        double xmin= 1e99, xmax= -1e99, sum= 0, sum_lim= 0;
        IdxT num= 0;
        for (size_t i=Xidx.GetNumUnknownsP1(); i<n; ++i)
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
        std::cerr << "extended pr: min/max/avg = " << xmin << ", " << xmax << ", " << sum/num << std::endl;
        std::cerr << "limited pr:  min/max/avg = " << lim_min << ", " << lim_max << ", " << sum_lim/num << std::endl;
    }

    L2ErrorPr( Stokes.p, lset, Stokes.prM.Data.GetFinest(), prJump, MG, Stokes.GetPrFE(), Stokes.GetXidx().GetFinest(), avg_ex);

    PostProcessPr( Stokes.p, new_pr, MG);

    // Initialize Ensight6 output
    std::string ensf( C.EnsDir + "/" + C.EnsCase);
    Ensight6OutCL ensight( C.EnsCase + ".case", C.num_steps + 1);
    ensight.Register( make_Ensight6Geom  ( MG, MG.GetLastLevel(),          "Cube",      ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),             "Levelset",  ensf + ".scl"));
    ensight.Register( make_Ensight6Scalar( Stokes.GetPrSolution( new_pr),  "Pressure",  ensf + ".pr"));
    ensight.Register( make_Ensight6Vector( Stokes.GetVelSolution(),        "Velocity",  ensf + ".vel"));
    ensight.Register( make_Ensight6Vector(  Stokes.GetVelSolution( curvForce),
                                                                           "Curvature", ensf + ".crv"));
    if (Stokes.UsesXFEM())
        ensight.Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.GetXidx().GetFinest(), Stokes.p.Data,
                                                                           "XPressure", ensf + ".pr"));

    if (C.EnsCase != "none") ensight.Write();

    std::cerr << std::endl;
}

} // end of namespace DROPS


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
        param.open( "prJump.param");
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    typedef DROPS::InstatStokes2PhaseP2P1CL<ZeroFlowCL>    MyStokesCL;

    const double L= 1; // Vol= 8*L*L*L;
    DROPS::Point3DCL orig(-L), e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 2*L;

    const int n= std::atoi( C.meshfile.c_str());
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, n, n, n);

    const DROPS::BndCondT bc[6]=
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel};

    DROPS::FiniteElementT prFE=DROPS::P1_FE;
    if (C.XFEMStab>=0) prFE=DROPS::P1X_FE;

    MyStokesCL prob(builder, ZeroFlowCL(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun),
        prFE, C.XFEMStab);

    DROPS::MultiGridCL& mg = prob.GetMG();
    EllipsoidCL::Init( C.Mitte, C.Radius );
    DROPS::AdapTriangCL adap( mg, C.ref_width, 0, C.ref_flevel);
    adap.MakeInitialTriang( EllipsoidCL::DistanceFct);

    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    Strategy( prob);  // do all the stuff

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
