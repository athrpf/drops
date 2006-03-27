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

// program for testing various FE pressure spaces 
// using a special constant surface force: 
//     \sigma \int_\Gamma v n ds
// => implies a constant pressure jump [p] = \sigma across \Gamma.

// rho*du/dt - mu*laplace u + Dp = f + rho*g - ovn
//                        -div u = 0
//                             u = u0, t=t0


class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleCL& c) 
      : rho( DROPS::JumpCL( c.rhoD, c.rhoF ), DROPS::H_sm, c.sm_eps),
         mu( DROPS::JumpCL( c.muD,  c.muF),   DROPS::H_sm, c.sm_eps),
        g( c.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const double L, U;
    const DROPS::SmoothedJumpCL rho, mu;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( double l, double u, const DROPS::ParamMesszelleCL& c) 
      : L( l), U( u),
        rho( DROPS::JumpCL( 1, c.rhoF/c.rhoD ), DROPS::H_sm, c.sm_eps/L),
         mu( DROPS::JumpCL( 1, c.muF/c.muD),   DROPS::H_sm, c.sm_eps/L),
        g( (L/U/U)*c.g)    {}
};


DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }
/*
double DistanceFct( const DROPS::Point3DCL& p)
{ // ball
    const DROPS::Point3DCL d= C.Mitte-p;
    return d.norm()-C.Radius;
}

double DistanceFct( const DROPS::Point3DCL& p)
{ // cube
    const DROPS::Point3DCL d= C.Mitte-p;
    double max=-1;
    for (int i=0; i<3; ++i)
        if (std::fabs(d[i])>max) max= std::fabs(d[i]);
    return max-C.Radius;
}
*/
double DistanceFct( const DROPS::Point3DCL& p)
{ // plane at z=z_0. z_0 determined by z-component of PosDrop.
    const DROPS::Point3DCL d= C.Mitte-p;
    return d[2]; // z-Plane
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


void InitPr( VecDescCL& p, double delta_p, const MultiGridCL& mg)
{
    const IdxDescCL& idx= *p.RowIdx;
    const Uint lvl= p.RowIdx->TriangLevel,
        idxnum= p.RowIdx->GetIdx();
    
    delta_p/= 2;
    if (idx.NumUnknownsVertex && idx.NumUnknownsTetra) // P01_FE
        for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl), 
            end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
        {
            p.Data[it->Unknowns(idxnum)]= 0;
        }    
    if (idx.NumUnknownsVertex && !idx.NumUnknownsTetra) // P1_FE
        for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl), 
            end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
        {
            const double dist= DistanceFct( it->GetCoord());
            p.Data[it->Unknowns(idxnum)]= dist > 0 ? -delta_p : delta_p;
        }    

    if (idx.NumUnknownsTetra)
        for( MultiGridCL::const_TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl), 
            end= mg.GetTriangTetraEnd(lvl); it != end; ++it)
        {
            const double dist= DistanceFct( GetBaryCenter( *it));
            p.Data[it->Unknowns(idxnum)]= dist > 0 ? -delta_p : delta_p;
        }    
}

double my_abs( double x) { return std::abs(x); }

void ErrorPr( const VecDescCL& p, const LevelsetP2CL& lset, double delta_p, const MultiGridCL& mg)
{
    if (p.RowIdx->NumUnknownsTetra) return;
    
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint lvl= p.GetLevel();
    IdxT prNumb[4];
    double L1Err= 0, L2Err= 0;
    JumpCL PressureJump( delta_p/2, -delta_p/2);
    LocalP2CL<double> locallset;
    Quad2CL<double> diff, exact;
    
    // ToDo: compute avg. of p for scaling
    double vol= 0, integr= 0, integr_ex= 0;
    for (MultiGridCL::const_TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(lvl),
      send= mg.GetTriangTetraEnd(lvl); sit != send; ++sit) 
    {
        vol+= sit->GetVolume();
        const double absdet= sit->GetVolume()*6.;
        GetLocalNumbP1NoBnd( prNumb, *sit, *p.RowIdx);
        double sum= 0;
        for (int i=0; i<4; ++i)
            sum+= diff[i]= p.Data[prNumb[i]];
        diff[4]= sum/4; // value in barycenter

        const double pr_ex= DistanceFct(GetBaryCenter( *sit)) > 0 ? -delta_p/2 : delta_p/2;
        for (int i= 0; i<5; ++i) exact[i]= pr_ex;

        integr+= diff.quad( absdet);
        integr_ex+= exact.quad( absdet);
    }
    const double avg= integr/vol, avg_ex= integr_ex/vol;
    std::cerr << "\navg pr = " << avg << "\tavg. exact pr = " << avg_ex << "\n";    

    for (MultiGridCL::const_TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(lvl),
      send= mg.GetTriangTetraEnd(lvl); sit != send; ++sit) {
        
        const double absdet= sit->GetVolume()*6.;

        const double pr_ex= DistanceFct(GetBaryCenter( *sit)) > 0 ? -delta_p/2 : delta_p/2;
        for (int i= 0; i<5; ++i) exact[i]= pr_ex - avg_ex;

        GetLocalNumbP1NoBnd( prNumb, *sit, *p.RowIdx);
        double sum= 0;
        for (int i=0; i<4; ++i)
            sum+= diff[i]= p.Data[prNumb[i]] - avg;
        diff[4]= sum/4; // value in barycenter
        


        diff-= exact;
//std::cerr << diff[0] << "\t" << diff[4] << "\n";        
        L2Err+= Quad2CL<>(diff*diff).quad( absdet);

        diff.apply( my_abs);
        L1Err+= diff.quad( absdet);
    }
    L2Err= std::sqrt( L2Err);
    std::cerr << "------- pressure --------\nL2 error:\t" << L2Err
              << "\nL1 error:\t" << L1Err << "\n---------------------------\n";
}

void ErrorPr( const VecDescCL& p, const MatrixCL& prM, double delta_p, const MultiGridCL& mg)
{
    const double min= p.Data.min(), max= p.Data.max();
    std::cerr << "pressure min/max/diff:\t" << min << "\t" << max << "\t" << (max-min-delta_p) << "\n";
    
    VectorCL ones( 1.0, p.Data.size());
    const double Vol= (prM*ones).sum()*C.muD; // note that prM is scaled by 1/mu !!
// std::cerr << "Vol = " << Vol << '\n';
    VecDescCL p_exakt( p.RowIdx);
    InitPr( p_exakt, delta_p, mg);
    const double p_avg= (prM*p.Data).sum()*C.muD/Vol; // note that prM is scaled by 1/mu !!
    VectorCL diff( p.Data - p_avg);
//    const double p0_avg= (prM*diff).sum()*C.muD/Vol;
    const double p_ex_avg= (prM*p_exakt.Data).sum()*C.muD/Vol;
// std::cerr << "average of pressure:\t" << p_avg << std::endl;
// std::cerr << "avg. of scaled pr:\t" << p0_avg << std::endl;
// std::cerr << "avg. of exact pr:\t" << p_ex_avg << std::endl;
    diff-= VectorCL( p_exakt.Data - p_ex_avg);
    const double err_sup= supnorm(diff), errM= std::sqrt( dot( prM*diff, diff));
    std::cerr << "max. pressure error:\t|| p ||_oo =  " << err_sup
              << "\napprox. L2-error:\t|| p ||_M = " << errM << "\n\n";
}

void PostProcessPr( const VecDescCL& p, VecDescCL& new_p, const MultiGridCL& mg)
// as Ensight output is for cont. data only, put values of P0 DoFs in tetras 
// into P1 DoFs in vertices (s.t. the average over each tetra gives
// the former P0 value)
{
    VectorCL num( new_p.Data.size()), sum( new_p.Data.size());
    
    const Uint lvl= p.RowIdx->TriangLevel,
        idxnum= p.RowIdx->GetIdx(),
        idxnum2= new_p.RowIdx->GetIdx();
    
    if (!p.RowIdx->NumUnknownsTetra) { new_p= p; return; }
    
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

    if (new_p.RowIdx->NumUnknownsVertex)
        for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl), 
            end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
        {
            const IdxT nr= it->Unknowns(idxnum),
                nr2= it->Unknowns(idxnum2);
            new_p.Data[nr2]= p.Data[nr] + sum[nr2]/num[nr2];
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
    // Levelset-Disc.: Crank-Nicholson
    LevelsetP2CL lset( MG, C.sigma, C.theta, C.lset_SD, C.RepDiff, C.lset_iter, C.lset_tol, C.CurvDiff); 
    lset.SetSurfaceForce( SF_Const);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    MatDescCL prM, prA;
    VecDescCL new_pr;  // for pressure output in Ensight

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
    new_pr.SetIdx( lidx);
    Stokes.InitVel( &Stokes.v, Null);
    Stokes.SetupPrMass(  &prM, lset);
//    Stokes.SetupPrStiff( &prA, lset); // makes no sense for P0

    PSchur_PCG_CL   schurSolver( prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
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

        InitPr( Stokes.p, C.sigma, MG);
        VectorCL surf( curv.Data), BTp( transp_mul( Stokes.B.Data, Stokes.p.Data));
        PrintNorm( "surf", surf);
        PrintNorm( "BT p", BTp);
        PrintNorm( "Diff.", VectorCL(surf-BTp));
        curvForce.Data= surf - BTp;
        std::cerr << "Solving velocity for exact pressure given...\n";
        PCG.Solve( Stokes.A.Data, Stokes.v.Data, VectorCL( Stokes.b.Data + curv.Data - transp_mul( Stokes.B.Data, Stokes.p.Data)) );
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
        std::cerr << "Discretizing Stokes/Curv for initial velocities took "<<time.GetTime()<<" sec.\n";

        DummyPcCL dpc;
        ISPreCL ispc( prA.Data, prM.Data, 0);
        SSORPCG_PreCL pcg( C.inner_iter, 0.2);
    //    typedef InexactUzawaCL<SSORPCG_PreCL, ISPreCL, APC_SYM> InexactUzawaT;
        typedef InexactUzawaCL<SSORPCG_PreCL, DummyPcCL, APC_SYM> InexactUzawaT;
        InexactUzawaT inexactUzawaSolver( pcg, dpc, C.outer_iter, C.outer_tol);

        time.Reset();
        schurSolver.Solve( Stokes.A.Data, Stokes.B.Data, 
//        inexactUzawaSolver.Solve( Stokes.A.Data, Stokes.B.Data, 
            Stokes.v.Data, Stokes.p.Data, VectorCL( Stokes.b.Data + curv.Data), Stokes.c.Data);
        time.Stop();
        std::cerr << "Solving Stokes for initial velocities took "<<time.GetTime()<<" sec.\n";
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        PCG.Solve( Stokes.M.Data, curvForce.Data, curv.Data);

      }
    }

    const VectorCL& u= Stokes.v.Data;       
    std::cerr << "\n----------------\n || u ||_oo = " << supnorm(u)
              << "\n || u ||_M  = " << std::sqrt( dot( Stokes.M.Data*u, u))
              << "\n || u ||_A  = " << std::sqrt( dot( Stokes.A.Data*u, u))
              << "\n----------------\n";

    ErrorPr( Stokes.p, prM.Data, C.sigma, MG);
    ErrorPr( Stokes.p, lset, C.sigma, MG);
    
    PostProcessPr( Stokes.p, new_pr, MG);

    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putVector( datcrv, Stokes.GetVelSolution( curvForce), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution( new_pr), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    DummyPcCL dpc;
    ISPreCL ispc( prA.Data, prM.Data, C.theta*C.dt);
    ISPSchur_PCG_CL ISPschurSolver( ispc,  C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
    ISPschurSolver.SetTol( C.outer_tol);
    
    SSORPCG_PreCL pcg( C.inner_iter, 0.2);
    typedef InexactUzawaCL<SSORPCG_PreCL, ISPreCL, APC_SYM> InexactUzawaT;
//    typedef InexactUzawaCL<SSORPCG_PreCL, DummyPcCL, APC_SYM> InexactUzawaT;
    InexactUzawaT inexactUzawaSolver( pcg, ispc, C.outer_iter, C.outer_tol);
//    InexactUzawaT inexactUzawaSolver( pcg, dpc, C.outer_iter, C.outer_tol);

//    CouplLevelsetStokes2PhaseCL<StokesProblemT, ISPSchur_PCG_CL> 
//        cpl( Stokes, lset, ISPschurSolver, C.theta);
    CouplLevelsetStokes2PhaseCL<StokesProblemT, InexactUzawaT> 
        cpl( Stokes, lset, inexactUzawaSolver, C.theta);
    cpl.SetTimeStep( C.dt);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( C.cpl_iter);
        const VectorCL& u= Stokes.v.Data;       
        std::cerr << "\n----------------\n || u ||_oo = " << supnorm(u)
                  << "\n || u ||_M  = " << std::sqrt( dot( Stokes.M.Data*u, u))
                  << "\n || u ||_A  = " << std::sqrt( dot( Stokes.A.Data*u, u))
                  << "\n----------------\n";
        curv.Clear();
        lset.AccumulateBndIntegral( curv);
        PCG.Solve( Stokes.A.Data, curvForce.Data, curv.Data);

        ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
        ensight.putVector( datcrv, Stokes.GetVelSolution( curvForce), step*C.dt);
        ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
        ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
        ensight.Commit();
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
/*            for (int i=0; i<4; ++i)
            {
                const double val= DistanceFct( It->GetVertex(i)->GetCoord() );
                if ( val<C.ref_width && val > -C.ref_width)
                    It->SetRegRefMark();   
            }
            const double val= DistanceFct( GetBaryCenter(*It));    
            if ( val<C.ref_width && val > -C.ref_width)
                It->SetRegRefMark();   
*/
            int neg= 0, zero= 0;
            for (int i=0; i<4; ++i)
            {
                const double val= DistanceFct( It->GetVertex(i)->GetCoord() );
                neg+= val<0; 
                zero+= std::fabs(val)<1e-9;
            }   
            const double val= DistanceFct( GetBaryCenter(*It));
            neg+= val<0; 
            zero+= std::fabs(val)<1e-9;
            
            if ( (neg!=0 && neg!=5) || zero) // change of sign or zero in tetra
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
    
    const double L= 2e-3; // Vol= 8*L*L*L;
    DROPS::Point3DCL orig(-L), e1, e2, e3;
    e1[0]= e2[1]= e3[2]= 2*L;
    
    const int n= std::atoi( C.meshfile.c_str());
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, n, n, n);
    
    const DROPS::BndCondT bc[6]= 
        { DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC, DROPS::WallBC};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Null}; 
        
    MyStokesCL prob(builder, ZeroFlowCL(C), DROPS::StokesBndDataCL( 6, bc, bnd_fun), 
        DROPS::P0_FE);
//               ^--------- FE type for pressure space
        
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
