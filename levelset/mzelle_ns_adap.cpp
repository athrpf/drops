//**************************************************************************
// File:    mzelle_instat.cpp                                              *
// Content: flow in drop cell                                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include <fstream>
#include <sstream>


DROPS::ParamMesszelleNsCL C;

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

namespace DROPS // for Strategy
{

template <class It>
  std::pair<double, double>
  h_interface (It begin, It end, const VecDescCL& ls)
{
    double hmin= 1e99, hmax= -1.0; //, hmean= 0.;
    // size_t num= 0;
    // EdgeCL* emax= 0;

    for (; begin != end; ++begin)
        if (   (ls.Data[begin->GetVertex( 0)->Unknowns( ls.RowIdx->GetIdx())]
                *ls.Data[begin->Unknowns( ls.RowIdx->GetIdx())] <= 0.)
            || (ls.Data[begin->GetVertex( 1)->Unknowns( ls.RowIdx->GetIdx())]
                *ls.Data[begin->Unknowns( ls.RowIdx->GetIdx())] <= 0.)) {
            const double h= (begin->GetVertex( 0)->GetCoord() - begin->GetVertex( 1)->GetCoord()).norm();
            hmin= std::min( hmin, h);
            // if (h > hmax) emax= &*begin;
            hmax= std::max( hmax, h);
            // hmean+= h;
            // ++num;
        }
    // std::cerr << "mean : " << hmean/num << '\n';
    // emax->DebugInfo( std::cerr);
    // emax->GetVertex( 0)->DebugInfo( std::cerr);
    // emax->GetVertex( 1)->DebugInfo( std::cerr);
    
    return std::make_pair( hmin, hmax);
}

template<class Coeff>
void WriteMatrices (InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, int i)
{
    std::string pfad( "matrices/");
    std::ostringstream n;
    n << pfad << "A" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    std::ofstream s( n.str().c_str());
    s << std::setprecision( 15);
    s << Stokes.A.Data;
    s.close();

    n.str( "");
    n << pfad << "B" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.B.Data;
    s.close();

    n.str( "");
    n << pfad << "M" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.M.Data;
    s.close();

    n.str( "");
    n << pfad << "prA" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.prA.Data;
    s.close();

    n.str( "");
    n << pfad << "prM" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.prM.Data;
    s.close();

    n.str( "");
    n << pfad << "N" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.N.Data;
    s.close();
}

template<class Coeff>
void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adap)
// flow control
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    sigma= Stokes.GetCoeff().SurfTens;
    LevelsetP2CL lset( MG, &sigmaf, &gsigma, C.lset_theta, C.lset_SD,
        -1, C.lset_iter, C.lset_tol, C.CurvDiff);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    VecDescCL cplN;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);
    lset.SetSurfaceForce( SF_ImprovedLB);
//    lset.SetSurfaceForce( SF_ImprovedLBVar);
//    IFInfo.Update( lset);
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);

    EnsightP2SolOutCL ensight( MG, lidx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datpr = filename+".pr" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo, true);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Pressure", datpr,  true);
    ensight.DescribeVector( "Velocity", datvec, true);

    MG.SizeInfo( std::cerr);
    Stokes.b.SetIdx( vidx);
    Stokes.v.SetIdx( vidx);
    cplN.SetIdx( vidx);
    Stokes.c.SetIdx( pidx);
    Stokes.p.SetIdx( pidx);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";
    Stokes.A.SetIdx(vidx, vidx);
    Stokes.B.SetIdx(pidx, vidx);
    Stokes.prM.SetIdx( pidx, pidx);
    Stokes.prA.SetIdx( pidx, pidx);
    Stokes.M.SetIdx(vidx, vidx);
    Stokes.N.SetIdx(vidx, vidx);

    Stokes.InitVel( &Stokes.v, ZeroVel);
    Stokes.SetupPrMass(  &Stokes.prM, lset);
    Stokes.SetupPrStiff( &Stokes.prA, lset);

    switch (C.IniCond)
    {
      case 1: case 2: // stationary flow with/without drop
      {
//           typedef CGSolverCL SPcSolverT;
//           SPcSolverT CGsolver( 50, 0.02, /*relative*/ true);
          typedef PCGSolverCL<JACPcCL> SPcSolverT;
          SPcSolverT CGsolver( JACPcCL(), 50, 0.02, /*relative*/ true);
          typedef ISNonlinearPreCL<SPcSolverT> SPcT;
          SPcT ispc( CGsolver, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 0., /*kM*/ 1.);
//          typedef ISPreCL SPcT;
//          SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 0., /*kM*/ 1.);

        // PC for A-Block-PC
//          typedef  DummyPcCL APcPcT;
          typedef JACPcCL  APcPcT;
//          typedef SSORPcCL APcPcT;
//          typedef GSPcCL   APcPcT;
          APcPcT Apcpc;

        // PC for A-block
          typedef PCGSolverCL<APcPcT> ASolverT;        // CG-based APcT
          ASolverT Asolver( Apcpc, 500, 0.02, /*relative=*/ true);
//          typedef GMResSolverCL<APcPcT>    ASolverT;        // GMRes-based APcT
//          ASolverT Asolver( Apcpc, 500, /*restart*/ 100, 0.02, /*relative=*/ true);
          typedef SolverAsPreCL<ASolverT> APcT;
          APcT Apc( Asolver);

        // PC for Oseen solver
//        typedef DummyPcCL OseenPcT;
//        OseenPcT oseenpc;
//        typedef DiagBlockPreCL<APcT, SPcT> OseenPcT;
//        OseenPcT oseenpc( Apc, ispc);

        // Oseen solver
          typedef InexactUzawaCL<APcT, SPcT, APC_SYM> OseenSolverT;
          OseenSolverT schurSolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.1);

//        PSchur_PCG_CL schurSolver( Stokes.prM.Data, C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);
//        Schur_GMRes_CL schurSolver( C.outer_iter, C.outer_tol, C.inner_iter, C.inner_tol);

        if (C.IniCond==2) // stationary flow without drop
            lset.Init( &One);

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
        std::cerr << "Solving Stokes for initial velocities took "<< time.GetTime()<<" sec.\n";

        if (C.IniCond==2)
            lset.Init( EllipsoidCL::DistanceFct);
      } break;

      case 3: // read from file
      { // FIXME: only the P1-part of P1X-elements can be read
        ReadEnsightP2SolCL reader( MG);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset, /*NumberingChanged*/ false);
        Stokes.p.SetIdx( pidx); // Zero-vector for now.
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr); // reads the P1-part of the pressure
      } break;

      default:
        lset.Init( EllipsoidCL::DistanceFct);
    }
    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

//     std::pair<double, double> h_int= h_interface( MG.GetTriangEdgeBegin(), MG.GetTriangEdgeEnd(), lset.Phi);
//     std::cerr << "h interface: " << h_int.first << '\t' << h_int.second << '\n';

    ensight.putGeom( datgeo, 0);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    if (C.scheme)
    {
        // PC for instat. Schur complement
//        typedef DummyPcCL SPcT;
//        SPcT ispc;
//         typedef ISPreCL SPcT;
//         SPcT ispc( Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));
        typedef ISBBTPreCL SPcT;
        SPcT ispc( Stokes.B.Data, Stokes.prM.Data, Stokes.M.Data,
            /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));
//        typedef PCGSolverCL<JACPcCL> SPcSolverT;
//        SPcSolverT JacCG( JACPcCL(), 50, 0.02, /*relative*/ true);
//        typedef ISNonlinearPreCL<SPcSolverT> SPcT;
//        SPcT ispc( JacCG, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));
//        typedef CGSolverCL SPcSolverT;
//        SPcSolverT CGsolver( 50, 0.02, /*relative*/ true);
//        typedef ISNonlinearPreCL<SPcSolverT> SPcT;
//        SPcT ispc( CGsolver, Stokes.prA.Data, Stokes.prM.Data, /*kA*/ 1./C.dt, /*kM=alpha*theta in FracStep*/ 3 - 2*std::sqrt(2.));

        // PC for A-Block-PC
//      typedef  DummyPcCL APcPcT;
        typedef JACPcCL  APcPcT;
//        typedef SSORPcCL APcPcT;
//        typedef GSPcCL   APcPcT;
        APcPcT Apcpc;

        // PC for A-block
//        typedef BiCGStabSolverCL<APcPcT> ASolverT;        // BiCGStab-based APcT
//        ASolverT Asolver( Apcpc, 500, 0.02, /*relative=*/ true);
        typedef GMResSolverCL<APcPcT>    ASolverT;        // GMRes-based APcT
        ASolverT Asolver( Apcpc, 500, /*restart*/ 100, 1e-5, /*relative=*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver);

        // PC for Oseen solver
//        typedef DummyPcCL OseenPcT;
//        OseenPcT oseenpc;
        typedef DiagBlockPreCL<APcT, SPcT> OseenPcT;
        OseenPcT oseenpc( Apc, ispc);

        // Oseen solver
//        typedef InexactUzawaCL<APcT, SPcT, APC_OTHER> OseenSolverT;
//        OseenSolverT oseensolver( Apc, ispc, C.outer_iter, C.outer_tol, 0.2);
       typedef GCRSolverCL<OseenPcT> OseenBaseSolverT;
       OseenBaseSolverT oseensolver0( oseenpc, /*truncate*/ C.outer_iter, C.outer_iter, C.outer_tol, /*relative*/ false);
       typedef BlockMatrixSolverCL<OseenBaseSolverT> OseenSolverT;
       OseenSolverT oseensolver( oseensolver0);

        // Navier-Stokes solver
        typedef AdaptFixedPtDefectCorrCL<StokesProblemT, OseenSolverT> NSSolverT;
        NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red);

        // Time-integration and coupling
//        typedef CouplLevelsetNavStokes2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//        CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear, C.cpl_stab);
        typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
        CouplingT cpl( Stokes, lset, nssolver, C.nonlinear, C.cpl_stab);

        cpl.SetTimeStep( C.dt);

        for (int step= 1; step<=C.num_steps; ++step)
        {
            std::cerr << "======================================================== Schritt " << step << ":\n";
//            IFInfo.Update( lset);
            cpl.DoStep( C.cpl_iter);
//            WriteMatrices( Stokes, step);
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr)
            {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }

            if (C.RepFreq && step%C.RepFreq==0) // reparam levelset function
            {
                lset.ReparamFastMarching( C.RepMethod);
                
                if (C.ref_freq != 0)
                     adap.UpdateTriang( Stokes, lset);
                if (adap.WasModified()) {
                    cpl.Update();
                    // don't forget to update the pr mass/stiff matrix for the schur compl. preconditioner!!
                    // This is now in cpl.Update().
                }
                
                std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                
                if (C.VolCorr)
                {
                    double dphi= lset.AdjustVolume( Vol, 1e-9);
                    std::cerr << "volume correction is " << dphi << std::endl;
                    lset.Phi.Data+= dphi;
                    std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
                }
            }
            ensight.putGeom( datgeo, step*C.dt);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), step*C.dt);
            ensight.putVector( datvec, Stokes.GetVelSolution(), step*C.dt);
            ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
            ensight.Commit();
        }
    } 
    // no Baensch scheme (FracStep + OpSplitting) 
    
    ensight.CaseEnd();
    std::cerr << std::endl;
}

} // end of namespace DROPS


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
        bnd_fun[i]= (bc[i]= builder.GetBC( i))==DROPS::DirBC ? &Inflow : &DROPS::ZeroVel;
        std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
    }

    EllipsoidCL::Init( C.Mitte, C.Radius );
    DROPS::AdapTriangCL adap( mg, C.ref_width, 0, C.ref_flevel);
//    adap.MakeInitialTriang( EllipsoidCL::DistanceFct);
    for (int i=0; i<C.ref_flevel; ++i)
    {
        DROPS::MarkInterface( EllipsoidCL::DistanceFct, C.ref_width, mg);
        mg.Refine();
    }
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;

    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL(
        num_bnd, bc, bnd_fun), DROPS::P1X_FE, C.XFEMStab);

    Strategy( prob, adap);    // do all the stuff

    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
