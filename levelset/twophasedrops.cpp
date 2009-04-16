/// \file
/// \brief flow in measurement cell or brick
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "poisson/transport2phase.h"
#include "num/stokessolverfactory.h"
#ifndef _PAR
#include "num/stokessolver.h"
#else
#include "num/parstokessolver.h"
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
#include <fstream>
#include <sstream>


DROPS::ParamMesszelleNsCL C;
// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0


//brickflow.cpp + brick_transp.cpp + brick_ns_adap.cpp
DROPS::SVectorCL<3> InflowBrick( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = p[0]*(2*C.r_inlet-p[0]) / (C.r_inlet*C.r_inlet),
                 z = p[2]*(2*C.r_inlet-p[2]) / (C.r_inlet*C.r_inlet);

    ret[1]= x * z * C.Anstroem * (1-C.inflow_ampl*std::cos(2*M_PI*C.inflow_freq*t));
    return ret;
}

//mzelle_ns_adap.cpp + mzelle_instat.cpp
DROPS::SVectorCL<3> InflowCell( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;
const DROPS::BndCondT c_bc[6]= {
    DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC,
    DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC
};
const c_bnd_val_fun c_bfun[6]= {0, 0, 0, 0, 0, 0};

double Initialcneg (const DROPS::Point3DCL& , double)
{
    return C.transp_cNeg;
}

double Initialcpos (const DROPS::Point3DCL& , double)
{
    return C.transp_cNeg;
}

namespace DROPS // for Strategy
{

template<class Coeff>
void WriteMatrices (InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, int i)
{
    std::string path( "matrices/");
    std::ostringstream suffix;
    suffix << std::setfill( '0') << std::setw( 4) << i << ".txt";
    WriteToFile( Stokes.A.Data.GetFinest(),   path + "A"   + suffix.str(), "A");
    WriteToFile( Stokes.B.Data.GetFinest(),   path + "B"   + suffix.str(), "B");
    WriteToFile( Stokes.M.Data.GetFinest(),   path + "M"   + suffix.str(), "M");
    WriteToFile( Stokes.prA.Data.GetFinest(), path + "prA" + suffix.str(), "prA");
    WriteToFile( Stokes.prM.Data.GetFinest(), path + "prM" + suffix.str(), "prM");
    WriteToFile( Stokes.N.Data.GetFinest(),   path + "N"   + suffix.str(), "N");

    WriteToFile( Stokes.v.Data, path + "v" + suffix.str(), "v");
    WriteToFile( Stokes.p.Data, path + "p" + suffix.str(), "p");
}

template< class StokesProblemT>
TimeDisc2PhaseCL<StokesProblemT>* CreateTimeDisc(StokesProblemT& Stokes, LevelsetP2CL& lset,
    NSSolverBaseCL<StokesProblemT>* solver, ParamMesszelleNsCL& C)
{
    if (C.num_steps == 0) return 0;
    switch (C.scheme)
    {
        case 1 :
            return (new LinThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_stab));
        break;
        case 3 :
            std::cout << "[WARNING] use of ThetaScheme2PhaseCL is deprecated using RecThetaScheme2PhaseCL instead\n";
        case 2 :
            return (new RecThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_proj, C.cpl_stab));
        break;
        case 4 :
            return (new OperatorSplitting2PhaseCL<StokesProblemT, StokesSolverBaseCL>
                        (Stokes, lset, solver->GetStokesSolver(), C.inner_iter, C.inner_tol, C.nonlinear));
        break;
        case 5 :
            return (new CrankNicolsonScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.nonlinear, C.cpl_proj, C.cpl_stab));
        break;
        default : throw DROPSErrCL("Unknown TimeDiscMethod");
    }
}


template <class Coeff>
void SolveStatProblem( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset,
                       NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL<Coeff> >& solver)
{
    TimerCL time;
    time.Reset();
    VelVecDescCL cplM, cplN;
    VecDescCL curv;
    cplM.SetIdx( &Stokes.vel_idx);
    cplN.SetIdx( &Stokes.vel_idx);
    curv.SetIdx( &Stokes.vel_idx);
    Stokes.SetIdx();
    Stokes.SetLevelSet( lset);
    lset.AccumulateBndIntegral( curv);
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &cplM, lset, Stokes.t);
    Stokes.SetupPrStiff( &Stokes.prA, lset);
    Stokes.SetupPrMass ( &Stokes.prM, lset);
    Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
    time.Stop();
    std::cout << "Discretizing took "<< time.GetTime() << " sec.\n";
    time.Reset();
    Stokes.b.Data += curv.Data;
    solver.Solve( Stokes.A.Data, Stokes.B.Data, Stokes.v, Stokes.p.Data, Stokes.b.Data, cplN, Stokes.c.Data, 1.0);
    time.Stop();
    std::cout << "Solving (Navier-)Stokes took "<< time.GetTime() << " sec.\n";
    std::cout << "iter: " << solver.GetIter() << "\tresid: " << solver.GetResid() << std::endl;
}

// For a two-level MG-solver: P2P1 -- P2P1X; canonical prolongations
void MakeP1P1XProlongation (size_t NumUnknownsVel, size_t NumUnknownsPr, size_t NumUnknownsPrP1,
    MatrixCL& PVel, MatrixCL& PPr)
{
    // finest level
    //P2-Prolongation (Id)
    PVel= MatrixCL( std::valarray<double>(  1.0, NumUnknownsVel));
    //P1-P1X-Prolongation
    VectorCL diag( 0., NumUnknownsPr);
    diag[std::slice(0, NumUnknownsPrP1, 1)]= 1.;
    PPr= MatrixCL( diag);
}

template<class Coeff>
void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adap)
// flow control
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    sigma= Stokes.GetCoeff().SurfTens;
    eps= C.st_jumpWidth;    lambda= C.st_relPos;    sigma_dirt_fac= C.st_red;
    instat_scalar_fun_ptr sigmap  = 0;
    instat_vector_fun_ptr gsigmap = 0;
    if (C.st_var)
    {
        sigmap  = &sigma_step;
        gsigmap = &gsigma_step;
    }
    else
    {
        sigmap  = &sigmaf;
        gsigmap = &gsigma;
    }
    LevelsetP2CL lset( MG, sigmap, gsigmap, C.theta, C.lset_SD,
        -1, C.lset_iter, C.lset_tol, C.CurvDiff);

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset);
    adap.push_back( &prrepair);

    IdxDescCL* lidx= &lset.idx;
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    if (C.st_var)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    if ( StokesSolverFactoryHelperCL<ParamMesszelleNsCL>().VelMGUsed(C))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL<ParamMesszelleNsCL>().PrMGUsed(C))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);
    // For a two-level MG-solver: P2P1 -- P2P1X; comment out the preceding CreateNumberings
//     Stokes.SetNumVelLvl ( 2);
//     Stokes.SetNumPrLvl  ( 2);
//     Stokes.vel_idx.GetCoarsest().CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
//     Stokes.vel_idx.GetFinest().  CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
//     Stokes.pr_idx.GetCoarsest(). GetXidx().SetBound( 1e99);
//     Stokes.pr_idx.GetCoarsest(). CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Pr, 0, &lset.Phi);
//     Stokes.pr_idx.GetFinest().   CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Pr, 0, &lset.Phi);

    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);
    switch (C.IniCond)
    {
#ifndef _PAR
      case -10: // read from ensight-file [deprecated]
      {
        std::cout << "[DEPRECATED] read from ensight-file [DEPRECATED]\n";
        ReadEnsightP2SolCL reader( MG);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        if (Stokes.UsesXFEM()) {
            VecDescCL pneg( pidx), ppos( pidx);
            reader.ReadScalar( C.IniData+".prNeg", pneg, Stokes.GetBndData().Pr);
            reader.ReadScalar( C.IniData+".prPos", ppos, Stokes.GetBndData().Pr);
            P1toP1X ( pidx->GetFinest(), Stokes.p.Data, pidx->GetFinest(), ppos.Data, pneg.Data, lset.Phi, MG);
        }
        else
            reader.ReadScalar( C.IniData+".pr", Stokes.p, Stokes.GetBndData().Pr);
      } break;
#endif
      case -1: // read from file
      {
        ReadFEFromFile( lset.Phi, MG, C.IniData+"levelset");
        ReadFEFromFile( Stokes.v, MG, C.IniData+"velocity");
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        ReadFEFromFile( Stokes.p, MG, C.IniData+"pressure", false, &lset.Phi); // pass also level set, as p may be extended
      } break;
      case 0: // zero initial condition
          lset.Init( EllipsoidCL::DistanceFct);
        break;
      case 1: // stationary flow
      {
        lset.Init( EllipsoidCL::DistanceFct);
#ifdef _PAR
        ParJac0CL jacpc( Stokes.vel_idx.GetFinest());
        typedef ParPCGSolverCL<ParJac0CL> PCGSolverT;
        typedef SolverAsPreCL<PCGSolverT> PCGPcT;
        PCGSolverT PCGSolver(200, 1e-2, Stokes.vel_idx.GetFinest(), jacpc, /*rel*/ true, /*acc*/ true);
        PCGPcT     apc(PCGSolver);
        ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), Stokes.vel_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
        ParInexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM> inexactuzawasolver( apc, bbtispc, Stokes.vel_idx.GetFinest(), Stokes.pr_idx.GetFinest(),
                                                                           C.outer_iter, C.outer_tol, 0.6, 50, &std::cout);
#else
        SSORPcCL ssorpc;
        PCG_SsorCL PCGsolver( ssorpc, 200, 1e-2, true);
        typedef SolverAsPreCL<PCG_SsorCL> PCGPcT;
        PCGPcT apc( PCGsolver);
        ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
        InexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM> inexactuzawasolver( apc, bbtispc, C.outer_iter, C.outer_tol, 0.6, 50);
#endif
        NSSolverBaseCL<StokesProblemT> stokessolver( Stokes, inexactuzawasolver);
        SolveStatProblem( Stokes, lset, stokessolver);
      } break;
      case  2: //flow without droplet
          lset.Init( &One);
      break;
      default : throw DROPSErrCL("Unknown initial condition");
    }

    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);

    const double Vol= EllipsoidCL::GetVolume();
    std::cout << "initial volume: " << lset.GetVolume()/Vol << std::endl;
    double dphi= lset.AdjustVolume( Vol, 1e-9);
    std::cout << "initial volume correction is " << dphi << std::endl;
    lset.Phi.Data+= dphi;
    std::cout << "new initial volume: " << lset.GetVolume()/Vol << std::endl;

    cBndDataCL Bnd_c( 6, c_bc, c_bfun);
    double D[2] = {C.transp_cPos, C.transp_cNeg};
    TransportP1CL massTransp( MG, Bnd_c, Stokes.GetBndData().Vel, C.transp_theta, D, C.transp_H, &Stokes.v, lset,
        /*t*/ 0., C.dt, C.transp_iter, C.transp_tol);
    TransportRepairCL transprepair(massTransp, MG);
    if (C.transp_do)
    {
        adap.push_back(&transprepair);
        MLIdxDescCL* cidx= &massTransp.idx;
        massTransp.CreateNumbering( MG.GetLastLevel(), cidx);
        massTransp.ct.SetIdx( cidx);
        if (C.IniCond != -1)
            massTransp.Init( &Initialcneg, &Initialcpos);
        else
        {
            ReadFEFromFile( massTransp.ct, MG, C.IniData+"concentrationTransf");
        }
        massTransp.Update();
        std::cout << massTransp.c.Data.size() << " concentration unknowns,\n";
    }
    // Stokes-Solver
    StokesSolverFactoryCL<StokesProblemT, ParamMesszelleNsCL> stokessolverfactory(Stokes, C);
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver();
//     StokesSolverAsPreCL pc (*stokessolver1, 1);
//     GCRSolverCL<StokesSolverAsPreCL> gcr(pc, C.outer_iter, C.outer_iter, C.outer_tol, /*rel*/ false);
//     BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> >* stokessolver =
//             new BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> > (gcr);

    // Navier-Stokes-Solver
    NSSolverBaseCL<StokesProblemT>* navstokessolver = 0;
    if (C.nonlinear==0.0)
        navstokessolver = new NSSolverBaseCL<StokesProblemT>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<StokesProblemT>(Stokes, *stokessolver, C.ns_iter, C.ns_tol, C.ns_red);

    // Time discretisation + coupling
    TimeDisc2PhaseCL<StokesProblemT>* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, C);
    if (C.num_steps != 0) timedisc->SetTimeStep( C.dt);

    if (C.nonlinear!=0.0 || C.num_steps == 0) {
        stokessolverfactory.SetMatrixA( &navstokessolver->GetAN()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( &navstokessolver->GetAN()->GetCoarsest(), &Stokes.B.Data.GetCoarsest(),
                                         &Stokes.M.Data.GetCoarsest(), &Stokes.prM.Data.GetCoarsest(), &Stokes.pr_idx.GetCoarsest());
    }
    else {
        stokessolverfactory.SetMatrixA( &timedisc->GetUpperLeftBlock()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( &timedisc->GetUpperLeftBlock()->GetCoarsest(), &Stokes.B.Data.GetCoarsest(),
                                         &Stokes.M.Data.GetCoarsest(), &Stokes.prM.Data.GetCoarsest(), &Stokes.pr_idx.GetCoarsest());
    }

    UpdateProlongationCL PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    // For a two-level MG-solver: P2P1 -- P2P1X;
//     MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
//         Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
//         stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());

    double lsetmaxGradPhi, lsetminGradPhi;
    std::ofstream* infofile = 0;
    IF_MASTER {
        infofile = new std::ofstream ((C.EnsCase+".info").c_str());
    }
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    if (C.num_steps == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // for serialization of geometry and numerical data
    TwoPhaseStoreCL<StokesProblemT> ser(MG, Stokes, lset, C.transp_do ? &massTransp : 0, C.ser_dir);

    // Initialize Ensight6 output
#ifndef _PAR
    std::string ensf( C.EnsDir + "/" + C.EnsCase);
    Ensight6OutCL ensight( C.EnsCase + ".case", (C.ensight ? C.num_steps/C.ensight+1 : 0), C.binary);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   C.geomName,      ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));
    ensight.Register( make_Ensight6Scalar    ( ScalarFunAsP2EvalCL( sigmap, 0., &MG, MG.GetLastLevel()),
                                                                        "Surfaceforce",  ensf + ".sf",  true));
    if (C.transp_do) {
        ensight.Register( make_Ensight6Scalar( massTransp.GetSolution(),"Concentration", ensf + ".c",   true));
        ensight.Register( make_Ensight6Scalar( massTransp.GetSolution( massTransp.ct),
                                                                        "TransConc",     ensf + ".ct",  true));
    }
    if (Stokes.UsesXFEM())
        ensight.Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));

    if (C.ensight) ensight.Write( 0.);

#else
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.EnsDir, C.EnsCase, C.geomName, /*adaptive=*/true,
                                  (C.ensight? C.num_steps/C.ensight+1 : 0), C.binary, C.masterOut);
    if (C.ensight) ensightwriter.write();
#endif

    // writer for vtk-format
    typedef TwoPhaseVTKCL<StokesProblemT, LevelsetP2CL> VTKWriterT;
    VTKWriterT vtkwriter( adap.GetMG(), Stokes, lset,  (C.vtk ? C.num_steps/C.vtk+1 : 0),
                          std::string(C.vtkDir + "/" + C.vtkName), C.vtkBinary);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cout << "============================================================ step " << step << std::endl;

        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        timedisc->DoStep( C.cpl_iter);
        if (C.transp_do) massTransp.DoStep( step*C.dt);

        // WriteMatrices( Stokes, step);
        std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

        bool forceVolCorr= false, forceUpdate= false,
             doReparam= C.RepFreq && step%C.RepFreq == 0,
             doGridMod= C.ref_freq && step%C.ref_freq == 0;

        // volume correction before reparam/grid modification
        if (C.VolCorr && (doReparam || doGridMod)) {
                dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cout << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cout << "new rel. volume: " << lset.GetVolume()/Vol << std::endl;
                forceUpdate = true; // volume correction modifies the level set
        }

        // reparam levelset function
        if (doReparam) {
            lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
            if (lsetmaxGradPhi > C.MaxGrad || lsetminGradPhi < C.MinGrad) {
                std::cout << "before reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
                lset.ReparamFastMarching( C.RepMethod, false, false, C.RepMethod==3);
                lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
                std::cout << "after  reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
                forceVolCorr = forceUpdate = true; // volume correction and update after reparam
            }
        }

        // grid modification
        if (doGridMod) {
            adap.UpdateTriang( lset);
            forceUpdate  |= adap.WasModified();
            forceVolCorr |= adap.WasModified();
            if (C.serialization)
                ser.Write();
        }

        // volume correction
        if (C.VolCorr && (step%C.VolCorr==0 || forceVolCorr)) {
            dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cout << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cout << "new rel. volume: " << lset.GetVolume()/Vol << std::endl;
            forceUpdate  = true;
        }

        // update
        if (forceUpdate) {
            timedisc->Update();
            if (C.transp_do) massTransp.Update();
        }

        if (C.ensight && step%C.ensight==0)
#ifndef _PAR
            ensight.Write( Stokes.t);
#else
            ensightwriter.write();
#endif
        if (C.vtk && step%C.vtk==0)
            vtkwriter.write();
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.t);
    std::cout << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
    if (infofile) delete infofile;
//     delete stokessolver1;
}

} // end of namespace DROPS

int main (int argc, char** argv)
{
  try
  {
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
#endif
    std::ifstream param;
    if (argc!=2)
    {
        std::cout << "Using default parameter file: risingdroplet.param\n";
        param.open( "risingdroplet.param");
    }
    else
        param.open( argv[1]);
    if (!param)
    {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cout << C << std::endl;

    typedef DROPS::ZeroFlowCL                             CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::MultiGridCL* mg= 0;
    DROPS::StokesBndDataCL* bnddata= 0;

    CreateGeom(mg, bnddata, C.GeomType == 0 ? InflowCell : InflowBrick, C.meshfile, C.GeomType, C.bnd_type, C.deserialization_file, C.r_inlet);
    DROPS::EllipsoidCL::Init( C.Mitte, C.Radius);
    DROPS::AdapTriangCL adap( *mg, C.ref_width, 0, C.ref_flevel, ((C.deserialization_file == "none") ? C.refineStrategy : -1));
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (C.deserialization_file == "none")
        adap.MakeInitialTriang( DROPS::EllipsoidCL::DistanceFct);

    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    if (DROPS::ProcCL::Check( CheckParMultiGrid( adap.GetPMG())))
        std::cout << "As far as I can tell the ParMultigridCl is sane\n";
#endif
    MyStokesCL prob( *mg, DROPS::ZeroFlowCL(C), *bnddata, C.XFEMStab<0 ? DROPS::P1_FE : DROPS::P1X_FE, C.XFEMStab);

    Strategy( prob, adap);    // do all the stuff

    delete mg;
    delete bnddata;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
