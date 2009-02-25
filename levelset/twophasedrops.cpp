/// \file
/// \brief flow in measurement cell or brick
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "out/output.h"
#include "out/ensightOut.h"
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

#ifdef _PAR
/// \brief Display a detailed list of unknowns
template <typename StokesT, typename LevelsetT>
  void DisplayUnks(const StokesT& Stokes, const LevelsetT& levelset, const DROPS::MultiGridCL& MG)
/** This functions write information about unknowns on the display. These
    informations are for the level-set-, pressure- and velocity-DOF:
    <ul>
     <li> global DOF
     <li> accumulated DOF
     <li> max and min DOF on a single processor (and the ratio)
     <li> max and min number of distributed DOF on a processor (and the ratio to the remaining DOF)
    </ul>
*/
{
    using namespace DROPS;
    const MLIdxDescCL* vidx = &Stokes.vel_idx,
                     * pidx = &Stokes.pr_idx;
    const IdxDescCL*   lidx = &levelset.idx;
    const ExchangeCL& ExV = Stokes.vel_idx.GetEx(),
                    & ExP = Stokes.pr_idx.GetEx(),
                    & ExL = levelset.idx.GetEx();

    // local number on unknowns
    Ulint Psize      = pidx->NumUnknowns();
    Ulint Vsize      = vidx->NumUnknowns();
    Ulint Lsize      = lidx->NumUnknowns();

    // global number of unknowns
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);
    Ulint GLsize     = lidx->GetGlobalNumUnknowns(MG);

    // accumulated size of unknwons
    Ulint Psize_acc = ProcCL::GlobalSum(Psize);
    Ulint Vsize_acc = ProcCL::GlobalSum(Vsize);
    Ulint Lsize_acc = ProcCL::GlobalSum(Lsize);

    // maximal and minimal number of unknowns
    Ulint P_min= ProcCL::GlobalMin(Psize); Ulint P_max= ProcCL::GlobalMax(Psize);
    Ulint V_min= ProcCL::GlobalMin(Vsize); Ulint V_max= ProcCL::GlobalMax(Vsize);
    Ulint L_min= ProcCL::GlobalMin(Lsize); Ulint L_max= ProcCL::GlobalMax(Lsize);

    // ratios between maximal number of unknowns/proc and minimal number
    double P_ratio   = (double)P_max/(double)P_min;
    double V_ratio   = (double)V_max/(double)V_min;
    double L_ratio   = (double)L_max/(double)L_min;

    // number on boundaries
    Ulint P_accmax=ProcCL::GlobalMax(ExP.AccDistIndex.size()), P_accmin=ProcCL::GlobalMin(ExP.AccDistIndex.size());
    Ulint V_accmax=ProcCL::GlobalMax(ExV.AccDistIndex.size()), V_accmin=ProcCL::GlobalMin(ExV.AccDistIndex.size());
    Ulint L_accmax=ProcCL::GlobalMax(ExL.AccDistIndex.size()), L_accmin=ProcCL::GlobalMin(ExL.AccDistIndex.size());

    // ratio of these unknowns
    double P_accratio= (double)P_accmax / (double)P_accmin;
    double V_accratio= (double)V_accmax / (double)V_accmin;
    double L_accratio= (double)L_accmax / (double)L_accmin;

    // output on screen
    std::cerr << "  + Number of DOF\n        "
                << std::setw(10)<<"global"<<std::setw(10)<<"accum"<<std::setw(10)
                << "max"<<std::setw(10)<<"min"<<std::setw(10)<<"ratio"<<"  |  "
                << std::setw(10)<<"max_acc" <<std::setw(10)<<"min_acc"<<std::setw(10)<<"ratio_acc"<<std::endl;

    std::cerr << "    "<<"pr  "
                << std::setw(10)<<GPsize<<std::setw(10)<<Psize_acc<<std::setw(10)<<P_max
                << std::setw(10)<<P_min<< std::setw(10)<<P_ratio<<"  |  "
                << std::setw(10)<<P_accmax<<std::setw(10)<<P_accmin<<std::setw(10)<<P_accratio<<std::endl;

    std::cerr << "    "<<"vel "
                << std::setw(10)<<GVsize<<std::setw(10)<<Vsize_acc<<std::setw(10)<<V_max
                << std::setw(10)<<V_min<< std::setw(10)<<V_ratio<<"  |  "
                << std::setw(10)<<V_accmax<<std::setw(10)<<V_accmin<<std::setw(10)<<V_accratio<<std::endl;

    std::cerr << "    "<<"scl "
                << std::setw(10)<<GLsize<<std::setw(10)<<Lsize_acc<<std::setw(10)<<L_max
                << std::setw(10)<<L_min<< std::setw(10)<<L_ratio<<"  |  "
                << std::setw(10)<<L_accmax<<std::setw(10)<<L_accmin<<std::setw(10)<<L_accratio<<std::endl;

    std::cerr << std::endl;
}

void DisplayDetailedGeom(DROPS::MultiGridCL& mg)
{
    const DROPS::Uint level=mg.GetLastLevel();
    DROPS::Uint *numTetrasAllProc=0;
    DROPS::Uint *numFacesAllProc=0;
    DROPS::Uint *numDistFaceAllProc=0;
    if (DROPS::ProcCL::IamMaster()){
        numTetrasAllProc  = new DROPS::Uint[DROPS::ProcCL::Size()];
        numFacesAllProc   = new DROPS::Uint[DROPS::ProcCL::Size()];
        numDistFaceAllProc= new DROPS::Uint[DROPS::ProcCL::Size()];
    }
    // Gather information about distribution on master processor
    DROPS::ProcCL::Gather(mg.GetNumTriangTetra(level),      numTetrasAllProc,   DROPS::ProcCL::Master());
    DROPS::ProcCL::Gather(mg.GetNumTriangFace(level),       numFacesAllProc,    DROPS::ProcCL::Master());
    DROPS::ProcCL::Gather(mg.GetNumDistributedFaces(level), numDistFaceAllProc, DROPS::ProcCL::Master());

    // Display information
    if (DROPS::ProcCL::IamMaster()){
        double ratioTetra       =  (double)*std::max_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size())
                                  /(double)*std::min_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size());
        DROPS::Uint allTetra    =  std::accumulate(numTetrasAllProc, numTetrasAllProc+DROPS::ProcCL::Size(), 0),
                    allFace     =  std::accumulate(numFacesAllProc, numFacesAllProc+DROPS::ProcCL::Size(), 0),
                    allDistFace =  std::accumulate(numDistFaceAllProc, numDistFaceAllProc+DROPS::ProcCL::Size(), 0);
        double      *ratioDistFace=new double[DROPS::ProcCL::Size()];

        // global information
        std::cerr << "Detailed information about the parallel multigrid:\n"
                  << "#(master tetras on finest level):    "<<allTetra<<'\n'
                  << "#(all Faces on finest level):        "<<allFace<<'\n'
                  << "#(distributed Faces on fines level): "<<allDistFace<<'\n';

        // local information for all processors
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            ratioDistFace[i]= ((double)numDistFaceAllProc[i]/(double)numFacesAllProc[i]*100.);

        double maxRatio= *std::max_element(ratioDistFace, ratioDistFace+DROPS::ProcCL::Size());
        std::cerr << "Ratio between max/min Tetra: "<<ratioTetra
                  <<" max Ratio DistFace/AllFace: "<<maxRatio<<std::endl;

        std::cerr << std::setw(6)  <<  "Proc"
                  << std::setw(8)  << "#Tetra"
                  << std::setw(8)  << "#Faces"
                  << std::setw(12) << "#DistFaces"
                  << std::setw(12) << "%DistFaces"
                  << '\n';
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            std::cerr << std::setw(6)  << i
                      << std::setw(8)  << numTetrasAllProc[i]
                      << std::setw(8)  << numFacesAllProc[i]
                      << std::setw(12) << numDistFaceAllProc[i]
                      << std::setw(12) << ratioDistFace[i] << std::endl;

        // free memory
        if (numTetrasAllProc)   delete[] numTetrasAllProc;
        if (numFacesAllProc)    delete[] numFacesAllProc;
        if (numDistFaceAllProc) delete[] numDistFaceAllProc;
        if (ratioDistFace)      delete[] ratioDistFace;
    }
}
#endif


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
        case 2 :
            return (new RecThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_proj, C.cpl_stab));
        break;
        case 3 :
            return (new ThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
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
    std::cerr << "Discretizing took "<< time.GetTime() << " sec.\n";
    time.Reset();
    Stokes.b.Data += curv.Data;
    solver.Solve( Stokes.A.Data, Stokes.B.Data, Stokes.v, Stokes.p.Data, Stokes.b.Data, cplN, Stokes.c.Data, 1.0);
    time.Stop();
    std::cerr << "Solving (Navier-)Stokes took "<< time.GetTime() << " sec.\n";
    std::cerr << "iter: " << solver.GetIter() << "\tresid: " << solver.GetResid() << std::endl;
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
    // For a two-level MG-solver: P2P1 -- P2P1X; comment out the preceeding CreateNumberings
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
      case -1: // read from file
      {
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
            reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
      } break;
      case 0: // zero initial condition
          lset.Init( EllipsoidCL::DistanceFct);
        break;
      case 1: // stationary flow
      {
        lset.Init( EllipsoidCL::DistanceFct);
        SSORPcCL ssorpc;
        PCG_SsorCL PCGsolver( ssorpc, 200, 1e-2, true);
        typedef SolverAsPreCL<PCG_SsorCL> PCGPcT;
        PCGPcT apc( PCGsolver);
#ifdef _PAR
        ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), Stokes.vel_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
#else
        ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
#endif
        InexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM> inexactuzawasolver( apc, bbtispc, C.outer_iter, C.outer_tol, 0.6);
        NSSolverBaseCL<StokesProblemT> stokessolver( Stokes, inexactuzawasolver);
        SolveStatProblem( Stokes, lset, stokessolver);
      } break;
      case  2: //flow without droplet
          lset.Init( &One);
      break;
      default : throw DROPSErrCL("Unknown initial condition");
    }

#ifndef _PAR
    MG.SizeInfo( std::cerr);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";
#else
    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);
#endif

    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "initial volume: " << lset.GetVolume()/Vol << std::endl;
    double dphi= lset.AdjustVolume( Vol, 1e-9);
    std::cerr << "initial volume correction is " << dphi << std::endl;
    lset.Phi.Data+= dphi;
    std::cerr << "new initial volume: " << lset.GetVolume()/Vol << std::endl;

    cBndDataCL Bnd_c( 6, c_bc, c_bfun);
    double D[2] = {C.transp_cPos, C.transp_cNeg};
    TransportP1CL c( MG, Bnd_c, Stokes.GetBndData().Vel, C.transp_theta, D, C.transp_H, &Stokes.v, lset,
        /*t*/ 0., C.dt, C.transp_iter, C.transp_tol);
    TransportRepairCL transprepair(c, MG);
    if (C.transp_do)
    {
        adap.push_back(&transprepair);
        MLIdxDescCL* cidx= &c.idx;
        c.CreateNumbering( MG.GetLastLevel(), cidx);
        c.ct.SetIdx( cidx);
        if (C.IniCond != -1)
            c.Init( &Initialcneg, &Initialcpos);
        else
        {
            ReadEnsightP2SolCL reader( MG);
            reader.ReadScalar( C.IniData+".ct", c.ct, c.GetBndData());
        }
        c.Update();
        std::cerr << c.c.Data.size() << " concentration unknowns,\n";
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

    bool second = false;
    double lsetmaxGradPhi, lsetminGradPhi;
    std::ofstream* infofile = 0;
    IF_MASTER {
        infofile = new std::ofstream ((C.EnsCase+".info").c_str());
        IFInfo.WriteHeader(*infofile);
    }

    if (C.num_steps == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // Initialize Ensight6 output
#ifndef _PAR
    std::string ensf( C.EnsDir + "/" + C.EnsCase);
    Ensight6OutCL ensight( C.EnsCase + ".case", C.num_steps + 1, C.binary);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   C.geomName,      ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
    ensight.Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));
    ensight.Register( make_Ensight6Scalar    ( ScalarFunAsP2EvalCL( sigmap, 0., &MG, MG.GetLastLevel()),
                                                                        "Surfaceforce",  ensf + ".sf",  true));
    if (C.transp_do) {
        ensight.Register( make_Ensight6Scalar( c.GetSolution(),         "Concentration", ensf + ".c",   true));
        ensight.Register( make_Ensight6Scalar( c.GetSolution( c.ct),    "TransConc",     ensf + ".ct",  true));
    }
    if (Stokes.UsesXFEM())
        ensight.Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));

    if (C.ensight) ensight.Write( 0.);

#else
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.EnsDir, C.EnsCase, C.geomName, /*adaptive=*/true,
                                  C.num_steps, C.binary, C.masterOut);
#endif
    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "============================================================ step " << step << std::endl;

        IFInfo.Update( lset, Stokes.GetVelSolution());

        IF_MASTER
            IFInfo.Write(Stokes.t, *infofile);

        timedisc->DoStep( C.cpl_iter);
        if (C.transp_do) c.DoStep( step*C.dt);

        // WriteMatrices( Stokes, step);
        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

        bool forceVolCorr= false, forceUpdate= false,
             doReparam= C.RepFreq && step%C.RepFreq == 0,
             doGridMod= C.ref_freq && step%C.ref_freq == 0;

        // volume correction before reparam/grid modification
        if (C.VolCorr && (doReparam || doGridMod)) {
                dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. volume: " << lset.GetVolume()/Vol << std::endl;
                forceUpdate = true; // volume correction modifies the level set
        }

        // reparam levelset function
        if (doReparam) {
            lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
            if (lsetmaxGradPhi > C.MaxGrad || lsetminGradPhi < C.MinGrad) {
                std::cerr << "before reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
                lset.ReparamFastMarching( C.RepMethod);
                lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
                std::cerr << "after  reparametrization: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n';
                forceVolCorr = forceUpdate = true; // volume correction and update after reparam
            }
        }

        // grid modification
        if (doGridMod) {
            adap.UpdateTriang( lset);
            forceUpdate  |= adap.WasModified();
            forceVolCorr |= adap.WasModified();
            if (C.serialization) {
                std::stringstream filename;
                filename << C.ser_dir;
                if (second) filename << "0";
                second = !second;
                MGSerializationCL ser( MG, filename.str().c_str());
                ser.WriteMG();
            }
        }

        // volume correction
        if (C.VolCorr && (step%C.VolCorr==0 || forceVolCorr)) {
            dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. volume: " << lset.GetVolume()/Vol << std::endl;
            forceUpdate  = true;
        }

        // update
        if (forceUpdate) {
            timedisc->Update();
            if (C.transp_do) c.Update();
        }

#ifndef _PAR
        if (C.ensight) ensight.Write( step*C.dt);
#else
        if (C.ensight) ensightwriter.write();
#endif
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IF_MASTER
        IFInfo.Write(Stokes.t, *infofile);
    std::cerr << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
//     delete stokessolver1;
}

} // end of namespace DROPS

void CreateGeom (DROPS::MultiGridCL* &mgp, DROPS::StokesBndDataCL* &bnddata)
{
    if (C.GeomType == 0) {
        std::ifstream meshfile( C.meshfile.c_str());
        if (!meshfile)
            throw DROPS::DROPSErrCL ("error while opening mesh file\n");

        DROPS::ReadMeshBuilderCL *mgb= 0;       // builder of the multigrid

        // read geometry information from a file and create the multigrid
        IF_MASTER
            mgb = new DROPS::ReadMeshBuilderCL( meshfile );
        IF_NOT_MASTER
            mgb = new DROPS::EmptyReadMeshBuilderCL( meshfile );
        // Create the multigrid
        if (C.deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( *mgb);
        else {
#ifdef _PAR
            throw DROPS::DROPSErrCL( "Sorry, no parallel deserialization yet");
#endif
            DROPS::FileBuilderCL filebuilder( C.deserialization_file, mgb);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        const DROPS::BoundaryCL& bnd= mgp->GetBnd();
        const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

        DROPS::BndCondT* bc = new DROPS::BndCondT[num_bnd];
        DROPS::StokesVelBndDataCL::bnd_val_fun* bnd_fun = new DROPS::StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= mgb->GetBC( i))==DROPS::DirBC ? &InflowCell : &DROPS::ZeroVel;
            std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
        }
        bnddata = new DROPS::StokesBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
        delete   mgb;
    }
    if (C.GeomType == 1) {
        int nx, ny, nz;
        double dx, dy, dz;
        std::string mesh( C.meshfile), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
        if (!brick_info)
            throw DROPS::DROPSErrCL("error while reading geometry information: " + mesh);
        C.r_inlet= dx/2;
        DROPS::Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;

        DROPS::BrickBuilderCL *mgb = 0;
        IF_MASTER
            mgb = new DROPS::BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
        IF_NOT_MASTER
            mgb = new DROPS::EmptyBrickBuilderCL(orig, px, py, pz);

        if (C.deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( *mgb);
        else {
#ifdef _PAR
            throw DROPS::DROPSErrCL( "Sorry, no parallel deserialization yet");
#endif
            DROPS::FileBuilderCL filebuilder( C.deserialization_file, mgb);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        DROPS::BndCondT bc[6]= { DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC };
        DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[6]=
            { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel };
        switch (C.bnd_type)
        {
            case 1 : // hom. Dirichlet bnd-data
              break;
            case 2 : // inhom. Dirichlet bnd-data for in-/outflow
            {
                bc[2]= bc[3]= DROPS::DirBC;
                bfun[2]= bfun[3]= &InflowBrick;
            } break;
            case 3 : // tube/canal
            {
                bc[3]= DROPS::DirBC;
                bc[2]= DROPS::NatBC; //Rohr
                //bc[2]= bc[4]= bc[5]= DROPS::NatBC;          //Kanal
                bfun[2]= &DROPS::ZeroVel;
                //bfun[2]=bfun[4]=bfun[5]= &DROPS::ZeroVel;   //Kanal
                bfun[3]= &InflowBrick;
            } break;
            default: throw DROPS::DROPSErrCL("Unknown boundary data type");
        }
        bnddata = new DROPS::StokesBndDataCL(6, bc, bfun);
        delete mgb;
    }
}

#ifdef _PAR
void DistributeGeom( DROPS::MultiGridCL& mg, DROPS::ParMultiGridCL* pmg, DROPS::LoadBalHandlerCL* &lb)
{
    // Create the multigrid and tell parallel multigrid about the geometry
    pmg->AttachTo( mg);

    // Create a load balancer and do initial distribution of the geometry
    lb = new DROPS::LoadBalHandlerCL( mg);
    lb->DoInitDistribution(DROPS::ProcCL::Master());
    int refineStrategy = 1;
    switch (refineStrategy) {
        case 0 : lb->SetStrategy(DROPS::NoMig);     break;
        case 1 : lb->SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb->SetStrategy(DROPS::Recursive); break;
    }
}
#endif

int main (int argc, char** argv)
{
  try
  {
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
#endif
    std::ifstream param;
    if (argc!=2)
    {
        std::cerr << "Using default parameter file: risingdroplet.param\n";
        param.open( "risingdroplet.param");
    }
    else
        param.open( argv[1]);
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

    DROPS::MultiGridCL* mg= 0;
    DROPS::StokesBndDataCL* bnddata= 0;
#ifdef _PAR
    DROPS::ParMultiGridCL *pmg= DROPS::ParMultiGridCL::InstancePtr();
    DROPS::LoadBalHandlerCL *lb=  0;
#endif
    CreateGeom(mg, bnddata);
    EllipsoidCL::Init( C.Mitte, C.Radius);
#ifdef _PAR
    DistributeGeom( *mg, pmg, lb);
    DROPS::AdapTriangCL adap( *pmg, *lb, C.ref_width, 0, C.ref_flevel);
#else
    DROPS::AdapTriangCL adap( *mg, C.ref_width, 0, C.ref_flevel);
#endif

    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (C.deserialization_file == "none")
        adap.MakeInitialTriang( EllipsoidCL::DistanceFct);

    bool mgok = false;
#ifdef _PAR
    DROPS::ProcCL::Check( CheckParMultiGrid(*pmg));
#endif
    std::cerr << DROPS::SanityMGOutCL(*mg) << std::endl;
    if (mgok)
        std::cerr << "As far as I can tell the ParMultigridCl is sane\n";
    MyStokesCL prob( *mg, ZeroFlowCL(C), *bnddata, C.XFEMStab<0 ? DROPS::P1_FE : DROPS::P1X_FE, C.XFEMStab);

    Strategy( prob, adap);    // do all the stuff

    delete mg;
    delete bnddata;
#ifdef _PAR
    if (pmg)     delete pmg;
    if (lb)      delete lb;
#endif
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
