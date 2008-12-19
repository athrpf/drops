/// \file
/// \brief flow in measurement cell or brick
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

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
#include "poisson/transport2phase.h"
#include "num/stokessolverfactory.h"
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
    LevelsetP2CL lset( MG, sigmap, gsigmap, C.lset_theta, C.lset_SD,
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
 //   Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
//   Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);
    // For a two-level MG-solver: P2P1 -- P2P1X; comment out the preceeding CreateNumberings
    Stokes.SetNumVelLvl ( 2);
    Stokes.SetNumPrLvl  ( 2);
    Stokes.vel_idx.GetCoarsest().CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
    Stokes.vel_idx.GetFinest().  CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
    Stokes.pr_idx.GetCoarsest(). GetXidx().SetBound( 1e99);
    Stokes.pr_idx.GetCoarsest(). CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Pr, 0, &lset.Phi);
    Stokes.pr_idx.GetFinest().   CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Pr, 0, &lset.Phi);

    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);
    switch (C.IniCond)
    {
      case  1: //flow without droplet
          lset.Init( &One);
      break;
      case -1: // read from file
      {
        ReadEnsightP2SolCL reader( MG);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        // reads the P1-part of the pressure
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        // reads the P1X-part of the pressure
        if (Stokes.UsesXFEM()) {
            std::ifstream fff  ( (C.IniData+".prx").c_str());
            std::ifstream fxidx( (C.IniData+".xidx").c_str());
            if (fff && fxidx) {
                size_t NumP1XUnknowns;
                fff >> NumP1XUnknowns;
                if (NumP1XUnknowns != (pidx->NumUnknowns() - Stokes.GetXidx().GetNumUnknownsStdFE()))
                    throw (DROPSErrCL("error while reading P1X unknowns"));
                for (Uint i=Stokes.GetXidx().GetNumUnknownsStdFE(); i < pidx->NumUnknowns(); ++i)
                    fff >> Stokes.p.Data[i];
                for (Uint k=0; k<Stokes.GetXidx().GetNumUnknownsStdFE(); ++k)
                    fxidx >> Stokes.GetXidx()[k];
            }
        }
      } break;
      default:
        lset.Init( EllipsoidCL::DistanceFct);
    }
    MG.SizeInfo( std::cerr);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";

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

    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

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
        navstokessolver = new AdaptFixedPtDefectCorrCL<StokesProblemT,DeltaSquaredPolicyCL>(Stokes, *stokessolver, C.ns_iter, C.ns_tol, C.ns_red);

    // Time discretisation + coupling
    TimeDisc2PhaseCL<StokesProblemT>* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, C);
    if (C.num_steps != 0) timedisc->SetTimeStep( C.dt);

    if (C.nonlinear!=0.0 || C.num_steps == 0) {
        stokessolverfactory.SetMatrixA( &navstokessolver->GetAN()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( &navstokessolver->GetAN()->GetCoarsest(), &Stokes.B.Data.GetCoarsest(),
                                         &Stokes.M.Data.GetCoarsest(), &Stokes.prM.Data.GetCoarsest());
    }
    else {
        stokessolverfactory.SetMatrixA( &timedisc->GetUpperLeftBlock()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( &timedisc->GetUpperLeftBlock()->GetCoarsest(), &Stokes.B.Data.GetCoarsest(),
                                         &Stokes.M.Data.GetCoarsest(), &Stokes.prM.Data.GetCoarsest());
    }        


    UpdateProlongationCL PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    // For a two-level MG-solver: P2P1 -- P2P1X;
    MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
        Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
        stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());
    stokessolverfactory.GetVankaSchurPc().SetAB( (C.nonlinear!=0.0 || C.num_steps == 0)
        ? &navstokessolver->GetAN()->GetFinest() : &timedisc->GetUpperLeftBlock()->GetFinest(),
        &Stokes.B.Data.GetFinest()
    );
    stokessolverfactory.GetVankaSmoother().SetRelaxation( 0.8);

    bool second = false;
    std::ofstream infofile((C.EnsCase+".info").c_str());
    double lsetmaxGradPhi, lsetminGradPhi;
    IFInfo.WriteHeader(infofile);
    if (C.num_steps == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // Initialize Ensight6 output
    std::string ensf( C.EnsDir + "/" + C.EnsCase);
    Ensight6OutCL ensight( C.EnsCase + ".case", C.num_steps + 1);
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   "Messzelle",     ensf + ".geo", true));
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

    if (C.EnsCase != "none") ensight.Write( 0.);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== step " << step << ":\n";

        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t, infofile);
        timedisc->DoStep( C.cpl_iter);
        if (C.transp_do) c.DoStep( step*C.dt);

        // WriteMatrices( Stokes, step);
        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

        bool forceVolCorr= false, forceUpdate= false,
             doReparam= C.RepFreq && step%C.RepFreq == 0,
             doGridMod= C.ref_freq && step%C.ref_freq == 0;

        // volume correction before reparam/grid modification
        if (C.VolCorr && (doReparam || doGridMod)) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
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
            if (C.serialization_file != "none") {
                std::stringstream filename;
                filename << C.serialization_file;
                if (second) filename << "0";
                second = !second;
                MGSerializationCL ser( MG, filename.str().c_str());
                ser.WriteMG();
            }
        }

        // volume correction
        if (C.VolCorr && (step%C.VolCorr==0 || forceVolCorr)) {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            forceUpdate  = true;
        }

        // update
        if (forceUpdate) {
            timedisc->Update();
            if (C.transp_do) c.Update();
        }

        if (C.EnsCase != "none") ensight.Write( step*C.dt);
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.t, infofile);
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

        DROPS::ReadMeshBuilderCL builder (meshfile);
        if (C.deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( builder);
        else {
            DROPS::FileBuilderCL filebuilder( C.deserialization_file, &builder);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        const DROPS::BoundaryCL& bnd= mgp->GetBnd();
        const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

        DROPS::BndCondT* bc = new DROPS::BndCondT[num_bnd];
        DROPS::StokesVelBndDataCL::bnd_val_fun* bnd_fun = new DROPS::StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= builder.GetBC( i))==DROPS::DirBC ? &InflowCell : &DROPS::ZeroVel;
            std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
        }
        bnddata = new DROPS::StokesBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
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
            DROPS::DROPSErrCL("error while reading geometry information: " + mesh);
        C.r_inlet= dx/2;
        DROPS::Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;
        DROPS::BrickBuilderCL builder ( orig, px, py, pz, nx, ny, nz);
        if (C.deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( builder);
        else {
            DROPS::FileBuilderCL filebuilder( C.deserialization_file, &builder);
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
                bfun[2]=
                bfun[3]= &InflowBrick;
            } break;
            default: throw DROPS::DROPSErrCL("Unknown boundary data type");
        }
        bnddata = new DROPS::StokesBndDataCL(6, bc, bfun);
    }
}

int main (int argc, char** argv)
{
  try
  {
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

    CreateGeom(mg, bnddata);
    EllipsoidCL::Init( C.Mitte, C.Radius);

    DROPS::AdapTriangCL adap( *mg, C.ref_width, 0, C.ref_flevel);

    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (C.deserialization_file == "none")
        adap.MakeInitialTriang( EllipsoidCL::DistanceFct);

    std::cerr << DROPS::SanityMGOutCL(*mg) << std::endl;
    MyStokesCL prob(*mg, ZeroFlowCL(C), *bnddata, C.XFEMStab<0 ? DROPS::P1_FE : DROPS::P1X_FE, C.XFEMStab);

    Strategy( prob, adap);    // do all the stuff

    delete mg;
    delete bnddata;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
