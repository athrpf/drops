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

/* Available stokes solver
 *   no | Solver         | APc       | SPc     | PC structure
 ------------------------------------------------------------
 *   11 | GCR            | MG        | BBT     | lower block
 *   12 | GCR            | MG        | MinComm | lower block
 *   13 | GCR            | GMRes     | BBT     | lower block
 *   14 | GCR            | GMRes     | MinComm | lower block
 ------------------------------------------------------------
 *   21 | Inexact Uzawa  | asymm. MG | BBT     |
 *   22 | Inexact Uzawa  | asymm. MG | MinComm |
 *  221 | Inexact Uzawa  | symm. MG  | BBT     |
 *  222 | Inexact Uzawa  | symm. MG  | MinComm |
 *   23 | Inexact Uzawa  | GMRes     | BBT     |
 *   24 | Inexact Uzawa  | GMRes     | MinComm |
 *   25 | Inexact Uzawa  | SSORPCG   | BBT     |
 *   26 | Inexact Uzawa  | SSORPCG   | MinComm |
 ------------------------------------------------------------
 *   31 | PMinRes        | MG        | BBT     | diag
 *   32 | PMinRes        | MG        | MinComm | diag
 *   33 | PMinRes        | SSORPCG   | BBT     | diag
 *   34 | PMinRes        | SSORPCG   | MinComm | diag*/

template <class StokesT>
class StokesSolverFactory
{
  private:
    StokesT& Stokes_;
    int    outer_iter_;
    double outer_tol_;
    double kA_, kM_;
    bool   mgused_;

// PC for instat. Schur complement
    ISBBTPreCL bbtispc_;

    MinCommPreCL mincommispc_;

// PC for A-block
    // MultiGrid
    MGDataCL velMG_;
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> MGSolver_;
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > MGPcT;
    MGPcT MGPc_;

    //JAC-GMRes
    JACPcCL  JACPc_;
    typedef GMResSolverCL<JACPcCL> GMResSolverT;        // GMRes-based APcT
    GMResSolverT GMResSolver_;
    typedef SolverAsPreCL<GMResSolverT> GMResPcT;
    GMResPcT GMResPc_;

    //PCG
    SSORPcCL SSORPc_;
    typedef PCGSolverCL<SSORPcCL> PCGSolverT;           // CG-based APcT
    PCGSolverT PCGSolver_;
    typedef SolverAsPreCL<PCGSolverT> PCGPcT;
    PCGPcT PCGPc_;

// PC for Oseen problem
    typedef BlockPreCL<MGPcT,    ISBBTPreCL>   DiagMGBBTOseenPcT;
    typedef BlockPreCL<MGPcT,    MinCommPreCL> DiagMGMinCommOseenPcT;
    typedef BlockPreCL<PCGPcT,   ISBBTPreCL>   DiagPCGBBTOseenPcT;
    typedef BlockPreCL<PCGPcT,   MinCommPreCL> DiagPCGMinCommOseenPcT;

    typedef BlockPreCL<MGPcT,    ISBBTPreCL,   LowerBlockPreCL> LBlockMGBBTOseenPcT;
    typedef BlockPreCL<MGPcT,    MinCommPreCL, LowerBlockPreCL> LBlockMGMinCommOseenPcT;
    typedef BlockPreCL<GMResPcT, ISBBTPreCL,   LowerBlockPreCL> LBlockGMResBBTOseenPcT;
    typedef BlockPreCL<GMResPcT, MinCommPreCL, LowerBlockPreCL> LBlockGMResMinCommOseenPcT;

    DiagMGBBTOseenPcT         DiagMGBBTOseenPc_;
    DiagMGMinCommOseenPcT     DiagMGMinCommOseenPc_;
    DiagPCGBBTOseenPcT        DiagPCGBBTOseenPc_;
    DiagPCGMinCommOseenPcT    DiagPCGMinCommOseenPc_;

    LBlockMGBBTOseenPcT        LBlockMGBBTOseenPc_;
    LBlockMGMinCommOseenPcT    LBlockMGMinCommOseenPc_;
    LBlockGMResBBTOseenPcT     LBlockGMResBBTOseenPc_;
    LBlockGMResMinCommOseenPcT LBlockGMResMinCommOseenPc_;

// GCR solver
    GCRSolverCL<LBlockMGBBTOseenPcT>        GCRMGBBT_;
    GCRSolverCL<LBlockMGMinCommOseenPcT>    GCRMGMinComm_;
    GCRSolverCL<LBlockGMResBBTOseenPcT>     GCRGMResBBT_;
    GCRSolverCL<LBlockGMResMinCommOseenPcT> GCRGMResMinComm_;

// Lanczos
    typedef PLanczosONBCL<BlockMatrixCL, VectorCL, DiagMGBBTOseenPcT>      Lanczos1T;
    typedef PLanczosONBCL<BlockMatrixCL, VectorCL, DiagMGMinCommOseenPcT>  Lanczos2T;
    typedef PLanczosONBCL<BlockMatrixCL, VectorCL, DiagPCGBBTOseenPcT>     Lanczos3T;
    typedef PLanczosONBCL<BlockMatrixCL, VectorCL, DiagPCGMinCommOseenPcT> Lanczos4T;

    Lanczos1T lanczos1_;
    Lanczos2T lanczos2_;
    Lanczos3T lanczos3_;
    Lanczos4T lanczos4_;

// MinRes solver
    PMResSolverCL<Lanczos1T> PMinResMGBBT_;
    PMResSolverCL<Lanczos2T> PMinResMGMinComm_;
    PMResSolverCL<Lanczos3T> PMinResPCGBBT_;
    PMResSolverCL<Lanczos4T> PMinResPCGMinComm_;

  public:
    StokesSolverFactory(StokesT& Stokes, int outer_iter, double outer_tol, /*1/dt*/ double kA=1, /*theta*/ double kM=1)
        : Stokes_(Stokes), outer_iter_(outer_iter), outer_tol_(outer_tol), kA_(kA), kM_(kM), mgused_(false),
            bbtispc_( Stokes_.B.Data, Stokes_.prM.Data, Stokes_.M.Data, kA_, kM_, 1e-4, 1e-4),
            mincommispc_( 0, Stokes_.B.Data, Stokes_.M.Data, Stokes_.prM.Data, 1e-4),
            smoother_( 1.0), coarsesolver_( SSORPcCL(1.0), 500, 1e-16),
            MGSolver_ ( velMG_, smoother_, coarsesolver_, 2, -1.0, false), MGPc_( MGSolver_),
            GMResSolver_( JACPc_, 500, /*restart*/ 100, 1e-2, /*relative=*/ true), GMResPc_( GMResSolver_),
            PCGSolver_( SSORPc_, 500, 0.02, true), PCGPc_( PCGSolver_),
            DiagMGBBTOseenPc_        ( MGPc_,    bbtispc_), DiagMGMinCommOseenPc_    ( MGPc_,    mincommispc_),
            DiagPCGBBTOseenPc_       ( PCGPc_,   bbtispc_), DiagPCGMinCommOseenPc_   ( PCGPc_,   mincommispc_),
            LBlockMGBBTOseenPc_      ( MGPc_,    bbtispc_), LBlockMGMinCommOseenPc_   ( MGPc_,    mincommispc_),
            LBlockGMResBBTOseenPc_   ( GMResPc_, bbtispc_), LBlockGMResMinCommOseenPc_( GMResPc_, mincommispc_),
            GCRMGBBT_        ( LBlockMGBBTOseenPc_,        /*trunc*/ outer_iter, outer_iter, outer_tol, /*rel*/ false),
            GCRMGMinComm_    ( LBlockMGMinCommOseenPc_,    /*trunc*/ outer_iter, outer_iter, outer_tol, /*rel*/ false),
            GCRGMResBBT_     ( LBlockGMResBBTOseenPc_,     /*trunc*/ outer_iter, outer_iter, outer_tol, /*rel*/ false),
            GCRGMResMinComm_ ( LBlockGMResMinCommOseenPc_, /*trunc*/ outer_iter, outer_iter, outer_tol, /*rel*/ false),
            lanczos1_ (DiagMGBBTOseenPc_),  lanczos2_ (DiagMGMinCommOseenPc_),
            lanczos3_ (DiagPCGBBTOseenPc_), lanczos4_ (DiagPCGMinCommOseenPc_),
            PMinResMGBBT_     ( lanczos1_, outer_iter, outer_tol, /*relative*/ false),
            PMinResMGMinComm_ ( lanczos2_, outer_iter, outer_tol, /*relative*/ false),
            PMinResPCGBBT_    ( lanczos3_, outer_iter, outer_tol, /*relative*/ false),
            PMinResPCGMinComm_( lanczos4_, outer_iter, outer_tol, /*relative*/ false)
            {}

    ~StokesSolverFactory() {}

    MGDataCL&  GetVelMG() {return velMG_;}
    void       SetMatrixA (const MatrixCL* A) {mincommispc_.SetMatrixA(A);}
    bool       MGUsed()   {return mgused_;}

    StokesSolverBaseCL* CreateStokesSolver(int StokesMethod)
    {
        StokesSolverBaseCL* stokessolver = 0;
        switch (StokesMethod)
        {
            case 11 :
                stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockMGBBTOseenPcT> >       ( GCRMGBBT_);
            break;
            case 12 :
                stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockMGMinCommOseenPcT> >   ( GCRMGMinComm_);
            break;
            case 13 :
                stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockGMResBBTOseenPcT> >    ( GCRGMResBBT_);
            break;
            case 14 :
                stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockGMResMinCommOseenPcT> >( GCRGMResMinComm_);
            break;
            case 211 :
                stokessolver = new InexactUzawaCL<MGPcT, ISBBTPreCL,      APC_SYM>
                            ( MGPc_, bbtispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 222 :
                stokessolver = new InexactUzawaCL<MGPcT, MinCommPreCL,    APC_SYM>
                            ( MGPc_, mincommispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 21 :
                stokessolver = new InexactUzawaCL<MGPcT, ISBBTPreCL,      APC_OTHER>
                            ( MGPc_, bbtispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 22 :
                stokessolver = new InexactUzawaCL<MGPcT, MinCommPreCL,    APC_OTHER>
                            ( MGPc_, mincommispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 23 :
                stokessolver = new InexactUzawaCL<GMResPcT, ISBBTPreCL,   APC_OTHER>
                            ( GMResPc_, bbtispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 24 :
                stokessolver = new InexactUzawaCL<GMResPcT, MinCommPreCL, APC_OTHER>
                            ( GMResPc_, mincommispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 25 :
                stokessolver = new InexactUzawaCL<PCGPcT, ISBBTPreCL,     APC_SYM>
                            ( PCGPc_, bbtispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 26 :
                stokessolver = new InexactUzawaCL<PCGPcT, MinCommPreCL,   APC_SYM>
                            ( PCGPc_, mincommispc_, outer_iter_, outer_tol_, 0.6);
            break;
            case 31 :
                stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos1T> >( PMinResMGBBT_);
            break;
            case 32 :
                stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos2T> >( PMinResMGMinComm_);
            break;
            case 33 :
                stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos3T> >( PMinResPCGBBT_);
            break;
            case 34 :
                stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos4T> >( PMinResPCGMinComm_);
            break;
            default: throw DROPSErrCL("Unknown StokesMethod");
        }
        mgused_ = (C.StokesMethod == 11 || C.StokesMethod == 12 || C.StokesMethod == 21 ||
                   C.StokesMethod == 22 || C.StokesMethod == 31 || C.StokesMethod == 32 ||
                   C.StokesMethod == 211 || C.StokesMethod == 222);
        return stokessolver;
    }
};

template< class StokesProblemT>
TimeDisc2PhaseCL<StokesProblemT>* CreateTimeDisc(StokesProblemT& Stokes, LevelsetP2CL& lset,
    NSSolverBaseCL<StokesProblemT>* solver, ParamMesszelleNsCL& C, bool usematMG, MGDataCL* matMG)
{
    switch (C.scheme)
    {
        case 1 : 
            return (new LinThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_stab, usematMG, matMG));
        break;
        case 2 :
            return (new RecThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_proj, C.cpl_stab, usematMG, matMG));
        break;
        case 3 :
            return (new ThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_proj, C.cpl_stab, usematMG, matMG));
        break;
        case 4 :
            return (new OperatorSplitting2PhaseCL<StokesProblemT, StokesSolverBaseCL>
                        (Stokes, lset, solver->GetStokesSolver(), C.inner_iter, C.inner_tol, C.nonlinear));
        break;
        case 5 :
            return (new CrankNicolsonScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.nonlinear, C.cpl_proj, C.cpl_stab, usematMG, matMG));
        break;
        default : throw DROPSErrCL("Unknown TimeDiscMethod");
    }
}

class FunAsP2EvalCL
{
  private:
    instat_scalar_fun_ptr f_;
  public:
    FunAsP2EvalCL( instat_scalar_fun_ptr f): f_(f)
    {}
    double val(const VertexCL& v) const {return f_(v.GetCoord(), 0.0);}
    double val(const EdgeCL& e, double) const {return f_(GetBaryCenter(e), 0.0);}
};


class EnsightWriterCL
{
  private:
    MultiGridCL& MG_;
    EnsightP2SolOutCL ensight_;
    std::string datgeo_,
                datpr_,
                datvec_,
                datscl_,
                datprx_,
                datsf_,
                datc_,
                datct_;
  public:
    EnsightWriterCL (MultiGridCL& MG, const IdxDescCL* idx, const ParamMesszelleNsCL& C);
    ~EnsightWriterCL();

    // To write at time t.
    template<class InstatNSCL>
    void
    WriteAtTime (const InstatNSCL& Stokes, const LevelsetP2CL& lset,
        instat_scalar_fun_ptr sigmap, const TransportP1CL& c, const double t);
};

EnsightWriterCL::EnsightWriterCL (MultiGridCL& MG, const IdxDescCL* idx, const ParamMesszelleNsCL& C)
    : MG_( MG), ensight_( MG, idx)
{
    if (C.EnsCase == "none") return; // no ensight output
    const std::string filename= C.EnsDir + "/" + C.EnsCase;
    datgeo_= filename+".geo";
    datpr_ = filename+".pr" ;
    datvec_= filename+".vel";
    datscl_= filename+".scl";
    datprx_= filename+".prx";
    datsf_ = filename+".sf";
    datc_  = filename+".c";
    datct_ = filename+".ct";
    ensight_.CaseBegin( std::string( C.EnsCase+".case").c_str(), C.num_steps + 1);
    ensight_.DescribeGeom  ( "Messzelle",     datgeo_, true);
    ensight_.DescribeScalar( "Levelset",      datscl_, true);
    ensight_.DescribeScalar( "Pressure",      datpr_,  true);
    ensight_.DescribeVector( "Velocity",      datvec_, true);
    ensight_.DescribeScalar( "Surfaceforce",  datsf_,  true);
    if (C.transp_do)
    {
        ensight_.DescribeScalar( "Concentration", datc_,   true);
        ensight_.DescribeScalar( "TransConc",     datct_, true);
    }
}

EnsightWriterCL::~EnsightWriterCL()
{
    if (C.EnsCase == "none") return;
    ensight_.CaseEnd();
}

template<class InstatNSCL>
void
EnsightWriterCL::WriteAtTime (const InstatNSCL& Stokes, const LevelsetP2CL& lset,
    instat_scalar_fun_ptr sigmap, const TransportP1CL& c, const double t)
{
    if (C.EnsCase == "none") return;
    ensight_.putGeom( datgeo_, t);
    ensight_.putVector( datvec_, Stokes.GetVelSolution(), t);
    ensight_.putScalar( datpr_,  Stokes.GetPrSolution(), t);
    ensight_.putScalar( datscl_, lset.GetSolution(), t);
    FunAsP2EvalCL sf( sigmap);
    ensight_.putScalar( datsf_,  sf, t);
    if (C.transp_do)
    {
        ensight_.putScalar( datc_,  c.GetSolution(), t);
        ensight_.putScalar( datct_, c.GetSolution( c.ct), t);
    }
    ensight_.Commit();
    if (Stokes.UsesXFEM()) {
        std::string datprxnow( datprx_);
        ensight_.AppendTimecode( datprxnow);
        std::ofstream fff( datprxnow.c_str());
        fff.precision( 16);
        size_t num_prx= Stokes.pr_idx.NumUnknowns - Stokes.GetXidx().GetNumUnknownsP1();
        out( fff, VectorCL( Stokes.p.Data[std::slice( Stokes.GetXidx().GetNumUnknownsP1(), num_prx, 1)]));
    }
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
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    if (C.st_var)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);

    EnsightWriterCL writer( MG, lidx, C);

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
        Stokes.UpdateXNumbering( pidx, lset, /*NumberingChanged*/ false);
        Stokes.p.SetIdx( pidx);
        // reads the P1-part of the pressure
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        // reads the P1X-part of the pressure
        if (Stokes.UsesXFEM()) {
            std::ifstream fff( (C.IniData+".prx").c_str());
            if (fff) {
                size_t NumP1XUnknowns;
                fff >> NumP1XUnknowns;
                if (NumP1XUnknowns != (pidx->NumUnknowns - Stokes.GetXidx().GetNumUnknownsP1()))
                    throw (DROPSErrCL("error while reading P1X unknowns"));
                for (Uint i=Stokes.GetXidx().GetNumUnknownsP1(); i < pidx->NumUnknowns; ++i)
                    fff >> Stokes.p.Data[i];
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
        IdxDescCL* cidx= &c.idx;
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

    writer.WriteAtTime( Stokes, lset, sigmap, c, 0.);

    // Stokes-Solver
    StokesSolverFactory<StokesProblemT> stokessolverfactory(Stokes, C.outer_iter, C.outer_tol, 1.0/C.dt, C.theta);
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver(C.StokesMethod);

    // Navier-Stokes-Solver
    NSSolverBaseCL<StokesProblemT>* navstokessolver = 0;
    if (C.nonlinear==0.0)
        navstokessolver = new NSSolverBaseCL<StokesProblemT>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<StokesProblemT>(Stokes, *stokessolver, C.ns_iter, C.ns_tol, C.ns_red);

    // Time discretisation + coupling
    MGDataCL& velMG = stokessolverfactory.GetVelMG();
    bool mgused = stokessolverfactory.MGUsed();
    TimeDisc2PhaseCL<StokesProblemT>* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, C, mgused, &velMG);
    timedisc->SetTimeStep( C.dt);

    stokessolverfactory.SetMatrixA( timedisc->GetUpperLeftBlock());
    bool second = false;
    std::ofstream infofile((C.EnsCase+".info").c_str());
    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== step " << step << ":\n";

        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(infofile);
        timedisc->DoStep( C.cpl_iter);
        if (C.transp_do) c.DoStep( step*C.dt);

//        WriteMatrices( Stokes, step);
        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        if (C.VolCorr) {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        }
        if (C.RepFreq && step%C.RepFreq==0) { // reparam levelset function
            lset.ReparamFastMarching( C.RepMethod);
            if (C.ref_freq != 0)
                adap.UpdateTriang( lset);
            if (adap.WasModified()) {
                timedisc->Update();
                if (C.transp_do) c.Update();
            }
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
            if (C.serialization_file != "none" && C.ref_freq != 0) {
                std::stringstream filename;
                filename << C.serialization_file;
                if (second) filename << "0";
                second = !second;
                MGSerializationCL ser( MG, filename.str().c_str());
                ser.WriteMG();
            }
        }
        writer.WriteAtTime( Stokes, lset, sigmap, c, step*C.dt);
    }
    std::cerr << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
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
        DROPS::StokesVelBndDataCL::bnd_val_fun* bnd_fun = new  DROPS::StokesVelBndDataCL::bnd_val_fun[num_bnd];
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

    EllipsoidCL::Init( C.Mitte, C.Radius );

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
