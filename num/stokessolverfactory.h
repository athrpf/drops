/// \file
/// \brief creates several standard Stokes-solver
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

#include "num/stokessolver.h"

/* Available stokes solver
 *   no | Solver         | APc       | SPc     | PC structure
 ------------------------------------------------------------
 *   11 | GCR            | MG        | BBT     | lower block
 *   12 | GCR            | MG        | MinComm | lower block
 *   13 | GCR            | GMRes     | BBT     | lower block
 *   14 | GCR            | GMRes     | MinComm | lower block
 *   15 | GCR            | BiCGStab  | BBT     | lower block
 *   16 | GCR            | BiCGStab  | MinComm | lower block
------------------------------------------------------------
 *   21 | Inexact Uzawa  | asymm. MG | BBT     |
 *   22 | Inexact Uzawa  | asymm. MG | MinComm |
 *  221 | Inexact Uzawa  | symm. MG  | BBT     |
 *  222 | Inexact Uzawa  | symm. MG  | MinComm |
 *   23 | Inexact Uzawa  | GMRes     | BBT     |
 *   24 | Inexact Uzawa  | GMRes     | MinComm |
 *   25 | Inexact Uzawa  | BiCGStab  | BBT     |
 *   26 | Inexact Uzawa  | BiCGStab  | MinComm |
 *  225 | Inexact Uzawa  | SSORPCG   | BBT     |
 *  226 | Inexact Uzawa  | SSORPCG   | MinComm |
 ------------------------------------------------------------
 *   31 | PMinRes        | MG        | BBT     | diag
 *   32 | PMinRes        | MG        | MinComm | diag
 *   33 | PMinRes        | SSORPCG   | BBT     | diag
 *   34 | PMinRes        | SSORPCG   | MinComm | diag
 ------------------------------------------------------------
 *   41 | GMRes          | MG        | BBT     | lower block
 *   42 | GMRes          | MG        | MinComm | lower block
 *   43 | GMRes          | GMRes     | BBT     | lower block
 *   44 | GMRes          | GMRes     | MinComm | lower block
 *   45 | GMRes          | BiCGStab  | BBT     | lower block
 *   46 | GMRes          | BiCGStab  | MinComm | lower block
 ------------------------------------------------------------
 *   51 | GMResR         | MG        | BBT     | lower block
 *   52 | GMResR         | MG        | MinComm | lower block
 *   53 | GMResR         | GMRes     | BBT     | lower block
 *   54 | GMResR         | GMRes     | MinComm | lower block
 *   55 | GMResR         | BiCGStab  | BBT     | lower block
 *   56 | GMResR         | BiCGStab  | MinComm | lower block
 ------------------------------------------------------------
 *   81 | StokesMGM      | PVanka-Smoother
 *   82 | StokesMGM      | Braess/Sarazin-Smoother
*/

namespace DROPS {

template <class StokesT, class ParamsT>
class StokesSolverFactoryCL
{
  private:
    StokesT& Stokes_;
    ParamsT& C_;
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

    //JAC-BiCGStab
    typedef BiCGStabSolverCL<JACPcCL> BiCGStabSolverT;        // BiCGStab-based APcT
    BiCGStabSolverT BiCGStabSolver_;
    typedef SolverAsPreCL<BiCGStabSolverT> BiCGStabPcT;
    BiCGStabPcT BiCGStabPc_;

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
    typedef BlockPreCL<BiCGStabPcT, ISBBTPreCL,   LowerBlockPreCL> LBlockBiCGBBTOseenPcT;
    typedef BlockPreCL<BiCGStabPcT, MinCommPreCL, LowerBlockPreCL> LBlockBiCGMinCommOseenPcT;

    DiagMGBBTOseenPcT         DiagMGBBTOseenPc_;
    DiagMGMinCommOseenPcT     DiagMGMinCommOseenPc_;
    DiagPCGBBTOseenPcT        DiagPCGBBTOseenPc_;
    DiagPCGMinCommOseenPcT    DiagPCGMinCommOseenPc_;

    LBlockMGBBTOseenPcT        LBlockMGBBTOseenPc_;
    LBlockMGMinCommOseenPcT    LBlockMGMinCommOseenPc_;
    LBlockGMResBBTOseenPcT     LBlockGMResBBTOseenPc_;
    LBlockGMResMinCommOseenPcT LBlockGMResMinCommOseenPc_;
    LBlockBiCGBBTOseenPcT      LBlockBiCGBBTOseenPc_;
    LBlockBiCGMinCommOseenPcT  LBlockBiCGMinCommOseenPc_;

//GCR solver
    GCRSolverCL<LBlockMGBBTOseenPcT>        GCRMGBBT_;
    GCRSolverCL<LBlockMGMinCommOseenPcT>    GCRMGMinComm_;
    GCRSolverCL<LBlockGMResBBTOseenPcT>     GCRGMResBBT_;
    GCRSolverCL<LBlockGMResMinCommOseenPcT> GCRGMResMinComm_;
    GCRSolverCL<LBlockBiCGBBTOseenPcT>      GCRBiCGStabBBT_;
    GCRSolverCL<LBlockBiCGMinCommOseenPcT>  GCRBiCGStabMinComm_;

//GMRes solver
    GMResSolverCL<LBlockMGBBTOseenPcT>        GMResMGBBT_;
    GMResSolverCL<LBlockMGMinCommOseenPcT>    GMResMGMinComm_;
    GMResSolverCL<LBlockGMResBBTOseenPcT>     GMResGMResBBT_;
    GMResSolverCL<LBlockGMResMinCommOseenPcT> GMResGMResMinComm_;
    GMResSolverCL<LBlockBiCGBBTOseenPcT>      GMResBiCGStabBBT_;
    GMResSolverCL<LBlockBiCGMinCommOseenPcT>  GMResBiCGStabMinComm_;

// GMResR solver
    GMResRSolverCL<LBlockMGBBTOseenPcT>        GMResRMGBBT_;
    GMResRSolverCL<LBlockMGMinCommOseenPcT>    GMResRMGMinComm_;
    GMResRSolverCL<LBlockGMResBBTOseenPcT>     GMResRGMResBBT_;
    GMResRSolverCL<LBlockGMResMinCommOseenPcT> GMResRGMResMinComm_;
    GMResRSolverCL<LBlockBiCGBBTOseenPcT>      GMResRBiCGStabBBT_;
    GMResRSolverCL<LBlockBiCGMinCommOseenPcT>  GMResRBiCGStabMinComm_;

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

    //MinRes solver
    LanczosONBCL<BlockMatrixCL, VectorCL> q;
    PMResSolverCL<LanczosONBCL<BlockMatrixCL, VectorCL> > minressolver;
    BlockMatrixSolverCL<PMResSolverCL<LanczosONBCL<BlockMatrixCL, VectorCL> > > blockminressolver;

    PVankaSmootherCL vankasmoother;
    BSSmootherCL bssmoother;

  public:
    StokesSolverFactoryCL(StokesT& Stokes, ParamsT& C);
    ~StokesSolverFactoryCL() {}

    MGDataCL&  GetVelMG() {return velMG_;}
    void       SetMatrixA (const MatrixCL* A) {mincommispc_.SetMatrixA(A);}
    bool       MGUsed()   {return mgused_;}

    StokesSolverBaseCL* CreateStokesSolver();
};

template <class StokesT, class ParamsT>
StokesSolverFactoryCL<StokesT, ParamsT>::StokesSolverFactoryCL(StokesT& Stokes, ParamsT& C)
    : Stokes_(Stokes), C_(C), kA_(1.0/C_.dt), kM_(C_.theta), mgused_(false),
        bbtispc_( Stokes_.B.Data, Stokes_.prM.Data, Stokes_.M.Data, kA_, kM_, C_.pcS_tol, C_.pcS_tol),
        mincommispc_( 0, Stokes_.B.Data, Stokes_.M.Data, Stokes_.prM.Data, C_.pcS_tol),
        smoother_( 1.0), coarsesolver_( SSORPcCL(1.0), 500, 1e-16),
        MGSolver_ ( velMG_, smoother_, coarsesolver_, C_.pcA_iter, C_.pcA_tol, false), MGPc_( MGSolver_),
        GMResSolver_( JACPc_, C_.pcA_iter, /*restart*/ 100, C_.pcA_tol, /*rel*/ true), GMResPc_( GMResSolver_),
        BiCGStabSolver_( JACPc_, C_.pcA_iter, C_.pcA_tol, /*rel*/ true),BiCGStabPc_( BiCGStabSolver_),
        PCGSolver_( SSORPc_, C_.pcA_iter, C_.pcA_tol, true), PCGPc_( PCGSolver_),
        DiagMGBBTOseenPc_        ( MGPc_,    bbtispc_), DiagMGMinCommOseenPc_     ( MGPc_,    mincommispc_),
        DiagPCGBBTOseenPc_       ( PCGPc_,   bbtispc_), DiagPCGMinCommOseenPc_    ( PCGPc_,   mincommispc_),
        LBlockMGBBTOseenPc_      ( MGPc_,    bbtispc_), LBlockMGMinCommOseenPc_   ( MGPc_,    mincommispc_),
        LBlockGMResBBTOseenPc_   ( GMResPc_, bbtispc_), LBlockGMResMinCommOseenPc_( GMResPc_, mincommispc_),
        LBlockBiCGBBTOseenPc_    ( BiCGStabPc_, bbtispc_), LBlockBiCGMinCommOseenPc_( BiCGStabPc_, mincommispc_),
        GCRMGBBT_           ( LBlockMGBBTOseenPc_,        C_.outer_iter, C_.outer_iter, C_.outer_tol, /*rel*/ false),
        GCRMGMinComm_       ( LBlockMGMinCommOseenPc_,    C_.outer_iter, C_.outer_iter, C_.outer_tol, /*rel*/ false),
        GCRGMResBBT_        ( LBlockGMResBBTOseenPc_,     C_.outer_iter, C_.outer_iter, C_.outer_tol, /*rel*/ false),
        GCRGMResMinComm_    ( LBlockGMResMinCommOseenPc_, C_.outer_iter, C_.outer_iter, C_.outer_tol, /*rel*/ false),
        GCRBiCGStabBBT_     ( LBlockBiCGBBTOseenPc_,      C_.outer_iter, C_.outer_iter, C_.outer_tol, /*rel*/ false),
        GCRBiCGStabMinComm_ ( LBlockBiCGMinCommOseenPc_,  C_.outer_iter, C_.outer_iter, C_.outer_tol, /*rel*/ false),
        GMResMGBBT_           ( LBlockMGBBTOseenPc_,        C_.outer_iter, C_.outer_iter, C_.outer_tol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResMGMinComm_       ( LBlockMGMinCommOseenPc_,    C_.outer_iter, C_.outer_iter, C_.outer_tol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResGMResBBT_        ( LBlockGMResBBTOseenPc_,     C_.outer_iter, C_.outer_iter, C_.outer_tol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResGMResMinComm_    ( LBlockGMResMinCommOseenPc_, C_.outer_iter, C_.outer_iter, C_.outer_tol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResBiCGStabBBT_     ( LBlockBiCGBBTOseenPc_,      C_.outer_iter, C_.outer_iter, C_.outer_tol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResBiCGStabMinComm_ ( LBlockBiCGMinCommOseenPc_,  C_.outer_iter, C_.outer_iter, C_.outer_tol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResRMGBBT_           ( LBlockMGBBTOseenPc_,        C_.outer_iter, C_.outer_iter, C_.inner_iter,
                                    C_.outer_tol, C_.inner_tol, /*rel*/ false),
        GMResRMGMinComm_       ( LBlockMGMinCommOseenPc_,    C_.outer_iter, C_.outer_iter, C_.inner_iter,
                                    C_.outer_tol, C_.inner_tol, /*rel*/ false),
        GMResRGMResBBT_        ( LBlockGMResBBTOseenPc_,     C_.outer_iter, C_.outer_iter, C_.inner_iter,
                                    C_.outer_tol, C_.inner_tol, /*rel*/ false),
        GMResRGMResMinComm_    ( LBlockGMResMinCommOseenPc_, C_.outer_iter, C_.outer_iter, C_.inner_iter,
                                    C_.outer_tol, C_.inner_tol, /*rel*/ false),
        GMResRBiCGStabBBT_     ( LBlockBiCGBBTOseenPc_,      C_.outer_iter, C_.outer_iter, C_.inner_iter,
                                    C_.outer_tol, C_.inner_tol, /*rel*/ false),
        GMResRBiCGStabMinComm_ ( LBlockBiCGMinCommOseenPc_,  C_.outer_iter, C_.outer_iter, C_.inner_iter,
                                    C_.outer_tol, C_.inner_tol, /*rel*/ false),
        lanczos1_ (DiagMGBBTOseenPc_),  lanczos2_ (DiagMGMinCommOseenPc_),
        lanczos3_ (DiagPCGBBTOseenPc_), lanczos4_ (DiagPCGMinCommOseenPc_),
        PMinResMGBBT_     ( lanczos1_, C_.outer_iter, C_.outer_tol, /*relative*/ false),
        PMinResMGMinComm_ ( lanczos2_, C_.outer_iter, C_.outer_tol, /*relative*/ false),
        PMinResPCGBBT_    ( lanczos3_, C_.outer_iter, C_.outer_tol, /*relative*/ false),
        PMinResPCGMinComm_( lanczos4_, C_.outer_iter, C_.outer_tol, /*relative*/ false),
        minressolver( q, C_.inner_iter, C_.inner_tol), blockminressolver(minressolver)
        {}

template <class StokesT, class ParamsT>
StokesSolverBaseCL* StokesSolverFactoryCL<StokesT, ParamsT>::CreateStokesSolver()
{
    StokesSolverBaseCL* stokessolver = 0;
    switch (C_.StokesMethod)
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
        case 15 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockBiCGBBTOseenPcT> >     ( GCRBiCGStabBBT_);
        break;
        case 16 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockBiCGMinCommOseenPcT> > ( GCRBiCGStabMinComm_);
        break;
        case 211 :
            stokessolver = new InexactUzawaCL<MGPcT, ISBBTPreCL,      APC_SYM>
                        ( MGPc_, bbtispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 222 :
            stokessolver = new InexactUzawaCL<MGPcT, MinCommPreCL,    APC_SYM>
                        ( MGPc_, mincommispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 21 :
            stokessolver = new InexactUzawaCL<MGPcT, ISBBTPreCL,      APC_OTHER>
                        ( MGPc_, bbtispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 22 :
            stokessolver = new InexactUzawaCL<MGPcT, MinCommPreCL,    APC_OTHER>
                        ( MGPc_, mincommispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 23 :
            stokessolver = new InexactUzawaCL<GMResPcT, ISBBTPreCL,   APC_OTHER>
                        ( GMResPc_, bbtispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 24 :
            stokessolver = new InexactUzawaCL<GMResPcT, MinCommPreCL, APC_OTHER>
                        ( GMResPc_, mincommispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 225 :
            stokessolver = new InexactUzawaCL<PCGPcT, ISBBTPreCL,     APC_SYM>
                        ( PCGPc_, bbtispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 226 :
            stokessolver = new InexactUzawaCL<PCGPcT, MinCommPreCL,   APC_SYM>
                        ( PCGPc_, mincommispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 25 :
            stokessolver = new InexactUzawaCL<BiCGStabPcT, ISBBTPreCL,     APC_SYM>
                        ( BiCGStabPc_, bbtispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
        break;
        case 26 :
            stokessolver = new InexactUzawaCL<BiCGStabPcT, MinCommPreCL,   APC_SYM>
                        ( BiCGStabPc_, mincommispc_, C_.outer_iter, C_.outer_tol, C_.inner_tol);
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
        case 41 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockMGBBTOseenPcT> >       ( GMResMGBBT_);
        break;
        case 42 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockMGMinCommOseenPcT> >   ( GMResMGMinComm_);
        break;
        case 43 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockGMResBBTOseenPcT> >    ( GMResGMResBBT_);
        break;
        case 44 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockGMResMinCommOseenPcT> >( GMResGMResMinComm_);
        break;
        case 45 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockBiCGBBTOseenPcT> >     ( GMResBiCGStabBBT_);
        break;
        case 46 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockBiCGMinCommOseenPcT> > ( GMResBiCGStabMinComm_);
        break;
        case 51 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockMGBBTOseenPcT> >       ( GMResRMGBBT_);
        break;
        case 52 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockMGMinCommOseenPcT> >   ( GMResRMGMinComm_);
        break;
        case 53 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockGMResBBTOseenPcT> >    ( GMResRGMResBBT_);
        break;
        case 54 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockGMResMinCommOseenPcT> >( GMResRGMResMinComm_);
        break;
        case 55 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockBiCGBBTOseenPcT> >     ( GMResRBiCGStabBBT_);
        break;
        case 56 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockBiCGMinCommOseenPcT> > ( GMResRBiCGStabMinComm_);
        break;
        case 81 : {
            if (C_.XFEMStab >= 0) // P1X
                throw DROPSErrCL("StokesMGM not implemented for P1X-elements");
            velMG_.SetStokesMG(true);
            vankasmoother.SetVankaMethod(2);
            stokessolver = new StokesMGSolverCL<PVankaSmootherCL>( velMG_, vankasmoother, blockminressolver, C_.outer_iter, C_.outer_tol, false, 1);
        }
        break;
        case 82 : {
            if (C_.XFEMStab >= 0) // P1X
                throw DROPSErrCL("StokesMGM not implemented for P1X-elements");
            velMG_.SetStokesMG(true);
            stokessolver = new StokesMGSolverCL<BSSmootherCL>( velMG_, bssmoother, blockminressolver, C_.outer_iter, C_.outer_tol, false, 4);
        }
        break;
        default: throw DROPSErrCL("Unknown StokesMethod");
    }
    mgused_ = (C_.StokesMethod == 11 || C_.StokesMethod == 12 || C_.StokesMethod == 21 ||
               C_.StokesMethod == 22 || C_.StokesMethod == 31 || C_.StokesMethod == 32 ||
               C_.StokesMethod == 41 || C_.StokesMethod == 42 || C_.StokesMethod == 51 ||
               C_.StokesMethod == 52 || C_.StokesMethod == 211 || C_.StokesMethod == 222 || C_.StokesMethod == 81  || C_.StokesMethod == 82);
    return stokessolver;
}

} // end of namespace DROPS
