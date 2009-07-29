/// \file
/// \brief creates several standard Stokes-solver
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

#include "num/stokessolver.h"

namespace DROPS {

/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y B a s e   C L            *
*******************************************************************/
/// \brief Creates a StokesSolverCL* and manages the preconditioner.
/// Interface for all Stokes solver factories.
/// Construction of an Oseen solver, e.g. Inexact Uzawa with GMRes and BBT preconditioner: 2*10000 + 4*100 + 2. Note: Not all combinations are implemented!
/**
    <table border="3">
    <tr><th> no </th><th> Oseen-Solver-Type </th><th> Type of Preconditioner for A-Block </th><th> Type of Preconditioner for S </th></tr>
    <tr><td>  1 </td><td> GCR               </td><td> MultiGrid V-cycle                  </td><td> ISBBTPreCL                   </td></tr>
    <tr><td>  2 </td><td> Inexact Uzawa     </td><td> symm. Multigrid V-cycle            </td><td> MinCommPreCL                 </td></tr>
    <tr><td>  3 </td><td> PMinRes           </td><td> PCG                                </td><td> ISPreCL                      </td></tr>
    <tr><td>  4 </td><td> GMRes             </td><td> GMRes                              </td><td> VankaSchurPreCL              </td></tr>
    <tr><td>  5 </td><td> GMResR            </td><td> BiCGStab                           </td><td>                              </td></tr>
    <tr><td>  6 </td><td>                   </td><td> VankaPre                           </td><td> VankaPre                     </td></tr>
    <tr><td>  7 </td><td>                   </td><td>                                    </td><td> ISMGPreCL                    </td></tr>
    <tr><td>  8 </td><td>                   </td><td>                                    </td><td>                              </td></tr>
    <tr><td>  9 </td><td>                   </td><td>                                    </td><td>                              </td></tr>
    <tr><td> 30 </td><td> StokesMGM         </td><td> PVankaSmootherCL                   </td><td> PVankaSmootherCL             </td></tr>
    <tr><td> 31 </td><td>                   </td><td> BSSmootherCL                       </td><td> BSSmootherCL                 </td></tr>
    </table>*/

template <class StokesT, class ParamsT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryBaseCL
{
  protected:
    StokesT& Stokes_;           ///< Stokes problem
    ParamsT& C_;                ///< Parameter for tolerances, iteration number, type of solver, ...
    int      SPc_,              ///< type of preconditioner for S
             APc_,              ///< type of preconditioner for A-block
             OseenSolver_;      ///< type of Oseen solver

  public:
    StokesSolverFactoryBaseCL( StokesT& Stokes, ParamsT& C) : Stokes_( Stokes), C_( C),
                               SPc_( C_.stk_StokesMethod % 100), APc_( (C_.stk_StokesMethod / 100) % 100),
                               OseenSolver_( (C_.stk_StokesMethod /10000) % 100) {}
    virtual ~StokesSolverFactoryBaseCL() {}

    /// Set the A-block in the minimal commutator
    virtual void       SetMatrixA ( const MatrixCL*) = 0;
    /// Set all matrices in Schur complement preconditioner
    virtual void       SetMatrices( const MatrixCL*, const MatrixCL*, const MatrixCL*, const MatrixCL*, const IdxDescCL* pr_idx) = 0;
    /// Returns pointer to prolongation for velocity
    virtual ProlongationVelT* GetPVel() = 0;
    /// Returns pointer to prolongation for pressure
    virtual ProlongationPT*   GetPPr() = 0;
    /// Returns a stokes solver with specifications from ParamsT C
    virtual StokesSolverBaseCL* CreateStokesSolver() = 0;
};

/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y H e l p e r  C L         *
********************************************************************/
template <class ParamsT>
class StokesSolverFactoryHelperCL
{
  public:
    bool VelMGUsed ( const ParamsT& C) const
    {
        const int APc = (C.stk_StokesMethod / 100) % 100;
        return (( APc == 1) || (APc == 2) || (APc == 30) || (APc == 31));
    }
    bool PrMGUsed  ( const ParamsT& C) const
    {
        const int APc = (C.stk_StokesMethod / 100) % 100,
            SPc = C.stk_StokesMethod % 100;
        return (( APc == 30) || ( APc == 31) || (SPc == 7));
    }
};

#ifndef _PAR
/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y  C L                     *
********************************************************************/
template <class StokesT, class ParamsT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryCL : public StokesSolverFactoryBaseCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>
{
  private:
    typedef StokesSolverFactoryBaseCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT> base_;
    using base_::Stokes_;
    using base_::C_;
    using base_::OseenSolver_;
    using base_::APc_;
    using base_::SPc_;
    double kA_, kM_;

// generic preconditioners
    JACPcCL  JACPc_;
    SSORPcCL SSORPc_;

// PC for instat. Schur complement
    ISBBTPreCL      bbtispc_;
    MinCommPreCL    mincommispc_;
    VankaSchurPreCL vankaschurpc_;
    ISPreCL         isprepc_;
    ISMGPreCL       ismgpre_;

// PC for A-block
    // MultiGrid symm.
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolversymm_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL, ProlongationVelT> MGSolversymm_;
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL, ProlongationVelT> > MGsymmPcT;
    MGsymmPcT MGPcsymm_;

    // Multigrid nonsymm.
    GMResSolverCL<JACPcCL> coarsesolver_;
    MGSolverCL<SSORsmoothCL, GMResSolverCL<JACPcCL>, ProlongationVelT > MGSolver_;
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, GMResSolverCL<JACPcCL>, ProlongationVelT> > MGPcT;
    MGPcT MGPc_;

    //JAC-GMRes
    typedef GMResSolverCL<JACPcCL> GMResSolverT;
    GMResSolverT GMResSolver_;
    typedef SolverAsPreCL<GMResSolverT> GMResPcT;
    GMResPcT GMResPc_;

    //JAC-BiCGStab
    typedef BiCGStabSolverCL<JACPcCL> BiCGStabSolverT;
    BiCGStabSolverT BiCGStabSolver_;
    typedef SolverAsPreCL<BiCGStabSolverT> BiCGStabPcT;
    BiCGStabPcT BiCGStabPc_;

    //PCG
    typedef PCGSolverCL<SSORPcCL> PCGSolverT;
    PCGSolverT PCGSolver_;
    typedef SolverAsPreCL<PCGSolverT> PCGPcT;
    PCGPcT PCGPc_;

// PC for Oseen problem
    typedef BlockPreCL<MGsymmPcT, ISBBTPreCL>   DiagMGBBTOseenPcT;
    typedef BlockPreCL<MGsymmPcT, MinCommPreCL> DiagMGMinCommOseenPcT;
    typedef BlockPreCL<MGsymmPcT, ISPreCL>      DiagMGISPreOseenPcT;
    typedef BlockPreCL<GMResPcT,  MinCommPreCL> DiagGMResMinCommOseenPcT;
    typedef BlockPreCL<PCGPcT,    ISBBTPreCL>   DiagPCGBBTOseenPcT;
    typedef BlockPreCL<PCGPcT,    MinCommPreCL> DiagPCGMinCommOseenPcT;
    typedef BlockPreCL<PCGPcT,    ISPreCL>      DiagPCGISPreOseenPcT;

    typedef BlockPreCL<MGPcT,       ISBBTPreCL,   LowerBlockPreCL> LBlockMGBBTOseenPcT;
    typedef BlockPreCL<MGPcT,       MinCommPreCL, LowerBlockPreCL> LBlockMGMinCommOseenPcT;
    typedef BlockPreCL<MGPcT,       ISPreCL,      LowerBlockPreCL> LBlockMGISPreOseenPcT;
    typedef BlockPreCL<MGPcT,       VankaSchurPreCL, LowerBlockPreCL> LBlockMGVankaOseenPcT;
    typedef BlockPreCL<GMResPcT,    ISBBTPreCL,   LowerBlockPreCL> LBlockGMResBBTOseenPcT;
    typedef BlockPreCL<GMResPcT,    MinCommPreCL, LowerBlockPreCL> LBlockGMResMinCommOseenPcT;
    typedef BlockPreCL<BiCGStabPcT, ISBBTPreCL,   LowerBlockPreCL> LBlockBiCGBBTOseenPcT;
    typedef BlockPreCL<BiCGStabPcT, MinCommPreCL, LowerBlockPreCL> LBlockBiCGMinCommOseenPcT;

    DiagMGBBTOseenPcT         DiagMGBBTOseenPc_;
    DiagMGMinCommOseenPcT     DiagMGMinCommOseenPc_;
    DiagMGISPreOseenPcT       DiagMGISPreOseenPc_;
    DiagGMResMinCommOseenPcT  DiagGMResMinCommPc_;
    DiagPCGBBTOseenPcT        DiagPCGBBTOseenPc_;
    DiagPCGMinCommOseenPcT    DiagPCGMinCommOseenPc_;
    DiagPCGISPreOseenPcT      DiagPCGISPreOseenPc_;

    LBlockMGBBTOseenPcT        LBlockMGBBTOseenPc_;
    LBlockMGMinCommOseenPcT    LBlockMGMinCommOseenPc_;
    LBlockMGISPreOseenPcT      LBlockMGISPreOseenPc_;
    LBlockMGVankaOseenPcT      LBlockMGVankaOseenPc_;
    LBlockGMResBBTOseenPcT     LBlockGMResBBTOseenPc_;
    LBlockGMResMinCommOseenPcT LBlockGMResMinCommOseenPc_;
    LBlockBiCGBBTOseenPcT      LBlockBiCGBBTOseenPc_;
    LBlockBiCGMinCommOseenPcT  LBlockBiCGMinCommOseenPc_;

    VankaPreCL vankapc_;

//GCR solver
    GCRSolverCL<LBlockMGBBTOseenPcT>        GCRMGBBT_;
    GCRSolverCL<LBlockMGMinCommOseenPcT>    GCRMGMinComm_;
    GCRSolverCL<LBlockMGISPreOseenPcT>      GCRMGISPre_;
    GCRSolverCL<LBlockMGVankaOseenPcT>      GCRMGVanka_;
    GCRSolverCL<LBlockGMResBBTOseenPcT>     GCRGMResBBT_;
    GCRSolverCL<LBlockGMResMinCommOseenPcT> GCRGMResMinComm_;
    GCRSolverCL<LBlockBiCGBBTOseenPcT>      GCRBiCGStabBBT_;
    GCRSolverCL<LBlockBiCGMinCommOseenPcT>  GCRBiCGStabMinComm_;
    GCRSolverCL<VankaPreCL>                 GCRVanka_;

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
    typedef PLanczosONBCL<VectorCL, DiagMGBBTOseenPcT>      Lanczos1T;
    typedef PLanczosONBCL<VectorCL, DiagMGMinCommOseenPcT>  Lanczos2T;
    typedef PLanczosONBCL<VectorCL, DiagPCGBBTOseenPcT>     Lanczos3T;
    typedef PLanczosONBCL<VectorCL, DiagPCGMinCommOseenPcT> Lanczos4T;
    typedef PLanczosONBCL<VectorCL, DiagMGISPreOseenPcT>    Lanczos5T;
    typedef PLanczosONBCL<VectorCL, DiagPCGISPreOseenPcT>   Lanczos6T;

    Lanczos1T lanczos1_;
    Lanczos2T lanczos2_;
    Lanczos3T lanczos3_;
    Lanczos4T lanczos4_;
    Lanczos5T lanczos5_;
    Lanczos6T lanczos6_;


// MinRes solver
    PMResSolverCL<Lanczos1T> PMinResMGBBT_;
    PMResSolverCL<Lanczos2T> PMinResMGMinComm_;
    PMResSolverCL<Lanczos3T> PMinResPCGBBT_;
    PMResSolverCL<Lanczos4T> PMinResPCGMinComm_;
    PMResSolverCL<Lanczos5T> PMinResMGISPre_;
    PMResSolverCL<Lanczos6T> PMinResPCGISPre_;

// coarse grid solver
    //MinRes solver
    PMResSolverCL<Lanczos3T> minressolver;
    BlockMatrixSolverCL<PMResSolverCL<Lanczos3T> > blockminressolver;

    //GCR solver
    GCRSolverCL<DiagGMResMinCommOseenPcT> gcrsolver;
    BlockMatrixSolverCL<GCRSolverCL<DiagGMResMinCommOseenPcT> > blockgcrsolver;

    PVankaSmootherCL vankasmoother;
    BSSmootherCL bssmoother;

    //StokesMGSolver
    StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>* mgvankasolversymm_;
    StokesMGSolverCL<BSSmootherCL,     ProlongationVelT, ProlongationPT>* mgbssolversymm_;
    StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>* mgvankasolver_;
    StokesMGSolverCL<BSSmootherCL,     ProlongationVelT, ProlongationPT>* mgbssolver_;


  public:
    StokesSolverFactoryCL(StokesT& Stokes, ParamsT& C);
    ~StokesSolverFactoryCL() {}

    /// Set the A-block in the minimal commutator
    void       SetMatrixA ( const MatrixCL* A) { mincommispc_.SetMatrixA(A); }
    /// Set all matrices in Schur complement preconditioner (only for StokesMGM)
    void       SetMatrices( const MatrixCL* A, const MatrixCL* B, const MatrixCL* Mvel, const MatrixCL* M, const IdxDescCL* pr_idx);
    /// Returns pointer to prolongation for velocity
    ProlongationVelT* GetPVel();
    /// Returns pointer to prolongation for pressure
    ProlongationPT*   GetPPr();
    /// Returns a stokes solver with specifications from ParamsT C
    StokesSolverBaseCL* CreateStokesSolver();

    PVankaSmootherCL&      GetVankaSmoother () { return vankasmoother; }
    VankaSchurPreCL&       GetVankaSchurPc ()  { return vankaschurpc_; }
};

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::
    StokesSolverFactoryCL(StokesT& Stokes, ParamsT& C)
    : base_(Stokes, C),
        kA_(C_.tm_NumSteps != 0 ? 1.0/C_.tm_StepSize : 0.0), // C_.tm_NumSteps == 0: stat. problem
        kM_(C_.stk_Theta),
        // schur complement preconditioner
        bbtispc_    ( &Stokes_.B.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), Stokes_.pr_idx.GetFinest(), kA_, kM_, C_.stk_PcSTol, C_.stk_PcSTol /* enable regularization: , 0.707*/),
        mincommispc_( 0, &Stokes_.B.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(),Stokes_.pr_idx.GetFinest(), C_.stk_PcSTol /* enable regularization: , 0.707*/),
        vankaschurpc_( &Stokes.pr_idx), isprepc_( Stokes.prA.Data, Stokes.prM.Data, kA_, kM_),
        ismgpre_( Stokes.prA.Data, Stokes.prM.Data, kA_, kM_), 
        // preconditioner for A
        smoother_( 1.0), coarsesolversymm_( SSORPc_, 500, 1e-6, true),
        MGSolversymm_ ( smoother_, coarsesolversymm_, C_.stk_PcAIter, C_.stk_PcATol, false),
        MGPcsymm_( MGSolversymm_),
        coarsesolver_( JACPc_, 500, 500, 1e-6, true),
        MGSolver_ ( smoother_, coarsesolver_, C_.stk_PcAIter, C_.stk_PcATol, false), MGPc_( MGSolver_),
        GMResSolver_( JACPc_, C_.stk_PcAIter, /*restart*/ 100, C_.stk_PcATol, /*rel*/ true), GMResPc_( GMResSolver_),
        BiCGStabSolver_( JACPc_, C_.stk_PcAIter, C_.stk_PcATol, /*rel*/ true),BiCGStabPc_( BiCGStabSolver_),
        PCGSolver_( SSORPc_, C_.stk_PcAIter, C_.stk_PcATol, true), PCGPc_( PCGSolver_),
        // block precondtioner
        DiagMGBBTOseenPc_        ( MGPcsymm_,   bbtispc_), DiagMGMinCommOseenPc_     ( MGPcsymm_,    mincommispc_),
        DiagMGISPreOseenPc_      ( MGPcsymm_,   isprepc_), DiagGMResMinCommPc_       ( GMResPc_,    mincommispc_),
        DiagPCGBBTOseenPc_       ( PCGPc_,      bbtispc_), DiagPCGMinCommOseenPc_    ( PCGPc_,   mincommispc_),
        DiagPCGISPreOseenPc_     ( PCGPc_,      isprepc_),
        LBlockMGBBTOseenPc_      ( MGPc_,       bbtispc_), LBlockMGMinCommOseenPc_   ( MGPc_,    mincommispc_),    LBlockMGISPreOseenPc_( MGPc_, isprepc_),
        LBlockMGVankaOseenPc_    ( MGPc_,       vankaschurpc_),
        LBlockGMResBBTOseenPc_   ( GMResPc_,    bbtispc_), LBlockGMResMinCommOseenPc_( GMResPc_, mincommispc_),
        LBlockBiCGBBTOseenPc_    ( BiCGStabPc_, bbtispc_), LBlockBiCGMinCommOseenPc_ ( BiCGStabPc_, mincommispc_),
        vankapc_                 ( &Stokes.pr_idx),
        // GCR solver
        GCRMGBBT_           ( LBlockMGBBTOseenPc_,        C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRMGMinComm_       ( LBlockMGMinCommOseenPc_,    C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRMGISPre_         ( LBlockMGISPreOseenPc_,      C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRMGVanka_         ( LBlockMGVankaOseenPc_,      C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRGMResBBT_        ( LBlockGMResBBTOseenPc_,     C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRGMResMinComm_    ( LBlockGMResMinCommOseenPc_, C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRBiCGStabBBT_     ( LBlockBiCGBBTOseenPc_,      C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRBiCGStabMinComm_ ( LBlockBiCGMinCommOseenPc_,  C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        GCRVanka_           ( vankapc_,                   C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false),
        // GMRes solver
        GMResMGBBT_           ( LBlockMGBBTOseenPc_,        C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResMGMinComm_       ( LBlockMGMinCommOseenPc_,    C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResGMResBBT_        ( LBlockGMResBBTOseenPc_,     C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResGMResMinComm_    ( LBlockGMResMinCommOseenPc_, C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResBiCGStabBBT_     ( LBlockBiCGBBTOseenPc_,      C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResBiCGStabMinComm_ ( LBlockBiCGMinCommOseenPc_,  C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol,
                                /*rel*/ false, false, RightPreconditioning),
        GMResRMGBBT_          ( LBlockMGBBTOseenPc_,        C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_InnerIter,
                                    C_.stk_OuterTol, C_.stk_InnerTol, /*rel*/ false),
        GMResRMGMinComm_      ( LBlockMGMinCommOseenPc_,    C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_InnerIter,
                                    C_.stk_OuterTol, C_.stk_InnerTol, /*rel*/ false),
        GMResRGMResBBT_       ( LBlockGMResBBTOseenPc_,     C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_InnerIter,
                                    C_.stk_OuterTol, C_.stk_InnerTol, /*rel*/ false),
        GMResRGMResMinComm_   ( LBlockGMResMinCommOseenPc_, C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_InnerIter,
                                    C_.stk_OuterTol, C_.stk_InnerTol, /*rel*/ false),
        GMResRBiCGStabBBT_    ( LBlockBiCGBBTOseenPc_,      C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_InnerIter,
                                    C_.stk_OuterTol, C_.stk_InnerTol, /*rel*/ false),
        GMResRBiCGStabMinComm_( LBlockBiCGMinCommOseenPc_,  C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_InnerIter,
                                    C_.stk_OuterTol, C_.stk_InnerTol, /*rel*/ false),
        // lanczos objects
        lanczos1_ (DiagMGBBTOseenPc_),   lanczos2_ (DiagMGMinCommOseenPc_),
        lanczos3_ (DiagPCGBBTOseenPc_),  lanczos4_ (DiagPCGMinCommOseenPc_),
        lanczos5_ (DiagMGISPreOseenPc_), lanczos6_ (DiagPCGISPreOseenPc_),
        // PMinRes solver
        PMinResMGBBT_     ( lanczos1_, C_.stk_OuterIter, C_.stk_OuterTol, /*relative*/ false),
        PMinResMGMinComm_ ( lanczos2_, C_.stk_OuterIter, C_.stk_OuterTol, /*relative*/ false),
        PMinResPCGBBT_    ( lanczos3_, C_.stk_OuterIter, C_.stk_OuterTol, /*relative*/ false),
        PMinResPCGMinComm_( lanczos4_, C_.stk_OuterIter, C_.stk_OuterTol, /*relative*/ false),
        PMinResMGISPre_   ( lanczos5_, C_.stk_OuterIter, C_.stk_OuterTol, /*relative*/ false),
        PMinResPCGISPre_  ( lanczos6_, C_.stk_OuterIter, C_.stk_OuterTol, /*relative*/ false),
        // coarse grid/direct solver for StokesMGM
        minressolver( lanczos3_, 500, 1e-6, true), blockminressolver(minressolver),
        gcrsolver( DiagGMResMinCommPc_, 500, 500, 1e-6, true), blockgcrsolver(gcrsolver),
        vankasmoother( 0, 0.8, &Stokes.pr_idx)
        {}

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
StokesSolverBaseCL* StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::CreateStokesSolver()
{
    if (Stokes_.UsesXFEM())
    { // check whether solver is well-defined for XFEM
        if (C_.stk_StokesMethod/10000 == 30)
            throw DROPSErrCL("StokesMGM not implemented for P1X-elements");
        if (C_.stk_StokesMethod%100 == 7)
            throw DROPSErrCL("ISMGPreCL not implemented for P1X-elements");
    }
    
    StokesSolverBaseCL* stokessolver = 0;
    switch (C_.stk_StokesMethod)
    {
        case 10101 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockMGBBTOseenPcT> >       ( GCRMGBBT_);
        break;
        case 10102 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockMGMinCommOseenPcT> >   ( GCRMGMinComm_);
        break;
        case 10103 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockMGISPreOseenPcT> >     ( GCRMGISPre_); 
        break;
        case 10401 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockGMResBBTOseenPcT> >    ( GCRGMResBBT_);
        break;
        case 10402 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockGMResMinCommOseenPcT> >( GCRGMResMinComm_);
        break;
        case 10501 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockBiCGBBTOseenPcT> >     ( GCRBiCGStabBBT_);
        break;
        case 10502 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockBiCGMinCommOseenPcT> > ( GCRBiCGStabMinComm_);
        break;
        case 10606 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<VankaPreCL> > ( GCRVanka_);
        break;
        case 10104 :
            stokessolver = new BlockMatrixSolverCL<GCRSolverCL<LBlockMGVankaOseenPcT> > ( GCRMGVanka_);
        break;
        case 20101 :
            stokessolver = new InexactUzawaCL<MGPcT, ISBBTPreCL, APC_OTHER>
                        ( MGPc_, bbtispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20102 :
            stokessolver = new InexactUzawaCL<MGPcT, MinCommPreCL, APC_OTHER>
                        ( MGPc_, mincommispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20103 :
            stokessolver = new InexactUzawaCL<MGPcT, ISPreCL, APC_OTHER>
                        ( MGPc_, isprepc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20107 :
            stokessolver = new InexactUzawaCL<MGPcT, ISMGPreCL, APC_OTHER>
                        ( MGPc_, ismgpre_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20201 :
            stokessolver = new InexactUzawaCL<MGsymmPcT, ISBBTPreCL, APC_SYM>
                        ( MGPcsymm_, bbtispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20202 :
            stokessolver = new InexactUzawaCL<MGsymmPcT, MinCommPreCL, APC_SYM>
                        ( MGPcsymm_, mincommispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20203 :
            stokessolver = new InexactUzawaCL<MGsymmPcT, ISPreCL, APC_SYM>
                        ( MGPcsymm_, isprepc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20207 :
            stokessolver = new InexactUzawaCL<MGsymmPcT, ISMGPreCL, APC_SYM>
                        ( MGPcsymm_, ismgpre_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20301 :
            stokessolver = new InexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM>
                        ( PCGPc_, bbtispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20302 :
            stokessolver = new InexactUzawaCL<PCGPcT, MinCommPreCL, APC_SYM>
                        ( PCGPc_, mincommispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20303 :
            stokessolver = new InexactUzawaCL<PCGPcT, ISPreCL, APC_SYM>
                        ( PCGPc_, isprepc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20401 :
            stokessolver = new InexactUzawaCL<GMResPcT, ISBBTPreCL, APC_OTHER>
                        ( GMResPc_, bbtispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20402 :
            stokessolver = new InexactUzawaCL<GMResPcT, MinCommPreCL, APC_OTHER>
                        ( GMResPc_, mincommispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20403 :
            stokessolver = new InexactUzawaCL<GMResPcT, ISPreCL, APC_OTHER>
                        ( GMResPc_, isprepc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20501 :
            stokessolver = new InexactUzawaCL<BiCGStabPcT, ISBBTPreCL, APC_OTHER>
                        ( BiCGStabPc_, bbtispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20502 :
            stokessolver = new InexactUzawaCL<BiCGStabPcT, MinCommPreCL, APC_OTHER>
                        ( BiCGStabPc_, mincommispc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20503 :
            stokessolver = new InexactUzawaCL<BiCGStabPcT, ISPreCL, APC_OTHER>
                        ( BiCGStabPc_, isprepc_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20507 :
            stokessolver = new InexactUzawaCL<BiCGStabPcT, ISMGPreCL, APC_OTHER>
                        ( BiCGStabPc_, ismgpre_, C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 30201 :
            stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos1T> >( PMinResMGBBT_);
        break;
        case 30202 :
            stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos2T> >( PMinResMGMinComm_);
        break;
        case 30203 :
            stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos5T> >( PMinResMGISPre_);
        break;
        case 30301 :
            stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos3T> >( PMinResPCGBBT_);
        break;
        case 30302 :
            stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos4T> >( PMinResPCGMinComm_);
        break;
        case 30303 :
            stokessolver = new BlockMatrixSolverCL<PMResSolverCL<Lanczos6T> >( PMinResPCGISPre_);
        break;
        case 40101 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockMGBBTOseenPcT> >        ( GMResMGBBT_);
        break;
        case 40102 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockMGMinCommOseenPcT> >    ( GMResMGMinComm_);
        break;
        case 40401 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockGMResBBTOseenPcT> >     ( GMResGMResBBT_);
        break;
        case 40402 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockGMResMinCommOseenPcT> > ( GMResGMResMinComm_);
        break;
        case 40501 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockBiCGBBTOseenPcT> >      ( GMResBiCGStabBBT_);
        break;
        case 40502 :
            stokessolver = new BlockMatrixSolverCL<GMResSolverCL<LBlockBiCGMinCommOseenPcT> >  ( GMResBiCGStabMinComm_);
        break;
        case 50101 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockMGBBTOseenPcT> >       ( GMResRMGBBT_);
        break;
        case 50102 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockMGMinCommOseenPcT> >   ( GMResRMGMinComm_);
        break;
        case 50401 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockGMResBBTOseenPcT> >    ( GMResRGMResBBT_);
        break;
        case 50402 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockGMResMinCommOseenPcT> >( GMResRGMResMinComm_);
        break;
        case 50501 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockBiCGBBTOseenPcT> >     ( GMResRBiCGStabBBT_);
        break;
        case 50502 :
            stokessolver = new BlockMatrixSolverCL<GMResRSolverCL<LBlockBiCGMinCommOseenPcT> > ( GMResRBiCGStabMinComm_);
        break;
        case 303030 : {
            if (C_.ns_Nonlinear==0.0){ // stokes
                mgvankasolversymm_ = new StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>
                           ( Stokes_.prM.Data, vankasmoother, blockminressolver, C_.stk_OuterIter, C_.stk_OuterTol, false, 2),
                stokessolver = mgvankasolversymm_;
            }
            else {
                mgvankasolver_ = new StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>
                           ( Stokes_.prM.Data, vankasmoother, blockgcrsolver, C_.stk_OuterIter, C_.stk_OuterTol, false, 2),
                stokessolver = mgvankasolver_;
            }
        }
        break;
        case 303131 : {
            if (C_.ns_Nonlinear==0.0){ // stokes
                mgbssolversymm_ = new StokesMGSolverCL<BSSmootherCL, ProlongationVelT, ProlongationPT>
                           ( Stokes_.prM.Data, bssmoother, blockminressolver, C_.stk_OuterIter, C_.stk_OuterTol, false, 2),
                stokessolver = mgbssolversymm_;
            }
            else {
                mgbssolversymm_ = new StokesMGSolverCL<BSSmootherCL, ProlongationVelT, ProlongationPT>
                           ( Stokes_.prM.Data, bssmoother, blockgcrsolver, C_.stk_OuterIter, C_.stk_OuterTol, false, 2),
                stokessolver = mgbssolver_;
            }
        }
        break;
        default: throw DROPSErrCL("StokesSolverFactoryCL: Unknown Stokes solver");
    }
    return stokessolver;
}

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::
    SetMatrices( const MatrixCL* A, const MatrixCL* B, const MatrixCL* Mvel, const MatrixCL* M, const IdxDescCL* pr_idx) {
    if ( APc_ == 30 || APc_ == 31) { //  Vanka or Braess Sarazin smoother
        mincommispc_.SetMatrices(A, B, Mvel, M, pr_idx);
        bbtispc_.SetMatrices(B, Mvel, M, pr_idx);
    }
    if ( SPc_ == 4) {              // VankaSchur
        vankaschurpc_.SetAB(A, B);
    }
}

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
ProlongationVelT* StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::GetPVel()
{
    switch ( APc_) {
        case  1 : return MGSolver_.GetProlongation();     break;  // general MG
        case  2 : return MGSolversymm_.GetProlongation(); break;  // symm. MG
        case 30 : return (C_.ns_Nonlinear == 0 ? mgvankasolversymm_->GetPVel() :  mgvankasolver_->GetPVel()); break;
        case 31 : return (C_.ns_Nonlinear == 0 ? mgbssolversymm_->GetPVel()    :  mgbssolver_->GetPVel());    break;
    }
    return 0;
}

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
ProlongationPT* StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::GetPPr()
{
    switch ( APc_) {
        case 30 : return (C_.ns_Nonlinear == 0 ? mgvankasolversymm_->GetPPr() :  mgvankasolver_->GetPPr()); break;
        case 31 : return (C_.ns_Nonlinear == 0 ? mgbssolversymm_->GetPPr()    :  mgbssolver_->GetPPr());    break;
    }
    if (SPc_ == 7 ) // ISMGPreCL
        return ismgpre_.GetProlongation();
    return 0;
}

#else // parallel part

/// \brief Structure that contains all neccessary parameter for the
///     ParStokesSolverFactoryCL
/** See documentation of parameter classes for detailed information*/
struct StokesSolverParamST
{
    int stk_StokesMethod;
    int tm_NumSteps;
    double tm_StepSize;
    int stk_OuterIter;
    double stk_OuterTol;
    int stk_InnerIter;
    double stk_InnerTol;
    int stk_PcAIter;
    double stk_PcATol;
    double stk_PcSTol;
    double stk_Theta;

    /// \brief Constructor which copies all values out of a parameter class
    ///   into this parameter class
    template <typename ParamT>
    StokesSolverParamST(const ParamT& C)
      : stk_StokesMethod(C.stk_StokesMethod), tm_NumSteps(C.tm_NumSteps), tm_StepSize(C.tm_StepSize),
        stk_OuterIter(C.stk_OuterIter), stk_OuterTol(C.stk_OuterTol),
        stk_InnerIter(C.stk_InnerIter), stk_InnerTol(C.stk_InnerTol),
        stk_PcAIter(C.stk_PcAIter), stk_PcATol(C.stk_PcATol), stk_PcSTol(C.stk_PcSTol),
        stk_Theta(C.stk_Theta)
    {}
};


/*************************************************************
*   S t o k e s S o l v e r F a c t o r y  C L               *
**************************************************************/
/// \brief Factory for producing parallel stokes solver
template <class StokesT, class ParamsT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryCL : public StokesSolverFactoryBaseCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>
{
  private:
    typedef StokesSolverFactoryBaseCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT> base_;
    using base_::Stokes_;
    using base_::C_;
    using base_::OseenSolver_;
    using base_::APc_;
    using base_::SPc_;
    double kA_, kM_;

// generic preconditioners
    ParDummyPcCL DummyPrPc_;
    ParDummyPcCL DummyVelPc_;
    ParJac0CL    JACPrPc_;
    ParJac0CL    JACVelPc_;

// PC for instat. Schur complement
    ISBBTPreCL      bbtispc_;

// PC for A-block
    //JAC-GMRes
    typedef ParPreGMResSolverCL<ParJac0CL> GMResSolverT;
    typedef SolverAsPreCL<GMResSolverT>    GMResPcT;
    GMResSolverT GMResSolver_;
    GMResPcT GMResPc_;

    //JAC-PCG
    typedef ParPCGSolverCL<ParJac0CL> PCGSolverT;
    typedef SolverAsPreCL<PCGSolverT> PCGPcT;
    PCGSolverT PCGSolver_;
    PCGPcT PCGPc_;

// BlockPC
    typedef BlockPreCL<GMResPcT, ISBBTPreCL, LowerBlockPreCL> LBlockGMResBBTOseenPcT;
    LBlockGMResBBTOseenPcT LBlockGMResBBTOseenPc_;

//GCR solver
    ParPreGCRSolverCL<LBlockGMResBBTOseenPcT> GCRGMResBBT_;

  public:
    StokesSolverFactoryCL( StokesT& Stokes, ParamsT& C);
    ~StokesSolverFactoryCL() {}

    /// Nothing is to be done in parallel, because special preconditioners does not exist
    void       SetMatrixA ( const MatrixCL*)  {};
    /// Nothing is to be done in parallel, because special preconditioners does not exist
    void       SetMatrices( const MatrixCL*, const MatrixCL*, const MatrixCL*, const MatrixCL*, const IdxDescCL*){}
    /// Nothing is to be done in parallel, because special preconditioners does not exist
    ProlongationVelT* GetPVel() { return 0; }
    /// Nothing is to be done in parallel, because special preconditioners does not exist
    ProlongationPT*   GetPPr()  { return 0; }
    /// Returns a stokes solver with specifications from ParamsT C
    StokesSolverBaseCL* CreateStokesSolver();
};

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
  StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::StokesSolverFactoryCL(StokesT& Stokes, ParamsT& C)
    : base_(Stokes, C),
      kA_(C_.tm_NumSteps != 0 ? 1.0/C_.tm_StepSize : 0.0), // C_.tm_NumSteps == 0: stat. problem
      kM_(C_.stk_Theta),
      DummyPrPc_( Stokes.pr_idx.GetFinest()), DummyVelPc_( Stokes.vel_idx.GetFinest()),
      JACPrPc_( Stokes.pr_idx.GetFinest()), JACVelPc_( Stokes.vel_idx.GetFinest()),
      bbtispc_ ( Stokes_.B.Data.GetFinestPtr(), Stokes_.prM.Data.GetFinestPtr(), Stokes_.M.Data.GetFinestPtr(),
                 Stokes.pr_idx.GetFinest(), Stokes.vel_idx.GetFinest(), kA_, kM_, C_.stk_PcSTol, C_.stk_PcSTol),
      GMResSolver_(/*restart*/ 100, C_.stk_PcAIter, C_.stk_PcATol, Stokes.vel_idx.GetFinest(), JACVelPc_,
                   /*rel*/ true, /*accure*/ true, /*ModGS*/ false),
      GMResPc_( GMResSolver_),
      PCGSolver_(C_.stk_PcAIter, C_.stk_PcATol, Stokes.vel_idx.GetFinest(), JACVelPc_,
                 /*rel*/ true, /*acc*/ true),
      PCGPc_(PCGSolver_),
      LBlockGMResBBTOseenPc_( GMResPc_, bbtispc_),
      GCRGMResBBT_( C.stk_OuterIter, C.stk_OuterIter, C.stk_OuterTol, LBlockGMResBBTOseenPc_, true, false, true, &std::cout)
    {}

template <class StokesT, class ParamsT, class ProlongationVelT, class ProlongationPT>
  StokesSolverBaseCL* StokesSolverFactoryCL<StokesT, ParamsT, ProlongationVelT, ProlongationPT>::CreateStokesSolver()
{
    StokesSolverBaseCL* stokessolver = 0;
    switch (C_.stk_StokesMethod)
    {
        case 20301 :
            stokessolver = new ParInexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM>
                        ( PCGPc_, bbtispc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                          C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 20400 :
            stokessolver = new ParInexactUzawaCL<GMResPcT, ParDummyPcCL, APC_OTHER>
                         ( GMResPc_, DummyPrPc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                           C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol, 500, &std::cout);
        break;
        case 20401 :
            stokessolver = new ParInexactUzawaCL<GMResPcT, ISBBTPreCL, APC_OTHER>
                        ( GMResPc_, bbtispc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                          C_.stk_OuterIter, C_.stk_OuterTol, C_.stk_InnerTol);
        break;
        case 10401 :
            stokessolver = new BlockMatrixSolverCL<ParPreGCRSolverCL<LBlockGMResBBTOseenPcT> >
                        ( GCRGMResBBT_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest());
        break;
        default: throw DROPSErrCL("Unknown StokesMethod");
    }
    return stokessolver;
}
#endif

} // end of namespace DROPS

