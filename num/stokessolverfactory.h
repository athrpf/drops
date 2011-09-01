/// \file stokessolverfactory.h
/// \brief creates several standard Stokes-solver
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen:

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef _PAR
#include "num/stokessolver.h"
#else
#include "num/parstokessolver.h"
#ifdef _HYPRE
#include "num/hypre.h"
#endif
#endif

#include "misc/params.h"

namespace DROPS {

/// codes for Oseen solvers
enum OseenSolverE {
    GCR_OS= 1, iUzawa_OS= 2, MinRes_OS= 3, GMRes_OS= 4, GMResR_OS= 5, IDRs_OS= 7, StokesMGM_OS= 30
};

/// codes for velocity preconditioners (also including smoothers for the StokesMGM_OS)
enum APcE {
    MG_APC= 1, MGsymm_APC= 2, PCG_APC= 3, GMRes_APC= 4, BiCGStab_APC= 5, VankaBlock_APC= 6, AMG_APC= 20, // preconditioners 
    PVanka_SM= 30, BraessSarazin_SM= 31, IDRs_APC=7                                                      // smoothers, nevertheless listed here
};

/// codes for the pressure Schur complement preconditioners
enum SPcE {
    ISBBT_SPC= 1, MinComm_SPC= 2, ISPre_SPC= 3, ISMG_SPC= 7, BDinvBT_SPC= 5, SIMPLER_SPC=8, MSIMPLER_SPC=9, VankaSchur_SPC= 4, VankaBlock_SPC=6, ISNonlinear_SPC=10
}; 

/// collects some information on the different Oseen solvers and preconditioners
struct StokesSolverInfoCL
{
    static std::string GetOseenSolverName( int solver) {
        switch(solver) {
            case GCR_OS:       return "GCR";
            case iUzawa_OS:    return "inexact Uzawa";
            case MinRes_OS:    return "PMinRes";
            case GMRes_OS:     return "GMRes";
            case GMResR_OS:    return "GMResR";
            case StokesMGM_OS: return "Stokes MG";
            case IDRs_OS:      return "IDR(s)";
            default:           return "unknown";
        }
    }
    static std::string GetVelPreName( int pre) {
        switch(pre) {
            case MG_APC:           return "multigrid V-cycle";
            case MGsymm_APC:       return "symm. multigrid V-cycle";
            case PCG_APC:          return "PCG iterations";
            case GMRes_APC:        return "GMRes iterations";
            case BiCGStab_APC:     return "BiCGStab iterations";
            case AMG_APC:          return "algebraic multigrid";
            case VankaBlock_APC:   return "block Vanka";
            case PVanka_SM:        return "Vanka smoother";
            case BraessSarazin_SM: return "Braess-Sarazin smoother";
            case IDRs_APC:         return "IDR(s) iterations";
            default:               return "unknown";
        }
    }
    static std::string GetSchurPreName( int pre) {
        switch(pre) {
            case ISBBT_SPC:        return "ISBBT (modified Cahouet-Chabard)";
            case MinComm_SPC:      return "MinComm (minimal commutator)";
            case ISPre_SPC:        return "ISPre (Cahouet-Chabard)";
            case ISNonlinear_SPC:  return "ISNonlinearPreCL (Cahouet-Chabard)";
            case ISMG_SPC:         return "ISMGPre (multigrid Cahouet-Chabard)";
            case BDinvBT_SPC:      return "B D^-1 B^T";
            case SIMPLER_SPC:      return "SIMPLER";
            case MSIMPLER_SPC:     return "MSIMPLER";
            case VankaSchur_SPC:   return "Vanka Schur";
            case VankaBlock_SPC:   return "block Vanka";
            case PVanka_SM:        return "Vanka smoother";
            case BraessSarazin_SM: return "Braess-Sarazin smoother";
            default:               return "unknown";
        }
    }
    static bool IsBlockPre( int pre) { return pre==VankaBlock_APC || pre==SIMPLER_SPC || pre==MSIMPLER_SPC; }
    static bool IsSmoother( int pre) { return pre==PVanka_SM || pre==BraessSarazin_SM; }
};

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
    <tr><td>  3 </td><td> MinRes            </td><td> PCG                                </td><td> ISPreCL                      </td></tr>
    <tr><td>  4 </td><td> GMRes             </td><td> GMRes                              </td><td> VankaSchurPreCL              </td></tr>
    <tr><td>  5 </td><td> GMResR            </td><td> BiCGStab                           </td><td> BD^{-1}BT                    </td></tr>
    <tr><td>  6 </td><td>                   </td><td> VankaPre                           </td><td> VankaPre                     </td></tr>
    <tr><td>  7 </td><td> IDR(s)            </td><td> IDR(s)                             </td><td> ISMGPreCL                    </td></tr>
    <tr><td>  8 </td><td>                   </td><td>                                    </td><td> SIMPLER                      </td></tr>
    <tr><td>  9 </td><td>                   </td><td>                                    </td><td> MSIMPLER                     </td></tr>
    <tr><td> 20 </td><td>                   </td><td> HYPRE-AMG                          </td><td>                              </td></tr>
    <tr><td> 30 </td><td> StokesMGM         </td><td> PVankaSmootherCL                   </td><td> PVankaSmootherCL             </td></tr>
    <tr><td> 31 </td><td>                   </td><td> BSSmootherCL                       </td><td> BSSmootherCL                 </td></tr>
    </table>*/
template <class StokesT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryBaseCL
{
  protected:
    StokesT& Stokes_;           ///< Stokes problem
    ParamCL& P_;                ///< Parameter for tolerances, iteration number, type of solver, ...
    int      SPc_,              ///< type of preconditioner for S
             APc_,              ///< type of preconditioner for A-block
             OseenSolver_;      ///< type of Oseen solver

  public:
    StokesSolverFactoryBaseCL( StokesT& Stokes, ParamCL& P) : Stokes_( Stokes), P_( P),
                               SPc_( P_.get<int>("Stokes.StokesMethod") % 100), APc_( (P_.get<int>("Stokes.StokesMethod") / 100) % 100),
                               OseenSolver_( (P_.get<int>("Stokes.StokesMethod") /10000) % 100) {}
    virtual ~StokesSolverFactoryBaseCL() {}

    /// print some infos about solver combination
    void PrintSolverInfo( std::ostream&) const;
    /// Set the A-block in the minimal commutator
    virtual void       SetMatrixA ( const MatrixCL*) = 0;
    /// Set all matrices in Schur complement preconditioner
    virtual void       SetMatrices( const MLMatrixCL*, const MLMatrixCL*, const MLMatrixCL*, const MLMatrixCL*, const MLIdxDescCL* pr_idx) = 0;
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
class StokesSolverFactoryHelperCL
{
  public:
    int GetOseenSolver( const ParamCL& P) const { return (P.get<int>("Stokes.StokesMethod") / 10000) % 100; }
    int GetAPc( const ParamCL& P) const { return (P.get<int>("Stokes.StokesMethod") / 100) % 100; }
    int GetSPc( const ParamCL& P) const { return P.get<int>("Stokes.StokesMethod") % 100; }
    bool VelMGUsed ( const ParamCL& P) const
    {
        const int APc = GetAPc( P);
        return (( APc == MG_APC) || (APc == MGsymm_APC) || (APc == PVanka_SM) || (APc == BraessSarazin_SM));
    }
    bool PrMGUsed  ( const ParamCL& P) const
    {
        const int APc = GetAPc( P),
            SPc = GetSPc( P);
        return (( APc == PVanka_SM) || ( APc == BraessSarazin_SM) || (SPc == ISMG_SPC));
    }
};

#ifndef _PAR
/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y  C L                     *
********************************************************************/
template <class StokesT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryCL : public StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT>
{
  private:
    typedef StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT> base_;
    using base_::Stokes_;
    using base_::P_;
    using base_::OseenSolver_;
    using base_::APc_;
    using base_::SPc_;
    using base_::PrintSolverInfo;
    double kA_, kM_;

// generic preconditioners
    JACPcCL  JACPc_;
    SSORPcCL SSORPc_;

// PC for instat. Schur complement
    SchurPreBaseCL  *spc_;
    ISBBTPreCL      bbtispc_;
    MinCommPreCL    mincommispc_;
    BDinvBTPreCL    bdinvbtispc_;
    VankaSchurPreCL vankaschurpc_;
    ISPreCL         isprepc_;
    ISMGPreCL       ismgpre_;
    
    PCG_SsorCL isnonlinearprepc_;
    ISNonlinearPreCL<PCG_SsorCL> isnonlinearpc_;

// PC for A-block
    PreBaseCL *apc_;
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

    //IDR(s)
    typedef IDRsSolverCL<SSORPcCL> IDRsSolverT;
    IDRsSolverT IDRsSolver_;
    typedef SolverAsPreCL<IDRsSolverT> IDRsPcT;
    IDRsPcT IDRsPc_;

// Block PC for Oseen problem
    typedef BlockPreCL<PreBaseCL, SchurPreBaseCL, DiagSpdBlockPreCL>  DiagBlockPcT;
    typedef BlockPreCL<PreBaseCL, SchurPreBaseCL, LowerBlockPreCL>    LowerBlockPcT;
    typedef BlockPreCL<PreBaseCL, BDinvBTPreCL, SIMPLERBlockPreCL>    SIMPLERBlockPcT;
    typedef BlockPreCL<PreBaseCL, SchurPreBaseCL, UpperBlockPreCL>    UpperBlockPcT;

    DiagBlockPcT    *DBlock_;
    LowerBlockPcT   *LBlock_;
    UpperBlockPcT   *UBlock_;
    SIMPLERBlockPcT *SBlock_;
    VankaPreCL      vankapc_;

//GCR solver
    typedef GCRSolverCL<LowerBlockPcT>   GCR_LBlockT;
    typedef GCRSolverCL<SIMPLERBlockPcT> GCR_SBlockT;
    typedef GCRSolverCL<VankaPreCL>      GCR_VankaT;

    GCR_LBlockT *GCRLBlock_;
    GCR_SBlockT *GCRSBlock_;
    GCR_VankaT  *GCRVanka_;

//GMRes solver
    typedef GMResSolverCL<LowerBlockPcT> GMRes_LBlockT;
    typedef GMResSolverCL<VankaPreCL>    GMRes_VankaT;

    GMRes_LBlockT *GMResLBlock_;
    GMRes_VankaT  *GMResVanka_;

// GMResR solver
    typedef GMResRSolverCL<LowerBlockPcT> GMResR_LBlockT;
    typedef GMResRSolverCL<VankaPreCL>    GMResR_VankaT;

    GMResR_LBlockT *GMResRLBlock_;
    GMResR_VankaT  *GMResRVanka_;

// Lanczos
    typedef PLanczosONBCL<VectorCL, DiagBlockPcT> LanczosT;

    LanczosT *lanczos_;

// MinRes solver
    typedef PMResSolverCL<LanczosT> MinResT;
    MinResT *MinRes_;

// IDR(s) solver
    typedef IDRsSolverCL<UpperBlockPcT> IDRs_UBlockT;
    IDRs_UBlockT *IDRsUBlock_;

// coarse grid solver
    DiagBlockPcT DiagPCGBBTOseenPc_, DiagGMResMinCommPc_;
    LanczosT lanczosPCGBBT_;
    MinResT minressolver_;
    BlockMatrixSolverCL<MinResT> coarse_blockminressolver_;

    //GCR solver
    GCRSolverCL<DiagBlockPcT> gcrsolver_;
    BlockMatrixSolverCL<GCRSolverCL<DiagBlockPcT> > coarse_blockgcrsolver_;

    PVankaSmootherCL vankasmoother_;
    BSSmootherCL bssmoother_;

    //StokesMGSolver
    StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>* mgvankasolver_;
    StokesMGSolverCL<BSSmootherCL,     ProlongationVelT, ProlongationPT>* mgbssolver_;

    PreBaseCL*      CreateAPc();
    SchurPreBaseCL* CreateSPc();

  public:
    StokesSolverFactoryCL(StokesT& Stokes, ParamCL& P);
    ~StokesSolverFactoryCL();

    // checks, whether the combination of Oseen solver, A and S preconditioner is valid.
    bool ValidSolverCombination( std::ostream* os= 0) const;
    /// Set the A-block in the minimal commutator
    void       SetMatrixA ( const MatrixCL* A) { mincommispc_.SetMatrixA(A); bdinvbtispc_.SetMatrixA(A); }
    /// Set all matrices in Schur complement preconditioner (only for StokesMGM)
    void       SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx);
    /// Returns pointer to prolongation for velocity
    ProlongationVelT* GetPVel();
    /// Returns pointer to prolongation for pressure
    ProlongationPT*   GetPPr();
    /// Returns a stokes solver with specifications from ParamsT C
    StokesSolverBaseCL* CreateStokesSolver();
    
    /// Returns a pointer to the used velocity preconditioner for the upper left block (aka A block)
    PreBaseCL*      GetVelPrePtr()   { return apc_; }
    /// Returns a pointer to the used pressure preconditioner for the schur complement
    SchurPreBaseCL* GetSchurPrePtr() { return spc_; }

    PVankaSmootherCL&      GetVankaSmoother () { return vankasmoother_; }
    VankaSchurPreCL&       GetVankaSchurPc ()  { return vankaschurpc_; }
};

template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    StokesSolverFactoryCL(StokesT& Stokes, ParamCL& P)
    : base_(Stokes, P),
        kA_(P.get<int>("Time.NumSteps") != 0 ? 1.0/P.get<double>("Time.StepSize") : 0.0), // P.get<int>("Time.NumSteps") == 0: stat. problem
        kM_(P.get<double>("Stokes.Theta")),
        // schur complement preconditioner
        bbtispc_    ( &Stokes_.B.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), Stokes_.pr_idx.GetFinest(), kA_, kM_, P.get<double>("Stokes.PcSTol"), P.get<double>("Stokes.PcSTol") /* enable regularization: , 0.707*/),
        mincommispc_( 0, &Stokes_.B.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(),Stokes_.pr_idx.GetFinest(), P.get<double>("Stokes.PcSTol") /* enable regularization: , 0.707*/),
        bdinvbtispc_( 0, &Stokes_.B.Data.GetFinest(), &Stokes_.M.Data.GetFinest(), &Stokes_.prM.Data.GetFinest(),Stokes_.pr_idx.GetFinest(), P.get<double>("Stokes.PcSTol") /* enable regularization: , 0.707*/),
        vankaschurpc_( &Stokes.pr_idx), isprepc_( Stokes.prA.Data, Stokes.prM.Data, kA_, kM_),
        ismgpre_( Stokes.prA.Data, Stokes.prM.Data, kA_, kM_),
        isnonlinearprepc_( SSORPc_, 100, P.get<double>("Stokes.PcSTol"), true),
        isnonlinearpc_( isnonlinearprepc_, Stokes_.prA.Data.GetFinest(), Stokes_.prM.Data.GetFinest(), kA_, kM_),
        // preconditioner for A
        smoother_( 1.0), coarsesolversymm_( SSORPc_, 500, 1e-6, true),
        MGSolversymm_ ( smoother_, coarsesolversymm_, P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), false),
        MGPcsymm_( MGSolversymm_),
        coarsesolver_( JACPc_, 500, 500, 1e-6, true),
        MGSolver_ ( smoother_, coarsesolver_, P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), false), MGPc_( MGSolver_),
        GMResSolver_( JACPc_, P.get<int>("Stokes.PcAIter"), /*restart*/ 100, P.get<double>("Stokes.PcATol"), /*rel*/ true), GMResPc_( GMResSolver_),
        BiCGStabSolver_( JACPc_, P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), /*rel*/ true),BiCGStabPc_( BiCGStabSolver_),
        PCGSolver_( SSORPc_, P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), true), PCGPc_( PCGSolver_),
        IDRsSolver_( SSORPc_, P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), true), IDRsPc_( IDRsSolver_),
        // block precondtioner
        DBlock_(0), LBlock_(0), SBlock_(0),
        vankapc_( &Stokes.pr_idx),
        // GCR solver
        GCRLBlock_(0), GCRSBlock_(0), GCRVanka_(0),
        // GMRes solver
        GMResLBlock_(0),  GMResVanka_(0),
        GMResRLBlock_(0), GMResRVanka_(0),
        // lanczos objects
        lanczos_ (0), 
        // PMinRes solver
        MinRes_(0),
        // IDRs solver
        IDRsUBlock_(0),
        // coarse grid/direct solver for StokesMGM
        DiagPCGBBTOseenPc_( PCGPc_, bbtispc_), DiagGMResMinCommPc_( GMResPc_, mincommispc_), lanczosPCGBBT_ (DiagPCGBBTOseenPc_),
        minressolver_( lanczosPCGBBT_, 500, 1e-6, true), coarse_blockminressolver_(minressolver_),
        gcrsolver_( DiagGMResMinCommPc_, 500, 500, 1e-6, true), coarse_blockgcrsolver_(gcrsolver_),
        vankasmoother_( 0, 0.8, &Stokes.pr_idx)
{
    apc_= CreateAPc();
    spc_= CreateSPc();
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    ~StokesSolverFactoryCL()
{
    delete MinRes_; delete lanczos_;
    delete GMResRVanka_; delete GMResRLBlock_;
    delete GMResVanka_; delete GMResLBlock_;
    delete GCRVanka_; delete GCRLBlock_; delete GCRSBlock_;
    delete SBlock_; delete LBlock_; delete DBlock_; delete IDRsUBlock_;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
bool StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::ValidSolverCombination( std::ostream* os) const
{
    std::string msg;
    bool ok= false;
    
    if (StokesSolverInfoCL::GetOseenSolverName(OseenSolver_)=="unknown")
        msg= "unknown Oseen solver";
    else if (StokesSolverInfoCL::GetVelPreName(APc_)=="unknown")
        msg= "unknown vel preconditioner";
    else if (StokesSolverInfoCL::GetSchurPreName(SPc_)=="unknown")
        msg= "unknown pr preconditioner";
    else if (OseenSolver_==StokesMGM_OS && (APc_!=SPc_ || !StokesSolverInfoCL::IsSmoother(APc_) )) 
        msg= "Stokes multigrid method requires smoother";
    else if ((StokesSolverInfoCL::IsSmoother(APc_) || StokesSolverInfoCL::IsSmoother(SPc_)) && OseenSolver_!=StokesMGM_OS)
        msg= "smoother makes no sense without multigrid solver";
    else if (OseenSolver_==iUzawa_OS && (StokesSolverInfoCL::IsBlockPre(APc_) || StokesSolverInfoCL::IsBlockPre(SPc_) ))
        msg= "block preconditioner not allowed for inexact Uzawa";
    else if (OseenSolver_==MinRes_OS && (StokesSolverInfoCL::IsBlockPre(APc_) || StokesSolverInfoCL::IsBlockPre(SPc_) ))
        msg= "MinRes requires diagonal block preconditioner";
    else if ((StokesSolverInfoCL::IsBlockPre(APc_) || StokesSolverInfoCL::IsBlockPre(SPc_)) && APc_!=SPc_ && SPc_!=SIMPLER_SPC && SPc_!=MSIMPLER_SPC)
        msg= "block preconditioner should be the same for vel and pr part";
    else // all tests passed successfully
        ok= true;

    if (os && !ok)
        (*os) << "invalid solver combination in StokesSolverFactoryCL:\t" << msg << std::endl;
    return ok;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT>::PrintSolverInfo( std::ostream& os) const
{
    os << "Oseen solver info:\t" << StokesSolverInfoCL::GetOseenSolverName( OseenSolver_)
       << "\n + vel precond.  :\t" << StokesSolverInfoCL::GetVelPreName( APc_)
       << "\n + pr  precond.  :\t" << StokesSolverInfoCL::GetSchurPreName( SPc_) << std::endl;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
PreBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateAPc()
{
    switch (APc_) {
        case MG_APC:       return &MGPc_;
        case MGsymm_APC:   return &MGPcsymm_;
        case PCG_APC:      return &PCGPc_;
        case GMRes_APC:    return &GMResPc_;
        case BiCGStab_APC: return &BiCGStabPc_;
        case IDRs_APC:     return &IDRsPc_;
        default:           return 0;
    }
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
SchurPreBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateSPc()
{
    switch (SPc_) {
        case ISBBT_SPC:      return &bbtispc_;
        case MinComm_SPC:    return &mincommispc_;
        case ISPre_SPC:      return &isprepc_;
        case ISMG_SPC:       return &ismgpre_;
        case SIMPLER_SPC:
        case MSIMPLER_SPC:
        case BDinvBT_SPC:    return &bdinvbtispc_;
        case VankaSchur_SPC: return &vankaschurpc_;
        case ISNonlinear_SPC:return &isnonlinearpc_;
        default:             return 0;
    }
}


template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateStokesSolver()
{
    PrintSolverInfo( std::cout);
    if (!ValidSolverCombination( &std::cout))
        throw DROPSErrCL("StokesSolverFactoryCL::CreateStokesSolver(): Invalid solver combination");
        
    if (Stokes_.UsesXFEM())
    { // check whether solver is well-defined for XFEM
        if (OseenSolver_ == StokesMGM_OS)
            throw DROPSErrCL("StokesMGM not implemented for P1X-elements");
        if (SPc_ == ISMG_SPC)
            throw DROPSErrCL("ISMGPreCL not implemented for P1X-elements");
    }

    StokesSolverBaseCL* stokessolver = 0;
    
    switch (OseenSolver_) {
        case iUzawa_OS: {
            if (APc_==MGsymm_APC) // symmetric A preconditionder -> use more efficient version of inexact Uzawa
                stokessolver= new InexactUzawaCL<PreBaseCL, SchurPreBaseCL, APC_SYM>  ( *apc_, *spc_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), P_.template get<int>("Stokes.InnerIter"));
            else
                stokessolver= new InexactUzawaCL<PreBaseCL, SchurPreBaseCL, APC_OTHER>( *apc_, *spc_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), P_.template get<int>("Stokes.InnerIter"));
        }
        break;
        
        case MinRes_OS: { // MinRes requires symmetric block preconditioner, hence we can only use diagonal block preconditioners
            DBlock_= new DiagBlockPcT( *apc_, *spc_);
            lanczos_= new LanczosT( *DBlock_);
            MinRes_= new MinResT( *lanczos_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), /*relative*/ false);
            stokessolver= new BlockMatrixSolverCL<MinResT>( *MinRes_);
        }      
        break;
        
        case GCR_OS: {
            if (APc_==VankaBlock_APC) {
                GCRVanka_= new GCR_VankaT( vankapc_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GCR_VankaT> ( *GCRVanka_);
            } else if (SPc_==SIMPLER_SPC || SPc_==MSIMPLER_SPC) {
                bdinvbtispc_.SetMassLumping( SPc_==MSIMPLER_SPC);
                SBlock_= new SIMPLERBlockPcT( *apc_, bdinvbtispc_);
                GCRSBlock_= new GCR_SBlockT( *SBlock_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GCR_SBlockT>( *GCRSBlock_);      
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                GCRLBlock_= new GCR_LBlockT( *LBlock_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GCR_LBlockT>( *GCRLBlock_);      
            }
        }
        break;
        
        case GMRes_OS: {
            if (APc_==VankaBlock_APC) {
                GMResVanka_= new GMRes_VankaT( vankapc_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), /*rel*/ false, false, RightPreconditioning);
                stokessolver= new BlockMatrixSolverCL<GMRes_VankaT> ( *GMResVanka_);
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                GMResLBlock_= new GMRes_LBlockT( *LBlock_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), /*rel*/ false, false, RightPreconditioning);
                stokessolver= new BlockMatrixSolverCL<GMRes_LBlockT>( *GMResLBlock_);      
            }
        }
        break;
        
        case GMResR_OS: {
            if (APc_==VankaBlock_APC) {
                GMResRVanka_= new GMResR_VankaT( vankapc_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.InnerIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GMResR_VankaT> ( *GMResRVanka_);
            } else {
                LBlock_= new LowerBlockPcT( *apc_, *spc_);
                GMResRLBlock_= new GMResR_LBlockT( *LBlock_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.InnerIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<GMResR_LBlockT>( *GMResRLBlock_);      
            }
        }
        break;
            
        case StokesMGM_OS: {
            if (APc_==PVanka_SM) {
                if (P_.template get<double>("NavierStokes.Nonlinear", 0.0)==0.0) // Stokes
                    mgvankasolver_ = new StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, vankasmoother_, coarse_blockminressolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), false, 2);
                else
                    mgvankasolver_ = new StokesMGSolverCL<PVankaSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, vankasmoother_, coarse_blockgcrsolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), false, 2);
                stokessolver = mgvankasolver_;
            }
            else if (APc_==BraessSarazin_SM) {
                if (P_.template get<double>("NavierStokes.Nonlinear", 0.0) ==0.0) // Stokes
                    mgbssolver_ = new StokesMGSolverCL<BSSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, bssmoother_, coarse_blockminressolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), false, 2);
                else
                    mgbssolver_ = new StokesMGSolverCL<BSSmootherCL, ProlongationVelT, ProlongationPT>
                               ( Stokes_.prM.Data, bssmoother_, coarse_blockgcrsolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), false, 2);
                stokessolver = mgbssolver_;
            }
        }
        break;
        
        case IDRs_OS: {
//            if (APc_==VankaBlock_APC) {
//                GCRVanka_= new GCR_VankaT( vankapc_,  C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false);
//                stokessolver= new BlockMatrixSolverCL<GCR_VankaT> ( *GCRVanka_);
//            } else if (SPc_==SIMPLER_SPC || SPc_==MSIMPLER_SPC) {
//                bdinvbtispc_.SetMassLumping( SPc_==MSIMPLER_SPC);
 //               SBlock_= new SIMPLERBlockPcT( *apc_, bdinvbtispc_);
  //              GCRSBlock_= new GCR_SBlockT( *SBlock_,  C_.stk_OuterIter, C_.stk_OuterIter, C_.stk_OuterTol, /*rel*/ false);
//                stokessolver= new BlockMatrixSolverCL<GCR_SBlockT>( *GCRSBlock_);
//            } else {
                UBlock_= new UpperBlockPcT( *apc_, *spc_);
                IDRsUBlock_= new IDRs_UBlockT( *UBlock_,  P_.template get<int>("Stokes.OuterIter"), P_.template get<int>("Stokes.OuterTol"), /*rel*/ false);
                stokessolver= new BlockMatrixSolverCL<IDRs_UBlockT>( *IDRsUBlock_);
//            }
        }
        break;

        default: throw DROPSErrCL("StokesSolverFactoryCL: Unknown Oseen solver");
    }
    if (stokessolver==0)
        throw DROPSErrCL("StokesSolverFactoryCL: Sorry, this solver combination is not implemented, yet");
    return stokessolver;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::
    SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx) {
    if ( APc_ == PVanka_SM || APc_ == BraessSarazin_SM) { //  Vanka or Braess Sarazin smoother
        bbtispc_.SetMatrices(B->GetCoarsestPtr(), Mvel->GetCoarsestPtr(), M->GetCoarsestPtr(), pr_idx->GetCoarsestPtr());
        return;
    }
    if ( SPc_ == VankaSchur_SPC) {              // VankaSchur
        vankaschurpc_.SetAB(A->GetCoarsestPtr(), B->GetCoarsestPtr());
    }
    else {
        mincommispc_.SetMatrices(A->GetFinestPtr(), B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
        bdinvbtispc_.SetMatrices(A->GetFinestPtr(), B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
        bbtispc_.SetMatrices(B->GetFinestPtr(), Mvel->GetFinestPtr(), M->GetFinestPtr(), pr_idx->GetFinestPtr());
    }
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ProlongationVelT* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::GetPVel()
{
    switch ( APc_) {
        case MG_APC           : return MGSolver_.GetProlongation();     break;  // general MG
        case MGsymm_APC       : return MGSolversymm_.GetProlongation(); break;  // symm. MG
        case PVanka_SM        : return mgvankasolver_->GetPVel(); break;
        case BraessSarazin_SM : return mgbssolver_->GetPVel();    break;
    }
    return 0;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ProlongationPT* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::GetPPr()
{
    switch ( APc_) {
        case PVanka_SM        : return mgvankasolver_->GetPPr(); break;
        case BraessSarazin_SM : return mgbssolver_->GetPPr();    break;
    }
    if (SPc_ == ISMG_SPC )
        return ismgpre_.GetProlongation();
    return 0;
}

#else // parallel part

/// \brief Structure that contains all neccessary parameter for the
///     ParStokesSolverFactoryCL
/** See documentation of parameter classes for detailed information*/
/*struct StokesSolverParamST
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
    StokesSolverParamST(const ParamCL& P)
      : stk_StokesMethod(P.get<int>("Stokes.StokesMethod")), tm_NumSteps(C.tm_NumSteps), tm_StepSize(C.tm_StepSize),
        stk_OuterIter(P.get<int>("Stokes.OuterIter")), stk_OuterTol(P.get<double>("Stokes.OuterTol")),
        stk_InnerIter(C.stk_InnerIter), stk_InnerTol(C.stk_InnerTol),
        stk_PcAIter(C.stk_PcAIter), stk_PcATol(C.stk_PcATol), stk_PcSTol(C.stk_PcSTol),
        stk_Theta(C.stk_Theta)
    {}
};*/


/*************************************************************
*   S t o k e s S o l v e r F a c t o r y  C L               *
**************************************************************/
/// \brief Factory for producing parallel stokes solver
template <class StokesT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryCL : public StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT>
{
  private:
    typedef StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT> base_;
    using base_::Stokes_;
    using base_::P_;
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

#ifdef _HYPRE
     //Algebraic MG solver
    HypreAMGSolverCL hypreAMG_;
    typedef SolverAsPreCL<HypreAMGSolverCL> AMGPcT;
    AMGPcT AMGPc_;
    typedef BlockPreCL<AMGPcT, ISBBTPreCL, LowerBlockPreCL> LBlockAMGBBTOseenPcT;
    LBlockAMGBBTOseenPcT LBlockAMGBBTOseenPc_;
    ParPreGCRSolverCL<LBlockAMGBBTOseenPcT> GCRAMGBBT_;
#endif

  public:
    StokesSolverFactoryCL( StokesT& Stokes, ParamCL& P);
    ~StokesSolverFactoryCL() {}

    /// Nothing is to be done in parallel, because special preconditioners does not exist
    void       SetMatrixA ( const MatrixCL*)  {};
    /// Nothing is to be done in parallel, because special preconditioners does not exist
    void       SetMatrices( const MLMatrixCL*, const MLMatrixCL*, const MLMatrixCL*, const MLMatrixCL*, const MLIdxDescCL*){}
    /// Nothing is to be done in parallel, because special preconditioners does not exist
    ProlongationVelT* GetPVel() { return 0; }
    /// Nothing is to be done in parallel, because special preconditioners does not exist
    ProlongationPT*   GetPPr()  { return 0; }
    /// Returns a stokes solver with specifications from ParamsT C
    StokesSolverBaseCL* CreateStokesSolver();
    /// Returns a pointer to the schur complement preconditioner
    SchurPreBaseCL* GetSchurPrePtr() { return &bbtispc_; }
};

template <class StokesT, class ProlongationVelT, class ProlongationPT>
  StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::StokesSolverFactoryCL(StokesT& Stokes, ParamCL& P)
    : base_(Stokes, P),
      kA_(P.get<int>("Time.NumSteps") != 0 ? 1.0/P.get<double>("Time.StepSize") : 0.0), // P.get<int>("Time.NumSteps") == 0: stat. problem
      kM_(P.get<double>("Stokes.Theta")),
      DummyPrPc_( Stokes.pr_idx.GetFinest()), DummyVelPc_( Stokes.vel_idx.GetFinest()),
      JACPrPc_( Stokes.pr_idx.GetFinest()), JACVelPc_( Stokes.vel_idx.GetFinest()),
      bbtispc_ ( Stokes_.B.Data.GetFinestPtr(), Stokes_.prM.Data.GetFinestPtr(), Stokes_.M.Data.GetFinestPtr(),
                 Stokes.pr_idx.GetFinest(), Stokes.vel_idx.GetFinest(), kA_, kM_, P.get<double>("Stokes.PcSTol"), P.get<double>("Stokes.PcSTol")),
      GMResSolver_(/*restart*/ 100, P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), Stokes.vel_idx.GetFinest(), JACVelPc_,
                   /*rel*/ true, /*accure*/ true, /*ModGS*/ false),
      GMResPc_( GMResSolver_),
      PCGSolver_(P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol"), Stokes.vel_idx.GetFinest(), JACVelPc_,
                 /*rel*/ true, /*acc*/ true),
      PCGPc_(PCGSolver_),
      LBlockGMResBBTOseenPc_( GMResPc_, bbtispc_),
      GCRGMResBBT_( P.get<int>("Stokes.OuterIter"), P.get<int>("Stokes.OuterIter"), P.get<double>("Stokes.OuterTol"), LBlockGMResBBTOseenPc_, true, false, true, &std::cout)
#ifdef _HYPRE
      , hypreAMG_( Stokes.vel_idx.GetFinest(), P.get<int>("Stokes.PcAIter"), P.get<double>("Stokes.PcATol")), AMGPc_(hypreAMG_),
      LBlockAMGBBTOseenPc_( AMGPc_, bbtispc_),
      GCRAMGBBT_( P.get<int>("Stokes.OuterIter"), P.get<int>("Stokes.OuterIter"), P.get<double>("Stokes.OuterTol"), LBlockAMGBBTOseenPc_, true, false, true, &std::cout)
#endif
    {}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
  StokesSolverBaseCL* StokesSolverFactoryCL<StokesT, ProlongationVelT, ProlongationPT>::CreateStokesSolver()
{
    StokesSolverBaseCL* stokessolver = 0;
    switch (P_.template get<int>("Stokes.StokesMethod"))
    {
        case 20301 :
            stokessolver = new ParInexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM>
                        ( PCGPc_, bbtispc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                          P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), P_.template get<int>("Stokes.InnerIter"), &std::cout);
        break;
        case 20400 :
            stokessolver = new ParInexactUzawaCL<GMResPcT, ParDummyPcCL, APC_OTHER>
                         ( GMResPc_, DummyPrPc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                           P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), P_.template get<int>("Stokes.InnerIter"), &std::cout);
        break;
        case 20401 :
            stokessolver = new ParInexactUzawaCL<GMResPcT, ISBBTPreCL, APC_OTHER>
                        ( GMResPc_, bbtispc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                          P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), P_.template get<int>("Stokes.InnerIter"), &std::cout);
        break;
        case 10401 :
            stokessolver = new BlockMatrixSolverCL<ParPreGCRSolverCL<LBlockGMResBBTOseenPcT> >
                        ( GCRGMResBBT_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest());
        break;
#ifdef _HYPRE
        case 22001 :
            stokessolver = new ParInexactUzawaCL<AMGPcT, ISBBTPreCL, APC_OTHER>
                        ( AMGPc_, bbtispc_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest(),
                          P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.InnerTol"), P_.template get<int>("Stokes.InnerIter"), &std::cout);
        break;template 
        case 12001 :
            stokessolver = new BlockMatrixSolverCL<ParPreGCRSolverCL<LBlockAMGBBTOseenPcT> >
                        ( GCRAMGBBT_, Stokes_.vel_idx.GetFinest(), Stokes_.pr_idx.GetFinest());
        break;
#endif
        default: throw DROPSErrCL("Unknown StokesMethod");
    }
    return stokessolver;
}
#endif


#ifndef _PAR
//coding for obsolete solvers!!!
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
    <tr><td> 50 00 00 </td><td> Uzawa+PCG               </td><td>                        </td><td>                              </td></tr>
    <tr><td> 50 01 00 </td><td> MG-Uzawa                </td><td>                        </td><td>                              </td></tr>
    <tr><td> 50 11 01 </td><td> Uzawa+ISMGPre+MGPCT     </td><td>                        </td><td>                              </td></tr>
    <tr><td> 51 00 00 </td><td> Schur+PCG               </td><td>                        </td><td>                              </td></tr>
    <tr><td> 51 01 00 </td><td> MG-Schur                </td><td>                        </td><td>                              </td></tr>
    <tr><td> 51 00 01 </td><td> Schur+PCG(GS)           </td><td>                        </td><td>                              </td></tr>
    <tr><td> 51 11 02 </td><td> PSchur-Full-MG          </td><td>                        </td><td>                              </td></tr>
    <tr><td> 51 00 03 </td><td> PSchur-PCGPre           </td><td>                        </td><td>                              </td></tr>
    <tr><td> 51 10 04 </td><td> PSchur-PCG-Pr-MG        </td><td>                        </td><td>                              </td></tr>
    <tr><td> 52 00 00 </td><td> Stokes-Minres           </td><td>                        </td><td>                              </td></tr>
    <tr><td> 53 00 00 </td><td> Stokes-PMinres          </td><td>                        </td><td>                              </td></tr>
    <tr><td> 53 11 01 </td><td> Stokes-PMinres_FullMG   </td><td>                        </td><td>                              </td></tr>
    <tr><td> 54 00 00 </td><td> Another Schur           </td><td>                        </td><td>                              </td></tr>
    <tr><td> 55 00 00 </td><td> Schur no PC             </td><td>                        </td><td>                              </td></tr>
    </table>*/

/*******************************************************************
*   S t o k e s S o l v e r F a c t o r y O b s o l e t e H e l p e r  C L         *
********************************************************************/

class StokesSolverFactoryObsoleteHelperCL
{
  public:
    bool VelMGUsed ( const ParamCL& P) const
    {
        const int APc = (P.get<int>("Stokes.StokesMethod") / 100) % 100;
        return (( APc == 1) || (APc == 11));
    }
    bool PrMGUsed  ( const ParamCL& P) const
    {
        const int SPc = (P.get<int>("Stokes.StokesMethod") / 1000) % 10;
        return (SPc == 1);
    }
};



/*******************************************************************
 *   S t o k e s S o l v e r F a c t o r y O b s o l e t e C L     *
 *******************************************************************/
template <class StokesT, class ProlongationVelT= MLMatrixCL, class ProlongationPT= MLMatrixCL>
class StokesSolverFactoryObsoleteCL : public StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT>
{
  private:
    typedef StokesSolverFactoryBaseCL<StokesT, ProlongationVelT, ProlongationPT> base_;
    using base_::Stokes_;
    using base_::P_;
    using base_::OseenSolver_;
    using base_::APc_;
    using base_::SPc_;
    double kA_, kM_;

    //solver
    SSORPcCL   ssor_;
    PCG_SsorCL PCGsolver_;
    CGSolverCL CGsolver_;
    SGSPcCL    sgs_;
    PCG_SgsCL  PCGsgssolver_;

    //Minres
    LanczosONBCL<VectorCL> q_;
    typedef PMResSolverCL<LanczosONBCL<VectorCL> > MinresSPT;
    MinresSPT minressolver_;

    //PMinres
    typedef SolverAsPreCL<PCG_SsorCL> PPcT;
    typedef BlockPreCL<PPcT, ISPreCL,DiagSpdBlockPreCL> BlockDiagPCGPreCL;
    typedef PMResSolverCL<PLanczosONBCL<VectorCL, BlockDiagPCGPreCL> > PMinresSP_DiagPCGT;

    PCG_SsorCL PPA_;
    PPcT PA_;
    ISPreCL PS_;
    BlockDiagPCGPreCL pre_;
    PLanczosONBCL<VectorCL, BlockDiagPCGPreCL> pq_;
    PMinresSP_DiagPCGT pminressolver_;

    //multigrid solver
    SSORPcCL ssorom_;
    SSORsmoothCL smoother_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> MGsolver_;

    //MG preconditioner
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > MGPCT;
    MGPCT MGpc_;
    ISMGPreCL ismgpcp_;

    //Pshur2-MG
    PCGSolverCL<ISMGPreCL> PCGMGPresolver_;

    //Pschur-PCG-Pre
    ISPreCL ispcp_;
    PCGSolverCL<ISPreCL> PCGPresolver_;

    //PMinresSP_FullMG
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > APcT;
    typedef BlockPreCL<APcT, ISMGPreCL,DiagSpdBlockPreCL> PcT;
    typedef PMResSolverCL<PLanczosONBCL<DROPS::VectorCL,PcT> > PMinresSP_FullMG;

    PMinresSP_FullMG pminresMGsolver_;
    PcT preMG_;
    APcT apc_;
    PLanczosONBCL<VectorCL, PcT> pqMG_;

  public:
    StokesSolverFactoryObsoleteCL(StokesT& Stokes, ParamCL& P);
    ~StokesSolverFactoryObsoleteCL() {}

    /// Set the A-block in the minimal commutator
    void       SetMatrixA ( const MatrixCL* /*A*/) { }
    /// Set all matrices in Schur complement preconditioner (only for StokesMGM)
    void       SetMatrices( const MLMatrixCL* A, const MLMatrixCL* B, const MLMatrixCL* Mvel, const MLMatrixCL* M, const MLIdxDescCL* pr_idx);
    /// Returns pointer to prolongation for velocity
    ProlongationVelT* GetPVel();
    /// Returns pointer to prolongation for pressure
    ProlongationPT*   GetPPr();
    /// Returns a stokes solver with specifications from ParamsT C
    StokesSolverBaseCL* CreateStokesSolver();

};

template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverFactoryObsoleteCL<StokesT, ProlongationVelT, ProlongationPT>::
    StokesSolverFactoryObsoleteCL(StokesT& Stokes, ParamCL& P)
    : base_( Stokes, P),
      kA_(P.get<int>("Time.NumSteps") != 0 ? 1.0/P.get<double>("Time.StepSize") : 0.0), // P.get<int>("Time.NumSteps") == 0: stat. problem
      kM_(P.get<double>("Stokes.Theta")),
      ssor_( P.get<double>("Stokes.Omega")), PCGsolver_( ssor_, P.get<int>("Stokes.InnerIter"), P.get<double>("Stokes.InnerTol")), CGsolver_( P.get<int>("Stokes.InnerIter"), P.get<double>("Stokes.InnerTol")),
      PCGsgssolver_(sgs_, P.get<int>("Stokes.InnerIter"), P.get<double>("Stokes.InnerTol")),
      q_(), minressolver_( q_, P.get<int>("Stokes.InnerIter"), P.get<double>("Stokes.OuterTol")),
      PPA_( ssor_, 8, 1e-20), PA_( PPA_), PS_( Stokes_.prM.Data.GetFinest(), Stokes_.prM.Data.GetFinest(), kA_, kM_, P.get<double>("Stokes.Omega")),
      pre_( PA_, PS_), pq_( pre_), pminressolver_( pq_, P.get<int>("Stokes.InnerIter"), P.get<double>("Stokes.OuterTol")),
      smoother_(1.0), coarsesolver_( ssorom_, 500, P.get<double>("Stokes.InnerTol")),
      MGsolver_( smoother_, coarsesolver_, P.get<int>("Stokes.InnerIter"), ( P.get<int>("Stokes.StokesMethod") == 500101)?-1:P.get<double>("Stokes.InnerTol")),
      MGpc_( MGsolver_), ismgpcp_( Stokes_.prA.Data, Stokes_.prM.Data, kA_, kM_),
      PCGMGPresolver_( ismgpcp_, P.get<int>("Stokes.OuterIter"), P.get<double>("Stokes.OuterTol")),
      ispcp_( Stokes_.prA.Data, Stokes_.prM.Data, kA_, kM_), PCGPresolver_( ispcp_, P.get<int>("Stokes.OuterIter"), P.get<double>("Stokes.OuterTol")),
      pminresMGsolver_( pqMG_, P.get<int>("Stokes.OuterIter"), P.get<double>("Stokes.OuterTol")), preMG_( apc_, ismgpcp_), apc_( MGsolver_), pqMG_( preMG_)
       {}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
StokesSolverBaseCL* StokesSolverFactoryObsoleteCL<StokesT, ProlongationVelT, ProlongationPT>::CreateStokesSolver()
{
    StokesSolverBaseCL* stokessolver = 0;
    switch (P_.template get<int>("Stokes.StokesMethod"))
    {
        case 500000 :
            stokessolver = new  UzawaSolverCL<PCG_SsorCL>( PCGsolver_, Stokes_.prM.Data.GetFinest(), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.Tau"));
            break;
        case 500100 :
            stokessolver = new  UzawaSolver2CL<PCG_SsorCL, MGSolverCL<SSORsmoothCL, PCG_SsorCL> >( PCGsolver_, MGsolver_, Stokes_.prM.Data.GetFinest(),
                P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.Tau"));
            break;
        case 501101 :
            stokessolver = new  UzawaSolver2ModifiedCL<ISMGPreCL, MGPCT>( ismgpcp_, MGpc_, Stokes_.prM.Data.GetFinest(), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"), P_.template get<double>("Stokes.Tau"));
            break;
        case 510000 :
            stokessolver = new  PSchurSolverCL<PCG_SsorCL>( PCGsolver_, Stokes_.prM.Data.GetFinest(), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        case 510100 :
            stokessolver = new  PSchurSolverCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> >( MGsolver_, Stokes_.prM.Data.GetFinest(), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol") );
            break;
        case 510001 :
            stokessolver = new  PSchurSolverCL<PCG_SgsCL>( PCGsgssolver_, Stokes_.prM.Data.GetFinest(), P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        case 511102 :
            stokessolver = new  PSchurSolver2CL<MGSolverCL<SSORsmoothCL, PCG_SsorCL>, PCGSolverCL<ISMGPreCL> >( MGsolver_, PCGMGPresolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        case 510003 :
            stokessolver = new  PSchurSolver2CL<PCG_SsorCL, PCGSolverCL<ISPreCL> >( PCGsolver_, PCGPresolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        case 511004 :
            stokessolver = new  PSchurSolver2CL<PCG_SsorCL, PCGSolverCL<ISMGPreCL> >( PCGsolver_, PCGMGPresolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        case 520000 :
            stokessolver = new  BlockMatrixSolverCL<MinresSPT>( minressolver_);
            break;
        case 530000 :
            stokessolver = new  BlockMatrixSolverCL<PMinresSP_DiagPCGT>( pminressolver_);
            break;
        case 531101 :
            stokessolver = new  BlockMatrixSolverCL<PMinresSP_FullMG>( pminresMGsolver_);
            break;
        case 540000 :
            stokessolver = new  SchurSolverCL<PCG_SsorCL>( PCGsolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        case 550000 :
            stokessolver = new  SchurNoPcSolverCL<CGSolverCL>( CGsolver_, P_.template get<int>("Stokes.OuterIter"), P_.template get<double>("Stokes.OuterTol"));
            break;
        default: throw DROPSErrCL("StokesSolverFactoryCL: Unknown Stokes solver");
    }
    return stokessolver;
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
void StokesSolverFactoryObsoleteCL<StokesT, ProlongationVelT, ProlongationPT>::
    SetMatrices( const MLMatrixCL* /*A*/, const MLMatrixCL* /*B*/, const MLMatrixCL* /*Mvel*/, const MLMatrixCL* /*M*/, const MLIdxDescCL* /*pr_idx*/) {
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ProlongationVelT* StokesSolverFactoryObsoleteCL<StokesT, ProlongationVelT, ProlongationPT>::GetPVel()
{
    return MGsolver_.GetProlongation();
}

template <class StokesT, class ProlongationVelT, class ProlongationPT>
ProlongationPT* StokesSolverFactoryObsoleteCL<StokesT, ProlongationVelT, ProlongationPT>::GetPPr()
{
    return ismgpcp_.GetProlongation();
}
#endif

} // end of namespace DROPS

