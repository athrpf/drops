/// \file
/// \brief levelset equation for two phase flow problems
/// \author Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen

#ifndef DROPS_LEVELSET_H
#define DROPS_LEVELSET_H

#include "num/spmat.h"
#include "num/discretize.h"
#include "num/solver.h"
#include "num/bndData.h"
#include "num/fe.h"
#include "levelset/mgobserve.h"
#include <vector>

namespace DROPS
{

enum SurfaceForceT
/// different types of surface forces
{
    SF_LB=0,             ///< Laplace-Beltrami discretization: \f$\mathcal O(h^{1/2})\f$
    SF_ImprovedLB=1,     ///< improved Laplace-Beltrami discretization: \f$\mathcal O(h)\f$
    SF_Const=2,          ///< surface force with constant curvature
    SF_ImprovedLBVar=3   ///< improved Laplace-Beltrami discretization with variable surface tension
};

class LevelsetP2CL
/// P2-discretization and solution of the level set equation for two phase flow problems.
{
  public:
    typedef BndDataCL<>    BndDataT;
    typedef P2EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    IdxDescCL             idx;
    VecDescCL             Phi;        ///< level set function
    instat_scalar_fun_ptr sigma;      ///< variable surface tension
    instat_vector_fun_ptr grad_sigma; ///< gradient of sigma

  private:
    MultiGridCL&        MG_;
    double              diff_,     ///< amount of diffusion in reparametrization
                        curvDiff_, ///< amount of diffusion in curvature calculation
                        SD_,       ///< streamline diffusion
                        theta_, dt_;
    MatrixCL            L_;
    BndDataT            Bnd_;
    mutable SSORPcCL    pc_;
    GMResSolverCL<SSORPcCL>  gm_;
    SurfaceForceT       SF_;

    void SetupReparamSystem( MatrixCL&, MatrixCL&, const VectorCL&, VectorCL&) const;
    void SetupSmoothSystem ( MatrixCL&, MatrixCL&)                             const;
    void SmoothPhi( VectorCL& SmPhi, double diff)                              const;

  public:
    MatrixCL            E, H;

    LevelsetP2CL( MultiGridCL& mg, instat_scalar_fun_ptr sig= 0,instat_vector_fun_ptr gsig= 0,
        double theta= 0.5, double SD= 0., double diff= 0., int iter= 1000, double tol= 1e-7,
        double curvDiff= -1.)
    : idx( P2_FE), sigma( sig), grad_sigma( gsig), MG_( mg), diff_(diff), curvDiff_( curvDiff), SD_( SD),
        theta_( theta), dt_( 0.), Bnd_( BndDataT(mg.GetBnd().GetNumBndSeg())),
        gm_( pc_, 100, iter, tol), SF_( SF_ImprovedLB)
    {}

    LevelsetP2CL( MultiGridCL& mg, const BndDataT& bnd, instat_scalar_fun_ptr sig= 0,
        instat_vector_fun_ptr gsig= 0, double theta= 0.5, double SD= 0, double diff= 0,
        Uint iter= 1000, double tol= 1e-7, double curvDiff= -1)
    : idx( P2_FE), sigma( sig), grad_sigma( gsig), MG_( mg), diff_(diff), curvDiff_( curvDiff), SD_( SD),
        theta_( theta), dt_( 0.), Bnd_( bnd), gm_( pc_, 100, iter, tol), SF_(SF_ImprovedLB)
    {}

    const MultiGridCL& GetMG() const { return MG_; }

    GMResSolverCL<SSORPcCL>& GetSolver() { return gm_; }

    const BndDataT& GetBndData() const { return Bnd_; }

    /// \name Numbering
    ///@{
    void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, MG_, Bnd_, match); }
    void DeleteNumbering( IdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    ///@}

    /// initialize level set function
    void Init( scalar_fun_ptr);

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetTimeStep( double dt, double theta=-1);
    /// \remarks call SetupSystem \em before calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&);
    /// perform one time step
    void DoStep();
    /// Reparametrization by solving evolution equation (not recommended).
    void Reparam( Uint steps, double dt);
    /// Reparametrization by Fast Marching method (recommended).
    void ReparamFastMarching( bool ModifyZero= true, bool Periodic= false, bool OnlyZeroLvl= false);

    /// tests whether level set function changes its sign on tetra \p t.
    bool   Intersects( const TetraCL&) const;
    /// returns information about level set function and interface.
    template<class DiscVelSolT>
    void   GetInfo( double& maxGradPhi, double& Volume, Point3DCL& bary, Point3DCL& vel, const DiscVelSolT& vel_sol, Point3DCL& minCoord, Point3DCL& maxCoord) const;
    /// returns the maximum and minimum of the gradient of phi
    void   GetMaxMinGradPhi(double& maxGradPhi, double& minGradPhi) const;
    /// returns approximate volume of domain where level set function is negative.
    double GetVolume( double translation= 0) const;
    /// volume correction to ensure no loss or gain of mass.
    double AdjustVolume( double vol, double tol, double surf= 0) const;
    /// Set type of surface force.
    void   SetSurfaceForce( SurfaceForceT SF) { SF_= SF; }
    /// Get type of surface force.
    SurfaceForceT GetSurfaceForce() const { return SF_; }
    /// Discretize surface force
    void   AccumulateBndIntegral( VecDescCL& f) const;
    /// Set surface tension and its gradient.
    void   SetSigma( instat_scalar_fun_ptr sig, instat_vector_fun_ptr gsig= 0) { sigma= sig; grad_sigma= gsig; }
    /// Clear all matrices, should be called after grid change to avoid reuse of matrix pattern
    void   ClearMat() { E.clear(); H.clear(); L_.clear(); }
    /// \name Evaluate Solution
    ///@{
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &Phi, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& MyPhi) const
        { return const_DiscSolCL( &MyPhi, &Bnd_, &MG_); }
    ///@}

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the coupling navstokes-levelset. They should not be called by a common user.
    /// Use LevelsetP2CL::DoStep() instead.
    ///@{
    void ComputeRhs( VectorCL&) const;
    const MatrixCL& GetL() const { return L_; }
    void DoStep    ( const VectorCL&);
    void DoLinStep ( const VectorCL&);
    ///@}
};


/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the Function ls.Phi.
///
/// The actual work is done in post_refine().
class LevelsetRepairCL : public MGObserverCL
{
  private:
    LevelsetP2CL& ls_;

  public:
    LevelsetRepairCL (LevelsetP2CL& ls)
        : ls_( ls) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
};


/// marks all tetrahedra in the band |\p DistFct(x)| < \p width for refinement
void MarkInterface (scalar_fun_ptr DistFct, double width, MultiGridCL&);
/// marks all tetrahedra in the band |\p lset(x)| < \p width for refinement
void MarkInterface ( const LevelsetP2CL::const_DiscSolCL& lset, double width, MultiGridCL& mg);

} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

