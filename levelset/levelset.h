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
    MatrixCL            E_, H_, L_;
    BndDataT            Bnd_;
    SSORPcCL            pc_;
    GMResSolverCL<SSORPcCL>  gm_;
    SurfaceForceT       SF_;

    void SetupReparamSystem( MatrixCL&, MatrixCL&, const VectorCL&, VectorCL&) const;
    void SetupSmoothSystem ( MatrixCL&, MatrixCL&)                             const;
    void SmoothPhi( VectorCL& SmPhi, double diff)                              const;

  public:
    LevelsetP2CL( MultiGridCL& mg, instat_scalar_fun_ptr sig= 0,instat_vector_fun_ptr gsig= 0,
        double theta= 0.5, double SD= 0., double diff= 0., int iter= 1000, double tol= 1e-7,
        double curvDiff= -1.)
    : idx( 1, 1), sigma( sig), grad_sigma( gsig), MG_( mg), diff_(diff), curvDiff_( curvDiff), SD_( SD),
        theta_( theta), dt_( 0.), Bnd_( BndDataT(mg.GetBnd().GetNumBndSeg())),
        gm_( pc_, 100, iter, tol), SF_( SF_ImprovedLB)
    {}

    LevelsetP2CL( MultiGridCL& mg, const BndDataT& bnd, instat_scalar_fun_ptr sig= 0,
        instat_vector_fun_ptr gsig= 0, double theta= 0.5, double SD= 0, double diff= 0,
        Uint iter= 1000, double tol= 1e-7, double curvDiff= -1)
    : idx( 1, 1), sigma( sig), grad_sigma( gsig), MG_( mg), diff_(diff), curvDiff_( curvDiff), SD_( SD),
        theta_( theta), dt_( 0.), Bnd_( bnd), gm_( pc_, 100, iter, tol), SF_(SF_ImprovedLB)
    {}

    const MultiGridCL& GetMG() const { return MG_; }

    GMResSolverCL<SSORPcCL>& GetSolver() { return gm_; }

    const BndDataT& GetBndData() const { return Bnd_; }

    /// \name Numbering
    ///@{
    void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, MG_, Bnd_, match); }
    void DeleteNumbering( IdxDescCL* idx)
        { DeleteNumb( *idx, MG_); }
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
    void   GetInfo( double& maxGradPhi, double& Volume, Point3DCL& bary, Point3DCL& minCoord, Point3DCL& maxCoord) const;
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
    void   ClearMat() { E_.clear(); H_.clear(); L_.clear(); }
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
    void DoStep    ( const VectorCL&);
    ///@}
};


class InterfacePatchCL
/// Computes approximation of interface.
/** Computes the planar interface patches, which are the intersection of a child T' of
 *  a tetrahedron T and the zero level of I(phi), where I(phi) is the linear interpolation
 *  of the level set function phi on T'.
 */
{
  private:
    static const double approxZero_;
    const RefRuleCL RegRef_;
    int             sign_[10], num_sign_[3];  // 0/1/2 = -/0/+
    int             intersec_, ch_, Edge_[4];
    int             numchildtriangles_; // The number of triangles in the intersection of a child with the interface.
    double          sqrtDetATA_;
    LocalP2CL<>     PhiLoc_;
    Point3DCL       PQRS_[4], Coord_[10], B_[3];
    BaryCoordCL     Bary_[4], BaryDoF_[10];
    Point2DCL       ab_;

    inline void Solve2x2( const double det, const SMatrixCL<2,2>& A, SVectorCL<2>& x, const SVectorCL<2>& b)
    { x[0]= (A(1,1)*b[0]-A(0,1)*b[1])/det;    x[1]= (A(0,0)*b[1]-A(1,0)*b[0])/det; }

  public:
    InterfacePatchCL();

    static int Sign( double phi) { return std::abs(phi)<approxZero_ ? 0 : (phi>0 ? 1 : -1); } ///< returns -1/0/1

    void Init( const TetraCL& t, const VecDescCL& ls, double translation= 0.);

    /// \name Use after Init
    /// \remarks The following functions are only valid, if Init(...) was called before!
    ///@{
    int    GetSign( Uint DoF)   const { return sign_[DoF]; }   ///< returns -1/0/1
    double GetPhi( Uint DoF)    const { return PhiLoc_[DoF]; } ///< returns value of level set function
    bool   Intersects()         const                          ///  returns wether patch exists (i.e. interface intersects tetra)
      { for(int i=1; i<10; ++i) if (sign_[0]!=sign_[i]) return true; return false; }
    bool   IntersectsInterior() const                          ///  returns wether patch exists, which is not subset of a face
      { for(int i=0; i<9; ++i) for (int j=i+1; j<10; ++j) if (sign_[i]*sign_[j]==-1) return true; return false; }
    bool   ComputeForChild( Uint ch);                          ///< returns true, if a patch exists for this child
    bool   ComputeCutForChild( Uint ch);                       ///< returns true, if a patch exists for this child
    ///@}

    /// \name Use after ComputeForChild
    /// \remarks The following functions are only valid, if ComputeForChild(...) was called before!
    ///@{
    int                GetNumTriangles()     const { return numchildtriangles_; } ///< Returns, how many triangles form the intersection of the child and the interface.
    bool               IsQuadrilateral()     const { return intersec_==4; }
    bool               EqualToFace()         const { return num_sign_[1]>=3; }   ///< returns true, if patch is shared by two tetras
    Uint               GetNumPoints()        const { return intersec_; }
    const Point3DCL&   GetPoint( Uint i)     const { return PQRS_[i]; }
    const BaryCoordCL& GetBary ( Uint i)     const { return Bary_[i]; } ///< The first three points are the vertices of the triangular patch; if the patch is quadrilateral, the last three points are the vertices of the second triangle.
    int                GetNumSign ( int sign) const { return num_sign_[sign+1]; } ///< returns number of child points with given sign, where sign is in {-1, 0, 1}
    double             GetFuncDet( int i= 0)  const { return sqrtDetATA_*(i==0 ? 1.0 : GetAreaFrac()); } ///< Returns the Determinant for surface integration for triangle i.
    double             GetAreaFrac()         const { return intersec_==4 ? ab_[0]+ab_[1]-1 : 0; }
    const Point3DCL&   GetGradId( Uint i)    const { return B_[i]; }
    const Point3DCL    ApplyProj( const Point3DCL& grad) const { return grad[0]*B_[0] + grad[1]*B_[1] + grad[2]*B_[2]; }

    void               WriteGeom( std::ostream&) const;                          ///< Geomview output for debugging
    void               DebugInfo( std::ostream&, bool InfoOnChild= false) const;
    ///@}

    /// \name Use after ComputeCutForChild
    /// \remarks The following functions are only valid, if ComputeCutForChild(...) was called before!
    ///@{
    template<class ValueT>
    ValueT quad( const LocalP2CL<ValueT>&, double absdet, bool posPart= true);   ///< integrate on pos./neg. part
    ///@}
};

/// marks all tetrahedra in the band |\p DistFct(x)| < \p width for refinement
void MarkInterface (scalar_fun_ptr DistFct, double width, MultiGridCL&);
/// marks all tetrahedra in the band |\p lset(x)| < \p width for refinement
void MarkInterface ( const LevelsetP2CL::const_DiscSolCL& lset, double width, MultiGridCL& mg);

} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

