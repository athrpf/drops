//**************************************************************************
// File:    levelset.h                                                     *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

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
{
  SF_CSF=0, SF_Const=1
};
  
class LevelsetP2CL
// P2-discretization and solution of the levelset equation for two phase
// flow problems.
{
  public:
    typedef BndDataCL<>    BndDataT;
    typedef P2EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    IdxDescCL           idx;
    VecDescCL           Phi;
    double              sigma;     // surface tension

  private:
    MultiGridCL&        MG_;
    double              diff_,     // amount of diffusion in reparametrization
                        curvDiff_, // amount of diffusion in curvature calculation
                        SD_,       // streamline diffusion
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
    LevelsetP2CL( MultiGridCL& mg, double sig= 0, double theta= 0.5, double SD= 0, 
                  double diff= 0, Uint iter=1000, double tol=1e-7, double curvDiff= -1)
      : idx( 1, 1), sigma( sig), MG_( mg), diff_(diff), curvDiff_( curvDiff), SD_( SD), 
        theta_( theta), dt_( 0.), Bnd_( BndDataT(mg.GetBnd().GetNumBndSeg()) ), 
        gm_( pc_, 100, iter, tol), SF_(SF_CSF)
    {}
    
    LevelsetP2CL( MultiGridCL& mg, const BndDataT& bnd, double sig= 0, double theta= 0.5, double SD= 0, 
                  double diff= 0, Uint iter=1000, double tol=1e-7, double curvDiff= -1)
      : idx( 1, 1), sigma( sig), MG_( mg), diff_(diff), curvDiff_( curvDiff), SD_( SD), 
        theta_( theta), dt_( 0.), Bnd_( bnd), gm_( pc_, 100, iter, tol), SF_(SF_CSF)
    {}
    
    GMResSolverCL<SSORPcCL>& GetSolver() { return gm_; }
    
    const BndDataT& GetBndData() const { return Bnd_; }
    
    void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, MG_, Bnd_, match); }
    void DeleteNumbering( IdxDescCL* idx)
        { DeleteNumb( *idx, MG_); }

    void Init( scalar_fun_ptr);

    void SetTimeStep( double dt, double theta=-1);
    // call SetupSystem *before* calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem( const DiscVelSolT&);
    void DoStep();
    void Reparam( Uint steps, double dt);
    void ReparamFastMarching( bool ModifyZero= true, bool Periodic= false, bool OnlyZeroLvl= false);
    
    bool   Intersects( const TetraCL&) const;
    double GetVolume( double translation= 0) const;
    double AdjustVolume( double vol, double tol, double surf= 0) const;
    void   SetSurfaceForce( SurfaceForceT SF) { SF_= SF; }
    void   AccumulateBndIntegral( VecDescCL& f) const;
    
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &Phi, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& MyPhi) const
        { return const_DiscSolCL( &MyPhi, &Bnd_, &MG_); }
        
    // the following member functions are added to enable an easier implementation
    // of the coupling navstokes-levelset. They should not be called by a common user.
    void ComputeRhs( VectorCL&) const;
    void DoStep    ( const VectorCL&);
};


class InterfacePatchCL
/// computes the planar interface patches, which are the intersection of a child T' of
/// a tetrahedron T and the zero level of I(phi), where I(phi) is the linear interpolation
/// of the level set function phi on T'.
{
  private:
    static const double approxZero_;
    const RefRuleCL RegRef_;
    int             sign_[10], num_sign_[3];  // 0/1/2 = -/0/+
    int             intersec_, ch_, Edge_[4];
    double          PhiLoc_[10], sqrtDetATA_;
    Point3DCL       PQRS_[4], Coord_[10], B_[3];
    BaryCoordCL     Bary_[4], BaryDoF_[10];
    Point2DCL       ab_;
  
    inline void Solve2x2( const double det, const SMatrixCL<2,2>& A, SVectorCL<2>& x, const SVectorCL<2>& b)
    { x[0]= (A(1,1)*b[0]-A(0,1)*b[1])/det;    x[1]= (A(0,0)*b[1]-A(1,0)*b[0])/det; }

  public:
    InterfacePatchCL();
    
    static int Sign( double phi) { return std::abs(phi)<approxZero_ ? 0 : (phi>0 ? 1 : -1); } //< returns -1/0/1
    
    void Init( const TetraCL& t, const VecDescCL& ls);

    // Remark: The following functions are only valid, if Init(...) was called before!
    int    GetSign( Uint DoF)   const { return sign_[DoF]; }   //< returns -1/0/1
    double GetPhi( Uint DoF)    const { return PhiLoc_[DoF]; } //< returns value of level set function
    bool   Intersects()         const                          //< returns wether patch exists (i.e. interface intersects tetra)
      { for(int i=1; i<10; ++i) if (sign_[0]!=sign_[i]) return true; return false; }
    bool   IntersectsInterior() const                          //< returns wether patch exists, which is not subset of a face
      { for(int i=0; i<9; ++i) for (int j=i+1; j<10; ++j) if (sign_[i]*sign_[j]==-1) return true; return false; }
    bool   ComputeForChild( Uint ch);                          //< returns true, if a patch exists for this child
    bool   ComputeCutForChild( Uint ch);                       //< returns true, if a patch exists for this child

    // Remark: The following functions are only valid, if ComputeForChild(...) was called before!
    bool               IsQuadrilateral()     const { return intersec_==4; }
    bool               EqualToFace()         const { return num_sign_[1]>=3; }   //< returns true, if patch is shared by two tetras
    Uint               GetNumPoints()        const { return intersec_; }
    const Point3DCL&   GetPoint( Uint i)     const { return PQRS_[i]; }
    const BaryCoordCL& GetBary ( Uint i)     const { return Bary_[i]; }
    int                GetNumSign( int sign) const { return num_sign_[sign+1]; } //< returns number of child points with given sign, where sign is in {-1, 0, 1}
    double             GetFuncDet()          const { return sqrtDetATA_; }
    double             GetAreaFrac()         const { return intersec_==4 ? ab_[0]+ab_[1]-1 : 0; }
    const Point3DCL&   GetGradId( Uint i)    const { return B_[i]; }

    void               WriteGeom( std::ostream&) const;                          //< Geomview output for debugging
    void               DebugInfo( std::ostream&, bool InfoOnChild= false) const;
    
    // Remark: The following functions are only valid, if ComputeCutForChild(...) was called before!
    template<class ValueT>
    ValueT quad( const LocalP2CL<ValueT>&, double absdet, bool posPart= true);   //< integrate on pos./neg. part
};

} // end of namespace DROPS

#include "levelset/levelset.tpp"

#endif

