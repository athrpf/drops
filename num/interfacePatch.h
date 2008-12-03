/// \file interfacePatch.h
/// \brief Computes 2D patches and 3D cuts of tetrahedra and interface
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

#ifndef DROPS_INTERFACEPATCH_H
#define DROPS_INTERFACEPATCH_H

#include "num/discretize.h"

namespace DROPS
{

class InterfacePatchCL
/// Computes approximation of interface.
/** Computes the planar interface patches, which are the intersection of a child T' of
 *  a tetrahedron T and the zero level of I(phi), where I(phi) is the linear interpolation
 *  of the level set function phi on T'. With LinearEdgeIntersection==false the planar interface
 *  patch on T' is determined by computing the roots of the quadratic level set function phi
 *  on each edge of T' where phi changes its sign.
 */
{
  public:
    typedef SArrayCL<BaryCoordCL,4> SubTetraT;

  private:
    static const double approxZero_;
    static const bool   LinearEdgeIntersection= true;
    const RefRuleCL RegRef_;
    int             sign_[10], num_sign_[3];  // 0/1/2 = -/0/+
    int             intersec_, ch_, Edge_[4], innersec_; // intersec_: # of all intersections, innersec_ = # of edge intersections
    int             numchildtriangles_; // The number of triangles in the intersection of a child with the interface.
    double          sqrtDetATA_;
    LocalP2CL<>     PhiLoc_;
    Point3DCL       PQRS_[4], Coord_[10], B_[3];
    BaryCoordCL     Bary_[4];
    static BaryCoordCL BaryDoF_[10], AllEdgeBaryCenter_[10][10];
    Point2DCL       ab_;
    std::vector<SubTetraT> posTetras, negTetras;

    inline void Solve2x2( const double det, const SMatrixCL<2,2>& A, SVectorCL<2>& x, const SVectorCL<2>& b)
        { x[0]= (A(1,1)*b[0]-A(0,1)*b[1])/det;    x[1]= (A(0,0)*b[1]-A(1,0)*b[0])/det; }
    void InsertSubTetra(SubTetraT& BaryCoords, bool pos);

  public:
    InterfacePatchCL();

    static int Sign( double phi) { return std::abs(phi)<approxZero_ ? 0 : (phi>0 ? 1 : -1); } ///< returns -1/0/1

    inline static double EdgeIntersection (Uint v0, Uint v1, LocalP2CL<>& philoc); ///< Compute the root of the LS-Function restricted to the edge (v0,v1) as barycentric coordinate on this edge.

    void Init( const TetraCL& t, const VecDescCL& ls, double translation= 0.);
    void Init( const TetraCL& t, const LocalP2CL<double>& ls, double translation= 0.);

    /// \name Use after Init
    /// \remarks The following functions are only valid, if Init(...) was called before!
    ///@{
    int    GetSign( Uint DoF)   const { return sign_[DoF]; }   ///< returns -1/0/1
    double GetPhi( Uint DoF)    const { return PhiLoc_[DoF]; } ///< returns value of level set function
    bool   Intersects()         const                          ///  returns whether patch exists (i.e. interface intersects tetra)
      { for(int i=1; i<10; ++i) if (sign_[0]!=sign_[i]) return true; return false; }
    bool   IntersectsInterior() const                          ///  returns whether patch exists, which is not subset of a face
      { for(int i=0; i<9; ++i) for (int j=i+1; j<10; ++j) if (sign_[i]*sign_[j]==-1) return true; return false; }
    bool   ComputeForChild( Uint ch);                          ///< returns true, if a patch exists for this child
    bool   ComputeCutForChild( Uint ch);                       ///< returns true, if a patch exists for this child
    void   ComputeSubTets();                                   ///< Computes a tetrahedralization of \f$\{\varphi<0\}\cap T\f$ and \f$\{\varphi>0\}\cap T\f$; the regular children of T are triangulated.
    ///@}

    /// \name Use after ComputeSubTets
    /// \remarks The following functions are only valid, if ComputeSubTets() was called before!
    ///@{
    const SubTetraT& GetTetra (Uint i)  const { return i < negTetras.size() ? negTetras[i] : posTetras[i-negTetras.size()];}
    Uint  GetNumTetra()         const {return negTetras.size() + posTetras.size();} ///< returns number of subtetras
    Uint  GetNumNegTetra()      const {return negTetras.size();}      ///< returns number of tetras with level set function < 0
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
    int                GetNumSign ( int sign)const { return num_sign_[sign+1]; } ///< returns number of child points with given sign, where sign is in {-1, 0, 1}
    double             GetFuncDet( Uint tri= 0) const { return sqrtDetATA_*(tri==0 ? 1.0 : GetAreaFrac()); } ///< Returns the Determinant for surface integration for triangle \p tri.
    double             GetAreaFrac()         const { return intersec_==4 ? ab_[0]+ab_[1]-1 : 0; }
    template<class ValueT>
    ValueT quad2D( const LocalP2CL<ValueT>&, Uint tri= 0) const;  ///< integrate on triangle \p tri, quadrature exact up to degree 2
    const Point3DCL&   GetGradId( Uint i)    const { return B_[i]; }
          Point3DCL    GetNormal ()          const; ///< Returns the unit normal to the linear approximation of \f$\Gamma\f$, that points from \f$\{\varphi<0\}\f$ to \f$\{\varphi<0\}\f$.
          Point3DCL    ApplyProj( const Point3DCL& grad) const { return grad[0]*B_[0] + grad[1]*B_[1] + grad[2]*B_[2]; }

    void               WriteGeom( std::ostream&) const;                          ///< Geomview output for debugging
    void               DebugInfo( std::ostream&, bool InfoOnChild= false) const;
    ///@}

    /// \name Use after ComputeCutForChild
    /// \remarks The following functions are only valid, if ComputeCutForChild(...) was called before!
    ///@{
    template<class ValueT>
    ValueT quad( const LocalP2CL<ValueT>&, double absdet, bool posPart= true);   ///< integrate on pos./neg. part
    template<class ValueT>
    void   quadBothParts( ValueT& int_pos, ValueT& int_neg, const LocalP2CL<ValueT>&, double absdet);   ///< integrate on pos. and neg. part
    ///@}
};

/// \brief Returns the length of the longest and shortest edge, which is cut by the interface.
///
/// \param ls VecDescCL that describes the levelset-function.
/// \param begin begin-iterator of triangulation-edges
/// \param end end-iterator of triangulation-edges
template <class It>
  std::pair<double, double>
  h_interface (It begin, It end, const VecDescCL& ls);

} // end of namespace DROPS

#include "num/interfacePatch.tpp"

#endif
