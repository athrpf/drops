//**************************************************************************
// File:    bndData.h                                                      *
// Content: boundary data                                                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

/// \file
/// \brief Classes for storing and handling boundary data.

/// \todo All BndVal-functions should have the signature BndValT func( const Point3DCL&, double).
/// Up to now this is not the case for the Poisson problems.

#ifndef DROPS_BNDDATA_H
#define DROPS_BNDDATA_H

#include "misc/utils.h"
#include "geom/multigrid.h"
#include "geom/boundary.h"

namespace DROPS
{

enum BndCondT
/// \brief enum for boundary conditions in DROPS.
///
/// Conventions:
/// - boundary condition with lower number is superior
/// - odd number for boundary with unknowns
/// - valid boundary conditions have numbers in the range 0..99
/// - interior simplices have no boundary data and return NoBC
///   in BndDataCL::GetBC
{
    Dir0BC= 0,                   ///< hom.   Dirichlet boundary conditions
    DirBC= 2,                    ///< inhom. Dirichlet boundary conditions
    Per1BC= 13,                  ///< periodic boundary conditions, where
    Per2BC= 11,                  ///< Per1BC and Per2BC denote corresponding boundaries
    Nat0BC= 21,                  ///< hom.   natural   boundary condition
    NatBC= 23,                   ///< inhom. natural   boundary conditions
    OutflowBC= 21,               ///< same as Nat0BC, for convenience
    WallBC= 0,                   ///< same as Dir0BC, for convenience

    NoBC= -1,                    ///< interior simplices
    UndefinedBC_= 99,            ///< ReadMeshBuilderCL: error, unknown bc
    MaxBC_= 100                  ///< upper bound for valid bc's
};


/// \brief Represents the boundary data of a single boundary segment for a certain variable.
///
/// This class stores the boundary condition and a function representing
/// the corresponding boundary values. It is only used as a part of the BndDataCL.
template<class BndValT = double>
class BndSegDataCL
{
  public:
    typedef BndValT (*bnd_val_fun)( const Point3DCL&, double);

  private:
    BndCondT     bc_;
    bnd_val_fun  bnd_val_;

  public:
    BndSegDataCL( BndCondT bc= Nat0BC, bnd_val_fun f= 0)
      : bc_(bc), bnd_val_(f)
    /// \param bc boundary condition for one segment
    /// \param f  boundary value function on this segment. f has to be specified if bc is non-homogeneous
    {
        if ( bc!=Nat0BC && bc!=Dir0BC && bc!=Per1BC && bc!=Per2BC && f==0)
            throw DROPSErrCL("BndSegDataCL: no boundary function for non-homogeneous condition specified!");
    }

    bool        WithUnknowns() const { return bc_ & 1; }
    bool        IsDirichlet()  const { return bc_==Dir0BC || bc_==DirBC; }
    bool        IsNatural()    const { return bc_==Nat0BC || bc_==NatBC; }
    bool        IsPeriodic()   const { return bc_==Per1BC || bc_==Per2BC ; }
    BndCondT    GetBC()        const { return bc_; }
    bnd_val_fun GetBndFun()    const { return bnd_val_; }
    BndValT     GetBndVal( const Point3DCL& p, double t= 0.0) const { return bnd_val_( p, t); }
};


/// \brief Contains the boundary data of all boundary segments for a certain variable.
///
/// For each sub-simplex on the boundary, the boundary condition and corresponding boundary values
/// can be accessed.
template<class BndValT = double>
class BndDataCL
{
  public:
    typedef BndSegDataCL<BndValT>             BndSegT;
    typedef BndValT                           bnd_type;
    typedef typename BndSegT::bnd_val_fun     bnd_val_fun;

  private:
    std::vector<BndSegT> BndData_;

  public:
    /// If \a bc and \a fun are given, they are assumed to be arrays of length \a numbndseg
    /// containing the boundary conditions and boundary data resp. of the boundary segments.
    /// If \a bc is omitted, hom. natural boundary conditions are imposed (Nat0BC) for all boundary segments.
    /// \a fun should only be omitted, if all boundary conditions given are homogenious
    /// and thus no boundary values have to be specified.
    BndDataCL( BndIdxT numbndseg, const BndCondT* bc= 0, const bnd_val_fun* fun= 0);
    /// Deprecated ctor, just for compatibility with older code
    BndDataCL( BndIdxT numbndseg, const bool* isneumann, const bnd_val_fun* fun); // deprecated ctor!

    /// \name boundary condition
    /// Returns superior boundary condition of sub-simplex
    /// \{
    inline BndCondT GetBC( const VertexCL&)    const;
    inline BndCondT GetBC( const EdgeCL&)      const;
    inline BndCondT GetBC( const FaceCL&)      const;
    inline BndCondT GetBC( const VertexCL&, BndIdxT&) const;
    inline BndCondT GetBC( const EdgeCL&, BndIdxT&)   const;
    inline BndCondT GetBC( const FaceCL&, BndIdxT&)   const;
    /// \}

    inline bool IsOnDirBnd( const VertexCL&) const;
    inline bool IsOnDirBnd( const EdgeCL&)   const;
    inline bool IsOnDirBnd( const FaceCL&)   const;
    inline bool IsOnNatBnd( const VertexCL&) const;
    inline bool IsOnNatBnd( const EdgeCL&)   const;
    inline bool IsOnNatBnd( const FaceCL&)   const;
    inline bool IsOnPerBnd( const VertexCL&) const;
    inline bool IsOnPerBnd( const EdgeCL&)   const;
    inline bool IsOnPerBnd( const FaceCL&)   const;

    /// \name boundary value
    /// Returns boundary value of sub-simplex
    /// \{
    inline BndValT GetDirBndValue( const VertexCL&, double) const;
    inline BndValT GetDirBndValue( const EdgeCL&, double)   const;
    inline BndValT GetDirBndValue( const FaceCL&, double)   const;
    inline BndValT GetNatBndValue( const FaceCL&, double)   const;

    inline BndValT GetDirBndValue( const VertexCL& v) const { return GetDirBndValue( v, 0); }
    inline BndValT GetDirBndValue( const EdgeCL& e)   const { return GetDirBndValue( e, 0); }
    inline BndValT GetDirBndValue( const FaceCL& f)   const { return GetDirBndValue( f, 0); }
    /// \}

    /// \name boundary segment data
    /// Returns boundary segment data with superior boundary condition of sub-simplex
    /// \{
    inline BndSegT GetBndSeg( const VertexCL&) const;
    inline BndSegT GetBndSeg( const EdgeCL&)   const;
    inline BndSegT GetBndSeg( const FaceCL&)   const;
           BndSegT GetBndSeg( BndIdxT idx)     const { return BndData_[idx]; }
    bnd_val_fun    GetBndFun( BndIdxT idx)     const { return BndData_[idx].GetBndFun(); }
    /// \}
};


template<class BndValT= double>
class NoBndDataCL
{
  public:
    typedef BndSegDataCL<BndValT>             BndSegT;
    typedef BndValT                           bnd_type;
    typedef typename BndSegT::bnd_val_fun     bnd_val_fun;

    // default ctor, dtor, whatever

    template<class SimplexT>
    static inline BndCondT GetBC  (const SimplexT&) { return NoBC; }
    template<class SimplexT>
    static inline bool IsOnDirBnd (const SimplexT&) { return false; }
    template<class SimplexT>
    static inline bool IsOnNatBnd (const SimplexT&) { return false; }
    template<class SimplexT>
    static inline bool IsOnPerBnd (const SimplexT&) { return false; }

    static inline BndValT GetDirBndValue (const VertexCL&, double= 0.0)
        { throw DROPSErrCL("NoBndDataCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
    static inline BndValT GetDirBndValue (const EdgeCL&, double= 0.0)
        { throw DROPSErrCL("NoBndDataCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on edge."); }
    static inline BndValT GetDirBndValue (const FaceCL&, double= 0.0)
        { throw DROPSErrCL("NoBndDataCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on face."); }
};


//======================================
//        inline functions
//======================================

/// \brief returns zero (commonly used boundary data function)
///
/// Commonly used boundary data function of type BndSegDataCL<double>::bnd_val_fun
/// returning zero
inline double Zero( const Point3DCL&, double) { return 0.; }

/// \brief returns zero vector (commonly used boundary data function)
///
/// Commonly used boundary data function of type BndSegDataCL<Point3DCL>::bnd_val_fun
/// returning zero velocities 
inline Point3DCL ZeroVel( const Point3DCL&, double) { return Point3DCL(0.); }


template<class BndValT>
inline BndDataCL<BndValT>::BndDataCL( BndIdxT numbndseg, const BndCondT* bc, const bnd_val_fun* fun)
{
    BndData_.reserve( numbndseg);
    for (Uint i=0; i<numbndseg; ++i)
        BndData_.push_back( BndSegT( bc ? bc[i] : Nat0BC, fun ? fun[i] : 0) );
}

template<class BndValT>
inline BndDataCL<BndValT>::BndDataCL( BndIdxT numbndseg, const bool* isneumann, const bnd_val_fun* fun) // deprecated ctor!
{
    BndData_.reserve( numbndseg);
    for (Uint i=0; i<numbndseg; ++i)
        BndData_.push_back( BndSegT( isneumann[i] ? NatBC : DirBC, fun[i]) );
}


//---------------------------------------
// functions for vertices

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnDirBnd( const VertexCL& v) const
{ // v is on dir bnd, iff it is on one or more dir bnd segments
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( BndData_[it->GetBndIdx()].IsDirichlet() )
            return true;
    return false;
}

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnNatBnd( const VertexCL& v) const
{ // v is on neu bnd, iff it is only on neu bnd segments
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !BndData_[it->GetBndIdx()].IsNatural() )
            return false;
    return true;
}

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnPerBnd( const VertexCL& v) const
{ // v is on per bnd, iff it is on one or more per bnd segments and not on a dir bnd segment
    if ( !v.IsOnBoundary() ) return false;
    bool HasPer= false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndData_[it->GetBndIdx()].GetBC();
        if (bc==DirBC || bc==Dir0BC)
            return false;
        if (bc==Per1BC || bc==Per2BC)
            HasPer= true;
    }
    return HasPer;
}

template<class BndValT>
inline BndCondT BndDataCL<BndValT>::GetBC( const VertexCL& v) const
{ /// Returns BC on v with lowest number (i.e. the superior BC on v)
    if ( !v.IsOnBoundary() ) return NoBC;
    BndCondT bc= MaxBC_;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        bc= std::min( bc, BndData_[it->GetBndIdx()].GetBC());
    return bc;
}

template<class BndValT>
inline BndCondT BndDataCL<BndValT>::GetBC( const VertexCL& v, BndIdxT& bidx) const
{ /// Returns BC on v with lowest number (i.e. the superior BC on v) and the number of the corresponding boundary segment
    if ( !v.IsOnBoundary() ) return NoBC;
    BndCondT bc= MaxBC_, tmp;
    BndIdxT tmp2= NoBndC;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ((tmp= BndData_[tmp2= it->GetBndIdx()].GetBC()) < bc) {
            bc= tmp;
            bidx= tmp2;
        }
    return bc;
}

template<class BndValT>
inline typename BndDataCL<BndValT>::BndSegT BndDataCL<BndValT>::GetBndSeg( const VertexCL& v) const
{ /// Returns bnd segment data on v with superior BC
    Assert( v.IsOnBoundary(), DROPSErrCL("BndDataCL::GetBndSeg(VertexCL): Not on boundary!"), ~0);
    BndCondT bc_min= MaxBC_;
    BndIdxT  idx= -1;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndData_[it->GetBndIdx()].GetBC();
        if (bc<bc_min) { bc_min= bc;    idx= it->GetBndIdx(); }
    }
    return BndData_[idx];
}

//---------------------------------------
// functions for edges

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnDirBnd( const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( BndData_[*it].IsDirichlet() )
            return true;
    return false;
}

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnNatBnd( const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !BndData_[*it].IsNatural() )
            return false;
    return true;
}

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnPerBnd( const EdgeCL& e) const
{ // e is on per bnd, iff it is on one or more per bnd segments and not on a dir bnd segment
    if ( !e.IsOnBoundary() ) return false;
    bool HasPer= false;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndData_[*it].GetBC();
        if (bc==DirBC || bc==Dir0BC)
            return false;
        if (bc==Per1BC || bc==Per2BC)
            HasPer= true;
    }
    return HasPer;
}

template<class BndValT>
inline BndCondT BndDataCL<BndValT>::GetBC( const EdgeCL& e) const
{ /// Returns BC on e with lowest number (i.e. the superior BC on e)
    if ( !e.IsOnBoundary() ) return NoBC;
    BndCondT bc= MaxBC_;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        bc= std::min( bc, BndData_[*it].GetBC());
    return bc;
}

template<class BndValT>
inline BndCondT BndDataCL<BndValT>::GetBC( const EdgeCL& e, BndIdxT& bidx) const
{ /// Returns BC on e with lowest number (i.e. the superior BC on e) and the number of the corresponding boundary segment
    if ( !e.IsOnBoundary() ) return NoBC;
    BndCondT bc= MaxBC_, tmp;
    BndIdxT tmp2= NoBndC;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ((tmp= BndData_[tmp2= *it].GetBC()) < bc) {
            bc= tmp;
            bidx= tmp2;
        }
    return bc;
}

template<class BndValT>
inline typename BndDataCL<BndValT>::BndSegT BndDataCL<BndValT>::GetBndSeg( const EdgeCL& e) const
{ /// Returns bnd segment data on e with superior BC
    Assert( e.IsOnBoundary(), DROPSErrCL("BndDataCL::GetBndSeg(EdgeCL): Not on boundary!"), ~0);
    BndCondT bc_min= MaxBC_;
    BndIdxT  idx= -1;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndData_[*it].GetBC();
        if (bc<bc_min) { bc_min= bc;    idx= *it; }
    }
    return BndData_[idx];
}

//---------------------------------------
// functions for faces

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnDirBnd( const FaceCL& f) const
{
    return f.IsOnBoundary() && BndData_[f.GetBndIdx()].IsDirichlet();
}

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnNatBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && BndData_[f.GetBndIdx()].IsNatural();
}

template<class BndValT>
inline bool BndDataCL<BndValT>::IsOnPerBnd( const FaceCL& f) const
{
    return f.IsOnBoundary() && BndData_[f.GetBndIdx()].IsPeriodic();
}

template<class BndValT>
inline BndCondT BndDataCL<BndValT>::GetBC( const FaceCL& f) const
{ // Returns BC on f
    if ( !f.IsOnBoundary() ) return NoBC;
    return BndData_[f.GetBndIdx()].GetBC();
}

template<class BndValT>
inline BndCondT BndDataCL<BndValT>::GetBC( const FaceCL& f, BndIdxT& bidx) const
{ // Returns BC and number of boundary segment on f
    if ( !f.IsOnBoundary() ) { bidx= NoBndC; return NoBC; }
    return BndData_[bidx= f.GetBndIdx()].GetBC();
}

template<class BndValT>
inline typename BndDataCL<BndValT>::BndSegT BndDataCL<BndValT>::GetBndSeg( const FaceCL& f) const
{
    Assert( f.IsOnBoundary(), DROPSErrCL("BndDataCL::GetBndSeg(FaceCL): Not on boundary!"), ~0);
    return BndData_[f.GetBndIdx()];
}


//---------------------------------------------------------
// definitions of BndDataCL<...>::GetXXXBndValue(...)

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetDirBndValue( const VertexCL& v, double t) const
/// Returns value of the Dirichlet boundary value.
/// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( BndData_[it->GetBndIdx()].IsDirichlet() )
            return BndData_[it->GetBndIdx()].GetBndFun() ? BndData_[it->GetBndIdx()].GetBndVal( v.GetCoord(), t) : BndValT();
    throw DROPSErrCL("GetDirBndValue(VertexCL): No Dirichlet Boundary Segment!");
}

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetDirBndValue( const EdgeCL& e, double t) const
/// Returns value of the Dirichlet boundary value.
/// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( BndData_[*it].IsDirichlet() )
            return BndData_[*it].GetBndFun() ? BndData_[*it].GetBndVal( GetBaryCenter(e), t) : BndValT();
    throw DROPSErrCL("GetDirBndValue(EdgeCL): No Dirichlet Boundary Segment!");
}

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetDirBndValue( const FaceCL& f, double t) const
/// Returns value of the Dirichlet boundary value.
/// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    Assert( BndData_[f.GetBndIdx()].IsDirichlet(), DROPSErrCL("GetDirBndValue(FaceCL): No Dirichlet Boundary Segment!"), ~0);
    return BndData_[f.GetBndIdx()].GetBndFun() ? BndData_[f.GetBndIdx()].GetBndVal( GetBaryCenter(f), t) : BndValT();
}

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetNatBndValue( const FaceCL& f, double t) const
/// Returns value of the Neumann boundary value.
/// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    Assert( BndData_[f.GetBndIdx()].IsNatural(), DROPSErrCL("GetNeuBndValue(FaceCL): No Neumann Boundary Segment!"), ~0);
    return BndData_[f.GetBndIdx()].GetBndVal( GetBaryCenter(f), t);
}


} //end of namespace DROPS

#endif
