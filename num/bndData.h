/// \file bndData.h
/// \brief Classes for storing and handling boundary data.
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_BNDDATA_H
#define DROPS_BNDDATA_H

#include "misc/utils.h"
#include "geom/multigrid.h"
#include "geom/boundary.h"

//#include <boost/function.hpp>

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

    NoBC= 98,                    ///< interior simplices
    UndefinedBC_= 99,            ///< ReadMeshBuilderCL: error, unknown bc
    MaxBC_= 100                  ///< upper bound for valid bc's
};


/// \brief Stores the boundary condition of a single boundary segment for a certain variable.
class BndCondInfoCL
{
  protected:
    BndCondT bc_;

  public:
    BndCondInfoCL( BndCondT bc= Nat0BC)
    /// \param bc boundary condition for one segment
      : bc_(bc) {}

    bool     WithUnknowns() const { return bc_ & 1; }
    bool     IsDirichlet()  const { return bc_==Dir0BC || bc_==DirBC; }
    bool     IsNatural()    const { return bc_==Nat0BC || bc_==NatBC; }
    bool     IsPeriodic()   const { return bc_==Per1BC || bc_==Per2BC ; }
    BndCondT GetBC()        const { return bc_; }
};

/// Prints a text-message describing the given boundary-condition.
void inline BndCondInfo (BndCondT bc, std::ostream& os)
/// \param bc Value of type BndCondT, which shall be described.
/// \param os Stream, to which the description is written.
{
    switch(bc)
    {
      case Dir0BC: /* WallBC has the same number */
                         os << "hom. Dirichlet BC / wall\n"; break;
      case DirBC:        os << "inhom. Dirichlet BC / inflow\n"; break;
      case Per1BC:       os << "periodic BC\n"; break;
      case Per2BC:       os << "periodic BC, correspondent\n"; break;
      case Nat0BC: /* OutflowBC has the same number */
                         os << "hom. Natural BC / outflow\n"; break;
      case NatBC:        os << "inhom. Natural BC\n"; break;
      case NoBC:         os << "no boundary\n"; break;
      case UndefinedBC_: os << "WARNING! unknown BC from ReadMeshBuilderCL\n"; break;
      default:           os << "WARNING! unknown BC\n";
    }
}


/// \brief Represents the boundary data of a single boundary segment for a certain variable.
///
/// This class stores the boundary condition and a function representing
/// the corresponding boundary values. It is only used as a part of the BndDataCL.
template<class BndValT = double>
class BndSegDataCL: public BndCondInfoCL
{
  public:
  typedef BndValT (*bnd_val_fun)( const Point3DCL&, double);
  
  private:
    bnd_val_fun  bnd_val_;

  public:
    BndSegDataCL( BndCondT bc= Nat0BC, bnd_val_fun f= 0)
      :  BndCondInfoCL(bc), bnd_val_(f)
    /// \param bc boundary condition for one segment
    /// \param f  boundary value function on this segment. f has to be specified if bc is non-homogeneous
    { CheckValid( bc, f); }

    static void CheckValid( BndCondT bc, bnd_val_fun f)
    /// check compatibility of boundary condition \a bc and boundary value function \a f.
    {
        if ( bc!=NoBC &&bc!=Nat0BC && bc!=Dir0BC && bc!=Per1BC && bc!=Per2BC && f==0)
           throw DROPSErrCL("BndSegDataCL: no boundary function for non-homogeneous condition specified!");
#ifdef _PAR
        if ( bc==Per1BC || bc==Per2BC )
            throw DROPSErrCL("BndSegDataCL::CheckValid: No periodic boundary conditions in parDROPS implemented, yet, sorry");
#endif
    }

    bnd_val_fun GetBndFun()    const { return bnd_val_; }
    BndValT     GetBndVal( const Point3DCL& p, double t= 0.0) const { return bnd_val_( p, t); }
};


/// \brief Contains the boundary conditions of all boundary segments for a certain variable.
///
/// For each sub-simplex on the boundary, the boundary condition can be accessed.
class BndCondCL
{
  protected:
    std::vector<BndCondInfoCL> BndCond_;

  public:
    /// If \a bc is given, it is assumed to be an array of length \a numbndseg
    /// containing the boundary conditions of the boundary segments.
    /// If \a bc is omitted, hom. natural boundary conditions are imposed (Nat0BC) for all boundary segments.
    /// For the special case \a numbndseg=0 we always have GetBC() = NoBC and IsOnXXXBnd(...) = false (aka NoBndCondCL)
    BndCondCL( BndIdxT numbndseg, const BndCondT* bc= 0)
    {
        BndCond_.resize( numbndseg);
        for (Uint i=0; i<numbndseg; ++i)
            BndCond_[i]= bc ? bc[i] : Nat0BC;
    }

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

    /// \name boundary segment
    /// Returns boundary segment with superior boundary condition of sub-simplex
    /// \{
    template<class SimplexT>
    inline const BndCondInfoCL& GetBndSeg( const SimplexT& s) const;
    const BndCondInfoCL& GetBndSeg( BndIdxT idx) const { return BndCond_[idx]; }
    /// \}
};


/// \brief Contains the boundary data of all boundary segments for a certain variable.
///
/// For each sub-simplex on the boundary, the boundary condition and corresponding boundary values
/// can be accessed.
template<class BndValT = double>
class BndDataCL: public BndCondCL
{
  public:
    typedef BndSegDataCL<BndValT>             BndSegT;
    typedef BndValT                           bnd_type;
    typedef typename BndSegT::bnd_val_fun     bnd_val_fun;

  private:
    std::vector<bnd_val_fun> BndFun_;

  public:
    /// If \a bc and \a fun are given, they are assumed to be arrays of length \a numbndseg
    /// containing the boundary conditions and boundary data resp. of the boundary segments.
    /// If \a bc is omitted, hom. natural boundary conditions are imposed (Nat0BC) for all boundary segments.
    /// \a fun should only be omitted, if all boundary conditions given are homogenious
    /// and thus no boundary values have to be specified.
    BndDataCL( BndIdxT numbndseg, const BndCondT* bc= 0, const bnd_val_fun* fun= 0);
    /// Deprecated ctor, just for compatibility with older code
    BndDataCL( BndIdxT numbndseg, const bool* isneumann, const bnd_val_fun* fun); // deprecated ctor!

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
           BndSegT GetBndSeg( BndIdxT idx)     const { return BndSegT( BndCond_[idx].GetBC(), BndFun_[idx]); }
    bnd_val_fun    GetBndFun( BndIdxT idx)     const { return BndFun_[idx]; }
    /// \}
};


class NoBndCondCL: public BndCondCL
{
  public:
         NoBndCondCL() : BndCondCL(0) {} // no bnd segments stored
     // default copyctor, dtor, whatever

    template<class SimplexT>
    static inline BndCondT GetBC  (const SimplexT&) { return NoBC; }
    template<class SimplexT>
    static inline bool IsOnDirBnd (const SimplexT&) { return false; }
    template<class SimplexT>
    static inline bool IsOnNatBnd (const SimplexT&) { return false; }
    template<class SimplexT>
    static inline bool IsOnPerBnd (const SimplexT&) { return false; }
};


template<class BndValT= double>
class NoBndDataCL: public NoBndCondCL
{
  public:
    typedef BndSegDataCL<BndValT>             BndSegT;
    typedef BndValT                           bnd_type;
    typedef typename BndSegT::bnd_val_fun     bnd_val_fun;

   // default ctor, dtor, whatever

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

template<class BndValT>
inline BndDataCL<BndValT>::BndDataCL( BndIdxT numbndseg, const BndCondT* bc, const bnd_val_fun* fun)
  : BndCondCL( numbndseg, bc)
{
    BndFun_.resize( numbndseg);
    for (Uint i=0; i<numbndseg; ++i) {
            BndSegDataCL<BndValT>::CheckValid( bc ? bc[i] : Nat0BC, fun ? fun[i] : 0);
        BndFun_[i]= fun ? fun[i] : 0;
    }
}

template<class BndValT>
inline BndDataCL<BndValT>::BndDataCL( BndIdxT numbndseg, const bool* isneumann, const bnd_val_fun* fun) // deprecated ctor!
  : BndCondCL( numbndseg)
{
    BndFun_.resize( numbndseg);
    for (Uint i=0; i<numbndseg; ++i) {
        BndCond_[i]= isneumann[i] ? NatBC : DirBC;
        BndFun_[i]= fun ? fun[i] : 0;
            BndSegDataCL<BndValT>::CheckValid( BndCond_[i].GetBC(), BndFun_[i]);
    }
}


//---------------------------------------
// functions for vertices

inline bool BndCondCL::IsOnDirBnd( const VertexCL& v) const
{ // v is on dir bnd, iff it is on one or more dir bnd segments
    if ( !v.IsOnBoundary() || !BndCond_.size()) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( BndCond_[it->GetBndIdx()].IsDirichlet() )
            return true;
    return false;
}

inline bool BndCondCL::IsOnNatBnd( const VertexCL& v) const
{ // v is on neu bnd, iff it is only on neu bnd segments
    if ( !v.IsOnBoundary() || !BndCond_.size()) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !BndCond_[it->GetBndIdx()].IsNatural() )
            return false;
    return true;
}

inline bool BndCondCL::IsOnPerBnd( const VertexCL& v) const
{ // v is on per bnd, iff it is on one or more per bnd segments and not on a dir bnd segment
    if ( !v.IsOnBoundary() || !BndCond_.size()) return false;
    bool HasPer= false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndCond_[it->GetBndIdx()].GetBC();
        if (bc==DirBC || bc==Dir0BC)
            return false;
        if (bc==Per1BC || bc==Per2BC)
            HasPer= true;
    }
    return HasPer;
}

inline BndCondT BndCondCL::GetBC( const VertexCL& v) const
/// Returns BC on vertex \a v with lowest number (i.e. the superior BC on \a v)
{
    if ( !v.IsOnBoundary() || !BndCond_.size()) return NoBC;
    BndCondT bc= MaxBC_;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        bc= std::min( bc, BndCond_[it->GetBndIdx()].GetBC());
    return bc;
}

inline BndCondT BndCondCL::GetBC( const VertexCL& v, BndIdxT& bidx) const
/// Returns BC on vertex \a v with lowest number (i.e. the superior BC on \a v) and the number of the corresponding boundary segment
{
    if ( !v.IsOnBoundary() || !BndCond_.size()) return NoBC;
    BndCondT bc= MaxBC_, tmp;
    BndIdxT tmp2= NoBndC;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ((tmp= BndCond_[tmp2= it->GetBndIdx()].GetBC()) < bc) {
            bc= tmp;
            bidx= tmp2;
        }
    return bc;
}

template<class SimplexT>
inline const BndCondInfoCL& BndCondCL::GetBndSeg( const SimplexT& s) const
/// Returns boundary segment with superior boundary condition of sub-simplex \a s
{
    BndIdxT  bidx;
    GetBC( s, bidx);
    return BndCond_[bidx];
}


template<class BndValT>
inline typename BndDataCL<BndValT>::BndSegT BndDataCL<BndValT>::GetBndSeg( const VertexCL& v) const
/// Returns bnd segment data on vertex \a v with superior BC
{
    Assert( v.IsOnBoundary() && BndCond_.size(), DROPSErrCL("BndDataCL::GetBndSeg(VertexCL): Not on boundary!"), ~0);
    BndCondT bc_min= MaxBC_;
    BndIdxT  idx= -1;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndCond_[it->GetBndIdx()];
        if (bc<bc_min) { bc_min= bc;    idx= it->GetBndIdx(); }
    }
    return BndSegDataCL<BndValT>( BndCond_[idx], BndFun_[idx]);
}

//---------------------------------------
// functions for edges

inline bool BndCondCL::IsOnDirBnd( const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() || !BndCond_.size()) return false;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( BndCond_[*it].IsDirichlet() )
            return true;
    return false;
}

inline bool BndCondCL::IsOnNatBnd( const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() || !BndCond_.size()) return false;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !BndCond_[*it].IsNatural() )
            return false;
    return true;
}

inline bool BndCondCL::IsOnPerBnd( const EdgeCL& e) const
{ // e is on per bnd, iff it is on one or more per bnd segments and not on a dir bnd segment
    if ( !e.IsOnBoundary() || !BndCond_.size()) return false;
    bool HasPer= false;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndCond_[*it].GetBC();
        if (bc==DirBC || bc==Dir0BC)
            return false;
        if (bc==Per1BC || bc==Per2BC)
            HasPer= true;
    }
    return HasPer;
}

inline BndCondT BndCondCL::GetBC( const EdgeCL& e) const
/// Returns BC on edge \a e with lowest number (i.e. the superior BC on \a e)
{
    if ( !e.IsOnBoundary() || !BndCond_.size()) return NoBC;
    BndCondT bc= MaxBC_;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        bc= std::min( bc, BndCond_[*it].GetBC());
    return bc;
}

inline BndCondT BndCondCL::GetBC( const EdgeCL& e, BndIdxT& bidx) const
/// Returns BC on edge \a e with lowest number (i.e. the superior BC on \a e) and the number of the corresponding boundary segment
{
    if ( !e.IsOnBoundary() || !BndCond_.size()) return NoBC;
    BndCondT bc= MaxBC_, tmp;
    BndIdxT tmp2= NoBndC;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ((tmp= BndCond_[tmp2= *it].GetBC()) < bc) {
            bc= tmp;
            bidx= tmp2;
        }
    return bc;
}

template<class BndValT>
inline typename BndDataCL<BndValT>::BndSegT BndDataCL<BndValT>::GetBndSeg( const EdgeCL& e) const
/// Returns bnd segment data on edge \a e with superior BC
{
    Assert( e.IsOnBoundary() && BndCond_.size(), DROPSErrCL("BndDataCL::GetBndSeg(EdgeCL): Not on boundary!"), ~0);
    BndCondT bc_min= MaxBC_;
    BndIdxT  idx= -1;
    for (const BndIdxT *it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
    {
        const BndCondT bc= BndCond_[*it].GetBC();
        if (bc<bc_min) { bc_min= bc;    idx= *it; }
    }
    return BndSegDataCL<BndValT>( BndCond_[idx], BndFun_[idx]);
}

//---------------------------------------
// functions for faces

inline bool BndCondCL::IsOnDirBnd( const FaceCL& f) const
{
    return f.IsOnBoundary() && BndCond_.size() && BndCond_[f.GetBndIdx()].IsDirichlet();
}

inline bool BndCondCL::IsOnNatBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && BndCond_.size() && BndCond_[f.GetBndIdx()].IsNatural();
}

inline bool BndCondCL::IsOnPerBnd( const FaceCL& f) const
{
    return f.IsOnBoundary() && BndCond_.size() && BndCond_[f.GetBndIdx()].IsPeriodic();
}

inline BndCondT BndCondCL::GetBC( const FaceCL& f) const
/// Returns BC on face \a f
{
    if ( !f.IsOnBoundary() || !BndCond_.size()) return NoBC;
    return BndCond_[f.GetBndIdx()].GetBC();
}

inline BndCondT BndCondCL::GetBC( const FaceCL& f, BndIdxT& bidx) const
/// Returns BC and number of boundary segment on face \a f
{
    if ( !f.IsOnBoundary() || !BndCond_.size()) { bidx= NoBndC; return NoBC; }
    return BndCond_[bidx= f.GetBndIdx()].GetBC();
}

template<class BndValT>
inline typename BndDataCL<BndValT>::BndSegT BndDataCL<BndValT>::GetBndSeg( const FaceCL& f) const
/// Returns bnd segment data on face \a f
{
    Assert( f.IsOnBoundary() && BndCond_.size(), DROPSErrCL("BndDataCL::GetBndSeg(FaceCL): Not on boundary!"), ~0);
    const BndIdxT idx= f.GetBndIdx();
    return BndSegDataCL<BndValT>( BndCond_[idx], BndFun_[idx]);
}


//---------------------------------------------------------
// definitions of BndDataCL<...>::GetXXXBndValue(...)

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetDirBndValue( const VertexCL& v, double t) const
/// Returns value of the Dirichlet boundary value.
/// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end && BndCond_.size(); ++it)
        if ( BndCond_[it->GetBndIdx()].IsDirichlet() )
            return BndFun_[it->GetBndIdx()] ? BndFun_[it->GetBndIdx()]( v.GetCoord(), t) : BndValT();
    throw DROPSErrCL("GetDirBndValue(VertexCL): No Dirichlet Boundary Segment!");
}

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetDirBndValue( const EdgeCL& e, double t) const
/// Returns value of the Dirichlet boundary value.
/// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end && BndCond_.size(); ++it)
        if ( BndCond_[*it].IsDirichlet() )
            return BndFun_[*it] ? BndFun_[*it]( GetBaryCenter(e), t) : BndValT();
    throw DROPSErrCL("GetDirBndValue(EdgeCL): No Dirichlet Boundary Segment!");
}

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetDirBndValue( const FaceCL& f, double t) const
/// Returns value of the Dirichlet boundary value.
/// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    Assert( BndCond_.size() && BndCond_[f.GetBndIdx()].IsDirichlet(), DROPSErrCL("GetDirBndValue(FaceCL): No Dirichlet Boundary Segment!"), ~0);
    return BndFun_[f.GetBndIdx()] ? BndFun_[f.GetBndIdx()]( GetBaryCenter(f), t) : BndValT();
}

template<class BndValT>
inline BndValT BndDataCL<BndValT>::GetNatBndValue( const FaceCL& f, double t) const
/// Returns value of the Neumann boundary value.
/// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    Assert( BndCond_.size() && BndCond_[f.GetBndIdx()].IsNatural(), DROPSErrCL("GetNeuBndValue(FaceCL): No Neumann Boundary Segment!"), ~0);
    return BndFun_[f.GetBndIdx()]( GetBaryCenter(f), t);
}


} //end of namespace DROPS

#endif
