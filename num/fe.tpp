/// \file
/// \brief Description of various finite-element functions

namespace DROPS
{


//**************************************************************************
// Class:   FE_P1CL                                                        *
//**************************************************************************
inline double
FE_P1CL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      default: throw DROPSErrCL("FE_P1CL::H: Invalid shape function.");
    };
}


//**************************************************************************
// Class:   FE_P1DCL                                                       *
//**************************************************************************
inline double
FE_P1DCL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      default: throw DROPSErrCL("FE_P1DCL::H: Invalid shape function.");
    };
}

inline double
FE_P1DCL::H(Uint dof, const BaryCoordCL& p)
{
    switch(dof) {
      case 0: return H0( p);
      case 1: return H1( p);
      case 2: return H2( p);
      case 3: return H3( p);
      default: throw DROPSErrCL("FE_P1DCL::H: Invalid shape function.");
    };
}


//**************************************************************************
// Class:   FE_P2CL                                                        *
//**************************************************************************
inline double
FE_P2CL::H(Uint dof, double v1, double v2, double v3)
{
    switch(dof) {
      case 0: return H0( v1, v2, v3);
      case 1: return H1( v1, v2, v3);
      case 2: return H2( v1, v2, v3);
      case 3: return H3( v1, v2, v3);
      case 4: return H4( v1, v2, v3);
      case 5: return H5( v1, v2, v3);
      case 6: return H6( v1, v2, v3);
      case 7: return H7( v1, v2, v3);
      case 8: return H8( v1, v2, v3);
      case 9: return H9( v1, v2, v3);
      default: throw DROPSErrCL("FE_P2CL::H: Invalid shape function.");
    };
}

inline double
FE_P2CL::H(Uint dof, const BaryCoordCL& p)
{
    switch(dof) {
      case 0: return H0( p);
      case 1: return H1( p);
      case 2: return H2( p);
      case 3: return H3( p);
      case 4: return H4( p);
      case 5: return H5( p);
      case 6: return H6( p);
      case 7: return H7( p);
      case 8: return H8( p);
      case 9: return H9( p);
      default: throw DROPSErrCL("FE_P2CL::H: Invalid shape function.");
    };
}

inline SVectorCL<3>
FE_P2CL::DH0Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret(4.*(v1 + v2 + v3) -3.);
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH1Ref(double v1, double, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v1 -1.; ret[1]= ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH2Ref(double, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= 0.; ret[1]= 4*v2 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH3Ref(double, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 0.; ret[2]= 4.*v3 -1.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH4Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*( 1. - (2.*v1 + v2 + v3) ); ret[1]= ret[2]= -4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH5Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[2]= -4.*v2; ret[1]= 4.*( 1. -(2.*v2 + v1 + v3) );
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH6Ref(double v1, double v2, double)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v2; ret[1]= 4.*v1; ret[2]= 0.;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH7Ref(double v1, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= ret[1]= -4.*v3; ret[2]= 4.*( 1. -(2.*v3 + v1 + v2) );
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH8Ref(double v1, double, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 4.*v3; ret[1]= 0.; ret[2]= 4.*v1;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DH9Ref(double, double v2, double v3)
{
    SVectorCL<3> ret;
    ret[0]= 0.; ret[1]= 4.*v3; ret[2]= 4.*v2;
    return ret;
}
inline SVectorCL<3>
FE_P2CL::DHRef(Uint dof, double v1, double v2, double v3)
{
    switch (dof)
    {
      case 0: return DH0Ref(v1, v2, v3);
      case 1: return DH1Ref(v1, v2, v3);
      case 2: return DH2Ref(v1, v2, v3);
      case 3: return DH3Ref(v1, v2, v3);
      case 4: return DH4Ref(v1, v2, v3);
      case 5: return DH5Ref(v1, v2, v3);
      case 6: return DH6Ref(v1, v2, v3);
      case 7: return DH7Ref(v1, v2, v3);
      case 8: return DH8Ref(v1, v2, v3);
      case 9: return DH9Ref(v1, v2, v3);
      default: throw DROPSErrCL("FE_P2CL::DHRef: Invalid shape function.");
    };
}

inline double
FE_P2CL::Laplace(Uint dof, const SMatrixCL<3,3>& M)
{
    double ret= 0.;
    for (Uint i=0; i<3; ++i)
      for(Uint j=0; j<3; ++j)
        for (Uint k=0; k<3; ++k)
        {
            ret+= M(i,j)*M(i,k)*D2HRef(dof, j, k);
        }
    return ret;
}


//**************************************************************************
// Class:   P1EvalCL                                                       *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1EvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1EvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF( s);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1EvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return  c[0];
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1EvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF( *s.GetVertex( 0));
    c[1]= GetDoF( *s.GetVertex( 1));
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P1CL::H0( v1)
          + GetDoF( *s.GetVertex( 1))*FE_P1CL::H1( v1);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P1EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0]*FE_P1CL::H0( v1) + c[1]*FE_P1CL::H1( v1);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF( *s.GetVertex( i));
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P1CL::H0( v1, v2, v3)
           +GetDoF( *s.GetVertex( 1))*FE_P1CL::H1( v1, v2, v3)
           +GetDoF( *s.GetVertex( 2))*FE_P1CL::H2( v1, v2, v3)
           +GetDoF( *s.GetVertex( 3))*FE_P1CL::H3( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
    inline Data
    P1EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, const BaryCoordCL& p) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P1CL::H0( p)
           +GetDoF( *s.GetVertex( 1))*FE_P1CL::H1( p)
           +GetDoF( *s.GetVertex( 2))*FE_P1CL::H2( p)
           +GetDoF( *s.GetVertex( 3))*FE_P1CL::H3( p);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P1EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return c[0] * FE_P1CL::H0(v1, v2, v3) + c[1] * FE_P1CL::H1(v1, v2, v3)
         + c[2] * FE_P1CL::H2(v1, v2, v3) + c[3] * FE_P1CL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline void
  P1EvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
        DROPSErrCL("P1EvalBaseCL::SetDoF: Trying to assign to"
        "Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()), d);
}


//**************************************************************************
// Class:   P1DEvalCL                                                      *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1DEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumFacesC; ++i)
        c[i]= GetDoF( *s.GetFace( i));
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P1DEvalCL<Data, _BndData, _VD>::GetDoFP1(const TetraCL& s, _Cont& c) const
{
        c[0]= val( s, 0., 0., 0.);
        c[1]= val( s, 1., 0., 0.);
        c[2]= val( s, 0., 1., 0.);
        c[3]= val( s, 0., 0., 1.);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P1DEvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetFace( 0))*FE_P1DCL::H0( v1, v2, v3)
           +GetDoF( *s.GetFace( 1))*FE_P1DCL::H1( v1, v2, v3)
           +GetDoF( *s.GetFace( 2))*FE_P1DCL::H2( v1, v2, v3)
           +GetDoF( *s.GetFace( 3))*FE_P1DCL::H3( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P1DEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return c[0] * FE_P1DCL::H0(v1, v2, v3) + c[1] * FE_P1DCL::H1(v1, v2, v3)
         + c[2] * FE_P1DCL::H2(v1, v2, v3) + c[3] * FE_P1DCL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline void
  P1DEvalCL<Data, _BndData, _VD>::SetDoF(const FaceCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
        DROPSErrCL("P1EvalBaseCL::SetDoF: Trying to assign to"
        "Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()), d);
}


//**************************************************************************
// Class:   P2EvalCL                                                       *
//**************************************************************************
template<class Data, class _BndData, class _VD>
  inline void
  P2EvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
            DROPSErrCL( "P2EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."),
            DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns( _sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return c[0];
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF( s);
}


template<class Data, class _BndData, class _VD>
  inline void
  P2EvalCL<Data, _BndData, _VD>::SetDoF(const EdgeCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
            DROPSErrCL( "P2EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-edge."),
            DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns( _sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF( *s.GetVertex( 0));
    c[1]= GetDoF( *s.GetVertex( 1));
    c[2]= GetDoF( s);
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0] * FE_P2CL::H0( v1) + c[1] * FE_P2CL::H1( v1) + c[2] * FE_P2CL::H2( v1);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF( *s.GetVertex(0))*FE_P2CL::H0( v1)
          + GetDoF( *s.GetVertex(1))*FE_P2CL::H1( v1)
          + GetDoF( s)*FE_P2CL::H2( v1);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const EdgeCL& s) const
{
    return GetDoF( s);
}


template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, Uint face, _Cont& c) const
{
    for(Uint i= 0; i < 3; ++i) {
        c[i]= GetDoF( *s.GetVertex( VertOfFace( face, i)));
        c[i+3]= GetDoF( *s.GetEdge( EdgeOfFace( face, i)));
    }
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2) const
{
    return c[0] * FE_P2CL::H0(v1, v2) + c[1] * FE_P2CL::H1(v1, v2)
         + c[2] * FE_P2CL::H2(v1, v2) + c[3] * FE_P2CL::H3(v1, v2)
         + c[4] * FE_P2CL::H4(v1, v2) + c[5] * FE_P2CL::H5(v1, v2);
}


template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, Uint face, double v1, double v2) const
{
    return  GetDoF( *s.GetVertex( VertOfFace( face, 0)))*FE_P2CL::H0( v1, v2)
           +GetDoF( *s.GetVertex( VertOfFace( face, 1)))*FE_P2CL::H1( v1, v2)
           +GetDoF( *s.GetVertex( VertOfFace( face, 2)))*FE_P2CL::H2( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 0)))*FE_P2CL::H3( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 1)))*FE_P2CL::H4( v1, v2)
           +GetDoF( *s.GetEdge( EdgeOfFace( face, 2)))*FE_P2CL::H5( v1, v2);
}


template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline void
    P2EvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for (Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF( *s.GetVertex( i));
    for (Uint i= 0; i < NumEdgesC; ++i)
        c[i+NumVertsC]= GetDoF( *s.GetEdge( i));
}

template<class Data, class _BndData, class _VD>
  template<class _Cont>
    inline Data
    P2EvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return  c[0] * FE_P2CL::H0( v1, v2, v3) + c[1] * FE_P2CL::H1( v1, v2, v3)
          + c[2] * FE_P2CL::H2( v1, v2, v3) + c[3] * FE_P2CL::H3 (v1, v2, v3)
          + c[4] * FE_P2CL::H4( v1, v2, v3) + c[5] * FE_P2CL::H5( v1, v2, v3)
          + c[6] * FE_P2CL::H6( v1, v2, v3) + c[7] * FE_P2CL::H7( v1, v2, v3)
          + c[8] * FE_P2CL::H8( v1, v2, v3) + c[9] * FE_P2CL::H9( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P2CL::H0( v1, v2, v3)
           +GetDoF( *s.GetVertex( 1))*FE_P2CL::H1( v1, v2, v3)
           +GetDoF( *s.GetVertex( 2))*FE_P2CL::H2( v1, v2, v3)
           +GetDoF( *s.GetVertex( 3))*FE_P2CL::H3( v1, v2, v3)
           +GetDoF( *s.GetEdge( 0))*FE_P2CL::H4( v1, v2, v3)
           +GetDoF( *s.GetEdge( 1))*FE_P2CL::H5( v1, v2, v3)
           +GetDoF( *s.GetEdge( 2))*FE_P2CL::H6( v1, v2, v3)
           +GetDoF( *s.GetEdge( 3))*FE_P2CL::H7( v1, v2, v3)
           +GetDoF( *s.GetEdge( 4))*FE_P2CL::H8( v1, v2, v3)
           +GetDoF( *s.GetEdge( 5))*FE_P2CL::H9( v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
  inline Data
  P2EvalCL<Data, _BndData, _VD>::val(const TetraCL& s, const BaryCoordCL& p) const
{
    return  GetDoF( *s.GetVertex( 0))*FE_P2CL::H0( p)
           +GetDoF( *s.GetVertex( 1))*FE_P2CL::H1( p)
           +GetDoF( *s.GetVertex( 2))*FE_P2CL::H2( p)
           +GetDoF( *s.GetVertex( 3))*FE_P2CL::H3( p)
           +GetDoF( *s.GetEdge( 0))*FE_P2CL::H4( p)
           +GetDoF( *s.GetEdge( 1))*FE_P2CL::H5( p)
           +GetDoF( *s.GetEdge( 2))*FE_P2CL::H6( p)
           +GetDoF( *s.GetEdge( 3))*FE_P2CL::H7( p)
           +GetDoF( *s.GetEdge( 4))*FE_P2CL::H8( p)
           +GetDoF( *s.GetEdge( 5))*FE_P2CL::H9( p);
}

template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::UnknownsMissing(const TetraCL& t) const
{
    const Uint idx= _sol->RowIdx->GetIdx();
    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd();  
         it!=end; ++it) 
        if ( !IsDefinedOn( **it)) return true;
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd();
         it!=end; ++it) 
        if ( !(_bnd->IsOnDirBnd( **it) || (*it)->Unknowns.Exist(idx) )) return true;
    return false;
}


template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(const VertexCL& v) const
{
    return _bnd->IsOnDirBnd( v)
           || (v.Unknowns.Exist() && v.Unknowns.Exist( _sol->RowIdx->GetIdx()));
}

template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(const EdgeCL& e) const
{
    return IsDefinedOn( *e.GetVertex( 0)) && IsDefinedOn( *e.GetVertex( 1))
           && (_bnd->IsOnDirBnd( e)
               || (e.Unknowns.Exist() && e.Unknowns.Exist( _sol->RowIdx->GetIdx())));
}

template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(
    const TetraCL& t, const Uint face) const
{
    const Uint idx= _sol->RowIdx->GetIdx();
    for (Uint i=0; i<3; ++i) {
        const VertexCL& v= *t.GetVertex( VertOfFace( face, i));
        if (!IsDefinedOn( v)) return false;
        const EdgeCL* const ep= t.GetEdge( EdgeOfFace( face, i));
        if (!(_bnd->IsOnDirBnd( *ep)
              || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx))))
            return false;
    }
    return true;
}


template<class Data, class _BndData, class _VD>
inline bool P2EvalCL<Data, _BndData, _VD>::IsDefinedOn(const TetraCL& t) const
{
    for (Uint i=0; i<NumVertsC; ++i)
        if (!IsDefinedOn( *t.GetVertex( i))) return false;
    const Uint idx= _sol->RowIdx->GetIdx();
    for (Uint i=0; i<NumEdgesC; ++i) {
        const EdgeCL* const ep= t.GetEdge( i);
        if (!(_bnd->IsOnDirBnd( *ep)
              || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx))))
            return false;
    }
    return true;
}


//**************************************************************************
// Class:   P1BubbleEvalCL                                                 *
//**************************************************************************
template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s) const
{
    return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s)
        : DoFHelperCL<Data,typename VecDescT::DataType>::get(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()));
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c[0]= GetDoF(s);
}


template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s) const
{
    return DoFHelperCL<Data,typename VecDescT::DataType>::get(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()));
}


template<class Data, class _BndData, class _VD>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert( !_bnd->IsOnDirBnd(s),
        DROPSErrCL("P1EvalCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."),
        DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(
        _sol->Data, s.Unknowns(_sol->RowIdx->GetIdx()), d);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return  c[0];
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF(s);
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c[0]= GetDoF(*s.GetVertex(0));
    c[1]= GetDoF(*s.GetVertex(1));
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    return c[0] * FE_P1BubbleCL::H0(v1) + c[1] * FE_P1BubbleCL::H1(v1);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1BubbleCL::H0(v1)
          + GetDoF(*s.GetVertex(1))*FE_P1BubbleCL::H1(v1);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    for(Uint i= 0; i < NumVertsC; ++i)
        c[i]= GetDoF(*s.GetVertex(0));
    c[4]= GetDoF(s);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    return c[0] * FE_P1BubbleCL::H0(v1, v2, v3) + c[1] * FE_P1BubbleCL::H1(v1, v2, v3)
         + c[2] * FE_P1BubbleCL::H2(v1, v2, v3) + c[3] * FE_P1BubbleCL::H3(v1, v2, v3)
         + c[4] * FE_P1BubbleCL::H4(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1BubbleCL::H0(v1, v2, v3)
           +GetDoF(*s.GetVertex(1))*FE_P1BubbleCL::H1(v1, v2, v3)
           +GetDoF(*s.GetVertex(2))*FE_P1BubbleCL::H2(v1, v2, v3)
           +GetDoF(*s.GetVertex(3))*FE_P1BubbleCL::H3(v1, v2, v3)
           +GetDoF(s)*FE_P1BubbleCL::H4(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::lin_val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1BubbleCL::H0(v1, v2, v3)
           +GetDoF(*s.GetVertex(1))*FE_P1BubbleCL::H1(v1, v2, v3)
           +GetDoF(*s.GetVertex(2))*FE_P1BubbleCL::H2(v1, v2, v3)
           +GetDoF(*s.GetVertex(3))*FE_P1BubbleCL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::bubble_val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(s)*FE_P1BubbleCL::H4(v1, v2, v3);
}


template<class Data, class _BndData, class _VD>
void
Interpolate(P1EvalCL<Data, _BndData, _VD>& sol, const P1EvalCL<Data, _BndData, const _VD>& old_sol)
// This only works, if Interpolate is called after every refinement of the multigrid.
// Take care, that x and old_x are on successive triangulations.
{
    typedef typename P1EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();

    const MultiGridCL& _MG= old_sol.GetMG();    
    const Uint old_level= old_sol.GetLevel();
    const Uint level= sol.GetLevel();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->GetIdx();
    Uint counter1= 0, counter2= 0;
    
//    Assert( level==old_level || old_level==level-1, DROPSErrCL("Interpolate: successive triangs are expected\n"));
    // Iterate over all edges, interpolate values on new mid vertices
    for (MultiGridCL::const_EdgeIterator sit= _MG.GetAllEdgeBegin(level), theend= _MG.GetAllEdgeEnd(level);
         sit!=theend; ++sit)
        if ( sit->IsRefined() && !_bnd->IsOnDirBnd(*sit->GetMidVertex())  ) // only new non-boundary vertices are interpolated
        {
            sol.SetDoF(*sit->GetMidVertex(), (old_sol.val(*sit->GetVertex(0)) + old_sol.val(*sit->GetVertex(1)))/2.);
            ++counter2;
        }
    // Iterate over the vertices of the old triangulation and copy values
    for (MultiGridCL::const_TriangVertexIteratorCL sit= _MG.GetTriangVertexBegin(old_level), theend= _MG.GetTriangVertexEnd(old_level);
         sit!=theend; ++sit)
        if ( !_bnd->IsOnDirBnd(*sit) )
        {
// @@@ Handle missing unknowns: Do it right here?
if (sit->Unknowns.Exist(old_idx))        
            sol.SetDoF(*sit, old_sol.val(*sit));
//else they were set in the for-loop before !
            ++counter1;
        }

        std::cerr << "Interpolate: " << counter1 << " vertex-Dof of " << old_sol.GetSolution()->Data.size() << " copied, " 
                                     << counter2 << " new Mid-vertex-DoF interpolated." << std::endl;
}


//**************************************************************************
// RepairAfterRefine: Repairs the P1-function old_f, which is possibly     *
//     damaged by a grid-refinement. This is done by copying to the new    *
//     VecDescBaseCL object vecdesc and interpolation of old values.       *
// Precondition: old_f is a damaged P1-function (through at most one       *
//     refinement step), vecdesc contains a valid IdxDescCL* to an index on*
//     the same level as old_f (If this level was deleted, vecdesc shall be*
//     defined on the next coarser level.); vecdesc.Data has the correct   *
//     size.                                                               *
// Postcondition: vecdesc, together with the boundary-data of old_f,       *
//     represents a P1-function on the triangulation tl. If old_f was      *
//     defined on the last level before refinement, which is then deleted, *
//     tl ==  old_f.GetLevel() -1; else tl is the level of old_f.          *
//**************************************************************************    
template <class P1T, class VecDesc>
  Uint
  RepairAfterRefineP1( const P1T& old_f, VecDesc& vecdesc)
{
    Uint tl= old_f.GetLevel();
    const MultiGridCL& MG= old_f.GetMG();
    const Uint maxlevel= MG.GetLastLevel();

    // If the level, on which f is defined, was completely deleted, old_f should
    // be lifted to the next coarser level. It is left to the caller to do that
    // -- use the return value of the current function.
    if (tl > maxlevel) {
        Assert( tl == maxlevel+1,
                "RepairAfterRefine (P1): old_f is defined on a level, "
                "which cannot have existed in the previous multigrid.",
                DebugNumericC);
        tl= maxlevel;
    }
    Assert( tl == vecdesc.GetLevel(),
            "RepairAfterRefine (P1): old and new function are "
            "defined on incompatible levels.",
            DebugNumericC);
    const Uint old_idx= old_f.GetSolution()->RowIdx->GetIdx();
    const Uint idx= vecdesc.RowIdx->GetIdx();
    typedef typename P1T::BndDataCL BndCL;
    BndCL* const bnd= old_f.GetBndData();
    typename P1T::modifiable_type f( &vecdesc, old_f.GetBndData(), &MG);
    f.SetTime( old_f.GetTime());
    Uint counter1= 0, counter2= 0;

    // Iterate over all edges on grids with smaller level than tl. If the
    // edge **is refined and its midvertex *does not have the old index,
    // *, has the new index, *is not on a Dirichlet-boundary, then
    // the new value on the midvertex can be interpolated from the vertices
    // of the edge. (These have an old value as the edge has existed in the old grid
    // and because of consistency so have its vertices.)
    // As new vertices in a grid are always midvertices of edges on the next coarser grid,
    // all new vertices in triangulation tl are initialized.
    // If f is defined on tl==0, there is no coarser level, thus this step is to be skipped.
    if (tl > 0) {
        for (MultiGridCL::const_EdgeIterator sit= MG.GetAllEdgeBegin(tl-1),
             theend= MG.GetAllEdgeEnd(tl-1); sit!=theend; ++sit)
            if ( sit->IsRefined()
                 && sit->GetMidVertex()->Unknowns.Exist()
                 && !sit->GetMidVertex()->Unknowns.Exist( old_idx)
                 && sit->GetMidVertex()->Unknowns.Exist( idx)
                 && !bnd->IsOnDirBnd( *sit->GetMidVertex())) {
                f.SetDoF( *sit->GetMidVertex(),
                          (old_f.val( *sit->GetVertex(0)) + old_f.val( *sit->GetVertex(1)))*0.5);
                ++counter2;
            }
    }
    // All vertices in tl, that **have the old index and **are not on the Dirichlet-boundary
    // hold an old value, which is copied to the new vector.
    // As the check for the Dirichlet-boundary is the first, sit->Unknowns.Exist(),
    // can be spared, since idx is required to exist on tl.
    for (MultiGridCL::const_TriangVertexIteratorCL sit= MG.GetTriangVertexBegin( tl),
         theend= MG.GetTriangVertexEnd( tl); sit!=theend; ++sit)
        if (!bnd->IsOnDirBnd( *sit) && sit->Unknowns.Exist( old_idx)) {
            f.SetDoF( *sit, old_f.val( *sit));
            ++counter1;
        }
    Comment( "RepairAfterRefine (P1): " << counter1 << " vertex-dof of "
             << old_f.GetSolution()->Data.size() << " copied, " << counter2
             << " new mid-vertex-doF interpolated." << std::endl,
             DebugNumericC);
    return tl;
}

template<class Data, class _BndData, class _VD>
void InterpolateChildren( const TetraCL& t, P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
{
    typedef typename P2EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();

    const double edgebary[3][2]=
        { {0.25, 0.25},
          {0.5 , 0.25},
          {0.25, 0.5}
        };
      
    // Hole des Tetraeders RefRule; gehe ueber alle Kinder ;-): durchlaufe die edges in
    // der childrule des kindes: falls IsSubEdge(edge): finde mit ParentEdge & NumOfSubEdge heraus,
    // von welcher Kante des tetras, (falls der wert bei edge != 0) und interpoliere ueber kante. 
    // sonst, falls IsSubInParFace(subedge): WhichEdgeInFace(subedge, face, pos) und interpoliere ueber face;
    // sonst, behandele raumdiagonale;
    const RefRuleCL& refrule= t.GetRefData();
    TetraCL::const_ChildPIterator child= t.GetChildBegin();
    const TetraCL::const_ChildPIterator childend= t.GetChildEnd();
    for (Uint childnum=0; child!=childend; ++childnum, ++child)
    {
        const ChildDataCL& childdata= GetChildData(refrule.Children[childnum]);
        for (Uint chedge=0; chedge<NumEdgesC; ++chedge)
        {
            const EdgeCL* const edgep= (*child)->GetEdge(chedge);
            const Uint chedgeinparent= childdata.Edges[chedge];
            if (!_bnd->IsOnDirBnd(*edgep))
            {
                if ( IsSubEdge(chedgeinparent) )
                {
                    const Uint paredge= ParentEdge(chedgeinparent);
                    const Uint num= NumOfSubEdge(chedgeinparent);
                    sol.SetDoF( *edgep, old_sol.val(*t.GetEdge(paredge), 0.25+0.5*num) );
                }
                else if ( IsSubInParFace(chedgeinparent) )
                {
                    Uint parface;
                    Uint pos;
                    WhichEdgeInFace(chedgeinparent, parface, pos);
                    sol.SetDoF( *edgep, old_sol.val(t, parface, edgebary[pos][0], edgebary[pos][1]) );
                }
                else
                {
                    sol.SetDoF( *edgep, old_sol.val(t, 0.25, 0.25, 0.25) );
                }
            }
        }
    }
}


template<class Data, class _BndData, class _VD>
inline void CopyValues( const TetraCL& t, P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
{
    typedef typename P2EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();
    
    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd();
        it!=end; ++it)
        if (!_bnd->IsOnDirBnd( **it) )
            sol.SetDoF( **it, old_sol.val( **it) );
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd();
        it!=end; ++it)
        if (!_bnd->IsOnDirBnd( **it) )
            sol.SetDoF( **it, old_sol.val( **it, 0.5) );
}

template<class Data, class _BndData, class _VD>
void Adapt( P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
{
// Adapt a solution on a triangulation of a certain level, that has changed during the refinement.
// Notation: T = old triang, T' = new triang. Both T and T' are of the same level.
// This change can be classified in several cases (not complete...): Fot tetra t in T:
//    a) t is missing not only in T' but in all new triang levels  
//       -> information is lost
//    b) t is missing in T', but is member of the new MultiGrid
//       -> t was refined, children of t are members of T', for these information is interpolated from t!
//    c) t and its brotherhood were replaced by other children
//       -> change of the parents' refinement rule, restore information in parent and interpolate!
// TODO: missing: handling of unrefined Tetras

// Adapt should be very robust in all occuring situations!!!

    const MultiGridCL& _MG= old_sol.GetMG();    
    const Uint level= sol.GetLevel();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->GetIdx();
    const Uint NumUnknowns=  old_sol.GetSolution()->RowIdx->NumUnknownsEdge;
    
    Assert( level==old_sol.GetLevel(),
        DROPSErrCL("Adapt: Same triang levels are expected\n"), ~0);

    // 1. Iterate tetras of triangulation: Interpolate missing unknowns
    // 2. Iterate tetras of triangulation: Copy known values
    // => known values override interpolated values

    for (MultiGridCL::const_TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(level), theend= _MG.GetTriangTetraEnd(level);
         sit!=theend; ++sit)
    {
         if ( old_sol.UnknownsMissing( *sit) )
         // Tetra war in T vorhanden -> neu in T'! Fall b)c)
             if ( sol.UnknownsMissing( *sit) ) 
             {    
                 // Ergaenze evtl. fehlenden Idx des Vaters auf Edges! (nur fuer c) noetig)
                 for (TetraCL::const_EdgePIterator it= sit->GetParent()->GetEdgesBegin(), end= sit->GetParent()->GetEdgesEnd();
                     it!=end; ++it)
                     if ((*it)->IsRefined() && (*it)->GetMidVertex()->Unknowns.Exist(old_idx) )
                     {
                         // evtl. UnknownIdx fuer old_idx anlegen
                         (*it)->Unknowns.Prepare( old_idx);
                         // Indexwerte von MidVertex kopieren
                         //for (Uint i=0; i<NumUnknowns; ++i)
                             (*it)->Unknowns(old_idx)= (*it)->GetMidVertex()->Unknowns(old_idx);
                     }
                 // Interpoliere Werte vom Vater
                 InterpolateChildren( *sit->GetParent(), sol, old_sol);
             }
    }
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(level), theend= _MG.GetTriangTetraEnd(level);
         sit!=theend; ++sit)
    {
        if ( !old_sol.UnknownsMissing(*sit) )
            CopyValues( *sit, sol, old_sol);
    }
}




template<class Data, class _BndData, class _VD>
void Interpolate(P2EvalCL<Data, _BndData, _VD>& sol, const P2EvalCL<Data, _BndData, const _VD>& old_sol)
// This only works, if Interpolate is called after every refinement of the multigrid.
// Take care, that x and old_x are on successive triangulations.
{
    typedef typename P2EvalCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->GetIdx();
      
    // All velocity-components use the same row-index and the same trianglevel
    const MultiGridCL& _MG= old_sol.GetMG();    
    const Uint old_level= old_sol.GetLevel();
    const Uint level= sol.GetLevel();
    Uint num_vert_copy= 0, num_edge_copy= 0, num_child_edge= 0;
    
//    Assert( level==old_level || old_level==level-1, DROPSErrCL("Interpolate: successive triangs are expected\n"));

    // Iterate over all tetras in old_level, interpolate edge-DoF on diagonal, if there is one;
    // interpolate edge-DoF's of the children's non-copied edges
    for (MultiGridCL::const_TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(old_level), theend= _MG.GetTriangTetraEnd(old_level);
         sit!=theend; ++sit)
    {
        // If *sit is unrefined, all interpolation will be done via copying DoF on edges and vertices later
        if ( !sit->IsUnrefined() && (*sit->GetChildBegin())->IsInTriang(level) )
            InterpolateChildren( *sit, sol, old_sol);
    }
    // Iterate over all edges, interpolate values on new mid-vertices and edge-copies (plain copying)
    for (MultiGridCL::const_EdgeIterator sit= _MG.GetAllEdgeBegin(level), theend= _MG.GetAllEdgeEnd(level);
         sit!=theend; ++sit)
    {
        if ( sit->IsRefined() && sit->IsInTriang(old_level)
             && sit->GetMidVertex()->IsInTriang(level) && !_bnd->IsOnDirBnd(*sit->GetMidVertex()) ) // only new non-boundary vertices are interpolated
        {
            if (!sit->Unknowns.Exist(old_idx)) continue;        
                sol.SetDoF( *sit->GetMidVertex(), old_sol.val(*sit, 0.5) );
                ++num_edge_copy;
        }
        else if ( sit->IsInTriang(old_level) && sit->IsInTriang(level) && !_bnd->IsOnDirBnd(*sit) )
            {
                if (!sit->Unknowns.Exist(old_idx)) continue;        
                    sol.SetDoF( *sit, old_sol.val(*sit, 0.5) );
                    ++num_edge_copy;
            }
    }
    // Iterate over the vertices of the old triangulation and copy values
    for (MultiGridCL::const_TriangVertexIteratorCL sit= _MG.GetTriangVertexBegin(old_level), theend= _MG.GetTriangVertexEnd(old_level);
         sit!=theend; ++sit)
        if ( !_bnd->IsOnDirBnd(*sit) )
        {
            if (!sit->Unknowns.Exist(old_idx)) continue;        
                sol.SetDoF(*sit, old_sol.val(*sit) );
                ++num_vert_copy;
        }

        std::cerr << "Interpolate: " << num_vert_copy << " vertex-DoF copied, " 
                                     << num_edge_copy << " edge-DoF copied, "
                                     << num_child_edge << " edge-DoF interpolated."
                                     << std::endl;
}

// Helper function for RepairOnChildren.
template <class P2T>
  bool
  TryGetValuesOnFace( const P2T& f,
      std::vector<typename P2T::DataT>& v, const TetraCL& t, Uint face)
{
    typedef typename P2T::BndDataCL BndCL;
    const BndCL* const bnd= f.GetBndData();
    const Uint idx= f.GetSolution()->RowIdx->GetIdx() ;

    for (Uint i=0; i<3; ++i) {
        const VertexCL* const vp= t.GetVertex( VertOfFace( face, i));
        if (bnd->IsOnDirBnd( *vp)
            || (vp->Unknowns.Exist() && vp->Unknowns.Exist( idx)))
            v[i]= f.val( *vp);
        else return false;
    }
    for (Uint i=0; i<3; ++i) {
        const EdgeCL* const ep= t.GetEdge( EdgeOfFace( face, i));
        if (bnd->IsOnDirBnd( *ep)
            || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx)))
            v[i+3]= f.val( *ep);
        else { // Never on dirichlet-boundary: this is handled by the preceeding if.
            if (ep->IsRefined()
                && ep->GetMidVertex()->Unknowns.Exist()
                && ep->GetMidVertex()->Unknowns.Exist( idx))
                v[i+3]= f.val( *ep->GetMidVertex());
            else return false;
        }
    }    
    return true;
}

// Helper function for RepairOnChildren.
template <class P2T>
  bool
  TryGetValuesOnTetra( const P2T& f,
      std::vector<typename P2T::DataT>& v, const TetraCL& t)
{
    typedef typename P2T::BndDataCL BndCL;
    const BndCL* const bnd= f.GetBndData();
    const Uint idx= f.GetSolution()->RowIdx->GetIdx() ;

    for (Uint i=0; i<4; ++i) {
        const VertexCL* const vp= t.GetVertex( i);
        if (bnd->IsOnDirBnd( *vp)
            || (vp->Unknowns.Exist() && vp->Unknowns.Exist( idx)))
            v[i]= f.val( *vp);
        else return false;
    }
    for (Uint i=0; i<6; ++i) {
        const EdgeCL* const ep= t.GetEdge( i);
        if (bnd->IsOnDirBnd( *ep)
            || (ep->Unknowns.Exist() && ep->Unknowns.Exist( idx)))
            v[i+4]= f.val( *ep);
        else { // Never on dirichlet-boundary: this is handled by the preceeding if.
            if (ep->IsRefined()
                && ep->GetMidVertex()->Unknowns.Exist()
                && ep->GetMidVertex()->Unknowns.Exist( idx))
                v[i+4]= f.val( *ep->GetMidVertex());
            else return false;
        }
    }    
    return true;
}


// Used in RepairAfterRefine for P2-elements. The all possible parent tetras of
// the new function f are walked over, and scanned for available indices of the
// old function old_f. These are used to define f by quadratic interpolation. If
// a former (ir)regularly refined tetra is now refined differently, it may be
// neccessary to interpolate only linearly because of lack of sufficient data.
template <class P2T, class VecDesc>
void
RepairOnChildren( const TetraCL& t,  const P2T& old_f, VecDesc& vecdesc)
{
    typedef typename P2T::BndDataCL BndCL;
    const BndCL* const bnd= old_f.GetBndData();
    const Uint old_idx= old_f.GetSolution()->RowIdx->GetIdx();
    const Uint idx= vecdesc.RowIdx->GetIdx() ;
    typename P2T::modifiable_type f( &vecdesc, old_f.GetBndData(), &old_f.GetMG(), old_f.GetTime());
    f.SetTime( old_f.GetTime());

    const double edgebary[3][2]= {
          {0.25, 0.25},
          {0.5 , 0.25},
          {0.25, 0.5}
        };
      
    const RefRuleCL& refrule= t.GetRefData();
    TetraCL::const_ChildPIterator child= t.GetChildBegin();
    const TetraCL::const_ChildPIterator childend= t.GetChildEnd();
    for (Uint childnum=0; child!=childend; ++childnum, ++child) {
        const ChildDataCL& childdata= GetChildData( refrule.Children[childnum]);
        for (Uint chedge=0; chedge<NumEdgesC; ++chedge) {
            const EdgeCL* const edgep= (*child)->GetEdge( chedge);
            if (!bnd->IsOnDirBnd( *edgep) 
                && edgep->Unknowns.Exist()
                && edgep->Unknowns.Exist( idx)
                && !edgep->Unknowns.Exist( old_idx)) { // not just a copy on the edge itself
                const Uint chedgeinparent= childdata.Edges[chedge];
                if (IsSubEdge( chedgeinparent)) { // Sub-Edges of parent edges
                    const Uint paredge= ParentEdge( chedgeinparent);
                    const EdgeCL* const paredgep= t.GetEdge( paredge);
                    if (paredgep->Unknowns.Exist()
                        && paredgep->Unknowns.Exist( old_idx)) { // refinement took place,
                                                                 // interpolate edge-dof
                        const Uint num= NumOfSubEdge( chedgeinparent);
                        f.SetDoF( *edgep, old_f.val( *paredgep, 0.25+0.5*num) );
                    }
                    else {  // this edge has been unrefined;
                            // There is not enough information available for 
                            // more than linear interpolation.
                        f.SetDoF( *edgep, 0.5*(old_f.val( *edgep->GetVertex( 0))
                                               + old_f.val( *edgep->GetVertex( 1))));
                    }
                }
                else if (IsSubInParFace( chedgeinparent)) { // Sub-edges in parent faces
                    Uint parface, pos;
                    WhichEdgeInFace(chedgeinparent, parface, pos);
                    std::vector<typename P2T::DataT> v( 6);
                    if (TryGetValuesOnFace( old_f, v, t, parface)) {
                        f.SetDoF( *edgep, old_f.val( v, edgebary[pos][0], edgebary[pos][1]));
                    }
                    else { // These are complicated; the parent-face was refined and
                           // is now refined (differently). Linear interpolation
                           // between the vertices is used for now.
                           // The values of f in the vertices of face are used
                           // (not old_f), because the vertices itself might be new.
                           // These have already been set in RepairAfterRefine.
                        f.SetDoF( *edgep, 0.5*(f.val( *edgep->GetVertex( 0))
                                               + f.val( *edgep->GetVertex( 1))));
                    }
                }
                // space diagonal
                else {
                    std::vector<typename P2T::DataT> v( 10);
                    if (TryGetValuesOnTetra( old_f, v, t)) {
                        f.SetDoF( *edgep, old_f.val( v, 0.25, 0.25, 0.25));
                    }
                    else {
                        f.SetDoF( *edgep, 0.5*(f.val( *edgep->GetVertex( 0))
                                               + f.val( *edgep->GetVertex( 1))));
                    }
                }
            }
        }
    }
}

//**************************************************************************
// RepairAfterRefine: Repairs the P2-function old_f, which is possibly     *
//     damaged by a grid-refinement. This is done by copying to the new    *
//     VecDescBaseCL object vecdesc and interpolation of old values.       *
// Precondition: old_f is a damaged P2-function (through at most one       *
//     refinement step), vecdesc contains a valid IdxDescCL* to an index on*
//     the same level as old_f (If this level was deleted, vecdesc shall be*
//     defined on the next coarser level.); vecdesc.Data has the correct   *
//     size.                                                               *
// Postcondition: vecdesc, together with the boundary-data of old_f,       *
//     represents a P2-function on the triangulation tl. If old_f was      *
//     defined on the last level before refinement, which is then deleted, *
//     tl ==  old_f.GetLevel -1; else tl is the level of old_f.            *
//**************************************************************************    
template <class P2T, class VecDesc>
Uint
RepairAfterRefineP2( const P2T& old_f, VecDesc& vecdesc)
{
    const Uint tl= vecdesc.GetLevel();
    const MultiGridCL& MG= old_f.GetMG();
    const Uint old_idx= old_f.GetSolution()->RowIdx->GetIdx();
    const Uint idx= vecdesc.RowIdx->GetIdx();
    typedef typename P2T::BndDataCL BndCL;
    BndCL* const bnd= old_f.GetBndData();
    typename P2T::modifiable_type f( &vecdesc, old_f.GetBndData(), &MG, old_f.GetTime());

    // The first two loops interpolate the values of all vertices (new ones as
    // mid-vertices and old ones as copies). This works similar to the P1 case.    

    // Iterate over all edges on grids with smaller level than tl. If the
    // edge **is refined and its midvertex *does not have the old index,
    // *, has the new index, *is not on a Dirichlet-boundary, then
    // the new value on the midvertex can be interpolated from the vertices
    // of the edge and the edge-dof. (These have an old value as the edge has
    // existed in the old grid and because of consistency so have its vertices.)
    // (Actually, the old vertices do not matter, as their dof are 0 in the midvertex.)
    // As new vertices in a grid are always midvertices of edges on the next coarser grid,
    // all new vertices in triangulation tl are initialized.
    // Also iterate over tl. This will not trigger the first if-statrement below as 
    // sit->GetMidVertex()->Unknowns.Exist( idx) will always be false; this is used to
    // copy edge-dofs, that did not change.
    Uint counter2= 0, counter3= 0;
    for (MultiGridCL::const_EdgeIterator sit= MG.GetAllEdgeBegin( tl),
         theend= MG.GetAllEdgeEnd( tl); sit!=theend; ++sit) {
        if (sit->IsRefined()
            && sit->GetMidVertex()->Unknowns.Exist()
            && !sit->GetMidVertex()->Unknowns.Exist( old_idx)
            && sit->GetMidVertex()->Unknowns.Exist( idx)) {
            f.SetDoF( *sit->GetMidVertex(), old_f.val( *sit));
        }
        else if (sit->Unknowns.Exist()
                 && sit->Unknowns.Exist( old_idx)
                 && sit->Unknowns.Exist( idx)) {
                f.SetDoF( *sit, old_f.val( *sit));   
                ++counter3;
            }
    }
    // All vertices in tl, that **have the old index and **are not on the Dirichlet-boundary
    // hold an old value, which is copied to the new vector.
    // As the check for the Dirichlet-boundary is the first, sit->Unknowns.Exist(),
    // can be spared, since idx is required to exist on tl.
    for (MultiGridCL::const_TriangVertexIteratorCL sit= MG.GetTriangVertexBegin( tl),
         theend= MG.GetTriangVertexEnd( tl); sit!=theend; ++sit)
        if (!bnd->IsOnDirBnd( *sit) && sit->Unknowns.Exist( old_idx)) {
            f.SetDoF( *sit, old_f.val( *sit));
            ++counter2;
        }
    // Coarsened edges on level 0 cannot be handled on child-tetras,
    // as level-0-tetras are never children. Therefore they are handled
    // here instead of in RepairOnChildren (see below). This is possible
    // as edges on level 0 never have parents, which would otherwise have to
    // be checked for refinement.
    // Linear interpolation.
    for (MultiGridCL::const_EdgeIterator sit= MG.GetAllEdgeBegin( 0),
         theend= MG.GetAllEdgeEnd( 0); sit!=theend; ++sit) {
        if (sit->Unknowns.Exist()
            && sit->Unknowns.Exist( idx)
            && !sit->Unknowns.Exist( old_idx)
            && !bnd->IsOnDirBnd( *sit)) {
            f.SetDoF( *sit, 0.5*(old_f.val( *sit->GetVertex( 0))
                                 + old_f.val( *sit->GetVertex( 1))));
        }
    }
    // Handle edges e that were refined: The condition that identifies such
    // an edge is given above -- unfortunately, there is no information
    // concerning the child edges available on the edge. Thus an iteration
    // over tetras is used.
    // The remaining cases, where green rules change, are handled here, too.
    if (tl>0)
        for (MultiGridCL::const_TetraIterator sit= MG.GetAllTetraBegin( tl-1),
             theend= MG.GetAllTetraEnd( tl-1); sit!=theend; ++sit) {
            // If *sit is unrefined, there are no diagonals and childedge-dof to interpolate.
            if ( !sit->IsUnrefined() )
                RepairOnChildren( *sit, old_f, vecdesc);
        }
    Comment( "RepairAfterRefine (P2): " << " Vertex-dof copies: " << counter2
             << " Edge-dof copies: " << counter3
             << " I'm too dumb to count the rest accurately",
             DebugNumericC);
    return tl;
}


template <class VecDescT, class BndDataT, class Cont>
void RestrictP2(const TetraCL& s, const VecDescT& vd, const BndDataT& bnd, Cont& c, double t)
{
    const Uint slvl= s.GetLevel();
    const Uint flvl= vd.GetLevel();
    Assert( slvl<=flvl, DROPSErrCL("RestrictP2: Tetra is on a finer level"
            "than the function."), ~0);

    typedef typename VecDescT::DataType VecT;
    typedef DoFHelperCL< typename BndDataT::bnd_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();
    for (Uint i= 0; i<NumVertsC; ++i)
        c[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), t);
    for (Uint i= 0; i<NumEdgesC; ++i) {
        const EdgeCL& e= *s.GetEdge( i);
        c[i+NumVertsC]= (slvl < flvl && e.IsRefined())
            ? ( !bnd.IsOnDirBnd( *e.GetMidVertex())
                ? DoFT::get( v, e.GetMidVertex()->Unknowns( idx))
                : bnd.GetDirBndValue( *e.GetMidVertex(), t))
            : ( !bnd.IsOnDirBnd( e)
                ? DoFT::get( v, e.Unknowns( idx))
                : bnd.GetDirBndValue( e, t));
    }
}


} // end of namespace DROPS
