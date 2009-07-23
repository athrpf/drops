//**************************************************************************
// File:    discretize.tpp                                                 *
// Content: discretizations for several PDEs and FE types                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Feb, 1 2005                                            *
//**************************************************************************

namespace DROPS
{

//**************************************************************************
// Class:   LocalP1CL                                                      *
//**************************************************************************
template<class T>
  inline LocalP1CL<T>&
  LocalP1CL<T>::assign(const TetraCL& s, instat_fun_ptr f, double t)
{
    for (Uint i= 0; i< NumVertsC; ++i)
        (*this)[i]= f( s.GetVertex( i)->GetCoord(), t);
    return *this;
}

template<class T>
  template<class BndDataT>
    inline LocalP1CL<T>&
    LocalP1CL<T>::assign(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd, double t)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint idx= vd.RowIdx->GetIdx();
    for (Uint i= 0; i< NumVertsC; ++i)
        (*this)[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
            ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
            : bnd.GetDirBndValue( *s.GetVertex( i), t);
    return *this;
}

template<class T>
  template<class P1FunT>
    inline LocalP1CL<T>&
    LocalP1CL<T>::assign(const TetraCL& s, const P1FunT& f, double t)
{
    const double tmp= f.GetTime();
    f.SetTime( t);
    f.GetDoF( s, *this);
    f.SetTime( tmp);
    return *this;
}


template<class T>
  LocalP1CL<T>::LocalP1CL(const TetraCL& s, instat_fun_ptr f , double t)
  : base_type( value_type(), FE_P1CL::NumDoFC)
{
    this->assign( s, f, t);
}

template<class T>
  template <class P1FunT>
    LocalP1CL<T>::LocalP1CL(const TetraCL& s, const P1FunT& f, double t)
    : base_type( value_type(), FE_P1CL::NumDoFC)
{
    this->assign( s, f, t);
}

template<class T>
  template<class BndDataT>
    LocalP1CL<T>::LocalP1CL(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd, double t)
    : base_type( value_type(), FE_P1CL::NumDoFC)
{
    this->assign( s, vd, bnd, t);
}

template<class T>
  inline typename LocalP1CL<T>::value_type
  LocalP1CL<T>::operator() (const BaryCoordCL& p) const
{
    return FE_P1CL::val( *this, p);
}


//**************************************************************************
// Class:   LocalP2CL                                                      *
//**************************************************************************
template<class T>
  inline LocalP2CL<T>&
  LocalP2CL<T>::assign(const TetraCL& s, instat_fun_ptr f, double t)
{
    for (Uint i= 0; i< NumVertsC; ++i)
        (*this)[i]= f( s.GetVertex( i)->GetCoord(), t);
    for (Uint i= 0; i< NumEdgesC; ++i)
        (*this)[i+NumVertsC]= f( GetBaryCenter( *s.GetEdge( i)), t);
    return *this;
}

template<class T>
  template<class BndDataT>
    inline LocalP2CL<T>&
    LocalP2CL<T>::assign(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd, double t)
{
    typedef VecDescCL::DataType VecT;
    typedef DoFHelperCL<value_type, VecT> DoFT;
    const VecT& v= vd.Data;
    const Uint tlvl= s.GetLevel();
    const Uint vlvl= vd.GetLevel();
    const Uint idx= vd.RowIdx->GetIdx();
    if (tlvl == vlvl) {
        for (Uint i= 0; i< NumVertsC; ++i)
            (*this)[i]= !bnd.IsOnDirBnd( *s.GetVertex( i))
                ? DoFT::get( v, s.GetVertex( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetVertex( i), t);
        for (Uint i= 0; i< NumEdgesC; ++i)
            (*this)[i+NumVertsC]= !bnd.IsOnDirBnd( *s.GetEdge( i))
                ? DoFT::get( v, s.GetEdge( i)->Unknowns( idx))
                : bnd.GetDirBndValue( *s.GetEdge( i), t);
    }
    else {
        if (tlvl < vlvl) RestrictP2( s, vd, bnd, *this);
        else throw DROPSErrCL( "LocalP2CL::Assign: Prolongation not implemented.\n");
    }
    return *this;
}

template<class T>
  template<class P2FunT>
    inline LocalP2CL<T>&
    LocalP2CL<T>::assign(const TetraCL& s, const P2FunT& f, double t)
{
    const Uint tlvl= s.GetLevel();
    const Uint flvl= f.GetLevel();
    const double tmp= f.GetTime();
    f.SetTime( t);
    if (tlvl == flvl)
        f.GetDoF( s, *this);
    else
        if (tlvl < flvl) RestrictP2( s, f, *this);
        else throw DROPSErrCL( "LocalP2CL::Assign: Prolongation not implemented.\n");
    f.SetTime( tmp);
    return *this;
}

template<class T>
  inline LocalP2CL<T>&
  LocalP2CL<T>::assign(const LocalP1CL<T>& p1, double)
{
    for (size_t i= 0; i < 4; ++i)
        (*this)[i]= p1[i];
    for (size_t i= 0; i < 6; ++i)
        (*this)[i + 4]= 0.5*(p1[VertOfEdge( i, 0)] + p1[VertOfEdge( i, 1)]);
    return *this;
}


template<class T>
  LocalP2CL<T>::LocalP2CL(const TetraCL& s, instat_fun_ptr f , double t)
  : base_type( value_type(), FE_P2CL::NumDoFC)
{
    this->assign( s, f, t);
}

template<class T>
  template <class P2FunT>
    LocalP2CL<T>::LocalP2CL(const TetraCL& s, const P2FunT& f, double t)
    : base_type( value_type(), FE_P2CL::NumDoFC)
{
    this->assign( s, f, t);
}

template<class T>
  template<class BndDataT>
    LocalP2CL<T>::LocalP2CL(const TetraCL& s,
        const VecDescCL& vd, const BndDataT& bnd, double t)
    : base_type( value_type(), FE_P2CL::NumDoFC)
{
    this->assign( s, vd, bnd, t);
}

template<class T>
  LocalP2CL<T>::LocalP2CL (const LocalP1CL<T>& p1, double t)
    : base_type( value_type(), FE_P2CL::NumDoFC)
{
    this->assign( p1, t);
}

template<class T>
  inline typename LocalP2CL<T>::value_type
  LocalP2CL<T>::operator() (const BaryCoordCL& p) const
{
    return FE_P2CL::val( *this, p);
}


//**************************************************************************
// Class: Quad2CL                                                          *
//**************************************************************************
template<class T>
const double Quad2CL<T>::Node[5][4]= {
    {1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}, {.25,.25,.25,.25}
};

template<class T>
const double Quad2CL<T>::Wght[5]= {
    1./120., 1./120., 1./120., 1./120., 2./15.
};

template<class T>
  BaryCoordCL*
  Quad2CL<T>::TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p)
{
    if (!p) p= new BaryCoordCL[NumNodesC];
    //tN[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
    std::memcpy( p, M.begin(), 4*sizeof( BaryCoordCL));
    for (Uint k= 0; k < 4; ++k)
            p[4][k]= 0.25*( M[0][k] + M[1][k] + M[2][k] + M[3][k]);
    return p;
}

template<class T>
  inline Quad2CL<T>&
  Quad2CL<T>::assign(const TetraCL& s, instat_fun_ptr f , double t)
{
    for (Uint i= 0; i<NumNodesC-1; ++i)
        (*this)[i]= f( s.GetVertex( i)->GetCoord(), t);
    (*this)[NumNodesC-1]= f( GetBaryCenter( s), t);
    return *this;
}

template<class T>
  inline Quad2CL<T>&
  Quad2CL<T>::assign(const LocalP1CL<value_type>& f)
{
    (*this)[std::slice( 0, 4, 1)]= f;
    (*this)[NumNodesC-1]= f( BaryCoordCL( 0.25));
    return *this;
}

template<class T>
  inline Quad2CL<T>&
  Quad2CL<T>::assign(const LocalP2CL<value_type>& f)
{
    (*this)[std::slice( 0, 4, 1)]= f[std::slice( 0, 4, 1)];
    (*this)[NumNodesC-1]= f( BaryCoordCL( 0.25));
    return *this;
}

template<class T>
  inline Quad2CL<T>&
  Quad2CL<T>::assign(const LocalP2CL<value_type>& f, const BaryCoordCL* const node)
{
    for (size_t i= 0; i < Quad2CL<T>::NumNodesC; ++i)
        (*this)[i]= f( node[i]);
    return *this;
}

template<class T>
  template <class P2FunT>
    inline Quad2CL<T>&
    Quad2CL<T>::assign(const TetraCL& s, const P2FunT& f, double t)
{
    const double oldt= f.GetTime();
    f.SetTime( t);
    for (Uint i= 0; i<NumNodesC-1; ++i)
        (*this)[i]= f.val( *s.GetVertex( i));
    (*this)[NumNodesC-1]= f.val( s, 0.25, 0.25, 0.25);
    f.SetTime( oldt);
    return *this;
}

template<class T>
  Quad2CL<T>::Quad2CL(const TetraCL& s,
      instat_fun_ptr f, double t)
  : base_type( value_type(), NumNodesC)
{
    this->assign( s, f, t);
}

template<class T>
  Quad2CL<T>::Quad2CL(const LocalP2CL<value_type>& f)
  : base_type( value_type(), NumNodesC)
{
    this->assign( f);
}

template<class T>
  Quad2CL<T>::Quad2CL(const LocalP2CL<value_type>& f, const BaryCoordCL* const node)
{
    this->assign( f, node);

}

template<class T>
  template <class PFunT>
    Quad2CL<T>::Quad2CL(const TetraCL& s, const PFunT& f, double t)
  : base_type( value_type(), NumNodesC)
{
    this->assign( s, f, t);
}


template<class T>
  T
  Quad2CL<T>::quadP1D(int i, double absdet) const
{
    const double* const tmp( FE_P1DCL::VertexVal( i));
    Quad2CL<T> tmp2( *this);
    for (Uint k= 0; k<4; ++k)
        tmp2[k]*= tmp[k];
    tmp2[4]*= 0.25; // The value in the barycenter, independent of i.
    return tmp2.quad( absdet);
}

template<class T>
  T
  Quad2CL<T>::quadP1D (int i, int j, double absdet) const
{
    const std::valarray<double> tmp(  FE_P1DCL::VertexVal( i), NumVertsC);
    const std::valarray<double> tmp3( FE_P1DCL::VertexVal( j), NumVertsC);
    Quad2CL<T> tmp2( *this);
    tmp2[std::slice(0, 4, 1)]*= tmp*tmp3;
    tmp2[4]*= 0.0625; // The value in the barycenter, independent of i, j.
    return tmp2.quad( absdet);

}

//**************************************************************************
// Class: Quad3CL                                                          *
//**************************************************************************

template<class T>
  inline Quad3CL<T>&
  Quad3CL<T>::assign(const TetraCL& s, instat_fun_ptr f , double t, const BaryCoordCL* const node)
{
    for (Uint i= 0; i < Quad3DataCL::NumNodesC; ++i)
        (*this)[i]= f( GetWorldCoord( s, node[i]), t);
    return *this;
}

template<class T>
  inline Quad3CL<T>&
  Quad3CL<T>::assign(const LocalP1CL<value_type>& f, const BaryCoordCL* const node)
{
    for (size_t i= 0; i < Quad3DataCL::NumNodesC; ++i)
        (*this)[i]= f( node[i]);
    return *this;
}

template<class T>
  inline Quad3CL<T>&
  Quad3CL<T>::assign(const LocalP2CL<value_type>& f)
{
    (*this)= f[0]*static_cast<GridFunctionCL<double> >( Quad3DataCL::P2_Val[0]);
    for (size_t i= 1; i < 10; ++i)
        (*this)+= f[i]*static_cast<GridFunctionCL<double> >( Quad3DataCL::P2_Val[i]);
    return *this;
}

template<class T>
  inline Quad3CL<T>&
  Quad3CL<T>::assign(const LocalP2CL<value_type>& f, const BaryCoordCL* const node)
{
    for (size_t i= 0; i < Quad3DataCL::NumNodesC; ++i)
        (*this)[i]= f( node[i]);
    return *this;
}

template<class T>
  template <class _BndData, class _VD>
      inline Quad3CL<T>&
      Quad3CL<T>::assign(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>& f, double t)
{
    const double oldt= f.GetTime();
    f.SetTime( t);
    value_type dof[10];
    f.GetDoF( s, dof);
    (*this)= T();
    for (size_t i= 0; i < Quad3DataCL::NumNodesC; ++i)
        for (size_t j= 0; j < 10; ++j)
            (*this)[i]+= dof[j]*Quad3DataCL::P2_Val[j][i];
    f.SetTime( oldt);
    return *this;
}

template<class T>
  template <class PFunT>
    inline Quad3CL<T>&
    Quad3CL<T>::assign(const TetraCL& s, const PFunT& f, double t)
{
    const double oldt= f.GetTime();
    f.SetTime( t);
    for (size_t i= 0; i < Quad3DataCL::NumNodesC; ++i)
        (*this)[i]= f.val( s, Quad3DataCL::Node[i]);
    f.SetTime( oldt);
    return *this;
}

template<class T>
  Quad3CL<T>::Quad3CL(const TetraCL& s,
                      instat_fun_ptr f, double t, const BaryCoordCL* const node)
  : base_type( value_type(), Quad3DataCL::NumNodesC)
{
    this->assign( s, f, t, node);
}

template<class T>
  Quad3CL<T>::Quad3CL(const LocalP1CL<value_type>& f, const BaryCoordCL* const node)
    : base_type( value_type(), Quad3DataCL::NumNodesC)
{
    this->assign( f, node);
}

template<class T>
  Quad3CL<T>::Quad3CL(const LocalP2CL<value_type>& f)
  : base_type( value_type(), Quad3DataCL::NumNodesC)
{
    this->assign( f);
}

template<class T>
  Quad3CL<T>::Quad3CL(const LocalP2CL<value_type>& f, const BaryCoordCL* const node)
    : base_type( value_type(), Quad3DataCL::NumNodesC)
{
    this->assign( f, node);
}

template<class T>
  template <class _BndData, class _VD>
    Quad3CL<T>::Quad3CL(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>& f, double t)
    : base_type( value_type(), Quad3DataCL::NumNodesC)
{
    this->assign( s, f, t);
}

template<class T>
  template <class PFunT>
    Quad3CL<T>::Quad3CL(const TetraCL& s, const PFunT& f, double t)
    : base_type( value_type(), Quad3DataCL::NumNodesC)
{
    this->assign( s, f, t);
}

template<class T>
  inline T Quad3CL<T>::quad (double absdet) const
{
    return absdet*(
          -2./15. *  (*this)[0] + 3./40. * ((*this)[1]+(*this)[2]+(*this)[3]+(*this)[4])
    );
}

    // Quadraturformel zur Annaeherung von \int f*phi, phi = P2-Hutfunktion
template<class T>
  inline T Quad3CL<T>::quadP2 (int i, double absdet) const
{
    typedef Quad3DataCL Q;
    return absdet*(
          -2./15.*Q::P2_Val[i][0]*(*this)[0]
        + 3./40.*(Q::P2_Val[i][1]*(*this)[1]+Q::P2_Val[i][2]*(*this)[2]+Q::P2_Val[i][3]*(*this)[3]+Q::P2_Val[i][4]*(*this)[4])
    );
}

//**************************************************************************
// Class: Quad5CL                                                          *
//**************************************************************************

template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const TetraCL& s, instat_fun_ptr f , double t, const BaryCoordCL* const node)
{
    for (Uint i= 0; i < Quad5DataCL::NumNodesC; ++i)
        (*this)[i]= f( GetWorldCoord( s, node[i]), t);
    return *this;
}

template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const LocalP1CL<value_type>& f, const BaryCoordCL* const node)
{
    for (size_t i= 0; i < Quad5DataCL::NumNodesC; ++i)
        (*this)[i]= f( node[i]);
    return *this;
}

template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const LocalP2CL<value_type>& f)
{
    (*this)= f[0]*static_cast<GridFunctionCL<double> >( Quad5DataCL::P2_Val[0]);
    for (size_t i= 1; i < 10; ++i)
        (*this)+= f[i]*static_cast<GridFunctionCL<double> >( Quad5DataCL::P2_Val[i]);
    return *this;
}

template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const LocalP2CL<value_type>& f, const BaryCoordCL* const node)
{
    for (size_t i= 0; i < Quad5DataCL::NumNodesC; ++i)
        (*this)[i]= f( node[i]);
    return *this;
}

template<class T>
  template <class _BndData, class _VD>
      inline Quad5CL<T>&
      Quad5CL<T>::assign(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>& f, double t)
{
    const double oldt= f.GetTime();
    f.SetTime( t);
    value_type dof[10];
    f.GetDoF( s, dof);
    (*this)= T();
    for (size_t i= 0; i < Quad5DataCL::NumNodesC; ++i)
        for (size_t j= 0; j < 10; ++j)
            (*this)[i]+= dof[j]*Quad5DataCL::P2_Val[j][i];
    f.SetTime( oldt);
    return *this;
}

template<class T>
  template <class PFunT>
    inline Quad5CL<T>&
    Quad5CL<T>::assign(const TetraCL& s, const PFunT& f, double t)
{
    const double oldt= f.GetTime();
    f.SetTime( t);
    for (size_t i= 0; i < Quad5DataCL::NumNodesC; ++i)
        (*this)[i]= f.val( s, Quad5DataCL::Node[i]);
    f.SetTime( oldt);
    return *this;
}

template<class T>
  Quad5CL<T>::Quad5CL(const TetraCL& s,
                      instat_fun_ptr f, double t, const BaryCoordCL* const node)
  : base_type( value_type(), Quad5DataCL::NumNodesC)
{
    this->assign( s, f, t, node);
}

template<class T>
  Quad5CL<T>::Quad5CL(const LocalP1CL<value_type>& f, const BaryCoordCL* const node)
    : base_type( value_type(), Quad5DataCL::NumNodesC)
{
    this->assign( f, node);
}

template<class T>
  Quad5CL<T>::Quad5CL(const LocalP2CL<value_type>& f)
  : base_type( value_type(), Quad5DataCL::NumNodesC)
{
    this->assign( f);
}

template<class T>
  Quad5CL<T>::Quad5CL(const LocalP2CL<value_type>& f, const BaryCoordCL* const node)
    : base_type( value_type(), Quad5DataCL::NumNodesC)
{
    this->assign( f, node);
}

template<class T>
  template <class _BndData, class _VD>
    Quad5CL<T>::Quad5CL(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>& f, double t)
    : base_type( value_type(), Quad5DataCL::NumNodesC)
{
    this->assign( s, f, t);
}

template<class T>
  template <class PFunT>
    Quad5CL<T>::Quad5CL(const TetraCL& s, const PFunT& f, double t)
    : base_type( value_type(), Quad5DataCL::NumNodesC)
{
    this->assign( s, f, t);
}

template<class T>
  inline T Quad5CL<T>::quad (double absdet) const
{
    return absdet*(
          8./405. *  (*this)[0]
        + (2665.0 + 14.0*std::sqrt( 15.0))/226800.0 * ((*this)[1]+(*this)[2]+(*this)[3]+(*this)[4])
        + (2665.0 - 14.0*std::sqrt( 15.0))/226800.0 * ((*this)[5]+(*this)[6]+(*this)[7]+(*this)[8])
        + 5./567.* ((*this)[9]+(*this)[10]+(*this)[11]+(*this)[12]+(*this)[13]+(*this)[14])
    );
}

    // Quadraturformel zur Annaeherung von \int f*phi, phi = P2-Hutfunktion
template<class T>
  inline T Quad5CL<T>::quadP2 (int i, double absdet) const
{
    typedef Quad5DataCL Q;
    return absdet*(
          8./405.*Q::P2_Val[i][0]*(*this)[0]
        + (2665.0 + 14.0*std::sqrt( 15.0))/226800.0*(Q::P2_Val[i][1]*(*this)[1]+Q::P2_Val[i][2]*(*this)[2]+Q::P2_Val[i][3]*(*this)[3]+Q::P2_Val[i][4]*(*this)[4])
        + (2665.0 - 14.0*std::sqrt( 15.0))/226800.0*(Q::P2_Val[i][5]*(*this)[5]+Q::P2_Val[i][6]*(*this)[6]+Q::P2_Val[i][7]*(*this)[7]+Q::P2_Val[i][8]*(*this)[8])
        + 5./567.*(Q::P2_Val[i][9]*(*this)[9]+Q::P2_Val[i][10]*(*this)[10]+Q::P2_Val[i][11]*(*this)[11]+Q::P2_Val[i][12]*(*this)[12]+Q::P2_Val[i][13]*(*this)[13]+Q::P2_Val[i][14]*(*this)[14])
    );
}


//**************************************************************************
// Class: Quad5_2DCL                                                       *
//**************************************************************************
template<class T>
  inline Quad5_2DCL<T>&
  Quad5_2DCL<T>::assign(const TetraCL& s, const BaryCoordCL* const p, instat_fun_ptr f , double t)
{
    Point3DCL v[3];
    for (Uint k= 0; k < 3; ++k)
        v[k]= GetWorldCoord(s, p[k]);
    for (Uint i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
        (*this)[i]= f( Quad5_2DDataCL::Node[i][0]*v[0]+Quad5_2DDataCL::Node[i][1]*v[1]+Quad5_2DDataCL::Node[i][2]*v[2], t) ;
    return *this;
}

template<class T>
  inline Quad5_2DCL<T>&
  Quad5_2DCL<T>::assign(const LocalP1CL<value_type>& f, const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    SetInterface( p, NodeInTetra);
    for (size_t i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
        (*this)[i]= f( NodeInTetra[i]);
    return *this;
}

template<class T>
  inline Quad5_2DCL<T>&
  Quad5_2DCL<T>::assign(const LocalP2CL<value_type>& f, const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    SetInterface( p, NodeInTetra);
    std::valarray<double> P2_Val[10];
    FE_P2CL::ApplyAll( Quad5_2DDataCL::NumNodesC, NodeInTetra, P2_Val);
    (*this)= f[0]*(*static_cast<GridFunctionCL<>* >( &P2_Val[0])); // XXX: So, oder mit Kopierkonstruktor?
    for (size_t i= 1; i < 10; ++i)
        (*this)+= f[i]*(*static_cast<GridFunctionCL<>* >( &P2_Val[i]));
    return *this;
}

template<class T>
  template <class _BndData, class _VD>
    inline Quad5_2DCL<T>&
    Quad5_2DCL<T>::assign(const TetraCL& s, const BaryCoordCL* const p,
        const P2EvalCL<T, _BndData, _VD>& f, double t)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    SetInterface( p, NodeInTetra);
    std::valarray<double> P2_Val[10];
    FE_P2CL::ApplyAll( Quad5_2DDataCL::NumNodesC, NodeInTetra, P2_Val);
    const double oldt= f.GetTime();
    f.SetTime( t);
    value_type dof[10];
    f.GetDoF( s, dof);
    (*this)= T();
    for (size_t i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
        for (size_t j= 0; j < 10; ++j)
            (*this)[i]+= dof[j]*P2_Val[j][i];
    f.SetTime( oldt);
    return *this;
}

template<class T>
  template <class PFunT>
    inline Quad5_2DCL<T>&
    Quad5_2DCL<T>::assign(const TetraCL& s, const BaryCoordCL* const p, const PFunT& f, double t)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    SetInterface( p, NodeInTetra);
    const double oldt= f.GetTime();
    f.SetTime( t);
    for (size_t i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
        (*this)[i]= f.val( s, NodeInTetra[i]);
    f.SetTime( oldt);
    return *this;
}

template<class T>
  Quad5_2DCL<T>::Quad5_2DCL(const TetraCL& s, const BaryCoordCL*const p, instat_fun_ptr f, double t)
  : base_type( value_type(), Quad5_2DDataCL::NumNodesC)
{
    this->assign( s, p, f, t);
}

template<class T>
  Quad5_2DCL<T>::Quad5_2DCL(const LocalP2CL<value_type>& f, const BaryCoordCL* const p)
  : base_type( value_type(), Quad5_2DDataCL::NumNodesC)
{
    this->assign( f, p);
}

template<class T>
  template <class _BndData, class _VD>
    Quad5_2DCL<T>::Quad5_2DCL(const TetraCL& s, const BaryCoordCL* const p,
        const P2EvalCL<T, _BndData, _VD>& f, double t)
    : base_type( value_type(), Quad5_2DDataCL::NumNodesC)
{
    this->assign( s, p, f, t);
}

template<class T>
  template <class PFunT>
    Quad5_2DCL<T>::Quad5_2DCL(const TetraCL& s, const BaryCoordCL* const p, const PFunT& f, double t)
    : base_type( value_type(), Quad5_2DDataCL::NumNodesC)
{
    this->assign( s, p, f, t);
}


} // end of namespace DROPS
