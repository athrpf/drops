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
    const Uint tlvl= s.GetLevel();
    const Uint vlvl= vd.GetLevel();
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
    const Uint tlvl= s.GetLevel();
    const Uint flvl= f.GetLevel();
    const double tmp= f.GetTime();
    const_cast<P1FunT&>( f).SetTime( t);
    f.GetDoF( s, *this);
    const_cast<P1FunT&>( f).SetTime( tmp);
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
    const_cast<P2FunT&>( f).SetTime( t);
    if (tlvl == flvl)
        f.GetDoF( s, *this);
    else
        if (tlvl < flvl) RestrictP2( s, f, *this);
        else throw DROPSErrCL( "LocalP2CL::Assign: Prolongation not implemented.\n");
    const_cast<P2FunT&>( f).SetTime( tmp);
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
  template <class P2FunT>
    inline Quad2CL<T>&
    Quad2CL<T>::assign(const TetraCL& s, const P2FunT& f, double t)
{
    const double oldt= f.GetTime();
    const_cast<P2FunT&>( f).SetTime( t);
    for (Uint i= 0; i<NumNodesC-1; ++i)
        (*this)[i]= f.val( *s.GetVertex( i));
    (*this)[NumNodesC-1]= f.val( s, 0.25, 0.25, 0.25);
    const_cast<P2FunT&>( f).SetTime( oldt);
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
// Class: Quad5CL                                                          *
//**************************************************************************

template<class T>
bool Quad5CL<T>::HaveNodes= false;

template <class T>
BaryCoordCL Quad5CL<T>::Node[NumNodesC];

template<class T>
const double Quad5CL<T>::Wght[4]= {
    8./405.,                                   /*Node[0]*/
    (2665.0 + 14.0*std::sqrt( 15.0))/226800.0, /*Node[1] bis Node[4]*/
    (2665.0 - 14.0*std::sqrt( 15.0))/226800.0, /*Node[5] bis Node[8]*/
    5./567.                                    /*Node[9] bis Node[14]*/
};

template<class T>
std::valarray<double> Quad5CL<T>::P2_Val[10]; // P2_Val[i] contains FE_P2CL::H_i( Node).

template<class T>
  void
  Quad5CL<T>::MaybeInitNodes() const
{
    if (HaveNodes) return;

    Node[0]= BaryCoordCL( 0.25);
    const double A1= (7.0 - std::sqrt( 15.0))/34.0,
                 B1= (13.0 + 3.0*std::sqrt( 15.0))/34.0;
    Node[1]= MakeBaryCoord( A1,A1,A1,B1);
    Node[2]= MakeBaryCoord( A1,A1,B1,A1);
    Node[3]= MakeBaryCoord( A1,B1,A1,A1);
    Node[4]= MakeBaryCoord( B1,A1,A1,A1);
    const double A2= (7.0 + std::sqrt( 15.0))/34.0,
                 B2= (13.0 - 3.0*std::sqrt( 15.0))/34.0;
    Node[5]= MakeBaryCoord( A2,A2,A2,B2);
    Node[6]= MakeBaryCoord( A2,A2,B2,A2);
    Node[7]= MakeBaryCoord( A2,B2,A2,A2);
    Node[8]= MakeBaryCoord( B2,A2,A2,A2);
    const double A3= (10.0 - 2.0*std::sqrt( 15.0))/40.0,
                 B3= (10.0 + 2.0*std::sqrt( 15.0))/40.0;
    Node[9] = MakeBaryCoord( A3,A3,B3,B3);
    Node[10]= MakeBaryCoord( A3,B3,A3,B3);
    Node[11]= MakeBaryCoord( A3,B3,B3,A3);
    Node[12]= MakeBaryCoord( B3,A3,A3,B3);
    Node[13]= MakeBaryCoord( B3,A3,B3,A3);
    Node[14]= MakeBaryCoord( B3,B3,A3,A3);

    FE_P2CL::ApplyAll( NumNodesC, Node, P2_Val);

    HaveNodes= true;
}


template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const TetraCL& s, instat_fun_ptr f , double t)
{
    (*this)[0]= f( GetBaryCenter( s), t);
    for (Uint i= 1; i < NumNodesC; ++i)
        (*this)[i]= f( GetWorldCoord( s, Node[i]), t);
    return *this;
}

template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const LocalP1CL<value_type>& f)
{
    for (size_t i= 0; i < NumNodesC; ++i)
        (*this)[i]= f( Node[i]);
    return *this;
}

template<class T>
  inline Quad5CL<T>&
  Quad5CL<T>::assign(const LocalP2CL<value_type>& f)
{
    (*this)= f[0]*P2_Val[0];
    for (size_t i= 1; i < 10; ++i)
        (*this)+= f[i]*P2_Val[i];
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
    for (size_t i= 0; i < NumNodesC; ++i)
        for (size_t j= 0; j < 10; ++j)
            (*this)[i]+= dof[j]*P2_Val[j][i];
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
    for (size_t i= 0; i < NumNodesC; ++i)
        (*this)[i]= f.val( s, Node[i]);
    f.SetTime( oldt);
    return *this;
}

template<class T>
  Quad5CL<T>::Quad5CL(const TetraCL& s,
      instat_fun_ptr f, double t)
  : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( s, f, t);
}

template<class T>
  Quad5CL<T>::Quad5CL(const LocalP2CL<value_type>& f)
  : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( f);
}

template<class T>
  template <class _BndData, class _VD>
    Quad5CL<T>::Quad5CL(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>& f, double t)
    : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( s, f, t);
}

template<class T>
  template <class PFunT>
    Quad5CL<T>::Quad5CL(const TetraCL& s, const PFunT& f, double t)
    : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( s, f, t);
}


//**************************************************************************
// Class: Quad5_2DCL                                                       *
//**************************************************************************
template<class T>
  bool Quad5_2DCL<T>::HaveNodes= false;

template <class T>
  Point3DCL Quad5_2DCL<T>::Node[NumNodesC];  // Barycentric coord for 2D

template<class T>
  const double Quad5_2DCL<T>::Wght[3]= {
      9./80.,                                   /*Node[0]*/
      (155. - std::sqrt( 15.0))/2400., /*Node[1] to Node[3]*/
      (155. + std::sqrt( 15.0))/2400., /*Node[4] to Node[6]*/
};

template<class T>
  void
  Quad5_2DCL<T>::MaybeInitNodes() const
{
    if (HaveNodes) return;

    Node[0]= Point3DCL( 1./3.);
    const double A1= (6.0 - std::sqrt( 15.0))/21.0,
    B1= (9.0 + 2.0*std::sqrt( 15.0))/21.0;
    Node[1]= MakePoint3D( A1,A1,B1);
    Node[2]= MakePoint3D( A1,B1,A1);
    Node[3]= MakePoint3D( B1,A1,A1);
    const double A2= (6.0 + std::sqrt( 15.0))/21.0,
    B2= (9.0 - 2.0*std::sqrt( 15.0))/21.0;
    Node[4]= MakePoint3D( A2,A2,B2);
    Node[5]= MakePoint3D( A2,B2,A2);
    Node[6]= MakePoint3D( B2,A2,A2);
    HaveNodes= true;
}

template<class T>
  void
  Quad5_2DCL<T>::SetInterface(const BaryCoordCL*const p, BaryCoordCL* NodeInTetra)
{
    for (Uint i= 0; i < NumNodesC; ++i)
        NodeInTetra[i]= Node[i][0]*p[0] + Node[i][1]*p[1] + Node[i][2]*p[2];
}

template<class T>
  inline Quad5_2DCL<T>&
  Quad5_2DCL<T>::assign(const TetraCL& s, const BaryCoordCL* const p, instat_fun_ptr f , double t)
{
    Point3DCL v[3];
    for (Uint k= 0; k < 3; ++k)
        v[k]= GetWorldCoord(s, p[k]);
    for (Uint i= 0; i < NumNodesC; ++i)
        (*this)[i]= f( Node[i][0]*v[0]+Node[i][1]*v[1]+Node[i][2]*v[2], t) ;
    return *this;
}

template<class T>
  inline Quad5_2DCL<T>&
  Quad5_2DCL<T>::assign(const LocalP1CL<value_type>& f, const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[NumNodesC];
    SetInterface( p, NodeInTetra);
    for (size_t i= 0; i < NumNodesC; ++i)
        (*this)[i]= f( NodeInTetra[i]);
    return *this;
}

template<class T>
  inline Quad5_2DCL<T>&
  Quad5_2DCL<T>::assign(const LocalP2CL<value_type>& f, const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[NumNodesC];
    SetInterface( p, NodeInTetra);
    std::valarray<double> P2_Val[10];
    FE_P2CL::ApplyAll( NumNodesC, NodeInTetra, P2_Val);
    (*this)= f[0]*P2_Val[0];
    for (size_t i= 1; i < 10; ++i)
        (*this)+= f[i]*P2_Val[i];
    return *this;
}

template<class T>
  template <class _BndData, class _VD>
    inline Quad5_2DCL<T>&
    Quad5_2DCL<T>::assign(const TetraCL& s, const BaryCoordCL* const p,
        const P2EvalCL<T, _BndData, _VD>& f, double t)
{
    BaryCoordCL NodeInTetra[NumNodesC];
    SetInterface( p, NodeInTetra);
    std::valarray<double> P2_Val[10];
    FE_P2CL::ApplyAll( NumNodesC, NodeInTetra, P2_Val);
    const double oldt= f.GetTime();
    f.SetTime( t);
    value_type dof[10];
    f.GetDoF( s, dof);
    (*this)= T();
    for (size_t i= 0; i < NumNodesC; ++i)
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
    BaryCoordCL NodeInTetra[NumNodesC];
    SetInterface( p, NodeInTetra);
    const double oldt= f.GetTime();
    f.SetTime( t);
    for (size_t i= 0; i < NumNodesC; ++i)
        (*this)[i]= f.val( s, NodeInTetra[i]);
    f.SetTime( oldt);
    return *this;
}

template<class T>
  Quad5_2DCL<T>::Quad5_2DCL(const TetraCL& s, const BaryCoordCL*const p, instat_fun_ptr f, double t)
  : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( s, p, f, t);
}

template<class T>
  Quad5_2DCL<T>::Quad5_2DCL(const LocalP2CL<value_type>& f, const BaryCoordCL* const p)
  : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( f, p);
}

template<class T>
  template <class _BndData, class _VD>
    Quad5_2DCL<T>::Quad5_2DCL(const TetraCL& s, const BaryCoordCL* const p,
        const P2EvalCL<T, _BndData, _VD>& f, double t)
    : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( s, p, f, t);
}

template<class T>
  template <class PFunT>
    Quad5_2DCL<T>::Quad5_2DCL(const TetraCL& s, const BaryCoordCL* const p, const PFunT& f, double t)
    : base_type( value_type(), NumNodesC)
{
    MaybeInitNodes();
    this->assign( s, p, f, t);
}


} // end of namespace DROPS
