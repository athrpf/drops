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
    return FE_P1DCL::val( *this, p);
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

} // end of namespace DROPS
