/// \file discretize.h
/// \brief discretizations for several PDEs and FE types
/// \author LNM RWTH Aachen: Patrick Esser, Sven Gross, Trung Hieu Nguyen, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_DISCRETIZE_H
#define DROPS_DISCRETIZE_H

#include <cstring>

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "misc/container.h"
#include "num/fe.h"

namespace DROPS
{

typedef double    (*scalar_fun_ptr)       (const Point3DCL&);
typedef Point3DCL (*vector_fun_ptr)       (const Point3DCL&);
typedef double    (*instat_scalar_fun_ptr)(const Point3DCL&, double);
typedef Point3DCL (*instat_vector_fun_ptr)(const Point3DCL&, double);
typedef bool      (*match_fun)        (const Point3DCL&, const Point3DCL&);

typedef double    (*SmoothFunT)           (double,double);


// SmoothedJumpCL for jumping coefficients

class JumpCL
{
  private:
    double Coeff[2];

  public:
    JumpCL (double In, double Out) { Coeff[0]=In; Coeff[1]=Out; }

    double operator() (bool b)   const { return Coeff[b]; }
    double operator() (double x) const { return (1-x)*Coeff[0]+x*Coeff[1]; }
};

class SmoothedJumpCL
{
  private:
    const JumpCL     jc;
    const SmoothFunT smf;
    const double     eps;

  public:
    SmoothedJumpCL (const JumpCL& myjc, const SmoothFunT f, double myeps)
      : jc(myjc), smf(f), eps(myeps) {}
    SmoothedJumpCL (double In, double Out, const SmoothedJumpCL& sjc)
      : jc(In,Out), smf(sjc.smf), eps(sjc.eps) {}

    double operator() (double x) const { return jc(smf(x,eps)); }
};

// smoothed Heaviside function
inline double H_sm( double s, double eps)
{
    if (s <= -eps) return 0;
    if (s >=  eps) return 1;
    // -eps < s < eps
    s/= eps;
    const double s2= s*s, s3= s2*s;
    return 0.5 + 1.40625*s - 1.5625*s3 + 0.65625*s2*s3;
}


//**************************************************************************
// Class:   GridFunctionCL                                                 *
// Template Parameter:                                                     *
//          T - The result-type of the finite-element-function             *
// Purpose: Base-class for LocalPiCL and Quadrature-classes. Defines       *
//          arithmetik via expression templates and a function-object bases*
//          apply.                                                         *
//          Defines mixed operations for Point3DCL and double as           *
//          value_types.                                                   *
//**************************************************************************
template<class T= double>
class GridFunctionCL: public std::valarray<T>
{
  public:
    typedef T value_type;
    typedef std::valarray<T> base_type;
    typedef value_type (*instat_fun_ptr)(const Point3DCL&, double);

  protected:
    typedef GridFunctionCL<T> self_;

  public:
    GridFunctionCL (value_type v, Uint s): base_type( v, s) {}
    GridFunctionCL (const value_type* p, Uint s): base_type( p, s) {}

DROPS_DEFINE_VALARRAY_DERIVATIVE(GridFunctionCL, T, base_type)

    template<typename FuncT>
      inline GridFunctionCL& apply(FuncT fun);

    template<typename U, typename MemberFuncT>
      inline GridFunctionCL& apply(U& obj, MemberFuncT fun);
};

template<class T>
  template<typename FuncT>
    inline GridFunctionCL<T>&
    GridFunctionCL<T>::apply(FuncT fun)
{
    for (size_t i= 0; i < this->size(); ++i)
        (*this)[i]= fun( (*this)[i]);
    return *this;
}

template<class T>
  template<typename U, typename MemberFuncT>
    inline GridFunctionCL<T>&
    GridFunctionCL<T>::apply(U& obj, MemberFuncT fun)
{
    for (size_t i= 0; i < this->size(); ++i)
        (*this)[i]= (obj.*fun)( (*this)[i]);
    return *this;
}

inline GridFunctionCL<Point3DCL>&
operator*=(GridFunctionCL<Point3DCL>& a, const GridFunctionCL<double>& b)
{
    for (size_t i= 0; i < b.size(); ++i)
        a[i]*= b[i];
    return a;
}

inline GridFunctionCL<Point3DCL>
operator*(const GridFunctionCL<Point3DCL>& a, const GridFunctionCL<double>& b)
{
    GridFunctionCL<Point3DCL> ret( a);
    return ret*= b;
}

inline GridFunctionCL<Point3DCL>
operator*(const GridFunctionCL<double>& a, const GridFunctionCL<Point3DCL>& b)
{
    return b*a;
}

inline GridFunctionCL<Point3DCL>&
operator*=(GridFunctionCL<Point3DCL>& a, double b)
{
    for (size_t i= 0; i < a.size(); ++i)
        a[i]*= b;
    return a;
}

inline GridFunctionCL<Point3DCL>
operator*( double a, const GridFunctionCL<Point3DCL>& b)
{

    GridFunctionCL<SVectorCL<3> > ret( b);
    return ret*= a;
}

inline GridFunctionCL<Point3DCL>
operator*( const GridFunctionCL<Point3DCL>& a, double b)
{

    return b*a;
}

inline GridFunctionCL<Point3DCL>
operator*( const Point3DCL& a, const GridFunctionCL<double>& b)
{

    GridFunctionCL<SVectorCL<3> > ret( a, b.size());
    for (size_t i= 0; i < b.size(); ++i)
        ret[i]*= b[i];
    return ret;
}

inline GridFunctionCL<Point3DCL>
operator*( const GridFunctionCL<double>& a, const Point3DCL& b)
{

    return b*a;
}

inline GridFunctionCL<double>
dot(const GridFunctionCL<Point3DCL>& a, const GridFunctionCL<Point3DCL>& b)
{
    GridFunctionCL<double> ret( 0.0, a.size());
    for (size_t i= 0; i<a.size(); ++i)
        ret[i]= inner_prod( a[i], b[i]);
    return ret;
}

inline GridFunctionCL<double>
dot(const Point3DCL& a, const GridFunctionCL<Point3DCL>& b)
{
    GridFunctionCL<double> ret( 0.0, b.size());
    for (size_t i= 0; i<b.size(); ++i)
        ret[i]= inner_prod( a, b[i]);
    return ret;
}

inline void ExtractComponent( const GridFunctionCL<Point3DCL>& src, GridFunctionCL<double>& target, int comp)
{
    Assert( target.size()==src.size(), DROPSErrCL("extractComponent: GridFunctionCL objects have different sizes"), DebugNumericC);
    for (size_t i= 0; i<src.size(); ++i)
        target[i]= src[i][comp];
}

inline GridFunctionCL< SMatrixCL<3,3> >
outer_product(const GridFunctionCL<Point3DCL>& a, const GridFunctionCL<Point3DCL>& b)
{
    GridFunctionCL< SMatrixCL<3,3> > ret( SMatrixCL<3,3>(), b.size());
    for (size_t i= 0; i<b.size(); ++i)
        ret[i]= outer_product( a[i], b[i]);
    return ret;
}

//**************************************************************************
// Class:   LocalP1CL                                                      *
// Template Parameter:                                                     *
//          T - The result-type of the finite-element-function             *
// Purpose: Evaluate a P1-function on a tetrahedron and calculate with such*
//          functions. As LocalP1CL is derived from valarray, arithmetic   *
//          operations are carried out efficiently.                        *
//          The valarray holds the values in the 4 degrees of freedom,     *
//          vertex_0,..., vertex_3.                                        *
//**************************************************************************
template<class T= double>
class LocalP1CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

  protected:
    typedef LocalP1CL<T> self_;

  public:
    LocalP1CL() : base_type( value_type(), FE_P1CL::NumDoFC) {}
    LocalP1CL(const value_type& t): base_type( t, FE_P1CL::NumDoFC) {}
    // Initialize from a given function
    LocalP1CL(const TetraCL&, instat_fun_ptr , double= 0.0);
    // Initialize from VecDescCL and boundary-data
    template<class BndDataT>
      LocalP1CL(const TetraCL&, const VecDescCL&, const BndDataT&);
    // Initialize from PiEvalCl
    template <class P1FunT>
      LocalP1CL(const TetraCL&, const P1FunT&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(LocalP1CL, T, base_type)

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const TetraCL&, instat_fun_ptr, double= 0.0);
    template<class BndDataT>
      inline self_&
      assign(const TetraCL&, const VecDescCL&, const BndDataT&);
    template <class P1FunT>
      inline self_&
      assign(const TetraCL&, const P1FunT&);

    // pointwise evaluation in barycentric coordinates
    inline value_type operator()(const BaryCoordCL&) const;
};

//**************************************************************************
// Class:   LocalP2CL                                                      *
// Template Parameter:                                                     *
//          T - The result-type of the finite-element-function             *
// Purpose: Evaluate a P2-function on a tetrahedron and calculate with such*
//          functions. As LocalP2CL is derived from valarray, arithmetic   *
//          operations are carried out efficiently.                        *
//          The valarray holds the values in the 10 degrees of freedom,    *
//          vertex_0,..., vertex_3, edge_0,..., edge5.                     *
//**************************************************************************
template<class T= double>
class LocalP2CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

  protected:
    typedef LocalP2CL<T> self_;

  public:
    LocalP2CL() : base_type( value_type(), FE_P2CL::NumDoFC) {}
    LocalP2CL(const value_type& t): base_type( t, FE_P2CL::NumDoFC) {}
    // Initialize from a given function
    LocalP2CL(const TetraCL&, instat_fun_ptr , double= 0.0);
    // Initialize from VecDescCL and boundary-data
    template<class BndDataT>
      LocalP2CL(const TetraCL&, const VecDescCL&, const BndDataT&);
    // Initialize from PiEvalCL
    template <class P2FunT>
      LocalP2CL(const TetraCL&, const P2FunT&);
    // Initialize from LocalP1CL
    LocalP2CL(const LocalP1CL<T>&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(LocalP2CL, T, base_type)

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const TetraCL&, instat_fun_ptr, double= 0.0);
    template<class BndDataT>
      inline self_&
      assign(const TetraCL&, const VecDescCL&, const BndDataT&);
    template <class P2FunT>
      inline self_&
      assign(const TetraCL&, const P2FunT&);
    inline self_&
    assign(const LocalP1CL<T>&);

    // pointwise evaluation in barycentric coordinates
    inline value_type operator()(const BaryCoordCL&) const;
};


/// \brief Extends/Interpolates P1 function on regular child to P1 function on parent.
///
/// \param isoP2 data interpreted as isoP2 function (i.e., P1 on each child) to be extended
/// \param child index 0,..,7 of regular child
/// \param P1onParent after function call, contains interpolated P1 values on whole parent extended from child \a ch
template<class T>
void ExtendP1onChild( const LocalP2CL<T>& isoP2, int child, LocalP2CL<T>& P1onParent);

// ===================================
//        Quadrature formulas
// ===================================
/// \brief Contains the nodes and weights of a positive quadrature rule on the reference tetrahedron. It uses 5 nodes an is exact up to degree 2.
///
/// The data is initialized exactly once on program-startup by the global object in num/discretize.cpp.
class Quad2DataCL
{
  public:
    Quad2DataCL ();

    enum { NumNodesC= 5 };

    static BaryCoordCL  Node[NumNodesC]; ///< quadrature nodes
    static const double Wght[2];         ///< quadrature weights
    static const double Weight[NumNodesC];///< quadrature weight for each node

    /// \param M contains the barycentric coordinates of a tetrahedron;
    /// \param p array to be used for the quadrature points for this tetrahedron.
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    static BaryCoordCL* TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p= 0);
};

template<class T=double>
class Quad2CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

    static const double Node[5][4]; // Stuetzstellen (NumNodesC*4 doubles)
    static const double Wght[5];    // Gewichte      (NumNodesC   doubles)

    /// \param M contains the barycentric coordinates of a tetrahedron;
    /// \param p array to be used for the quadrature points for this tetrahedron.
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    static BaryCoordCL* TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p= 0);

  protected:
    typedef Quad2CL<T> self_;

  public:
    Quad2CL(): base_type( value_type(), Quad2DataCL::NumNodesC) {}
    Quad2CL(const value_type& t): base_type( t, Quad2DataCL::NumNodesC) {}

    Quad2CL(const TetraCL&, instat_fun_ptr, double= 0.0);
    Quad2CL(const LocalP2CL<value_type>&);
    Quad2CL(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class PFunT>
      Quad2CL(const TetraCL&, const PFunT&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(Quad2CL, T, base_type)

    inline self_&
    assign(const TetraCL&, instat_fun_ptr , double= 0.0);
    inline self_&
    assign(const LocalP1CL<value_type>&);
    inline self_&
    assign(const LocalP2CL<value_type>&);
    inline self_&
    assign(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class P2FunT>
      inline self_&
      assign(const TetraCL&, const P2FunT&);

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const
    {
        value_type sum= this->sum()/120.;
        return (sum + 0.125*(*this)[Quad2DataCL::NumNodesC-1])*absdet;
    }

    // Folgende Spezialformeln nutzen die spezielle Lage der Stuetzstellen aus
    // zur Annaeherung von \int f*phi,    phi = P1-/P1D-/P2-Hutfunktion
    T quadP1 (int i, double absdet) const
      { return ((1./120.)*(*this)[i] + (1./30.)*(*this)[4])*absdet; }
    T quadP1 (int i, int j, double absdet) const
      { return (i!=j ? (1./720.)*((*this)[i]+(*this)[j]) + (1./180.)*(*this)[4]
                     : (1./180.)*(*this)[i] + (1./90.)*(*this)[4]  )*absdet;}
    // Die P1D-Formeln sind nur exakt bis Grad 1 bzw. 0.
    T quadP1D (int i, double absdet) const;
    T quadP1D (int i, int j, double absdet) const;

    T quadP2 (int i, double absdet) const
    {
        return (i<4 ? (1./360.)*(*this)[i] - (1./90.)*(*this)[4]
                    : (1./180.)*((*this)[VertOfEdge(i-4,0)]+(*this)[VertOfEdge(i-4,1)]) + (1./45.)*(*this)[4]
               )*absdet;
    }

    T quadP2 (int i, int j, double absdet) const
    {
        const double valBary= (i<4 ? -0.125 : 0.25)*(j<4 ? -0.125 : 0.25);
        return ((i!=j || i>=4) ? Wght[4]*(*this)[4]*valBary
                               : Wght[4]*(*this)[4]*valBary + Wght[i]*(*this)[i]
               )*absdet;
    }
};

/// \brief Contains the nodes and weights of a quadrature rule on the reference tetrahedron. It uses 5 nodes an is exact up to degree 3.
///
/// The data is initialized exactly once on program-startup by the global object in num/discretize.cpp.
class Quad3DataCL
{
  public:
    Quad3DataCL ();

    enum { NumNodesC= 5};

    static BaryCoordCL           Node[NumNodesC]; ///< quadrature nodes
    static const double          Wght[2];         ///< quadrature weights
    static const double          Weight[NumNodesC];///< quadrature weight for each node
    static std::valarray<double> P2_Val[10];      ///< P2_Val[i] contains FE_P2CL::H_i( Node).

    /// \param M contains the barycentric coordinates of a tetrahedron;
    /// \param p array to be used for the quadrature points for this tetrahedron.
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    static BaryCoordCL* TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p= 0);
};


/// \brief Quadrature rule on a tetrahedron of degree 3 with 5 nodes.
///
/// The numerical data for the actual quadrature is contained in Quad3DataCL -- it does not depend on the template parameter.
template<class T=double>
class Quad3CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

  protected:
    typedef Quad3CL<T> self_;

  public:
    Quad3CL(): base_type( value_type(), Quad3DataCL::NumNodesC) {}
    Quad3CL(const value_type& t): base_type( t, Quad3DataCL::NumNodesC) {}

    Quad3CL(const TetraCL&, instat_fun_ptr, double= 0.0, const BaryCoordCL* const= Quad3DataCL::Node);
    Quad3CL(const LocalP1CL<value_type>&, const BaryCoordCL* const= Quad3DataCL::Node);
    Quad3CL(const LocalP2CL<value_type>&);
    Quad3CL(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class _BndData, class _VD>
      Quad3CL(const TetraCL&, const P2EvalCL<T, _BndData, _VD>&);
    template <class PFunT>
      Quad3CL(const TetraCL&, const PFunT&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(Quad3CL, T, base_type)

    inline self_&
    assign(const TetraCL&, instat_fun_ptr , double= 0.0, const BaryCoordCL* const= Quad3DataCL::Node);
    inline self_&
    assign(const LocalP1CL<value_type>&, const BaryCoordCL* const= Quad3DataCL::Node);
    inline self_&
    assign(const LocalP2CL<value_type>&);
    inline self_&
    assign(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class _BndData, class _VD>
      inline self_&
      assign(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>&);
    template <class PFunT>
      inline self_&
      assign(const TetraCL&, const PFunT&);

    /// \param M contains the barycentric coordinates of a tetrahedron;
    /// \param p array to be used for the quadrature points for this tetrahedron.
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    static BaryCoordCL* TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p= 0) {
        return Quad3DataCL::TransformNodes( M, p);
    }

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    inline T quad (double absdet) const;

    // Quadraturformel zur Annaeherung von \int f*phi, phi = P2-Hutfunktion
    inline T quadP2 (int i, double absdet) const;
};

/// \brief Contains the nodes and weights of a positive quadrature rule on the reference tetrahedron. It uses 15 nodes an is exact up to degree 5.
///
/// The data is initialized exactly once on program-startup by the global object in num/discretize.cpp.
class Quad5DataCL
{
  public:
    Quad5DataCL ();

    enum { NumNodesC= 15 };

    static BaryCoordCL           Node[NumNodesC]; ///< quadrature nodes
    static const double          Wght[4];         ///< quadrature weights
    static const double          Weight[NumNodesC];///< quadrature weight for each node
    static std::valarray<double> P2_Val[10];      ///< P2_Val[i] contains FE_P2CL::H_i( Node).

    /// \param M contains the barycentric coordinates of a tetrahedron;
    /// \param p array to be used for the quadrature points for this tetrahedron.
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    static BaryCoordCL* TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p= 0);
};


/// \brief Positive quadrature rule on a tetrahedron of degree 5 with 15 nodes.
///
/// The numerical data for the actual quadrature is contained in Quad5DataCL -- it does not depend on the template parameter.
template<class T=double>
class Quad5CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

  protected:
    typedef Quad5CL<T> self_;

  public:
    Quad5CL(): base_type( value_type(), Quad5DataCL::NumNodesC) {}
    Quad5CL(const value_type& t): base_type( t, Quad5DataCL::NumNodesC) {}

    Quad5CL(const TetraCL&, instat_fun_ptr, double= 0.0, const BaryCoordCL* const= Quad5DataCL::Node);
    Quad5CL(const LocalP1CL<value_type>&, const BaryCoordCL* const= Quad5DataCL::Node);
    Quad5CL(const LocalP2CL<value_type>&);
    Quad5CL(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class _BndData, class _VD>
      Quad5CL(const TetraCL&, const P2EvalCL<T, _BndData, _VD>&);
    template <class PFunT>
      Quad5CL(const TetraCL&, const PFunT&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(Quad5CL, T, base_type)

    inline self_&
    assign(const TetraCL&, instat_fun_ptr , double= 0.0, const BaryCoordCL* const= Quad5DataCL::Node);
    inline self_&
    assign(const LocalP1CL<value_type>&, const BaryCoordCL* const= Quad5DataCL::Node);
    inline self_&
    assign(const LocalP2CL<value_type>&);
    inline self_&
    assign(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class _BndData, class _VD>
      inline self_&
      assign(const TetraCL& s, const P2EvalCL<T, _BndData, _VD>&);
    template <class PFunT>
      inline self_&
      assign(const TetraCL&, const PFunT&);

    /// \param M contains the barycentric coordinates of a tetrahedron;
    /// \param p array to be used for the quadrature points for this tetrahedron.
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    static BaryCoordCL* TransformNodes (const SArrayCL<BaryCoordCL,4>& M, BaryCoordCL* p= 0) {
        return Quad5DataCL::TransformNodes( M, p);
    }

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    inline T quad (double absdet) const;

    // Quadraturformel zur Annaeherung von \int f*phi, phi = P2-Hutfunktion
    inline T quadP2 (int i, double absdet) const;
};

/// \brief Contains the nodes and weights of a positive quadrature rule on the reference triangle. It uses 7 nodes an is exact up to degree 5.
///
/// The data is initialized exactly once on program-startup by the global object in num/discretize.cpp.
class Quad5_2DDataCL
{
  public:
    Quad5_2DDataCL ();

    enum { NumNodesC= 7 };

    static Point3DCL           Node[NumNodesC];   ///< quadrature nodes
    static const double        Wght[3];           ///< quadrature weights
    static const double        Weight[NumNodesC]; ///< quadrature weights for each node

    /// Calculates the barycentric coordinates of the quadrature points
    /// of the triangle given by the 1st argument with respect to the
    /// tetrahedron and stores them in the 2nd argument.
    /// RAIterT is a random-access-iterator to a sequence of BaryCoordCL
    template <class RAIterT>
    static void SetInterface (const BaryCoordCL* const, RAIterT);
};

template<class T=double>
class Quad5_2DCL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;

  protected:
    typedef Quad5_2DCL<T> self_;

  public:
    Quad5_2DCL(): base_type( value_type(), Quad5_2DDataCL::NumNodesC) {}
    Quad5_2DCL(const value_type& t): base_type( t, Quad5_2DDataCL::NumNodesC) {}
    Quad5_2DCL(const TetraCL&, const BaryCoordCL* const, instat_fun_ptr, double= 0.0);
    Quad5_2DCL(const LocalP2CL<value_type>&, const BaryCoordCL*const );
    template <class _BndData, class _VD>
      Quad5_2DCL(const TetraCL&, const BaryCoordCL* const, const P2EvalCL<T, _BndData, _VD>&);
    template <class PFunT>
      Quad5_2DCL(const TetraCL&, const BaryCoordCL* const, const PFunT&);

DROPS_DEFINE_VALARRAY_DERIVATIVE(Quad5_2DCL, T, base_type)

    inline self_&
    assign(const TetraCL&, const BaryCoordCL* const, instat_fun_ptr , double= 0.0);
    inline self_&
    assign(const LocalP1CL<value_type>&, const BaryCoordCL* const);
    inline self_&
    assign(const LocalP2CL<value_type>&, const BaryCoordCL* const);
    template <class _BndData, class _VD>
        inline self_&
        assign(const TetraCL& s, const BaryCoordCL* const, const P2EvalCL<T, _BndData, _VD>&);
    template <class PFunT>
      inline self_&
      assign(const TetraCL&, const BaryCoordCL* const, const PFunT&);

    // Calculates the barycentric coordinates of the quadrature points
    // of the triangle given by the 1st argument with respect to the
    // tetrahedron and stores them in the 2nd argument.
    static void SetInterface (const BaryCoordCL* const p, BaryCoordCL* NodeInTetra) {
        Quad5_2DDataCL::SetInterface ( p, NodeInTetra);
    }

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const {
      return absdet*(9./80.*(*this)[0] + (155. - std::sqrt( 15.0))/2400.*((*this)[1]+(*this)[2]+(*this)[3]) + (155. + std::sqrt( 15.0))/2400.*((*this)[4]+(*this)[5]+(*this)[6]));
    }
};


class Quad3PosWeightsCL
// contains cubatur on reference-tetra, that is exact up to degree 3, positive,
// and uses only 8 points.
// Do not forget to multiply the result of Quad() by the determinant of the affine trafo
// from the reference tetra to the given tetra.
{
  private:
    static const double _points[8][3];

  public:
    static Uint GetNumPoints() { return 8; }
    static const Point3DCL* GetPoints() { return reinterpret_cast<const Point3DCL*>(_points[0]); }

    static inline double Quad(const TetraCL&, instat_scalar_fun_ptr, double);
    static inline double Quad(const TetraCL&, scalar_fun_ptr);
    static inline SVectorCL<3> Quad(const TetraCL&, instat_vector_fun_ptr, double);
    static inline SVectorCL<3> Quad(const TetraCL&, vector_fun_ptr);

    static inline double Quad(const double*);
    static inline SVectorCL<3> Quad(const SVectorCL<3>*);
    template<class IteratorT, class ValueT>
      static inline void
      Quad(IteratorT, ValueT*const);
};

inline double Quad3PosWeightsCL::Quad(const TetraCL& t, instat_scalar_fun_ptr f, double tt)
{
    const Point3DCL* pts= GetPoints();
    return ( f(GetWorldCoord(t, pts[0]), tt) + f(GetWorldCoord(t, pts[1]), tt)
            +f(GetWorldCoord(t, pts[2]), tt) + f(GetWorldCoord(t, pts[3]), tt) )/240.
          +( f(GetWorldCoord(t, pts[4]), tt) + f(GetWorldCoord(t, pts[5]), tt)
            +f(GetWorldCoord(t, pts[6]), tt) + f(GetWorldCoord(t, pts[7]), tt) )*3./80.;

}
inline double Quad3PosWeightsCL::Quad(const TetraCL& t, scalar_fun_ptr f)
{
    const Point3DCL* pts= GetPoints();
    return ( f(GetWorldCoord(t, pts[0])) + f(GetWorldCoord(t, pts[1]))
            +f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )/240.
          +( f(GetWorldCoord(t, pts[4])) + f(GetWorldCoord(t, pts[5]))
            +f(GetWorldCoord(t, pts[6])) + f(GetWorldCoord(t, pts[7])) )*3./80.;

}
inline SVectorCL<3> Quad3PosWeightsCL::Quad(const TetraCL& t, instat_vector_fun_ptr f, double tt)
{
    const Point3DCL* pts= GetPoints();
    return ( f(GetWorldCoord(t, pts[0]), tt) + f(GetWorldCoord(t, pts[1]), tt)
            +f(GetWorldCoord(t, pts[2]), tt) + f(GetWorldCoord(t, pts[3]), tt) )/240.
          +( f(GetWorldCoord(t, pts[4]), tt) + f(GetWorldCoord(t, pts[5]), tt)
            +f(GetWorldCoord(t, pts[6]), tt) + f(GetWorldCoord(t, pts[7]), tt) )*3./80.;

}
inline SVectorCL<3> Quad3PosWeightsCL::Quad(const TetraCL& t, vector_fun_ptr f)
{
    const Point3DCL* pts= GetPoints();
    return ( f(GetWorldCoord(t, pts[0])) + f(GetWorldCoord(t, pts[1]))
            +f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )/240.
          +( f(GetWorldCoord(t, pts[4])) + f(GetWorldCoord(t, pts[5]))
            +f(GetWorldCoord(t, pts[6])) + f(GetWorldCoord(t, pts[7])) )*3./80.;

}

inline double Quad3PosWeightsCL::Quad(const double* vals)
{
    return (vals[0] + vals[1] + vals[2] + vals[3])/240.
          +(vals[4] + vals[5] + vals[6] + vals[7])*3./80.;
}

inline SVectorCL<3> Quad3PosWeightsCL::Quad(const SVectorCL<3>* vals)
{
    return (vals[0] + vals[1] + vals[2] + vals[3])/240.
          +(vals[4] + vals[5] + vals[6] + vals[7])*3./80.;
}

template<class IteratorT, class ValueT>
inline void Quad3PosWeightsCL::Quad(IteratorT beg, ValueT*const ret)
{
    ValueT tmp0= *beg++; tmp0+= *beg++; tmp0+= *beg++; tmp0+= *beg++;
    tmp0/=240.;
    ValueT tmp1= *beg++; tmp1+= *beg++; tmp1+= *beg++; tmp1+= *beg;
    tmp1*=3./240.;
    *ret= tmp0 + tmp1;
}


class FaceQuad2CL
// contains cubatur on reference-face, that is exact up to degree 2, positive,
// and uses 4 points.
{
  private:
    static const double _points[4][2];

  public:
    static Uint GetNumPoints() { return 4; }
    static const Point2DCL* GetPoints() { return reinterpret_cast<const Point2DCL*>(_points[0]); }

    static inline double Quad(const double*);
    static inline SVectorCL<3> Quad(const SVectorCL<3>*);
    template<class IteratorT, class ValueT>
      static inline void
      Quad(IteratorT, ValueT*const);
};

inline double FaceQuad2CL::Quad(const double* vals)
{
    return (vals[0] + vals[1] + vals[2])/24. + vals[3]*3./8.;
}

inline SVectorCL<3> FaceQuad2CL::Quad(const SVectorCL<3>* vals)
{
    return (vals[0] + vals[1] + vals[2])/24.  +vals[3]*3./8.;
}

template<class IteratorT, class ValueT>
inline void FaceQuad2CL::Quad(IteratorT beg, ValueT*const ret)
{
    ValueT tmp= *beg++; tmp+= *beg++; tmp+= *beg++;
    tmp/=24.;
    *ret= tmp + (*beg)*3./8.;
}


//=========================================
//    Finite Elements: P1, P1Bubble, P2
//=========================================

class P1DiscCL
// contains cubatur etc. for linear FE
{
  private:
    template<class SVecT, Uint D>
    static inline SVecT lin_comb (double a0, const SVecT& m0, double a1, const SVecT& m1, double a2, const SVecT& m2, double a3, const SVecT& m3, double a4, const SVecT& m4);

  public:
    // cubatur formula for int f(x) dx, exact up to degree 2
    static inline double Quad(const TetraCL&, scalar_fun_ptr);
    static inline double Quad(const TetraCL&, instat_scalar_fun_ptr, double= 0.0);
    static inline SVectorCL<3> Quad(const TetraCL&, instat_vector_fun_ptr, double= 0.0);
    static inline SMatrixCL<3,3> Quad( const LocalP2CL< SMatrixCL<3,3> >& f, BaryCoordCL** bp);
    static inline SVectorCL<3> Quad( const LocalP2CL< SVectorCL<3> >& f, BaryCoordCL** bp);
    template<class ValueT>
    static inline ValueT Quad( const LocalP2CL<ValueT>&, BaryCoordCL**);
    template<class ValueT>
    static inline ValueT Quad( const LocalP2CL<ValueT>&, const BaryCoordCL[4]);
    template<class ValueT>
    static inline ValueT Quad( const LocalP2CL<ValueT>& f);
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 1
    static inline double Quad(const TetraCL&, scalar_fun_ptr, Uint);
    static inline double Quad(const TetraCL&, instat_scalar_fun_ptr, Uint, double= 0.0);
    static inline SVectorCL<3> Quad(const TetraCL&, vector_fun_ptr, Uint);
    // cubatur formula for int f(x)*phi_i*phi_j dx, exact up to degree 1
    static inline double Quad(const TetraCL&, scalar_fun_ptr, Uint, Uint);
    static inline double Quad(const TetraCL&, instat_scalar_fun_ptr, Uint, Uint, double= 0.0);
    // cubatur formula for int f(x)*phi_i over face, exact up to degree 1
    static inline double Quad2D(const TetraCL&, Uint face, instat_scalar_fun_ptr, Uint, double= 0.0);
    static inline SVectorCL<3> Quad2D(const TetraCL&, Uint face, instat_vector_fun_ptr, Uint, double= 0.0);
    // computes the square of the L2-norm of a given function f:
    // f^2 is integrated exact up to degree 2
    static inline double norm_L2_sq(const TetraCL&, scalar_fun_ptr);
    static inline double norm_L2_sq(const TetraCL&, instat_scalar_fun_ptr, double= 0.0);
    // returns int phi_i*phi_j dx
    static inline double GetMass( int i, int j) { return i!=j ? 1./120. : 1./60.; }
    // returns int phi_i dx
    static inline double GetLumpedMass( int) { return 1./24.; }
    // the gradient of hat function i is in column i of H
    static inline void   GetGradients( SMatrixCL<3,4>& H, double& det, const TetraCL& t);
    static inline void   GetGradients( Point3DCL H[4],    double& det, const TetraCL& t);
    static inline void   GetGradients( SMatrixCL<3,4>& H, double& det, const Point3DCL pt[4]);
    static inline void   GetGradients( Point3DCL H[4], const SMatrixCL<3,3>& T);
    static void GetP1Basis( Quad5_2DCL<> p1[4], const BaryCoordCL* const p);
};

class P1DDiscCL
// contains cubatur etc. for linear, discontinuous FE
{
  public:
    // the gradient of hat function i is in column i of H
    static inline void   GetGradients( SMatrixCL<3,4>& H, double& det, const TetraCL& t);
};

class P1BubbleDiscCL
// contains cubatur etc. for linear FE with bubble function for the barycenter
{
  private:
    static const double _points1[26][3];

  public:
    static Uint GetNumPoints1(Uint i) { return i<4 ? 4 : 10; }
    static const Point3DCL* GetPoints1(Uint i)
        { return reinterpret_cast<const Point3DCL*>(_points1[i*4]); }

    // cubatur formula for int f(x)*phi_i dx, exact up to degree 2
    // the last two forms take an array that contains the values of
    // the function to be evaluated in the points obtained by GetPoints1(i)
    static inline double Quad(const TetraCL&, scalar_fun_ptr, Uint);
    static inline SVectorCL<3> Quad(const TetraCL&, vector_fun_ptr, Uint);
    static inline double Quad(const double*, Uint);
    static inline SVectorCL<3> Quad(const SVectorCL<3>*, Uint);
};

class P2DiscCL
{
  public:
    // gradients on reference tetra
    static void GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[10]);
    static void GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10]);
    static void GetGradientsOnRef( Quad5CL<Point3DCL> GRef[10]);
    // The 2nd arg points to 3 vertices of the triangle
    static void GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[10], const BaryCoordCL* const);
    // p2[i] contains a Quad5_2DCL-object that is initialized with FE_P2CL::Hi
    static void GetP2Basis( Quad5_2DCL<> p2[10], const BaryCoordCL* const p);
    // compute gradients
    static void GetGradients( LocalP1CL<Point3DCL> G[10], LocalP1CL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<4; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad2CL<Point3DCL> G[10], Quad2CL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<5; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad5CL<Point3DCL> G[10], Quad5CL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<Quad5DataCL::NumNodesC; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradients( Quad5_2DCL<Point3DCL> G[10], Quad5_2DCL<Point3DCL> GRef[10], const SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<Quad5_2DDataCL::NumNodesC; ++j) G[i][j]= T*GRef[i][j]; }
    static void GetGradient( Quad2CL<Point3DCL> &G, Quad2CL<Point3DCL> &GRef, const SMatrixCL<3,3> &T)
    { for (int j=0; j<5; ++j) G[j]= T*GRef[j]; }
    static void GetGradient( Quad5CL<Point3DCL> &G, Quad5CL<Point3DCL> &GRef, const SMatrixCL<3,3> &T)
    { for (int j=0; j<Quad5DataCL::NumNodesC; ++j) G[j]= T*GRef[j]; }
    /// compute gradient of a function provided as LocalP2CL<double> object
    template<class GradT>
    static void GetFuncGradient( GradT& gradF, const LocalP2CL<>& F, const GradT G[10])
    { gradF= F[0]*G[0]; for (int i=1; i<10; ++i) gradF+= F[i]*G[i]; }
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 1
    static inline SVectorCL<3> Quad( const TetraCL& tetra, instat_vector_fun_ptr, Uint, double= 0.0);
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 2
    template<class valT>
    static inline valT Quad( valT f[10], int i);
    // returns int phi_i phi_j dx
    static inline double GetMass( int i, int j);
    // returns int phi_i dx
    static inline double GetLumpedMass( int i) { return i<4 ? -1./120. : 1./30.; }
};

class P2RidgeDiscCL
/// \brief contains helper functions for the XFEM discretization based on ridge enrichment.
///
/// The P2 ridge XFEM space is the direct sum of the  P2 FE space and P1*F_R (extended part)
/// with the ridge function \f$ F_R = \sum |\phi_i| v_i - |\phi|\f$
/// as defined in Moes et al., Comput. Methods Appl. Mech Engrg. 192 (2003), pp. 3163-3177,
/// where \f$\phi=\sum \phi_i v_i\f$ is the level set function,
{
  public:
    static void GetEnrichmentFunction( LocalP2CL<>& ridgeFunc_p,    LocalP2CL<>& ridgeFunc_n,    const LocalP2CL<>& lset);
    static void GetExtBasisOnChildren( LocalP2CL<> p1ridge_p[4][8], LocalP2CL<> p1ridge_n[4][8], const LocalP2CL<>& lset);
    static void GetExtBasisPointwise ( LocalP2CL<> p1ridge_p[4],    LocalP2CL<> p1ridge_n[4],    const LocalP2CL<>& lset);
};

/// splits a p2r-VectorCL into two p2-VectorCL (should work for both P2R_FE and vecP2R_FE)
void P2RtoP2( const IdxDescCL& p2ridx, const VectorCL& p2r, const IdxDescCL& p2idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const BndDataCL<>& lsetbnd, const MultiGridCL& mg);


inline double FuncDet2D( const Point3DCL& p, const Point3DCL& q)
{
    const double d0= p[1]*q[2] - p[2]*q[1];
    const double d1= p[2]*q[0] - p[0]*q[2];
    const double d2= p[0]*q[1] - p[1]*q[0];
    return std::sqrt(d0*d0 + d1*d1 + d2*d2);
}


/********************************************************************************
*
*        definition of   i n l i n e   f u n c t i o n s
*
********************************************************************************/

inline double VolFrac(BaryCoordCL** bp)
{
    double M[3][3];
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            (M[i][j]= (*bp[j+1])[i+1]-(*bp[0])[i+1]);
    return std::abs( M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
                   - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
                   + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]) );
}

inline double VolFrac(const SArrayCL<BaryCoordCL,4>& bp)
{
    double M[3][3];
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            (M[i][j]= (bp[j+1])[i+1]-(bp[0])[i+1]);
    return std::abs( M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
            - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
            + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]) );
}

inline double VolFrac (const SMatrixCL<4, 4>& A)
{
    double M[3][3];
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            M[i][j]= A(i+1, j+1) - A(i+1, 0);
    return std::fabs( M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
           - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
           + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]));
}

inline double P1DiscCL::Quad(const TetraCL& s, instat_scalar_fun_ptr coeff, double t)
{
    return ( coeff(s.GetVertex(0)->GetCoord(), t)
            +coeff(s.GetVertex(1)->GetCoord(), t)
            +coeff(s.GetVertex(2)->GetCoord(), t)
            +coeff(s.GetVertex(3)->GetCoord(), t))/120.
            + 2./15.*coeff(GetBaryCenter(s), t);
}

inline double P1DiscCL::Quad(const TetraCL& t, scalar_fun_ptr coeff)
{
    return ( coeff(t.GetVertex(0)->GetCoord())
            +coeff(t.GetVertex(1)->GetCoord())
            +coeff(t.GetVertex(2)->GetCoord())
            +coeff(t.GetVertex(3)->GetCoord()))/120.
            + 2./15.*coeff(GetBaryCenter(t));
}

inline SVectorCL<3> P1DiscCL::Quad(const TetraCL& s, instat_vector_fun_ptr coeff, double t)
{
    return ( coeff(s.GetVertex(0)->GetCoord(), t)
            +coeff(s.GetVertex(1)->GetCoord(), t)
            +coeff(s.GetVertex(2)->GetCoord(), t)
            +coeff(s.GetVertex(3)->GetCoord(), t))/120.
            + 2./15.*coeff(GetBaryCenter(s), t);
}

inline double P1DiscCL::Quad( const TetraCL& s, instat_scalar_fun_ptr coeff, Uint i, double t)
{
    double f_Vert_i= coeff( s.GetVertex(i)->GetCoord(), t ),
           f_Bary  = coeff( GetBaryCenter(s), t ),
           f_Other = 0;

    for (Uint k=0; k<4; ++k)
        if (k!=i) f_Other+= coeff( s.GetVertex(k)->GetCoord(), t );
    return f_Vert_i/108. + f_Other/1080. + 4./135.*f_Bary;
}

template<class SVecT, Uint D>
inline SVecT
P1DiscCL::lin_comb (double a0, const SVecT& m0, double a1, const SVecT& m1, double a2, const SVecT& m2, double a3, const SVecT& m3, double a4, const SVecT& m4)
{
    SVecT ret( Uninitialized);
    for (Uint i= 0; i < D; ++i)
        ret[i]= a0*m0[i] + a1*m1[i] + a2*m2[i] + a3*m3[i] +a4*m4[i];
    return ret;
}

inline SMatrixCL<3,3> P1DiscCL::Quad( const LocalP2CL< SMatrixCL<3,3> >& f, BaryCoordCL** bp)
{
    const BaryCoordCL bc= 0.25*(*bp[0]+*bp[1]+*bp[2]+*bp[3]);
    return P1DiscCL::lin_comb<SMatrixCL<3,3>, 9>( 1./120, f(*bp[0]), 1./120, f(*bp[1]), 1./120, f(*bp[2]), 1./120, f(*bp[3]), 2./15, f(bc));
}

inline SVectorCL<3> P1DiscCL::Quad( const LocalP2CL< SVectorCL<3> >& f, BaryCoordCL** bp)
{
    const BaryCoordCL bc= 0.25*(*bp[0]+*bp[1]+*bp[2]+*bp[3]);
    return P1DiscCL::lin_comb<SVectorCL<3>, 3>( 1./120, f(*bp[0]), 1./120, f(*bp[1]), 1./120, f(*bp[2]), 1./120, f(*bp[3]), 2./15, f(bc));
}

template<class ValueT>
inline ValueT P1DiscCL::Quad( const LocalP2CL<ValueT>& f, BaryCoordCL** bp)
{
    const BaryCoordCL bc= 0.25*(*bp[0]+*bp[1]+*bp[2]+*bp[3]);
    return ( f(*bp[0]) + f(*bp[1]) + f(*bp[2]) + f(*bp[3]) )/120. + 2./15. * f(bc);
}

template<class ValueT>
inline ValueT P1DiscCL::Quad( const LocalP2CL<ValueT>& f, const BaryCoordCL bp[4])
{
    const BaryCoordCL bc= 0.25*(bp[0]+bp[1]+bp[2]+bp[3]);
    return ( f(bp[0]) + f(bp[1]) + f(bp[2]) + f(bp[3]) )/120. + 2./15. * f(bc);
}

template<class ValueT>
inline ValueT P1DiscCL::Quad( const LocalP2CL<ValueT>& f)
{
    return 1./30.*(f[4] + f[5] + f[6] + f[7] + f[8] + f[9]) - 1./120.*(f[0] + f[1] + f[2] + f[3]);
}


inline double P1DiscCL::Quad( const TetraCL& t, scalar_fun_ptr coeff, Uint i)
{
    double f_Vert_i= coeff( t.GetVertex(i)->GetCoord() ),
           f_Bary  = coeff( GetBaryCenter(t) ),
           f_Other = 0;

    for (Uint k=0; k<4; ++k)
        if (k!=i) f_Other+= coeff( t.GetVertex(k)->GetCoord() );
    return f_Vert_i/108. + f_Other/1080. + 4./135.*f_Bary;
}
inline SVectorCL<3> P1DiscCL::Quad( const TetraCL& t, vector_fun_ptr coeff, Uint i)
{
    SVectorCL<3> f_Vert_i= coeff( t.GetVertex(i)->GetCoord() ),
                 f_Bary  = coeff( GetBaryCenter(t) ),
                 f_Other(0.0);;

    for (Uint k=0; k<4; ++k)
        if (k!=i) f_Other+= coeff( t.GetVertex(k)->GetCoord() );
    return f_Vert_i/108. + f_Other/1080. + 4./135.*f_Bary;
}

inline double P1DiscCL::Quad( const TetraCL& s, instat_scalar_fun_ptr coeff, Uint i, Uint j, double t)
{
    double f_Vert_ij= coeff( s.GetVertex(i)->GetCoord(), t ),
           f_Bary  = coeff( GetBaryCenter(s), t ),
           f_Other = 0;

    if (i==j)
    {
        for (Uint k=0; k<4; ++k)
            if (k!=i) f_Other+= coeff( s.GetVertex(k)->GetCoord(), t );
        return 43./7560.*f_Vert_ij + f_Other/7560. + 2./189.*f_Bary;
    }
    else
    {
        f_Vert_ij+= coeff( s.GetVertex(j)->GetCoord(), t );
        for (Uint k=0; k<4; ++k)
            if (k!=i && k!=j) f_Other+= coeff( s.GetVertex(k)->GetCoord(), t );
        return 11./7560.*f_Vert_ij + f_Other/15120. + f_Bary/189.;
    }
}

inline double P1DiscCL::Quad( const TetraCL& t, scalar_fun_ptr coeff, Uint i, Uint j)
{
    double f_Vert_ij= coeff( t.GetVertex(i)->GetCoord() ),
           f_Bary  = coeff( GetBaryCenter(t) ),
           f_Other = 0;

    if (i==j)
    {
        for (Uint k=0; k<4; ++k)
            if (k!=i) f_Other+= coeff( t.GetVertex(k)->GetCoord() );
        return 43./7560.*f_Vert_ij + f_Other/7560. + 2./189.*f_Bary;
    }
    else
    {
        f_Vert_ij+= coeff( t.GetVertex(j)->GetCoord() );
        for (Uint k=0; k<4; ++k)
            if (k!=i && k!=j) f_Other+= coeff( t.GetVertex(k)->GetCoord() );
        return 11./7560.*f_Vert_ij + f_Other/15120. + f_Bary/189.;
    }
}

inline double P1DiscCL::Quad2D(const TetraCL& t, Uint face, instat_scalar_fun_ptr bfun, Uint vert, double time)
// Integrate neu_val() * phi_vert over face
{
    Point3DCL v[3];

    if (face==OppFace(vert)) return 0.; // vert is not on face

    v[0]= t.GetVertex(vert)->GetCoord();
    for (Uint i=0, k=1; i<3; ++i)
    {
        if (VertOfFace(face,i)!=vert)
            v[k++]= t.GetVertex( VertOfFace(face,i) )->GetCoord();
    }
    const double f0= bfun(v[0], time);
    const double f1= bfun(v[1], time) +  bfun( v[2], time);
    const double f2= bfun(1./3.*(v[0] + v[1] + v[2]), time);    //Barycenter of Face
    const double absdet= FuncDet2D(v[1] - v[0], v[2] - v[0]);
    return (11./240.*f0 + 1./240.*f1 + 9./80.*f2) * absdet;
}

inline SVectorCL<3> P1DiscCL::Quad2D(const TetraCL& t, Uint face, instat_vector_fun_ptr bfun, Uint vert, double time)
// Integrate neu_val() * phi_vert over face
{
    Point3DCL v[3];

    if (face==OppFace(vert)) return SVectorCL<3>(0.); // vert is not on face

    v[0]= t.GetVertex(vert)->GetCoord();
    for (Uint i=0, k=1; i<3; ++i)
    {
        if (VertOfFace(face,i)!=vert)
            v[k++]= t.GetVertex( VertOfFace(face,i) )->GetCoord();
    }
    const SVectorCL<3> f0= bfun(v[0], time);
    const SVectorCL<3> f1= bfun(v[1], time) +  bfun( v[2], time);
    const SVectorCL<3> f2= bfun(1./3.*(v[0] + v[1] + v[2]), time);    //Barycenter of Face
    const double absdet= FuncDet2D(v[1] - v[0], v[2] - v[0]);
    return (11./240.*f0 + 1./240.*f1 + 9./80.*f2) * absdet;
}

inline double P1DiscCL::norm_L2_sq(const TetraCL& s, instat_scalar_fun_ptr coeff, double t)
{
    const double f0= coeff(s.GetVertex(0)->GetCoord(), t);
    const double f1= coeff(s.GetVertex(1)->GetCoord(), t);
    const double f2= coeff(s.GetVertex(2)->GetCoord(), t);
    const double f3= coeff(s.GetVertex(3)->GetCoord(), t);
    const double fb= coeff(GetBaryCenter(s), t);
    return (f0*f0 + f1*f1 + f2*f2 + f3*f3)/120. + 2./15.*fb*fb;

}

inline double P1DiscCL::norm_L2_sq(const TetraCL& t, scalar_fun_ptr coeff)
{
    const double f0= coeff(t.GetVertex(0)->GetCoord());
    const double f1= coeff(t.GetVertex(1)->GetCoord());
    const double f2= coeff(t.GetVertex(2)->GetCoord());
    const double f3= coeff(t.GetVertex(3)->GetCoord());
    const double fb= coeff(GetBaryCenter(t));
    return (f0*f0 + f1*f1 + f2*f2 + f3*f3)/120. + 2./15.*fb*fb;

}


inline void P1DiscCL::GetGradients( SMatrixCL<3,4>& H, double& det, const TetraCL& t)
{
    double M[3][3];
    const Point3DCL& pt0= t.GetVertex(0)->GetCoord();
    for(Uint i=0; i<3; ++i)
        for(Uint j=0; j<3; ++j)
            M[j][i]= t.GetVertex(i+1)->GetCoord()[j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    H(0,1)= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    H(0,2)= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    H(0,3)= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    H(1,1)= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    H(1,2)= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    H(1,3)= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    H(2,1)= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    H(2,2)= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    H(2,3)= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
    // in H(1:3,0:2) steht jetzt die Adjunkte von M ...
    H(0,0)= -H(0,1)-H(0,2)-H(0,3);
    H(1,0)= -H(1,1)-H(1,2)-H(1,3);
    H(2,0)= -H(2,1)-H(2,2)-H(2,3);
}

inline void P1DiscCL::GetGradients( Point3DCL H[4], double& det, const TetraCL& t)
{
    double M[3][3];
    const Point3DCL& pt0= t.GetVertex(0)->GetCoord();
    for(Uint i=0; i<3; ++i)
        for(Uint j=0; j<3; ++j)
            M[j][i]= t.GetVertex(i+1)->GetCoord()[j] - pt0[j];
    det=   M[0][0] * (M[1][1]*M[2][2] - M[1][2]*M[2][1])
         - M[0][1] * (M[1][0]*M[2][2] - M[1][2]*M[2][0])
         + M[0][2] * (M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    H[1][0]= (M[1][1]*M[2][2] - M[1][2]*M[2][1])/det;
    H[1][1]= (M[2][1]*M[0][2] - M[0][1]*M[2][2])/det;
    H[1][2]= (M[0][1]*M[1][2] - M[1][1]*M[0][2])/det;
    H[2][0]= (M[2][0]*M[1][2] - M[1][0]*M[2][2])/det;
    H[2][1]= (M[0][0]*M[2][2] - M[2][0]*M[0][2])/det;
    H[2][2]= (M[1][0]*M[0][2] - M[0][0]*M[1][2])/det;
    H[3][0]= (M[1][0]*M[2][1] - M[2][0]*M[1][1])/det;
    H[3][1]= (M[2][0]*M[0][1] - M[0][0]*M[2][1])/det;
    H[3][2]= (M[0][0]*M[1][1] - M[1][0]*M[0][1])/det;
    // in H[1:3][0:2] steht jetzt die Adjunkte von M ...
    H[0]= -H[1]-H[2]-H[3];
}

inline void P1DiscCL::GetGradients( Point3DCL H[4], const SMatrixCL<3,3>& T)
{
    for(Uint i=0; i<4; ++i)
        H[i]= T*FE_P1CL::DHRef( i);
}

inline void P1DiscCL::GetGradients (SMatrixCL<3,4>& H, double& det, const Point3DCL pt[4])
{
    SMatrixCL<3 ,3> M;
    GetTrafoTr( M, det, pt);
    for (Uint i= 0; i<4; ++i) {
        const Point3DCL tmp( M*FE_P1CL::DHRef( i));
        H(0, i)= tmp[0];
        H(1, i)= tmp[1];
        H(2, i)= tmp[2];
    }
}

inline void P1DDiscCL::GetGradients( SMatrixCL<3,4>& H, double& det, const TetraCL& t)
{
    SMatrixCL<3 ,3> M;
    GetTrafoTr( M, det, t);
    for (Uint i= 0; i<4; ++i) {
        const Point3DCL tmp( M*FE_P1DCL::DHRef( i));
        H(0, i)= tmp[0];
        H(1, i)= tmp[1];
        H(2, i)= tmp[2];
    }
}


inline double P1BubbleDiscCL::Quad(const TetraCL& t, scalar_fun_ptr f, Uint i)
{
    const Point3DCL* pts= GetPoints1(i);
    if (i<4)
        return   f(GetWorldCoord(t, pts[0]))/240.
              +( f(GetWorldCoord(t, pts[1])) + f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )/80.;
    return -4./945.*( f(GetWorldCoord(t, pts[0])) + f(GetWorldCoord(t, pts[1])) + f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )
           +32./2835.*( f(GetWorldCoord(t, pts[4])) + f(GetWorldCoord(t, pts[5])) + f(GetWorldCoord(t, pts[6]))
                       +f(GetWorldCoord(t, pts[7])) + f(GetWorldCoord(t, pts[8])) + f(GetWorldCoord(t, pts[9])) );
}

inline SVectorCL<3> P1BubbleDiscCL::Quad(const TetraCL& t, vector_fun_ptr f, Uint i)
{
    const Point3DCL* pts= GetPoints1(i);
    if (i<4)
        return   f(GetWorldCoord(t, pts[0]))/240.
              +( f(GetWorldCoord(t, pts[1])) + f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )/80.;
    return -4./945.*( f(GetWorldCoord(t, pts[0])) + f(GetWorldCoord(t, pts[1])) + f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )
           +32./2835.*( f(GetWorldCoord(t, pts[4])) + f(GetWorldCoord(t, pts[5])) + f(GetWorldCoord(t, pts[6]))
                       +f(GetWorldCoord(t, pts[7])) + f(GetWorldCoord(t, pts[8])) + f(GetWorldCoord(t, pts[9])) );
}

inline double P1BubbleDiscCL::Quad(const double* vals, Uint i)
{
    if (i<4)
        return vals[0]/240. + (vals[1] + vals[2] + vals[3])/80.;
    return -4./945.*( vals[0] + vals[1] + vals[2] + vals[3] )
           +32./2835.*(vals[4] + vals[5] + vals[6] + vals[7] + vals[8] + vals[9]);
}

inline SVectorCL<3> P1BubbleDiscCL::Quad(const SVectorCL<3>* vals, Uint i)
{
    if (i<4)
        return vals[0]/240. + (vals[1] + vals[2] + vals[3])/80.;
    return -4./945.*( vals[0] + vals[1] + vals[2] + vals[3] )
           +32./2835.*(vals[4] + vals[5] + vals[6] + vals[7] + vals[8] + vals[9]);
}

inline SVectorCL<3> P2DiscCL::Quad( const TetraCL& tetra, instat_vector_fun_ptr coeff, Uint i, double t)
{
    SVectorCL<3> f[5];

    if (i<4) // hat function on vert
    {
        f[0]= coeff( tetra.GetVertex(i)->GetCoord(), t);
        for (Uint k=0, l=1; k<4; ++k)
            if (k!=i) f[l++]= coeff( tetra.GetVertex(k)->GetCoord(), t);
        f[4]= coeff( GetBaryCenter(tetra), t);
        return f[0]/504. - (f[1] + f[2] + f[3])/1260. - f[4]/126.;
    }
    else  // hat function on edge
    {
        const double ve= 4./945.,  // coeff for verts of edge
                     vn= -1./756.,  // coeff for other verts
                     vs= 26./945.;   // coeff for barycenter
        double a[4];
        a[VertOfEdge(i-4,0)]= a[VertOfEdge(i-4,1)]= ve;
        a[VertOfEdge(OppEdge(i-4),0)]= a[VertOfEdge(OppEdge(i-4),1)]= vn;

        SVectorCL<3> sum= vs * coeff( GetBaryCenter(tetra), t);
        for(Uint k=0; k<4; ++k)
            sum+= a[k] * coeff( tetra.GetVertex(k)->GetCoord(), t);

        return sum;
    }
}

template<class valT>
inline valT P2DiscCL::Quad( valT f[10], int i)
{
    valT sum= valT(), result;
    if (i<4) // hat function on vert
    {
        // Q = sum c[i]*f[i]
        // Gewichte c[i] = 1/420    fuer vert i
        //                 1/2520   fuer uebrige verts
        //                -1/630    fuer an vert i anliegende edges
        //                -1/420    fuer uebrige drei edges
        result= f[i]*(1/420.-1./2520.);
        for (int k=0; k<4; ++k)
            sum+= f[k];
        result+= sum/2520.;
        sum= valT();
        for (int k=0; k<3; ++k)
            sum+= f[EdgeByVert(i,k)+4];
        result+= -sum/630.;
        sum= valT();
        const int oppF=OppFace(i);
        for (int k=0; k<3; ++k)
            sum+= f[EdgeOfFace(oppF, k)+4];
        result+= -sum/420.;
        return result;
    }
    else  // hat function on edge
    {
        i-= 4;
        // Q = sum c[i]*f[i]
        // Gewichte c[i] = 4/315    fuer egde i
        //                 1/315    fuer opposite edge
        //                 2/315    fuer uebrige edges
        //                -1/630    fuer an edge i anliegende verts
        //                -1/420    fuer uebrige zwei verts
        result=  f[i+4]*4./315.;
        const int opp= OppEdge(i);
        result+= f[opp+4]/315.;
        for(int k=0; k<6; ++k)
            if (k!=i && k!=opp)
                sum+= f[k+4];
        result+= sum*2./315.;
        sum= f[VertOfEdge(i,0)] + f[VertOfEdge(i,1)];
        result+= -sum/630.;
        sum= f[VertOfEdge(OppEdge(i),0)] + f[VertOfEdge(OppEdge(i),1)];
        result+= -sum/420.;

        return result;
    }
}


inline double P2DiscCL::GetMass( int i, int j)
{
// zur Erlaeuterung der zurueckgegebenen Gewichte
// beachte Kommentare zu Funktion P2DiscCL::Quad obendrueber!

    if (i>j) { int h=i; i=j; j=h; } // swap such that i<=j holds
    if (j<4) // i,j are vertex-indices
    {
        if (i==j)
            return 1./420.;
        else
            return 1./2520.;
    }
    else if (i>=4) // i,j are edge-indices
    {
        if (i==j)
            return 4./315.;
        else
        {
            if (i-4==OppEdge(j-4))
                return 1./315.;
            else
                return 2./315.;
        }
    }
    else // i<4, j>=4: Sonderfall...
    {
        if (i==VertOfEdge(j-4,0) || i==VertOfEdge(j-4,1))
            return -1./630.;
        else
            return -1./420.;
    }
}

} // end of namespace DROPS

#include "num/discretize.tpp"

#endif
