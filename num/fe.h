//**************************************************************************
// File:    fe.h                                                           *
// Content: description of various finite-element functions                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - May, 2 2001                                            *
//**************************************************************************
// TODO: Faces for P1/2EvalCL!!!

#ifndef _FE_H_
#define _FE_H_


#include "misc/container.h"
#include "geom/multigrid.h"
#include "misc/problem.h"


namespace DROPS
{

//**************************************************************************
// Class:   FE_P1CL                                                        *
// Purpose: Shape functions and their gradients for piecewise linear,      *
//          continuous finite elements on the reference tetrahedron        *
//          The number of the H-functions refer to the number of the       *
//          vertex in the tetrahedron as defined in topo.h, where the      *
//          degree of freedom sits.                                        *
//**************************************************************************
class FE_P1CL
{
private:
    static const double _gradient[4][3];
    
public:
    // default ctor, copy-ctor, assignment-op, dtor

    static bool IsC0() { return true; }
    static bool IsC1() { return false; }

    // restriction of the shape functions to reference edge
    static double H0(double v1) { return 1. -v1; }
    static double H1(double v1) { return v1; }

    // restriction of the shape functions to reference face
    static double H0(double v1, double v2) { return 1. -v1 -v2; }
    static double H1(double v1, double)    { return v1; }
    static double H2(double, double v2)    { return v2; }

    // shape functions on the reference tetrahedron
    static double H0(double v1, double v2, double v3) { return 1. -v1 -v2 -v3; }
    static double H1(double v1, double, double)       { return v1; }
    static double H2(double, double v2, double)       { return v2; }
    static double H3(double, double, double v3)       { return v3; }

    // gradients of the shape functions on the reference tetrahedron.
    // To obtain the gradient on tetra T: Let f:x->A*x+b be the affine
    // transformation that carries the reference tetrahefron to T.
    // Then DHi(y) == inverse(transpose(A))*DHiRef(x), where y=f(x).
    static SVectorCL<3> DH0Ref() { return _gradient[0]; }
    static SVectorCL<3> DH1Ref() { return _gradient[1]; }
    static SVectorCL<3> DH2Ref() { return _gradient[2]; }
    static SVectorCL<3> DH3Ref() { return _gradient[3]; }

    // DHRef(i, ...) == DHiRef(...)
    static SVectorCL<3> DHRef(Uint i) { return _gradient[i]; }
};


//**************************************************************************
// Class:   FE_P2CL                                                        *
// Purpose: Shape functions and their gradients for piecewise quadratic,   *
//          continuous finite elements on the reference tetrahedron        *
//          The number of the H-functions refer to the number of the       *
//          (mid-) vertex in the tetrahedron as defined in topo.h, where   *
//          the degree of freedom lives.                                   *
//**************************************************************************
class FE_P2CL
{
  private:
    static const double _D2H[10][3][3];

  public:
    // default ctor, copy-ctor, assignment-op, dtor

    static bool IsC0() { return true; }
    static bool IsC1() { return false; }

    // restriction of the shape functions to reference edge
    static double H0(double v1) { return 1. +v1*(2.*v1 -3.); }
    static double H1(double v1) { return v1*(2.*v1 -1.); }
    static double H2(double v1) { return 4.*v1*(1. -v1); }

    // restriction of the shape functions to reference face
    static double H0(double v1, double v2) { const double sum= v1 + v2; return 1. +sum*(2.*sum -3.); }
    static double H1(double v1, double)    { return v1*(2.*v1 -1.); }
    static double H2(double, double v2)    { return v2*(2.*v2 -1.); }
    static double H3(double v1, double v2) { return 4.*v1*( 1. -(v1 + v2) ); }
    static double H4(double v1, double v2) { return 4.*v2*( 1. -(v1 + v2) ); }
    static double H5(double v1, double v2) { return 4.*v1*v2; }

    // restriction of the shape functions to reference tetrahedron
    static double H0(double v1, double v2, double v3) { const double sum= v1 + v2 + v3; return 1. +sum*(2.*sum -3.); }
    static double H1(double v1, double, double)       { return v1*(2.*v1 -1.); }
    static double H2(double, double v2, double)       { return v2*(2.*v2 -1.); }
    static double H3(double, double, double v3)       { return v3*(2.*v3 -1.); }
    static double H4(double v1, double v2, double v3) { return 4.*v1*( 1. -(v1 + v2 + v3) ); }
    static double H5(double v1, double v2, double v3) { return 4.*v2*( 1. -(v1 + v2 + v3) ); }
    static double H6(double v1, double v2, double)    { return 4.*v1*v2; }
    static double H7(double v1, double v2, double v3) { return 4.*v3*( 1. -(v1 + v2 + v3) ); }
    static double H8(double v1, double, double v3)    { return 4.*v1*v3; }
    static double H9(double, double v2, double v3)    { return 4.*v2*v3; }

    // gradients of the shape functions on the reference tetrahedron.
    // To obtain the gradient on tetra T: See comments in FE_P1CL.
    static inline SVectorCL<3> DH0Ref(double, double, double);
    static inline SVectorCL<3> DH1Ref(double, double, double);
    static inline SVectorCL<3> DH2Ref(double, double, double);
    static inline SVectorCL<3> DH3Ref(double, double, double);
    static inline SVectorCL<3> DH4Ref(double, double, double);
    static inline SVectorCL<3> DH5Ref(double, double, double);
    static inline SVectorCL<3> DH6Ref(double, double, double);
    static inline SVectorCL<3> DH7Ref(double, double, double);
    static inline SVectorCL<3> DH8Ref(double, double, double);
    static inline SVectorCL<3> DH9Ref(double, double, double);

    // DHREef(i, ...) == DHiRef(...)
    static inline SVectorCL<3> DHRef(Uint dof, double v1, double v2, double v3);

    // D2HRef(i) == second derivative of H_i; this is constant
    // diff( H_d(x), k,i) = sum( M_ij*M_kl*D2HRef(d, j, l), j,l=0..2)
    static inline double D2HRef(Uint dof, Uint r, Uint s)
        { return _D2H[dof][r][s]; }

    // Laplace(H_d) on a tetrahedron T; M:= transpose(inverse(A)), where A is the matrix
    // of the affine transformation that maps the reference tetrahedron onto T.
    static inline double Laplace(Uint dof, const SMatrixCL<3,3>& M);
};


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
// Class:   FE_P1BubbleCL                                                  *
// Purpose: Shape functions and their gradients for the mini element.      *
//          The number of the H-functions refer to the number of the       *
//          vertex in the tetrahedron as defined in topo.h, where the      *
//          degree of freedom sits. The last shape function is the bubble  *
//          function on the tetrahedron.                                   *
//**************************************************************************
class FE_P1BubbleCL
{
private:
    static const FE_P1CL _LinFun;
    
public:
    // default ctor, copy-ctor, assignment-op, dtor

    static bool IsC0() { return true; }
    static bool IsC1() { return false; }

    // restriction of the shape functions to reference edge
    static double H0(double v1) { return 1. -v1; }
    static double H1(double v1) { return v1; }

    // restriction of the shape functions to reference face
    static double H0(double v1, double v2) { return 1. -v1 -v2; }
    static double H1(double v1, double)    { return v1; }
    static double H2(double, double v2)    { return v2; }

    // shape functions on the reference tetrahedron
    static double H0(double v1, double v2, double v3) { return 1. -v1 -v2 -v3; }
    static double H1(double v1, double, double)       { return v1; }
    static double H2(double, double v2, double)       { return v2; }
    static double H3(double, double, double v3)       { return v3; }
    static double H4(double v1, double v2, double v3) { return 256.*v1*v2*v3*(1 -v1 -v2 -v3); }

    // gradients of the shape functions on the reference tetrahedron.
    // To obtain the gradient on tetra T: Let f:x->A*x+b be the affine
    // transformation that carries the reference tetrahefron to T.
    // Then DHi(y) == inverse(transpose(A))*DHiRef(x), where y=f(x).
    static SVectorCL<3> DH0Ref() { return _LinFun.DH0Ref(); }
    static SVectorCL<3> DH1Ref() { return _LinFun.DH1Ref(); }
    static SVectorCL<3> DH2Ref() { return _LinFun.DH2Ref(); }
    static SVectorCL<3> DH3Ref() { return _LinFun.DH3Ref(); }
    static SVectorCL<3> DH4Ref(double v1, double v2, double v3)
        { const double tmp1= 1. -v1 -v2 -v3; SVectorCL<3> ret(-v1*v2*v3);
          ret[0]+= tmp1*v2*v3; ret[1]+= tmp1*v1*v3; ret[2]+= tmp1*v1*v2;
          return 256.*ret; }

    // DHREef(i, ...) == DHiRef(...)
    static SVectorCL<3> DHRef(Uint i) { return _LinFun.DHRef(i); }
    static SVectorCL<3> DHRef(Uint i, double v1, double v2, double v3)
        { return i<4 ? DHRef(i) : DH4Ref(v1, v2, v3); }
};


//**************************************************************************
// Class:   DoFHelperCL                                                    *
// Template Parameter:                                                     *
//          Data - The result-type of the finite-element-function on the   *
//                 multigrid                                               *
//          _Vec - The type of the object, that is used for actually       *
//                 storing numerical data, which is referenced by          *
//                 UnknownIdxCL. It is assumed to provide operator[] to    *
//                 access numerical data.                                  *
// Purpose: This class provides an interface for the PiEvalCL to translate *
//          numerical data in _Vec-objects and Data-objects. It must be    *
//          specialized for each Data - type that is used to represent     *
//          numerical data in PDE-solvers, discretizers, etc., as it       *
//          defines the mapping between the "PDE-data-format" (Data) and   *
//          the format for storing a collection of objects of this type    *
//          (_Vec).                                                        *
//          As there is no default-mapping, the general template will throw*
//          an exception upon usage.                                       *
//**************************************************************************
template<class Data, class _Vec>
struct DoFHelperCL
{
    // read the Data-object referenced by IdxT in _Vec
    static inline Data
    get(const _Vec&, const IdxT*)
    { throw DROPSErrCL("DoFHelperCL::get: Unknown data-type"); return Data(); }
    // write the Data-object at position IdxT to the _Vec-object 
    static inline void
    set(_Vec&, const IdxT*, const Data&)
    { throw DROPSErrCL("DoFHelperCL::set: Unknown data-type"); }
};

//**************************************************************************
// Class:   DoFHelperCL                                                    *
// Template Parameter:                                                     *
//          Data - double, partial specialization of the general template  *                                                  *
//          _Vec - The type of the object, that is used for actually       *
//                 storing numerical data, which is referenced by          *
//                 UnknownIdxCL. It is assumed to provide operator[] to    *
//                 access numerical data.                                  *
// Purpose: reads and writes double-objects from / to _Vec objects         *
//          interpreting IdxT as index in _Vec.                            *
//**************************************************************************
template<class _Vec>
struct DoFHelperCL<double, _Vec>
{
    static inline double
    get(const _Vec& sol, IdxT* idx)
    { return sol[*idx]; }
    static inline void
    set(_Vec& sol, IdxT* idx, double val)
    { sol[*idx]= val; }
};

//**************************************************************************
// Class:   DoFHelperCL                                                    *
// Template Parameter:                                                     *
//          _Len - Length of SVectorCL used as Data                        *
//          Data - SVector<_Len>, partial specialization of the general    *
//                 template                                                *
//          _Vec - The type of the object, that is used for actually       *
//                 storing numerical data, which is referenced by          *
//                 UnknownIdxCL. It is assumed to provide operator[] to    *
//                 access numerical data.                                  *
// Purpose: reads and writes SVectoCL<_Len>-objects from / to _Vec objects *
//          interpreting IdxT as an array of _Len indices in _Vec.         *
//**************************************************************************
template<Uint _Len, class _Vec>
struct DoFHelperCL<SVectorCL<_Len>, _Vec>
{
    static inline SVectorCL<_Len>
    get(const _Vec& sol, const IdxT* idx)
    { SVectorCL<_Len> ret; for (IdxT i=0; i<_Len; ++i) ret[i]= sol[idx[i]]; return ret;}
    static inline void
    set(_Vec& sol, const IdxT* idx, const SVectorCL<_Len>& val)
    { for (IdxT i=0; i<_Len; ++i) sol[idx[i]]= val[i]; }
};


//**************************************************************************
// Class:   P1EvalBaseCL                                                   *
// Template Parameter:                                                     *
//          Data     - The result-type of this finite-element-function on  *
//                     the multigrid                                       *
//          _BndData - Class-type that contains functions that describe,   *
//                     how to handle values in boundary-simplices:         *
//                     bool IsOnDirBnd(VertexCL) - iff true, we use        *
//                     Data GetDirBndValue(T) to obtain the function value *
//          _VD      - (const-) VecDescBaseCL<> like type - type of the    *
//                     container that stores numerical data of this finite-*
//                     element function. The class must contain the typedef*
//                     DataType representing the type used for storing     *
//                     numerical data.                                     *
// Purpose: Abstraction that represents boundary-data and VecDescCL-objects*
//          as a function on the multigrid, that can be evaluated on       *
//          vertices, edges, faces and tetrahedrons via val()-functions and*
//          coordinates in the reference tetrahedron.                      *
//          Degrees of freedom can be set, via SetDoF.                     *
//          We provide GetDoF(S) for simplices S to store all relevant     *
//          numerical data via push_back() in a container. This container  *
//          can be passed to a special val() function and allows for faster*
//          evaluation of the FE-function, if several evaluations on the   *
//          same simpex are necessary.                                     *
//          Generally, evaluations on lower-dimensional simplices are      *
//          faster as only a smaller amount of shape-functions has to be   *
//          evaluated.                                                     *
//**************************************************************************
template<class Data, class _BndData, class _VD>
class P1EvalBaseCL
{
public:
    typedef Data     DataT;
    typedef _BndData BndDataCL;
    typedef _VD      VecDescT;
    
protected:
    // numerical data
    VecDescT*          _sol;
    // boundary-data
    BndDataCL*         _bnd;
    // the multigrid
    const MultiGridCL* _MG;

    virtual inline DataT // helper-function to evaluate the degree of freedom in a vertex,
    GetDoF(const VertexCL&) const= 0;// use val() instead
       
public:
    P1EvalBaseCL() :_sol(0), _MG(0) {}
    P1EvalBaseCL(_VD* sol, _BndData* bnd, const MultiGridCL* MG)
        :_sol(sol), _bnd(bnd), _MG(MG) {}
    // default copy-ctor, dtor, assignment-op
    // copying P1EvalBaseCL-objects is safe - it is a flat copy, which is fine,
    // as P1EvalBaseCL does not take possession of the pointed to _sol, _bnd and _MG.

    void // set / get the container of numerical data
    SetSolution(VecDescT* sol)
        { _sol= sol; }
    VecDescT*
    GetSolution() const
        { return _sol; }
    void // set / get the container of boundary-data
    SetBndData(BndDataCL* bnd)
        { _bnd= bnd; }
    BndDataCL*
    GetBndData() const
        { return _bnd; }
    const MultiGridCL&
    GetMG() const // the multigrid we refer to
        { return *_MG; }

    // evaluation on vertices    
    template<class _Cont>
      inline void
      GetDoF(const VertexCL&, _Cont&) const;
    inline void // set the degree of freedom in this vertex; fails for Dirichlet-vertices
    SetDoF(const VertexCL&, const DataT&);
    template<class _Cont>
      inline DataT
      val(const _Cont&) const;
    inline DataT
    val(const VertexCL&) const;
    
    // evaluation on edges
    template<class _Cont>
      inline void
      GetDoF(const EdgeCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double) const;
    inline DataT
    val(const EdgeCL&, double) const;

    // evaluation on the tetrahedron
    template<class _Cont>
      inline void
      GetDoF(const TetraCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double, double, double) const;
    inline DataT
    val(const TetraCL&, double, double, double) const;
};

template<class DataT, class BndDataT, class VecDescT>
class P1EvalCL: public P1EvalBaseCL<DataT, BndDataT, VecDescT>
{
  // without using-instruction helper-function will shadow other 
  // GetDoF member-functions
  using P1EvalBaseCL<DataT, BndDataT, VecDescT>::GetDoF;
  
  private:
    typedef P1EvalBaseCL<DataT, BndDataT, VecDescT> _base;
    
    inline DataT // helper-function to evaluate on a vertex; use val() instead
    GetDoF(const VertexCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s)
                                   : DoFHelperCL<DataT,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
    }
    
public:
    P1EvalCL() : _base() {}
    P1EvalCL(VecDescT* sol, BndDataT* bnd, const MultiGridCL* MG)
        : _base( sol, bnd, MG) {}
    //default copy-ctor, dtor, assignment-op
};

template<class DataT, class BndDataT, class VecDescT>
class InstatP1EvalCL: public P1EvalBaseCL<DataT, BndDataT, VecDescT>
{
  private:
    typedef P1EvalBaseCL<DataT, BndDataT, VecDescT> _base;
    
    double _t;
    
    inline DataT // helper-function to evaluate on a vertex; use val() instead
    GetDoF(const VertexCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s, _t)
                                   : DoFHelperCL<DataT,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
    }
    
public:
    InstatP1EvalCL() : _base(), _t(0) {}
    InstatP1EvalCL(VecDescT* sol, BndDataT* bnd, const MultiGridCL* MG, double t)
        : _base( sol, bnd, MG), _t(t) {}
    //default copy-ctor, dtor, assignment-op
};

template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1EvalBaseCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c.push_back( GetDoF(s) );
}

template<class Data, class _BndData, class _VD>
inline void
P1EvalBaseCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert(!_bnd->IsOnDirBnd(s), DROPSErrCL("P1EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(_sol->Data, s.Unknowns(_sol->RowIdx->Idx), d);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return  *c.begin();
}

template<class Data, class _BndData, class _VD>
inline Data
P1EvalBaseCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF(s);
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1EvalBaseCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c.push_back( GetDoF(*s.GetVertex(0)) );
    c.push_back( GetDoF(*s.GetVertex(1)) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    typename _Cont::const_iterator it= c.begin();
    DataT ret= *it++ * FE_P1CL::H0(v1);
    return ret + *it * FE_P1CL::H1(v1);
}

template<class Data, class _BndData, class _VD>
inline Data
P1EvalBaseCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1CL::H0(v1)
          + GetDoF(*s.GetVertex(1))*FE_P1CL::H1(v1);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1EvalBaseCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    c.push_back( GetDoF(*s.GetVertex(0)) );
    c.push_back( GetDoF(*s.GetVertex(1)) );
    c.push_back( GetDoF(*s.GetVertex(2)) );
    c.push_back( GetDoF(*s.GetVertex(3)) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    typename _Cont::const_iterator it= c.begin();
    DataT ret= *it++ * FE_P1CL::H0(v1, v2, v3);
    ret+= *it++ * FE_P1CL::H1(v1, v2, v3);
    ret+= *it++ * FE_P1CL::H2(v1, v2, v3);
    return ret + *it * FE_P1CL::H3(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P1EvalBaseCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P1CL::H0(v1, v2, v3)
           +GetDoF(*s.GetVertex(1))*FE_P1CL::H1(v1, v2, v3)
           +GetDoF(*s.GetVertex(2))*FE_P1CL::H2(v1, v2, v3)
           +GetDoF(*s.GetVertex(3))*FE_P1CL::H3(v1, v2, v3);
}



//**************************************************************************
// Class:   P2EvalBaseCL                                                   *
// Template Parameter:                                                     *
//          Data     - The result-type of this finite-element-function on  *
//                     the multigrid                                       *
//          _BndData - Class-type that contains functions that describe,   *
//                     how to handle values in boundary-simplices:         *
//                     bool IsOnDirBnd(VertexCL) - iff true, we use        *
//                     Data GetDirBndValue(T) to obtain the function value *
//          _VD      - (const-) VecDescBaseCL<> like type - type of the    *
//                     container that stores numerical data of this finite-*
//                     element function. The class must contain the typedef*
//                     DataType representing the type used for storing     *
//                     numerical data.                                     *
// Purpose: Abstraction that represents boundary-data and VecDescCL-objects*
//          as a function on the multigrid, that can be evaluated on       *
//          vertices, edges, faces and tetrahedrons via val()-functions and*
//          coordinates in the reference tetrahedron.                      *
//          Degrees of freedom can be set, via SetDoF.                     *
//          We provide GetDoF(S) for simplices S to store all relevant     *
//          numerical data via push_back() in a container. This container  *
//          can be passed to a special val() function and allows for faster*
//          evaluation of the FE-function, if several evaluations on the   *
//          same simpex are necessary.                                     *
//          Generally, evaluations on lower-dimensional simplices are      *
//          faster as only a smaller amount of shape-functions has to be   *
//          evaluated.                                                     *
//**************************************************************************
template<class Data, class _BndData, class _VD>
class P2EvalBaseCL
{
public:
    typedef Data     DataT;
    typedef _BndData BndDataCL;
    typedef _VD      VecDescT;

protected:
    // numerical data
    VecDescT*          _sol;
    // boundary-data
    BndDataCL*         _bnd;
    // the multigrid
    const MultiGridCL* _MG;

    // TODO: Ouch!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    virtual inline DataT // helper-function to evaluate on a vertex; use val() instead
    GetDoF(const VertexCL&) const= 0;
    virtual inline DataT // helper-function to evaluate on an edge; use val() instead
    GetDoF(const EdgeCL&) const= 0;
    
public:
    P2EvalBaseCL() :_sol(0), _bnd(0), _MG(0) {}
    P2EvalBaseCL(VecDescT* sol, BndDataCL* bnd, const MultiGridCL* MG)
        :_sol(sol), _bnd(bnd), _MG(MG) {}
    //default copy-ctor, dtor, assignment-op
    // copying P1EvalBaseCL-objects is safe - it is a flat copy, which is fine,
    // as P1EvalBaseCL does not take possession of the pointed to _sol, _bnd and _MG.

    void // set / get the container of numerical data
    SetSolution(VecDescT* sol)
        { _sol= sol; }
    VecDescT*
    GetSolution() const
        { return _sol; }
    void // set / get the container of boundary-data
    SetBndData(BndDataCL* bnd)
        { _bnd= bnd; }
    BndDataCL*
    GetBndData() const
        { return _bnd; }
    const MultiGridCL&
    GetMG() const // the multigrid we refer to
        { return *_MG; }

    // evaluation on vertices    
    inline void // set the degree of freedom in the vertex; fails if, we are on a Dirichlet-boundary
    SetDoF(const VertexCL&, const DataT&);
    template<class _Cont>
      inline void
      GetDoF(const VertexCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&) const;
    inline DataT
    val(const VertexCL&) const;
    
    // evaluation on edges
    inline void // set the degree of freedom on the edge; fails if, we are on a Dirichlet-boundary
    SetDoF(const EdgeCL&, const DataT&);
    template<class _Cont>
      inline void
      GetDoF(const EdgeCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double) const;
    inline DataT
    val(const EdgeCL&, double) const;
    inline DataT // for your convenience: val(edge) == val(edge, 0.5)
    val(const EdgeCL&) const;

    // FIXME: evaluation on faces are currently bogus
    template<class _Cont>
      inline void
      GetDoF(const TetraCL&, Uint, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double, double) const;
    inline DataT
    val(const TetraCL&, Uint, double, double) const;

    // evaluation on a tetrahedron
    template<class _Cont>
      inline void
      GetDoF(const TetraCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double, double, double) const;
    inline DataT
    val(const TetraCL&, double, double, double) const;

    inline bool UnknownsMissing(const TetraCL&) const;
};

template<class DataT, class BndDataT, class VecDescT>
class P2EvalCL: public P2EvalBaseCL<DataT, BndDataT, VecDescT>
{
  private:
    typedef P2EvalBaseCL<DataT, BndDataT, VecDescT> _base;
    using _base::_bnd;
    using _base::_sol;
    
    inline DataT // helper-function to evaluate on a vertex; use val() instead
    GetDoF(const VertexCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s)
                                   : DoFHelperCL<DataT,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
    }
    inline DataT // helper-function to evaluate on an edge; use val() instead
    GetDoF(const EdgeCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s)
                                   : DoFHelperCL<DataT,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
    }
    
public:
    using _base::GetDoF;

    P2EvalCL() : _base() {}
    P2EvalCL(VecDescT* sol, BndDataT* bnd, const MultiGridCL* MG)
        : _base( sol, bnd, MG) {}
    //default copy-ctor, dtor, assignment-op
};

template<class DataT, class BndDataT, class VecDescT>
class InstatP2EvalCL: public P2EvalBaseCL<DataT, BndDataT, VecDescT>
{
  private:
    typedef P2EvalBaseCL<DataT, BndDataT, VecDescT> _base;
    using _base::_bnd;
    using _base::_sol;
    
    double _t;
    
    inline DataT // helper-function to evaluate on a vertex; use val() instead
    GetDoF(const VertexCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s, _t)
                                   : DoFHelperCL<DataT,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
    }
    inline DataT // helper-function to evaluate on an edge; use val() instead
    GetDoF(const EdgeCL& s) const
    {
        return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s, _t)
                                   : DoFHelperCL<DataT,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
    }
    
public:
    InstatP2EvalCL() : _base(), _t(0) {}
    InstatP2EvalCL(VecDescT* sol, BndDataT* bnd, const MultiGridCL* MG, double t)
        : _base( sol, bnd, MG), _t(t) {}
    //default copy-ctor, dtor, assignment-op
};


template<class Data, class _BndData, class _VD>
inline void
P2EvalBaseCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert(!_bnd->IsOnDirBnd(s), DROPSErrCL("P2EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(_sol->Data, s.Unknowns(_sol->RowIdx->Idx), d);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P2EvalBaseCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c.push_back( GetDoF(s) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return *c.begin();
}

template<class Data, class _BndData, class _VD>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const VertexCL& s) const
{
    return GetDoF(s);
}


template<class Data, class _BndData, class _VD>
inline void
P2EvalBaseCL<Data, _BndData, _VD>::SetDoF(const EdgeCL& s, const Data& d)
{
    Assert(!_bnd->IsOnDirBnd(s), DROPSErrCL("P2EvalBaseCL::SetDoF: Trying to assign to Dirichlet-boundary-edge."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(_sol->Data, s.Unknowns(_sol->RowIdx->Idx), d);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P2EvalBaseCL<Data, _BndData, _VD>::GetDoF(const EdgeCL& s, _Cont& c) const
{
    c.push_back( GetDoF(*s.GetVertex(0)) );
    c.push_back( GetDoF(*s.GetVertex(1)) );
    c.push_back( GetDoF(s) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    typename _Cont::const_iterator it= c.begin();
    DataT ret= *it++ * FE_P2CL::H0(v1);
    ret+= *it++ * FE_P2CL::H1(v1);
    return ret + *it * FE_P2CL::H2(v1);
}

template<class Data, class _BndData, class _VD>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const EdgeCL& s, double v1) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P2CL::H0(v1)
          + GetDoF(*s.GetVertex(1))*FE_P2CL::H1(v1)
          + GetDoF(s)*FE_P2CL::H2(v1);
}

template<class Data, class _BndData, class _VD>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const EdgeCL& s) const
{
    return GetDoF(s);
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P2EvalBaseCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, Uint face, _Cont& c) const
{
    typename _Cont::iterator it= c.begin();
    *it++= GetDoF(*s.GetVertex(VertOfFace(face, 0)));
    *it++= GetDoF(*s.GetVertex(VertOfFace(face, 1)));
    *it++= GetDoF(*s.GetVertex(VertOfFace(face, 2)));
    *it++= GetDoF(*s.GetEdge(EdgeOfFace(face, 0)));
    *it++= GetDoF(*s.GetEdge(EdgeOfFace(face, 1)));
    *it= GetDoF(*s.GetEdge(EdgeOfFace(face, 2)));
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2) const
{
    typename _Cont::const_iterator it= c.begin();

    DataT ret= *it++ * FE_P2CL::H0(v1, v2);
    ret+= *it++ * FE_P2CL::H1(v1, v2);
    ret+= *it++ * FE_P2CL::H2(v1, v2);
    ret+= *it++ * FE_P2CL::H3(v1, v2);
    ret+= *it++ * FE_P2CL::H4(v1, v2);
    return ret + *it * FE_P2CL::H5(v1, v2);
}


template<class Data, class _BndData, class _VD>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const TetraCL& s, Uint face, double v1, double v2) const
{
    return  GetDoF(*s.GetVertex(VertOfFace(face, 0)))*FE_P2CL::H0(v1, v2)
           +GetDoF(*s.GetVertex(VertOfFace(face, 1)))*FE_P2CL::H1(v1, v2)
           +GetDoF(*s.GetVertex(VertOfFace(face, 2)))*FE_P2CL::H2(v1, v2)
           +GetDoF(*s.GetEdge(EdgeOfFace(face, 0)))*FE_P2CL::H3(v1, v2)
           +GetDoF(*s.GetEdge(EdgeOfFace(face, 1)))*FE_P2CL::H4(v1, v2)
           +GetDoF(*s.GetEdge(EdgeOfFace(face, 2)))*FE_P2CL::H5(v1, v2);
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P2EvalBaseCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s, _Cont& c) const
{
    c.push_back( GetDoF(*s.GetVertex(0)) );
    c.push_back( GetDoF(*s.GetVertex(1)) );
    c.push_back( GetDoF(*s.GetVertex(2)) );
    c.push_back( GetDoF(*s.GetVertex(3)) );
    c.push_back( GetDoF(*s.GetEdge(0)) );
    c.push_back( GetDoF(*s.GetEdge(1)) );
    c.push_back( GetDoF(*s.GetEdge(2)) );
    c.push_back( GetDoF(*s.GetEdge(3)) );
    c.push_back( GetDoF(*s.GetEdge(4)) );
    c.push_back( GetDoF(*s.GetEdge(5)) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    typename _Cont::const_iterator it= c.begin();

    DataT ret= *it++ * FE_P2CL::H0(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H1(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H2(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H3(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H4(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H5(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H6(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H7(v1, v2, v3);
    ret+= *it++ * FE_P2CL::H8(v1, v2, v3);
    return ret + *it * FE_P2CL::H9(v1, v2, v3);
}

template<class Data, class _BndData, class _VD>
inline Data
P2EvalBaseCL<Data, _BndData, _VD>::val(const TetraCL& s, double v1, double v2, double v3) const
{
    return  GetDoF(*s.GetVertex(0))*FE_P2CL::H0(v1, v2, v3)
           +GetDoF(*s.GetVertex(1))*FE_P2CL::H1(v1, v2, v3)
           +GetDoF(*s.GetVertex(2))*FE_P2CL::H2(v1, v2, v3)
           +GetDoF(*s.GetVertex(3))*FE_P2CL::H3(v1, v2, v3)
           +GetDoF(*s.GetEdge(0))*FE_P2CL::H4(v1, v2, v3)
           +GetDoF(*s.GetEdge(1))*FE_P2CL::H5(v1, v2, v3)
           +GetDoF(*s.GetEdge(2))*FE_P2CL::H6(v1, v2, v3)
           +GetDoF(*s.GetEdge(3))*FE_P2CL::H7(v1, v2, v3)
           +GetDoF(*s.GetEdge(4))*FE_P2CL::H8(v1, v2, v3)
           +GetDoF(*s.GetEdge(5))*FE_P2CL::H9(v1, v2, v3);
}


template<class Data, class _BndData, class _VD>
inline bool P2EvalBaseCL<Data, _BndData, _VD>::UnknownsMissing(const TetraCL& t) const
{
    const Uint idx= _sol->RowIdx->Idx;
    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd();
        it!=end; ++it)
        if ( !(_bnd->IsOnDirBnd( **it) || (*it)->Unknowns.Exist(idx) ))
            return true;
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd();
        it!=end; ++it)
        if ( !(_bnd->IsOnDirBnd( **it) || (*it)->Unknowns.Exist(idx) ))
            return true;
    return false;
}


//**************************************************************************
// Class:   P1BubbleEvalCL                                                 *
// Template Parameter:                                                     *
//          Data     - The result-type of this finite-element-function on  *
//                     the multigrid                                       *
//          _BndData - Class-type that contains functions that describe,   *
//                     how to handle values in boundary-simplices:         *
//                     bool IsOnDirBnd(VertexCL) - iff true, we use        *
//                     Data GetDirBndValue(T) to obtain the function value *
//          _VD      - (const-) VecDescBaseCL<> like type - type of the    *
//                     container that stores numerical data of this finite-*
//                     element function. The class must contain the typedef*
//                     DataType representing the type used for storing     *
//                     numerical data.                                     *
// Purpose: Abstraction that represents boundary-data and VecDescCL-objects*
//          as a function on the multigrid, that can be evaluated on       *
//          vertices, edges, faces and tetrahedrons via val()-functions and*
//          coordinates in the reference tetrahedron.                      *
//          Degrees of freedom can be set, via SetDoF.                     *
//          We provide GetDoF(S) for simplices S to store all relevant     *
//          numerical data via push_back() in a container. This container  *
//          can be passed to a special val() function and allows for faster*
//          evaluation of the FE-function, if several evaluations on the   *
//          same simpex are necessary.                                     *
//          Generally, evaluations on lower-dimensional simplices are      *
//          faster as only a smaller amount of shape-functions has to be   *
//          evaluated.                                                     *
//**************************************************************************
template<class Data, class _BndData, class _VD>
class P1BubbleEvalCL
{
public:
    typedef Data     DataT;
    typedef _BndData BndDataCL;
    typedef _VD      VecDescT;
    
private:
    // numerical data
    VecDescT*          _sol;
    // boundary-data
    BndDataCL*         _bnd;
    // the multigrid
    const MultiGridCL* _MG;

    inline DataT // helper-function to evaluate the degree of freedom in a vertex,
    GetDoF(const VertexCL&) const;// use val() instead
    inline DataT // helper-function to evaluate the degree of freedom in a tetra
    GetDoF(const TetraCL&) const;
       
public:
    P1BubbleEvalCL() :_sol(0), _MG(0) {}
    P1BubbleEvalCL(_VD* sol, _BndData* bnd, const MultiGridCL* MG)
        :_sol(sol), _bnd(bnd), _MG(MG) {}
    // default copy-ctor, dtor, assignment-op
    // copying P1BubbleEvalCL-objects is safe - it is a flat copy, which is fine,
    // as P1BubbleEvalCL does not take possession of the pointed to _sol, _bnd and _MG.

    void // set / get the container of numerical data
    SetSolution(VecDescT* sol)
        { _sol= sol; }
    VecDescT*
    GetSolution() const
        { return _sol; }
    void // set / get the container of boundary-data
    SetBndData(BndDataCL* bnd)
        { _bnd= bnd; }
    BndDataCL*
    GetBndData() const
        { return _bnd; }
    const MultiGridCL&
    GetMG() const // the multigrid we refer to
        { return *_MG; }

    // evaluation on vertices    
    template<class _Cont>
      inline void
      GetDoF(const VertexCL&, _Cont&) const;
    inline void // set the degree of freedom in this vertex; fails for Dirichlet-vertices
    SetDoF(const VertexCL&, const DataT&);
    template<class _Cont>
      inline DataT
      val(const _Cont&) const;
    inline DataT
    val(const VertexCL&) const;
    
    // evaluation on edges
    template<class _Cont>
      inline void
      GetDoF(const EdgeCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double) const;
    inline DataT
    val(const EdgeCL&, double) const;

    // evaluation on the tetrahedron
    template<class _Cont>
      inline void
      GetDoF(const TetraCL&, _Cont&) const;
    template<class _Cont>
      inline DataT
      val(const _Cont&, double, double, double) const;
    inline DataT
    val(const TetraCL&, double, double, double) const;
    inline DataT
    lin_val(const TetraCL&, double, double, double) const;
    inline DataT
    bubble_val(const TetraCL&, double, double, double) const;
};


template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s) const
{
    return _bnd->IsOnDirBnd(s) ? _bnd->GetDirBndValue(s)
                               : DoFHelperCL<Data,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
}


template<class Data, class _BndData, class _VD> template<class _Cont>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const VertexCL& s, _Cont& c) const
{
    c.push_back( GetDoF(s) );
}


template<class Data, class _BndData, class _VD>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::GetDoF(const TetraCL& s) const
{
    return DoFHelperCL<Data,typename VecDescT::DataType>::get(_sol->Data, s.Unknowns(_sol->RowIdx->Idx));
}


template<class Data, class _BndData, class _VD>
inline void
P1BubbleEvalCL<Data, _BndData, _VD>::SetDoF(const VertexCL& s, const Data& d)
{
    Assert(!_bnd->IsOnDirBnd(s), DROPSErrCL("P1EvalCL::SetDoF: Trying to assign to Dirichlet-boundary-vertex."), DebugNumericC);
    DoFHelperCL<Data, typename VecDescT::DataType>::set(_sol->Data, s.Unknowns(_sol->RowIdx->Idx), d);
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c) const
{
    return  *c.begin();
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
    c.push_back( GetDoF(*s.GetVertex(0)) );
    c.push_back( GetDoF(*s.GetVertex(1)) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1) const
{
    typename _Cont::const_iterator it= c.begin();
    DataT ret= *it++ * FE_P1BubbleCL::H0(v1);
    return ret + *it * FE_P1BubbleCL::H1(v1);
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
    c.push_back( GetDoF(*s.GetVertex(0)) );
    c.push_back( GetDoF(*s.GetVertex(1)) );
    c.push_back( GetDoF(*s.GetVertex(2)) );
    c.push_back( GetDoF(*s.GetVertex(3)) );
    c.push_back( GetDoF(s) );
}

template<class Data, class _BndData, class _VD> template<class _Cont>
inline Data
P1BubbleEvalCL<Data, _BndData, _VD>::val(const _Cont& c, double v1, double v2, double v3) const
{
    typename _Cont::const_iterator it= c.begin();
    DataT ret= *it++ * FE_P1BubbleCL::H0(v1, v2, v3);
    ret+= *it++ * FE_P1BubbleCL::H1(v1, v2, v3);
    ret+= *it++ * FE_P1BubbleCL::H2(v1, v2, v3);
    ret+= *it++ * FE_P1BubbleCL::H3(v1, v2, v3);
    return ret + *it * FE_P1BubbleCL::H4(v1, v2, v3);
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
    const Uint old_level= old_sol.GetSolution()->RowIdx->TriangLevel;
    const Uint level= sol.GetSolution()->RowIdx->TriangLevel;
    const Uint old_idx= old_sol.GetSolution()->RowIdx->Idx;
    Uint counter1= 0, counter2= 0;
    
//    Assert( level==old_level || old_level==level-1, DROPSErrCL("Interpolate: successive triangs are expected\n"));
    // Iterate over all edges, interpolate values on new mid vertices
    for (MultiGridCL::const_AllEdgeIteratorCL sit= _MG.GetAllEdgeBegin(level), theend= _MG.GetAllEdgeEnd(level);
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


inline void
GetEdgeInFace(Uint chedge, Uint& parface, Uint& pos)
// Returns, on which face of the parent chedge lies and also determines,
// whether the barycenter of chedge lies at v0 (0), v1 (1) or v2 (2) of the face.
// As to why why this works, take a look at topo.h, specifically the sorting
// of VertOfEdgeAr. 
{
    parface= chedge<30 ? chedge/3 -6 : chedge/3 -10;
    pos= chedge%3;
}

template<class Data, class _BndData, class _VD>
void InterpolateChildren( const TetraCL& t, P2EvalBaseCL<Data, _BndData, _VD>& sol, const P2EvalBaseCL<Data, _BndData, const _VD>& old_sol)
{
    typedef typename P2EvalBaseCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();

    const double edgebary[3][2]=
        { {0.25, 0.25},
          {0.5 , 0.25},
          {0.25, 0.5}
        };
      
    // Hole des Tetraeders RefRule; gehe ueber alle Kinder ;-): durchlaufe die edges in
    // der childrule des kindes: falls IsSubEdge(edge): finde mit ParentEdge & NumOfSubEdge heraus,
    // von welcher Kante des tetras, (falls der wert bei edge != 0) und interpoliere ueber kante. 
    // sonst, falls IsSubInParFace(subedge): GetEdgeInFace(subedge, face, pos) und interpoliere ueber face;
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
                    GetEdgeInFace(chedgeinparent, parface, pos);
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
inline void CopyValues( const TetraCL& t, P2EvalBaseCL<Data, _BndData, _VD>& sol, const P2EvalBaseCL<Data, _BndData, const _VD>& old_sol)
{
    typedef typename P2EvalBaseCL<Data, _BndData, _VD>::BndDataCL BndCL;
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
void Adapt( P2EvalBaseCL<Data, _BndData, _VD>& sol, const P2EvalBaseCL<Data, _BndData, const _VD>& old_sol)
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
    const Uint level= sol.GetSolution()->RowIdx->TriangLevel;
    const Uint old_idx= old_sol.GetSolution()->RowIdx->Idx;
    const Uint NumUnknowns=  old_sol.GetSolution()->RowIdx->NumUnknownsEdge;
    
    Assert( level==old_sol.GetSolution()->RowIdx->TriangLevel, 
        DROPSErrCL("Adapt: Same triang levels are expected\n"), -1u);

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
                         if (!(*it)->Unknowns.Exist()) (*it)->Unknowns.Init(old_idx+1);
                         else if ( !(old_idx < (*it)->Unknowns.Get()->GetNumSystems()) )
                             (*it)->Unknowns.Get()->resize(old_idx+1); 
                         if ( !(*it)->Unknowns.Exist(old_idx) )
                             (*it)->Unknowns.Get()->Alloc(old_idx, NumUnknowns);
                         // Indexwerte von MidVertex kopieren
                         for (Uint i=0; i<NumUnknowns; ++i)
                             (*it)->Unknowns(old_idx)[i]= (*it)->GetMidVertex()->Unknowns(old_idx)[i];
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
void Interpolate(P2EvalBaseCL<Data, _BndData, _VD>& sol, const P2EvalBaseCL<Data, _BndData, const _VD>& old_sol)
// This only works, if Interpolate is called after every refinement of the multigrid.
// Take care, that x and old_x are on successive triangulations.
{
    typedef typename P2EvalBaseCL<Data, _BndData, _VD>::BndDataCL BndCL;
    const BndCL* const _bnd= old_sol.GetBndData();
    const Uint old_idx= old_sol.GetSolution()->RowIdx->Idx;
      
    // All velocity-components use the same row-index and the same trianglevel
    const MultiGridCL& _MG= old_sol.GetMG();    
    const Uint old_level= old_sol.GetSolution()->RowIdx->TriangLevel;
    const Uint level= sol.GetSolution()->RowIdx->TriangLevel;
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
    for (MultiGridCL::const_AllEdgeIteratorCL sit= _MG.GetAllEdgeBegin(level), theend= _MG.GetAllEdgeEnd(level);
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
                                     << num_child_edge << " edge-DoF interpolated." << std::endl;
}


} // end of namespace DROPS

#endif
