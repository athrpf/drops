//**************************************************************************
// File:    discretize.h                                                   *
// Content: discretizations for several PDEs and FE types                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Juli, 10 2001                                          *
//**************************************************************************

#ifndef DROPS_DISCRETIZE_H
#define DROPS_DISCRETIZE_H

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "misc/container.h"

namespace DROPS
{

typedef double (*scalar_fun_ptr)(const Point3DCL&);
typedef Point3DCL (*vector_fun_ptr)(const Point3DCL&);
typedef double (*SmoothFunT) (double,double);


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



// ===================================
//        Quadrature formulas
// ===================================

template<class T=double>
class QuadBaseCL
{
  protected:
    QuadBaseCL( size_t size) : val( size) {}
    QuadBaseCL( size_t size, const T& t) : val( t, size) {}
    
  public:
    std::valarray<T> val;
    
    size_t size() const 
        { return val.size(); }
        
    // Arithmetik
    template<typename FuncT>
    QuadBaseCL& apply ( FuncT fun)
      { for (size_t i=0; i<val.size(); ++i) val[i]=fun(val[i]); return *this; }

    QuadBaseCL& operator+= (const QuadBaseCL &q)
      { val+=q.val; return *this; }
    QuadBaseCL& operator+= (const T &t)
      { val+= t; return *this; }
    QuadBaseCL& operator*= (const T &t)
      { val*= t; return *this; }
    template <typename U> QuadBaseCL& operator*= (const QuadBaseCL<U>& q)
      { for (size_t i=0; i<val.size(); ++i) val[i]*=q.val[i]; return *this; }
//    template <typename U> QuadBaseCL& operator*= (const U& c)
//      { val*=c; return *this; }

    QuadBaseCL operator+ (const QuadBaseCL &q) const
      { return QuadBaseCL(*this)+= q; }
    QuadBaseCL operator+ (const T &t) const
      { return QuadBaseCL(*this)+= t; }
    QuadBaseCL operator* (const T &t) const
      { return QuadBaseCL(*this)*= t; }
    template <typename U> QuadBaseCL operator* (const QuadBaseCL<U> &q) const
      { return QuadBaseCL(*this)*= q; }
//    template <typename U> QuadBaseCL operator* (const U& c) const
//      { QuadBaseCL ret(*this); ret.val*=c; return ret; }
    friend QuadBaseCL<double> dot (const QuadBaseCL<Point3DCL>&, const QuadBaseCL<Point3DCL>&);
};

inline QuadBaseCL<double> dot (const QuadBaseCL<Point3DCL> &q1, const QuadBaseCL<Point3DCL> &q2)
{ 
    QuadBaseCL<double> res( q1.size()); 
    for (size_t i=0; i<q1.size(); ++i) 
        res.val[i]= inner_prod( q1.val[i], q2.val[i]); 
    return res; 
}


template<class T=double>
class Quad2CL: public QuadBaseCL<T>
{
  public:
    using QuadBaseCL<T>::size;
    using QuadBaseCL<T>::val;

    static const BaryCoordCL Node[5]; // Stuetzstellen
    static const double      Wght[5]; // Gewichte
    
    Quad2CL() : QuadBaseCL<T>( 5) {}
    Quad2CL( const T& t) : QuadBaseCL<T>( 5, t) {}
    Quad2CL( const QuadBaseCL<T>& q) : QuadBaseCL<T>(q) {}
    
    // Initialisiere die Knotenwerte
/*    template <ElemT et> 
    void set( const ElemBaseCL<et,T>& e)
      { for (size_t i=0; i<size(); ++i) val[i]= e.eval( Node[i]); }
*/
    // Werte Quadraturformel aus
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const
      { T sum= T(); for (size_t i=0; i<size(); ++i) sum+= Wght[i]*val[i]; return sum*absdet; }
    // Folgende Spezialformeln nutzen die spezielle Lage der Stuetzstellen aus
    // zur Annaeherung von \int f*phi,    phi = P1-/P2-Hutfunktion
    T quadP1 (int i, double absdet) const
      { return ((1./120.)*val[i] + (1./30.)*val[4])*absdet; }
    T quadP1 (int i, int j, double absdet) const
      { return (i!=j ? (1./720.)*(val[i]+val[j]) + (1./180.)*val[4]
                     : (1./180.)*val[i]          + (1./90.)*val[4]  )*absdet;}
    T quadP2 (int i, double absdet) const
    { 
        return (i<4 ? (1./360.)*val[i] - (1./90.)*val[4]
                    : (1./180.)*(val[VertOfEdge(i-4,0)]+val[VertOfEdge(i-4,1)]) + (1./45.)*val[4]
               )*absdet;
    }

    T quadP2 (int i, int j, double absdet) const
    { 
        const double valBary= (i<4 ? -0.125 : 0.25)*(j<4 ? -0.125 : 0.25);
        return ((i!=j || i>=4) ? Wght[4]*val[4]*valBary
                               : Wght[4]*val[4]*valBary + Wght[i]*val[i]
               )*absdet;
    }
};


class Quad3CL
// contains cubatur on reference-tetra, that is exact up to degree 3, positive,
// and uses only 8 points.
// Do not forget to multiply the result of Quad() by the determinant of the affine trafo
// from the reference tetra to the given tetra.
{
  private:
    static const double _points[8][3];

  public:
    static Uint GetNumPoints() { return 8; }
    static const Point3DCL* GetPoints() { return reinterpret_cast<const Point3DCL*>(_points); }

    static inline double Quad(const TetraCL&, scalar_fun_ptr);
    static inline SVectorCL<3> Quad(const TetraCL&, vector_fun_ptr);

    static inline double Quad(const double*);
    static inline SVectorCL<3> Quad(const SVectorCL<3>*);
    template<class IteratorT, class ValueT>
      static inline void
      Quad(IteratorT, ValueT*const);
};

inline double Quad3CL::Quad(const TetraCL& t, scalar_fun_ptr f)
{
    const Point3DCL* pts= GetPoints();
    return ( f(GetWorldCoord(t, pts[0])) + f(GetWorldCoord(t, pts[1]))
            +f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )/240.
          +( f(GetWorldCoord(t, pts[4])) + f(GetWorldCoord(t, pts[5]))
            +f(GetWorldCoord(t, pts[6])) + f(GetWorldCoord(t, pts[7])) )*3./80.;

}

inline SVectorCL<3> Quad3CL::Quad(const TetraCL& t, vector_fun_ptr f)
{
    const Point3DCL* pts= GetPoints();
    return ( f(GetWorldCoord(t, pts[0])) + f(GetWorldCoord(t, pts[1]))
            +f(GetWorldCoord(t, pts[2])) + f(GetWorldCoord(t, pts[3])) )/240.
          +( f(GetWorldCoord(t, pts[4])) + f(GetWorldCoord(t, pts[5]))
            +f(GetWorldCoord(t, pts[6])) + f(GetWorldCoord(t, pts[7])) )*3./80.;

}

inline double Quad3CL::Quad(const double* vals)
{
    return (vals[0] + vals[1] + vals[2] + vals[3])/240.
          +(vals[4] + vals[5] + vals[6] + vals[7])*3./80.;
}

inline SVectorCL<3> Quad3CL::Quad(const SVectorCL<3>* vals)
{
    return (vals[0] + vals[1] + vals[2] + vals[3])/240.
          +(vals[4] + vals[5] + vals[6] + vals[7])*3./80.;
}

template<class IteratorT, class ValueT>
inline void Quad3CL::Quad(IteratorT beg, ValueT*const ret)
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
    static const Point2DCL* GetPoints() { return reinterpret_cast<const Point2DCL*>(_points); }

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

/*
template <class Fun>
struct FunTraitsCL
{
    typedef void argument_type;
    typedef void return_type;
};
template <>
struct FunTraitsCL<scalar_fun_ptr>
{
    typedef Point3DCL argument_type;
    typedef double    return_type;
};
template <>
struct FunTraitsCL<vector_fun_ptr>
{
    typedef Point3DCL    argument_type;
    typedef SVectorCL<3> return_type;
};
*/


//=========================================
//    Finite Elements: P1, P1Bubble, P2
//=========================================

class P1DiscCL
// contains cubatur etc. for linear FE
{
  public:
    // cubatur formula for int f(x) dx, exact up to degree 2
    static inline double Quad(const TetraCL&, scalar_fun_ptr);
    // cubatur formula for int f(x)*phi_i dx, exact up to degree 1
    static inline double Quad(const TetraCL&, scalar_fun_ptr, Uint);
    static inline SVectorCL<3> Quad(const TetraCL&, vector_fun_ptr, Uint);
    // cubatur formula for int f(x)*phi_i*phi_j dx, exact up to degree 1
    static inline double Quad(const TetraCL&, scalar_fun_ptr, Uint, Uint);
    // computes the square of the L2-norm of a given function f: 
    // f^2 is integrated exact up to degree 2
    static inline double norm_L2_sq(const TetraCL&, scalar_fun_ptr);

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
    static void GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10]);
    // compute gradients
    static void GetGradients( Quad2CL<Point3DCL> G[10], Quad2CL<Point3DCL> GRef[10], SMatrixCL<3,3> &T)
    { for (int i=0; i<10; ++i) for (int j=0; j<5; ++j) G[i].val[j]= T*GRef[i].val[j]; }
    static void GetGradient( Quad2CL<Point3DCL> &G, Quad2CL<Point3DCL> &GRef, SMatrixCL<3,3> &T)
    { for (int j=0; j<5; ++j) G.val[j]= T*GRef.val[j]; }
};


inline double FuncDet2D( const Point3DCL& p, const Point3DCL& q)
{
    const double d0= p[1]*q[2] - p[2]*q[1];
    const double d1= p[2]*q[0] - p[0]*q[2];
    const double d2= p[0]*q[1] - p[1]*q[0];
    return sqrt(d0*d0 + d1*d1 + d2*d2);
}

/*
template<class T=double>
class P1ElemCL
{
  private:
    valarray<T> val;
    
  public:
    P1ElemCL()
      : val(4) {}
    
    const T& operator[] (int i) const { return val[i]; }
    
    T eval( const BaryCoordCL& c) const
    {
        T sum= T();
        for (int i=0; i<4; ++i) sum+= val[i]*c[i];
        return sum;
    }
};
*/


/********************************************************************************
*
*        definition of   i n l i n e   f u n c t i o n s
*
********************************************************************************/


inline double P1DiscCL::Quad(const TetraCL& t, scalar_fun_ptr coeff)
{
    return ( coeff(t.GetVertex(0)->GetCoord())
            +coeff(t.GetVertex(1)->GetCoord())
            +coeff(t.GetVertex(2)->GetCoord())
            +coeff(t.GetVertex(3)->GetCoord()))/120. 
            + 2./15.*coeff(GetBaryCenter(t));
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

template<class T>
const BaryCoordCL Quad2CL<T>::Node[5]= {
    {1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}, {.25,.25,.25,.25}
    }; 
template<class T>
const double Quad2CL<T>::Wght[5]= { 1./120., 1./120., 1./120., 1./120., 2./15.};

} // end of namespace DROPS

#endif
