//**************************************************************************
// File:    discretize.h                                                   *
// Content: discretizations for several PDEs and FE types                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Juli, 10 2001                                          *
//**************************************************************************

#ifndef _DISCRETIZE_H_
#define _DISCRETIZE_H_

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "misc/container.h"

namespace DROPS
{

typedef double (*scalar_fun_ptr)(const Point3DCL&);
typedef SVectorCL<3> (*vector_fun_ptr)(const Point3DCL&);


class Quad3CL
// contains cubatur on reference-tetra, that is exact up to degree 3, positive,
// and uses only 8 points.
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


inline double FuncDet2D( const Point3DCL& p, const Point3DCL& q)
{
    const double d0= p[1]*q[2] - p[2]*q[1];
    const double d1= p[2]*q[0] - p[0]*q[2];
    const double d2= p[0]*q[1] - p[1]*q[0];
    return sqrt(d0*d0 + d1*d1 + d2*d2);
}


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

} // end of namespace DROPS

#endif
