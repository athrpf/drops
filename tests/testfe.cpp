#include "num/fe.h"
#include "misc/container.h"
#include "geom/topo.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace DROPS;

/*****************************************************************************************************
* formulas for   n u m e r i c   i n t e g r a t i o n   on the reference tetrahedron
*****************************************************************************************************/
typedef double (*scalar_coeff_ptr)(const Point3DCL&);

inline double Quad(const TetraCL& t, scalar_coeff_ptr coeff)
// exact up to degree 2
{
    return ( coeff(t.GetVertex(0)->GetCoord())
            +coeff(t.GetVertex(1)->GetCoord())
            +coeff(t.GetVertex(2)->GetCoord())
            +coeff(t.GetVertex(3)->GetCoord()))/120. 
            + 2./15.*coeff(GetBaryCenter(t));
}

/*inline StokesBndCL::bnd_type Quad(const TetraCL& t, vector_coeff_ptr coeff)
// exact up to degree 2
{
    return ( coeff(t.GetVertex(0)->GetCoord())
            +coeff(t.GetVertex(1)->GetCoord())
            +coeff(t.GetVertex(2)->GetCoord())
            +coeff(t.GetVertex(3)->GetCoord()))/120. 
            + 2./15.*coeff(GetBaryCenter(t));
}
*/
inline double QuadGrad(const SMatrixCL<3,5>* G, int i, int j)
// computes int( grad(phi_i) * grad(phi_j) ) for P2-elements on ref. tetra
{
    SVectorCL<5> tmp(0.0);
    
    for(int k=0; k<5; ++k)
        for(int l=0; l<3; ++l)
            tmp[k]+= G[i](l,k) * G[j](l,k);
            
    return ( tmp[0] + tmp[1] + tmp[2] + tmp[3] )/120. + 2./15.*tmp[4];
}

//inline SVectorCL<3> Quad( const TetraCL& t, vector_coeff_ptr coeff, int i)
inline double Quad( const TetraCL& t, scalar_coeff_ptr coeff, int i)
{
    double f[5];
    
    if (i<4) // hat function on vert
    {
        f[0]= coeff( t.GetVertex(i)->GetCoord() );
        for (int k=0, l=1; k<4; ++k)
            if (k!=i) f[l++]= coeff( t.GetVertex(k)->GetCoord() );
        f[4]= coeff( GetBaryCenter(t) );
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

        double sum= vs * coeff( GetBaryCenter(t) );
        for(int k=0; k<4; ++k)
            sum+= a[k] * coeff( t.GetVertex(k)->GetCoord() );

        return sum;
    }
}


// TODO: folgende Quadraturformeln anpassen an P2

inline double Quad( const TetraCL& t, scalar_coeff_ptr f, int i, int j)
{
    double a[5];
    if (i>j) std::swap(i,j);
    switch(i*10+j)
    {
      case  0: a[0]= 1./1260.; a[1]= a[2]= a[3]= 0.; a[4]= 1./630.; break;
      case  1: a[0]= a[1]= -1./8505.; a[2]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case  2: a[0]= a[2]= -1./8505.; a[1]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case  3: a[0]= a[3]= -1./8505.; a[1]= a[2]= 11./136080.; a[4]= 4./8505; break;
      case  4: a[0]= 1./2520.; a[1]= -1./2520.; a[2]= a[3]= 0.; a[4]= -1./630.; break;
      case  5: a[0]= 1./2520.; a[2]= -1./2520.; a[1]= a[3]= 0.; a[4]= -1./630.; break;
      case  6: a[0]= a[3]= 1./8505.; a[1]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case  7: a[0]= 1./2520.; a[3]= -1./2520.; a[1]= a[2]= 0.; a[4]= -1./630.; break;
      case  8: a[0]= a[2]= 1./8505.; a[1]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case  9: a[0]= a[1]= 1./8505.; a[2]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 11: a[1]= 1./1260.; a[0]= a[2]= a[3]= 0.; a[4]= 1./630.; break;
      case 12: a[1]= a[2]= -1./8505.; a[0]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case 13: a[1]= a[3]= -1./8505.; a[0]= a[2]= 11./136080.; a[4]= 4./8505; break;
      case 14: a[1]= 1./2520.; a[0]= -1./2520.; a[2]= a[3]= 0.; a[4]= -1./630.; break;
      case 15: a[1]= a[3]= 1./8505.; a[0]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 16: a[1]= 1./2520.; a[2]= -1./2520.; a[0]= a[3]= 0.; a[4]= -1./630.; break;
      case 17: a[1]= a[2]= 1./8505.; a[0]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 18: a[1]= 1./2520.; a[3]= -1./2520.; a[0]= a[2]= 0.; a[4]= -1./630.; break;
      case 19: a[0]= a[1]= 1./8505.; a[2]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 22: a[2]= 1./1260.; a[0]= a[1]= a[3]= 0.; a[4]= 1./630.; break;
      case 23: a[2]= a[3]= -1./8505.; a[0]= a[1]= 11./136080.; a[4]= 4./8505; break;
      case 24: a[2]= a[3]= 1./8505.; a[0]= a[1]= -19./68040.; a[4]= -1./486.; break;
      case 25: a[2]= 1./2520.; a[0]= -1./2520.; a[1]= a[3]= 0.; a[4]= -1./630.; break;
      case 26: a[2]= 1./2520.; a[1]= -1./2520.; a[0]= a[3]= 0.; a[4]= -1./630.; break;
      case 27: a[2]= a[1]= 1./8505.; a[0]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 28: a[2]= a[0]= 1./8505.; a[1]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 29: a[2]= 1./2520.; a[3]= -1./2520.; a[0]= a[1]= 0.; a[4]= -1./630.; break;
      case 33: a[3]= 1./1260.; a[0]= a[1]= a[2]= 0.; a[4]= 1./630.; break;
      case 34: a[3]= a[2]= 1./8505.; a[0]= a[1]= -19./68040.; a[4]= -1./486.; break;
      case 35: a[3]= a[1]= 1./8505.; a[0]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 36: a[3]= a[0]= 1./8505.; a[1]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 37: a[3]= 1./2520.; a[0]= -1./2520.; a[1]= a[2]= 0.; a[4]= -1./630.; break;
      case 38: a[3]= 1./2520.; a[1]= -1./2520.; a[0]= a[2]= 0.; a[4]= -1./630.; break;
      case 39: a[3]= 1./2520.; a[2]= -1./2520.; a[0]= a[1]= 0.; a[4]= -1./630.; break;
      case 44: a[0]= a[1]= 37./17010.; a[2]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 45: a[0]= 1./972.; a[1]= a[2]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 46: a[1]= 1./972.; a[0]= a[2]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 47: a[0]= 1./972.; a[1]= a[3]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 48: a[1]= 1./972.; a[0]= a[3]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 49: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 55: a[0]= a[2]= 37./17010.; a[1]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 56: a[2]= 1./972.; a[0]= a[1]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 57: a[0]= 1./972.; a[2]= a[3]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 58: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 59: a[2]= 1./972.; a[0]= a[3]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 66: a[1]= a[2]= 37./17010.; a[0]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 67: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 68: a[1]= 1./972.; a[2]= a[3]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 69: a[2]= 1./972.; a[1]= a[3]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 77: a[0]= a[3]= 37./17010.; a[1]= a[2]= -17./17010.; a[4]= 88./8505.; break;
      case 78: a[3]= 1./972.; a[0]= a[1]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 79: a[3]= 1./972.; a[0]= a[2]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 88: a[1]= a[3]= 37./17010.; a[0]= a[2]= -17./17010.; a[4]= 88./8505.; break;
      case 89: a[3]= 1./972.; a[1]= a[2]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 99: a[2]= a[3]= 37./17010.; a[0]= a[1]= -17./17010.; a[4]= 88./8505.; break;
      default: throw DROPSErrCL("Quad(i,j): no such shape function");
    }
    double sum= a[4]*f(GetBaryCenter(t));
    for(Uint i=0; i<4; ++i)
        sum+= a[i]*f(t.GetVertex(i)->GetCoord());
    return sum;
}


double f0(const Point3DCL&)   { return 1; }
double f1(const Point3DCL& p) { return p[0]; }
double f2(const Point3DCL& p) { return p[1]; }
double f3(const Point3DCL& p) { return p[2]; }
double f4(const Point3DCL& p) { return p[0]*p[0] + p[1]*p[1] + p[2]*p[2]
                                     + p[0]*p[1] + p[0]*p[2] + p[1]*p[2]; }

scalar_coeff_ptr testFct(int i)
{
    switch (i)
    {
    case 0: return &f0;
    case 1: return &f1;
    case 2: return &f2;
    case 3: return &f3;
    case 4: return &f4;
    }
    return 0;
}

void testHat( int i)
{
    double (*H)(double, double, double)=0;
    switch(i)
    {
      case 0: H= &FE_P2CL::H0; break;
      case 1: H= &FE_P2CL::H1; break;
      case 2: H= &FE_P2CL::H2; break;
      case 3: H= &FE_P2CL::H3; break;
      case 4: H= &FE_P2CL::H4; break;
      case 5: H= &FE_P2CL::H5; break;
      case 6: H= &FE_P2CL::H6; break;
      case 7: H= &FE_P2CL::H7; break;
      case 8: H= &FE_P2CL::H8; break;
      case 9: H= &FE_P2CL::H9; break;
    }

    std::cout << "Values on Verts 0 1 2 3: "
              << H(0,0,0)<<", "<<H(1,0,0)<<", "<<H(0,1,0)<<", "<<H(0,0,1)<<std::endl;
    std::cout << "Values on Edges 01 02 03 12 13 23: "
              << H(0.5,0,0)<<", "<<H(0,0.5,0)<<", "<<H(0,0,0.5)<<", "
              << H(0.5,0.5,0)<<", "<<H(0.5,0,0.5)<<", "<<H(0,0.5,0.5)<<", "<<std::endl;
}         

void testGrad( int i)
{
    std::cout << "Values on Verts 0 1 2 3: "
              << FE_P2CL::DHRef(i,0,0,0)<<", "<<FE_P2CL::DHRef(i,1,0,0)<<", "
              << FE_P2CL::DHRef(i,0,1,0)<<", "<<FE_P2CL::DHRef(i,0,0,1)<<std::endl;
}         

void testQuad(const TetraCL& t, int i)
{
    std::cout << "Integrating f*phi_"<<i<<" with f = 1, x, y, z, x^2+y^2+z^2+xy+xz+yz:\n";
    for(int j=0; j<5; ++j)
        std::cout << Quad(t,testFct(j),i) << std::endl;
}

void testQuad2(const TetraCL& t)
{
    for(int k=0; k<5; ++k)
    {
        std::cout << "Checking with Test Fct "<<k<<std::endl;
        std::ostringstream name;  name<<"data"<<k<<".dat";
        std::ofstream fil(name.str().c_str());
        for(int i=0; i<10; ++i)
        {
            for(int j=0; j<10; ++j)
                fil << Quad(t, testFct(k),i,j) << " ";
            fil << std::endl;
        }
        fil.close();
    }
}

int main()
{
    for(int i=0; i<10; ++i)
    {
        std::cout << "Checking Hat Fct "<<i<<std::endl;
        testHat(i);
        std::cout << std::endl;
    }
    for(int i=0; i<10; ++i)
    {
        std::cout << "Checking Gradient of Hat Fct "<<i<<std::endl;
        testGrad(i);
        std::cout << std::endl;
    }
    Point3DCL p0(0.), p1(0.), p2(0.), p3(0.);
    p1[0]= p2[1]= p3[2]= 1;
    VertexCL v0(p0,1), v1(p1,1),
             v2(p2,1), v3(p3,1);
    TetraCL tetra( &v0, &v1, &v2, &v3, 0);
    for(int i=0; i<10; ++i)
    {
        std::cout << "Checking numeric integration of Hat Fct "<<i<<std::endl;
        testQuad(tetra,i);
        std::cout << std::endl;
    }
    std::cout << "Checking numeric integration of two Hat Fcts:"<<std::endl;
    testQuad2(tetra);
    return 0;
}                      
