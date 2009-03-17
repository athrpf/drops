//**************************************************************************
// File:    discretize.cpp                                                 *
// Content: discretizations for several PDEs and FE types                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - October, 2001                                          *
//**************************************************************************

#include "num/discretize.h"
#include "num/fe.h"

namespace DROPS
{
//**************************************************************************
// Class: Quad3DataCL                                                      *
//**************************************************************************

BaryCoordCL Quad3DataCL::Node[NumNodesC];

const double Quad3DataCL::Wght[2]= {
    -2./15., /* -(n+1)^2/[4(n+2)] /6       Node[0]*/
    3./40.,  /* (n+3)^2/[4(n+1)(n+2)] /6 , Node[1] bis Node[4]*/
};

std::valarray<double> Quad3DataCL::P2_Val[10]; // P2_Val[i] contains FE_P2CL::H_i( Node).

Quad3DataCL::Quad3DataCL()
{
    Node[0]= BaryCoordCL( 0.25);
    const double A= 1./6.,
                 B= 0.5;
    Node[1]= MakeBaryCoord( A,A,A,B);
    Node[2]= MakeBaryCoord( A,A,B,A);
    Node[3]= MakeBaryCoord( A,B,A,A);
    Node[4]= MakeBaryCoord( B,A,A,A);

    FE_P2CL::ApplyAll( NumNodesC, Node, P2_Val);
}

BaryCoordCL* Quad3DataCL::TransformNodes (const SArrayCL<BaryCoordCL,4>& M)
{
    BaryCoordCL* tN = new BaryCoordCL[NumNodesC];
    for (Uint i=0; i< NumNodesC; ++i)
        //tN[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
        for (Uint k=0; k<4; ++k)
            tN[i][k]= M[0][k]*Node[i][0] + M[1][k]*Node[i][1]
                    + M[2][k]*Node[i][2] + M[3][k]*Node[i][3];
    return tN;
}

namespace {
    Quad3DataCL theQuad3DataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace

//**************************************************************************
// Class: Quad5DataCL                                                      *
//**************************************************************************

BaryCoordCL Quad5DataCL::Node[NumNodesC];

const double Quad5DataCL::Wght[4]= {
    0.01975308641975308641975308641975308641975, /* 8./405.,                                   Node[0]*/
    0.01198951396316977000173064248499538621936, /* (2665.0 + 14.0*std::sqrt( 15.0))/226800.0, Node[1] bis Node[4]*/
    0.01151136787104539754677023934921978132914, /* (2665.0 - 14.0*std::sqrt( 15.0))/226800.0, Node[5] bis Node[8]*/
    0.008818342151675485008818342151675485008818 /* 5./567.,                                   Node[9] bis Node[14]*/
};

std::valarray<double> Quad5DataCL::P2_Val[10]; // P2_Val[i] contains FE_P2CL::H_i( Node).

Quad5DataCL::Quad5DataCL()
{
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
}

BaryCoordCL* Quad5DataCL::TransformNodes (const SArrayCL<BaryCoordCL,4>& M)
{
    BaryCoordCL* tN = new BaryCoordCL[NumNodesC];
    for (Uint i=0; i< NumNodesC; ++i)
        //tN[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
        for (Uint k=0; k<4; ++k)
            tN[i][k]= M[0][k]*Node[i][0] + M[1][k]*Node[i][1]
                    + M[2][k]*Node[i][2] + M[3][k]*Node[i][3];
    return tN;
}

namespace {
    Quad5DataCL theQuad5DataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace

//**************************************************************************
// Class: Quad5_2DDataCL                                                   *
//**************************************************************************
Point3DCL Quad5_2DDataCL::Node[NumNodesC];  // Barycentric coord for 2D

const double Quad5_2DDataCL::Wght[3]= {
      9./80.,                          /*Node[0]*/
      (155. - std::sqrt( 15.0))/2400., /*Node[1] to Node[3]*/
      (155. + std::sqrt( 15.0))/2400., /*Node[4] to Node[6]*/
};

Quad5_2DDataCL::Quad5_2DDataCL ()
{
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
}

void
Quad5_2DDataCL::SetInterface (const BaryCoordCL*const p, BaryCoordCL* NodeInTetra)
{
    for (Uint i= 0; i < NumNodesC; ++i)
        NodeInTetra[i]= Node[i][0]*p[0] + Node[i][1]*p[1] + Node[i][2]*p[2];
}

namespace {
    Quad5_2DDataCL theQuad52DDataInitializer_; // The constructor sets up the static arrays
} // end of anonymous namespace


const double Quad3PosWeightsCL::_points[8][3]= {
    {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.},
    {1./3.,1./3.,0.}, {1./3.,0.,1./3.}, {0.,1./3.,1./3.},
    {1./3.,1./3.,1./3.}
    };

const double FaceQuad2CL::_points[4][2]= {
    {0., 0.}, {1., 0.}, {0., 1.}, {1./3., 1./3.} };

const double P1BubbleDiscCL::_points1[26][3]= {
    {0.,0.,0.}, {1./3.,1./3.,0.}, {1./3.,0.,1./3.}, {0.,1./3.,1./3.},
    {1.,0.,0.}, {1./3.,1./3.,0.}, {1./3.,0.,1./3.}, {1./3.,1./3.,1./3.},
    {0.,1.,0.}, {1./3.,1./3.,0.}, {0.,1./3.,1./3.}, {1./3.,1./3.,1./3.},
    {0.,0.,1.}, {1./3.,0.,1./3.}, {0.,1./3.,1./3.}, {1./3.,1./3.,1./3.},
    {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}, {.5,0.,0.}, {0.,.5,0.}, {0.,0.,.5}, {.5,.5,0.}, {.5,0.,.5}, {0.,.5,.5}
    };

void P2DiscCL::GetGradientsOnRef( LocalP1CL<Point3DCL> GRef[10])
{
    for (int i= 0; i < 10; ++i)
    {
        GRef[i][0]= FE_P2CL::DHRef( i, 0,0,0);
        GRef[i][1]= FE_P2CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P2CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P2CL::DHRef( i, 0,0,1);
    }
}

void P2DiscCL::GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
    {
        GRef[i][0]= FE_P2CL::DHRef( i, 0,0,0);
        GRef[i][1]= FE_P2CL::DHRef( i, 1,0,0);
        GRef[i][2]= FE_P2CL::DHRef( i, 0,1,0);
        GRef[i][3]= FE_P2CL::DHRef( i, 0,0,1);
        GRef[i][4]= FE_P2CL::DHRef( i, 0.25,0.25,0.25);
    }
}

void P2DiscCL::GetGradientsOnRef( Quad5CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
        for (int j=0; j<Quad5DataCL::NumNodesC; ++j)
        {
            const BaryCoordCL& Node= Quad5DataCL::Node[j];
            GRef[i][j]= FE_P2CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P2DiscCL::GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[10],
    const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<Point3DCL>::SetInterface( p, NodeInTetra);
    for (int i= 0; i < 10; ++i)
        for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
            const BaryCoordCL& Node= NodeInTetra[j];
            GRef[i][j]= FE_P2CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P2DiscCL::GetP2Basis( Quad5_2DCL<> p2[10], const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DDataCL::NumNodesC];
    Quad5_2DCL<>::SetInterface( p, NodeInTetra);
    for (int j= 0; j < Quad5_2DDataCL::NumNodesC; ++j) {
        const BaryCoordCL& Node= NodeInTetra[j];
        p2[0][j]= FE_P2CL::H0( Node);
        p2[1][j]= FE_P2CL::H1( Node);
        p2[2][j]= FE_P2CL::H2( Node);
        p2[3][j]= FE_P2CL::H3( Node);
        p2[4][j]= FE_P2CL::H4( Node);
        p2[5][j]= FE_P2CL::H5( Node);
        p2[6][j]= FE_P2CL::H6( Node);
        p2[7][j]= FE_P2CL::H7( Node);
        p2[8][j]= FE_P2CL::H8( Node);
        p2[9][j]= FE_P2CL::H9( Node);
    }
}

} // end of namespace DROPS
