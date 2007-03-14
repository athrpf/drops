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

const double Quad3CL::_points[8][3]= {
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
        for (int j=0; j<Quad5CL<Point3DCL>::NumNodesC; ++j)
        {
            const BaryCoordCL& Node= Quad5CL<Point3DCL>::Node[j];
            GRef[i][j]= FE_P2CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P2DiscCL::GetGradientsOnRef( Quad5_2DCL<Point3DCL> GRef[10],
    const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DCL<Point3DCL>::NumNodesC];
    Quad5_2DCL<Point3DCL>::SetInterface( p, NodeInTetra);
    for (int i= 0; i < 10; ++i)
        for (int j= 0; j < Quad5_2DCL<Point3DCL>::NumNodesC; ++j) {
            const BaryCoordCL& Node= NodeInTetra[j];
            GRef[i][j]= FE_P2CL::DHRef( i, Node[1], Node[2], Node[3]);
        }
}

void P2DiscCL::GetP2Basis( Quad5_2DCL<> p2[10], const BaryCoordCL* const p)
{
    BaryCoordCL NodeInTetra[Quad5_2DCL<Point3DCL>::NumNodesC];
    Quad5_2DCL<>::SetInterface( p, NodeInTetra);
    for (int j= 0; j < Quad5_2DCL<Point3DCL>::NumNodesC; ++j) {
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
