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

void P2DiscCL::GetGradientsOnRef( Quad2CL<Point3DCL> GRef[10])
{
    for (int i=0; i<10; ++i)
    {
        GRef[i].val[0]= FE_P2CL::DHRef( i, 0,0,0);
        GRef[i].val[1]= FE_P2CL::DHRef( i, 1,0,0);
        GRef[i].val[2]= FE_P2CL::DHRef( i, 0,1,0);
        GRef[i].val[3]= FE_P2CL::DHRef( i, 0,0,1);
        GRef[i].val[4]= FE_P2CL::DHRef( i, 0.25,0.25,0.25);
    }
}

} // end of namespace DROPS
