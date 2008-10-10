//**************************************************************************
// File:    builder.cpp                                                    *
// Content: MGBuilderCL objects for some domains                           *
// Author:  Joerg Peters, Volker Reichelt, Patrick Esser, IGPM RWTH Aachen *
// Version: 0.1                                                            *
// History: begin - October, 3 2000                                        *
//                                                                         *
// Remarks: We should use the const-qualifier to make it difficult to      *
//          accidentally change the multigrid structure from anywhere      *
//          outside of the multigrid algorithms.                           *
//          Thus the pointer to user data structures should probably be    *
//          a pointer to mutable.                                          *
//**************************************************************************

#include <vector>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "geom/builder.h"

namespace DROPS
{


BrickBuilderCL::BrickBuilderCL(const Point3DCL& origin,
                               const Point3DCL& e1,
                               const Point3DCL& e2,
                               const Point3DCL& e3,
                               Uint n1, Uint n2, Uint n3)
    :_orig(origin), _e1(e1), _e2(e2), _e3(e3), _n1(n1), _n2(n2), _n3(n3)
{}

void BrickBuilderCL::buildBoundary (MultiGridCL* mgp) const
{
    BoundaryCL::SegPtrCont& Bnd= GetBnd(mgp);

    Bnd.push_back( new AffineSquareCL(_orig,     _orig    +_e2, _orig    +_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e1, _orig+_e1+_e2, _orig+_e1+_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig,     _orig    +_e1, _orig    +_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e2, _orig+_e2+_e1, _orig+_e2+_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig,     _orig    +_e1, _orig    +_e2) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e3, _orig+_e3+_e1, _orig+_e3+_e2) ); // e1-e2-plane
}

void BrickBuilderCL::build (MultiGridCL* mgp) const
{
    // Check, if the parallelepiped spanned by e1,e2,e3 is degenerated

    AppendLevel(mgp);

    // Create boundary
    buildBoundary(mgp);

    // Create vertices
    MultiGridCL::VertexLevelCont& verts= GetVertices(mgp)[0];
    std::vector<VertexCL*> va( (_n3+1)*(_n2+1)*(_n1+1) );
    const Point3DCL off1= 1.0/static_cast<double>(_n1) * _e1;
    const Point3DCL off2= 1.0/static_cast<double>(_n2) * _e2;
    const Point3DCL off3= 1.0/static_cast<double>(_n3) * _e3;
    Point2DCL e1_2D(0.0);
    Point2DCL e2_2D(0.0);
    e2_2D[1]= e1_2D[0]= 1.0;

    for (Uint i3=0; i3<=_n3; ++i3)
        for (Uint i2=0; i2<=_n2; ++i2)
            for (Uint i1=0; i1<=_n1; ++i1)
            {
                verts.push_back( VertexCL(_orig + static_cast<double>(i1)*off1
                                                + static_cast<double>(i2)*off2
                                                + static_cast<double>(i3)*off3, 0) );
                va[v_idx(i3,i2,i1)]= &verts.back();
                if (i1 == 0)    // y-z-plane
                    verts.back().AddBnd( BndPointCL(0, static_cast<double>(i2)/static_cast<double>(_n2)*e1_2D
                                                      +static_cast<double>(i3)/static_cast<double>(_n3)*e2_2D) );
                if (i1 == _n1)    // y-z-plane
                    verts.back().AddBnd( BndPointCL(1, static_cast<double>(i2)/static_cast<double>(_n2)*e1_2D
                                                      +static_cast<double>(i3)/static_cast<double>(_n3)*e2_2D) );
                if (i2 == 0)    // x-z-plane
                    verts.back().AddBnd( BndPointCL(2, static_cast<double>(i1)/static_cast<double>(_n1)*e1_2D
                                                      +static_cast<double>(i3)/static_cast<double>(_n3)*e2_2D) );
                if (i2 == _n2)    // x-z-plane
                    verts.back().AddBnd( BndPointCL(3, static_cast<double>(i1)/static_cast<double>(_n1)*e1_2D
                                                      +static_cast<double>(i3)/static_cast<double>(_n3)*e2_2D) );
                if (i3 == 0)    // x-y-plane
                    verts.back().AddBnd( BndPointCL(4, static_cast<double>(i1)/static_cast<double>(_n1)*e1_2D
                                                      +static_cast<double>(i2)/static_cast<double>(_n2)*e2_2D) );
                if (i3 == _n3)    // x-y-plane
                    verts.back().AddBnd( BndPointCL(5, static_cast<double>(i1)/static_cast<double>(_n1)*e1_2D
                                                      +static_cast<double>(i2)/static_cast<double>(_n2)*e2_2D) );
                if ( verts.back().IsOnBoundary() ) verts.back().BndSort();
            }

    // Create edges by calling BuildEdges() an BuildFaces() for every new tetrahedron;
    // this will search for all the ones needed and add missing edges automatically;
    // Create tetras
    MultiGridCL::EdgeLevelCont& edges= GetEdges(mgp)[0];
    MultiGridCL::FaceLevelCont& faces= GetFaces(mgp)[0];
    MultiGridCL::TetraLevelCont& tetras= GetTetras(mgp)[0];
    std::vector<TetraCL*> ta(_n3*_n2*_n1*6);

    for (Uint i3=0; i3<_n3; ++i3)
        for (Uint i2=0; i2<_n2; ++i2)
            for (Uint i1=0; i1<_n1; ++i1)
            {   // Add tetrahedrons in one mini-brick; less-than ordering of indices of e_i used
                tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2,i1+1)], va[v_idx(i3,i2+1,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                ta[t_idx(i3,i2,i1,0)]= &tetras.back();
                tetras.back().BuildEdges(edges);
                tetras.back().BuildAndLinkFaces(faces);
                tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2,i1+1)], va[v_idx(i3+1,i2,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                ta[t_idx(i3,i2,i1,1)]= &tetras.back();
                tetras.back().BuildEdges(edges);
                tetras.back().BuildAndLinkFaces(faces);
                tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2+1,i1)], va[v_idx(i3,i2+1,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                ta[t_idx(i3,i2,i1,2)]= &tetras.back();
                tetras.back().BuildEdges(edges);
                tetras.back().BuildAndLinkFaces(faces);
                tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2+1,i1)], va[v_idx(i3+1,i2+1,i1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                ta[t_idx(i3,i2,i1,3)]= &tetras.back();
                tetras.back().BuildEdges(edges);
                tetras.back().BuildAndLinkFaces(faces);
                tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3+1,i2,i1)], va[v_idx(i3+1,i2,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                ta[t_idx(i3,i2,i1,4)]= &tetras.back();
                tetras.back().BuildEdges(edges);
                tetras.back().BuildAndLinkFaces(faces);
                tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3+1,i2,i1)], va[v_idx(i3+1,i2+1,i1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                ta[t_idx(i3,i2,i1,5)]= &tetras.back();
                tetras.back().BuildEdges(edges);
                tetras.back().BuildAndLinkFaces(faces);
            }
    std::for_each( va.begin(), va.end(), std::mem_fun( &VertexCL::DestroyRecycleBin ) );
}


LBuilderCL::LBuilderCL(const Point3DCL& origin,
                       const Point3DCL& e1,
                       const Point3DCL& e2,
                       const Point3DCL& e3,
                       Uint n1, Uint n2, Uint n3,
                       Uint b1, Uint b2)
    :_orig(origin), _e1(e1), _e2(e2), _e3(e3), _n1(n1), _n2(n2), _n3(n3), _b1(b1), _b2(b2)
{}

void LBuilderCL::buildBoundary(MultiGridCL* mgp) const
{
    const double _dn1= static_cast<double>(_n1);
    const double _dn2= static_cast<double>(_n2);
    const double _db1= static_cast<double>(_b1);
    const double _db2= static_cast<double>(_b2);

    const double b1= _db1/_dn1;
    const double b2= _db2/_dn2;

    BoundaryCL::SegPtrCont& Bnd= GetBnd(mgp);

    Bnd.push_back( new AffineSquareCL(_orig, _orig+b2*_e2, _orig+_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2, _orig+_e2, _orig+b2*_e2+_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e1, _orig+_e1+b2*_e2, _orig+_e1+_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b2*_e2, _orig+b1*_e1+_e2, _orig+b1*_e1+b2*_e2+_e3) ); // e2-e3-plane

    Bnd.push_back( new AffineSquareCL(_orig, _orig+b1*_e1, _orig+_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1, _orig+_e1, _orig+b1*_e1+_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e2, _orig+b1*_e1+_e2, _orig+_e2+_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b2*_e2, _orig+_e1+b2*_e2, _orig+b1*_e1+b2*_e2+_e3) ); // e1-e3-plane

    Bnd.push_back( new AffineSquareCL(_orig, _orig+b1*_e1, _orig+b2*_e2) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2, _orig+b1*_e1+b2*_e2, _orig+_e2) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1, _orig+_e1, _orig+b1*_e1+b2*_e2) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e3, _orig+b1*_e1+_e3, _orig+b2*_e2+_e3) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2+_e3, _orig+b1*_e1+b2*_e2+_e3, _orig+_e2+_e3) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+_e3, _orig+_e1+_e3, _orig+b1*_e1+b2*_e2+_e3) ); // e1-e2-plane
}

void LBuilderCL::build (MultiGridCL* mgp) const
{
    const double _dn1= static_cast<double>(_n1);
    const double _dn2= static_cast<double>(_n2);
    const double _dn3= static_cast<double>(_n3);
    const double _db1= static_cast<double>(_b1);
    const double _db2= static_cast<double>(_b2);

    // Check, if the parallelepiped spanned by e1,e2,e3 is degenerated

    AppendLevel(mgp);

    // Create boundary
    buildBoundary(mgp);

    // Create vertices
    MultiGridCL::VertexLevelCont& verts= GetVertices(mgp)[0];
    std::vector<VertexCL*> va( (_n3+1)*(_n2+1)*(_n1+1) );
    const Point3DCL off1= 1.0/_dn1 * _e1;
    const Point3DCL off2= 1.0/_dn2 * _e2;
    const Point3DCL off3= 1.0/_dn3 * _e3;
    Point2DCL e1_2D(0.0);
    Point2DCL e2_2D(0.0);
    e2_2D[1]= e1_2D[0]= 1.0;

    for (Uint i3=0; i3<=_n3; ++i3)
        for (Uint i2=0; i2<=_n2; ++i2)
            for (Uint i1=0; i1<=_n1; ++i1)
                if (i1<=_b1 || (i1>_b1 && i2<=_b2)) // stay in the L-shaped form
                {
                    verts.push_back( VertexCL(_orig + static_cast<double>(i1)*off1
                                                    + static_cast<double>(i2)*off2
                                                    + static_cast<double>(i3)*off3,
                                              static_cast<Uint>(0)) );
                    va[v_idx(i3,i2,i1)]= &verts.back();
                    if (i1 == 0 && i2 <= _b2)    // y-z-plane
                        verts.back().AddBnd( BndPointCL(0, static_cast<double>(i2)/_db2*e1_2D
                                                          +static_cast<double>(i3)/_dn3*e2_2D) );
                    if (i1 == 0 && i2 >= _b2)    // y-z-plane
                        verts.back().AddBnd( BndPointCL(1, static_cast<double>(i2-_b2)/(_dn2-_db2)*e1_2D
                                                          +static_cast<double>(i3)    /_dn3       *e2_2D) );
                    if (i1 == _n1 && i2<=_b2)    // y-z-plane
                        verts.back().AddBnd( BndPointCL(2, static_cast<double>(i2)/_db2*e1_2D
                                                          +static_cast<double>(i3)/_dn3*e2_2D) );
                    if (i1 == _b1 && i2>=_b2)    // y-z-plane
                        verts.back().AddBnd( BndPointCL(3, static_cast<double>(i2-_b2)/(_dn2-_db2)*e1_2D
                                                          +static_cast<double>(i3)    /_dn3       *e2_2D) );
                    if (i2 == 0 && i1 <= _b1)    // x-z-plane
                        verts.back().AddBnd( BndPointCL(4, static_cast<double>(i1)/_db1*e1_2D
                                                          +static_cast<double>(i3)/_dn3*e2_2D) );
                    if (i2 == 0 && i1 >= _b1)    // x-z-plane
                        verts.back().AddBnd( BndPointCL(5, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                          +static_cast<double>(i3)    /_dn3       *e2_2D) );
                    if (i2 == _n2 && i1<=_b1)    // x-z-plane
                        verts.back().AddBnd( BndPointCL(6, static_cast<double>(i1)/_db1*e1_2D
                                                          +static_cast<double>(i3)/_dn3*e2_2D) );
                    if (i2 == _b2 && i1>=_b1)    // x-z-plane
                        verts.back().AddBnd( BndPointCL(7, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                          +static_cast<double>(i3)    /_dn3       *e2_2D) );
                    if (i3 == 0 && i1<=_b1 && i2<=_b2)    // x-y-plane
                        verts.back().AddBnd( BndPointCL(8, static_cast<double>(i1)/_db1*e1_2D
                                                          +static_cast<double>(i2)/_db2*e2_2D) );
                    if (i3 == 0 && i1<=_b1 && i2>=_b2)    // x-y-plane
                        verts.back().AddBnd( BndPointCL(9, static_cast<double>(i1)    /_db1       *e1_2D
                                                          +static_cast<double>(i2-_b2)/(_dn2-_db2)*e2_2D) );
                    if (i3 == 0 && i1>=_b1 && i2<=_b2)    // x-y-plane
                        verts.back().AddBnd( BndPointCL(10, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i2)    /_db2       *e2_2D) );
                    if (i3 == _n3 && i1<=_b1 && i2<=_b2)    // x-y-plane
                        verts.back().AddBnd( BndPointCL(11, static_cast<double>(i1)/_db1*e1_2D
                                                           +static_cast<double>(i2)/_db2*e2_2D) );
                    if (i3 == _n3 && i1<=_b1 && i2>=_b2)    // x-y-plane
                        verts.back().AddBnd( BndPointCL(12, static_cast<double>(i1)    /_db1       *e1_2D
                                                           +static_cast<double>(i2-_b2)/(_dn2-_db2)*e2_2D) );
                    if (i3 == _n3 && i1>=_b1 && i2<=_b2)    // x-y-plane
                        verts.back().AddBnd( BndPointCL(13, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i2)    /_db2       *e2_2D) );
                    if ( verts.back().IsOnBoundary() ) verts.back().BndSort();
                }

    // Create edges/faces by calling BuildEdges() and BuildAndLinkFaces(() for every new tetrahedron;
    // this will search for all the ones needed and add missing edges automatically;
    // Create tetras
    MultiGridCL::EdgeLevelCont& edges= GetEdges(mgp)[0];
    MultiGridCL::FaceLevelCont& faces= GetFaces(mgp)[0];
    MultiGridCL::TetraLevelCont& tetras= GetTetras(mgp)[0];
    std::vector<TetraCL*> ta(_n3*_n2*_n1*6);

    for (Uint i3=0; i3<_n3; ++i3)
        for (Uint i2=0; i2<_n2; ++i2)
            for (Uint i1=0; i1<_n1; ++i1)
            // Add tetrahedrons in one mini-brick; less-than ordering of indices of e_i used
                if (i1<_b1 || (i1>=_b1 && i2<_b2))
                {
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2,i1+1)], va[v_idx(i3,i2+1,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,0)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2,i1+1)], va[v_idx(i3+1,i2,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,1)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2+1,i1)], va[v_idx(i3,i2+1,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,2)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2+1,i1)], va[v_idx(i3+1,i2+1,i1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,3)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3+1,i2,i1)], va[v_idx(i3+1,i2,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,4)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3+1,i2,i1)], va[v_idx(i3+1,i2+1,i1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,5)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                }
    std::for_each( va.begin(), va.end(), std::mem_fun( &VertexCL::DestroyRecycleBin ) );
}


BBuilderCL::BBuilderCL(const Point3DCL& origin,
                       const Point3DCL& e1,
                       const Point3DCL& e2,
                       const Point3DCL& e3,
                       Uint n1, Uint n2, Uint n3,
                       Uint b1, Uint b2, Uint b3)
    :_orig(origin), _e1(e1), _e2(e2), _e3(e3), _n1(n1), _n2(n2), _n3(n3), _b1(b1), _b2(b2), _b3(b3)
{}

void BBuilderCL::buildBoundary(MultiGridCL* mgp) const
{
    const double _dn1= static_cast<double>(_n1);
    const double _dn2= static_cast<double>(_n2);
    const double _dn3= static_cast<double>(_n3);
    const double _db1= static_cast<double>(_b1);
    const double _db2= static_cast<double>(_b2);
    const double _db3= static_cast<double>(_b3);

    const double b1= _db1/_dn1;
    const double b2= _db2/_dn2;
    const double b3= _db3/_dn3;

    BoundaryCL::SegPtrCont& Bnd= GetBnd(mgp);
    // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig, _orig+b2*_e2, _orig+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2, _orig+_e2, _orig+b2*_e2+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b3*_e3, _orig+b2*_e2+b3*_e3, _orig+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2+b3*_e3, _orig+_e2+b3*_e3, _orig+b2*_e2+_e3) );

    Bnd.push_back( new AffineSquareCL(_orig+_e1, _orig+_e1+b2*_e2, _orig+_e1+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+_e1+b2*_e2, _orig+_e1+_e2, _orig+_e1+b2*_e2+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+_e1+b3*_e3, _orig+_e1+b2*_e2+b3*_e3, _orig+_e1+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b2*_e2+b3*_e3, _orig+b1*_e1+_e2+b3*_e3, _orig+b1*_e1+b2*_e2+_e3) );
    // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig, _orig+b1*_e1, _orig+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1, _orig+_e1, _orig+b1*_e1+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b3*_e3, _orig+b1*_e1+b3*_e3, _orig+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b3*_e3, _orig+_e1+b3*_e3, _orig+b1*_e1+_e3) );

    Bnd.push_back( new AffineSquareCL(_orig+_e2, _orig+b1*_e1+_e2, _orig+_e2+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+_e2, _orig+_e1+_e2, _orig+b1*_e1+_e2+b3*_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+_e2+b3*_e3, _orig+b1*_e1+_e2+b3*_e3, _orig+_e2+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b2*_e2+b3*_e3, _orig+_e1+b2*_e2+b3*_e3, _orig+b1*_e1+b2*_e2+_e3) );
    // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig, _orig+b1*_e1, _orig+b2*_e2) );
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2, _orig+b1*_e1+b2*_e2, _orig+_e2) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1, _orig+_e1, _orig+b1*_e1+b2*_e2) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b2*_e2, _orig+_e1+b2*_e2, _orig+b1*_e1+_e2) );

    Bnd.push_back( new AffineSquareCL(_orig+_e3, _orig+b1*_e1+_e3, _orig+b2*_e2+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b2*_e2+_e3, _orig+b1*_e1+b2*_e2+_e3, _orig+_e2+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+_e3, _orig+_e1+_e3, _orig+b1*_e1+b2*_e2+_e3) );
    Bnd.push_back( new AffineSquareCL(_orig+b1*_e1+b2*_e2+b3*_e3, _orig+_e1+b2*_e2+b3*_e3, _orig+b1*_e1+_e2+b3*_e3) );
    
}


void BBuilderCL::build (MultiGridCL* mgp) const
{
    const double _dn1= static_cast<double>(_n1);
    const double _dn2= static_cast<double>(_n2);
    const double _dn3= static_cast<double>(_n3);
    const double _db1= static_cast<double>(_b1);
    const double _db2= static_cast<double>(_b2);
    const double _db3= static_cast<double>(_b3);

//    const double b1= _db1/_dn1;
//    const double b2= _db2/_dn2;
//    const double b3= _db3/_dn3;

    // Check, if the parallelepiped spanned by e1,e2,e3 is degenerated

    AppendLevel(mgp);

    // Create boundary
    buildBoundary(mgp);
    
    // Create vertices
    MultiGridCL::VertexLevelCont& verts= GetVertices(mgp)[0];
    std::vector<VertexCL*> va( (_n3+1)*(_n2+1)*(_n1+1) );
    const Point3DCL off1= 1.0/_dn1 * _e1;
    const Point3DCL off2= 1.0/_dn2 * _e2;
    const Point3DCL off3= 1.0/_dn3 * _e3;
    Point2DCL e1_2D(0.0);
    Point2DCL e2_2D(0.0);
    e2_2D[1]= e1_2D[0]= 1.0;
    for (Uint i3=0; i3<=_n3; ++i3)
        for (Uint i2=0; i2<=_n2; ++i2)
            for (Uint i1=0; i1<=_n1; ++i1)
                if (i3<=_b3 || i2<=_b2 || i1<=_b1 ) // stay in the domain
                {
                    verts.push_back( VertexCL(_orig + static_cast<double>(i1)*off1
                                                    + static_cast<double>(i2)*off2
                                                    + static_cast<double>(i3)*off3,
                                              static_cast<Uint>(0)) );
                    va[v_idx(i3,i2,i1)]= &verts.back();

                    // y-z-plane
                    if (i1 == 0 && i2 <= _b2 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(0, static_cast<double>(i2)/_db2*e1_2D
                                                          +static_cast<double>(i3)/_db3*e2_2D) );
                    if (i1 == 0 && i2 >= _b2 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(1, static_cast<double>(i2-_b2)/(_dn2-_db2)*e1_2D
                                                          +static_cast<double>(i3)    /_db3       *e2_2D) );
                    if (i1 == 0 && i2 <= _b2 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(2, static_cast<double>(i2)    /_db2       *e1_2D
                                                          +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );
                    if (i1 == 0 && i2 >= _b2 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(3, static_cast<double>(i2-_b2)/(_dn2-_db2)*e1_2D
                                                          +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );
                    if (i1 == _n1 && i2 <= _b2 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(4, static_cast<double>(i2)/_db2*e1_2D
                                                          +static_cast<double>(i3)/_db3*e2_2D) );
                    if (i1 == _n1 && i2 >= _b2 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(5, static_cast<double>(i2-_b2)/(_dn2-_db2)*e1_2D
                                                          +static_cast<double>(i3)    /_db3       *e2_2D) );
                    if (i1 == _n1 && i2 <= _b2 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(6, static_cast<double>(i2)    /_db2       *e1_2D
                                                          +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );
                    if (i1 == _b1 && i2 >= _b2 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(7, static_cast<double>(i2-_b2)/(_dn2-_db2)*e1_2D
                                                          +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );

                    // x-z-plane
                    if (i2 == 0 && i1 <= _b1 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(8, static_cast<double>(i1)/_db1*e1_2D
                                                          +static_cast<double>(i3)/_db3*e2_2D) );
                    if (i2 == 0 && i1 >= _b1 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(9, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                          +static_cast<double>(i3)    /_db3       *e2_2D) );
                    if (i2 == 0 && i1 <= _b1 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(10, static_cast<double>(i1)    /_db1       *e1_2D
                                                           +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );
                    if (i2 == 0 && i1 >= _b1 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(11, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );
                    if (i2 == _n2 && i1 <= _b1 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(12, static_cast<double>(i1)/_db1*e1_2D
                                                           +static_cast<double>(i3)/_db3*e2_2D) );
                    if (i2 == _n2 && i1 >= _b1 && i3<=_b3)
                        verts.back().AddBnd( BndPointCL(13, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i3)    /_db3       *e2_2D) );
                    if (i2 == _n2 && i1 <= _b1 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(14, static_cast<double>(i1)    /_db1       *e1_2D
                                                           +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );
                    if (i2 == _b2 && i1 >= _b1 && i3>=_b3)
                        verts.back().AddBnd( BndPointCL(15, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i3-_b3)/(_dn3-_db3)*e2_2D) );

                    // x-y-plane
                    if (i3 == 0 && i1<=_b1 && i2<=_b2)
                        verts.back().AddBnd( BndPointCL(16, static_cast<double>(i1)/_db1*e1_2D
                                                           +static_cast<double>(i2)/_db2*e2_2D) );
                    if (i3 == 0 && i1<=_b1 && i2>=_b2)
                        verts.back().AddBnd( BndPointCL(17, static_cast<double>(i1)    /_db1       *e1_2D
                                                           +static_cast<double>(i2-_b2)/(_dn2-_db2)*e2_2D) );
                    if (i3 == 0 && i1>=_b1 && i2<=_b2)
                        verts.back().AddBnd( BndPointCL(18, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i2)    /_db2       *e2_2D) );
                    if (i3 == 0 && i1>=_b1 && i2>=_b2)
                        verts.back().AddBnd( BndPointCL(19, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i2-_b2)/(_dn2-_db2)*e2_2D) );
                    if (i3 == _n3 && i1<=_b1 && i2<=_b2)
                        verts.back().AddBnd( BndPointCL(20, static_cast<double>(i1)/_db1*e1_2D
                                                           +static_cast<double>(i2)/_db2*e2_2D) );
                    if (i3 == _n3 && i1<=_b1 && i2>=_b2)
                        verts.back().AddBnd( BndPointCL(21, static_cast<double>(i1)    /_db1       *e1_2D
                                                           +static_cast<double>(i2-_b2)/(_dn2-_db2)*e2_2D) );
                    if (i3 == _n3 && i1>=_b1 && i2<=_b2)
                        verts.back().AddBnd( BndPointCL(22, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i2)    /_db2       *e2_2D) );
                    if (i3 == _b3 && i1>=_b1 && i2>=_b2)
                        verts.back().AddBnd( BndPointCL(23, static_cast<double>(i1-_b1)/(_dn1-_db1)*e1_2D
                                                           +static_cast<double>(i2-_b2)/(_dn2-_db2)*e2_2D) );
                    if ( verts.back().IsOnBoundary() ) verts.back().BndSort();
                }

    // Create edges by calling BuildEdges() for every new tetrahedron; this will search for all the ones
    // needed and add missing edges automatically;
    // Create tetras
    MultiGridCL::EdgeLevelCont& edges= GetEdges(mgp)[0];
    MultiGridCL::FaceLevelCont& faces= GetFaces(mgp)[0];
    MultiGridCL::TetraLevelCont& tetras= GetTetras(mgp)[0];
    std::vector<TetraCL*> ta(_n3*_n2*_n1*6);
    for (Uint i3=0; i3<_n3; ++i3)
        for (Uint i2=0; i2<_n2; ++i2)
            for (Uint i1=0; i1<_n1; ++i1)
            // Add tetrahedrons in one mini-brick; less-than ordering of indices of e_i used
                if (i1<_b1 || i2<_b2 || i3<_b3 )
                {
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2,i1+1)], va[v_idx(i3,i2+1,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,0)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2,i1+1)], va[v_idx(i3+1,i2,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,1)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2+1,i1)], va[v_idx(i3,i2+1,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,2)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3,i2+1,i1)], va[v_idx(i3+1,i2+1,i1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,3)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3+1,i2,i1)], va[v_idx(i3+1,i2,i1+1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,4)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                    tetras.push_back( TetraCL(va[v_idx(i3,i2,i1)], va[v_idx(i3+1,i2,i1)], va[v_idx(i3+1,i2+1,i1)], va[v_idx(i3+1,i2+1,i1+1)], 0) );
                    ta[t_idx(i3,i2,i1,5)]= &tetras.back();
                    tetras.back().BuildEdges(edges);
                    tetras.back().BuildAndLinkFaces(faces);
                }
    std::for_each( verts.begin(), verts.end(), std::mem_fun_ref( &VertexCL::DestroyRecycleBin ) );
}


TetraBuilderCL::TetraBuilderCL(Ubyte rule)
    :rule_(rule), p0_( std_basis<3>(0)), p1_( std_basis<3>(1)),
                  p2_( std_basis<3>(2)), p3_( std_basis<3>(3))
{}

TetraBuilderCL::TetraBuilderCL(Ubyte rule, const Point3DCL& p0, const Point3DCL& p1,
                                           const Point3DCL& p2, const Point3DCL& p3)
    :rule_(rule), p0_( p0), p1_( p1), p2_( p2), p3_( p3)
{}


// This is absolutely evil.
// One can force the refinement algorithm to go from the current rule to rule.
// This function only works for multigrids built by the current class.
void TetraBuilderCL::BogoReMark(DROPS::MultiGridCL& mg, DROPS::Uint rule)
{
    TetraCL& t= *mg.GetTetrasBegin( 0);

    if (rule == t.GetRefRule()) return;
    switch (rule) {
      case NoRefMarkC:
        if (t.IsRegularlyRef()) {
            for (Uint i=0; i<8; ++i)
                t.GetChild( i)->SetRemoveMark();
        }
        else {
            for (Uint i=0; i<6; ++i)
                if (t.GetRefRule() & (1<<i))
                    const_cast<EdgeCL*>( t.GetEdge( i))->DecMarkForRef();
        }
        break;
      case RegRefMarkC:
        if (t.IsUnrefined()) t.SetRegRefMark();
        else {
            for (Uint i=0; i<6; ++i)
                if (t.GetRefRule() & (1<<i))
                    const_cast<EdgeCL*>( t.GetEdge( i))->DecMarkForRef();
            (*t.GetChildBegin())->SetRegRefMark();
        }
        break;
      case GreenRegRefMarkC:
        if (t.IsRegularlyRef()) {
            for (Uint i=0; i<8; ++i)
                t.GetChild( i)->SetRemoveMark();
            for (Uint i=0; i<6; ++i)
                const_cast<EdgeCL*>( t.GetEdge( i))->IncMarkForRef();
        }
        else if (t.IsUnrefined())
                for (Uint i=0; i<6; ++i)
                    const_cast<EdgeCL*>( t.GetEdge( i))->IncMarkForRef();
            else
                for (Uint i=0; i<6; ++i)
                    if ( !(t.GetRefRule() & (1<<i)))
                        const_cast<EdgeCL*>( t.GetEdge( i))->IncMarkForRef();
        break;
      default:
        if (t.IsRegularlyRef()) {
            for (Uint i=0; i<8; ++i)
                t.GetChild( i)->SetRemoveMark();
        }
        else {
            for (Uint i=0; i<6; ++i)
                if (t.GetRefRule() & (1<<i))
                    const_cast<EdgeCL*>( t.GetEdge( i))->DecMarkForRef();
        }
        for (Uint i=0; i<6; ++i)
            if (rule & (1<<i))
                const_cast<EdgeCL*>( t.GetEdge( i))->IncMarkForRef();
    }
    mg.Refine();
}

void TetraBuilderCL::buildBoundary(MultiGridCL* mgp) const
{
    BoundaryCL::SegPtrCont& Bnd= GetBnd( mgp);
    Bnd.push_back( new AffineTriangleCL( p0_, p1_, p2_)); // e1-e2-plane
    Bnd.push_back( new AffineTriangleCL( p0_ ,p1_, p3_)); // e1-e3-plane
    Bnd.push_back( new AffineTriangleCL( p0_, p2_, p3_)); // e1-e2-plane
    Bnd.push_back( new AffineTriangleCL( p1_, p2_, p3_)); // lid
}

void TetraBuilderCL::build(MultiGridCL* mgp) const
{
    AppendLevel( mgp);

    // Create boundary
    buildBoundary(mgp);
    
    // Create vertices
    MultiGridCL::VertexLevelCont& verts= GetVertices(mgp)[0];
    std::vector<VertexCL*> va( 4);
    // origin
    verts.push_back( VertexCL( p0_, 0));
    va[0]= &verts.back();
    verts.back().AddBnd( BndPointCL( 0, std_basis<2>( 0)));
    verts.back().AddBnd( BndPointCL( 1, std_basis<2>( 0)));
    verts.back().AddBnd( BndPointCL( 2, std_basis<2>( 0)));
    verts.back().BndSort();
    // e1
    verts.push_back( VertexCL( p1_, 0));
    va[1]= &verts.back();
    verts.back().AddBnd( BndPointCL( 0, std_basis<2>( 1)));
    verts.back().AddBnd( BndPointCL( 1, std_basis<2>( 1)));
    verts.back().AddBnd( BndPointCL( 3, std_basis<2>( 0)));
    verts.back().BndSort();
    // e2
    verts.push_back( VertexCL( p2_, 0));
    va[2]= &verts.back();
    verts.back().AddBnd( BndPointCL( 0, std_basis<2>( 2)));
    verts.back().AddBnd( BndPointCL( 2, std_basis<2>( 1)));
    verts.back().AddBnd( BndPointCL( 3, std_basis<2>( 1)));
    verts.back().BndSort();
    // e3
    verts.push_back( VertexCL( p3_, 0));
    va[3]= &verts.back();
    verts.back().AddBnd( BndPointCL( 1, std_basis<2>( 2)));
    verts.back().AddBnd( BndPointCL( 2, std_basis<2>( 2)));
    verts.back().AddBnd( BndPointCL( 3, std_basis<2>( 2)));
    verts.back().BndSort();

    // Create edges by calling BuildEdges() an BuildFaces() for every new tetrahedron;
    // this will search for all the ones needed and add missing edges automatically;
    // Create tetras
    MultiGridCL::EdgeLevelCont& edges= GetEdges(mgp)[0];
    MultiGridCL::FaceLevelCont& faces= GetFaces(mgp)[0];
    MultiGridCL::TetraLevelCont& tetras= GetTetras(mgp)[0];
    TetraCL* tp;

    // Add tetrahedron
    tetras.push_back( TetraCL( va[0], va[1], va[2], va[3], 0));
    tp= &tetras.back();
    tetras.back().BuildEdges( edges);
    tetras.back().BuildAndLinkFaces( faces);

    // Clean up recycle bins, that are used by BuildEdges and BuildAndLinkFaces.
    std::for_each( va.begin(), va.end(), std::mem_fun( &VertexCL::DestroyRecycleBin));

    // The preceeding part was routine. Now, we artificially mark the edges
    // of the desired refinement rule *regularly* and refine once.
    if ( rule_ == RegRefMarkC) tp->SetRegRefMark();
    else
        for (Uint i= 0; i < 6; ++i)
            if (rule_ & (1<<i)) const_cast<EdgeCL*>(tp->GetEdge( i))->IncMarkForRef();
    // Builder are called between PrepareModify and FinalizeModify.
    // PrepareModify and FinalizeModify may only be called in turns and building
    // is finished anyway.
    FinalizeModify( mgp);
    mgp->Refine();
    // Needed due to FinalizeModify() in MultiGridCL-constructor.
    PrepareModify( mgp);
}


//--------------------------------------------------------------------
// Mesh-file-parser
//--------------------------------------------------------------------

MeshStringCL&
operator>>(std::istream& is, MeshStringCL& m)
{
    m.isp_= &is;
    return m;
}

std::istream&
MeshStringCL::operator>>(std::string& s)
{
    Assert( isp_!= 0, "std::istream& MeshStringCL::operator>>(std::string&): Which stream?", ~0u);
    char c;
    *isp_ >> c;
    if (c != '\"')
        throw DROPSErrCL( "std::istream& MeshStringCL::operator>>(std::string&): Not a mesh string.");
    std::getline( *isp_, s, '\"'); // Reads the closing ",
                                   // but does not enter it into s.
                                   // Embedded \" (quoted ") kill us.
    return *isp_;
}

std::istream&
MeshStringCL::SkipMeshString(std::istream& is, bool skipleadingquote)
{
    if (skipleadingquote) {
        char c;
        is >> c;
    }
    std::string s;
    std::getline( is, s, '\"');
    return is;
}

void
MeshNodeCL::Check(std::ostream* msg)
{
    if (section.size() == 1) {
        if (section[0].headerinfo[1] != 1)
            throw DROPSErrCL( "MeshNodeCL::Check(): Node-indices must begin with 1.\n");
        if (num_expected != section[0].point.size()) {
            num_expected= section[0].point.size();
            if (msg)
                *msg << "MeshNodeCL::Check(): Unexpected number of nodes. Corrected.\n";
        }
        return;
    }
    else
        throw DROPSErrCL("MeshNodeCL::Check(): Multiple node sections are currently not supported.\n");
}

void
MeshFaceCL::Check(std::ostream* msg)
{
    std::sort( section.begin(), section.end(), FirstIndexLessCL());
    Uint num= section[0].mface.size();
    if (section[0].headerinfo[1] != 1)
        throw DROPSErrCL( "MeshFaceCL::Check(): Face-indices must begin with 1.\n");
    for (Uint i= 0; i<section.size()-1; ++i) {
        num+= section[i+1].mface.size();
        if (section[i].headerinfo[2] + 1 != section[i+1].headerinfo[1]) {
            throw DROPSErrCL( "MeshFaceCL::Check(): Error: The enumeration of faces has holes.\n");
        }
    }
    if (num != num_expected && msg)
        *msg << "MeshFaceCL::Check(): Unexpected number of faces. Corrected.\n";
    num_expected= num;
}

MFaceCL
MeshFaceCL::operator[](Uint i) const
{
    Assert( i<=num_expected, "MeshFaceCL::operator[](Uint): Index out of bounds.\n", DebugRefineEasyC);
    Uint s= 0;
    for (; i>section[s].headerinfo[2]; ++s) ;
    return section[s].mface[i-section[s].headerinfo[1]];
}

void
MeshCellCL::Check(std::ostream*)
{
    if (section.size() == 1) {
        return;
    }
    else
        throw DROPSErrCL("MeshCellCL::Check(): Multiple cell sections are currently not supported.\n");
}

std::vector<Uint>
HybridCellCL::Vertices()
{
    std::vector<Uint> v;
    for (Uint i= 0; i<3; ++i) {
        v.push_back( mf[0][i]);
        v.push_back( mf[1][i]);
    }
    std::sort( v.begin(), v.end());
    v.erase( std::unique( v.begin(), v.end()), v.end());
    if (v.size() != 4)
        DROPSErrCL( "HybridCellCL::Vertices(): Cell does not have 4 nodes.\n");
    return v;
}

std::vector<FaceCL*>
HybridCellCL::Faces()
{
    return fp;
}

FaceCL*
HybridCellCL::Face(Uint i, Uint j, Uint k)
{
    for (Uint l= 0; l<4; ++l) {
        std::vector<Uint> face( mf[l].begin(), mf[l].begin() + 3);
        if (is_in( face.begin(), face.end(), i)
            &&  is_in( face.begin(), face.end(), j)
            &&  is_in( face.begin(), face.end(), k))
            return fp[l];
    }
    throw( "HybridCellCL::Face(Uint, Uint, Uint): Face not found.\n");
}

FaceCL*
HybridCellCL::face(Uint i)
{
    return fp[i];
}

void
HybridCellCL::Check()
{
    if (fp.size() != 4)
        throw DROPSErrCL( "HybridCellCL::Check(): Cell does not have 4 faces.\n");
    std::vector<MFaceCL> mf_= mf;
    std::sort( mf_.begin(), mf_.end());
    if ( mf_.end() != std::unique( mf_.begin(), mf_.end()))
        throw DROPSErrCL( "HybridCellCL::Check(): Cell-faces coincide.\n");
}

const char *ReadMeshBuilderCL::SymbolicName_[]= {
    "XF_COMMENT",
    "XF_HEADER",
    "XF_DIMENSION",
    "XF_NODE",
    "XF_CELL",
    "XF_FACE"
    };

bool
ReadMeshBuilderCL::NextSection() const
{
    bool done= false;
    do {
        switch(f_.peek()) {
          case '\"':
            MeshStringCL::SkipMeshString( f_);
            break;
          case '(':
            done= true;
            break;
          default:
            f_.get();
            break;
        }
    } while(f_ && !done);
    return f_ && done;
}

bool
ReadMeshBuilderCL::SkipSection(bool skipleadingparenthesis) const
{
    if (skipleadingparenthesis) { // Also, read whitespace until (.
        char c;
        f_ >> c;
    }
    Uint nesting_level= 0;
    bool done= false;
    do {
        switch(f_.peek()) {
          case '\"':
            MeshStringCL::SkipMeshString( f_);
            break;
          case '(':
            ++nesting_level;
            f_.get();
            break;
          case ')':
            if (nesting_level == 0) done= true;
            else --nesting_level;
            f_.get();
            break;
          default:
            f_.get();
            break;
        }
    } while(f_ && !done);
    return f_ && done;
}

Uint
ReadMeshBuilderCL::ReadId()
{
    char c= f_.get(); // Read (.
    if (msg_ && c != '(')
        *msg_ << "Uint ReadMeshBuilderCL::ReadId(): Not at the start of a section.\n";
    Uint id;
    f_ >> id;
    id_history_.push_back( id);
    return id;
}

HeaderInfoCL
ReadMeshBuilderCL::ReadHeaderInfoHex()
{
    std::vector<Uint> ret;
    ret.reserve( 5); // Most Headers have 5 entries;
    char c;
    f_ >> c; // Eat Whitespace and (.
    f_ >> std::hex;
    std::copy( std::istream_iterator<Uint>( f_), std::istream_iterator<Uint>(),
               std::back_inserter( ret));

    f_.clear();
    f_ >> std::dec; // TODO: Preserve Stream-State by copying.
    f_ >> c; // Eat whitespace and ).
    return ret;
}

void
ReadMeshBuilderCL::ReadNode()
{
    if (msg_)
       *msg_ << "XF_NODE: " << std::endl;
    HeaderInfoCL hi= ReadHeaderInfoHex();
    if (hi.empty())
        throw DROPSErrCL( "ReadMeshBuilderCL::ReadNode(): Empty Header.\n");
    switch (hi[0]) {
      case 0: // Declare total number of nodes, no vertices defined.
        nodes_.num_expected= hi[2];
        if (hi[1]!= 1 && msg_)
            *msg_ << "ReadMeshBuilderCL::ReadNode(): Header out of spec. Ignored.\n";
        SkipSection( false);
        break;
      default:
        NextSection();
        char ch; f_ >> ch; // Eat (;
        nodes_.section.push_back( NodeSectionCL());
        NodeSectionCL& ns= nodes_.section.back();
        ns.headerinfo= hi;
        ns.point.reserve( hi[2] - hi[1] + 1);
        Point3DCL p;
        while(true) {
            f_ >> p[0] >> p[1] >> p[2];
            if (f_)
                ns.point.push_back( p);
            else {
                f_.clear();
                break;
            }
        }
        if (ns.point.size() != hi[2] - hi[1] + 1) {
            if (msg_)
                *msg_ << "ReadMeshBuilderCL::ReadNode(): Wrong number of nodes. Adjusting last-index for this section.\n";
            ns.headerinfo[2]= ns.headerinfo[1] + ns.point.size();
        }
        f_ >> ch; // Eat );
        SkipSection( false);
        break;
    }
}

void
ReadMeshBuilderCL::ReadFace()
{
    if (msg_)
       *msg_ << "XF_FACE: " << std::endl;
    HeaderInfoCL hi= ReadHeaderInfoHex();
    if (hi.empty())
        throw DROPSErrCL( "ReadMeshBuilderCL::ReadFace(): Empty Header.\n");
    switch (hi[0]) {
      case 0: // Declare total number of faces, no cells defined.
        mfaces_.num_expected= hi[2];
        SkipSection( false);
        break;
      default:
        if (hi[4] != 3)
            throw DROPSErrCL("ReadMeshBuilderCL::ReadFace(): Error: Faces are not triangular.\n");
        NextSection();
        char ch; f_ >> ch; // Eat (;
        mfaces_.section.push_back( MFaceSectionCL());
        MFaceSectionCL& s= mfaces_.section.back();
        s.headerinfo= hi;
        s.mface.reserve( hi[2] - hi[1] + 1);
        MFaceCL mf;
        f_ >> std::hex;
        while(true) {
            f_ >> mf[0] >> mf[1] >> mf[2] >> mf[3] >> mf[4];
            if (f_)
                s.mface.push_back( mf);
            else {
                f_.clear();
                f_ >> std::dec; // TODO: Preserve Stream-State by copying.
                break;
            }
        }
        if (s.mface.size() != hi[2] - hi[1] + 1) {
            if (msg_)
                *msg_ << "ReadMeshBuilderCL::ReadFace(): Wrong number of faces. Adjusting last-index for this section.\n";
            s.headerinfo[2]= s.headerinfo[1] + s.mface.size();
        }
        f_ >> ch; // Eat );
        SkipSection( false);
        break;
    }
}


void
ReadMeshBuilderCL::ReadCell()
{
    if (msg_)
       *msg_ << "XF_CELL: " << std::endl;
    HeaderInfoCL hi= ReadHeaderInfoHex();
    if (hi.empty())
        throw DROPSErrCL( "ReadMeshBuilderCL::ReadCell(): Empty Header.\n");
    switch (hi[0]) {
      case 0: // Declare total number of cells, no cells defined.
        cells_.num_expected= hi[2];
        if (hi[2] == 0)
            throw DROPSErrCL("ReadMeshBuilderCL::ReadCell(): Error: No volume-mesh in file.\n");
        SkipSection( false);
        break;
      default:
        if (hi[3] == 0) { // dead-zone, ignore
            if (msg_)
                *msg_ << "ReadMeshBuilderCL::ReadCell(): Skipping dead zone " << hi[0]
                      << '\n';
            SkipSection( false);
            break;
        }
        if (hi[4] != 2)
            throw DROPSErrCL("ReadMeshBuilderCL::ReadCell(): Error: Cells are not tetrahedral.\n");
        cells_.section.push_back( CellSectionCL());
        CellSectionCL& cs= cells_.section.back();
        cs.headerinfo= hi;
        SkipSection( false);
        break;
    }
}

void
ReadMeshBuilderCL::ReadFile()
{
    while(NextSection()) {
        switch (ReadId()) {
          case 0: // XF_COMMENT
            if (msg_) {
                std::string s;
                MeshStringCL ms;
                f_ >> ms >> s;
                *msg_ << std::string("XF_COMMENT: ") + s << std::endl;
            }
            SkipSection( false);
            break;
          case 1: // XF_HEADER
            if (msg_) {
                std::string s;
                MeshStringCL ms;
                f_ >> ms >> s;
                *msg_ << std::string("XF_HEADER: ") + s << std::endl;
            }
            SkipSection( false);
            break;
          case 2: // XF_DIMENSION
            if (msg_)
                *msg_ << "XF_DIMENSION: ignored" << std::endl;
            SkipSection( false);
            break;
          case 10: // XF_NODE
            ReadNode();
            break;
          case 12: // XF_CELL
            ReadCell();
            break;
          case 13: // XF_FACE
            ReadFace();
            break;
          default:
            if (msg_)
                *msg_ << "ReadMeshBuilderCL::ReadFile(): Skipping unknown section.\n";
            SkipSection( false);
            break;
        }
    }
}

void
ReadMeshBuilderCL::AddVertexBndDescription(VertexCL* vp, Uint bndidx) const
{
    BndPointCL bpt( bndidx, std_basis<2>( 0));
    if (vp->IsOnBoundary() && std::binary_search( vp->GetBndVertBegin(), vp->GetBndVertEnd(),
                                                  bpt, BndPointSegLessCL()))
        return;
    vp->AddBnd( bpt);
    vp->BndSort();
}

void
ReadMeshBuilderCL::CreateUpdateBndEdge(MultiGridCL::EdgeLevelCont& edges,
                                       VertexCL* vp0, VertexCL* vp1,
                                       Uint bndidx) const
{
    if (vp1->GetId() < vp0->GetId() ) // Take care: We have to look for the edge at its first vertex.
        std::swap( vp0, vp1); // Swap pointer, not pointee.
    EdgeCL* ep;
    if ((ep= vp0->FindEdge( vp1)) != 0) {
        if (!is_in( ep->GetBndIdxBegin(), ep->GetBndIdxEnd(), bndidx))
            ep->AddBndIdx( bndidx);
    }
    else {
        edges.push_back( EdgeCL( vp0, vp1, 0, bndidx));
        ep= &edges.back();
        ep->SortVertices(); // Shoud be a nop due to the beginning of CreateUpdateEdge.
        ep->RecycleMe();
    }
}

void
ReadMeshBuilderCL::Clear() const
{
    std::vector<NodeSectionCL> tmp1;
    nodes_.section.swap( tmp1);
    std::vector<MFaceSectionCL> tmp2;
    mfaces_.section.swap( tmp2);
    std::vector<CellSectionCL> tmp3;
    cells_.section.swap( tmp3);
}

ReadMeshBuilderCL::ReadMeshBuilderCL(std::istream& f, std::ostream* msg)
        :f_( f), msg_( msg) {}

void 
ReadMeshBuilderCL::buildBoundaryImp(MultiGridCL* mgp) const
{
    BoundaryCL::SegPtrCont& Bnd= GetBnd( mgp);
    for (Uint i= 0; i<mfaces_.section.size(); ++i) {
        MFaceSectionCL& section= mfaces_.section[i];
        switch(section.headerinfo[3]) { // switch on bondary-condition.
          case 2: break; // interior faces
          default:
            Bnd.push_back( new MeshBoundaryCL( section.headerinfo[0], // zone-id
                                               section.headerinfo[3])); // the bc-type; see Mesh-File-Format C.8
            BC_.push_back( MapBC(section.headerinfo[3]));
            zone_id2bndidx_[section.headerinfo[0]]= Bnd.size()-1;
            break;
        }
    }
}

void 
ReadMeshBuilderCL::buildBoundary(MultiGridCL* mgp) const
{
    const_cast<ReadMeshBuilderCL*>( this)->ReadFile();
    buildBoundaryImp(mgp);
    Clear();
}

void
ReadMeshBuilderCL::build(MultiGridCL* mgp) const
{
    // Read the mesh file.
    const_cast<ReadMeshBuilderCL*>( this)->ReadFile(); // It is not useful that build is a
                                                       // const member-function by inheritance.
    nodes_.Check();
    mfaces_.Check();
    cells_.Check();

    AppendLevel( mgp);

    // Create boundary
    buildBoundaryImp(mgp);
    // Enter vertices into mgp;
    // We assume that the nodes are numbered consecutively, starting at 1.
    IdCL<VertexCL>::ResetCounter( 1);
    MultiGridCL::VertexLevelCont& verts= GetVertices( mgp)[0];
    std::vector<VertexCL*> va;
    va.reserve( nodes_.num_expected);
    for (Uint i= 0; i<nodes_.num_expected; ++i) {
        verts.push_back( VertexCL( nodes_.section[0].point[i], 0));
        va.push_back( &verts.back());
    }

    // We create the faces. On the fly, we gather tetra-definitions and add boundary descriptions
    // to the appropriate vertices.
    // To obtain correct boundary-descriptions, boundary-edges are created from boundary-faces.
    // We assume that the MFace-Sections are sorted by ascending index, starting at 1. (Calling .Check() does this.)
    IdCL<FaceCL>::ResetCounter( 1);
    MultiGridCL::FaceLevelCont& faces= GetFaces(mgp)[0];
    std::vector<FaceCL*> fa;
    std::vector<HybridCellCL> thecells( cells_.num_expected);
    MultiGridCL::EdgeLevelCont& edges= GetEdges(mgp)[0];
    for (Uint s= 0; s<mfaces_.section.size(); ++s) {
        MFaceSectionCL& section= mfaces_.section[s];
        switch(section.headerinfo[3]) { // switch on the boundary condition.
          case 2: // Inner faces
            for (Uint i= 0; i<section.mface.size(); ++i) {
                faces.push_back( FaceCL( 0)); // Default is no boundary segment.
                fa.push_back( &faces.back());
                thecells[section.mface[i][3]-1].push_back( &faces.back(), section.mface[i]);
                thecells[section.mface[i][4]-1].push_back( &faces.back(), section.mface[i]);
            }
            break;
          default: // boundary-faces
            for (Uint i= 0; i<section.mface.size(); ++i) {
                faces.push_back( FaceCL( 0, zone_id2bndidx_[section.headerinfo[0]])); // Default is no boundary segment.
                fa.push_back( &faces.back());
                if (section.mface[i][3] != 0) // A cell on the right
                    thecells[section.mface[i][3]-1].push_back( &faces.back(), section.mface[i]);
                if (section.mface[i][4] != 0) // A cell on the left
                    thecells[section.mface[i][4]-1].push_back( &faces.back(), section.mface[i]);
                AddVertexBndDescription( va[section.mface[i][0]-1], zone_id2bndidx_[section.headerinfo[0]]);
                AddVertexBndDescription( va[section.mface[i][1]-1], zone_id2bndidx_[section.headerinfo[0]]);
                AddVertexBndDescription( va[section.mface[i][2]-1], zone_id2bndidx_[section.headerinfo[0]]);
                CreateUpdateBndEdge( edges, va[section.mface[i][0]-1],
                                     va[section.mface[i][1]-1],
                                     zone_id2bndidx_[section.headerinfo[0]]);
                CreateUpdateBndEdge( edges, va[section.mface[i][0]-1],
                                     va[section.mface[i][2]-1],
                                     zone_id2bndidx_[section.headerinfo[0]]);
                CreateUpdateBndEdge( edges, va[section.mface[i][1]-1],
                                     va[section.mface[i][2]-1],
                                     zone_id2bndidx_[section.headerinfo[0]]);

            }
            break;
        }
    }

    // Finally, build the tetras. This step also creates interior edges.
    IdCL<TetraCL>::ResetCounter( 1);
    MultiGridCL::TetraLevelCont& tetras= GetTetras( mgp)[0];
    std::vector<TetraCL*> ta;
    ta.reserve( cells_.num_expected);
    for (Uint i= 0; i<thecells.size(); ++i) {
        thecells[i].Check();
        std::vector<Uint> vi= thecells[i].Vertices();
        tetras.push_back( TetraCL( va[vi[0]-1], va[vi[1]-1], va[vi[2]-1], va[vi[3]-1], 0));
        ta.push_back( &tetras.back());
        tetras.back().BuildEdges( edges);
        for (Uint f= 0; f<4; ++f) {
            thecells[i].face( f)->LinkTetra( &tetras.back());
            tetras.back().SetFace( f, thecells[i].Face( vi[VertOfFace( f, 0)],
                                                        vi[VertOfFace( f, 1)],
                                                        vi[VertOfFace( f, 2)]));
        }
    }
    // Empty recycle-bins used to build edges.
    std::for_each( va.begin(), va.end(), std::mem_fun( &VertexCL::DestroyRecycleBin));
    // Should be a nop by construction in the builder, but one never knows...
    mgp->MakeConsistentNumbering();
    Clear(); // save memory.
}

const char*
ReadMeshBuilderCL::Symbolic(const Uint id)
{
    switch (id) {
      case 0:   return SymbolicName_[0];
      case 1:   return SymbolicName_[1];
      case 2:   return SymbolicName_[2];
      case 10:  return SymbolicName_[3];
      case 12:  return SymbolicName_[4];
      case 13:  return SymbolicName_[5];
      default:  return "UNKNOWN_TO_DROPS";
    }
}

BndCondT ReadMeshBuilderCL::MapBC( Uint gambit_bc)
{
    switch(gambit_bc)
    { // gambit code: drops bc          gambit name
      case  3: return WallBC;           // wall
      case  4: return NatBC;            // pressure-inlet
      case 10: return DirBC;            // velocity-inlet
      case 36: return OutflowBC;        // outflow
      case 12: return Per1BC;           // periodic
      case  8: return Per2BC;           // periodic-shadow
      default: return UndefinedBC_;
    }
}


/*******************************************************************
*   F I L E B U I L D E R  C L                                    *
*******************************************************************/

void FileBuilderCL::BuildVerts(MultiGridCL* mgp) const
{
    std::ifstream vertex_file((path_+"Vertices").c_str());
    std::ifstream bndvtx_file((path_+"BoundaryVertices").c_str());
    CheckFile(vertex_file);
    CheckFile(bndvtx_file);
    MultiGridCL::VertexCont& verts= GetVertices(mgp);
    size_t idx=0;

    // Input-Variables
    size_t id=0, max_id= 0;
    Point3DCL point;
    Uint level=0, oldlevel=0;
    bool rmmark;

    // Boundary
    size_t bndidx=0;
    bndvtx_file >> bndidx;
    Point2DCL p2d;
    size_t bidx;

    while (!vertex_file.eof())
    {
        idx++;
        vertex_file >> id 
                    >> point[0] >> point[1] >> point[2]
                    >> level    >> rmmark;
        if (level>oldlevel) {
            AppendLevel(mgp);  // append Level in ALL lists
            oldlevel=level;
        }

        max_id= std::max( max_id, id);
        verts[level].push_back( VertexCL ( point, level, IdCL<VertexCL>( id)));
        if (rmmark) verts[level].back().SetRemoveMark();

        // BndVerts
        while (idx==bndidx)
        {
            bndvtx_file >> bidx 
                        >> p2d[0] >> p2d[1];
            verts[level].back().AddBnd( BndPointCL(bidx, p2d));
            if (!bndvtx_file.eof()) bndvtx_file >> bndidx; else bndidx=0;
        }
        vertexAddressMap[idx]= &verts[level].back();
    }
    IdCL<VertexCL>::ResetCounter( max_id + 1);
    CheckFile(vertex_file);
    CheckFile(bndvtx_file);
}

void FileBuilderCL::BuildEdges(MultiGridCL* mgp) const
{
    std::ifstream edge_file((path_+"Edges").c_str());
    CheckFile(edge_file);
    MultiGridCL::EdgeCont& edges= GetEdges(mgp);
    size_t idx=0;

    size_t vert0, vert1, midvert;
    Uint level=0;
    BndIdxT bnd0, bnd1;
    short int mfr=0;
    bool rmmark;

    while (!edge_file.eof())
    {
        idx++;
        edge_file >> vert0   >> vert1
                  >> midvert >> bnd0
                  >> bnd1    >> mfr
                  >> level   >> rmmark;

        VertexCL* vertex0   (vertexAddressMap[vert0]);
        VertexCL* vertex1   (vertexAddressMap[vert1]);
        VertexCL* midvertex (vertexAddressMap[midvert]);
        edges[level].push_back( EdgeCL (vertex0, vertex1, level, bnd0, bnd1, mfr));
        edges[level].back().SetMidVertex (midvertex);
        if (rmmark) edges[level].back().SetRemoveMark();
        edgeAddressMap[idx]= &edges[level].back();
    }
    CheckFile(edge_file);
}

void FileBuilderCL::BuildFacesI(MultiGridCL* mgp) const
{
    std::ifstream face_file((path_+"Faces").c_str());
    CheckFile(face_file);
    MultiGridCL::FaceCont& faces= GetFaces(mgp);
    size_t idx=0;
    size_t tmp;

    Uint level=0;
    BndIdxT bnd;
    bool rmmark;

    while (!face_file.eof())
    {
        idx++;
        for (int i=0; i<4; ++i) 
            face_file >> tmp; // 4 Neighbors
        face_file >> bnd >> level
                  >> rmmark;

        faces[level].push_back(FaceCL (level, bnd));
        if (rmmark) faces[level].back().SetRemoveMark();
        faceAddressMap[idx]= &faces[level].back();
    }
    CheckFile(face_file);
}

void FileBuilderCL::BuildTetras(MultiGridCL* mgp) const
{
    std::ifstream tetra_file((path_+"Tetras").c_str());
    CheckFile(tetra_file);
    std::string buffer;
    MultiGridCL::TetraCont& tetras= GetTetras(mgp);
    size_t idx=0;
    size_t id=0, max_id=0;
    Uint refrule, refmark;
    Uint level=0;

    VertexCL* verts[4];
    size_t    vertaddr[4];
    size_t    edgeaddr[6];
    size_t    faceaddr[4];
    size_t    parent;

    while (!tetra_file.eof())
    {
        idx++;
        tetra_file >> id          >> level
                   >> refrule     >> refmark
                   >> vertaddr[0] >> vertaddr[1]
                   >> vertaddr[2] >> vertaddr[3]
                   >> edgeaddr[0] >> edgeaddr[1]
                   >> edgeaddr[2] >> edgeaddr[3]
                   >> edgeaddr[4] >> edgeaddr[5]
                   >> faceaddr[0] >> faceaddr[1]
                   >> faceaddr[2] >> faceaddr[3]
                   >> parent;

        max_id= std::max( max_id, id);
        for (Uint i=0; i<4; ++i) verts[i]=vertexAddressMap[vertaddr[i]];
        TetraCL* par = tetraAddressMap[parent];
        tetras[level].push_back( TetraCL (verts[0], verts[1], verts[2], verts[3], par, IdCL<TetraCL>( id)));
        tetras[level].back().SetRefRule(refrule);
        tetras[level].back().SetRefMark (refmark);

        for (Uint i=0; i<6; ++i) {
            tetras[level].back().SetEdge(i, edgeAddressMap[edgeaddr[i]]);
        }

        for (Uint i=0; i<4; ++i) {
            tetras[level].back().SetFace(i, faceAddressMap[faceaddr[i]]);
        }

        tetraAddressMap[idx]= &tetras[level].back();
    }
    IdCL<TetraCL>::ResetCounter( max_id + 1);
    CheckFile(tetra_file);
}

void FileBuilderCL::BuildFacesII(MultiGridCL* mgp) const
{
    std::ifstream face_file((path_+"Faces").c_str());
    CheckFile(face_file);
    size_t idx=0;
    size_t neighbor[4];

    Uint level=0;
    BndIdxT bnd;
    bool rmmark;
    MultiGridCL::FaceIterator it= mgp->GetAllFaceBegin();

    while (!face_file.eof())
    {
        idx++;
        for (int i=0; i<4; ++i) 
            face_file >> neighbor[i]; // 4 Neighbors
        face_file >> bnd >> level
                  >> rmmark;
        for (Uint i=0; i<4; ++i) {
            it->SetNeighbor(i, tetraAddressMap[neighbor[i]]);
        }
        ++it;
    }
    CheckFile(face_file);
}

void FileBuilderCL::AddChildren() const
{
    std::ifstream child_file((path_+"Children").c_str());
    CheckFile(child_file);
    size_t tetra=0;
    size_t child;

    while (!child_file.eof())
    {
        child_file >> tetra;
        for (Uint i=0; (i<MaxChildrenC) && (tetra!=0); ++i) {
            child_file >> child;
            if (child != 0) tetraAddressMap[tetra]->SetChild(i, tetraAddressMap[child]);
        }
    }
    CheckFile(child_file);
}

void FileBuilderCL::build(MultiGridCL* mgp) const
{
    AppendLevel(mgp);

    // Create vertices
    std::cout << "Building Vertices ";
    BuildVerts(mgp);
    std::cout << "--> success\n";

    // Create Edges
    std::cout << "Building Edges ";
    BuildEdges(mgp);
    std::cout << "--> success\n";

    // Create Faces (without Neighbors)
    std::cout << "Building Faces Part I ";
    BuildFacesI (mgp);
    std::cout << "--> success\n";

    // Create Tetras
    std::cout << "Building Tetras ";
    BuildTetras (mgp);
    AddChildren ();
    std::cout << "--> success\n";

    FinalizeModify(mgp);

    // Link Tetras to Faces
    std::cout << "Building Faces Part II ";
    BuildFacesII(mgp);
    std::cout << "--> success\n";

    // Build Boundary
    std::cout << "Building Boundary ";
    buildBoundary(mgp);
    std::cout << "--> success\n";
}

void FileBuilderCL::CheckFile( const std::ifstream& is) const
{
    if (!is) throw DROPSErrCL( "FileBuilderCL: error while opening file!");
}

/*******************************************************************
*   M G S E R I A L I Z A T I O N   C L                           *
*******************************************************************/

template<class itT, class T>
void MGSerializationCL::GetAddr (itT b, itT e, std::map<T, size_t> &m)
{
    int i=1;
    for (itT it=b; it != e; ++it, ++i)
        m[&*it]=i;
}

void MGSerializationCL::WriteEdges()
{
    std::ofstream edge_file((path_+"Edges").c_str());
    CheckFile(edge_file);
    int i=0;
    for (MultiGridCL::EdgeIterator p=mg_.GetAllEdgeBegin(); p!=mg_.GetAllEdgeEnd(); ++p, ++i) {
        if (i!=0) edge_file << '\n';
        edge_file << vertexAddressMap[p->GetVertex(0)]   << " " << vertexAddressMap[p->GetVertex(1)] << " "
                  << vertexAddressMap[p->GetMidVertex()] << " " << *p->GetBndIdxBegin()        << " " 
                  << *(p->GetBndIdxBegin() + 1)          << " " << p->GetMFR()                 << " " 
                  << p->GetLevel()                       << " " << p->IsMarkedForRemovement();
    }
    CheckFile(edge_file);
}

void MGSerializationCL::WriteFaces()
{
    std::ofstream face_file((path_+"Faces").c_str());
    CheckFile(face_file);
    int i=0;
    for (MultiGridCL::FaceIterator p=mg_.GetAllFaceBegin(); p!=mg_.GetAllFaceEnd(); ++p, ++i) {
        if (i!=0) face_file << '\n';
        face_file << tetraAddressMap[p->GetNeighbor(0)]   << " " << tetraAddressMap[p->GetNeighbor(1)] << " "
                  << tetraAddressMap[p->GetNeighbor(2)]   << " " << tetraAddressMap[p->GetNeighbor(3)] << " "
                  << p->GetBndIdx()                       << " " << p->GetLevel()               << " "
                  << p->IsMarkedForRemovement();
    }
    CheckFile(face_file);
}

void MGSerializationCL::WriteVertices()
{
    std::ofstream vertex_file((path_+"Vertices").c_str());
    std::ofstream bndvtx_file((path_+"BoundaryVertices").c_str()); 
    CheckFile(vertex_file);
    CheckFile(bndvtx_file);
    int i=0, j=0;
    for (MultiGridCL::VertexIterator p=mg_.GetAllVertexBegin(); p!=mg_.GetAllVertexEnd(); ++p, ++i) {
        if (i!=0) vertex_file << '\n';
        vertex_file << p->GetId().GetIdent() << " " << std::scientific << std::setprecision(16) 
                    << p->GetCoord() //<< " "
                    << p->GetLevel()         << " " << p->IsMarkedForRemovement();// <<'\n';
        if (p->IsOnBoundary()) {
            for (VertexCL::const_BndVertIt it= p->GetBndVertBegin(); it != p->GetBndVertEnd(); ++it, ++j) {
                if (j!=0) bndvtx_file << '\n';
                bndvtx_file << vertexAddressMap[&*p] << " " << it->GetBndIdx() << " " 
                            << std::scientific << std::setprecision(16) << it->GetCoord2D()[0] << " " << it->GetCoord2D()[1];
            }
        }
    }
    CheckFile(vertex_file);
    CheckFile(bndvtx_file);
}

void MGSerializationCL::WriteTetras()
{
    std::ofstream tetra_file ((path_+"Tetras").c_str());
    std::ofstream child_file ((path_+"Children").c_str());
    CheckFile(tetra_file);
    CheckFile(child_file);
    bool start=true, child_start=true;
    for (MultiGridCL::TetraIterator p=mg_.GetAllTetraBegin(); p!=mg_.GetAllTetraEnd(); ++p) {
        if (!start) tetra_file << '\n';
        tetra_file << p->GetId().GetIdent()             << " " << p->GetLevel() << " " 
                   << p->GetRefRule()                   << " " << p->GetRefMark() << " "
                   << vertexAddressMap[p->GetVertex(0)] << " " << vertexAddressMap[p->GetVertex(1)] << " "
                   << vertexAddressMap[p->GetVertex(2)] << " " << vertexAddressMap[p->GetVertex(3)] << " "
                   << edgeAddressMap[p->GetEdge(0)]     << " " << edgeAddressMap[p->GetEdge(1)] << " "
                   << edgeAddressMap[p->GetEdge(2)]     << " " << edgeAddressMap[p->GetEdge(3)] << " "
                   << edgeAddressMap[p->GetEdge(4)]     << " " << edgeAddressMap[p->GetEdge(5)] << " "
                   << faceAddressMap[p->GetFace(0)]     << " " << faceAddressMap[p->GetFace(1)] << " "
                   << faceAddressMap[p->GetFace(2)]     << " " << faceAddressMap[p->GetFace(3)] << " "
                   << tetraAddressMap[p->GetParent()];
        if (!p->IsUnrefined()) {
            if (!child_start) child_file << '\n';
            else child_start=false;
            child_file << tetraAddressMap[&*p] << " ";
            for (Uint i=0; i<MaxChildrenC - 1; ++i) {
                child_file << tetraAddressMap[p->GetChild(i)] << " ";
            }
            child_file << tetraAddressMap[p->GetChild(MaxChildrenC - 1)];
        }
        start=false;
    }
    CheckFile(tetra_file);
    CheckFile(child_file);
}

void MGSerializationCL::CreateAddrMaps()
{
    GetAddr (mg_.GetAllEdgeBegin(),   mg_.GetAllEdgeEnd(),     edgeAddressMap);
    GetAddr (mg_.GetAllVertexBegin(), mg_.GetAllVertexEnd(), vertexAddressMap);
    GetAddr (mg_.GetAllFaceBegin(),   mg_.GetAllFaceEnd(),     faceAddressMap);
    GetAddr (mg_.GetAllTetraBegin(),  mg_.GetAllTetraEnd(),   tetraAddressMap);
}

void MGSerializationCL::WriteMG()
{
    CreateAddrMaps();

    // Write vertices
    std::cout << "Writing Vertices ";
    WriteVertices();
    std::cout << "--> success\n";

    // Write Edges
    std::cout << "Writing Edges ";
    WriteEdges();
    std::cout << "--> success\n";

    // Write Tetras
    std::cout << "Writing Tetras ";
    WriteTetras();
    std::cout << "--> success\n";    

    // Write Faces
    std::cout << "Writing Faces ";
    WriteFaces();
    std::cout << "--> success\n";
}

void MGSerializationCL::CheckFile( const std::ofstream& os) const
{
    if (!os) throw DROPSErrCL( "MGSerializationCL: error while opening file!");
}

} //end of namespace DROPS
