//**************************************************************************
// File:    builder.cpp                                                    *
// Content: MGBuilderCL objects for some domains                           *
// Author:  Joerg Peters, Volker Reichelt, IGPM RWTH Aachen                *
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


void BrickBuilderCL::build (MultiGridCL* mgp) const
{
    // Check, if the parallelepiped spanned by e1,e2,e3 is degenerated

    AppendLevel(mgp);

    // Create boundary
    BoundaryCL::SegPtrCont& Bnd= GetBnd(mgp);

    Bnd.push_back( new AffineSquareCL(_orig,     _orig    +_e2, _orig    +_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e1, _orig+_e1+_e2, _orig+_e1+_e3) ); // e2-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig,     _orig    +_e1, _orig    +_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e2, _orig+_e2+_e1, _orig+_e2+_e3) ); // e1-e3-plane
    Bnd.push_back( new AffineSquareCL(_orig,     _orig    +_e1, _orig    +_e2) ); // e1-e2-plane
    Bnd.push_back( new AffineSquareCL(_orig+_e3, _orig+_e3+_e1, _orig+_e3+_e2) ); // e1-e2-plane

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


void LBuilderCL::build (MultiGridCL* mgp) const
{
    const double _dn1= static_cast<double>(_n1);
    const double _dn2= static_cast<double>(_n2);
    const double _dn3= static_cast<double>(_n3);
    const double _db1= static_cast<double>(_b1);
    const double _db2= static_cast<double>(_b2);
    
    const double b1= _db1/_dn1;
    const double b2= _db2/_dn2;

    // Check, if the parallelepiped spanned by e1,e2,e3 is degenerated

    AppendLevel(mgp);

    // Create boundary
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


void BBuilderCL::build (MultiGridCL* mgp) const
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

    // Check, if the parallelepiped spanned by e1,e2,e3 is degenerated

    AppendLevel(mgp);

    // Create boundary
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


TetraBuilderCL::TetraBuilderCL(const Ubyte rule)
    :rule_(rule), p0_( std_basis<3>(0)), p1_( std_basis<3>(1)),
                  p2_( std_basis<3>(2)), p3_( std_basis<3>(3))
{}

TetraBuilderCL::TetraBuilderCL(const Ubyte rule,  const Point3DCL& p0, const Point3DCL& p1,
                                                  const Point3DCL& p2, const Point3DCL& p3)
    :rule_(rule), p0_( p0), p1_( p1), p2_( p2), p3_( p3)
{}


void TetraBuilderCL::build(MultiGridCL* mgp) const
{
    AppendLevel( mgp);

    // Create boundary
    BoundaryCL::SegPtrCont& Bnd= GetBnd( mgp);
    Bnd.push_back( new AffineTriangleCL( p0_, p1_, p2_)); // e1-e2-plane
    Bnd.push_back( new AffineTriangleCL( p0_ ,p1_, p3_)); // e1-e3-plane
    Bnd.push_back( new AffineTriangleCL( p0_, p2_, p3_)); // e1-e2-plane
    Bnd.push_back( new AffineTriangleCL( p1_, p2_, p3_)); // lid

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
    mgp->Refine();
}


} //end of namespace DROPS
