//**************************************************************************
// File:     boundary.cpp                                                  *
// Content:  Boundary description                                          *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.2                                                           *
// History:  begin - 09.10.2000                                            *
//                                                                         *
// Remarks:                                                                *
//**************************************************************************

#include <numeric>
#include "geom/boundary.h"

namespace DROPS
{

//**************************************************************************
// Class:    AffineSquareCL                                                *
// Base:     BndSegCL                                                      *
// Purpose:  affine image of the 2D-unit-square (0,0), (1,0), (1,1), (0,1) *
// Remarks:                                                                *
//**************************************************************************


AffineSquareCL::AffineSquareCL (const Point3DCL& p00, const Point3DCL& p10, const Point3DCL& p01)
    : BndSegCL(false), _Orig(p00), _D0(p10-p00), _D1(p01-p00),
      _D0D1(std::inner_product(_D0.begin(), _D0.end(), _D1.begin(), 0.0)),
      _D0sq(std::inner_product(_D0.begin(), _D0.end(), _D0.begin(), 0.0)),
      _D1sq(std::inner_product(_D1.begin(), _D1.end(), _D1.begin(), 0.0)),
      _Det(_D0sq*_D1sq-_D0D1*_D0D1) {}


bool AffineSquareCL::IsInBounds (const Point2DCL& p2D) const
{
    return p2D[0]>=0.0 && p2D[0]<=1.0 && p2D[1]>=0.0 && p2D[1]<=1.0;
}


Point3DCL AffineSquareCL::Map (const Point2DCL& p2D) const
{
    return _Orig + p2D[0]*_D0 + p2D[1]*_D1;
}


Point2DCL AffineSquareCL::ProjectRaw (const Point3DCL& p3D, const Point2DCL*) const
{
    Point3DCL Diff(p3D-_Orig);
    double    s0(std::inner_product(_D0.begin(), _D0.end(), Diff.begin(), 0.0));
    double    s1(std::inner_product(_D1.begin(), _D1.end(), Diff.begin(), 0.0));

    Point2DCL Result;
    Result[0] = (s0*_D1sq - s1*_D0D1) / _Det;
    Result[1] = (s1*_D0sq - s0*_D0D1) / _Det;
    return Result;
}


Point2DCL AffineSquareCL::Project (const Point3DCL& p3D, const Point2DCL*) const
{
    Point2DCL Result(ProjectRaw(p3D));
    // TODO: Auf den Rand projizieren
    return Result;
}


BndPairCL AffineSquareCL::MidProject (const BndPointCL& bp0, const BndPointCL& bp1) const
// We don't need the edge-pointer in the affine case
{
    Point2DCL bc = BaryCenter( bp0.GetCoord2D(), bp1.GetCoord2D() );
    return std::make_pair( bc, Map(bc) );
}


//**************************************************************************
// Class:    AffineTriangleCL                                              *
// Base:     BndSegCL                                                      *
// Purpose:  affine image of the reference-triangle (0,0), (1,0), (0,1)    *
//**************************************************************************

AffineTriangleCL::AffineTriangleCL (const Point3DCL& p00, const Point3DCL& p10, const Point3DCL& p01)
    : BndSegCL(false), _Orig(p00), _D0(p10-p00), _D1(p01-p00),
      _D0D1(std::inner_product(_D0.begin(), _D0.end(), _D1.begin(), 0.0)),
      _D0sq(std::inner_product(_D0.begin(), _D0.end(), _D0.begin(), 0.0)),
      _D1sq(std::inner_product(_D1.begin(), _D1.end(), _D1.begin(), 0.0)),
      _Det(_D0sq*_D1sq-_D0D1*_D0D1) {}


bool AffineTriangleCL::IsInBounds (const Point2DCL& p2D) const
{
    return p2D[0]>=0.0 && p2D[1]>=0.0 && p2D[0]+p2D[1]<=1.0;
}


Point3DCL AffineTriangleCL::Map (const Point2DCL& p2D) const
{
    return _Orig + p2D[0]*_D0 + p2D[1]*_D1;
}


Point2DCL AffineTriangleCL::ProjectRaw (const Point3DCL& p3D, const Point2DCL*) const
{
    Point3DCL Diff(p3D-_Orig);
    double    s0(std::inner_product(_D0.begin(), _D0.end(), Diff.begin(), 0.0));
    double    s1(std::inner_product(_D1.begin(), _D1.end(), Diff.begin(), 0.0));

    Point2DCL Result;
    Result[0] = (s0*_D1sq - s1*_D0D1) / _Det;
    Result[1] = (s1*_D0sq - s0*_D0D1) / _Det;
    return Result;
}


Point2DCL AffineTriangleCL::Project (const Point3DCL& p3D, const Point2DCL*) const
{
    Point2DCL Result(ProjectRaw(p3D));
    // TODO: Auf den Rand projizieren
    return Result;
}


BndPairCL AffineTriangleCL::MidProject (const BndPointCL& bp0, const BndPointCL& bp1) const
// We don't need the edge-pointer in the affine case
{
    Point2DCL bc = BaryCenter( bp0.GetCoord2D(), bp1.GetCoord2D() );
    return std::make_pair( bc, Map(bc) );
}

} // end of namespace DROPS
