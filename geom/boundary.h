//**************************************************************************
// File:     boundary.h                                                    *
// Content:  Boundary description                                          *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.2                                                           *
// History:  begin - 09.10.2000                                            *
//                                                                         *
// Remarks:                                                                *
//**************************************************************************

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_


#include "geom/parallel.h"
#include "misc/container.h"


namespace DROPS
{


class EdgeCL;
class BndPointCL;

typedef std::pair<Point2DCL, Point3DCL> BndPairCL;
typedef Usint BndIdxT;
const BndIdxT NoBndC= static_cast<BndIdxT>(-1);

//**************************************************************************
// Class:    BndSegCL                                                      *
// Purpose:  base class for the boundary descriptions used by DROPS        *
// Remarks:  To define your own boundary segment descriptions inherit from *
//           this class.                                                   *
//**************************************************************************

class BndSegCL
{
  private:
    bool           _IsNonPlanar;

  public:
    BndSegCL (bool IsNonPl)
        :_IsNonPlanar(IsNonPl) {}
    virtual ~BndSegCL () {}
    // Default copy-ctor, assignment-op.

    bool                  IsNonPlanar () const  { return _IsNonPlanar; }

    // true, iff p2D is in the domain of the boundary segment
    virtual bool      IsInBounds (const Point2DCL&) const = 0;
    // If IsInBounds(p2D)==true, the result is the image of p2D, else undefined
    virtual Point3DCL Map        (const Point2DCL&) const = 0;
    // Project p3D into to parameter space. It is not checked, whether p2D is in
    // the actual domain of the boundary segment. The last argument can be a hint,
    // where p2D lies - it is implementation dependant, whether it is used.
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = NULL) const = 0;
    // Map p3D into the domain
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = NULL) const = 0;
    // Calculate the midvertex coordinates of an edge, the endpoints of which have boundary-descriptions
    // bp1 & bp2. The midvertex shall lie on the boundary.
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const = 0;
};


//**************************************************************************
// Class:    BndPointCL                                                    *
// Purpose:  stores the boundary-description of a vertex                   *
// Remarks:  _Coord2D are the parametric coordinates of the vertex in the  *
//           boundary segment *_BndIdx.                                    *
//**************************************************************************

class BndPointCL
{
  private:
    BndIdxT   _BndIdx;
    Point2DCL _Coord2D;

  public:
    BndPointCL (BndIdxT BndIdx, const Point2DCL& p2D)
        : _BndIdx(BndIdx), _Coord2D(p2D) {}
    // Default copy-ctor, dtor, assignnment-op.

    // Returns the BndSegCL-object, to which these 2D-coordinates belong
    BndIdxT GetBndIdx() const { return _BndIdx; }
    // Returns the 2D-coordinates of the vertex in BndSegCL *GetBndIdx()
    const Point2DCL& GetCoord2D () const { return _Coord2D; }
};


//**************************************************************************
// Class:    BndPointSegLessCL                                             *
// Purpose:  orders BndPointCL-objects by their BndSegCL-object-Id         *
// Remarks:  The boundary descriptions of a vertex are stored in this order*
//           to allow for the use of set_intersection in BuildMidVertex.   *
//**************************************************************************

class BndPointSegLessCL : public std::binary_function<BndPointCL,BndPointCL,bool>
{
  public:
    bool operator () (const BndPointCL& bp0, const BndPointCL& bp1) const
        { return bp0.GetBndIdx() < bp1.GetBndIdx(); }
};


//**************************************************************************
// Class:    BndPointSegEqCL                                               *
// Purpose:  Compare BndPointCL-objects by their BndSegCL-object-Id        *
// Remarks:                                                                *
//**************************************************************************

class BndPointSegEqCL : public std::unary_function<BndPointCL,bool>
{
  private:
    BndIdxT _bidx;
    
  public:
    BndPointSegEqCL(BndIdxT bidx) :_bidx(bidx) {}
    
    bool operator () (const BndPointCL& bp) const
        { return bp.GetBndIdx() == _bidx; }
};


//**************************************************************************
// Class:    AffineSquareCL                                                *
// Base:     BndSegCL                                                      *
// Purpose:  affine image of the 2D-unit-square (0,0), (1,0), (1,1), (0,1) *
// Remarks:                                                                *
//**************************************************************************

class AffineSquareCL : public BndSegCL
{
  private:
    Point3DCL _Orig, _D0, _D1;
    double    _D0D1, _D0sq, _D1sq, _Det;

  public:
    // Images of          (0,0),            (1,0),            (0,1).
    AffineSquareCL (const Point3DCL&, const Point3DCL&, const Point3DCL&);
    // Default copy-ctor, assignment-op.

    virtual bool      IsInBounds (const Point2DCL&) const;
    virtual Point3DCL Map        (const Point2DCL&) const;
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const;
};


//**************************************************************************
// Class:    AffineTriangleCL                                              *
// Base:     BndSegCL                                                      *
// Purpose:  affine image of the reference-triangle (0,0), (1,0), (0,1)    *
//**************************************************************************

class AffineTriangleCL : public BndSegCL
{
  private:
    Point3DCL _Orig, _D0, _D1;
    double    _D0D1, _D0sq, _D1sq, _Det;

  public:
    // Images of            (0,0),            (1,0),            (0,1).
    AffineTriangleCL (const Point3DCL&, const Point3DCL&, const Point3DCL&);
    // Default copy-ctor, assignment-op.

    virtual bool      IsInBounds (const Point2DCL&) const;
    virtual Point3DCL Map        (const Point2DCL&) const;
    virtual Point2DCL ProjectRaw (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual Point2DCL Project    (const Point3DCL&, const Point2DCL* = NULL) const;
    virtual BndPairCL MidProject (const BndPointCL&, const BndPointCL&) const;
};

} // end of namespace DROPS

#endif
