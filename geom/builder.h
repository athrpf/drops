//**************************************************************************
// File:    builder.h                                                      *
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


#ifndef DROPS_BUILDER_H
#define DROPS_BUILDER_H


#include "geom/multigrid.h"


namespace DROPS
{


class BrickBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    // for vertices:
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    // for tetras:
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

  public:
    BrickBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&, Uint, Uint, Uint);

    virtual void
    build(MultiGridCL*) const;
};


class LBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    const Uint _b1;
    const Uint _b2;
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

  public:
    LBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&,
               Uint, Uint, Uint, Uint, Uint);

    virtual void
    build(MultiGridCL*) const;
};


class BBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    const Uint _b1;
    const Uint _b2;
    const Uint _b3;
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

  public:
    BBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&,
               Uint, Uint, Uint, Uint, Uint, Uint);

    virtual void
    build(MultiGridCL*) const;
};


class TetraBuilderCL : public MGBuilderCL
{
  private:
    const Ubyte rule_;

    const Point3DCL p0_;
    const Point3DCL p1_;
    const Point3DCL p2_;
    const Point3DCL p3_;
    
  public:

    TetraBuilderCL(Ubyte rule);
    TetraBuilderCL(Ubyte rule, const Point3DCL& p0, const Point3DCL& p1,
                               const Point3DCL& p2, const Point3DCL& p3);

    virtual void
    build(MultiGridCL*) const;
};

} //end of namespace DROPS

#endif
