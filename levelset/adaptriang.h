//**************************************************************************
// File:    adaptriang.h                                                   *
// Content: adaptive triangulation based on position of the interface      *
//          provided by the levelset function                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_ADAPTRIANG_H
#define DROPS_ADAPTRIANG_H

#include "levelset/levelset.h"

namespace DROPS
{

class AdapTriangCL
{
  private:
    MultiGridCL& mg_;
    double width_;
    Uint c_level_, f_level_;
    bool modified_;
    
    template <class DistFctT>
    double GetValue( DistFctT& dist, const VertexCL& v)      { return dist.val( v); }
    template <class DistFctT>
    double GetValue( DistFctT& dist, const EdgeCL& e)        { return dist.val( e); }
    template <class DistFctT>
    double GetValue( DistFctT& dist, const TetraCL& t)       { return dist.val( t, 0.25, 0.25, 0.25); }
    double GetValue( scalar_fun_ptr dist, const VertexCL& v) { return dist( v.GetCoord() ); }
    double GetValue( scalar_fun_ptr dist, const EdgeCL& e)   { return dist( GetBaryCenter( e) ); }
    double GetValue( scalar_fun_ptr dist, const TetraCL& t)  { return dist( GetBaryCenter( t) ); }

    template <class DistFctT>
    bool ModifyGridStep( DistFctT&);
    // One step of grid change; returns true if modifications were necessary,
    // false, if nothing changed.

  public:
    AdapTriangCL( MultiGridCL& mg, double width, Uint c_level, Uint f_level)
      : mg_(mg), width_(width), c_level_(c_level), f_level_(f_level), modified_(false) 
      { Assert( 0<=c_level && c_level<=f_level, "AdapTriangCL: Levels are cheesy.\n", ~0); }
    
    template <class DistFctT>
    void MakeInitialTriang( DistFctT&);

    template <class StokesT>
    void UpdateTriang( StokesT&, LevelsetP2CL&);

    bool WasModified() const { return modified_; }
};

#include "levelset/adaptriang.tpp"

} // end of namespace DROPS

#endif

