//**************************************************************************
// File:    output.h                                                       *
// Content: geometry and solution output in various formats                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 26 2001                                           *
//**************************************************************************

#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "geom/multigrid.h"
#include "misc/problem.h"

namespace DROPS
{


class MGOutCL
// virtual base class for output of geometry (+ solution)
{
  protected:
    const MultiGridCL* _MG;

  public:
    MGOutCL(const MultiGridCL* MG) : _MG(MG) {}

    virtual std::ostream& put(std::ostream&) const = 0;
};

class ColorMapperCL
{
  public:
    typedef SArrayCL<double, 4> RGBAType;
    
    virtual RGBAType map(double val) const
        { RGBAType rgba; rgba[3]= 1.0; rgba[2]= val; rgba[1]= val; rgba[0]= val; return rgba; }
    // returns an RGBA-quadruple for input between 0.0 and 1.0
};


class RBColorMapperCL: public ColorMapperCL
{
  public:
    virtual RGBAType map(double val) const
    {
        RGBAType rgba;
        rgba[3]= 1.0;
        rgba[0]= val;
        rgba[1]= 0.0;
        rgba[2]= 1.0 - val;

        return rgba;
    }
};


class DumpMGCL : public MGOutCL
{
  public:
    DumpMGCL(const MultiGridCL& MG) : MGOutCL(&MG) {}

    virtual std::ostream& put (std::ostream&) const;
};


class SanityMGOutCL : public MGOutCL
{
  public:
    SanityMGOutCL(const MultiGridCL& MG) : MGOutCL(&MG) {}

    virtual std::ostream& put (std::ostream&) const;
};



class GeomMGOutCL : public MGOutCL
// output of geometry in GeomView format
{
  private:
    Uint   _level;
    bool   _onlyBnd;
    double _explode;

  public:
    GeomMGOutCL (const MultiGridCL& MG, int TriLevel=-1, bool onlyBnd=false, double explode=0.5)
        : MGOutCL(&MG), _level( TriLevel<0 ? MG.GetLastLevel() : TriLevel ),
          _onlyBnd(onlyBnd), _explode(explode) {}

    void   SetExplode (double explode) { _explode = explode; }
    double GetExplode () const         { return _explode; }

    virtual std::ostream& put (std::ostream&) const;
};



template <class DiscSol>
class GeomSolOutCL : public MGOutCL
// output of solution in GeomView format
{
  private:
    Uint   _level;
    bool   _onlyBnd;
    double _explode;
    double _min;
    double _max;
    const ColorMapperCL* _color;
    DiscSol _discsol;

  public:
    GeomSolOutCL (const MultiGridCL& MG, const DiscSol& discsol, const ColorMapperCL* colmap, int TriLevel=-1, bool onlyBnd=false,
                  double explode=0.5, double min=0., double max=1.)
        : MGOutCL(&MG), _level( TriLevel<0 ? MG.GetLastLevel() : TriLevel ),
          _onlyBnd(onlyBnd), _explode(explode), _min(min), _max(max), _color(colmap), _discsol(discsol) {}

    void   SetExplode (double explode) { _explode = explode; }
    double GetExplode () const         { return _explode; }
    void   SetMinMax  (double min, double max) { _min= min; _max= max; }
    void   SetColorMap (const ColorMapperCL* colmap) { _color= colmap; }
    virtual std::ostream& put (std::ostream&) const;
};




template <class DiscVelT, class DiscPrT>
class TecPlot2DSolOutCL: public MGOutCL
// output of solution in TecPlot format
// on a x-/y-/z-cutplane
{
  private:
    Uint _idx;
    Uint _level;
    Uint _cut;
    double _cutplane;

    DiscVelT _v;
    DiscPrT  _p;
    
    bool IsInPlane( const Point3DCL& p) const { return p[_cut]==_cutplane; }
    
  public:
    TecPlot2DSolOutCL( const MultiGridCL& mg, const DiscVelT& v, const DiscPrT& p, 
                       const IdxDescCL& freeVertIdx, int lvl, Uint cutvar, double cutplane)
      : MGOutCL( &mg), _idx( freeVertIdx.Idx), _level( lvl<0 ? mg.GetLastLevel() : lvl), 
        _cut( cutvar), _cutplane( cutplane), _v( v), _p( p) 
    {}
      
    std::ostream& put( std::ostream&) const;
};


template <class DiscVelT, class DiscPrT>
class TecPlotSolOutCL: public MGOutCL
// output of solution in TecPlot format
{
  private:
    Uint _level;

    DiscVelT _v;
    DiscPrT  _p;
    
  public:
    TecPlotSolOutCL( const MultiGridCL& mg, const DiscVelT& v, const DiscPrT& p, int lvl= -1)
      : MGOutCL( &mg), _level( lvl<0 ? mg.GetLastLevel() : lvl), _v( v), _p( p) 
    {}
      
    std::ostream& put( std::ostream&) const;
};



//=====================================================
//              inline definitions
//=====================================================

inline std::ostream& operator << (std::ostream& os, const MGOutCL& mgout)
{
    return mgout.put(os);
}


//=====================================================
//              template definitions
//=====================================================

template <class DiscSol>
inline std::ostream& operator << (std::ostream& os, const GeomSolOutCL<DiscSol>& solout)
{
    return solout.put(os);
}


template <class DiscSol>
std::ostream&
GeomSolOutCL<DiscSol>::put(std::ostream &os) const
{
    const double val_diff= _max-_min;
    ColorMapperCL::RGBAType rgba;
//    std::ios_base::fmtflags my_format= std::ios_base::fixed|std::ios_base::showpoint;
//    std::ios_base::fmtflags old_format= os.flags(my_format);
    std::ios::fmtflags my_format= std::ios::fixed|std::ios::showpoint;
    std::ios::fmtflags old_format= os.flags(my_format);
    Assert(_level==_discsol.GetSolution()->RowIdx->TriangLevel, DROPSErrCL("GeomSolOutCL::put: wrong level"), ~0);
    os << "appearance {\n-concave\nshading smooth\n}\n";
    os << "LIST {\n";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level);
          tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)) )
            continue;
//        if (GetBaryCenter(*tit)[2]>0.55) continue;

        Point3DCL Offset(0.0);
        std::vector<double> val;
        val.reserve(4);
        _discsol.GetDoF(*tit, val);

        for (TetraCL::const_VertexPIterator it= tit->GetVertBegin(), vend= tit->GetVertEnd(); it!=vend; ++it)
            Offset+= (*it)->GetCoord();
        os << "geom { COFF 4 4 6\n";
        for ( int i=0; i<4; i++ )
        {
            for (int j=0; j<3; ++j)
                os << _explode*Offset[j]/4+tit->GetVertex(i)->GetCoord()[j] << ' ';
            rgba= _color->map( (val[i]-_min)/val_diff );
            os << rgba[0] << ' ' << rgba[1] << ' ' << rgba[2] << ' ' << rgba[3] << '\n';

        }
        os <<   "3 1 2 3"
           << "\n3 0 2 3"
           << "\n3 0 1 3"
           << "\n3 0 1 2"
           << "\n}\n";
    }
    os.flags(old_format);
    return os << '}' << std::endl;
}


template<class DiscVelT, class DiscPrT>
std::ostream& TecPlot2DSolOutCL<DiscVelT,DiscPrT>::put( std::ostream& os) const
{
    const char xyz[]= "XYZ",
               uvw[]= "UVW";
               
    const MultiGridCL::const_TriangTetraIteratorCL tbegin= _MG->GetTriangTetraBegin( _level),
                                                   tend=   _MG->GetTriangTetraEnd( _level);
               
    const MultiGridCL::TriangVertexIteratorCL vbegin= const_cast<MultiGridCL*>(_MG)->GetTriangVertexBegin( _level), 
                                              vend=   const_cast<MultiGridCL*>(_MG)->GetTriangVertexEnd( _level); 
    Uint numv= 0, numel= 0;
    for (MultiGridCL::TriangVertexIteratorCL it=vbegin; it!=vend; ++it)
        if (IsInPlane( it->GetCoord() ))
            it->Unknowns( _idx)[0]= ++numv;

    Uint VertsInPlane, NotInPlane= 0;
    for (MultiGridCL::const_TriangTetraIteratorCL it= tbegin; it!=tend; ++it)
    {
        VertsInPlane= 0;
        for (Uint v=0; v<NumVertsC; ++v)
            if (IsInPlane( it->GetVertex(v)->GetCoord()))
                ++VertsInPlane;
        if (VertsInPlane==3) // found triangle in plane!
            ++numel;
    }
            
    os << "# output written by DROPS::TecPlot2DSolOutCL,\n";
    os << "# cutplane is " << xyz[_cut] << " = " << _cutplane << "\n\n";
    os << "TITLE = " << '"' << "cutplane " << xyz[_cut] << " = " << _cutplane << '"' << '\n';
    
    os << "VARIABLES = ";
    for (Uint i=0; i<3; ++i)
        if (i!=_cut) os << '"' << xyz[i] << '"' << " , ";
    for (Uint i=0; i<3; ++i)
        if (i!=_cut) os << '"' << uvw[i] << '"' << " , ";
    os << '"' << 'P' << '"' << '\n';

    os << "ZONE F=FEPOINT, ET=TRIANGLE, E=" << numel << ", N=" << numv << "\n\n";

    // write numerical data on verts
    Point3DCL data;
    for (MultiGridCL::TriangVertexIteratorCL it=vbegin; it!=vend; ++it)
        if (IsInPlane( it->GetCoord() ))
        {
            data= it->GetCoord();  // XYZ
            for (Uint i=0; i<3; ++i)
                if (i!=_cut) os << data[i] << ' ';
                
            data= _v.val( *it);    // UVW
            for (Uint i=0; i<3; ++i)
                if (i!=_cut) os << data[i] << ' ';
            
            os << _p.val( *it) << '\n';
                
        }
    os << '\n';
    // write connectivity of the verts -> triangles
    for (MultiGridCL::const_TriangTetraIteratorCL it= tbegin; it!=tend; ++it)
    {
        VertsInPlane= 0;
        for (Uint v=0; v<NumVertsC; ++v)
            if (IsInPlane( it->GetVertex(v)->GetCoord() ))
                ++VertsInPlane;
            else
                NotInPlane= v;
        if (VertsInPlane==3) // found triangle in plane!
        {
            for (Uint v=0; v<NumVertsC; ++v)
                if (v!=NotInPlane)
                    os << it->GetVertex(v)->Unknowns(_idx)[0] << '\t';
            os << '\n';
        }
    }

    return os;
}
    
    
template<class DiscVelT, class DiscPrT>
std::ostream& TecPlotSolOutCL<DiscVelT,DiscPrT>::put( std::ostream& os) const
{
    const char xyz[]= "XYZ",
               uvw[]= "UVW";
               
    const MultiGridCL::const_TriangTetraIteratorCL tbegin= _MG->GetTriangTetraBegin( _level),
                                                   tend=   _MG->GetTriangTetraEnd( _level);
               
    const MultiGridCL::TriangVertexIteratorCL vbegin= const_cast<MultiGridCL*>(_MG)->GetTriangVertexBegin( _level), 
                                              vend=   const_cast<MultiGridCL*>(_MG)->GetTriangVertexEnd( _level); 
    Uint numv= std::distance( vbegin, vend), 
         numel= std::distance( tbegin, tend),
         pidx= _p.GetSolution()->RowIdx->Idx;
            
    os << "# output written by DROPS::TecPlotSolOutCL,\n";
    os << "TITLE = " << '"' << "3D-Solution" << '"' << '\n';
    
    os << "VARIABLES = ";
    for (Uint i=0; i<3; ++i)
        os << '"' << xyz[i] << '"' << " , ";
    for (Uint i=0; i<3; ++i)
        os << '"' << uvw[i] << '"' << " , ";
    os << '"' << 'P' << '"' << '\n';

    os << "ZONE F=FEPOINT, ET=TETRAHEDRON, E=" << numel << ", N=" << numv << "\n\n";

    // write numerical data on verts
    Point3DCL data;
    for (MultiGridCL::TriangVertexIteratorCL it=vbegin; it!=vend; ++it)
    {
        data= it->GetCoord();  // XYZ
        for (Uint i=0; i<3; ++i)
            os << data[i] << ' ';

        data= _v.val( *it);    // UVW
        for (Uint i=0; i<3; ++i)
            os << data[i] << ' ';

        os << _p.val( *it) << '\n';

    }
    os << '\n';
    // write connectivity of the verts -> tetras
    for (MultiGridCL::const_TriangTetraIteratorCL it= tbegin; it!=tend; ++it)
    {
        for (Uint v=0; v<NumVertsC; ++v)
                os << (it->GetVertex(v)->Unknowns(pidx)[0]+1) << '\t';
        os << '\n';
    }

    return os;
}
    
    
} // end of namespace DROPS

// I don't want to mess around too much.
#include "out/mapleout.h"

#endif
