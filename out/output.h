/// \file output.h
/// \brief Classes for output to geomview, maple, text-files, etc.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#ifndef DROPS_OUTPUT_H
#define DROPS_OUTPUT_H

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

    virtual ~MGOutCL() {}
};

class ColorMapperCL
{
  public:
    typedef SArrayCL<double, 4> RGBAType;

    virtual RGBAType map(double val) const
        { RGBAType rgba; rgba[3]= 1.0; rgba[2]= val; rgba[1]= val; rgba[0]= val; return rgba; }
    // returns an RGBA-quadruple for input between 0.0 and 1.0

    virtual ~ColorMapperCL() {}
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
    double _procExplode;

  public:
    GeomMGOutCL (const MultiGridCL& MG, int TriLevel=-1, bool onlyBnd=false,
                 double explode=0, double procExplode=0.5)
        : MGOutCL(&MG), _level( TriLevel<0 ? MG.GetLastLevel() : TriLevel ),
          _onlyBnd(onlyBnd), _explode(explode), _procExplode(procExplode) {}

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
      : MGOutCL( &mg), _idx( freeVertIdx.GetIdx()), _level( lvl<0 ? mg.GetLastLevel() : lvl),
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


// Plane in normal-form: n*x - r = 0
struct PlaneCL
{
  public:
    SVectorCL<3> n;
    double       r;

    PlaneCL(const SVectorCL<3>& normal=std_basis<3>(1), double rr=.5)
      : n(normal), r(rr) {}
};

// contains coordinates of 2D-polygon-vertices (in 3D-space) expressed in
// the coordinate-system attached to the tetra
typedef std::vector<SVectorCL<3> > TetraSectionT;


// Options for Maple-plots
class Maple3DOptionCL
{
  public:
    // Name of the PLOT-variable in the Maple-output.
    std::string name;
    // Title
    std::string title;
    // Labels for axes
    std::string xlabel, ylabel, zlabel;
    // how to plot the axes
    std::string axesstyle;
    // constrained or unconstrained plot
    std::string scaling;

    Maple3DOptionCL(const std::string& n= "pl", const std::string& t= "DROPS-plot",
                    const std::string& x= "X", const std::string& y= "Y", const std::string& z= "Z",
                    const std::string& axst= "BOXED", const std::string& scal= "CONSTRAINED")
        : name(n), title(t), xlabel(x), ylabel(y),zlabel(z), axesstyle(axst), scaling(scal) {}

    bool // Write title, label and axesstyle into a PLOT; assumes, that it shall not write a leading ','.
    WriteGlobalOptions(std::ostream&) const;
};


// Write the geometry or a cut thereof in Maples format
class MapleMGOutCL : public MGOutCL
{
  private:
    Uint            _level;
    bool            _onlyBnd;
    bool            _cut;
    PlaneCL         _plane;
    Maple3DOptionCL _options;

    bool
    _Write_Intersection(const TetraCL&, std::ostream&, bool) const;
    bool
    _Write_Tetra(const TetraCL&, std::ostream&, bool) const;
    void // Writes the intersection POLYGON, but does not close the last parenthesis
    _Intersect_Plot() const;

  public:
    MapleMGOutCL(const MultiGridCL& MG, int TriLevel=-1, bool onlyBnd=false,
                 bool cut= true, const PlaneCL& pl= PlaneCL() , Maple3DOptionCL opt=Maple3DOptionCL())
      : MGOutCL(&MG), _level( TriLevel<0 ? MG.GetLastLevel() : TriLevel ),
        _onlyBnd(onlyBnd), _cut(cut), _plane(pl), _options(opt) {}
    virtual ~MapleMGOutCL() {}

    void // Writes the intersection POLYGON, but does not close the last parenthesis
    _Intersect_Plot(std::ostream&) const;

    std::ostream& put (std::ostream&) const;
};


template <class DiscVelT, class DiscPrT>
class MapleSolOutCL: public MGOutCL
// output of solution in Maples format
{
  private:
    Uint            _level;
    bool            _onlyBnd;
    PlaneCL         _plane;
    Maple3DOptionCL _options;
    double          _scalevec;

    DiscVelT _v;
    DiscPrT  _p;

    bool
    _Write_Intersection(const TetraCL&, std::ostream&, bool) const;
    bool
    _Write_Vel_Intersection(const TetraCL&, std::ostream&, bool) const;

  public:
    int gridlines;  // Maple linestyle for the grid; -1 == no grid

    MapleSolOutCL(const MultiGridCL& mg, const DiscVelT& v, const DiscPrT& p,
                  int TriLevel= -1, const PlaneCL& pl= PlaneCL() ,
                  Maple3DOptionCL opt=Maple3DOptionCL(), double scale= 1.0, int gr= 0)
      : MGOutCL(&mg), _level(TriLevel<0 ? mg.GetLastLevel() : TriLevel),
        _plane(pl), _options(opt), _scalevec(scale), _v( v), _p( p), gridlines(gr)
    {}
    virtual ~MapleSolOutCL() {}

    std::ostream& put(std::ostream&) const;
};

/// \brief Write finite element function, stored in \a v, in a file, named \a filename
void WriteFEToFile( const VecDescCL& v, MultiGridCL& mg, std::string filename, bool binary=false, const VecDescCL* lsetp=0);

/// Read a serialized finite element function from a file
/// \pre CreateNumbering of v.RowIdx must have been called before
void ReadFEFromFile( VecDescCL& v, MultiGridCL& mg, std::string filename, bool binary=false, const VecDescCL* lsetp=0);

/// \brief Write the permutation p (of an IdxDescCL), in a file, named \a filename
/// The empty permutation is treated as identity.
void WritePermutationToFile (const PermutationT& p, std::string filename);

/// \brief Read the permutation , stored in \a idx, from a file, named \a filename
/// The empty permutation is treated as identity.
void ReadPermutationFromFile (PermutationT& p, std::string filename);

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
    Assert(_level==_discsol.GetLevel(), DROPSErrCL("GeomSolOutCL::put: wrong level"), ~0);
    os << "appearance {\n-concave\nshading smooth\n}\n";
    os << "LIST {\n";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level);
          tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)) )
            continue;
//        if (GetBaryCenter(*tit)[2]>0.55) continue;

        Point3DCL Offset(0.0);
        std::vector<double> val( NumVertsC);
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
            it->Unknowns( _idx)= ++numv;

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
                    os << it->GetVertex(v)->Unknowns(_idx) << '\t';
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
         pidx= _p.GetSolution()->RowIdx->GetIdx();

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
                os << (it->GetVertex(v)->Unknowns(pidx)+1) << '\t';
        os << '\n';
    }

    return os;
}


inline double
CutEdge(const SVectorCL<3>& v0, const SVectorCL<3>& v1, const PlaneCL& pl)
// Calculates the parameter, for which the convex-combination of v0, v1 lies in pl.
// It is only called for existing, nondegenerated cuts.
{
    return (pl.r - inner_prod(pl.n, v0))/inner_prod(v1 - v0, pl.n);
}

// Intersects t with pl; possibly empty
TetraSectionT
intersect(const TetraCL& t, const PlaneCL& pl);

template <class  DiscVelT, class DiscPrT>
bool
MapleSolOutCL<DiscVelT, DiscPrT>::_Write_Vel_Intersection(const TetraCL& t, std::ostream& os, bool have_previous) const
{
    TetraSectionT section= intersect(t, this->_plane);
    if (section.empty()) return have_previous;

    const SVectorCL<3> bary= std::accumulate(section.begin(), section.end(), SVectorCL<3>() )/section.size();
    const SVectorCL<3> vel= _v.val(t, bary[0], bary[1], bary[2]);
// Maple does not like too short vectors.
    if (vel.norm()< 1.e-7) return have_previous;
    const SVectorCL<3> base= GetWorldCoord(t, bary);
    if (have_previous) os << ", ";
    os << "[[" << base[0] << ", " << base[1] << ", " << base[2] << "], ["
               <<  vel[0] << ", " <<  vel[1] << ", " <<  vel[2] << "]]";
    return true;
}

template <class  DiscVelT, class DiscPrT>
std::ostream&
MapleSolOutCL<DiscVelT, DiscPrT>::put(std::ostream& os) const
{
    os << _options.name << ":= PLOT3D(";
    if ( _options.WriteGlobalOptions(os) ) os << ", ";

    if (gridlines >= 0)
    {
        MapleMGOutCL grid(*_MG, _level, false, true, _plane, _options);
        grid._Intersect_Plot(os);
        os << ", LINESTYLE(" << gridlines << ")) ):" << std::endl;
    }

    bool have_previous= false;
    os << _options.name << "_v:= arrow([ ";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level); tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        have_previous= _Write_Vel_Intersection(*tit, os, have_previous);
    }
    os << " ]," << std::endl
       << "shape=arrow, difference=false, length=[" << _scalevec << ", relative=true]):" << std::endl;
    return os;
}


} // end of namespace DROPS

#endif
