#include <string>
#include <numeric>

namespace DROPS
{

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


inline double
CutEdge(const SVectorCL<3>& v0, const SVectorCL<3>& v1, const PlaneCL& pl)
// Calculates the parameter, for which the convex-combination of v0, v1 lies in pl.
{
    return (pl.r - inner_prod(pl.n, v0))/inner_prod(v1 - v0, pl.n);
}

// Intersects t with pl; possibly empty
TetraSectionT
intersect(const TetraCL& t, const PlaneCL& pl);


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

    Maple3DOptionCL(const std::string& n= "DROPSMGPLOT", const std::string& t= "DROPS Geometry-plot",
                    const std::string& x= "X", const std::string& y= "Y", const std::string& z= "Z",
                    const std::string& axst= "BOXED")
        : name(n), title(t), xlabel(x), ylabel(y),zlabel(z), axesstyle(axst) {}

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
//    SMatrixCL<3,3>  _arrow;

    DiscVelT _v;
    DiscPrT  _p;

    bool
    _Write_Intersection(const TetraCL&, std::ostream&, bool) const;
    bool
    _Write_Vel_Intersection(const TetraCL&, std::ostream&, bool) const;
    
  public:
    int gridlines;  // Maple linestyle for the grid; -1 == no grid

    MapleSolOutCL(const MultiGridCL& mg, const DiscVelT& v, const DiscPrT& p, 
                  int TriLevel= -1, const PlaneCL& pl= PlaneCL() , Maple3DOptionCL opt=Maple3DOptionCL(), int gr= 0)
      : MGOutCL(&mg), _level(TriLevel<0 ? mg.GetLastLevel() : TriLevel), 
        _plane(pl), _options(opt), _v( v), _p( p), gridlines(gr)
    {}
    virtual ~MapleSolOutCL() {}
      
    std::ostream& put(std::ostream&) const;
};

template <class  DiscVelT, class DiscPrT>
bool
MapleSolOutCL<DiscVelT, DiscPrT>::_Write_Vel_Intersection(const TetraCL& t, std::ostream& os, bool have_previous) const
{
    TetraSectionT section= intersect(t, this->_plane);
    if (section.empty()) return have_previous;
    if (have_previous) os << ", ";

    const SVectorCL<3> bary= std::accumulate(section.begin(), section.end(), SVectorCL<3>() )/section.size();
    const SVectorCL<3> vel= _v.val(t, bary[0], bary[1], bary[2]);
    const SVectorCL<3> base= GetWorldCoord(t, bary);
    const SVectorCL<3> end= base+vel;
    os << "[[" << base[0] << ", " << base[1] << ", " << base[2] << "], ["
               << end[0] << ", " << end[1] << ", " << end[2] << "]]";
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
        os << ", LINESTYLE(" << gridlines << ")), ";
    }
    
    bool have_previous= false;
    os << "CURVES(";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level); tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        have_previous= _Write_Vel_Intersection(*tit, os, have_previous);
    }
    return os << ")):" << std::endl;
}


} // end of namespace DROPS
