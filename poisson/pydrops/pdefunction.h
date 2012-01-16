#ifndef PDEFUNCTION_H
#define PDEFUNCTION_H

#include "num/discretize.h"

#include <boost/shared_ptr.hpp>

class PdeFunction {
public:
  typedef boost::shared_ptr<PdeFunction> Ptr;
  typedef boost::shared_ptr<const PdeFunction> ConstPtr;
  /// Return the function value at the given grid point.
  virtual double operator()(int nx, int ny, int nz, int nt) const=0;
  /** Set inputs to number of grid points in this PdeFunction.

      Returns true if Nx, Ny,... that were passed in match the number of gridpoints underlying this function. */
  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const=0;
};

typedef std::pair<double, double> d_pair;
typedef std::pair<double, d_pair> cmp_key;
typedef std::map<cmp_key, DROPS::FaceCL*>  FACE_MAP;
typedef std::map<cmp_key, DROPS::TetraCL*> TETRA_MAP;
typedef FACE_MAP::const_iterator face_it;
typedef TETRA_MAP::const_iterator tetra_it;

inline int rd( double d) { return static_cast<int>( d+0.5); }                   // rounding

double rnd(double d) {// rounding four digits after comma
  int i_d = (int)(d*10000);
  double d_d = (double)i_d/10000;
  return d_d;
}

class GridFunction {
public:
  typedef boost::shared_ptr<GridFunction> Ptr;
  typedef boost::shared_ptr<const GridFunction> ConstPtr;
  GridFunction(bool adjoint_) : adjoint(adjoint_) {}
  virtual bool get_indices(const DROPS::Point3DCL& p,double t,int& ix,int& iy,int& iz,int& it)const =0;
  virtual void get_barycenter_indices(const DROPS::Point3DCL& p,double t,int& ix,int& iy,int& iz,int& it,int k) const=0;
protected:
  bool adjoint;
private:
  GridFunction();
};

class VolumeGridFunction : public GridFunction {
 public:
  /** Constructor
   *
   *  Nt: Number of grid points.
   *  adjoint: if true, the time steps will be inverted.
   */
 VolumeGridFunction(int Nt_, double dx,double dy,double dz,double dt, const TETRA_MAP* tetra_map_, bool adjoint_) : GridFunction(adjoint_), Nt(Nt_), tetra_map(tetra_map_), dx_(dx), dy_(dy), dz_(dz), dt_(dt) { if (dt_==0) dt_=1.;}

  void GetNum(const DROPS::Point3DCL& p, double t, int& ix, int& iy, int& iz, int& it) const {
    ix = rd(p[0]/dx_);
    iy = rd(p[1]/dy_);
    iz = rd(p[2]/dz_);
    if (!adjoint)
      it = rd(t/dt_);
    else
      it = Nt-rd(t/dt_)-1;
  }

  /// Returns true if p,t is a non-barycentric (a true gridpoint).
  /// If true, ix, iy, ... are set to the indices of that grid point
  virtual bool get_indices(const DROPS::Point3DCL& p,double t,int& ix,int& iy,int& iz,int& it) const {
    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    tetra_it tetra= tetra_map->find(key);
    if (tetra == tetra_map->end()) {//non-barycenter
      GetNum(p,t,ix,iy,iz,it);
      return true;
    }
    return false;
  }

  /// Returns the indices for the barycentric coordinate p.
  virtual void get_barycenter_indices(const DROPS::Point3DCL& p,double t,int& ix,int& iy,int& iz,int& it,int k) const {
    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    tetra_it tetra= tetra_map->find(key);
    assert(tetra!=tetra_map->end());
    GetNum(tetra->second->GetVertex(k)->GetCoord(),t,ix,iy,iz,it);
  }
private:
  VolumeGridFunction();
  int Nt;
  const TETRA_MAP* tetra_map;
  double dx_, dy_, dz_, dt_;
};

class SurfaceGridFunction : public GridFunction {
public:
  SurfaceGridFunction(int Nt_, double dx, double dy, double dz, double dt, const FACE_MAP* face_map_, int surface_index_, bool adjoint_) :
    GridFunction(adjoint_), Nt(Nt_),
      face_map(face_map_), dx_(dx), dy_(dy), dz_(dz), dt_(dt), surface_index(surface_index_) {if(dt_==0) dt_=1.;}

  void GetNum(const DROPS::Point3DCL& p, double t, int& ix, int& iy, int& iz, int& it) const {
    ix = rd(p[0]/dx_);
    iy = rd(p[1]/dy_);
    iz = rd(p[2]/dz_);
    it = rd(t/dt_);
    if (surface_index==0 || surface_index==1) {ix=0;}
    else if (surface_index==2 || surface_index==3) {iy=0;}
    else if (surface_index==4 || surface_index==5) {iz=0;}
  }

  virtual bool get_indices(const DROPS::Point3DCL& p, double t, int& ix, int& iy, int& iz, int& it) const {
    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    face_it face = face_map->find(key);
    if (face == face_map->end()) {//non-barycenter
      GetNum(p,t,ix,iy,iz,it);
      return true;
    }
    return false;
  }

  virtual void get_barycenter_indices(const DROPS::Point3DCL&, double, int&, int&, int&, int&, int) const {
    assert(false); // all boundary functions are dirichlet...
  }
private:
  SurfaceGridFunction();
  int Nt;
  const FACE_MAP* face_map;
  double dx_, dy_, dz_, dt_;
  int surface_index;
};

class DropsFunction {
 public:
  typedef boost::shared_ptr<DropsFunction> Ptr;
  typedef boost::shared_ptr<const DropsFunction> ConstPtr;
  DropsFunction(PdeFunction::ConstPtr f_, GridFunction::ConstPtr g_, int n_) : f(f_), g(g_), n(n_){}

  double operator()(const DROPS::Point3DCL& p, double t)
  {
    int ix, iy, iz, it;
    bool non_barycenter = g->get_indices(p,t,ix,iy,iz,it);
    if (non_barycenter) {
      return (*f)(ix,iy,iz,it);
    } else {
      double retval = 0.;
      for (int k=0; k<n; ++k) {
	g->get_barycenter_indices(p,t,ix,iy,iz,it,k);
	retval += 1./n*(*f)(ix,iy,iz,it);
      }
      return retval;
    }
  }

  PdeFunction::ConstPtr get_pdefunction() const { return f;}

 private:
  PdeFunction::ConstPtr f;
  GridFunction::ConstPtr g;
  /// 4 for volumetric grid point, 3 for surface point
  int n;
};

/**
 *  Class for interpolating a function for scalar products.
 *
 *  The main difference to DropsFunction is that this class assumes that the
 *  function is only evaluated on grids.
 */
class DropsScalarProdFunction {
public:
  typedef boost::shared_ptr<DropsScalarProdFunction> Ptr;
  typedef boost::shared_ptr<const DropsScalarProdFunction> ConstPtr;
  DropsScalarProdFunction(PdeFunction::ConstPtr pdefun_, double dx_, double dy_, double dz_, double dt_) : pdefun(pdefun_), dx(dx_), dy(dy_), dz(dz_), dt(dt_) {}
  void getnum(const DROPS::Point3DCL& p, double t, int& ix, int& iy, int& iz, int& it) const
  {
    ix=rd(p[0]/dx); iy=rd(p[1]/dy); iz=rd(p[2]/dz); it=rd(t/dt);
  }
  double operator()(const DROPS::Point3DCL& p, double t) const
  {
    int ix, iy, iz, it;
    getnum(p, t, ix, iy, iz, it);
    //std::cout << "dx: " << dx << "dy: " << dy << "dz: " << dz << "dt: " << dt << std::endl;
    //std::cout << "ix: " << ix << "iy: " << iy << "iz: " << iz << "it: " << it << std::endl;
    return pdefun->operator()(ix, iy, iz, it);
  }
private:
  DropsScalarProdFunction();
  PdeFunction::ConstPtr pdefun;
  double dx, dy, dz, dt;
};


//For testing use, construct from function
class TestPdeFunction : public PdeFunction {

private:
  int nx_, ny_, nz_, nt_;
  double dx_, dy_, dz_, dt_;
  DROPS::instat_scalar_fun_ptr f_;
public:
  TestPdeFunction( DROPS::ParamCL P, DROPS::instat_scalar_fun_ptr f) : f_(f)
  {


    int refinesteps;
    double lx, ly, lz;
    int nx, ny, nz;
    refinesteps= P.get<int>("DomainCond.RefineSteps");
    std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
    size_t idx_;
    while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx_]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> lx >> ly >> lz >> nx >> ny >> nz;
    nx_ = nx * pow (2, refinesteps)+1;
    ny_ = ny * pow (2, refinesteps)+1;
    nz_ = nz * pow (2, refinesteps)+1;
    dx_= lx/(nx_-1);
    dy_= ly/(ny_-1);
    dz_= lz/(nz_-1);
    dt_ = P.get<double>("Time.StepSize");
    nt_ = P.get<double>("Time.NumSteps")+1;
  }

  virtual double operator()(int ix, int iy, int iz, int it) const {
    assert (ix>=0 && ix<nx_);
    assert (iy>=0 && iy<ny_);
    assert (iz>=0 && iz<nz_);
    assert (it>=0 && it<nt_);
    DROPS::Point3DCL p;
    double t;
    t    = it * dt_;
    p[0] = ix * dx_;
    p[1] = iy * dy_;
    p[2] = iz * dz_;
    return f_(p, t);
  }
  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const {
    bool retval = false;
    if (Nx==nx_ && Ny==ny_ && Nz==nz_ && Nt==nt_) {
      retval = true;
    }
    Nx = nx_; Ny = ny_; Nz = nz_; Nt = nt_;
    return retval;
  }
};

#endif
