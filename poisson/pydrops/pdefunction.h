#ifndef PDEFUNCTION_H
#define PDEFUNCTION_H

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include "misc/params.h"
#include <boost/shared_ptr.hpp>

#include "pyconnect.h"
#include "convection_diffusion.cpp"



class PdeFunction {
public:
  /// Return the function value at the given grid point.
  virtual double operator()(int nx, int ny, int nz, int nt) const=0;
  /** Set inputs to values in this PdeFunction. Return true, if the values that were passed in match those of the grid underlying this function. */
  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const=0;
};


//
class PyPdeFunction : public PdeFunction {
private:
  int nx, ny, nz, nt;
  boost::python::numeric::array data;
public:

  PyPdeFunction(boost::python::numeric::array& data_) : data(data_)
  {
    using namespace boost::python;
    tuple t = extract<tuple>(data.attr("shape"));
    assert(len(t)==4);
    nx = boost::python::extract<int>(t[0]);
    ny = boost::python::extract<int>(t[1]);
    nz = boost::python::extract<int>(t[2]);
    nt = boost::python::extract<int>(t[3]);
  }


  virtual ~PyPdeFunction(){}

  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const {
    bool retval = false;
    if (Nx==nx && Ny==ny && Nz==nz && Nt==nt) {
      retval = true;
    }
    Nx = nx; Ny = ny; Nz = nz; Nt = nt;
    return retval;
  }

  virtual double operator()(int ix, int iy, int iz, int it) const {
    assert (ix>=0 && ix<nx);
    assert (iy>=0 && iy<ny);
    assert (iz>=0 && iz<nz);
    assert (it>=0 && it<nt);
    using namespace boost::python;
    tuple t = make_tuple<int,int,int,int>(ix,iy,iz,it);
    return extract<double>(data[t]);
  }

  double at(int ix, int iy, int iz, int it) const {
    assert (ix>=0 && ix<nx);
    assert (iy>=0 && iy<ny);
    assert (iz>=0 && iz<nz);
    assert (it>=0 && it<nt);
    using namespace boost::python;
    tuple t = make_tuple<int,int,int,int>(ix,iy,iz,it);
    return extract<double>(data[t]);
  }
};


//For testing use, construct from function
class TestPdeFunction : public PdeFunction {
private:
  int nx_, ny_, nz_, nt_;
  double dx_, dy_, dz_, dt_;
  instat_scalar_fun_ptr f_;
public:
  TestPdeFunction( DROPS::ParamCL p, instat_scalar_fun_ptr f) : f_(f)
  {
      nx_ = P.get<int>("DomainCond.nx") * pow(2, P.get<int>("RefineSteps")) +1;
      ny_ = P.get<int>("DomainCond.ny") * pow(2, P.get<int>("RefineSteps")) +1;
      nz_ = P.get<int>("DomainCond.nz") * pow(2, P.get<int>("RefineSteps")) +1;
      nt_ = P.get<int>("Time.NumSteps") +1;
      
      dx_ = P.get<double>("DomainCond.lx")/(nx_ -1);
      dy_ = P.get<double>("DomainCond.ly")/(ny_ -1);
      dz_ = P.get<double>("DomainCond.lz")/(nz_ -1);
      dt_ = P.get<double>("Time.StepSize");
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
    return f(p, t);
  }
  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const {
    bool retval = false;
    if (Nx==nx && Ny==ny && Nz==nz && Nt==nt) {
      retval = true;
    }
    Nx = nx; Ny = ny; Nz = nz; Nt = nt;
    return retval;
  }    
}

#endif