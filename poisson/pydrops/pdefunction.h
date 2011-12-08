#ifndef PDEFUNCTION_H
#define PDEFUNCTION_H

#include "num/discretize.h"

class PdeFunction {
public:
  /// Return the function value at the given grid point.
  virtual double operator()(int nx, int ny, int nz, int nt) const=0;
  /** Set inputs to values in this PdeFunction. Return true, if the values that were passed in match those of the grid underlying this function. */
  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const=0;
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
      nx_ = P.get<int>("DomainCond.nx") * pow(2, P.get<int>("DomainCond.RefineSteps")) +1;
      ny_ = P.get<int>("DomainCond.ny") * pow(2, P.get<int>("DomainCond.RefineSteps")) +1;
      nz_ = P.get<int>("DomainCond.nz") * pow(2, P.get<int>("DomainCond.RefineSteps")) +1;
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
