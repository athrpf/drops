
#ifndef __PY_PARAM_HPP__
#define __PY_PARAM_HPP__
class ParamCL
{
public:
  double lx_, ly_, lz_;
  int    nx_, ny_, nz_;                                                      // n=number of intervals
  int flag_pr_; // adjoint problem?

  /*
  ParamCL()
    : lx_(180), ly_(0.8), lz_(0.3), nx_(8), ny_(2), nz_(2), nt_(50),                 // in mm
      dt_(0.02), hflux_(350), rho_(912), cp_(1540), a_mol_(0.118/1540/912), a_(1.0),// in SI
      flag_pr_(false)
      {}*/
};

class InstatParamCL : virtual public ParamCL
{
public:
  double dt_;
  int nt_;
};

class DiffusionParamCL : virtual public ParamCL
{
public:
  double a_mol_;
};

class SourceParamCL : virtual public InstatParamCL, virtual public DiffusionParamCL
{
public:
  double uN_;
};

class CoeffParamCL : virtual public DiffusionParamCL
{

};

#endif
