#ifndef __FUNCTORS_HPP__
#define __FUNCTORS_HPP__

#include "misc/container.h"
#include "poisson/poisson.h"


class Velocity {
private:
  Velocity();
  double h_;
  double uN_;
public:
  Velocity(double height,  double uN) : h_(height), uN_(uN) {}
  DROPS::Point3DCL operator()(const DROPS::Point3DCL& p, double) const
  {
    DROPS::Point3DCL ret;
    const double d= p[1]/h_;
    ret[0] = uN_*(2-d)*d; // Nusselt
    return ret;
  }
};

/// Stabilization coefficient
class Stability_Coeff {
private:
  DROPS::instat_scalar_fun_ptr alpha;
  Velocity Vel;
  int nx_;
  double lx_;
  Stability_Coeff();
public:
  Stability_Coeff(DROPS::instat_scalar_fun_ptr alpha_in, Velocity v, int nx, double lx) : alpha(alpha_in), Vel(v), nx_(nx), lx_(lx){}

  double operator()(const DROPS::Point3DCL& p, double t) const
  {
    if (PecNum(p,t)<=1)
      return 0.0;
    else
      return h_Value()/(2.*fabs(Vel(p, t)[0]))*(1.-1./PecNum(p, t));
  }

  double PecNum(const DROPS::Point3DCL& p, double t) const
  {//Peclet Number
    double Pec=0.;
    Pec=fabs(Vel(p, t)[0])*h_Value()/(2.*alpha(p,t));
    return Pec;
  }

  double h_Value() const
  {//mesh size in flow direction
    double h = lx_/nx_;
    return h;
  }
};

// Functor that always returns a fixed value
class FixedValue {
private:
  double val_;
  FixedValue();
public:

  FixedValue(double val) : val_(val)
  {
  }

  double operator()(const DROPS::Point3DCL& p, double t) const
  {
    return val_;
  }
};

/// Functor that always returns zero (this could be replaced by a function was here mainly for testing)
class Zero {
public:
  double operator()(const DROPS::Point3DCL& p, double t) const { return 0.0; }
};

class AddFunctor {
private:
  /** Functors to be added */
  DROPS::instat_scalar_fun_ptr f1_;
  DROPS::instat_scalar_fun_ptr f2_;

  /** Linear coefficients */
  double alpha1_, alpha2_;
  /** Default constructor is not implemented. */
  AddFunctor();
public:
  AddFunctor(double alpha1, DROPS::instat_scalar_fun_ptr f1, double alpha2, DROPS::instat_scalar_fun_ptr f2)
    : f1_(f1), f2_(f2),
      alpha1_(alpha1), alpha2_(alpha2)
  {}
  double operator()(const DROPS::Point3DCL& p, double t) const { return alpha1_*f1_(p,t) + alpha2_*f2_(p,t); }
};

class PoissonCoeffCL
{
private:
  PoissonCoeffCL();
public:
  DROPS::instat_scalar_fun_ptr alpha;
  Velocity Vel;
  DROPS::instat_scalar_fun_ptr f;
  DROPS::instat_scalar_fun_ptr Sta_Coeff;

  PoissonCoeffCL(DROPS::instat_scalar_fun_ptr alpha_in, DROPS::instat_scalar_fun_ptr f_in, Velocity v, DROPS::instat_scalar_fun_ptr Sta_Coeff_in)
    :
    alpha(alpha_in),
    Vel(v),
    f(f_in),
    Sta_Coeff(Sta_Coeff_in)
  {
  }
};

class SolutionContainer {
public:
  virtual void set_solution(const DROPS::PoissonP1CL<PoissonCoeffCL>& sol, double t) =0;
};

#endif
