
#ifndef __PY_COEFF_GRADIENT_HPP__
#define __PY_COEFF_GRADIENT_HPP__

#include "functors.hpp"
/*
PyObject* drops_gradient_stat(PyObject *self, PyObject *args);

PyObject* drops_sensitivity_l2_stat(PyObject *self, PyObject *args);
*/

/// Coefficient class for the stationary coefficient inverse problem
class CoeffStatGradientCoeff
{
private:
  CoeffStatGradientCoeff();
public:
  FixedValue alpha;
  Velocity Vel;
  Zero Sta_Coeff;
  Zero f;

  CoeffStatGradientCoeff(bool h1_gradient) : alpha(h1_gradient ? 1.0 : 0.0), Vel(1.,0.){}
};

#endif
