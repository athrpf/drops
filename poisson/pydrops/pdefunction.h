#ifndef PDEFUNCTION_H
#define PDEFUNCTION_H

#include "misc/params.h"

class PdeFunction {
public:
  /// Return the function value at the given grid point.
  virtual double operator()(int nx, int ny, int nz, int nt) const=0;
  /** Set inputs to values in this PdeFunction. Return true, if the values that were passed in match those of the grid underlying this function. */
  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const=0;
};




#endif
