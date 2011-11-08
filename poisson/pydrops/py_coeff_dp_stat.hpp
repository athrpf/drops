
#ifndef __PY_COEFF_DP_STAT_HPP__
#define __PY_COEFF_DP_STAT_HPP__

#include "functors.hpp"
#include "drops_utils.hpp"
#include "py_journalist.hpp"
namespace DROPS {
  void coeff_stat(Journalist& jnlst, MassTransferBrick& brick, PoissonCoeffCL& pcl, SolutionContainer* sol_container, double tol, int maxiter);
}
#endif
