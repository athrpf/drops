
#ifndef __PY_KDELTA_PSI_HPP__
#define __PY_KDELTA_PSI_HPP__

namespace DROPS
{
  double compute_kdelta(Journalist& jnlst, MassTransferBrick& brick, PoissonCoeffCL& pcl, instat_scalar_fun_ptr dalpha, instat_scalar_fun_ptr u, double tol, int maxiter);
}
#endif
