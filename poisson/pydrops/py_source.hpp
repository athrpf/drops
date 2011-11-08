
#ifndef __PY_SOURCE_HPP__
#define __PY_SOURCE_HPP__

#include "functors.hpp"
#include "drops_utils.hpp"
#include "py_journalist.hpp"

void source_dp(Journalist& jnlst, MassTransferBrick& brick, PoissonCoeffCL& pcl, SolutionContainer* sol_container, DROPS::instat_scalar_fun_ptr initial, int nt, double dt, double theta, double tol, int iter, int Flag_pr, int Flag_bc, int Flag_SUPG);


#endif
