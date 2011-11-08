
#include "functors.hpp"
#include "poisson/poisson.h"
#include "py_coeff_dp_stat.hpp"
namespace DROPS {
  class KDeltaSolution : public SolutionContainer
  {
  private:

    instat_scalar_fun_ptr dalpha_;
    instat_scalar_fun_ptr u_;
    double scalar_prod_;

  public:
    KDeltaSolution(instat_scalar_fun_ptr dalpha,
		   instat_scalar_fun_ptr u)
      : dalpha_(dalpha), u_(u)
    {}

    double get_scalarprod() const {return scalar_prod_;}

    virtual void set_solution(const DROPS::InstatPoissonP1CL<PoissonCoeffCL>& poisson, double t)
    {
      const VecDescCL& psi = poisson.x;
      VecDescCL kdelta;
      IdxDescCL idx= poisson.idx;
      kdelta.SetIdx(&idx);
      poisson.SetupGradSrc(kdelta, u_, dalpha_, NULL);
      scalar_prod_ = dot(kdelta.Data, psi.Data);
    }
  };

  double compute_kdelta(Journalist& jnlst, MassTransferBrick& brick, PoissonCoeffCL& pcl, instat_scalar_fun_ptr dalpha, instat_scalar_fun_ptr u, double tol, int maxiter)
  {
    KDeltaSolution* sol = new KDeltaSolution(dalpha, u);
    coeff_stat(jnlst, brick, pcl, sol, tol, maxiter);
    double retval = sol->get_scalarprod();
    delete sol;
    return retval;
  }
}
