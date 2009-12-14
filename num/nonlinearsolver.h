/// \file
/// \brief nonlinear solvers, e.g. for coupling classes

#ifndef NONLINEARSOLVER_H_
#define NONLINEARSOLVER_H_

#include "num/solver.h"

namespace DROPS {

// What every iterative nonlinear solver should have
template <class FunctionT>
class NonlinearSolverBaseCL: public SolverBaseCL
{
  protected:
    FunctionT& Function_;
  public:
    NonlinearSolverBaseCL<FunctionT> (FunctionT &Function, int maxiter, double tol, bool rel= false, std::ostream* output= 0)
      : SolverBaseCL(maxiter, tol, rel, output), Function_(Function) {}

  virtual void Solve( VectorCL& v) = 0;
};

template <class FunctionT>
class FixedPointSolverCL: public NonlinearSolverBaseCL<FunctionT>
{
  private:
    typedef NonlinearSolverBaseCL<FunctionT> base_;
    using base_::Function_;
    using base_::_maxiter;
    using base_::_tol;
    using base_::_res;
    using base_::output_;

  public:
    FixedPointSolverCL<FunctionT> (FunctionT &Function, int maxiter, double tol, bool rel= false, std::ostream* output= 0)
      : base_(Function, maxiter, tol, rel, output) {}

    void Solve( VectorCL& v) {
        VectorCL v0(v);
        for (int i=0; i<_maxiter; ++i)
        {
            (*output_) << "~~~~~~~~~~~~~~~~ FP-Iter " << i+1 << '\n';
            v-= Function_.eval( v);
            if ( Function_.NoChange() || norm(VectorCL( v0-v)) < _tol)
            {
                (*output_) << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
                _res = norm(VectorCL( v0-v));
                break;
            }
        }
    }
};

template <class FunctionT>
class DampedFixedPointSolverCL: public NonlinearSolverBaseCL<FunctionT>
{
  private:
    typedef NonlinearSolverBaseCL<FunctionT> base_;
    using base_::Function_;
    using base_::_maxiter;
    using base_::_tol;
    using base_::_res;
    using base_::output_;

  public:
    DampedFixedPointSolverCL<FunctionT> (FunctionT &Function, int maxiter, double tol, bool rel= false, std::ostream* output= 0)
      : base_(Function, maxiter, tol, rel, output) {}

    void Solve( VectorCL& v) {
        VectorCL v0(v);
        double omega = 1.0;
        VectorCL Fold (v.size());
        VectorCL diff(v.size());
        for (int i=0; i<_maxiter; ++i)
        {
            (*output_) << "~~~~~~~~~~~~~~~~ FP-Iter " << i+1 << '\n';
            VectorCL tmp( Function_.eval( v));
            if (i > 0) {
                diff = tmp - Fold;
                omega *= -dot( diff, Fold)/ norm_sq( diff);
            }
            (*output_) << "DampedFixedPointSolverCL: omega = " << omega << std::endl;
            v-= omega*tmp;
            Fold = tmp;
            if ( Function_.NoChange() || norm(VectorCL( v0-v)) < _tol)
            {
                (*output_) << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
                _res = norm(VectorCL( v0-v));
                break;
            }
        }
    }
};

template <class FunctionT>
class BroydenSolverCL: public NonlinearSolverBaseCL<FunctionT>
{
  private:
    typedef NonlinearSolverBaseCL<FunctionT> base_;
    using base_::Function_;
    using base_::_maxiter;
    using base_::_tol;
    using base_::_res;
    using base_::output_;

  public:
    BroydenSolverCL<FunctionT> (FunctionT &Function, int maxiter, double tol, bool rel= false, std::ostream* output= 0)
      : base_(Function, maxiter, tol, rel, output) {}

    void Solve( VectorCL& x);
};

template <class FunctionT>
void BroydenSolverCL<FunctionT>::Solve( VectorCL& x){

    double thetamax_ = 0.45;
    double kappamax_ = 10000;

    VectorCL v0(x);

    typedef std::vector<VectorCL> VecValueT;

    VecValueT F, deltaF;
    std::vector<double> gamma;

    // eval F(x^{0})
    F.push_back(Function_.eval( x));

    double sigma0, sigma;

    sigma0 = norm_sq(F.back());

    VectorCL deltax( -F.back());

    double kappa = 1.0;

    for (int i=0; i<_maxiter; ++i)
    {
        (*output_) << "~~~~~~~~~~~~~~~~ FP-Iter " << i+1 << '\n';

        x += deltax;

        F.push_back(Function_.eval( x));

        deltaF.push_back( VectorCL( F.back() - F[i]));

        sigma = norm_sq(F.back());

        if (sigma < _tol*_tol) {
            (*output_) << "Solution found\n";
            break;
        }
        const double theta = std::sqrt(sigma/sigma0);
        if (theta >= thetamax_) {
            (*output_) << "No convergence: theta = " << theta << std::endl;
            //break;
        }
        sigma0 = sigma;

        VectorCL w( deltaF.back());
        gamma.push_back( norm_sq(w));

        kappa /= (1.0-2.0*theta);

        if (kappa >= kappamax_){
            (*output_) << "ill-conditioned update: kappa = "<< kappa << std::endl;
            //break;
        }

        VectorCL v( F.back());
        const double factor = 1.0 - dot( w, v) / gamma[i];
        v*= factor;

        for (int j=i-1; j>=0; --j) {
            const double beta = ( dot (deltaF[j], v)) / gamma[j];
            v -= beta*F[j+1];
        }
        deltax = -v;

        if ( Function_.NoChange() || norm(VectorCL( v0-v)) < _tol)
        {
            (*output_) << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            _res = norm(VectorCL( v0-v));
            break;
        }

    }
}


} // end of namespace DROPS

#endif /* NONLINEARSOLVER_H_ */
