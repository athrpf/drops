/// \file
/// \brief nonlinear solvers for the Navier-Stokes equation
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser
// History: begin - Nov, 20 2001

#include "num/stokessolver.h"
#include "misc/problem.h"

#ifndef DROPS_NSSOLVER_H
#define DROPS_NSSOLVER_H

namespace DROPS
{

/// \brief Base class for Navier-Stokes solver. The base class version forwards all operations to the Stokes solver.
template<class NavStokesT>
class NSSolverBaseCL : public SolverBaseCL
{
  protected:
    NavStokesT& NS_;
    StokesSolverBaseCL& solver_;
    using SolverBaseCL::_iter;
    using SolverBaseCL::_maxiter;
    using SolverBaseCL::_tol;
    using SolverBaseCL::_res;

  public:
    NSSolverBaseCL (NavStokesT& NS, StokesSolverBaseCL& solver, int maxiter= -1, double tol= -1.0)
        : SolverBaseCL(maxiter, tol), NS_( NS), solver_( solver) {}

    virtual ~NSSolverBaseCL() {}

    virtual double   GetResid ()           const { return solver_.GetResid(); }
    virtual int      GetIter  ()           const { return solver_.GetIter(); }
    StokesSolverBaseCL& GetStokesSolver () const { return solver_; }
    virtual const MLMatrixCL* GetAN()            { return &NS_.A.Data; }

    /// solves the system   A v + BT p = b
    ///                     B v        = c
    virtual void Solve (const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
        VectorCL& b, VecDescCL& cplN, VectorCL& c, double)
    {
        solver_.Solve( A, B, v.Data, p, b, c);
        cplN.Data= 0.;
    }
    virtual void Solve (const MLMatrixCL& A, const MLMatrixCL& B, VecDescCL& v, VectorCL& p,
        VectorCL& b, VecDescCL& cplN, VectorCL& c, double)
    {
        solver_.Solve( A, B, v.Data, p, b, c);
        cplN.Data= 0.;
    }
};

// forward declaration
class LineSearchPolicyCL;

/// \brief adaptive fixed-point defect correction (TUREK p. 187f) for the Navier-Stokes equation.
///
/// The NS problem is of type NavStokesT, the inner problems of Stokes-type
/// are solved via a StokesSolverBaseCL-solver.
/// After the run, the NS class contains the nonlinear part N / cplN belonging
/// to the iterated solution.
/// RelaxationPolicyT defines the method, by which the relaxation factor omega is computed.
/// The default LineSearchPolicyCL corresponds to the method in Turek's book.
template <class NavStokesT, class RelaxationPolicyT= LineSearchPolicyCL>
class AdaptFixedPtDefectCorrCL : public NSSolverBaseCL<NavStokesT>
{
  private:
    typedef NSSolverBaseCL<NavStokesT> base_;
    using base_::NS_;
    using base_::solver_;
    using base_::_iter;
    using base_::_maxiter;
    using base_::_tol;
    using base_::_res;

    MLMatrixCL* AN_;

    double      red_;
    bool        adap_;

  public:
    AdaptFixedPtDefectCorrCL( NavStokesT& NS, StokesSolverBaseCL& solver, int maxiter,
                              double tol, double reduction= 0.1, bool adap=true)
        : base_( NS, solver, maxiter, tol), AN_( new MLMatrixCL( NS.vel_idx.size()) ), red_( reduction), adap_( adap) { }

    ~AdaptFixedPtDefectCorrCL() { delete AN_; }

    void SetReduction( double red) { red_= red; }
    double   GetResid ()         const { return _res; }
    int      GetIter  ()         const { return _iter; }

    const MLMatrixCL* GetAN()          { return AN_; }

    /// solves the system   [A + alpha*N] v + BT p = b + alpha*cplN
    ///                                 B v        = c
    /// (param. alpha is used for time integr. schemes)
    void Solve( const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
                VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha= 1.);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VecDescCL& v, VectorCL& p,
                VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha= 1.);
};


/// \brief Computes the relaxation factor in AdapFixedPtDefectCorrCL by line search.
class LineSearchPolicyCL
{
  private:
    double omega_;
    VectorCL d_,
             e_;

  public:
    LineSearchPolicyCL (size_t vsize, size_t psize)
        : omega_( 1.0), d_( vsize), e_( psize) {}

    template<class NavStokesT>
      void
      Update (NavStokesT&, const MatrixCL&, const MatrixCL&,
        const VecDescCL&, VectorCL&, VectorCL&, VecDescCL&, VectorCL&,
        const VectorCL&, const VectorCL&, double);

    double RelaxFactor () const { return omega_; }
};

/// \brief Always returns 1 as relaxation factor in AdapFixedPtDefectCorrCL.
class FixedPolicyCL
{
  public:
    FixedPolicyCL (size_t, size_t) {}

    template<class NavStokesT>
      void
      Update (NavStokesT&, const MatrixCL&, const MatrixCL&,
        const VecDescCL&, VectorCL&, VectorCL&, VecDescCL&, VectorCL&,
        const VectorCL&, const VectorCL&, double) {}

    double RelaxFactor () const { return 1.0; }
};

/// \brief Compute the relaxation factor in AdapFixedPtDefectCorrCL by Aitken's delta-squared method.
///
/// This vector version of classical delta-squared concergence-acceleration computes the
/// relaxation factor in span{ (w, q)^T}.
class DeltaSquaredPolicyCL
{
  private:
    bool firststep_;
    double omega_;
    VectorCL w_old_,
             q_old_,
             w_diff_,
             q_diff_;

  public:
    DeltaSquaredPolicyCL (size_t vsize, size_t psize)
        : firststep_( true), omega_( 1.0),  w_old_( vsize), q_old_( psize),
          w_diff_( vsize), q_diff_( psize) {}

    template<class NavStokesT>
      void
      Update (NavStokesT& ns, const MatrixCL& A, const MatrixCL& B,
        const VecDescCL& v, VectorCL& p, VectorCL& b, VecDescCL& cplN, VectorCL& c,
        const VectorCL& w, const VectorCL& q, double alpha);

    double RelaxFactor () const { return omega_; }
};

//=================================
//     template definitions
//=================================
template<class NavStokesT>
  inline void
  LineSearchPolicyCL::Update (NavStokesT& ns, const MatrixCL& A, const MatrixCL& B,
    const VecDescCL& v, VectorCL& p, VectorCL& b, VecDescCL& cplN, VectorCL& c,
    const VectorCL& w, const VectorCL& q, double alpha)
{
    VecDescCL v_omw( v.RowIdx);
    v_omw.Data= v.Data - omega_*w;
    ns.SetupNonlinear( ns.N.Data.GetFinest(), &v_omw, &cplN, ns.N.RowIdx->GetFinest());

    d_= A*w + alpha*(ns.N.Data.GetFinest()*w) + transp_mul( B, q);
    e_= B*w;
    omega_= dot( d_, VectorCL( A*v.Data + alpha*(ns.N.Data.GetFinest()*v.Data)
    + transp_mul( B, p) - b - alpha*cplN.Data))
    + dot( e_, VectorCL( B*v.Data - c));
    omega_/= norm_sq( d_) + norm_sq( e_);
}

template<class NavStokesT>
  inline void
  DeltaSquaredPolicyCL::Update (NavStokesT&, const MatrixCL&, const MatrixCL&,
    const VecDescCL&, VectorCL&, VectorCL&, VecDescCL&, VectorCL&,
    const VectorCL& w, const VectorCL& q, double)
{
    if (firststep_) {
        w_old_= w; q_old_= q;
        firststep_ = false;
        return;
    }
    w_diff_=  w - w_old_; q_diff_= q - q_old_;
    omega_*= -(dot( w_diff_, w_old_) + dot( q_diff_, q_old_))
              / (norm_sq( w_diff_) + norm_sq( q_diff_));

    w_old_= w; q_old_= q;
}

template<class NavStokesT, class RelaxationPolicyT>
void
AdaptFixedPtDefectCorrCL<NavStokesT, RelaxationPolicyT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
    VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha)
{
    VectorCL d( v.Data.size()), e( p.size()),
             w( v.Data.size()), q( p.size());
    RelaxationPolicyT relax( v.Data.size(), p.size());

    _iter= 0;
    for(;;++_iter) { // ever
        NS_.SetupNonlinear(&NS_.N, &v, &cplN);
        //std::cerr << "sup_norm : N: " << supnorm( _NS.N.Data) << std::endl;
        AN_->GetFinest().LinComb( 1., A, alpha, NS_.N.Data.GetFinest());

        // calculate defect:
        d= *AN_*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
        e= B*v.Data - c;
        std::cerr << _iter << ": res = " << (_res= std::sqrt( norm_sq( d) + norm_sq( e) ) ) << std::endl;
        //if (_iter == 0) std::cerr << "new tol: " << (_tol= std::min( 0.1*_res, 5e-10)) << '\n';
        if (_res < _tol || _iter>=_maxiter)
            break;

        // solve correction:
        double outer_tol= _res*red_;
        if (outer_tol < 0.5*_tol) outer_tol= 0.5*_tol;
        w= 0.0; q= 0.0;
        solver_.SetTol( outer_tol);
        solver_.Solve( AN_->GetFinest(), B, w, q, d, e); // solver_ should use a relative termination criterion.

        // calculate step length omega:
        relax.Update( NS_, A,  B, v, p, b, cplN, c, w,  q, alpha);

        // update solution:
        const double omega( relax.RelaxFactor());
        std::cerr << "omega = " << omega << std::endl;
        v.Data-= omega*w;
        p     -= omega*q;
    }
}

template<class NavStokesT, class RelaxationPolicyT>
void
AdaptFixedPtDefectCorrCL<NavStokesT, RelaxationPolicyT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VecDescCL& v, VectorCL& p,
    VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha)
{
    VectorCL d( v.Data.size()), e( p.size()),
             w( v.Data.size()), q( p.size());
    RelaxationPolicyT relax( v.Data.size(), p.size());

    _iter= 0;
    for(;;++_iter) { // ever
        NS_.SetupNonlinear(&NS_.N, &v, &cplN);
        //std::cerr << "sup_norm : N: " << supnorm( _NS.N.Data) << std::endl;
        AN_->LinComb( 1., A, alpha, NS_.N.Data);
        // calculate defect:
        d= *AN_*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
        e= B*v.Data - c;
        std::cerr << _iter << ": res = " << (_res= std::sqrt( norm_sq( d) + norm_sq( e) ) ) << std::endl;
        //if (_iter == 0) std::cerr << "new tol: " << (_tol= std::min( 0.1*_res, 5e-10)) << '\n';
        if (_res < _tol || _iter>=_maxiter)
            break;

        // solve correction:
        double outer_tol= _res*red_;
        if (outer_tol < 0.5*_tol) outer_tol= 0.5*_tol;
        w= 0.0; q= 0.0;
        solver_.SetTol( outer_tol);
        solver_.Solve( *AN_, B, w, q, d, e); // solver_ should use a relative termination criterion.

        // calculate step length omega:
        relax.Update( NS_, A.GetFinest(),  B.GetFinest(), v, p, b, cplN, c, w,  q, alpha);

        // update solution:
        const double omega( relax.RelaxFactor());
        std::cerr << "omega = " << omega << std::endl;
        v.Data-= omega*w;
        p     -= omega*q;
    }
}

}    // end of namespace DROPS

#endif
