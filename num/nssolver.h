/// \file nssolver.h
/// \brief nonlinear solvers for the Navier-Stokes equation
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "num/stokessolver.h"
#include "misc/problem.h"
#ifdef _PAR
#  include "parallel/exchange.h"
#  include "parallel/logger.h"
#endif

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
    using SolverBaseCL::rel_;

  public:
    NSSolverBaseCL (NavStokesT& NS, StokesSolverBaseCL& solver, int maxiter= -1, double tol= -1.0, bool rel= false)
        : SolverBaseCL(maxiter, tol, rel), NS_( NS), solver_( solver) {}

    virtual ~NSSolverBaseCL() {}

    virtual double   GetResid ()           const { return solver_.GetResid(); }
    virtual int      GetIter  ()           const { return solver_.GetIter(); }
    StokesSolverBaseCL& GetStokesSolver () const { return solver_; }
    virtual const MLMatrixCL* GetAN()            { return &NS_.A.Data; }

    /// solves the system   A v + BT p = b
    ///                     B v        = c
    virtual void Solve (const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
        const VectorCL& b, VecDescCL& cplN, const VectorCL& c, double)
    {
        solver_.Solve( A, B, v.Data, p, b, c);
        cplN.Data= 0.;
    }
    virtual void Solve (const MLMatrixCL& A, const MLMatrixCL& B, VecDescCL& v, VectorCL& p,
        const VectorCL& b, VecDescCL& cplN, const VectorCL& c, double)
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
                const VectorCL& b, VecDescCL& cplN, const VectorCL& c, double alpha= 1.);
    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VecDescCL& v, VectorCL& p,
                const VectorCL& b, VecDescCL& cplN, const VectorCL& c, double alpha= 1.);
};


/// \brief Computes the relaxation factor in AdapFixedPtDefectCorrCL by line search.
class LineSearchPolicyCL
{
  private:
    double omega_;
    VectorCL d_,
             e_;
#ifdef _PAR
    VectorCL d_acc_,
             e_acc_;
#endif
  public:
    LineSearchPolicyCL (size_t vsize, size_t psize)
        : omega_( 1.0), d_( vsize), e_( psize)
#ifdef _PAR
          , d_acc_( vsize), e_acc_( psize)
#endif
        {}

    /// \todo (merge) More "const" in Update make it easier to see, that no parameter vectors changes ...
    template<class NavStokesT>
      void
      Update (NavStokesT&, const MatrixCL&, const MatrixCL&,
        const VecDescCL&, const VectorCL&, const VectorCL&, VecDescCL&, const VectorCL&,
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
        const VecDescCL&, const VectorCL&, const VectorCL&, VecDescCL&, const VectorCL&,
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
        const VecDescCL& v, const VectorCL& p, const VectorCL& b, VecDescCL& cplN, const VectorCL& c,
        const VectorCL& w, const VectorCL& q, double alpha);

    double RelaxFactor () const { return omega_; }
};

//=================================
//     template definitions
//=================================
template<class NavStokesT>
  inline void
  LineSearchPolicyCL::Update (NavStokesT& ns, const MatrixCL& A, const MatrixCL& B,
    const VecDescCL& v, const VectorCL& p, const VectorCL& b, VecDescCL& cplN, const VectorCL& c,
    const VectorCL& w, const VectorCL& q, double alpha)
// accumulated and non accumulated vectors:
// v - acc, p - acc, b - non-acc, cplN - non-acc, c - non-acc, w - acc, q - acc
{
#ifdef _PAR
    ExchangeCL& ExVel= ns.vel_idx.GetEx();
    ExchangeCL& ExPr = ns.pr_idx.GetEx();
    const bool useAccur=true;
#endif
    VecDescCL v_omw( v.RowIdx);
    v_omw.t = v.t;
    v_omw.Data= v.Data - omega_*w;
    ns.SetupNonlinear( ns.N.Data.GetFinest(), &v_omw, &cplN, ns.N.RowIdx->GetFinest());

    d_= A*w + alpha*(ns.N.Data.GetFinest()*w) + transp_mul( B, q);
    e_= B*w;
#ifndef _PAR
    omega_= dot( d_, VectorCL( A*v.Data + alpha*(ns.N.Data.GetFinest()*v.Data)
                     + transp_mul( B, p) - b - alpha*cplN.Data))
            + dot( e_, VectorCL( B*v.Data - c));
    omega_/= norm_sq( d_) + norm_sq( e_);
#else
    omega_ = ExVel.ParDot(d_, false,
                            VectorCL( A*v.Data + alpha*(ns.N.Data.GetFinest()*v.Data)
                                      + transp_mul( B, p) - b - alpha*cplN.Data),
                            false, useAccur, &d_acc_);
    omega_+= ExPr.ParDot(e_, false, VectorCL( B*v.Data - c), false, useAccur, &e_acc_);
    omega_/= ExVel.Norm_sq(d_acc_, true) + ExPr.Norm_sq(e_acc_, true);
#endif
}

template<class NavStokesT>
  inline void
  DeltaSquaredPolicyCL::Update (__UNUSED__ NavStokesT& ns, const MatrixCL&, const MatrixCL&,
    const VecDescCL&, const VectorCL&, const VectorCL&, VecDescCL&, const VectorCL&,
    const VectorCL& w, const VectorCL& q, double)
{
    if (firststep_) {
        w_old_= w; q_old_= q;
        firststep_ = false;
        return;
    }
    w_diff_=  w - w_old_; q_diff_= q - q_old_;
#ifndef _PAR
    omega_*= -(dot( w_diff_, w_old_) + dot( q_diff_, q_old_))
              / (norm_sq( w_diff_) + norm_sq( q_diff_));
#else
    ExchangeCL& ExVel  = ns.vel_idx.GetEx();
    ExchangeCL& ExPr   = ns.pr_idx.GetEx();
    const bool useAccur= true;
    omega_*= -(ExVel.ParDot( w_diff_, true, w_old_, true, useAccur)
               + ExPr.ParDot( q_diff_, true, q_old_, true, useAccur))
              / (ExVel.Norm_sq( w_diff_, true, useAccur) + ExPr.Norm_sq( q_diff_, true, useAccur));
#endif

    w_old_= w; q_old_= q;
}

template<class NavStokesT, class RelaxationPolicyT>
void
AdaptFixedPtDefectCorrCL<NavStokesT, RelaxationPolicyT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
    const VectorCL& b, VecDescCL& cplN, const VectorCL& c, double alpha)
{
    VectorCL d( v.Data.size()), e( p.size()),
             w( v.Data.size()), q( p.size());
    RelaxationPolicyT relax( v.Data.size(), p.size());
#ifdef _PAR
    VectorCL d_acc( v.Data.size()), e_acc( p.size()),
             w_acc( v.Data.size()), q_acc( p.size());
    ExchangeCL& ExVel= NS_.vel_idx.GetEx();
    ExchangeCL& ExPr = NS_.pr_idx.GetEx();
    const bool useAccur=true;
#endif
    double res0= 1.;
    int oseenIter= 0;
    _iter= 0;
    for(;;++_iter) { // ever
        NS_.SetupNonlinear(&NS_.N, &v, &cplN);
        //std::cout << "sup_norm : N: " << supnorm( _NS.N.Data) << std::endl;
        AN_->GetFinest().LinComb( 1., A, alpha, NS_.N.Data.GetFinest());

        // calculate defect:
        d= *AN_*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
        e= B*v.Data - c;
#ifndef _PAR
        _res= std::sqrt( norm_sq( d) + norm_sq( e) );
#else
        _res= std::sqrt( ExVel.Norm_sq(d, false, useAccur, &d_acc) + ExPr.Norm_sq(e, false, useAccur, &e_acc) );
#endif
        /// \todo(merge) Do we need this output? Or should/could we use the (*output_)?
        std::cout << _iter << ": res = " << _res << " reltol: " << this->GetRelError() << std::endl;
        if (this->GetRelError() == true && _iter == 0)
            res0= _res;
        if (_res < _tol*res0 || _iter>=_maxiter) // if absolute errors are required, res0==1.
            break;

        // solve correction:
        double outer_tol= _res*red_;
        if (outer_tol < 0.5*_tol && this->GetRelError() == false)
            outer_tol= 0.5*_tol;
        solver_.SetTol( outer_tol);
        w= 0.0; q= 0.0;
        solver_.Solve( AN_->GetFinest(), B, w, q, d, e); // solver_ should use a relative termination criterion.
        oseenIter+= solver_.GetIter();

        // calculate step length omega:
        relax.Update( NS_, A,  B, v, p, b, cplN, c, w,  q, alpha);

        // update solution:
        const double omega( relax.RelaxFactor());
        std::cout << "omega = " << omega << std::endl;
        v.Data-= omega*w;
        p     -= omega*q;
    }
      std::cout << "overall iterations of Oseen solver:\t" << oseenIter << std::endl;
}

template<class NavStokesT, class RelaxationPolicyT>
void
AdaptFixedPtDefectCorrCL<NavStokesT, RelaxationPolicyT>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B, VecDescCL& v, VectorCL& p,
    const VectorCL& b, VecDescCL& cplN, const VectorCL& c, double alpha)
{
    VectorCL d( v.Data.size()), e( p.size()),
             w( v.Data.size()), q( p.size());
    RelaxationPolicyT relax( v.Data.size(), p.size());
#ifdef _PAR
    VectorCL d_acc( v.Data.size()), e_acc( p.size()),
             w_acc( v.Data.size()), q_acc( p.size());
    ExchangeCL& ExVel= NS_.vel_idx.GetEx();
    ExchangeCL& ExPr = NS_.pr_idx.GetEx();
    const bool useAccur=true;
#endif
    double res0= 1.;
    int oseenIter= 0;

    _iter= 0;
    for(;;++_iter) { // ever
        NS_.SetupNonlinear(&NS_.N, &v, &cplN);
        //std::cout << "sup_norm : N: " << supnorm( _NS.N.Data) << std::endl;
        AN_->LinComb( 1., A, alpha, NS_.N.Data);
        // calculate defect:
        d= *AN_*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
        e= B*v.Data - c;
#ifndef _PAR
        _res= std::sqrt( norm_sq( d) + norm_sq( e) );
#else
        _res= std::sqrt( ExVel.Norm_sq(d, false, useAccur, &d_acc) + ExPr.Norm_sq(e, false, useAccur, &e_acc) );
#endif
        std::cout << _iter << ": res = " << _res << " reltol: " << this->GetRelError() << std::endl;
        if (this->GetRelError() == true && _iter == 0)
            res0= _res;
        if (_res < _tol*res0 || _iter>=_maxiter) // if absolute errors are required, res0==1.
            break;

        // solve correction:
        double outer_tol= _res*red_;
        if (outer_tol < 0.5*_tol && this->GetRelError() == false)
            outer_tol= 0.5*_tol;
        w= 0.0; q= 0.0;
        solver_.SetTol( outer_tol);
        solver_.Solve( *AN_, B, w, q, d, e); // solver_ should use a relative termination criterion.
        oseenIter+= solver_.GetIter();

        // calculate step length omega:
        relax.Update( NS_, A.GetFinest(),  B.GetFinest(), v, p, b, cplN, c, w,  q, alpha);

        // update solution:
        const double omega( relax.RelaxFactor());
        std::cout << "omega = " << omega << std::endl;
        v.Data-= omega*w;
        p     -= omega*q;
    }
    std::cout << "overall iterations of Oseen solver:\t" << oseenIter << std::endl;
}

}    // end of namespace DROPS

#endif
