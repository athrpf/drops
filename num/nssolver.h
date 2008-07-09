/// \file
/// \brief nonlinear solvers for the Navier-Stokes equation

//**************************************************************************
// Author:  Sven Gross, Joerg Grande, Volker Reichelt, Patrick Esser       *
//          IGPM RWTH Aachen                                               *
// Version: 0.2                                                            *
// History: begin - Nov, 20 2001                                           *
//**************************************************************************

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

    virtual double   GetResid ()         const { return solver_.GetResid(); }
    virtual int      GetIter  ()         const { return solver_.GetIter(); }
    StokesSolverBaseCL& GetStokesSolver () const { return solver_; }
    virtual MatrixCL& GetAN()                  { return NS_.A.Data; }

    /// solves the system   A v + BT p = b
    ///                     B v        = c
    virtual void Solve (const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
        VectorCL& b, VecDescCL& cplN, VectorCL& c, double)
    {
        solver_.Solve( A, B, v.Data, p, b, c);
        cplN.Data= 0.;
    }
};


/// \brief adaptive fixedpoint defect correction (TUREK p. 187f) for the Navier-Stokes equation.
///
/// The NS problem is of type NavStokesT, the inner problems of Stokes-type
/// are solved via a StokesSolverBaseCL-solver.
/// After the run, the NS class contains the nonlinear part N / cplN belonging
/// to the iterated solution.
template <class NavStokesT>
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

    MatrixCL AN_;
    const VectorCL* basevel_;

    double      red_;
    bool        adap_;

  public:
    AdaptFixedPtDefectCorrCL( NavStokesT& NS, StokesSolverBaseCL& solver, int maxiter,
                              double tol, double reduction= 0.1, bool adap=true)
        : base_( NS, solver, maxiter, tol), basevel_( 0), red_( reduction), adap_( adap) {}

    ~AdaptFixedPtDefectCorrCL() {}

    void SetReduction( double red) { red_= red; }

    void SetBaseVel( const VectorCL* bv) { basevel_= bv; }

    MatrixCL& GetAN()          { return AN_; }

    /// solves the system   [A + alpha*N] v + BT p = b + alpha*cplN
    ///                                 B v        = c
    /// (param. alpha is used for time integr. schemes)
    void Solve( const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
                VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha= 1.);
};

//=================================
//     template definitions
//=================================

template<class NavStokesT>
void
AdaptFixedPtDefectCorrCL<NavStokesT>::Solve(
    const MatrixCL& A, const MatrixCL& B, VecDescCL& v, VectorCL& p,
    VectorCL& b, VecDescCL& cplN, VectorCL& c, double alpha)
{
        VecDescCL v_omw( v.RowIdx);
        VectorCL d( v.Data.size()), e( p.size()),
                 w( v.Data.size()), q( p.size());
        double omega= 1.; // initial value (no damping)
        _iter= 0;
        for(;;++_iter) // ever
        {
            if (basevel_ != 0) v.Data= *basevel_ - v.Data;
            NS_.SetupNonlinear(&NS_.N, &v, &cplN);
            if (basevel_ != 0) v.Data= *basevel_ - v.Data;
            //std::cerr << "sup_norm : N: " << supnorm( _NS.N.Data) << std::endl;
            AN_.LinComb( 1., A, alpha, NS_.N.Data);

            // calculate defect:
            d= AN_*v.Data + transp_mul( B, p) - b - alpha*cplN.Data;
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
            solver_.Solve( AN_, B, w, q, d, e); // solver_ should use a relative termination criterion.

            // calculate step length omega:
            if (adap_) {
                v_omw.Data= v.Data - omega*w;
                if (basevel_ != 0) v_omw.Data= *basevel_ - v_omw.Data;
                NS_.SetupNonlinear( &NS_.N, &v_omw, &cplN);
                if (basevel_ != 0) v_omw.Data= *basevel_ - v_omw.Data;

                d= A*w + alpha*(NS_.N.Data*w) + transp_mul( B, q);
                e= B*w;
                omega= dot( d, VectorCL( A*v.Data + NS_.N.Data*VectorCL( alpha*v.Data)
                    + transp_mul( B, p) - b - alpha*cplN.Data))
                    + dot( e, VectorCL( B*v.Data - c));
                omega/= norm_sq( d) + norm_sq( e);
                std::cerr << "omega = " << omega << std::endl;
            }

            // update solution:
            v.Data-= omega*w;
            p     -= omega*q;
        }
}

}    // end of namespace DROPS

#endif
