//**************************************************************************
// File:    integrTime.h                                                   *
// Content: classes that perform time-integration steps                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 19 2001                                           *
//**************************************************************************

#ifndef DROPS_STO_INTEGRTIME_H
#define DROPS_STO_INTEGRTIME_H

#include "num/stokessolver.h"
#include "stokes/instatstokes.h"

namespace DROPS
{

template <class StokesT, class SolverT>
class InstatStokesThetaSchemeCL
//**************************************************************************
//  for solving the instationary Stokes equation of type StokesT with a
//  1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
//                       theta=1/2 -> Crank-Nicholson (Trapezregel)
//
//  Inner stationary Stokes-type problems are solved with a SolverT-solver.
//  The matrices A, B, M and the rhs b, c of the Stokes class have to be set
//  properly! After construction, SetTimeStep has to be called once. Then
//  every DoStep performs one step in time. Changing time steps require 
//  further calls to SetTimeStep.
//**************************************************************************
{
  private:
    StokesT& _Stokes;
    SolverT& _solver;
    
    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MatrixCL      _mat;               // M + theta*dt*A
    
    double _theta, _dt;
    
  public:
    InstatStokesThetaSchemeCL( StokesT& Stokes, SolverT& solver, double theta= 0.5)
        : _Stokes( Stokes), _solver( solver), _b( &Stokes.b), _old_b( new VelVecDescCL),
          _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns), 
          _theta( theta)
    { 
        _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
        _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    }

    ~InstatStokesThetaSchemeCL()
    {
        if (_old_b == &_Stokes.b)
            delete _b;
        else
            delete _old_b; 
        delete _cplM; delete _old_cplM;
    }
    
    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
        _dt= dt;
        _mat.LinComb( 1., _Stokes.M.Data, _theta*dt, _Stokes.A.Data);
    }
       
    void DoStep( VectorCL& v, VectorCL& p);
};


//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// A Poisson-problem with natural boundary-conditions for the pressure is
// solved via a SSOR-PCG-solver, a problem with the mass-matrix aswell.
// The constant k_ has to be chosen according to h and dt, see Theorem 4.1
// of the above paper.
//
// A_ is the pressure-Poisson-Matrix for natural boundary-conditions, M_ the
// pressure-mass-matrix.
//**************************************************************************
class ISPreCL
{
  private:
    DROPS::PCG_SsorCL solver_;
    DROPS::MatrixCL& A_;
    DROPS::MatrixCL& M_;
    double k_;

  public:
    ISPreCL(DROPS::MatrixCL& A_pr, DROPS::MatrixCL& M_pr,
            double k_pc, DROPS::Uint max_iter)
        :solver_( DROPS::SSORPcCL( 1.), max_iter, 1e-12),
         A_( A_pr), M_( M_pr), k_( k_pc)
    {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
};

//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// Confer ISPreCL for details. This preconditioner uses multigrid-solvers.
//**************************************************************************
class ISMGPreCL
{
  private:
    DROPS::MGDataCL& A_;
    DROPS::MGDataCL& M_;
    DROPS::Uint max_iter_;
    double k_;

  public:
    ISMGPreCL(DROPS::MGDataCL& A_pr, DROPS::MGDataCL& M_pr,
              double k_pc, DROPS::Uint max_iter)
        :A_( A_pr), M_( M_pr), max_iter_( max_iter), k_( k_pc)
    {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
};


//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// Confer ISPreCL for details regarding preconditioning of S. This
// preconditioner uses multigrid-solvers.
// It is a block-diagonal-preconditioner for Minres-solvers.
//**************************************************************************
class ISMinresMGPreCL
{
  private:
    DROPS::MGDataCL& Avel_;
    DROPS::MGDataCL& Apr_;
    DROPS::MGDataCL& Mpr_;
    DROPS::Uint iter_vel_;
    DROPS::Uint iter_prA_;
    DROPS::Uint iter_prM_;
    double tol_prA_;
    double k_;
//    DROPS::VectorCL onesMpr;
//    double vol_;

  public:
    ISMinresMGPreCL(DROPS::MGDataCL& A_vel,
                    DROPS::MGDataCL& A_pr, DROPS::MGDataCL& M_pr,
                    double k_pc, DROPS::Uint iter_vel, DROPS::Uint iter_prA,
                    DROPS::Uint iter_prM, double tol_prA)
        :Avel_( A_vel), Apr_( A_pr), Mpr_( M_pr), iter_vel_( iter_vel),
         iter_prA_( iter_prA), iter_prM_( iter_prM), tol_prA_( tol_prA), k_( k_pc)//,
//         onesMpr ( Mpr_.back().A.Data.num_rows())
    {
//        // Compute projection on constant pressure function and volume only once.
//        DROPS::VectorCL ones( 1.0, Mpr_.back().A.Data.num_rows());
//        onesMpr= Mpr_.back().A.Data*ones;
//        vol_= onesMpr*ones;
    }

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& /*A*/, const Mat& /*B*/, Vec& v, Vec& p, const Vec& b, const Vec& c) const;
};


//=================================
//     template definitions
//=================================

template <class StokesT, class SolverT>
void InstatStokesThetaSchemeCL<StokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * v;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*v + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*_old_b->Data);

    p*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);
    p/= _dt;

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
}


template <typename Mat, typename Vec>
void
ISPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
    const_cast<ISPreCL*const>( this)->solver_.Solve( A_, p, c);
//    std::cerr << "ISPreCL: A: iterations: " << solver_.GetIter()
//              << "\tresidual: " << solver_.GetResid() << '\t';
    DROPS::VectorCL p2( 0.0, p.size());
    const_cast<ISPreCL*const>( this)->solver_.SetTol(solver_.GetTol()/k_);
    const_cast<ISPreCL*const>( this)->solver_.Solve( M_, p2, c);
//    std::cerr << "M: iterations: " << solver_.GetIter()
//              << "\tresidual: " << solver_.GetResid() << '\n';
    const_cast<ISPreCL*const>( this)->solver_.SetTol(solver_.GetTol()*k_);
//    p+= k_*p2;
    axpy( k_, p2, p);
}


template <typename Mat, typename Vec>
void
ISMGPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
    DROPS::Uint sm=  2; // how many smoothing steps?
    int lvl= -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
    DROPS::SORsmoothCL smoother( omega);  // Gauss-Seidel with over-relaxation
    DROPS::SSORPcCL directpc;
    DROPS::PCG_SsorCL solver( directpc, 200, 1e-12);
    for (DROPS::Uint i=0; i<max_iter_; ++i)
        DROPS::MGM( A_.begin(), --A_.end(), p, c, smoother, sm, solver, lvl, -1);
    std::cerr << "IsMGPcCL: iterations: " << max_iter_ << '\t'
              << "residual: " << (A_.back().A.Data*p - c).norm() << '\t';
    DROPS::VectorCL p2( 0.0, p.size());
    for (DROPS::Uint i=0; i<2; ++i)
        DROPS::MGM( M_.begin(), --M_.end(), p2, c, smoother, sm, solver, lvl, -1);
    std::cerr << "M: iterations: " << 2 << '\t'
              << "residual: " << (M_.back().A.Data*p2 - c).norm() << '\n';
//    p+= k_*p2;
    axpy( k_, p2, p);
}

template <typename Mat, typename Vec>
void
ISMinresMGPreCL::Apply(const Mat& /*A*/, const Mat& /*B*/, Vec& v, Vec& p, const Vec& b, const Vec& c) const
{
    DROPS::Uint sm= 2; // how many smoothing steps?
    int lvl= -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
//    DROPS::SORsmoothCL smoother( omega);  // Gauss-Seidel with over-relaxation
    DROPS::SSORsmoothCL smoother( omega);  // Symmetric-Gauss-Seidel with over-relaxation
    DROPS::SSORPcCL directpc; DROPS::PCG_SsorCL solver( directpc, 200, 1e-12);

    for (DROPS::Uint i=0; i<iter_vel_; ++i)
        DROPS::MGM( Avel_.begin(), --Avel_.end(), v, b, smoother, sm, solver, lvl, -1);
    std::cerr << "ISMinresMGPreCL: Velocity: iterations: " << iter_vel_ << '\t';
//              << " residual: " <<  (Avel_.back().A.Data*v - b).norm() << '\t';

//    p= 0.;
//    for (DROPS::Uint i=0; i<iter_prA_; ++i) {
//        DROPS::MGM( Apr_.begin(), --Apr_.end(), p, c, smoother, sm, solver, lvl, -1);
//    }
//    std::cerr << "Pressure: iterations: " << iter_prA_ <<'\t'
//              << " residual: " <<  (Apr_.back().A.Data*p - c).norm() << '\t';
    DROPS::PCG_SsorCL solver_( DROPS::SSORPcCL( 1.), 1000, tol_prA_);
    solver_.Solve( Apr_.back().A.Data, p, c);
    std::cerr << "Pressure: iterations: " << solver_.GetIter() <<'\t'
              << " residual: " <<  solver_.GetResid() << '\t';

    DROPS::VectorCL p2( 0.0, p.size());
    for (DROPS::Uint i=0; i<iter_prM_; ++i)
        DROPS::MGM( Mpr_.begin(), --Mpr_.end(), p2, c, smoother, sm, solver, lvl, -1);
    std::cerr << "Mass: iterations: " << iter_prM_ << '\t'
              << " residual: " <<  (Mpr_.back().A.Data*p2 - c).norm() << '\n';
    
//    p+= k_*p2;
    axpy( k_, p2, p);
}

}    // end of namespace DROPS

#endif
