/// \file integrTime.h
/// \brief classes that perform time-integration steps
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_STO_INTEGRTIME_H
#define DROPS_STO_INTEGRTIME_H

#include "stokes/stokes.h"
#include "num/MGsolver.h"
#include "num/spblockmat.h"
#ifdef _PAR
# include "num/parprecond.h"
# include "num/parsolver.h"
# include "misc/problem.h"
#endif

namespace DROPS
{

template< class StokesT, class SolverT>
class TimeDiscStokesCL
{
  protected:
    StokesT& _Stokes;
    SolverT& _solver;

    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MLMatrixCL    _mat;               // M + theta*dt*A

    double _theta, _dt;

  public:
    TimeDiscStokesCL( StokesT& Stokes, SolverT& solver, double theta= 0.5)
        : _Stokes( Stokes), _solver( solver), _b( &Stokes.b), _old_b( new VelVecDescCL),
        _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns()),
        _theta( theta)
        {
            _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
            _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
        };

    virtual ~TimeDiscStokesCL()
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

    /// \brief Get constant reference on solver
    const SolverT& GetSolver() const { return _solver; }
    /// \brief Get reference on solver
    SolverT& GetSolver()       { return _solver; }

    virtual void SetTimeStep( double dt)= 0;

    virtual void SetTimeStep( double dt, double theta)= 0;

    virtual void DoStep( VectorCL& v, VectorCL& p)= 0;
};

template < class StokesT, class SolverT>
class InstatStokesThetaSchemeCL : public TimeDiscStokesCL< StokesT, SolverT>
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
	typedef TimeDiscStokesCL< StokesT, SolverT> base_;
    using base_:: _Stokes;
    using base_:: _solver;

    using base_:: _b;        using base_::_old_b;        // rhs + couplings with poisson matrix A
    using base_:: _cplM;     using base_::_old_cplM;  // couplings with mass matrix M
    using base_:: _rhs;
    using base_:: _mat;               // M + theta*dt*A

    using base_:: _theta;    using base_:: _dt;

  public:
    InstatStokesThetaSchemeCL( StokesT& Stokes, SolverT& solver, double theta= 0.5)
        :base_(Stokes, solver, theta){};

    ~InstatStokesThetaSchemeCL(){}

    void SetTimeStep( double dt)
    {
        _dt= dt;
        _mat.LinComb( 1./dt, _Stokes.M.Data, _theta, _Stokes.A.Data);
    }

    void SetTimeStep( double dt, double theta)
    {
        _dt= dt;
        _mat.LinComb( 1./dt, _Stokes.M.Data, _theta, _Stokes.A.Data);
         _theta= theta;
    }

    void DoStep( VectorCL& v, VectorCL& p);
};

template< template<class, class> class BaseMethod, class StokesT, class SolverT>
class StokesFracStepSchemeCL : public BaseMethod<StokesT, SolverT>
{
  private:
    static const double facdt_[3];
    static const double theta_[3];

    typedef BaseMethod<StokesT, SolverT> base_;

    double dt3_;
    int step_;



  public:
    StokesFracStepSchemeCL( StokesT& Stokes, SolverT& solver, int step = -1)
    	   : base_( Stokes, solver), step_((step >= 0) ? step%3 : 0) {}
    double GetSubTimeStep() const { return facdt_[step_]*dt3_; }
    double GetSubTheta()    const { return theta_[step_]; }
    int    GetSubStep()     const { return step_; }

    void SetTimeStep (double dt) { // overwrites baseclass-version
        dt3_= dt;
    }

    void DoSubStep( VectorCL& v, VectorCL& p) {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fractional Step Method: Substep " << step_ << '\n';
        base_::SetTimeStep( GetSubTimeStep(), GetSubTheta());
        base_::DoStep( v, p);
        step_= (step_ + 1)%3;
    }

    void DoStep( VectorCL& v, VectorCL& p) {
        DoSubStep( v, p);
        DoSubStep( v, p);
        DoSubStep( v, p);
    }
};

template < template<class, class> class BaseMethod, class StokesT, class SolverT>
const double StokesFracStepSchemeCL<BaseMethod, StokesT, SolverT>::facdt_[3]
  = { 1.0 - std::sqrt( 0.5), std::sqrt( 2.0) - 1.0, 1.0 - std::sqrt( 0.5) };

template < template<class, class> class BaseMethod, class StokesT, class SolverT>
const double StokesFracStepSchemeCL<BaseMethod, StokesT, SolverT>::theta_[3]
  = { 2.0 - std::sqrt( 2.0), std::sqrt( 2.0) - 1.0, 2.0 - std::sqrt( 2.0) };


class SchurPreBaseCL
{
  protected:
    double kA_,   ///< scaling factor for pressure stiffness matrix or equivalent
           kM_;   ///< scaling factor for pressure mass matrix
           
  public:
    SchurPreBaseCL( double kA, double kM) : kA_( kA), kM_( kM) {}
    void SetWeights( double kA, double kM) { kA_ = kA; kM_ = kM; }
};

//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// A Poisson-problem with natural boundary-conditions for the pressure is
// solved via 1 SSOR-step, a problem with the mass-matrix aswell.
// The constants kA_, kM_ have to be chosen according to h and dt, see Theorem 4.1
// of the above paper.
// kA_ = theta/Re and kM_ = 1/dt will do a good job,
// where Re is proportional to the ratio density/viscosity.
//
// A_ is the pressure-Poisson-Matrix for natural boundary-conditions, M_ the
// pressure-mass-matrix.
//**************************************************************************
class ISPreCL : public SchurPreBaseCL
{
  private:
    MatrixCL& A_;
    MatrixCL& M_;
    SSORPcCL  ssor_;

  public:
    ISPreCL( MatrixCL& A_pr, MatrixCL& M_pr,
        double kA= 0., double kM= 1., double om= 1.)
        : SchurPreBaseCL( kA, kM), A_( A_pr), M_( M_pr), ssor_( om)  {}
    ISPreCL( MLMatrixCL& A_pr, MLMatrixCL& M_pr,
             double kA= 0., double kM= 1., double om= 1.)
    : SchurPreBaseCL( kA, kM), A_( A_pr.GetFinest()), M_( M_pr.GetFinest()), ssor_( om)  {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
};


//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// Confer ISPreCL for details. This preconditioner uses a few CG-iterations
// to solve the linear systems.
// It is well suited for InexactUzawa-Solvers.
//**************************************************************************
#ifndef _PAR
template <class SolverT>
class ISNonlinearPreCL : public SchurPreBaseCL
{
  private:
    MatrixCL&  A_;
    MatrixCL&  M_;
    SolverT&   solver_;

  public:
    ISNonlinearPreCL(SolverT& solver, MatrixCL& A_pr, MatrixCL& M_pr,
        double kA= 0., double kM= 1.)
        : SchurPreBaseCL( kA, kM), A_( A_pr), M_( M_pr),
          solver_( solver)  {}

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;
};
#else
template <typename ASolverT, typename MSolverT>
class ISNonlinearPreCL : public SchurPreBaseCL
{
  private:
    MatrixCL&  A_;
    MatrixCL&  M_;
    mutable    ASolverT& Asolver_;
    mutable    MSolverT& Msolver_;
    mutable typename ASolverT::PrecondT PcA_;
    mutable typename MSolverT::PrecondT PcM_;

  public:
    ISNonlinearPreCL(ASolverT& Asolver, MSolverT& Msolver, MatrixCL& A_pr, MatrixCL& M_pr,
        double kA= 0., double kM= 1.)
        : SchurPreBaseCL( kA, kM), A_( A_pr), M_( M_pr),
          Asolver_(Asolver), Msolver_(Msolver), PcA_(Asolver_.GetPC()), PcM_(Msolver_.GetPC())  {}

    /// \brief Apply preconditioner
    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;

    /// \brief preconditionied vector is accumulated after "Apply"
    inline bool RetAcc() const {
        return true;
    }

    /// \brief Check if preconditioners needs diagonal of the matrices
    inline bool NeedDiag() const {
        return Asolver_.GetPC().NeedDiag() && Msolver_.GetPC().NeedDiag();
    }

    /// \brief Set diagonal of the preconditioner of the solvers (the matrices are known by this class)
    template <typename Mat>
    void SetDiag(const Mat&)
    {
        PcA_.SetDiag(A_);
        PcM_.SetDiag(M_);
    }


};
#endif

//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// cf. "Iterative Techniques For Time Dependent Stokes Problems",
//     James H. Bramble, Joseph E. Pasciak, January 1994
//
// Confer ISPreCL for details regarding preconditioning of S. This
// preconditioner uses multigrid-solvers.
//
//**************************************************************************
class ISMGPreCL : public SchurPreBaseCL
{
  private:
    const Uint sm; // how many smoothing steps?
    const int lvl; // how many levels? (-1=all)
    const double omega; // relaxation parameter for smoother
    SSORsmoothCL smoother;  // Symmetric-Gauss-Seidel with over-relaxation
    SSORPcCL directpc;
    mutable PCG_SsorCL solver;

    DROPS::MLMatrixCL& Apr_;
    DROPS::MLMatrixCL& Mpr_;
    DROPS::MLMatrixCL  P_;
    DROPS::Uint iter_prA_;
    DROPS::Uint iter_prM_;
    mutable std::vector<DROPS::VectorCL> ones_;

    void MaybeInitOnes() const;

  public:
    ISMGPreCL(DROPS::MLMatrixCL& A_pr, DROPS::MLMatrixCL& M_pr,
                    double kA, double kM, DROPS::Uint iter_prA=1,
                    DROPS::Uint iter_prM = 1)
        : SchurPreBaseCL( kA, kM), sm( 1), lvl( -1), omega( 1.0), smoother( omega), solver( directpc, 200, 1e-12),
          Apr_( A_pr), Mpr_( M_pr), iter_prA_( iter_prA), iter_prM_( iter_prM), ones_(0)
    {}

    template <typename Mat, typename Vec>
    void
    Apply(const Mat& /*A*/, Vec& p, const Vec& c) const;

    MLMatrixCL* GetProlongation() { return &P_; }
};


//**************************************************************************
// Preconditioner for the instationary Stokes-equations.
// It uses BB^T instead of a Laplacian on the pressure space and is
// therefore suited for P1X-elements;
// cf. "Uniform Preconditioners for a Parameter Dependent Saddle Point
// Problem with Application to Generalized Stokes Interface Equations",
// Olshanskii, Peters, Reusken, 2005
//**************************************************************************
class ISBBTPreCL : public SchurPreBaseCL
{
  private:
    const MatrixCL*  B_;
    mutable MatrixCL*  Bs_;                                     ///< scaled Matrix B
    mutable size_t Bversion_;
    const MatrixCL*  M_, *Mvel_;

    double     tolA_, tolM_;                                    ///< tolerances of the solvers
    mutable VectorCL Dprsqrtinv_;                               ///< diag(M)^{-1/2}
#ifndef _PAR
    typedef NEGSPcCL SPcT_;
    SPcT_            spc_;
    JACPcCL          jacpc_;
    mutable PCGNESolverCL<SPcT_> solver_;
    mutable PCGSolverCL<JACPcCL> solver2_;
#else
    mutable CompositeMatrixCL BBT_;
    typedef ParJacNEG0CL    PCSolver1T;                         ///< type of the preconditioner for solver 1
    typedef ParJac0CL       PCSolver2T;                         ///< type of the preconditioner for solver 2
    PCSolver1T PCsolver1_;
    PCSolver2T PCsolver2_;
    mutable ParPCGSolverCL<PCSolver1T> solver_;                 ///< solver for BB^T
    mutable ParPCGSolverCL<PCSolver2T> solver2_;                ///< solver for M

    const IdxDescCL* vel_idx_;                                  ///< Accessing ExchangeCL for velocity
#endif
    const IdxDescCL* pr_idx_;                                   ///< Accessing ExchangeCL for pressure; also used to determine, how to represent the kernel of BB^T in case of pure Dirichlet-BCs.
    double regularize_;                                         ///< If regularize_==0. no regularization is performed. Otherwise, a column is attached to Bs.
    void Update () const;                                       ///< Updating the diagonal matrices D and Dprsqrtinv

  public:
#ifndef _PAR
    ISBBTPreCL (const MatrixCL* B, const MatrixCL* M_pr, const MatrixCL* Mvel,
        const IdxDescCL& pr_idx,
        double kA= 0., double kM= 1., double tolA= 1e-2, double tolM= 1e-2, double regularize= 0.)
        : SchurPreBaseCL( kA, kM), B_( B), Bs_( 0), Bversion_( 0),
          M_( M_pr), Mvel_( Mvel), tolA_(tolA), tolM_(tolM),
          solver_( spc_, 500, tolA_, /*relative*/ true),
          solver2_( jacpc_, 500, tolM_, /*relative*/ true),
          pr_idx_( &pr_idx), regularize_( regularize) {}

    ISBBTPreCL (const ISBBTPreCL& pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), B_( pc.B_), Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Bversion_( pc.Bversion_),
          M_( pc.M_), Mvel_( pc.Mvel_),
          tolA_(pc.tolA_), tolM_(pc.tolM_),
          Dprsqrtinv_( pc.Dprsqrtinv_),
          spc_( pc.spc_),
          solver_( spc_, 500, tolA_, /*relative*/ true),
          solver2_( jacpc_, 500, tolM_, /*relative*/ true),
          pr_idx_( pc.pr_idx_), regularize_( pc.regularize_) {}
#else
    ISBBTPreCL (const MatrixCL* B, const MatrixCL* M_pr, const MatrixCL* Mvel,
        const IdxDescCL& pr_idx, const IdxDescCL& vel_idx,
        double kA= 0., double kM= 1., double tolA= 1e-2, double tolM= 1e-2, double regularize= 0.)
        : SchurPreBaseCL( kA, kM), B_( B), Bs_( 0), Bversion_( 0),
          M_( M_pr), Mvel_( Mvel), tolA_(tolA), tolM_(tolM),
          BBT_( 0, TRANSP_MUL, 0, MUL, vel_idx, pr_idx),
          PCsolver1_( pr_idx), PCsolver2_(pr_idx),
          solver_( 800, tolA_, pr_idx, PCsolver1_, /*relative*/ true, /*accure*/ true),
          solver2_( 500, tolM_, pr_idx, PCsolver2_, /*relative*/ true),
          vel_idx_( &vel_idx), pr_idx_( &pr_idx), regularize_( regularize) {}
    ISBBTPreCL (const ISBBTPreCL& pc)
        : SchurPreBaseCL( pc.kA_, pc.kM_), B_( pc.B_), Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Bversion_( pc.Bversion_),
          M_( pc.M_), Mvel_( pc.Mvel_),
          tolA_(pc.tolA_), tolM_(pc.tolM_),
          Dprsqrtinv_( pc.Dprsqrtinv_),
          BBT_( Bs_, TRANSP_MUL, Bs_, MUL, *pc.vel_idx_, *pc.pr_idx_),
          PCsolver1_( *pc.pr_idx_), PCsolver2_( *pc.pr_idx_),
          solver_( 800, tolA_, *pc.pr_idx_, PCsolver1_, /*relative*/ true, /*accure*/ true),
          solver2_( 500, tolM_, *pc.pr_idx_, PCsolver2_, /*relative*/ true),
          vel_idx_( pc.vel_idx_), pr_idx_( pc.pr_idx_), regularize_( pc.regularize_){}

    /// \name Parallel preconditioner setup ...
    //@{
    bool NeedDiag() const { return false; }
    void SetDiag(const VectorCL&) {}        // just for consistency
    template<typename Mat>
    void SetDiag(const Mat&) {}             // just for consistency
    bool RetAcc()   const { return true; }
    //@}
#endif

    ISBBTPreCL& operator= (const ISBBTPreCL&) {
        throw DROPSErrCL( "ISBBTPreCL::operator= is not permitted.\n");
    }

    ~ISBBTPreCL () { delete Bs_; }

    template <typename Mat, typename Vec>
    void Apply(const Mat&, Vec& p, const Vec& c) const;

    void SetMatrices (const MatrixCL* B, const MatrixCL* Mvel, const MatrixCL* M, const IdxDescCL* pr_idx) {
        B_= B;
        Mvel_= Mvel;
        M_= M;
        pr_idx_= pr_idx;
        Bversion_ = 0;
    }
};

template <typename Mat, typename Vec>
void ISBBTPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
    if (B_->Version() != Bversion_)
        Update();

    p= 0.0;
    if (kA_ != 0.0) {
#ifndef _PAR
        solver_.Solve( *Bs_, p, VectorCL( Dprsqrtinv_*c));
#else
        solver_.Solve( BBT_, p, VectorCL( Dprsqrtinv_*c));
#endif
//        IF_MASTER
//            std::cout << "ISBBTPreCL p: iterations: " << solver_.GetIter()
//                       << "\tresidual: " <<  solver_.GetResid();
        if (solver_.GetIter() == solver_.GetMaxIter()){
          IF_MASTER
            std::cout << "ISBBTPreCL::Apply: BBT-solve: " << solver_.GetIter()
                    << '\t' << solver_.GetResid() << '\n';
        }
        p= kA_*(Dprsqrtinv_*p);
    }
    if (kM_ != 0.0) {
        Vec p2_( c.size());
        solver2_.Solve( *M_, p2_, c);
//        IF_MASTER
//            std::cout << "\tISBBTPreCL p2: iterations: " << solver2_.GetIter()
//                       << "\tresidual: " <<  solver2_.GetResid()
//                       << '\n';
        if (solver2_.GetIter() == solver2_.GetMaxIter()){
          IF_MASTER
            std::cout << "ISBBTPreCL::Apply: M-solve: " << solver2_.GetIter()
                    << '\t' << solver2_.GetResid() << '\n';
        }

        p+= kM_*p2_;
    }
}

//**************************************************************************
// Preconditioner for the instationary (Navier-) Stokes-equations.
// It is a scaled version of the Min-Commutator-PC of Elman and can be used
// with P1X-elements.
//**************************************************************************
#ifndef _PAR
class MinCommPreCL
{
  private:
    const MatrixCL* A_, *B_, *Mvel_, *M_;
    mutable MatrixCL* Bs_;
    mutable size_t Aversion_, Bversion_, Mvelversion_, Mversion_;
    mutable VectorCL Dprsqrtinv_, Dvelsqrtinv_;
    double  tol_;

    typedef NEGSPcCL SPcT_;
    SPcT_ spc_;
    mutable PCGNESolverCL<SPcT_> solver_;
    const IdxDescCL* pr_idx_;                                   ///< Used to determine, how to represent the kernel of BB^T in case of pure Dirichlet-BCs.
    double regularize_;

    void Update () const;

  public:
    MinCommPreCL (const MatrixCL* A, MatrixCL* B, MatrixCL* Mvel, MatrixCL* M_pr, const IdxDescCL& pr_idx,
                  double tol=1e-2, double regularize= 0.0)
        : A_( A), B_( B), Mvel_( Mvel), M_( M_pr), Bs_( 0),
          Aversion_( 0), Bversion_( 0), Mvelversion_( 0), Mversion_( 0),
          tol_(tol),
          spc_( /*symmetric GS*/ true), solver_( spc_, 200, tol_, /*relative*/ true),
          pr_idx_( &pr_idx), regularize_( regularize) {}

    MinCommPreCL (const MinCommPreCL & pc)
        : A_( pc.A_), B_( pc.B_), Mvel_( pc.Mvel_), M_( pc.M_),
          Bs_( pc.Bs_ == 0 ? 0 : new MatrixCL( *pc.Bs_)),
          Aversion_( pc.Aversion_), Bversion_( pc.Bversion_), Mvelversion_( pc.Mvelversion_),
          Mversion_( pc.Mversion_),
          Dprsqrtinv_( pc.Dprsqrtinv_), Dvelsqrtinv_( pc.Dvelsqrtinv_), tol_(pc.tol_),
          spc_( pc.spc_), solver_( spc_, 200, tol_, /*relative*/ true),
          pr_idx_( pc.pr_idx_), regularize_( pc.regularize_) {}

    MinCommPreCL& operator= (const MinCommPreCL&) {
        throw DROPSErrCL( "MinCommPreCL::operator= is not permitted.\n");
    }

    ~MinCommPreCL () { delete Bs_; }

    template <typename Mat, typename Vec>
    void Apply (const Mat&, Vec& x, const Vec& b) const;

    void SetMatrixA  (const MatrixCL* A) { A_= A; Aversion_= 0; }
    void SetMatrices (const MatrixCL* A, const MatrixCL* B, const MatrixCL* Mvel, const MatrixCL* M,
                      const IdxDescCL* pr_idx) {
        A_= A;
        B_= B;
        Mvel_= Mvel;
        M_= M;
        pr_idx_= pr_idx;
        Aversion_ = Bversion_ = Mvelversion_ = Mversion_ = 0;
    }
};

template <typename Mat, typename Vec>
  void
  MinCommPreCL::Apply (const Mat&, Vec& x, const Vec& b) const
{
    if ((A_->Version() != Aversion_) || (Mvel_->Version() != Mvelversion_) || (B_->Version() != Bversion_))
        Update();

    VectorCL y( b.size());
    solver_.Solve( *Bs_, y, VectorCL( Dprsqrtinv_*b));
    if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "MinCommPreCL::Apply: 1st BBT-solve: " << solver_.GetIter()
                  << '\t' << solver_.GetResid() << '\n';
    y*= Dprsqrtinv_;
    VectorCL z( Dprsqrtinv_*((*B_)*VectorCL( Dvelsqrtinv_*Dvelsqrtinv_*
        ( (*A_)*VectorCL( Dvelsqrtinv_*Dvelsqrtinv_*transp_mul( *B_, y)) ))));
    VectorCL t( b.size());
    solver_.Solve( *Bs_, t, z);
    if (solver_.GetIter() == solver_.GetMaxIter())
        std::cout << "MinCommPreCL::Apply: 2nd BBT-solve: " << solver_.GetIter()
                  << '\t' << solver_.GetResid() << '\n';
    x= Dprsqrtinv_*t;
}
#endif

//=================================
//     template definitions
//=================================

template <class StokesT, class SolverT>
void InstatStokesThetaSchemeCL<StokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * v;
    _rhs*= (_theta-1.);
    _rhs+= (1./_dt)*(_Stokes.M.Data*v + _cplM->Data - _old_cplM->Data)
         +  _theta*_b->Data + (1.-_theta)*_old_b->Data;

    _solver.Solve( _mat, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
}


template <typename Mat, typename Vec>
void ISPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
//    double new_res;
//    double old_res= norm( c);
    ssor_.Apply( A_, p, c);
//    std::cout << " residual: " <<  (new_res= norm( A_*p - c)) << '\t';
//    std::cout << " reduction: " << new_res/old_res << '\t';
    p*= kA_;
//    double mnew_res;
//    double mold_res= norm( c);
    Vec p2_( c.size());
    ssor_.Apply( M_, p2_, c);
//    std::cout << " residual: " <<  (mnew_res= norm( M_*p2_ - c)) << '\t';
//    std::cout << " reduction: " << mnew_res/mold_res << '\n';
    p+= kM_*p2_;
}

#ifndef _PAR
template <class SolverT>
template <typename Mat, typename Vec>
void ISNonlinearPreCL<SolverT>::Apply(const Mat&, Vec& p, const Vec& c) const
{
    p= 0.0;
//    std::cout << "ISNonlinearPreCL";
    if (kA_ != 0.0) {
        solver_.Solve( A_, p, c);
//        std::cout << "\tp: iterations: " << solver_.GetIter()
//                  << "\tresidual: " <<  solver_.GetResid();
        p*= kA_;
    }
    if ( kM_ != 0.0) {
        Vec p2_( c.size());
        solver_.Solve( M_, p2_, c);
//        std::cout << "\t p2: iterations: " << solver_.GetIter()
//                  << "\tresidual: " <<  solver_.GetResid()
//                  << '\n';
        p+= kM_*p2_;
    }
}
#else
template <typename ASolverT, typename MSolverT>
template <typename Mat, typename Vec>
void ISNonlinearPreCL<ASolverT, MSolverT>::Apply(const Mat&, Vec& p, const Vec& c) const
{
    p= 0.0;
    if (kA_ != 0.0) {
        Asolver_.Solve( A_, p, c);
//     IF_MASTER
//       std::cout << "ISNonlinearPreCL p: iterations: " << solver_.GetIter()
//                 << "\tresidual: " <<  solver_.GetResid()<<std::endl;
        p*= kA_;
    }
    if ( kM_ != 0.0) {
        Vec p2_( c.size());
        Msolver_.Solve( M_, p2_, c);
//         IF_MASTER
//           std::cout << "\t p2: iterations: " << solver_.GetIter()
//                     << "\tresidual: " <<  solver_.GetResid()
//                     << '\n';
        p+= kM_*p2_;
    }
}
#endif

/* with orthogonalization
template <typename Mat, typename Vec>
void ISNonlinearPreCL::Apply(const Mat&, Vec& p, const Vec& c) const
{
    VectorCL e( 1., p.size());
    const double ee= p.size(); // = dot( e, e);
    VectorCL c2( (dot(c, e)/ee)*e);
    VectorCL c3( c - c2);
std::cout << "norm( c): " << norm( c) << "\tnorm( e): " << norm( e)
          << "\tdot(c, e)/norm( c)/norm( e): " << dot(c, e)/norm( c)/ee << '\n';
    p= 0.0;
    solver_.Solve( A_, p, c3);
    p+= c2;
    std::cout << "ISNonlinearPreCL p: iterations: " << solver_.GetIter()
              << "\tresidual: " <<  solver_.GetResid();
    p*= kA_;
    if ( kM_ != 0.0) {
        Vec p2_( c.size());
        solver_.Solve( M_, p2_, c);
std::cout << "norm( p2): " << norm( p2_);
        std::cout << "\t p2: iterations: " << solver_.GetIter()
                  << "\tresidual: " <<  solver_.GetResid()
                  << '\n';
        p+= kM_*p2_;
    }
}
*/

template<class SmootherCL, class DirectSolverCL>
void
MGMPr(const std::vector<VectorCL>::const_iterator& ones,
      const MLMatrixCL::const_iterator& begin, const MLMatrixCL::const_iterator& fine,
      MLMatrixCL::const_iterator P, VectorCL& x, const VectorCL& b,
      const SmootherCL& Smoother, const Uint smoothSteps,
      DirectSolverCL& Solver, const int numLevel, const int numUnknDirect)
// Multigrid method, V-cycle. If numLevel==0 or #Unknowns <= numUnknDirect,
// the direct solver Solver is used.
// If one of the parameters is -1, it will be neglected.
// If MLMatrixCL.begin() has been reached, the direct solver is used too.
// Concerning the stabilization see Hackbusch Multigrid-Methods and Applications;
// Basically we project on the orthogonal complement of the kernel of A before
// the coarse-grid correction.
{
    MLMatrixCL::const_iterator coarse = fine;
    MLMatrixCL::const_iterator coarseP= P;
    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) )
       || fine==begin)
    { // use direct solver
        Solver.Solve( *fine, x, b);
        x-= dot( *ones, x);
        return;
    }
    --coarse;
    --coarseP;
    // presmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine, x, b);
    // restriction of defect
    VectorCL d( transp_mul( *P, VectorCL( b - (*fine)*x)));
    d-= dot( *(ones-1), d);
    VectorCL e( d.size());
    // calculate coarse grid correction
    MGMPr( ones-1, begin, coarse, coarseP, e, d, Smoother, smoothSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    // add coarse grid correction
    x+= (*P) * e;
    // postsmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine, x, b);
    // This projection could probably be avoided, but it is cheap and with it,
    // we are on the safe side.
    x-= dot( *ones, x);
}


template <typename Mat, typename Vec>
void
ISMGPreCL::Apply(const Mat& /*A*/, Vec& p, const Vec& c) const
{
    MaybeInitOnes();
    p= 0.0;
    const Vec c2_( c - dot( ones_.back(), c));
    MLMatrixCL::const_iterator finestP = --P_.end();
//    double new_res= (Apr_.back().A.Data*p - c).norm();
//    double old_res;
//    std::cout << "Pressure: iterations: " << iter_prA_ <<'\t';
    for (DROPS::Uint i=0; i<iter_prA_; ++i) {
        DROPS::MGMPr( ones_.end()-1, Apr_.begin(), --Apr_.end(), finestP, p, c2_, smoother, sm, solver, lvl, -1);
//        old_res= new_res;
//        std::cout << " residual: " <<  (new_res= (Apr_.back().A.Data*p - c).norm()) << '\t';
//        std::cout << " reduction: " << new_res/old_res << '\n';
    }
    p*= kA_;
//    std::cout << " residual: " <<  (Apr_.back().A.Data*p - c).norm() << '\t';

    Vec p2( p.size());
    for (DROPS::Uint i=0; i<iter_prM_; ++i)
        DROPS::MGM( Mpr_.begin(), --Mpr_.end(), finestP, p2, c, smoother, sm, solver, lvl, -1);
//    std::cout << "Mass: iterations: " << iter_prM_ << '\t'
//              << " residual: " <<  (Mpr_.back().A.Data*p2 - c).norm() << '\n';

    p+= kM_*p2;
}

}    // end of namespace DROPS

#endif
