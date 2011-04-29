/// \file stokes.h
/// \brief classes that constitute the stokes-problem
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_STOKES_H
#define DROPS_STOKES_H

#include <vector>
#include "misc/problem.h"
#include "num/solver.h"
#include "num/fe.h"
#include "num/discretize.h"

#ifdef _PAR
#  include "parallel/exchange.h"
#endif

namespace DROPS
{

class StokesBndDataCL
{
  public:
    typedef BndDataCL<Point3DCL> VelBndDataCL;
    typedef BndDataCL<double>    PrBndDataCL;

    StokesBndDataCL( Uint numbndseg, const BndCondT* bc_vel, const VelBndDataCL::bnd_val_fun* fun, const BndCondT* bc_pr= 0)
        : Pr( numbndseg, bc_pr), Vel( numbndseg, bc_vel, fun) {}
    StokesBndDataCL(Uint numbndseg, const bool* isneumann, const VelBndDataCL::bnd_val_fun* fun)
        : Pr( numbndseg), Vel(numbndseg, isneumann, fun) {} // deprecated
    StokesBndDataCL(const VelBndDataCL & aVel, const PrBndDataCL & aPr):Pr(aPr),Vel(aVel){}
    const PrBndDataCL  Pr;
    const VelBndDataCL Vel;
    typedef VelBndDataCL::bnd_val_fun bnd_val_fun;
};


typedef StokesBndDataCL::VelBndDataCL StokesVelBndDataCL;
typedef StokesBndDataCL::PrBndDataCL  StokesPrBndDataCL;

typedef SMatrixCL<3, 3> (*jacobi_fun_ptr)(const Point3DCL&);

typedef VecDescBaseCL<VectorCL> VelVecDescCL;


template <class Coeff>
class StokesP2P1CL : public ProblemCL<Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, StokesBndDataCL> base_;
    typedef typename base_::CoeffCL                 CoeffCL;
    typedef typename base_::BndDataCL               BndDataCL;
    using base_::MG_;
    using base_::Coeff_;
    using base_::BndData_;
    using base_::GetBndData;
    using base_::GetMG;

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>           DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>           const_DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double);

    MLIdxDescCL  vel_idx;  // for velocity unknowns
    MLIdxDescCL  pr_idx;   // for pressure unknowns
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MLMatDescCL  A,
                 B,
                 M,
                 prA,
                 prM;

    StokesP2P1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : base_( mgb, coeff, bdata), vel_idx( vecP2_FE), pr_idx( P1_FE){}
    StokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : base_( mg, coeff, bdata), vel_idx( vecP2_FE), pr_idx( P1_FE) {}

    /// \name Create and delete numbering of unknowns
    //@{
    /// Within parallel these functions also create the Exchange classes.
    void CreateNumberingVel( Uint level, MLIdxDescCL* idx, match_fun match= 0);
    void CreateNumberingPr ( Uint level, MLIdxDescCL* idx, match_fun match= 0);
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    void SetNumVelLvl( size_t n);
    void SetNumPrLvl ( size_t n);
    //@}

    /// \brief Set up matrices and complete rhs
    void SetupSystem(MLMatDescCL*, VelVecDescCL*, MLMatDescCL*, VelVecDescCL*, double = 0.0) const;
    /// \brief  Set up only A.
    void SetupStiffnessMatrix(MLMatDescCL*) const;
    /// Set up the stiffness matrix for the pressure, scaled by \f$\rho^{-1}\f$.
    void SetupPrStiff(MLMatDescCL* ) const;
    /// \brief  Set up mass-matrix for pressure-unknowns (P1)
    void SetupPrMass(MLMatDescCL*) const;
    /// \brief  Set up mass-matrix for velocity-unknowns (P2) -- needed for MG-Theta-scheme,
    /// Time-independent
    void SetupMassMatrix(MLMatDescCL* matI) const;

    /// \brief  Setup time independent part of system
    void SetupInstatSystem( MLMatDescCL* A, MLMatDescCL* B, MLMatDescCL* M) const;
    /// \brief  Setup time dependent parts: couplings with bnd unknowns, coefficient f(t)
    /** If the function is called with the same vector for some arguments (out of 1, 2, 4),
        the vector will contain the sum of the results after the call*/
    void SetupInstatRhs( VelVecDescCL* vA, VelVecDescCL* vB, VelVecDescCL* vI, double tA, VelVecDescCL* vf, double tf) const;
    /// \brief  Set initial value for velocities
    void InitVel( VelVecDescCL*, instat_vector_fun_ptr, double t0= 0.) const;

    /// \brief  Check system and computed solution
    void GetDiscError (instat_vector_fun_ptr LsgVel, instat_scalar_fun_ptr LsgPr, double t= 0.0) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, instat_vector_fun_ptr, instat_matrix_fun_ptr, instat_scalar_fun_ptr, bool) const;
    
    /// \brief  estimation a la Verfuerth
    static double ResidualErrEstimator(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double t=0.0);

    //@{
    /// Get solutions as FE-functions
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }

    const_DiscPrSolCL GetPrSolution( const VecDescCL& pr) const
        { return const_DiscPrSolCL( &pr, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution( const VelVecDescCL& vel) const
        { return const_DiscVelSolCL( &vel, &GetBndData().Vel, &GetMG()); }
    //@}
    bool UsesXFEM() { return false; } // just for consistency
};

#ifndef _PAR
template <class Coeff>
class StokesP1BubbleP1CL : public ProblemCL<Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, StokesBndDataCL> base_;
    typedef typename base_::CoeffCL           CoeffCL;
    typedef typename base_::BndDataCL         BndDataCL;
    using                                     base_::MG_;
    using                                     base_::Coeff_;
    using                                     base_::BndData_;
    using                                     base_::GetBndData;
    using                                     base_::GetMG;

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>                 DiscPrSolCL;
    typedef P1BubbleEvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>                 const_DiscPrSolCL;
    typedef P1BubbleEvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double);

    MLIdxDescCL  vel_idx;  // for velocity unknowns
    MLIdxDescCL  pr_idx;   // for pressure unknowns
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MLMatDescCL  A,        // Bad, bad comma.
                 B;

    StokesP1BubbleP1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : base_( mgb, coeff, bdata), vel_idx( vecP1Bubble_FE), pr_idx( P1_FE) {}

    // Create and delete numbering of unknowns
    void CreateNumberingVel( Uint level, MLIdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, MG_, BndData_.Vel, match); }
    void CreateNumberingPr ( Uint level, MLIdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, MG_, BndData_.Pr, match); }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    void SetNumVelLvl( size_t n);
    void SetNumPrLvl ( size_t n);

    // Set up matrices and rhs
    void SetupSystem(MLMatDescCL*, VelVecDescCL*, MLMatDescCL*, VelVecDescCL*) const;
    void SetupPrMass(MLMatDescCL*) const;

    // Check system and computed solution
    void GetDiscError (instat_vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, instat_vector_fun_ptr, scalar_fun_ptr) const;

    // estimation a la Verfuerth
    static double ResidualErrEstimator(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double= 0.0);

    //@{
    /// Get solutions as FE-functions
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }

    const_DiscPrSolCL GetPrSolution( const VecDescCL& pr) const
        { return const_DiscPrSolCL( &pr, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution( const VelVecDescCL& vel) const
        { return const_DiscVelSolCL( &vel, &GetBndData().Vel, &GetMG()); }
    //@}
};


// Works only for the stationary Stokes-equations (id est, always t==0.0)
// because Estimate calls the estimation-function always with time==0.0.
template <class _TetraEst, class _ProblemCL>
class StokesDoerflerMarkCL
{
  private:
    double        _InitGlobErr;
    double        _RelReduction;
    double        _Threshold;
    double        _Meas;
    double        _ActGlobErr;
    _TetraEst     _Estimator;
    _ProblemCL&   _Problem;
    Uint          _NumLastMarkedForRef;
    Uint          _NumLastMarkedForDel;
    bool          _DoMark;
    std::ostream* _outp;

  public:
    typedef typename _ProblemCL::BndDataCL BndDataCL;
    typedef typename BndDataCL::PrBndDataCL PrBndDataCL;
    typedef typename BndDataCL::VelBndDataCL VelBndDataCL;
    typedef typename _ProblemCL::const_DiscPrSolCL const_DiscPrSolCL;
    typedef typename _ProblemCL::const_DiscVelSolCL const_DiscVelSolCL;

  // the tetras are sorted: T_1 with biggest error, last T_n with smallest
  // a tetra T_i is marked for refinement, iff (a) or (b) holds:
  // (a) it is among "min_ratio" % of the tetras with biggest errors,
  // (b) the sum err_1+..+err_i accounts for less than "Threshold" % of the global error

      StokesDoerflerMarkCL(double RelReduction, double Threshold, double Meas, bool DoMark, _TetraEst est,
                           _ProblemCL& problem, std::ostream* osp= &std::cout)
        : _InitGlobErr(0), _RelReduction(RelReduction), _Threshold(Threshold), _Meas(Meas), _ActGlobErr(-1), _Estimator(est),
          _Problem(problem), _NumLastMarkedForRef(0), _NumLastMarkedForDel(0), _DoMark(DoMark), _outp(osp)
        {}
    // default assignment-op, copy-ctor, dtor

    void Init(const const_DiscPrSolCL&, const const_DiscVelSolCL&);

    double GetRelRed() { return _RelReduction; }
    void   SetRelRed(double newred) { _RelReduction= newred; }
    double GetThreshold() { return _Threshold; }
    void   SetThreshold(double newThreshold) { _Threshold= newThreshold; }
    bool   DoesMark() { return _DoMark; }
    void   SwitchMark() { _DoMark= _DoMark ? false : true; }
    bool Estimate(const const_DiscPrSolCL&, const const_DiscVelSolCL&);
};
#endif // end of ifndef _PAR

//======================================
//        inline functions
//======================================


} // end of namespace DROPS

#include "stokes/stokes.tpp"

#endif
