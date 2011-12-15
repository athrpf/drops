/// \file poisson.h
/// \brief classes that constitute the poisson-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_POISSON_H
#define DROPS_POISSON_H

#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/bndData.h"
#include <deque>
#include <iostream>
#include <numeric>


namespace DROPS
{

typedef BndSegDataCL<> PoissonBndSegDataCL;
typedef BndDataCL<> PoissonBndDataCL;

double Py_product(MultiGridCL& mg, MLIdxDescCL& Idx, MLMatDescCL& A, MLMatDescCL& M,
		   instat_scalar_fun_ptr f1, instat_scalar_fun_ptr f2, double t, bool H1);

class StripTimeCL
// converts time dependent function to one, that is time independent.
// i.e.:  scalar_instat_fun_ptr  -->  scalar_fun_ptr
// useful where one wants to use all the stuff written for the stationary problems
{
  private:
    static instat_scalar_fun_ptr _func;
    static double                _t;
  public:
    StripTimeCL( instat_scalar_fun_ptr func, double t)
      { _func= func; _t= t; }
    void SetTime( double t) { _t= t; }

    static double GetFunc( const Point3DCL& x)
      { return _func( x, _t); }
};

template <class Coeff>
class PoissonP1CL : public ProblemCL<Coeff, PoissonBndDataCL>
{
  private:
    bool adjoint_;

  public:
    typedef ProblemCL<Coeff, PoissonBndDataCL> base_;
    typedef typename base_::BndDataCL          BndDataCL;
    typedef typename base_::CoeffCL            CoeffCL;
    using                                      base_::MG_;
    using                                      base_::Coeff_;
    using                                      base_::BndData_;
    using                                      base_::GetBndData;
    using                                      base_::GetMG;

    typedef P1EvalCL<double, const BndDataCL, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataCL, const VecDescCL> const_DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    MLIdxDescCL idx;
    VecDescCL   x;
    VecDescCL   b;
    VecDescCL   vU;
    MLMatDescCL A;
    MLMatDescCL M;
    MLMatDescCL U;

    PoissonP1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata, bool adj=false)
        : base_( mgb, coeff, bdata), adjoint_( adj), idx( P1_FE) {}

    PoissonP1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata, bool adj=false)
        : base_( mg, coeff, bdata), adjoint_( adj), idx( P1_FE) {}
    // numbering of unknowns
    void CreateNumbering( Uint level, MLIdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, MG_, BndData_, match); }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    void SetNumLvl( size_t n);

    // set up matrices and rhs
    void SetupSystem         (MLMatDescCL&, VecDescCL&, bool SUPG=false) const;
    ///  \brief set up matrices (M is time independent)
    void SetupInstatSystem( MLMatDescCL& A, MLMatDescCL& M, double t, bool SUPG=false) const;
    /// \brief set up matrix and couplings with bnd unknowns for convection term
    void SetupConvection( MLMatDescCL& U, VecDescCL& vU, double t) const;

    /// \brief Setup time dependent parts
    ///
    /// couplings with bnd unknowns, coefficient f(t)
    /// If the function is called with the same vector for some arguments,
    /// the vector will contain the sum of the results after the call
    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf, bool SUPG=false) const;
    /// \brief Setup special source term including the gradient of a given P1 function
    void SetupGradSrc( VecDescCL& src, instat_scalar_fun_ptr T, instat_scalar_fun_ptr dalpha, double t= 0.) const;

    void SetupL2ProjGrad(VecDescCL& r, instat_scalar_fun_ptr T, instat_scalar_fun_ptr Psi, instat_scalar_fun_ptr flux, double t= 0.) const;
    
    /// \brief Set initial value
    void Init( VecDescCL&, instat_scalar_fun_ptr, double t0= 0.) const;

    /// \brief check computed solution etc.
    double CheckSolution( const VecDescCL&, instat_scalar_fun_ptr) const;
    double CheckSolution( instat_scalar_fun_ptr Lsg) const { return CheckSolution(x, Lsg); }
    double CheckSolution( const VecDescCL&, instat_scalar_fun_ptr, double) const;
    double CheckSolution( instat_scalar_fun_ptr Lsg, double& t) const { return CheckSolution(x, Lsg, t); }

    void GetDiscError   ( const MLMatDescCL&, instat_scalar_fun_ptr) const;
    void GetDiscError   ( scalar_fun_ptr Lsg) const { GetDiscError(A, Lsg); }

    bool          EstimateError         ( const VecDescCL&, const double, double&, est_fun);
    static double ResidualErrEstimator  ( const TetraCL&, const VecDescCL&, const BndDataCL&);
    static double ResidualErrEstimatorL2( const TetraCL&, const VecDescCL&, const BndDataCL&);

    DiscSolCL GetSolution()
        { return DiscSolCL(&x, &GetBndData(), &GetMG()); }
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL(&x, &GetBndData(), &GetMG()); }
};

template <class Coeff>
class PoissonP2CL : public ProblemCL<Coeff, PoissonBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, PoissonBndDataCL> base_;
    typedef typename base_::BndDataCL               BndDataCL;
    typedef typename base_::CoeffCL                 CoeffCL;
    using                                           base_::MG_;
    using                                           base_::Coeff_;
    using                                           base_::BndData_;
    using                                           base_::GetBndData;
    using                                           base_::GetMG;

    typedef P2EvalCL<double, const BndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const BndDataCL, const VecDescCL> const_DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    // new fields for the matrix A, the rhs b and the solution x
    MLIdxDescCL idx;
    VecDescCL   x;
    VecDescCL   b;
    MLMatDescCL A;

    MLMatDescCL U;  //Convection matrix
    MLMatDescCL M;  //Mass matrix
    VecDescCL   vU; //Coupling with convection matrix

    //create an element of the class
    PoissonP2CL(const MGBuilderCL& mgb, const CoeffCL& coeff,
                const BndDataCL& bdata) : base_(mgb, coeff, bdata), idx(P2_FE) {}
    PoissonP2CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : base_( mg, coeff, bdata), idx( P2_FE) {}
    // numbering of unknowns
    void CreateNumbering( Uint level, MLIdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, MG_, BndData_, match); }
    void DeleteNumbering( MLIdxDescCL* idx)
        { idx->DeleteNumbering( MG_); }
    void SetNumLvl( size_t n);

    // set up matrices and rhs
    void SetupSystem         ( MLMatDescCL&, VecDescCL&, bool SUPG=false) const;

    ///  \brief set up matrices for instatProblem
    void SetupInstatSystem( MLMatDescCL& A, MLMatDescCL& M, double t, bool SUPG=false) const;

    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf, bool SUPG=false) const;

    //Set up convection
    void SetupConvection( MLMatDescCL&, VecDescCL&, double) const;

    //Set up initial value
    void Init( VecDescCL&, instat_scalar_fun_ptr, double t0= 0.) const;

    // check computed solution, etc.
    double CheckSolution( const VecDescCL&, instat_scalar_fun_ptr) const;
    double CheckSolution( instat_scalar_fun_ptr Lsg) const { return CheckSolution(x, Lsg); }

    double CheckSolution( const VecDescCL&, instat_scalar_fun_ptr, double) const;
    double CheckSolution( instat_scalar_fun_ptr Lsg, double& t) const { return CheckSolution(x, Lsg, t); }

    DiscSolCL GetSolution()
        { return DiscSolCL(&x, &GetBndData(), &GetMG()); }
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL(&x, &GetBndData(), &GetMG()); }
};


double SimpleGradEstimator ( const TetraCL& t, const VecDescCL& lsg, const PoissonBndDataCL&);

//==============================================================
//                  Marker classes
//==============================================================

template <class _TetraEst, class _ProblemCL>
class PoissonErrEstCL
{
  private:
    double    _InitGlobErr;
    double    _RelReduction;
    double    _ConvExponent;
    double    _Meas;
    double    _ActGlobErr;
    _TetraEst _Estimator;
    _ProblemCL& _Problem;
    Uint      _NumLastMarkedForRef;
    Uint      _NumLastMarkedForDel;
    bool      _DoMark;
    std::ostream*  _outp;

  public:
      PoissonErrEstCL(double RelReduction, double ConvExp, double Meas, bool DoMark, _TetraEst est,
                      _ProblemCL& problem, std::ostream* osp= &std::cout)
        : _InitGlobErr(0), _RelReduction(RelReduction), _ConvExponent(ConvExp), _Meas(Meas), _ActGlobErr(-1), _Estimator(est),
          _Problem(problem), _NumLastMarkedForRef(0), _NumLastMarkedForDel(0), _DoMark(DoMark), _outp(osp)
        {
#ifdef _PAR
            throw DROPSErrCL("This class has not been parallelized yes, sorry");
#endif
        }
    // default assignment-op, copy-ctor, dtor

    template <class BndData_, class _VD>
      void Init(const P1EvalCL<double, BndData_, _VD>&);

    double GetRelRed() { return _RelReduction; }
    void SetRelRed(double newred) { _RelReduction= newred; }
    bool DoesMark() { return _DoMark; }
    void SwitchMark() { _DoMark= _DoMark ? false : true; }
    template <class BndData_, class _VD>
      bool Estimate(const P1EvalCL<double, BndData_, const _VD>&);
};


template <class _TetraEst, class _ProblemCL>
class DoerflerMarkCL
{
  private:
    double        _InitGlobErr;
    double        _RelReduction;
    double        _min_tetra_ratio;
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
  // the tetras are sorted: T_1 with biggest error, last T_n with smallest
  // a tetra T_i is marked for refinement, iff (a) or (b) holds:
  // (a) it is among "min_ratio" % of the tetras with biggest errors,
  // (b) the sum err_1+..+err_i accounts for less than "Threshold" % of the global error

      DoerflerMarkCL(double RelReduction, double min_ratio, double Threshold, double Meas, bool DoMark, _TetraEst est,
                     _ProblemCL& problem, std::ostream* osp= &std::cout)
        : _InitGlobErr(0), _RelReduction(RelReduction), _min_tetra_ratio(min_ratio), _Threshold(Threshold), _Meas(Meas), _ActGlobErr(-1), _Estimator(est),
          _Problem(problem), _NumLastMarkedForRef(0), _NumLastMarkedForDel(0), _DoMark(DoMark), _outp(osp)
        {
#ifdef _PAR
            throw DROPSErrCL("This class has not been parallelized yes, sorry");
#endif
        }
    // default assignment-op, copy-ctor, dtor

    template <class BndData_, class _VD>
      void Init(const P1EvalCL<double, BndData_, _VD>&);

    double GetRelRed() { return _RelReduction; }
    void   SetRelRed(double newred) { _RelReduction= newred; }
    double GetThreshold() { return _Threshold; }
    void   SetThreshold(double newThreshold) { _Threshold= newThreshold; }
    bool   DoesMark() { return _DoMark; }
    void   SwitchMark() { _DoMark= _DoMark ? false : true; }
    template <class BndData_, class _VD>
      bool Estimate(const P1EvalCL<double, BndData_, const _VD>&);
};


} // end of namespace DROPS

#include "poisson/poisson.tpp"

#endif
