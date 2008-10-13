//**************************************************************************
// File:    poisson.h                                                      *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

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


template <class Coeff>
class PoissonP1CL : public ProblemCL<Coeff, PoissonBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, PoissonBndDataCL> _base;
    typedef typename _base::BndDataCL          BndDataCL;
    typedef typename _base::CoeffCL            CoeffCL;
    using                                      _base::_MG;
    using                                      _base::_Coeff;
    using                                      _base::_BndData;
    using                                      _base::GetBndData;
    using                                      _base::GetMG;

    typedef P1EvalCL<double, const BndDataCL, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataCL, const VecDescCL> const_DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    IdxDescCL idx;
    VecDescCL x;
    VecDescCL b;
    MatDescCL A;

    PoissonP1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mgb, coeff, bdata), idx( P1_FE) {}
    // numbering of unknowns
    void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, _MG, _BndData, match); }
    void DeleteNumbering( IdxDescCL* idx)
        { idx->DeleteNumbering( _MG); }
    // set up matrices and rhs
    void SetupSystem         (MatDescCL&, VecDescCL&) const;
    void SetupStiffnessMatrix(MatDescCL&) const;
    void SetupProlongation   (MatDescCL& P, IdxDescCL* cIdx, IdxDescCL* fIdx) const;
    // check computed solution etc.
    double CheckSolution(const VecDescCL&, instat_scalar_fun_ptr) const;
    double CheckSolution(instat_scalar_fun_ptr Lsg) const { return CheckSolution(x, Lsg); }
    void GetDiscError (const MatDescCL&, instat_scalar_fun_ptr) const;
    void GetDiscError (scalar_fun_ptr Lsg) const { GetDiscError(A, Lsg); }

    bool          EstimateError         (const VecDescCL&, const double, double&, est_fun);
    static double ResidualErrEstimator  (const TetraCL&, const VecDescCL&, const BndDataCL&);
    static double ResidualErrEstimatorL2(const TetraCL&, const VecDescCL&, const BndDataCL&);

    DiscSolCL GetSolution()
        { return DiscSolCL(&x, &GetBndData(), &GetMG()); }
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL(&x, &GetBndData(), &GetMG()); }
};

template <class Coeff>
class PoissonP2CL : public ProblemCL<Coeff, PoissonBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, PoissonBndDataCL> _base;
    typedef typename _base::BndDataCL               BndDataCL;
    typedef typename _base::CoeffCL                 CoeffCL;
    using                                           _base::_MG;
    using                                           _base::_Coeff;
    using                                           _base::_BndData;
    using                                           _base::GetBndData;
    using                                           _base::GetMG;

    typedef P2EvalCL<double, const BndDataCL, VecDescCL>       DiscSolCL;
    typedef P2EvalCL<double, const BndDataCL, const VecDescCL> const_DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    // new fields for the matrix A, the rhs b and the solution x
    IdxDescCL idx;
    VecDescCL x;
    VecDescCL b;
    MatDescCL A;

    //create an element of the class
    PoissonP2CL(const MGBuilderCL& mgb, const CoeffCL& coeff,
                const BndDataCL& bdata) : _base(mgb, coeff, bdata), idx(P2_FE) {}

    // numbering of unknowns
    void CreateNumbering( Uint level, IdxDescCL* idx, match_fun match= 0)
        { idx->CreateNumbering( level, _MG, _BndData, match); }
    void DeleteNumbering( IdxDescCL* idx)
        { idx->DeleteNumbering( _MG); }

    // set up matrices and rhs
    void SetupSystem         (MatDescCL&, VecDescCL&) const;
    void SetupStiffnessMatrix(MatDescCL&) const;

    // check computed solution, etc.
    double CheckSolution(const VecDescCL&, instat_scalar_fun_ptr) const;
    double CheckSolution(instat_scalar_fun_ptr Lsg) const { return CheckSolution(x, Lsg); }
};


double SimpleGradEstimator (const TetraCL& t, const VecDescCL& lsg, const PoissonBndDataCL&);

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
        {}
    // default assignment-op, copy-ctor, dtor

    template <class _BndData, class _VD>
      void Init(const P1EvalCL<double, _BndData, _VD>&);

    double GetRelRed() { return _RelReduction; }
    void SetRelRed(double newred) { _RelReduction= newred; }
    bool DoesMark() { return _DoMark; }
    void SwitchMark() { _DoMark= _DoMark ? false : true; }
    template <class _BndData, class _VD>
      bool Estimate(const P1EvalCL<double, _BndData, const _VD>&);
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
        {}
    // default assignment-op, copy-ctor, dtor

    template <class _BndData, class _VD>
      void Init(const P1EvalCL<double, _BndData, _VD>&);

    double GetRelRed() { return _RelReduction; }
    void   SetRelRed(double newred) { _RelReduction= newred; }
    double GetThreshold() { return _Threshold; }
    void   SetThreshold(double newThreshold) { _Threshold= newThreshold; }
    bool   DoesMark() { return _DoMark; }
    void   SwitchMark() { _DoMark= _DoMark ? false : true; }
    template <class _BndData, class _VD>
      bool Estimate(const P1EvalCL<double, _BndData, const _VD>&);
};


} // end of namespace DROPS

#include "poisson/poisson.tpp"

#endif
