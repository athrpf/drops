//**************************************************************************
// File:    instatpoisson.h                                                *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// History: begin - November, 11 2002                                      *
//**************************************************************************

#ifndef DROPS_INSTATPOISSON_H
#define DROPS_INSTATPOISSON_H

#include <vector>
#include "misc/problem.h"
#include "num/solver.h"
#include "num/fe.h"
#ifdef _PAR
# include "parallel/exchange.h"
# include "parallel/parallel.h"
#endif


namespace DROPS
{

typedef BndSegDataCL<double> InstatPoissonBndSegDataCL;
typedef BndDataCL<double>    InstatPoissonBndDataCL;
typedef double  (*scalar_instat_fun_ptr)( const Point3DCL&, double);
typedef double  (*scalar_fun_ptr)( const Point3DCL&);

class StripTimeCL
// converts time dependent function to one, that is time independent.
// i.e.:  scalar_instat_fun_ptr  -->  scalar_fun_ptr
// useful where one wants to use all the stuff written for the stationary problems
{
  private:
    static scalar_instat_fun_ptr _func;
    static double                _t;
  public:
    StripTimeCL( scalar_instat_fun_ptr func, double t)
      { _func= func; _t= t; }
    void SetTime( double t) { _t= t; }

    static double GetFunc( const Point3DCL& x)
      { return _func( x, _t); }
};


template <class Coeff>
#ifndef _PAR
class InstatPoissonP1CL : public ProblemCL<Coeff, InstatPoissonBndDataCL>
#else
class InstatPoissonP1CL : public ProblemCL<Coeff, InstatPoissonBndDataCL, ExchangeCL>
#endif
{
  private:
    bool adjoint_;

  public:
#ifndef _PAR
    typedef ProblemCL<Coeff, InstatPoissonBndDataCL>       _base;
#else
    typedef ProblemCL<Coeff, InstatPoissonBndDataCL, ExchangeCL> _base;
#endif
    typedef typename _base::BndDataCL                       BndDataCL;
    typedef typename _base::CoeffCL                         CoeffCL;
    using _base::GetBndData;
    using _base::GetMG;
    using _base::_BndData;
    using _base::_MG;
    using _base::_Coeff;
#ifdef _PAR
    using _base::ex_;
    using _base::GetEx;
#endif

    typedef P1EvalCL<double, const BndDataCL, const VecDescCL> DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    double    t;        // time

    MLIdxDescCL idx;
    VecDescCL   x;
    VecDescCL   b;
    MLMatDescCL A;
    MLMatDescCL M;
    MLMatDescCL U;


    InstatPoissonP1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata, bool adj=false)
        : _base( mgb, coeff, bdata), adjoint_( adj), t( 0.), idx( P1_FE) {}

    // numbering of unknowns
    void CreateNumbering( Uint level, MLIdxDescCL* idx, match_fun match= 0);
    void DeleteNumbering( MLIdxDescCL* idx)
    { idx->DeleteNumbering( _MG); }
    void SetNumLvl( size_t n);

    ///  \brief set up matrices (M is time independent)
    void SetupInstatSystem( MLMatDescCL& A, MLMatDescCL& M, double tA) const;
    /// \brief set up matrix and couplings with bnd unknowns for convection term
    void SetupConvection( MLMatDescCL& U, VecDescCL& vU, double t) const;

    /// \brief Setup time dependent parts
    ///
    /// couplings with bnd unknowns, coefficient f(t)
    /// If the function is called with the same vector for some arguments,
    /// the vector will contain the sum of the results after the call
    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf) const;
    /// \brief Setup special source term including the gradient of a given P1 function
    void SetupGradSrc( VecDescCL& src, scalar_instat_fun_ptr T, scalar_instat_fun_ptr dalpha, double t= 0.) const;

    /// \brief Set initial value
    void Init( VecDescCL&, scalar_instat_fun_ptr, double t0= 0.) const;

    /// \brief check computed solution etc.
    void CheckSolution( const VecDescCL&, scalar_instat_fun_ptr, double) const;
    void CheckSolution( scalar_instat_fun_ptr Lsg, double t) const { CheckSolution(x, Lsg, t); }

    // void GetDiscError ( const VecDescCL&, scalar_instat_fun_ptr, double) const;
    // void GetDiscError ( scalar_instat_fun_ptr Lsg, double t) const { GetDiscError(x, Lsg, t); }


    DiscSolCL GetSolution() const
      { return DiscSolCL(&x, &GetBndData(), &GetMG(), t); }
};

} // end of namespace DROPS

#include "poisson/instatpoisson.tpp"

#endif
