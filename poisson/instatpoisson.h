//**************************************************************************
// File:    instatpoisson.h                                                *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - November, 11 2002                                      *
//**************************************************************************

#ifndef DROPS_INSTATPOISSON_H
#define DROPS_INSTATPOISSON_H

#include <vector>
#include "misc/problem.h"
#include "num/solver.h"
#include "num/fe.h"


namespace DROPS
{

typedef BndSegDataCL<double> InstatPoissonBndSegDataCL;
typedef BndDataCL<double>    InstatPoissonBndDataCL;
typedef double	(*scalar_instat_fun_ptr)( const Point3DCL&, double);
typedef double 	(*scalar_fun_ptr)( const Point3DCL&);

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
class InstatPoissonP1CL : public ProblemCL<Coeff, InstatPoissonBndDataCL>
{
  private:
    bool adjoint_;

  public:
    typedef ProblemCL<Coeff, InstatPoissonBndDataCL> _base;
    typedef typename _base::BndDataCL                BndDataCL;
    typedef typename _base::CoeffCL                  CoeffCL;
    using _base::GetBndData;
    using _base::GetMG;
    using _base::_BndData;
    using _base::_MG;
    using _base::_Coeff;

    typedef P1EvalCL<double, const BndDataCL, const VecDescCL> DiscSolCL;
    typedef double (*est_fun)(const TetraCL&, const VecDescCL&, const BndDataCL&);

    double     	t;        // time

    IdxDescCL idx;
    VecDescCL x;
    VecDescCL b;
    MatDescCL A;
    MatDescCL M;
    MatDescCL U;
	

    InstatPoissonP1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata, bool adj=false)
        : _base( mgb, coeff, bdata), adjoint_( adj), t( 0.), idx( 1) {}
				
    // numbering of unknowns
    void CreateNumbering( Uint level, IdxDescCL* idx);
    void DeleteNumbering( IdxDescCL* idx)
        { DeleteNumb( *idx, _MG); }
		
    // set up matrices (M is time independent)
    void SetupInstatSystem( MatDescCL& A, MatDescCL& M, double tA) const;
    // set up matrix and couplings with bnd unknowns for convection term
    void SetupConvection( MatDescCL& U, VecDescCL& vU, double t) const;
		
    // Setup time dependent parts: couplings with bnd unknowns, coefficient f(t)
    // If the function is called with the same vector for some arguments,
    // the vector will contain the sum of the results after the call
    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf) const;
    // Setup special source term including the gradient of a given P1 function
    void SetupGradSrc( VecDescCL& src, scalar_instat_fun_ptr T, scalar_instat_fun_ptr dalpha, double t= 0.) const;

    // create prolongation
    void SetupProlongation( MatDescCL& P, IdxDescCL* cIdx, IdxDescCL* fIdx) const;

    // Set initial value
    void Init( VecDescCL&, scalar_instat_fun_ptr, double t0= 0.) const;
		
    // check computed solution etc.
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
