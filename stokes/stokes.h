//**************************************************************************
// File:    stokes.h                                                       *
// Content: classes that constitute the stokes-problem                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 30 2001                                         *
//**************************************************************************

#ifndef DROPS_STOKES_H
#define DROPS_STOKES_H

#include <vector>
#include "misc/problem.h"
#include "num/solver.h"
#include "num/fe.h"
#include "num/discretize.h"


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
    
    const PrBndDataCL  Pr;
    const VelBndDataCL Vel;   
};


typedef StokesBndDataCL::VelBndDataCL StokesVelBndDataCL;
typedef StokesBndDataCL::PrBndDataCL  StokesPrBndDataCL;

typedef SMatrixCL<3, 3> (*jacobi_fun_ptr)(const Point3DCL&);

typedef VecDescBaseCL<VectorCL> VelVecDescCL;

template <class Coeff>
class StokesP2P1CL : public ProblemCL<Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, StokesBndDataCL> _base;
    typedef typename _base::CoeffCL           CoeffCL;
    typedef typename _base::BndDataCL         BndDataCL;
    using                                     _base::_MG;
    using                                     _base::_Coeff;
    using                                     _base::_BndData;
    using                                     _base::GetBndData;
    using                                     _base::GetMG;

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>           DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>           const_DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double);
    
    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    double       t;        // time, 0.0 in the stationary case
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A, 
                 B,
                 M;
    
    StokesP2P1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mgb, coeff, bdata), vel_idx( 3, 3), pr_idx( 1), t( 0.0) {}  
    StokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mg, coeff, bdata), vel_idx( 3, 3), pr_idx( 1), t( 0.0) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Vel, match); }
    void CreateNumberingPr ( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Pr, match); }
    void DeleteNumberingVel(IdxDescCL*);
    void DeleteNumberingPr (IdxDescCL*);
    
    // Set up matrices and complete rhs
    void SetupSystem(MatDescCL*, VelVecDescCL*, MatDescCL*, VelVecDescCL*, double= 0.0) const;
    // Set up only A.
    void SetupStiffnessMatrix(MatDescCL*) const;
    // Set up mass-matrix for pressure-unknowns (P1)
    void SetupPrMass(MatDescCL*) const;
    // Set up mass-matrix for velocity-unknowns (P2) -- needed for MG-Theta-scheme
    // Time-independent
    void SetupMassMatrix(MatDescCL* matI) const;

    // Setup time independent part of system
    void SetupInstatSystem( MatDescCL* A, MatDescCL* B, MatDescCL* M) const;
    // Setup time dependent parts: couplings with bnd unknowns, coefficient f(t)
    // If the function is called with the same vector for some arguments (out of 1, 2, 4), 
    // the vector will contain the sum of the results after the call
    void SetupInstatRhs( VelVecDescCL* vA, VelVecDescCL* vB, VelVecDescCL* vI, double tA, VelVecDescCL* vf, double tf) const;
    // Set initial value for velocities
    void InitVel( VelVecDescCL*, instat_vector_fun_ptr, double t0= 0.) const;

    // Check system and computed solution
    void GetDiscError (instat_vector_fun_ptr LsgVel, instat_scalar_fun_ptr LsgPr, double t= 0.0) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, instat_vector_fun_ptr, jacobi_fun_ptr, scalar_fun_ptr) const;
    // XXX: merge stationary and instationary version
    void CheckSolution(const VelVecDescCL*, const VecDescCL*,
                       instat_vector_fun_ptr, instat_scalar_fun_ptr, double t) const;

    // estimation a la Verfuerth
    static double ResidualErrEstimator(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double t=0.0);

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution()
        { return DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution()
        { return DiscVelSolCL(&v, &GetBndData().Vel, &GetMG(), t); }
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL(&v, &GetBndData().Vel, &GetMG(), t); }

};

template <class Coeff>
class StokesP1BubbleP1CL : public ProblemCL<Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, StokesBndDataCL> _base;
    typedef typename _base::CoeffCL           CoeffCL;
    typedef typename _base::BndDataCL         BndDataCL;
    using                                     _base::_MG;
    using                                     _base::_Coeff;
    using                                     _base::_BndData;
    using                                     _base::GetBndData;
    using                                     _base::GetMG;

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>                 DiscPrSolCL;
    typedef P1BubbleEvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, VelVecDescCL> DiscVelSolCL;
    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>                 const_DiscPrSolCL;
    typedef P1BubbleEvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, const VelVecDescCL> const_DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double);
    
    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A,        // Bad, bad comma.
                 B;
    
    StokesP1BubbleP1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mgb, coeff, bdata), vel_idx( 3, 0, 0, 3), pr_idx( 1) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel(Uint, IdxDescCL*);
    void CreateNumberingPr (Uint, IdxDescCL*);
    void DeleteNumberingVel(IdxDescCL*);
    void DeleteNumberingPr (IdxDescCL*);
    
    // Set up matrices and rhs
    void SetupSystem(MatDescCL*, VelVecDescCL*, MatDescCL*, VelVecDescCL*) const;
    void SetupMass(MatDescCL*) const;

    // Check system and computed solution
    void GetDiscError (instat_vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, instat_vector_fun_ptr, scalar_fun_ptr) const;

    // estimation a la Verfuerth
    static double ResidualErrEstimator(const TetraCL&, const const_DiscPrSolCL&, const const_DiscVelSolCL&, double= 0.0);

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution()
        { return DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution()
        { return DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }
    const_DiscPrSolCL GetPrSolution() const
        { return const_DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    const_DiscVelSolCL GetVelSolution() const
        { return const_DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }

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

//======================================
//        inline functions 
//======================================


} // end of namespace DROPS

#include "stokes/stokes.tpp"

#endif
