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
    typedef NoBndDataCL<>        PrBndDataCL;
    typedef BndDataCL<Point3DCL> VelBndDataCL;
    
    StokesBndDataCL(Uint numbndseg, const bool* isneumann, const VelBndDataCL::bnd_val_fun* fun)
        : Pr(), Vel(numbndseg, isneumann, fun) {}
    StokesBndDataCL( Uint numbndseg, const BndCondT* bc, const VelBndDataCL::bnd_val_fun* fun)
        : Pr(), Vel( numbndseg, bc, fun) {}
    
    const PrBndDataCL  Pr;
    const VelBndDataCL Vel;   
};

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

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL> DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, const VelVecDescCL> DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);
    
    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A, 
                 B;
    
    StokesP2P1CL(const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mgb, coeff, bdata), vel_idx( 3, 3), pr_idx( 1) {}  
    StokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mg, coeff, bdata), vel_idx( 3, 3), pr_idx( 1) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel(Uint, IdxDescCL*);
    void CreateNumberingPr (Uint, IdxDescCL*);
    void DeleteNumberingVel(IdxDescCL*);
    void DeleteNumberingPr (IdxDescCL*);
    
    // Set up matrices and rhs
    void SetupSystem(MatDescCL*, VelVecDescCL*, MatDescCL*, VelVecDescCL*) const;
    void SetupStiffnessMatrix(MatDescCL*) const;
    void SetupMass(MatDescCL*) const;

    // Check system and computed solution
    void GetDiscError (instat_vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, instat_vector_fun_ptr, jacobi_fun_ptr, scalar_fun_ptr) const;

    // work of Joerg :-)
    static double ResidualErrEstimator(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution() const
        { return DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution() const
        { return DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }

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

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>                 DiscPrSolCL;
    typedef P1BubbleEvalCL<SVectorCL<3>, const StokesBndDataCL::VelBndDataCL, const VelVecDescCL> DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);
    
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

    // work of Joerg :-)  Very well then: Let the games begin!
    static double ResidualErrEstimator(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution() const
        { return DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution() const
        { return DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }

};

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
    typedef typename _ProblemCL::DiscPrSolCL DiscPrSolCL;
    typedef typename _ProblemCL::DiscVelSolCL DiscVelSolCL;
    

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

    void Init(const DiscPrSolCL&, const DiscVelSolCL&);

    double GetRelRed() { return _RelReduction; }
    void   SetRelRed(double newred) { _RelReduction= newred; }
    double GetThreshold() { return _Threshold; }
    void   SetThreshold(double newThreshold) { _Threshold= newThreshold; }
    bool   DoesMark() { return _DoMark; }
    void   SwitchMark() { _DoMark= _DoMark ? false : true; }
    bool Estimate(const DiscPrSolCL&, const DiscVelSolCL&);
};

//======================================
//        inline functions 
//======================================


} // end of namespace DROPS

#include "stokes/stokes.tpp"

#endif
