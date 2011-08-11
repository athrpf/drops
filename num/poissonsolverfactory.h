/// \file Poissonsolverfactory.h
/// \brief Solver for Poisson problems functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Marcus Soemers, Yuanjun Zhang, Thorolf Schulte; SC RWTH Aachen: Oliver Fortmeier

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


#ifndef POISSONSOLVERFACTORY_H_
#define POISSONSOLVERFACTORY_H_

#include "num/solver.h"
#include "misc/params.h"
#ifdef _HYPRE
#include "num/hypre.h"
#endif

namespace DROPS
{


class PoissonSolverFactoryHelperCL
{
  public:
    bool MGUsed ( ParamCL& P)
    {
        const int PM = P.get<int>(std::string("Poisson.Method"));
        return ( PM / 100 == 1 || PM % 10 == 1);
    }
};


class PoissonSolverBaseCL : public SolverBaseCL
{
  public:
    PoissonSolverBaseCL (int maxiter, double tol, bool rel= false, std::ostream* output= 0)
        : SolverBaseCL(maxiter, tol, rel, output){}
    virtual void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b) = 0;
    virtual void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b) = 0;
};

template <class SolverT>
class PoissonSolverCL : public PoissonSolverBaseCL
{
  private:
    SolverT& solver_;

  public:
    PoissonSolverCL( SolverT& solver) : PoissonSolverBaseCL( -1, -1.0), solver_( solver) {}
    void Solve( const MatrixCL& A, VectorCL& x, const VectorCL& b) { solver_.Solve( A, x, b); }
    void Solve( const MLMatrixCL& A, VectorCL& x, const VectorCL& b) { solver_.Solve( A, x, b); }
	// We overwrite these functions.
    void   SetTol     (double tol) { solver_.SetTol( tol); }
    void   SetMaxIter (int iter)   { solver_.SetMaxIter( iter); }
    void   SetRelError(bool rel)   { solver_.SetRelError( rel); }

    double GetTol     () const { return solver_.GetTol(); }
    int    GetMaxIter () const { return solver_.GetMaxIter(); }
    double GetResid   () const { return solver_.GetResid(); }
    int    GetIter    () const { return solver_.GetIter(); }
    bool   GetRelError() const { return solver_.GetRelError(); }
};

/// \brief Create a Poisson-Solver (Design-Pattern: Factory class)
/// Construction of a Poisson solver, e.g. CG and SSOR preconditioner: 2*100 + 3. Note: Not all combinations are implemented!
/**
    <table border="3">
    <tr><th> no </th><th> Poisson-Solver-Type </th><th> Type of PC/smoother  </th></tr>
    <tr><td>  0 </td><td>                     </td><td> DummyPc              </td></tr>
    <tr><td>  1 </td><td> MultiGrid V-cycle   </td><td>                      </td></tr>
    <tr><td>  2 </td><td> Preconditioned CG   </td><td> JOR                  </td></tr>
    <tr><td>  3 </td><td> GMRes               </td><td> SSOR                 </td></tr>
    <tr><td>  4 </td><td> Hypre-AMG           </td><td> GS                   </td></tr>
    <tr><td>  5 </td><td>                     </td><td> SGS                  </td></tr>
    <tr><td>  6 </td><td>                     </td><td> SOR                  </td></tr>
    <tr><td>  7 </td><td>                     </td><td>                      </td></tr>
    </table>*/
#ifndef _PAR
template <class ProlongationT= MLMatrixCL>
class PoissonSolverFactoryCL
{
  private:
    ParamCL& P_;
    MLIdxDescCL& idx_;
    ProlongationT* prolongptr_;

// generic preconditioners
    JACPcCL  JACPc_;
    SSORPcCL SSORPc_;

    // MultiGrid symm.
    JORsmoothCL  jorsmoother_;   // Jacobi
    GSsmoothCL   gssmoother_;    // Gauss-Seidel
    SGSsmoothCL  sgssmoother_;   // symmetric Gauss-Seidel
    SORsmoothCL  sorsmoother_;   // Gauss-Seidel with over-relaxation
    SSORsmoothCL ssorsmoother_;  // symmetric Gauss-Seidel with over-relaxation
    PCG_SsorCL   coarsesolversymm_;
    typedef MGSolverCL<JORsmoothCL, PCG_SsorCL, ProlongationT> MGSolversymmJORT;
    MGSolversymmJORT MGSolversymmJOR_;
    typedef MGSolverCL<GSsmoothCL, PCG_SsorCL, ProlongationT> MGSolversymmGST;
    MGSolversymmGST MGSolversymmGS_;
    typedef MGSolverCL<SGSsmoothCL, PCG_SsorCL, ProlongationT> MGSolversymmSGST;
    MGSolversymmSGST MGSolversymmSGS_;
    typedef MGSolverCL<SORsmoothCL, PCG_SsorCL, ProlongationT> MGSolversymmSORT;
    MGSolversymmSORT MGSolversymmSOR_;
    typedef MGSolverCL<SSORsmoothCL, PCG_SsorCL, ProlongationT> MGSolversymmSSORT;
    MGSolversymmSSORT MGSolversymmSSOR_;

    //JAC-GMRes
    typedef GMResSolverCL<JACPcCL> GMResSolverT;
    GMResSolverT GMResSolver_;
    typedef GMResSolverCL<SSORPcCL> GMResSolverSSORT;
    GMResSolverSSORT GMResSolverSSOR_;

    //PCG
    typedef PCGSolverCL<SSORPcCL> PCGSolverT;
    PCGSolverT PCGSolver_;

  public:
    PoissonSolverFactoryCL( ParamCL& P, MLIdxDescCL& idx);
    ~PoissonSolverFactoryCL() {}

    /// Returns pointer to prolongation for velocity
    ProlongationT* GetProlongation();
    PoissonSolverBaseCL* CreatePoissonSolver();

};


template <class ProlongationT>
PoissonSolverFactoryCL<ProlongationT>::
    PoissonSolverFactoryCL(ParamCL& P, MLIdxDescCL& idx)
    : P_(P), idx_(idx), prolongptr_( 0), JACPc_( P.get<double>("Poisson.Relax")), SSORPc_( P.get<double>("Poisson.Relax")),
        jorsmoother_( P.get<double>("Poisson.Relax")), gssmoother_( P.get<double>("Poisson.Relax")), sgssmoother_( P.get<double>("Poisson.Relax")), sorsmoother_( P.get<double>("Poisson.Relax")), ssorsmoother_( P.get<double>("Poisson.Relax")),
        coarsesolversymm_( SSORPc_, 500, 1e-6, true),
        MGSolversymmJOR_( jorsmoother_, coarsesolversymm_, P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), false, P.get<int>("Poisson.SmoothingSteps"), P.get<int>("Poisson.NumLvl")),
        MGSolversymmGS_( gssmoother_, coarsesolversymm_, P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), false, P.get<int>("Poisson.SmoothingSteps"), P.get<int>("Poisson.NumLvl")),
        MGSolversymmSGS_( sgssmoother_, coarsesolversymm_, P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), false, P.get<int>("Poisson.SmoothingSteps"), P.get<int>("Poisson.NumLvl")),
        MGSolversymmSOR_( sorsmoother_, coarsesolversymm_, P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), P.get<double>("Poisson.RelativeErr"), P.get<int>("Poisson.SmoothingSteps"), P.get<int>("Poisson.NumLvl")),
        MGSolversymmSSOR_( ssorsmoother_, coarsesolversymm_, P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), P.get<double>("Poisson.RelativeErr"), P.get<int>("Poisson.SmoothingSteps"), P.get<int>("Poisson.NumLvl")),
        GMResSolver_( JACPc_, P.get<int>("Poisson.Restart"), P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), P.get<double>("Poisson.RelativeErr")),
        GMResSolverSSOR_( SSORPc_, P.get<int>("Poisson.Restart"), P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), P.get<double>("Poisson.RelativeErr")),
        PCGSolver_( SSORPc_, P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), P.get<double>("Poisson.RelativeErr"))
        {}

template <class ProlongationT>
PoissonSolverBaseCL* PoissonSolverFactoryCL<ProlongationT>::CreatePoissonSolver()
{
    PoissonSolverBaseCL* Poissonsolver = 0;
    switch (P_.get<int>("Poisson.Method"))
    {
        case  102 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmJORT>( MGSolversymmJOR_);
            prolongptr_ = MGSolversymmJOR_.GetProlongation();
        } break;
        case  103 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmSSORT>( MGSolversymmSSOR_);
            prolongptr_ = MGSolversymmSSOR_.GetProlongation();
        } break;
        case  104 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmGST>( MGSolversymmGS_);
            prolongptr_ = MGSolversymmGS_.GetProlongation();
        } break;
        case  105 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmSGST>( MGSolversymmSGS_);
            prolongptr_ = MGSolversymmSGS_.GetProlongation();
        } break;
        case  106 : {
            Poissonsolver = new PoissonSolverCL<MGSolversymmSORT>( MGSolversymmSOR_);
            prolongptr_ = MGSolversymmSOR_.GetProlongation();
        } break;
        case  302 : Poissonsolver = new PoissonSolverCL<GMResSolverT>( GMResSolver_);  break;
        case  303 : Poissonsolver = new PoissonSolverCL<GMResSolverSSORT>( GMResSolverSSOR_);  break;
        case  203 : Poissonsolver = new PoissonSolverCL<PCGSolverT>( PCGSolver_); break;
        default: throw DROPSErrCL("PoissonSolverFactoryCL: Unknown Poisson solver");
    }
    return Poissonsolver;
}

template <class ProlongationT>
ProlongationT* PoissonSolverFactoryCL<ProlongationT>::GetProlongation()
{
    return prolongptr_;
}

#else
template <class ProlongationT= MLMatrixCL>
class PoissonSolverFactoryCL
{
  private:
    MLIdxDescCL & idx_;
    ParamCL& P_;

    // generic preconditioners
    ParJac0CL  JACPc_;
    ParDummyPcCL DummyPC_;

    //JAC-PCG
    typedef ParPCGSolverCL<ParJac0CL> JacPCGSolverT;
    JacPCGSolverT JacPCGSolver_;
    //CG without precondition
    typedef ParCGSolverCL CGSolverT;
    CGSolverT CGSolver_;

    //JAC-GMRes
    typedef ParPreGMResSolverCL<ParJac0CL> JacGMResSolverT;
    JacGMResSolverT JacGMResSolver_;
    //Dummy-GMRes
    typedef ParPreGMResSolverCL<ParDummyPcCL> DummyGMResSolverT;
    DummyGMResSolverT DummyGMResSolver_;

#ifdef _HYPRE
     //Algebraic MG solver
    typedef HypreAMGSolverCL AMGSolverT;
    AMGSolverT   hypreAMG_;
#endif

  public:
    PoissonSolverFactoryCL( ParamCL& P, MLIdxDescCL& idx);
    ~PoissonSolverFactoryCL() {}

    /// Returns pointer to prolongation for velocity
    ProlongationT* GetProlongation() {return 0;}
    PoissonSolverBaseCL* CreatePoissonSolver();

};

template <class ProlongationT>
PoissonSolverFactoryCL<ProlongationT>::
    PoissonSolverFactoryCL(ParamCL& P, MLIdxDescCL& idx)
    : idx_(idx), P_(P), JACPc_( idx_.GetFinest(), P.get<double>("Poisson.Relax")), DummyPC_(idx_.GetFinest()),
      JacPCGSolver_( P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), idx_.GetFinest(), JACPc_, P.get<double>("Poisson.RelativeErr")),
      CGSolver_( P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), idx_.GetFinest(), P.get<double>("Poisson.RelativeErr")),
      JacGMResSolver_( P.get<int>("Poisson.Restart"), P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), idx_.GetFinest(), JACPc_, P.get<double>("Poisson.RelativeErr")),
      DummyGMResSolver_( P.get<int>("Poisson.Restart"), P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"), idx_.GetFinest(), DummyPC_, P.get<double>("Poisson.RelativeErr"))
#ifdef _HYPRE
      , hypreAMG_( idx.getFinest(), P.get<int>("Poisson.Iter"), P.get<double>("Poisson.Tol"))
#endif
        {}

template <class ProlongationT>
PoissonSolverBaseCL* PoissonSolverFactoryCL<ProlongationT>::CreatePoissonSolver()
{
    PoissonSolverBaseCL* Poissonsolver = 0;
    switch (P_.get<int>("Poisson.Method"))
    {
        case 200 : Poissonsolver = new PoissonSolverCL<CGSolverT>( CGSolver_); break;
        case 202 : Poissonsolver = new PoissonSolverCL<JacPCGSolverT>( JacPCGSolver_); break;
        case 300 : Poissonsolver = new PoissonSolverCL<DummyGMResSolverT>( DummyGMResSolver_); break;
        case 302 : Poissonsolver = new PoissonSolverCL<JacGMResSolverT>( JacGMResSolver_);  break;
        case 400 : 
#ifdef _HYPRE
            Poissonsolver = new PoissonSolverCL<AMGSolverT>(   hypreAMG_); break;
#else
            throw DROPSErrCL("PoissonSolverFactoryCL::CreatePoissonSolver: Hypre not found, see the Wiki system for help"); break;
#endif
        default: throw DROPSErrCL("PoissonSolverFactoryCL: Unknown Poisson solver");
    }
    return Poissonsolver;
}

#endif
} //end of namespace DROPS

#endif /* PoissonSOLVERFACTORY_H_ */
