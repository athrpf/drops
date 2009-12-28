/// \file poissonsolverfactory.h
/// \brief Solver for Poisson problem with P2 functions
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Marcus Soemers, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

/** We solve \f$ -\Delta u = f\;\mbox{in}\; \Omega:=[0,1]^3 \f$ for the given
    solution \f$ u(x,y,z):= 64 \cdot xyz (1-x) (1-y) (1-z) \f$, i.e. homogeneous
    Dirichlet conditions are used. A uniform tetrahedral grid is applied as
    a triangulation of \f$ \Omega \f$. GMRES is used as a linear solver for the
    discretized linear equation system. Note, that CG-type methods can be used
    as well because the resulting linear equation system is s.p.d. However,
    since this program acts as a base performance test, GMRES is used here.
*/

#ifndef POISSONSOLVERFACTORY_H_
#define POISSONSOLVERFACTORY_H_

namespace DROPS
{

template <class ParamsT>
class PoissonSolverFactoryHelperCL
{
  public:
    bool MGUsed ( const ParamsT& C) const
    {
        return ( C.pos_Method / 100 == 1 || C.pos_Method % 10 == 1);
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

#ifndef _PAR
template <class ParamsT, class ProlongationT= MLMatrixCL>
class PoissonSolverFactoryCL
{
  private:
	MLIdxDescCL& idx_;
    ParamsT& C_;
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
    PoissonSolverFactoryCL( ParamsT& C, MLIdxDescCL& idx);
    ~PoissonSolverFactoryCL() {}

    /// Returns pointer to prolongation for velocity
    ProlongationT* GetProlongation();
    PoissonSolverBaseCL* CreatePoissonSolver();

};

template <class ParamsT, class ProlongationT>
PoissonSolverFactoryCL<ParamsT, ProlongationT>::
    PoissonSolverFactoryCL(ParamsT& C, MLIdxDescCL& idx)
    : idx_(idx), C_( C), prolongptr_( 0), JACPc_( C_.pos_Relax), SSORPc_( C_.pos_Relax),
        jorsmoother_( C_.pos_Relax), gssmoother_( C_.pos_Relax), sgssmoother_( C_.pos_Relax), sorsmoother_( C_.pos_Relax), ssorsmoother_( C_.pos_Relax),
        coarsesolversymm_( SSORPc_, 500, 1e-6, true),
        MGSolversymmJOR_( jorsmoother_, coarsesolversymm_, C_.pos_Iter, C_.pos_Tol, false, C_.pos_SmoothingSteps, C_.pos_NumLvl),
        MGSolversymmGS_( gssmoother_, coarsesolversymm_, C_.pos_Iter, C_.pos_Tol, false, C_.pos_SmoothingSteps, C_.pos_NumLvl),
        MGSolversymmSGS_( sgssmoother_, coarsesolversymm_, C_.pos_Iter, C_.pos_Tol, false, C_.pos_SmoothingSteps, C_.pos_NumLvl),
        MGSolversymmSOR_( sorsmoother_, coarsesolversymm_, C_.pos_Iter, C_.pos_Tol, C_.pos_RelativeErr, C_.pos_SmoothingSteps, C_.pos_NumLvl),
        MGSolversymmSSOR_( ssorsmoother_, coarsesolversymm_, C_.pos_Iter, C_.pos_Tol, C_.pos_RelativeErr, C_.pos_SmoothingSteps, C_.pos_NumLvl),
        GMResSolver_( JACPc_, C_.pos_Restart, C_.pos_Iter, C_.pos_Tol, C_.pos_RelativeErr),
        GMResSolverSSOR_( SSORPc_, C_.pos_Restart, C_.pos_Iter, C_.pos_Tol, C_.pos_RelativeErr),
        PCGSolver_( SSORPc_, C_.pos_Iter, C_.pos_Tol, C_.pos_RelativeErr)
        {}

template <class ParamsT, class ProlongationT>
PoissonSolverBaseCL* PoissonSolverFactoryCL<ParamsT, ProlongationT>::CreatePoissonSolver()
{
    PoissonSolverBaseCL* poissonsolver = 0;
    switch (C_.pos_Method)
    {
		case  102 : {
            poissonsolver = new PoissonSolverCL<MGSolversymmJORT>( MGSolversymmJOR_);
            prolongptr_ = MGSolversymmJOR_.GetProlongation();
		} break;
		case  103 : {
            poissonsolver = new PoissonSolverCL<MGSolversymmSSORT>( MGSolversymmSSOR_);
            prolongptr_ = MGSolversymmSSOR_.GetProlongation();
		} break;
		case  104 : {
            poissonsolver = new PoissonSolverCL<MGSolversymmGST>( MGSolversymmGS_);
            prolongptr_ = MGSolversymmGS_.GetProlongation();
		} break;
		case  105 : {
            poissonsolver = new PoissonSolverCL<MGSolversymmSGST>( MGSolversymmSGS_);
            prolongptr_ = MGSolversymmSGS_.GetProlongation();
		} break;
		case  106 : {
            poissonsolver = new PoissonSolverCL<MGSolversymmSORT>( MGSolversymmSOR_);
            prolongptr_ = MGSolversymmSOR_.GetProlongation();
		} break;
		case  302 : poissonsolver = new PoissonSolverCL<GMResSolverT>( GMResSolver_);  break;
		case  303 : poissonsolver = new PoissonSolverCL<GMResSolverSSORT>( GMResSolverSSOR_);  break;
		case  203 : poissonsolver = new PoissonSolverCL<PCGSolverT>( PCGSolver_); break;
        default: throw DROPSErrCL("PoissonSolverFactoryCL: Unknown Poisson solver");
    }
    return poissonsolver;
}

template <class ParamsT, class ProlongationT>
ProlongationT* PoissonSolverFactoryCL<ParamsT, ProlongationT>::GetProlongation()
{
    return prolongptr_;
}

#else
template <class ParamsT, class ProlongationT= MLMatrixCL>
class PoissonSolverFactoryCL
{
  private:
	MLIdxDescCL & idx_;
    ParamsT& C_;

// generic preconditioners
    ParJac0CL  JACPc_;

     //JAC-GMRes
    typedef ParPreGMResSolverCL<ParJac0CL> GMResSolverT;
    GMResSolverT GMResSolver_;

  public:
    PoissonSolverFactoryCL( ParamsT& C, MLIdxDescCL& idx);
    ~PoissonSolverFactoryCL() {}

    /// Returns pointer to prolongation for velocity
    ProlongationT* GetProlongation() {return 0;}
    PoissonSolverBaseCL* CreatePoissonSolver();

};

template <class ParamsT, class ProlongationT>
PoissonSolverFactoryCL<ParamsT, ProlongationT>::
    PoissonSolverFactoryCL(ParamsT& C, MLIdxDescCL& idx)
    : idx_(idx), C_( C), JACPc_( idx_.GetFinest(), C_.pos_Relax),
        GMResSolver_( C_.pos_Restart, C_.pos_Iter, C_.pos_Tol, idx_.GetFinest(), JACPc_, C_.pos_RelativeErr)
        {}

template <class ParamsT, class ProlongationT>
PoissonSolverBaseCL* PoissonSolverFactoryCL<ParamsT, ProlongationT>::CreatePoissonSolver()
{
    PoissonSolverBaseCL* poissonsolver = 0;
    switch (C_.pos_Method)
    {
		case  302 : poissonsolver = new PoissonSolverCL<GMResSolverT>( GMResSolver_);  break;
        default: throw DROPSErrCL("PoissonSolverFactoryCL: Unknown Poisson solver");
    }
    return poissonsolver;
}

#endif
} //end of namespace DROPS

#endif /* POISSONSOLVERFACTORY_H_ */
