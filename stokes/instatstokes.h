//**************************************************************************
// File:    instatstokes.h                                                 *
// Content: classes that constitute the stokes-problem                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Sep, 13 2001                                           *
//**************************************************************************

#ifndef DROPS_INSTATSTOKES_H
#define DROPS_INSTATSTOKES_H

#include <vector>
#include "misc/problem.h"
#include "num/stokessolver.h"
#include "stokes/stokes.h"


namespace DROPS
{

typedef BndDataCL<Point3DCL> InstatStokesVelBndDataCL;
typedef BndDataCL<double>    InstatStokesPrBndDataCL;

class InstatStokesBndDataCL
{
  public:
    InstatStokesBndDataCL( Uint numbndseg, const BndCondT* bc_vel, const InstatStokesVelBndDataCL::bnd_val_fun* fun, const BndCondT* bc_pr= 0)
        : Pr( numbndseg, bc_pr), Vel( numbndseg, bc_vel, fun) {}
    InstatStokesBndDataCL( Uint numbndseg, const bool* isneumann, const InstatStokesVelBndDataCL::bnd_val_fun* fun)
        : Pr( numbndseg), Vel( numbndseg, isneumann, fun) {} // deprecated ctor!
    
    const InstatStokesPrBndDataCL     Pr;
    const InstatStokesVelBndDataCL    Vel;   
};

typedef double       (*scalar_instat_fun_ptr)( const Point3DCL&, double);
typedef SVectorCL<3> (*vector_instat_fun_ptr)( const Point3DCL&, double);

template <class Coeff>
class InstatStokesP2P1CL : public ProblemCL<Coeff, InstatStokesBndDataCL>
{
  public:
    typedef ProblemCL<Coeff, InstatStokesBndDataCL> _base;
    typedef typename _base::CoeffCL                 CoeffCL;
    typedef typename _base::BndDataCL               BndDataCL;
    using                                           _base::_MG;
    using                                           _base::_Coeff;
    using                                           _base::_BndData;
    using                                           _base::GetBndData;
    using                                           _base::GetMG;

    typedef P1EvalCL<double, const InstatStokesPrBndDataCL, const VecDescCL>   DiscPrSolCL;
    typedef InstatP2EvalCL<SVectorCL<3>, const InstatStokesVelBndDataCL, 
                                                     const VelVecDescCL> DiscVelSolCL;

    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    double       t;        // time
    VelVecDescCL v;        // velocity
    VecDescCL    p;        // pressure
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A, 
                 B,
                 M;
    
    InstatStokesP2P1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mgb, coeff, bdata), vel_idx( 3, 3), pr_idx( 1), t( 0.) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Vel, match); }
    void CreateNumberingPr ( Uint level, IdxDescCL* idx, match_fun match= 0)
        { CreateNumb( level, *idx, _MG, _BndData.Vel, match); }
    void DeleteNumberingVel( IdxDescCL*);
    void DeleteNumberingPr ( IdxDescCL*);
    
    // Set up matrices and complete rhs
    void SetupSystem( MatDescCL*, VelVecDescCL*, MatDescCL*, VelVecDescCL*, double) const;
    // Set up mass-matrix for pressure-unknowns (P1) -- needed for Uzawa solver
    // Time-independent
    void SetupPrMass(MatDescCL* matM) const;
    // Setup time independent part of system
    void SetupInstatSystem( MatDescCL* A, MatDescCL* B, MatDescCL* M) const;
    // Setup time dependent parts: couplings with bnd unknowns, coefficient f(t)
    // If the function is called with the same vector for some arguments (out of 1, 2, 4), 
    // the vector will contain the sum of the results after the call
    void SetupInstatRhs( VelVecDescCL* vA, VelVecDescCL* vB, VelVecDescCL* vI, double tA, VelVecDescCL* vf, double tf) const;
    // Set up only A.
    void SetupStiffnessMatrix(MatDescCL*) const;
    // Set up mass-matrix for velocity-unknowns (P2) -- needed for MG-Theta-scheme
    // Time-independent
    void SetupMassMatrix(MatDescCL* matI) const;
    // Set initial value for velocities
    void InitVel( VelVecDescCL*, vector_instat_fun_ptr, double t0= 0.) const;

    // Check system and computed solution
    void GetDiscError ( vector_instat_fun_ptr LsgVel, scalar_instat_fun_ptr LsgPr, double t) const;
    void CheckSolution( const VelVecDescCL*, const VecDescCL*,
                        vector_instat_fun_ptr, scalar_instat_fun_ptr, double t) const;

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution() const
        { return DiscPrSolCL( &p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution() const
        { return DiscVelSolCL( &v, &GetBndData().Vel, &GetMG(), t); }

};



class FracStepMatrixCL
{
  private:
    const MatrixCL& _matA;
    const MatrixCL& _matI;
    double          _coeff;
    
  public:
    FracStepMatrixCL( const MatrixCL& I, const MatrixCL& A, double coeff)
        : _matA( A), _matI( I), _coeff( coeff) {}
    
    Uint num_cols() const { return _matA.num_cols(); }
    
    VectorCL operator* (const VectorCL& v) const
    {
        return _coeff*(_matA*v) + _matI*v;
    }
};

class SchurComplNoPcMatrixCL
{
  private:
    const FracStepMatrixCL& _matA;
    const MatrixCL&         _matB;
    double    _tol;
    
  public:
    SchurComplNoPcMatrixCL( const FracStepMatrixCL& A, const MatrixCL& B, double tol)
        : _matA(A), _matB(B), _tol(tol) {}
    friend VectorCL operator*( const SchurComplNoPcMatrixCL& M, const VectorCL& v);
};

} // end of namespace DROPS

#include "stokes/instatstokes.tpp"

#endif
