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

class InstatVelBndSegDataCL
{
  public:
    typedef SVectorCL<3>    bnd_type;
    typedef InstatVelBndSegDataCL self;
    typedef bnd_type (*bnd_val_fun)( const Point3DCL&, double);
    
  private:
    bool         _Neumann;
    bnd_val_fun  _bnd_val;

  public:
    InstatVelBndSegDataCL( bool neumann, bnd_val_fun f)
    : _Neumann(neumann), _bnd_val(f) {}
    
    bool        IsNeumann() const { return _Neumann; }
    bnd_val_fun GetBndFun() const { return _bnd_val; }
    bnd_type    GetBndVal( const Point3DCL& p, double t) const { return _bnd_val( p, t); }
};

class InstatStokesVelBndDataCL
{
  private:
    std::vector<InstatVelBndSegDataCL> _BndData;
    
  public:
    typedef InstatVelBndSegDataCL::bnd_type    bnd_type;
    typedef InstatVelBndSegDataCL::bnd_val_fun bnd_val_fun;
    
    InstatStokesVelBndDataCL(Uint, const bool*, const bnd_val_fun*);

    inline bool IsOnDirBnd( const VertexCL&) const;
    inline bool IsOnNeuBnd( const VertexCL&) const;
    inline bool IsOnDirBnd( const EdgeCL&) const;
    inline bool IsOnNeuBnd( const EdgeCL&) const;
    inline bool IsOnDirBnd( const FaceCL&) const;
    inline bool IsOnNeuBnd( const FaceCL&) const;
    
    inline bnd_type GetDirBndValue( const VertexCL&, double) const;
    inline bnd_type GetNeuBndValue( const VertexCL&, double) const;
    inline bnd_type GetDirBndValue( const EdgeCL&, double) const;
    inline bnd_type GetNeuBndValue( const EdgeCL&, double) const;
    inline bnd_type GetDirBndValue( const FaceCL&, double) const;
    inline bnd_type GetNeuBndValue( const FaceCL&, double) const;

    const InstatVelBndSegDataCL& GetSegData(BndIdxT idx) const { return _BndData[idx]; }
};

inline InstatStokesVelBndDataCL::InstatStokesVelBndDataCL( Uint numbndseg,
                                                           const bool* isneumann,
                                                           const bnd_val_fun* fun)
{
    _BndData.reserve(numbndseg);
    for (Uint i=0; i<numbndseg; ++i)
        _BndData.push_back( InstatVelBndSegDataCL( isneumann[i], fun[i]) );
}



class InstatStokesBndDataCL
{
  public:
    InstatStokesBndDataCL( Uint numbndseg, const bool* isneumann, const InstatVelBndSegDataCL::bnd_val_fun* fun)
        : Pr(), Vel( numbndseg, isneumann, fun) {}
    
    const StokesBndDataCL::PrBndDataCL Pr;
    const InstatStokesVelBndDataCL     Vel;   
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

    typedef P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>   DiscPrSolCL;
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
    void CreateNumberingVel( Uint, IdxDescCL*);
    void CreateNumberingPr ( Uint, IdxDescCL*);
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


//======================================
//        inline functions 
//======================================

inline bool InstatStokesVelBndDataCL::IsOnDirBnd( const VertexCL& v) const
{
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return true;
    return false; 
}        

inline bool InstatStokesVelBndDataCL::IsOnNeuBnd( const VertexCL& v) const
{
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return false;
    return true;
}        

inline bool InstatStokesVelBndDataCL::IsOnDirBnd( const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return true;
    return false;
}

inline bool InstatStokesVelBndDataCL::IsOnNeuBnd( const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return false;
    return true;
}        

inline bool InstatStokesVelBndDataCL::IsOnDirBnd( const FaceCL& f) const
{
    return f.IsOnBoundary() && !_BndData[f.GetBndIdx()].IsNeumann();
}

inline bool InstatStokesVelBndDataCL::IsOnNeuBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && _BndData[f.GetBndIdx()].IsNeumann();
}        




inline InstatStokesVelBndDataCL::bnd_type InstatStokesVelBndDataCL::GetDirBndValue( const VertexCL& v, double t) const
// Returns value of the Dirichlet boundary value. 
// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return _BndData[it->GetBndIdx()].GetBndVal( v.GetCoord(), t);
    throw DROPSErrCL("GetDirBndValue(VertexCL): No Dirichlet Boundary Segment!");
}


inline InstatStokesVelBndDataCL::bnd_type InstatStokesVelBndDataCL::GetNeuBndValue( const VertexCL& v, double t) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( _BndData[it->GetBndIdx()].IsNeumann() )
            return _BndData[it->GetBndIdx()].GetBndVal( v.GetCoord(), t);
    throw DROPSErrCL("GetNeuBndValue(VertexCL): No Neumann Boundary Segment!");
}

inline InstatStokesVelBndDataCL::bnd_type InstatStokesVelBndDataCL::GetDirBndValue( const EdgeCL& e, double t) const
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return _BndData[*it].GetBndVal( GetBaryCenter(e), t);
    throw DROPSErrCL("GetDirBndValue(EdgeCL): No Dirichlet Boundary Segment!");
}

inline InstatStokesVelBndDataCL::bnd_type InstatStokesVelBndDataCL::GetNeuBndValue( const EdgeCL& e, double t) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( _BndData[*it].IsNeumann() )
            return _BndData[*it].GetBndVal( GetBaryCenter(e), t);
    throw DROPSErrCL("GetNeuBndValue(EdgeCL): No Neumann Boundary Segment!");
}

inline InstatStokesVelBndDataCL::bnd_type InstatStokesVelBndDataCL::GetDirBndValue( const FaceCL& f, double t) const
{
    Assert( !_BndData[f.GetBndIdx()].IsNeumann(), DROPSErrCL("GetDirBndValue(FaceCL): No Dirichlet Boundary Segment!"), ~0);
    return _BndData[f.GetBndIdx()].GetBndVal( GetBaryCenter(f), t);
}

inline InstatStokesVelBndDataCL::bnd_type InstatStokesVelBndDataCL::GetNeuBndValue( const FaceCL& f, double t) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    Assert( _BndData[f.GetBndIdx()].IsNeumann(), DROPSErrCL("GetNeuBndValue(FaceCL): No Neumann Boundary Segment!"), ~0);
    return _BndData[f.GetBndIdx()].GetBndVal( GetBaryCenter(f), t);
}


} // end of namespace DROPS

#include "stokes/instatstokes.tpp"

#endif
