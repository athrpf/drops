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


namespace DROPS
{

class VelBndSegDataCL
{
  public:
    typedef SVectorCL<3>    bnd_type;
    typedef VelBndSegDataCL self;
    typedef bnd_type (*bnd_val_fun)(const Point3DCL&);
    
  private:
    bool         _Neumann;
    bnd_val_fun  _bnd_val;

  public:
    VelBndSegDataCL(bool neumann, bnd_val_fun f)
    : _Neumann(neumann), _bnd_val(f) {}
    
    bool        IsNeumann() const { return _Neumann; }
    bnd_val_fun GetBndFun() const { return _bnd_val; }
    bnd_type    GetBndVal(const Point3DCL& p) const { return _bnd_val(p); }
};

class StokesVelBndDataCL
{
  private:
    std::vector<VelBndSegDataCL> _BndData;

    
  public:
    typedef VelBndSegDataCL::bnd_type    bnd_type;
    typedef VelBndSegDataCL::bnd_val_fun bnd_val_fun;
    
    StokesVelBndDataCL(Uint, const bool*, const bnd_val_fun*);

    inline bool IsOnDirBnd(const VertexCL&) const;
    inline bool IsOnNeuBnd(const VertexCL&) const;
    inline bool IsOnDirBnd(const EdgeCL&) const;
    inline bool IsOnNeuBnd(const EdgeCL&) const;
    inline bool IsOnDirBnd(const FaceCL&) const;
    inline bool IsOnNeuBnd(const FaceCL&) const;
    
    inline bnd_type GetDirBndValue(const VertexCL&) const;
    inline bnd_type GetNeuBndValue(const VertexCL&) const;
    inline bnd_type GetDirBndValue(const EdgeCL&) const;
    inline bnd_type GetNeuBndValue(const EdgeCL&) const;
    inline bnd_type GetDirBndValue(const FaceCL&) const;
    inline bnd_type GetNeuBndValue(const FaceCL&) const;

    const VelBndSegDataCL& GetSegData(BndIdxT idx) const { return _BndData[idx]; }
};

inline StokesVelBndDataCL::StokesVelBndDataCL(Uint numbndseg,
                                              const bool* isneumann,
                                              const bnd_val_fun* fun)
{
    _BndData.reserve(numbndseg);
    for (Uint i=0; i<numbndseg; ++i)
        _BndData.push_back( VelBndSegDataCL(isneumann[i], fun[i]) );
}

// This class is more or less a hack to describe boundary conditions for
// pressure
class StokesPrBndDataCL
{
  public:
    // default ctor, dtor, whatever
    typedef double bnd_type;

    static inline bool IsOnDirBnd (const VertexCL&) { return false; }
    static inline bool IsOnNeuBnd (const VertexCL&) { return false; }
    
    static inline bnd_type GetDirBndValue (const VertexCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
};


class StokesBndDataCL
{
  public:
    StokesBndDataCL(Uint numbndseg, const bool* isneumann, const VelBndSegDataCL::bnd_val_fun* fun)
        : Pr(), Vel(numbndseg, isneumann, fun) {}
    
    typedef StokesPrBndDataCL PrBndDataCL;
    typedef StokesVelBndDataCL VelBndDataCL;

    const StokesPrBndDataCL  Pr;
    const StokesVelBndDataCL Vel;   
};

typedef double (*scalar_fun_ptr)(const Point3DCL&);
typedef SVectorCL<3> (*vector_fun_ptr)(const Point3DCL&);

template <class MGB, class Coeff>
class StokesP2P1CL : public ProblemCL<MGB, Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<MGB, Coeff, StokesBndDataCL> _base;
    typedef typename _base::MultiGridBuilderCL     MultiGridBuilderCL;
    typedef typename _base::CoeffCL                CoeffCL;
    typedef typename _base::BndDataCL              BndDataCL;
    using                                          _base::_MG;
    using                                          _base::_Coeff;
    using                                          _base::_BndData;
    using                                          _base::GetBndData;
    using                                          _base::GetMG;

    typedef VecDescBaseCL<VectorCL> VelVecDescCL;
    typedef P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>           DiscPrSolCL;
    typedef P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);
    
    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A, 
                 B;
    
    StokesP2P1CL(const MultiGridBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mgb, coeff, bdata) {}  
    StokesP2P1CL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mg, coeff, bdata) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel(Uint, IdxDescCL*);
    void CreateNumberingPr (Uint, IdxDescCL*);
    void DeleteNumberingVel(IdxDescCL*);
    void DeleteNumberingPr (IdxDescCL*);
    
    // Set up matrices and rhs
    void SetupSystem(MatDescCL*, VelVecDescCL*, MatDescCL*, VelVecDescCL*) const;
    void SetupMass(MatDescCL*) const;

    // Check system and computed solution
    void GetDiscError (vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, vector_fun_ptr, scalar_fun_ptr) const;

    // work of Joerg :-)
    static double ResidualErrEstimator(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);

    // Get solutions as FE-functions
    DiscPrSolCL GetPrSolution() const
        { return DiscPrSolCL(&p, &GetBndData().Pr, &GetMG()); }
    DiscVelSolCL GetVelSolution() const
        { return DiscVelSolCL(&v, &GetBndData().Vel, &GetMG()); }

};

template <class MGB, class Coeff>
class StokesP1BubbleP1CL : public ProblemCL<MGB, Coeff, StokesBndDataCL>
{
  public:
    typedef ProblemCL<MGB, Coeff, StokesBndDataCL> _base;
    typedef typename _base::MultiGridBuilderCL     MultiGridBuilderCL;
    typedef typename _base::CoeffCL                CoeffCL;
    typedef typename _base::BndDataCL              BndDataCL;
    using                                          _base::_MG;
    using                                          _base::_Coeff;
    using                                          _base::_BndData;
    using                                          _base::GetBndData;
    using                                          _base::GetMG;

    typedef VecDescBaseCL<VectorCL> VelVecDescCL;
    typedef P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>                 DiscPrSolCL;
    typedef P1BubbleEvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> DiscVelSolCL;

    typedef double (*est_fun)(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);
    
    IdxDescCL    vel_idx;  // for velocity unknowns
    IdxDescCL    pr_idx;   // for pressure unknowns
    VelVecDescCL v;
    VecDescCL    p;
    VelVecDescCL b;
    VecDescCL    c;
    MatDescCL    A,        // Bad, bad comma.
                 B;
    
    StokesP1BubbleP1CL(const MultiGridBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base(mgb, coeff, bdata) {}  

    // Create and delete numbering of unknowns
    void CreateNumberingVel(Uint, IdxDescCL*);
    void CreateNumberingPr (Uint, IdxDescCL*);
    void DeleteNumberingVel(IdxDescCL*);
    void DeleteNumberingPr (IdxDescCL*);
    
    // Set up matrices and rhs
    void SetupSystem(MatDescCL*, VelVecDescCL*, MatDescCL*, VelVecDescCL*) const;
    void SetupMass(MatDescCL*) const;

    // Check system and computed solution
    void GetDiscError (vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const;
    void CheckSolution(const VelVecDescCL*, const VecDescCL*, vector_fun_ptr, scalar_fun_ptr) const;

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
    typedef VecDescBaseCL<VectorCL> VelVecDescCL;
//    typedef P1EvalCL<double, const PrBndDataCL, const VecDescCL>                 DiscPrSolCL;
// TODO: New template parameter???
//    typedef P1BubbleEvalCL<SVectorCL<3>, const VelBndDataCL, const VelVecDescCL> DiscVelSolCL;
//    typedef P2EvalCL<SVectorCL<3>, const VelBndDataCL, const VelVecDescCL> DiscVelSolCL;
//  DONE:
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

class SchurComplMatrixCL
{
  private:
    const MatrixCL& _matA;
    const MatrixCL& _matB;
    double          _tol;
    SSORPcCL        _pc;

  public:
    SchurComplMatrixCL(const MatrixCL& A, const MatrixCL& B, double tol, double om)
        : _matA(A), _matB(B), _tol(tol), _pc(om) {}
    friend VectorCL operator* (const SchurComplMatrixCL& M, const VectorCL& v);
};

void Uzawa(const MatrixCL& A, const MatrixCL& B, const MatrixCL& I, VectorCL& x, VectorCL& y, const VectorCL& f, const VectorCL& g, 
           double tau, int& max_iter, double& tol, Uint inner_iter, double inner_iter_tol);

//double GradPrEstimator(const TetraCL&, const DiscPrSolCL&, const DiscVelSolCL&);



//======================================
//        inline functions 
//======================================

inline bool StokesVelBndDataCL::IsOnDirBnd(const VertexCL& v) const
{
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return true;
    return false; 
}        

inline bool StokesVelBndDataCL::IsOnNeuBnd(const VertexCL& v) const
{
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return false;
    return true;
}        

inline bool StokesVelBndDataCL::IsOnDirBnd(const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return true;
    return false;
}

inline bool StokesVelBndDataCL::IsOnNeuBnd(const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return false;
    return true;
}        

inline bool StokesVelBndDataCL::IsOnDirBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && !_BndData[f.GetBndIdx()].IsNeumann();
}

inline bool StokesVelBndDataCL::IsOnNeuBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && _BndData[f.GetBndIdx()].IsNeumann();
}        




inline StokesVelBndDataCL::bnd_type StokesVelBndDataCL::GetDirBndValue(const VertexCL& v) const
// Returns value of the Dirichlet boundary value. 
// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return _BndData[it->GetBndIdx()].GetBndVal( v.GetCoord() );
    throw DROPSErrCL("GetDirBndValue(VertexCL): No Dirichlet Boundary Segment!");
}


inline StokesVelBndDataCL::bnd_type StokesVelBndDataCL::GetNeuBndValue(const VertexCL& v) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( _BndData[it->GetBndIdx()].IsNeumann() )
            return _BndData[it->GetBndIdx()].GetBndVal( v.GetCoord() );
    throw DROPSErrCL("GetNeuBndValue(VertexCL): No Neumann Boundary Segment!");
}

inline StokesVelBndDataCL::bnd_type StokesVelBndDataCL::GetDirBndValue(const EdgeCL& e) const
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return _BndData[*it].GetBndVal( GetBaryCenter(e) );
    throw DROPSErrCL("GetDirBndValue(EdgeCL): No Dirichlet Boundary Segment!");
}

inline StokesVelBndDataCL::bnd_type StokesVelBndDataCL::GetNeuBndValue(const EdgeCL& e) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( _BndData[*it].IsNeumann() )
            return _BndData[*it].GetBndVal( GetBaryCenter(e) );
    throw DROPSErrCL("GetNeuBndValue(EdgeCL): No Neumann Boundary Segment!");
}

inline StokesVelBndDataCL::bnd_type StokesVelBndDataCL::GetDirBndValue(const FaceCL& f) const
{
    Assert( !_BndData[f.GetBndIdx()].IsNeumann(), DROPSErrCL("GetDirBndValue(FaceCL): No Dirichlet Boundary Segment!"), ~0);
    return _BndData[f.GetBndIdx()].GetBndVal( GetBaryCenter(f) );
}

inline StokesVelBndDataCL::bnd_type StokesVelBndDataCL::GetNeuBndValue(const FaceCL& f) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    Assert( _BndData[f.GetBndIdx()].IsNeumann(), DROPSErrCL("GetNeuBndValue(FaceCL): No Neumann Boundary Segment!"), ~0);
    return _BndData[f.GetBndIdx()].GetBndVal( GetBaryCenter(f) );
}


} // end of namespace DROPS

#include "stokes/stokes.tpp"

#endif
