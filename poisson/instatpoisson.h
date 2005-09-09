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


class InstatPoissonBndSegDataCL
{
  public:
    typedef double bnd_type;
    typedef bnd_type (*bnd_val_fun)( const Point2DCL&, double);
    
  private:
    bool        _Neumann;
    bnd_val_fun _bnd_val;
    
  public:
    InstatPoissonBndSegDataCL( bool neumann, bnd_val_fun bnd_val)
        : _Neumann(neumann), _bnd_val(bnd_val) {}

    bool        IsNeumann()                    const { return _Neumann; }
    bnd_val_fun GetBndFun()                    const { return _bnd_val; }
    bnd_type    GetBndVal( const Point2DCL& p, double t) const { return _bnd_val(p,t); }
};


class InstatPoissonBndDataCL
{
  private:
    std::vector<InstatPoissonBndSegDataCL> _BndData;

  public:
    typedef InstatPoissonBndSegDataCL::bnd_type    bnd_type;
    typedef InstatPoissonBndSegDataCL::bnd_val_fun bnd_val_fun;
    typedef InstatPoissonBndSegDataCL              BndSegDataCL;

    InstatPoissonBndDataCL( Uint, const bool*, const bnd_val_fun*);

    inline bool IsOnDirBnd( const VertexCL&) const;
    inline bool IsOnNeuBnd( const VertexCL&) const;
    inline bool IsOnDirBnd( const EdgeCL&) const;
    inline bool IsOnNeuBnd( const EdgeCL&) const;
    inline bool IsOnDirBnd( const FaceCL&) const;
    inline bool IsOnNeuBnd( const FaceCL&) const;
    
    inline bnd_type GetDirBndValue( const VertexCL&, double) const;
//    inline bnd_type GetDirBndValue(const EdgeCL&, double) const;
//    inline bnd_type GetDirBndValue(const FaceCL&, double) const;
    inline bnd_type GetNeuBndValue( const VertexCL&, double) const;

    bnd_val_fun         GetBndFun( BndIdxT i)  const { return _BndData[i].GetBndFun(); }
    const BndSegDataCL& GetSegData( BndIdxT i) const { return _BndData[i]; }
};

inline InstatPoissonBndDataCL::InstatPoissonBndDataCL( Uint numbndseg,
                                          const bool* isneumann,
                                          const bnd_val_fun* fun)
{
    _BndData.reserve( numbndseg);
    for (Uint i=0; i<numbndseg; ++i)
        _BndData.push_back( InstatPoissonBndSegDataCL( isneumann[i], fun[i]) );
}


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
	
    
    InstatPoissonP1CL( const MGBuilderCL& mgb, const CoeffCL& coeff, const BndDataCL& bdata)
        : _base( mgb, coeff, bdata), t( 0.), idx( 1) {}  
				
    // numbering of unknowns
    void CreateNumbering( Uint level, IdxDescCL* idx);
    void DeleteNumbering( IdxDescCL* idx)
        { DeleteNumb( *idx, _MG); }
		
    // set up matrices (which are time independent expressions)
    void SetupInstatSystem( MatDescCL& A, MatDescCL& M) const;
    // set up matrix and couplings with bnd unknowns for convection term
    void SetupConvection( MatDescCL& U, VecDescCL& vU, double t) const;
		
    // Setup time dependent parts: couplings with bnd unknowns, coefficient f(t)
    // If the function is called with the same vector for some arguments, 
    // the vector will contain the sum of the results after the call
    void SetupInstatRhs( VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf) const;
  
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



//======================================
//   declaration of inline functions 
//======================================

inline bool InstatPoissonBndDataCL::IsOnDirBnd(const VertexCL& v) const
{
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return true;
    return false; 
}        

inline bool InstatPoissonBndDataCL::IsOnNeuBnd(const VertexCL& v) const
{
    if ( !v.IsOnBoundary() ) return false;
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return false;
    return true;
}        

inline bool InstatPoissonBndDataCL::IsOnDirBnd(const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return true;
    return false;
}

inline bool InstatPoissonBndDataCL::IsOnNeuBnd(const EdgeCL& e) const
{
    if ( !e.IsOnBoundary() ) return false;
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return false;
    return true;
}        

inline bool InstatPoissonBndDataCL::IsOnDirBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && !_BndData[f.GetBndIdx()].IsNeumann();
}

inline bool InstatPoissonBndDataCL::IsOnNeuBnd(const FaceCL& f) const
{
    return f.IsOnBoundary() && _BndData[f.GetBndIdx()].IsNeumann();
}        


inline InstatPoissonBndDataCL::bnd_type InstatPoissonBndDataCL::GetDirBndValue( const VertexCL& v, double t) const
// Returns value of the Dirichlet boundary value. 
// Expects, that there is any Dirichlet boundary ( IsOnDirBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( !_BndData[it->GetBndIdx()].IsNeumann() )
            return _BndData[it->GetBndIdx()].GetBndVal( it->GetCoord2D(), t);
    throw DROPSErrCL("GetDirBndValue(VertexCL): No Dirichlet Boundary Segment!");
}


inline InstatPoissonBndDataCL::bnd_type InstatPoissonBndDataCL::GetNeuBndValue(const VertexCL& v, double t) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    for (VertexCL::const_BndVertIt it= v.GetBndVertBegin(), end= v.GetBndVertEnd(); it!=end; ++it)
        if ( _BndData[it->GetBndIdx()].IsNeumann() )
            return _BndData[it->GetBndIdx()].GetBndVal( it->GetCoord2D(), t);
    throw DROPSErrCL("GetNeuBndValue(VertexCL): No Neumann Boundary Segment!");
}

/*
inline InstatPoissonBndDataCL::bnd_type InstatPoissonBndDataCL::GetDirBndValue(const EdgeCL& e, double t) const
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( !_BndData[*it].IsNeumann() )
            return _BndData[*it].GetBndVal( GetBaryCenter(e), t);
    throw DROPSErrCL("GetDirBndValue(EdgeCL): No Dirichlet Boundary Segment!");
}

inline InstatPoissonBndDataCL::bnd_type InstatPoissonBndDataCL::GetNeuBndValue(const EdgeCL& e, double t) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    for (const BndIdxT* it= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); it!=end; ++it)
        if ( _BndData[*it].IsNeumann() )
            return _BndData[*it].GetBndVal( GetBaryCenter(e), t);
    throw DROPSErrCL("GetNeuBndValue(EdgeCL): No Neumann Boundary Segment!");
}

inline InstatPoissonBndDataCL::bnd_type InstatPoissonBndDataCL::GetDirBndValue( const FaceCL& f, double t) const
{
    Assert( !_BndData[f.GetBndIdx()].IsNeumann(), DROPSErrCL("GetDirBndValue(FaceCL): No Dirichlet Boundary Segment!"), ~0);
    return _BndData[f.GetBndIdx()].GetBndVal( GetBaryCenter(f), t);
}

inline InstatPoissonBndDataCL::bnd_type InstatPoissonBndDataCL::GetNeuBndValue( const FaceCL& f, double t) const
// Returns value of the Neumann boundary value. 
// Expects, that there is any Neumann boundary ( IsOnNeuBnd(...) == true )
{
    Assert( _BndData[f.GetBndIdx()].IsNeumann(), DROPSErrCL("GetNeuBndValue(FaceCL): No Neumann Boundary Segment!"), ~0);
    return _BndData[f.GetBndIdx()].GetBndVal( GetBaryCenter(f), t);
}
*/

} // end of namespace DROPS

#include "poisson/instatpoisson.tpp"

#endif
