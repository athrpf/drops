#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/discretize.h"
#include "num/fe.h"
#include "misc/problem.h"

using namespace DROPS;

typedef double (*fun_ptr)(const SVectorCL<3>&);
typedef SVectorCL<3> (*vfun_ptr)(const SVectorCL<3>&);

enum  OutputModeT { SILENT, NOISY };


double f(const SVectorCL<3>& p)
{ return p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2] +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]; }

double g(const SVectorCL<3>& p)
{  return p[0] +10.*p[1] +100.*p[2]+1000.; }

double h(const SVectorCL<3>& p)
{  return sin(M_PI*p[0])*sin(M_PI*p[1])*sin(M_PI*p[2]); }

double g2(const DROPS::Point3DCL& p)
{
    return (-1.)*p[0];
}

Point3DCL fv(const SVectorCL<3>& p)
{
    return (p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2]
        +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]); }

void MarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte( 0.5);

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It( mg.GetTriangTetraBegin(maxLevel)),
             ItEnd( mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It) {
        if ( (GetBaryCenter( *It)-Mitte).norm() <= std::max( 0.2, 1.5*std::pow( It->GetVolume(), 1.0/3.0)) )
            It->SetRemoveMark();
    }
}


typedef NoBndDataCL<double> BndCL;
BndCL Bnd;
typedef NoBndDataCL<Point3DCL> VBndCL;
VBndCL VBnd;

/*
template<class T=double>
class NewQuadBaseCL: public std::valarray<T>
{
  public:
    typedef T value_type;
    typedef std::valarray<value_type> base_type;

  private:
    NewQuadBaseCL()                      : base_type()       {}
  protected:
#ifdef VALARRAY_BUG
    NewQuadBaseCL (size_t s)             : base_type( T(),s) {}
#else
    NewQuadBaseCL (size_t s)             : base_type( s)     {}
#endif
    NewQuadBaseCL (size_t s, T c)        : base_type( c, s)  {}
    NewQuadBaseCL (const T* tp, size_t s): base_type( tp, s) {}
    template <class X> // For valarray expression-templates
      NewQuadBaseCL (const X& x)         : base_type( x)     {}
    
  public:

DROPS_ASSIGNMENT_OPS_FOR_VALARRAY_DERIVATIVE(NewQuadBaseCL, T, base_type)

//    template<typename FuncT>
//    NewQuadBaseCL& apply ( FuncT fun)
//      { for (size_t i=0; i<val.size(); ++i) val[i]=fun(val[i]); return *this; }

    friend NewQuadBaseCL<double> dot (const NewQuadBaseCL<Point3DCL>&, const NewQuadBaseCL<Point3DCL>&);
};
*/

template<class T=double>
class NewQuad2CL: public std::valarray<T>
{
  public:
    typedef T value_type;
    typedef std::valarray<T> base_type;
    typedef value_type (*instat_fun_ptr)(const Point3DCL&, double);

    static const Uint NumNodesC= 5;
    static const double Node[NumNodesC][4]; // Stuetzstellen 5*4 doubles
    static const double Wght[NumNodesC];    // Gewichte

    static inline BaryCoordCL // Das kopiert leider.
    GetNode( Uint i) { return Node[i]; }

  protected:
    typedef NewQuad2CL<T> self_;

  public:
    NewQuad2CL(): base_type( value_type(), NumNodesC) {}
    NewQuad2CL(const value_type& t): base_type( t, NumNodesC) {}
    template <class X> // For valarray expression-templates
      NewQuad2CL(const X& x): base_type( x) {}

    NewQuad2CL(const TetraCL&, instat_fun_ptr, double= 0.0);
    NewQuad2CL(const LocalP2CL<value_type>&);
    template <class PFunT> 
      NewQuad2CL(const TetraCL&, const PFunT&, double= 0.0);
    
DROPS_ASSIGNMENT_OPS_FOR_VALARRAY_DERIVATIVE(NewQuad2CL, T, base_type)

    inline self_&
    assign(const TetraCL&, value_type (*)(const Point3DCL&, double) , double= 0.0);
    inline self_&
    assign(const LocalP2CL<value_type>&);
    template <class P2FunT> 
      inline self_&
      assign(const TetraCL&, const P2FunT&, double= 0.0);

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const
    {
        value_type sum= this->sum()/120.;
        return (sum + 0.125*(*this)[NumNodesC-1])*absdet;
    }

    // Folgende Spezialformeln nutzen die spezielle Lage der Stuetzstellen aus
    // zur Annaeherung von \int f*phi,    phi = P1-/P2-Hutfunktion
    T quadP1 (int i, double absdet) const
      { return ((1./120.)*(*this)[i] + (1./30.)*(*this)[4])*absdet; }
    T quadP1 (int i, int j, double absdet) const
      { return (i!=j ? (1./720.)*((*this)[i]+(*this)[j]) + (1./180.)*(*this)[4]
                     : (1./180.)*(*this)[i] + (1./90.)*(*this)[4]  )*absdet;}
    T quadP2 (int i, double absdet) const
    { 
        return (i<4 ? (1./360.)*(*this)[i] - (1./90.)*(*this)[4]
                    : (1./180.)*((*this)[VertOfEdge(i-4,0)]+(*this)[VertOfEdge(i-4,1)]) + (1./45.)*(*this)[4]
               )*absdet;
    }

    T quadP2 (int i, int j, double absdet) const
    { 
        const double valBary= (i<4 ? -0.125 : 0.25)*(j<4 ? -0.125 : 0.25);
        return ((i!=j || i>=4) ? Wght[4]*(*this)[4]*valBary
                               : Wght[4]*(*this)[4]*valBary + Wght[i]*(*this)[i]
               )*absdet;
    }
};

template<class T>
const double NewQuad2CL<T>::Node[NewQuad2CL<T>::NumNodesC][4]= {
    {1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}, {.25,.25,.25,.25}
}; 

template<class T>
const double NewQuad2CL<T>::Wght[NewQuad2CL<T>::NumNodesC]= {
    1./120., 1./120., 1./120., 1./120., 2./15.
};

template<class T>
  inline NewQuad2CL<T>&
  NewQuad2CL<T>::assign(const TetraCL& s, value_type (*f)(const Point3DCL&, double) , double t)
{
    for (Uint i= 0; i<NumNodesC-1; ++i)
        (*this)[i]= f( s.GetVertex( i)->GetCoord(), t);
    (*this)[NumNodesC-1]= f( GetBaryCenter( s), t);
    return *this;
}

template<class T>
  inline NewQuad2CL<T>&
  NewQuad2CL<T>::assign(const LocalP2CL<value_type>& f)
{
    (*this)[std::slice( 0, 4, 1)]= f[std::slice( 0, 4, 1)];
    (*this)[NumNodesC-1]= f( BaryCoordCL( 0.25));
    return *this;
}

template<class T>
  template <class P2FunT> 
    inline NewQuad2CL<T>&
    NewQuad2CL<T>::assign(const TetraCL& s, const P2FunT& f, double t)
{
    const double oldt= f.GetTime();
    const_cast<P2FunT&>( f).SetTime( t);
    for (Uint i= 0; i<NumNodesC-1; ++i)
        (*this)[i]= f.val( *s.GetVertex( i));
    (*this)[NumNodesC-1]= f.val( s, 0.25, 0.25, 0.25);
    const_cast<P2FunT&>( f).SetTime( oldt);
    return *this;
}

template<class T>
  NewQuad2CL<T>::NewQuad2CL(const TetraCL& s,
      value_type (*f)(const Point3DCL&, double), double t)
  : base_type( value_type(), NumNodesC)  
{
    this->assign( s, f, t);
}

template<class T>
  NewQuad2CL<T>::NewQuad2CL(const LocalP2CL<value_type>& f)
  : base_type( value_type(), NumNodesC)  
{
    this->assign( f);
}

template<class T>
  template <class PFunT> 
    NewQuad2CL<T>::NewQuad2CL(const TetraCL& s, const PFunT& f, double t)
  : base_type( value_type(), NumNodesC)  
{
    this->assign( s, f, t);
}


inline std::valarray<double>
dot (const std::valarray<Point3DCL>& a, const std::valarray<Point3DCL>& b)
{
    std::valarray<double> ret; 
    for (size_t i=0; i<a.size(); ++i) 
        ret[i]= dot( a[i], b[i]); 
    return ret; 
}


void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns);
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &Bnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord()));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5));
    }
}

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, vfun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns);
    P2EvalCL<Point3DCL, VBndCL,VecDescBaseCL<VectorCL> > fun( &vd, &VBnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord()));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5));
    }
}

typedef P2EvalCL<double, BndCL, VecDescCL> P2FuncT;
typedef Quad2CL<double> OldQuadT;
typedef NewQuad2CL<double> NewQuadT;
typedef NewQuad2CL<Point3DCL> NewVQuadT;

BndCL theBnd;

double Quadrature( DROPS::MultiGridCL& mg, VecDescCL& vd0, VecDescCL& vd1,
    VecDescCL& vd2)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nOld Quadrature:\n";
    double ret= 0.;
    P2FuncT f( &vd0, &theBnd, &mg);
    P2FuncT g( &vd1, &theBnd, &mg);
    P2FuncT h( &vd2, &theBnd, &mg);
    OldQuadT quad;
    OldQuadT quad0;
    OldQuadT quad1;
    OldQuadT quad2;
    for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); it != end; ++it) {
        const double absdet= 6.0*it->GetVolume();

        quad0.val[0]= f.val( *it, 0.0, 0.0, 0.0);
        quad0.val[1]= f.val( *it, 1.0, 0.0, 0.0);
        quad0.val[2]= f.val( *it, 0.0, 1.0, 0.0);
        quad0.val[3]= f.val( *it, 0.0, 0.0, 1.0);
        quad0.val[4]= f.val( *it, 0.25, 0.25, 0.25);

        quad1.val[0]= g.val( *it, 0.0, 0.0, 0.0);
        quad1.val[1]= g.val( *it, 1.0, 0.0, 0.0);
        quad1.val[2]= g.val( *it, 0.0, 1.0, 0.0);
        quad1.val[3]= g.val( *it, 0.0, 0.0, 1.0);
        quad1.val[4]= g.val( *it, 0.25, 0.25, 0.25);

        quad2.val[0]= g.val( *it, 0.0, 0.0, 0.0);
        quad2.val[1]= g.val( *it, 1.0, 0.0, 0.0);
        quad2.val[2]= g.val( *it, 0.0, 1.0, 0.0);
        quad2.val[3]= g.val( *it, 0.0, 0.0, 1.0);
        quad2.val[4]= g.val( *it, 0.25, 0.25, 0.25);

        quad1*= 0.1;
        quad0+= quad1;
        quad2*= 1./3.;
        quad= quad0 + quad2;
        ret+= quad.quad( absdet);
    }
    return ret;
}


double NewQuadrature(DROPS::MultiGridCL& mg, VecDescCL& vd0, VecDescCL& vd1,
    VecDescCL& vd2)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nNew Quadrature:\n";
    double ret= 0.;
    P2FuncT f( &vd0, &theBnd, &mg);
    P2FuncT g( &vd1, &theBnd, &mg);
    P2FuncT h( &vd2, &theBnd, &mg);
    LocalP2CL<> localg;
    NewQuadT quad;
    NewQuadT quad0;
    NewQuadT quad1;
    NewQuadT quad2;
    for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); it != end; ++it) {
        const double absdet= 6.0*it->GetVolume();

        quad0.assign( *it, f);
        
        localg.assign( *it, g);
        quad1.assign( localg);

        quad2.assign( *it, g);

        quad1*= 0.1;
        quad0+= quad1;
        quad2*= 1./3.;
        quad= quad0 + quad2;
        ret+= quad.quad( absdet);
    }
    return ret;

}


int main ()
{
  try {
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
				DROPS::std_basis<3>(2),
				DROPS::std_basis<3>(3),
				40, 40, 40);
    DROPS::MultiGridCL mg( brick);
    DROPS::IdxDescCL idx( 1,1,0,0);
    idx.TriangLevel= 0;
    idx.NumUnknowns= 0;
    DROPS::CreateNumbOnVertex( idx.GetIdx(), idx.NumUnknowns, 1,
        mg.GetTriangVertexBegin( idx.TriangLevel),
        mg.GetTriangVertexEnd( idx.TriangLevel),
        Bnd);
    DROPS::CreateNumbOnEdge( idx.GetIdx(), idx.NumUnknowns, 1,
        mg.GetTriangEdgeBegin( idx.TriangLevel),
        mg.GetTriangEdgeEnd( idx.TriangLevel),
        Bnd);
    DROPS::VecDescCL vd0( &idx);
    SetFun( vd0, mg, f);
    DROPS::VecDescCL vd1( &idx);
    SetFun( vd1, mg, g);
    DROPS::VecDescCL vd2( &idx);
    SetFun( vd2, mg, h);
    TimerCL time;
    time.Start();
    double q0= Quadrature( mg, vd0, vd1, vd2);
    time.Stop();
    std::cout << "integral: " << q0 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    time.Start();
    double q1= NewQuadrature( mg, vd0, vd1, vd2);
    time.Stop();
    std::cout << "integral: " << q1 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    DROPS::IdxDescCL vidx( 3,3,0,0);
    vidx.TriangLevel= 0;
    vidx.NumUnknowns= 0;
    DROPS::CreateNumbOnVertex( vidx.GetIdx(), vidx.NumUnknowns, 3,
        mg.GetTriangVertexBegin( idx.TriangLevel),
        mg.GetTriangVertexEnd( idx.TriangLevel),
        Bnd);
    DROPS::CreateNumbOnEdge( vidx.GetIdx(), vidx.NumUnknowns, 3,
        mg.GetTriangEdgeBegin( vidx.TriangLevel),
        mg.GetTriangEdgeEnd( vidx.TriangLevel),
        Bnd);
    DROPS::VecDescCL vd3( &vidx);
    SetFun( vd3, mg, fv);
    NewVQuadT vq1, vq2;
    std::cout << "dot: "<< dot( vq1, vq2)[1] << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
