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
{  return std::sin(M_PI*p[0])*std::sin(M_PI*p[1])*std::sin(M_PI*p[2]); }

double g2(const DROPS::Point3DCL& p)
{
    return (-1.)*p[0];
}

Point3DCL fv(const SVectorCL<3>& p)
{
    return Point3DCL(p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2]
                    +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]);
}

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


template<class T=double>
class OldQuadBaseCL
{
  protected:
    OldQuadBaseCL( size_t size) : val( size) {}
    OldQuadBaseCL( size_t size, const T& t) : val( t, size) {}

  public:
    std::valarray<T> val;

    size_t size() const
        { return val.size(); }

    // Arithmetik
    template<typename FuncT>
    OldQuadBaseCL& apply ( FuncT fun)
      { for (size_t i=0; i<val.size(); ++i) val[i]=fun(val[i]); return *this; }

    OldQuadBaseCL& operator+= (const OldQuadBaseCL &q)
      { val+=q.val; return *this; }
    OldQuadBaseCL& operator+= (const T &t)
      { val+= t; return *this; }
    OldQuadBaseCL& operator*= (const T &t)
      { val*= t; return *this; }
    template <typename U> OldQuadBaseCL& operator*= (const OldQuadBaseCL<U>& q)
      { for (size_t i=0; i<val.size(); ++i) val[i]*=q.val[i]; return *this; }
//    template <typename U> OldQuadBaseCL& operator*= (const U& c)
//      { val*=c; return *this; }

    OldQuadBaseCL operator+ (const OldQuadBaseCL &q) const
      { return OldQuadBaseCL(*this)+= q; }
    OldQuadBaseCL operator+ (const T &t) const
      { return OldQuadBaseCL(*this)+= t; }
    OldQuadBaseCL operator* (const T &t) const
      { return OldQuadBaseCL(*this)*= t; }
    template <typename U> OldQuadBaseCL operator* (const OldQuadBaseCL<U> &q) const
      { return OldQuadBaseCL(*this)*= q; }
//    template <typename U> OldQuadBaseCL operator* (const U& c) const
//      { OldQuadBaseCL ret(*this); ret.val*=c; return ret; }
    friend OldQuadBaseCL<double> dot (const OldQuadBaseCL<Point3DCL>&, const OldQuadBaseCL<Point3DCL>&);
};

inline OldQuadBaseCL<double> dot (const OldQuadBaseCL<Point3DCL> &q1, const OldQuadBaseCL<Point3DCL> &q2)
{
    OldQuadBaseCL<double> res( q1.size());
    for (size_t i=0; i<q1.size(); ++i)
        res.val[i]= inner_prod( q1.val[i], q2.val[i]);
    return res;
}


template<class T=double>
class OldQuad2CL: public OldQuadBaseCL<T>
{
  public:
    using OldQuadBaseCL<T>::size;
    using OldQuadBaseCL<T>::val;

    static const BaryCoordCL Node[5]; // Stuetzstellen
    static const double      Wght[5]; // Gewichte

    OldQuad2CL() : OldQuadBaseCL<T>( 5) {}
    OldQuad2CL( const T& t) : OldQuadBaseCL<T>( 5, t) {}
    OldQuad2CL( const OldQuadBaseCL<T>& q) : OldQuadBaseCL<T>(q) {}

    // Initialisiere die Knotenwerte
/*    template <ElemT et>
    void set( const ElemBaseCL<et,T>& e)
      { for (size_t i=0; i<size(); ++i) val[i]= e.eval( Node[i]); }
*/
    // Werte Quadraturformel aus
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const
      { T sum= T(); for (size_t i=0; i<size(); ++i) sum+= Wght[i]*val[i]; return sum*absdet; }
    // Folgende Spezialformeln nutzen die spezielle Lage der Stuetzstellen aus
    // zur Annaeherung von \int f*phi,    phi = P1-/P2-Hutfunktion
    T quadP1 (int i, double absdet) const
      { return ((1./120.)*val[i] + (1./30.)*val[4])*absdet; }
    T quadP1 (int i, int j, double absdet) const
      { return (i!=j ? (1./720.)*(val[i]+val[j]) + (1./180.)*val[4]
                     : (1./180.)*val[i]          + (1./90.)*val[4]  )*absdet;}
    T quadP2 (int i, double absdet) const
    {
        return (i<4 ? (1./360.)*val[i] - (1./90.)*val[4]
                    : (1./180.)*(val[VertOfEdge(i-4,0)]+val[VertOfEdge(i-4,1)]) + (1./45.)*val[4]
               )*absdet;
    }

    T quadP2 (int i, int j, double absdet) const
    {
        const double valBary= (i<4 ? -0.125 : 0.25)*(j<4 ? -0.125 : 0.25);
        return ((i!=j || i>=4) ? Wght[4]*val[4]*valBary
                               : Wght[4]*val[4]*valBary + Wght[i]*val[i]
               )*absdet;
    }
};

template<class T>
const BaryCoordCL OldQuad2CL<T>::Node[5]= {
    {1.,0.,0.,0.}, {0.,1.,0.,0.}, {0.,0.,1.,0.}, {0.,0.,0.,1.}, {.25,.25,.25,.25}
    };
template<class T>
const double OldQuad2CL<T>::Wght[5]= { 1./120., 1./120., 1./120., 1./120., 2./15.};


/*
inline std::valarray<double>
dot (const std::valarray<Point3DCL>& a, const std::valarray<Point3DCL>& b)
{
    std::valarray<double> ret;
    for (size_t i=0; i<a.size(); ++i)
        ret[i]= dot( a[i], b[i]);
    return ret;
}
*/

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
typedef OldQuad2CL<double> OldQuadT;
typedef Quad2CL<double> NewQuadT;
typedef Quad2CL<Point3DCL> NewVQuadT;

BndCL theBnd;

double Quadrature( DROPS::MultiGridCL& mg, VecDescCL& vd0, VecDescCL& vd1,
    VecDescCL& /*vd2*/)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nOld Quadrature:\n";
    double ret= 0.;
    P2FuncT f( &vd0, &theBnd, &mg);
    P2FuncT g( &vd1, &theBnd, &mg);
//    P2FuncT h( &vd2, &theBnd, &mg);
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
    VecDescCL& /*vd2*/)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nNew Quadrature:\n";
    double ret= 0.;
    P2FuncT f( &vd0, &theBnd, &mg);
    P2FuncT g( &vd1, &theBnd, &mg);
//    P2FuncT h( &vd2, &theBnd, &mg);
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
    DROPS::IdxDescCL idx( P2_FE);
    DROPS::CreateNumb( 0, idx, mg, Bnd);
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
    DROPS::IdxDescCL vidx( vecP2_FE);
    DROPS::CreateNumb( 0, vidx, mg, Bnd);
    DROPS::VecDescCL vd3( &vidx);
    SetFun( vd3, mg, fv);
    NewVQuadT vq1, vq2;
    std::cout << "dot: "<< dot( vq1, vq2)[1] << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
