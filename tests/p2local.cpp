#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/discretize.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "num/bndData.h"

using namespace DROPS;

typedef double (*fun_ptr)(const SVectorCL<3>&);
typedef SVectorCL<3> (*v_inst_fun_ptr)(const SVectorCL<3>&, double);

enum  OutputModeT { SILENT, NOISY };


double f(const SVectorCL<3>& p)
{ return p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2] +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]; }

double g(const SVectorCL<3>& p)
{  return p[0] +10.*p[1] +100.*p[2]+1000.; }

double h(const SVectorCL<3>& p)
{  return sin(M_PI*p[0])*sin(M_PI*p[1])*sin(M_PI*p[2]); }

double h2(const SVectorCL<3>& p)
{  return cos(M_PI*p[0])*cos(M_PI*p[1])*cos(M_PI*p[2]); }

double g2(const DROPS::Point3DCL& p)
{
    return (-1.)*p[0];
}

Point3DCL fv(const SVectorCL<3>& p, double t= 0.0)
{
    return (1.-t)*(p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2]
        +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]) + t*(-2.); }

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
typedef NoBndDataCL<Point3DCL> VBndCL;

BndCL theBnd;
VBndCL theVBnd;

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns);
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &theBnd, &mg);
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

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, v_inst_fun_ptr f, double t)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns);
    P2EvalCL<Point3DCL, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &theBnd, &mg, t);
    const Uint lvl= vd.RowIdx->TriangLevel;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord(), t));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5, t));
    }
}


typedef P2EvalCL<double, BndCL, VecDescCL> P2FuncT;

double therho( const Point3DCL& p) { return norm( 1e-2*p); }
double themu( const Point3DCL& p) { return norm( 1e-2*std_basis<3>( 1) + p); }
double theRe= 1e1;
Point3DCL theg= 9.81*std_basis<3>( 3);
const BaryCoordCL bary( 0.25);
const Point3DCL bary3d( 0.25);

double Quadrature( DROPS::MultiGridCL& mg, VecDescCL& vd0, VecDescCL& /*vd1*/,
    VecDescCL& vd2)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nOld Setup:\n";
//    const IdxT num_unks_vel= vd0.RowIdx->NumUnknowns;
//    const Uint lvl         = vd0.GetLevel();
    const Uint vidx        = vd0.RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    double ret= 0.;
    Quad2CL<double> rho, mu_Re, Phi, kreuzterm;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    P2FuncT ls( &vd2, &theBnd, &mg);
    SMatrixCL<3,3> T;
    double coupA[10][10], coupM[10][10];
    double det, absdet;

    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
    
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            if(!(IsOnDirBnd[i]= theBnd.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            rhs[i]= f( sit->GetVertex(i)->GetCoord());
            Phi[i]= ls.val( *sit->GetVertex(i));
        }
        for (int i=0; i<6; ++i) {
            if (!(IsOnDirBnd[i+4]= theBnd.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        rhs[4]= h( GetBaryCenter( *sit));
        Phi[4]= ls.val( *sit, 0.25, 0.25, 0.25);

        // rho = rho( Phi),    mu_Re= mu( Phi)/Re
        rho=   Phi;     rho.apply( therho);
        mu_Re= Phi;     mu_Re.apply( themu);     mu_Re*= 1./theRe;

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( theg)*rho;
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j) {
                // dot-product of the gradients
                const Quad2CL<double> dotGrad( dot( Grad[i], Grad[j]) * mu_Re);
                coupA[i][j]= coupA[j][i]= dotGrad.quad( absdet);
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j, absdet);
                ret+= coupM[i][j] - coupA[i][j];
            }
    }
    return ret;
}


double NewQuadrature(DROPS::MultiGridCL& mg, VecDescCL& vd0, VecDescCL& /*vd1*/,
    VecDescCL& vd2)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nNew Setup:\n";
//    const IdxT num_unks_vel= vd0.RowIdx->NumUnknowns;
//    const Uint lvl         = vd0.GetLevel();
      const Uint vidx        = vd0.RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    double ret= 0.;
    Quad2CL<double> rho, mu_Re, Phi, kreuzterm;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    SMatrixCL<3,3> T;
    double coupA[10][10], coupM[10][10];
    double det, absdet;
    LocalP2CL<> ls;
    
    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        ls.assign( *sit, vd2, theBnd);
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            if(!(IsOnDirBnd[i]= theBnd.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            rhs[i]= f( sit->GetVertex(i)->GetCoord());
            Phi[i]= ls[i];
        }
        for (int i=0; i<6; ++i) {
            if (!(IsOnDirBnd[i+4]= theBnd.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        rhs[4]= h( GetBaryCenter( *sit));
        Phi[4]= ls( bary);

        // rho = rho( Phi),    mu_Re= mu( Phi)/Re
        rho=   Phi;     rho.apply( therho);
        mu_Re= Phi;     mu_Re.apply( themu);     mu_Re*= 1./theRe;

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( theg)*rho;
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j) {
                // dot-product of the gradients
                const Quad2CL<double> dotGrad( dot( Grad[i], Grad[j]) * mu_Re);
                coupA[i][j]= coupA[j][i]= dotGrad.quad( absdet);
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j, absdet);
                ret+= coupM[i][j] - coupA[i][j];
            }
    }
    return ret;

}

double Tests(DROPS::MultiGridCL& mg, VecDescCL& vd0, VecDescCL& vd1)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTests:\n";
    Point3DCL ret( 0.);
    Point3DCL c( 0.2);

    P2FuncT lsfun( &vd1, &theBnd, &mg);
    LocalP2CL<> ls2( *mg.GetTriangTetraBegin(), lsfun);
    ret+= ls2[4];

    LocalP2CL<Point3DCL> f1;
    LocalP2CL<Point3DCL> f2( *mg.GetTriangTetraBegin(), fv, 0.3);

    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        f1.assign( *sit, vd0, theVBnd, 0.5);
        ret+= (c*f1 + f1)[4] + fv( sit->GetVertex( 0)->GetCoord(), 0.1);
        ret-= f1( bary)[2];
    }
    
    return ret[0] + ret[1] + ret[2] + f2( bary)[2];    

}

int main ()
{
  try {
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
				DROPS::std_basis<3>(2),
				DROPS::std_basis<3>(3),
				10, 20, 30);
    DROPS::MultiGridCL mg( brick);
    DROPS::IdxDescCL idx( 1,1,0,0);
    idx.TriangLevel= 0;
    idx.NumUnknowns= 0;
    DROPS::CreateNumbOnVertex( idx.GetIdx(), idx.NumUnknowns, 1,
        mg.GetTriangVertexBegin( idx.TriangLevel),
        mg.GetTriangVertexEnd( idx.TriangLevel),
        theBnd);
    DROPS::CreateNumbOnEdge( idx.GetIdx(), idx.NumUnknowns, 1,
        mg.GetTriangEdgeBegin( idx.TriangLevel),
        mg.GetTriangEdgeEnd( idx.TriangLevel),
        theBnd);
    DROPS::IdxDescCL vidx( 3,3,0,0);
    vidx.TriangLevel= 0;
    vidx.NumUnknowns= 0;
    DROPS::CreateNumbOnVertex( vidx.GetIdx(), vidx.NumUnknowns, 3,
        mg.GetTriangVertexBegin( vidx.TriangLevel),
        mg.GetTriangVertexEnd( vidx.TriangLevel),
        theVBnd);
    DROPS::CreateNumbOnEdge( vidx.GetIdx(), vidx.NumUnknowns, 3,
        mg.GetTriangEdgeBegin( vidx.TriangLevel),
        mg.GetTriangEdgeEnd( vidx.TriangLevel),
        theVBnd);
    DROPS::VecDescCL vd0( &idx);
    SetFun( vd0, mg, f);
    DROPS::VecDescCL vd1( &idx);
    SetFun( vd1, mg, g);
    DROPS::VecDescCL vd2( &idx);
    SetFun( vd2, mg, h2);
    DROPS::VecDescCL vd3( &vidx);
    SetFun( vd3, mg, fv, 0.5);
    DROPS::VecDescCL vd4( &vidx);
    SetFun( vd4, mg, fv, 0.2);
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
    time.Start();
    double q2= Tests( mg, vd3, vd2);
    time.Stop();
    std::cout << "Tests: " << q2 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();

    MarkAll( mg);
    mg.Refine();
    DROPS::IdxDescCL idx1( 1,1,0,0);
    idx1.TriangLevel= 1;
    idx1.NumUnknowns= 0;
    DROPS::CreateNumbOnVertex( idx1.GetIdx(), idx1.NumUnknowns, 1,
        mg.GetTriangVertexBegin( idx1.TriangLevel),
        mg.GetTriangVertexEnd( idx1.TriangLevel),
        theBnd);
    DROPS::CreateNumbOnEdge( idx1.GetIdx(), idx1.NumUnknowns, 1,
        mg.GetTriangEdgeBegin( idx1.TriangLevel),
        mg.GetTriangEdgeEnd( idx1.TriangLevel),
        theBnd);
    DROPS::VecDescCL vd5( &idx1);
    SetFun( vd5, mg, f);
    LocalP2CL<> restr1( *mg.GetTriangTetraBegin( 0), vd5, theBnd, 0.3);
    double maxdiff= 0.0, maxdiff2= 0.0;
    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin( 0),
         end=mg.GetTriangTetraEnd( 0); sit != end; ++sit) {
        TetraCL& tet= *sit;
        restr1.assign( tet, vd5, theBnd, 0.3);
        maxdiff= std::max( maxdiff, std::abs( restr1( bary) - f( tet.GetVertex( 0)->GetCoord()*0.25
            + tet.GetVertex( 1)->GetCoord()*0.25 + tet.GetVertex( 2)->GetCoord()*0.25
            + tet.GetVertex( 3)->GetCoord()*0.25)));
        for (int i=0; i<4; ++i) {
            maxdiff2= std::max( maxdiff2, std::abs( restr1( std_basis<4>( i+1))
                - f( tet.GetVertex( i)->GetCoord()))); 
        }
        for (int i=0; i<6; ++i) {
            maxdiff2= std::max( maxdiff2, std::abs( restr1( 0.5*(std_basis<4>( VertOfEdge( i, 0)+1) 
                + std_basis<4>( VertOfEdge( i, 1)+1)))
                - f( GetBaryCenter( *tet.GetEdge( i))))); 
        }
        
    }
    std::cout << "max. Differenz im Schwerpunkt: " << maxdiff
        << "\tin den Freiheitsgraden: " << maxdiff2 << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
