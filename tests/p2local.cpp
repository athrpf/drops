/// \file p2local.cpp
/// \brief tests implementation of LocalP2CL
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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
{  return std::sin(M_PI*p[0])*std::sin(M_PI*p[1])*std::sin(M_PI*p[2]); }

double h2(const SVectorCL<3>& p)
{  return std::cos(M_PI*p[0])*std::cos(M_PI*p[1])*std::cos(M_PI*p[2]); }

double g2(const DROPS::Point3DCL& p)
{
    return (-1.)*p[0];
}

Point3DCL fv(const SVectorCL<3>& p, double t= 0.0)
{
    return Point3DCL((1.-t)*(p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2]
                     +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2])
                     + t*(-2.));
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
typedef NoBndDataCL<Point3DCL> VBndCL;

BndCL theBnd;
VBndCL theVBnd;

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &theBnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel();
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
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    vd.t = t;
    P2EvalCL<Point3DCL, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &theBnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel();
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

double therho( double x) { return std::fabs(x)+0.1; }
double themu( double x) { return x*x; }
Point3DCL theg= 9.81*std_basis<3>( 3);
const BaryCoordCL bary( 0.25);
const Point3DCL bary3d( 0.25);


double Quadrature( DROPS::MultiGridCL& mg, VecDescCL& /*vd0*/, VecDescCL& /*vd1*/,
    VecDescCL& vd2)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nOld Setup:\n";
//    const IdxT num_unks_vel= vd0.RowIdx->NumUnknowns;
//    const Uint lvl         = vd0.GetLevel();

    double ret= 0.;
    Quad2CL<double> rho, mu, Phi, kreuzterm;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    P2FuncT ls( &vd2, &theBnd, &mg);
    SMatrixCL<3,3> T;
    double det, absdet;

    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            rhs[i]= Point3DCL (h( sit->GetVertex(i)->GetCoord()));
            Phi[i]= ls.val( *sit->GetVertex(i));
        }
        rhs[4]= Point3DCL (h( GetBaryCenter( *sit)));
        Phi[4]= ls.val( *sit, 0.25, 0.25, 0.25);

        // rho = rho( Phi),    mu= mu( Phi)
        rho= Phi;     rho.apply( therho);
        mu=  Phi;     mu.apply( themu);

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( theg)*rho;
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j) {
                // dot-product of the gradients
                const double cA= Quad2CL<>(dot( Grad[i], Grad[j]) * mu).quad( absdet);
                const double cM= rho.quadP2(i,j, absdet);
                ret+= cM - cA;
            }
    }
    return ret;
}


double NewQuadrature(DROPS::MultiGridCL& mg, VecDescCL& /*vd0*/, VecDescCL& /*vd1*/,
    VecDescCL& vd2)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nNew Setup:\n";
//    const IdxT num_unks_vel= vd0.RowIdx->NumUnknowns;
//    const Uint lvl         = vd0.GetLevel();

    double ret= 0.;
    Quad2CL<double> rho, mu, Phi, kreuzterm;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls;

    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        ls.assign( *sit, vd2, theBnd);
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            rhs[i]= Point3DCL (h( sit->GetVertex(i)->GetCoord()));
            Phi[i]= ls[i];
        }
        rhs[4]= Point3DCL (h( GetBaryCenter( *sit)));
        Phi[4]= ls( bary);

        // rho = rho( Phi),    mu= mu( Phi)
        rho= Phi;     rho.apply( therho);
        mu=  Phi;     mu.apply( themu);

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( theg)*rho;
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j) {
                // dot-product of the gradients
                const double cA= Quad2CL<>(dot( Grad[i], Grad[j]) * mu).quad( absdet);
                const double cM= rho.quadP2(i,j, absdet);
                ret+= cM - cA;
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
    ret+= Point3DCL( ls2[4]);

    LocalP2CL<Point3DCL> f1;
    LocalP2CL<Point3DCL> f2( *mg.GetTriangTetraBegin(), fv, 0.3);

    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        vd0.t = 0.5;
        f1.assign( *sit, vd0, theVBnd);
        ret+= (c*f1 + f1)[4] + fv( sit->GetVertex( 0)->GetCoord(), 0.1);
        ret-= Point3DCL( f1( bary)[2]);
    }

    return ret[0] + ret[1] + ret[2] + f2( bary)[2];

}


class MulCL
{
  public:
    int i;
    MulCL () : i( 1) {}
    double operator() (double x) { return i++*x;}
};

void MemberApplyTest()
{
    LocalP2CL<> l;
    l= 1.0;
    MulCL m;
    l.apply( m);
    std::cout << l[0] << '\n' << l[1] << '\n' << l[2] << '\n'
        << l[3] << '\n' << l[4] << '\n' << l[5] << '\n'
        << l[6] << '\n' << l[7] << '\n' << l[8] << '\n' << l[9] << std::endl;
}


int main ()
{
  try {
    MemberApplyTest();

    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                10, 20, 30);
    DROPS::MultiGridCL mg( brick);
    DROPS::IdxDescCL idx( P2_FE, theBnd);
    idx.CreateNumbering( 0, mg);
    DROPS::IdxDescCL vidx( vecP2_FE, theVBnd);
    vidx.CreateNumbering( 0, mg);
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
    DROPS::IdxDescCL idx1( P2_FE, theBnd);
    idx1.CreateNumbering( 1, mg);
    DROPS::VecDescCL vd5( &idx1);
    SetFun( vd5, mg, f);
    vd5.t = 0.3;
    LocalP2CL<> restr1( *mg.GetTriangTetraBegin( 0), vd5, theBnd);
    double maxdiff= 0.0, maxdiff2= 0.0;
    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin( 0),
         end=mg.GetTriangTetraEnd( 0); sit != end; ++sit) {
        TetraCL& tet= *sit;
        vd5.t = 0.3;
        restr1.assign( tet, vd5, theBnd);
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
