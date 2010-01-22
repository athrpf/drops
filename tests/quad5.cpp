/// \file quad5.cpp
/// \brief tests 3D quadrature of order 5
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
#include <iomanip>

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

double hsq (const Point3DCL& p, double)
{
    return pow(p.norm_sq(),1); // p[0];
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

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &Bnd, &mg);
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

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, vfun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    P2EvalCL<Point3DCL, VBndCL,VecDescBaseCL<VectorCL> > fun( &vd, &VBnd, &mg);
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

typedef P2EvalCL<double, BndCL, VecDescCL> P2FuncT;
typedef Quad5CL<double> NewQuadT;
typedef Quad5CL<Point3DCL> NewVQuadT;

BndCL theBnd;


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

int degx, degy, degz;

double f(const Point3DCL& p, double)
{
    return  (degx==0 ? 1. : std::pow( p[0], degx))
           *(degy==0 ? 1. : std::pow( p[1], degy))
           *(degz==0 ? 1. : std::pow( p[2], degz));
}

double exactint[56] = { // with maple
  0.1666666666666666666666667, 0.04166666666666666666666667,
  0.01666666666666666666666667, 0.008333333333333333333333333,
  0.004761904761904761904761905, 0.002976190476190476190476190,
  0.04166666666666666666666667, 0.008333333333333333333333333,
  0.002777777777777777777777778, 0.001190476190476190476190476,
  0.0005952380952380952380952381, 0.01666666666666666666666667,
  0.002777777777777777777777778, 0.0007936507936507936507936508,
  0.0002976190476190476190476190, 0.008333333333333333333333333,
  0.001190476190476190476190476, 0.0002976190476190476190476190,
  0.004761904761904761904761905, 0.0005952380952380952380952381,
  0.002976190476190476190476190, 0.04166666666666666666666667,
  0.008333333333333333333333333, 0.002777777777777777777777778,
  0.001190476190476190476190476, 0.0005952380952380952380952381,
  0.008333333333333333333333333, 0.001388888888888888888888889,
  0.0003968253968253968253968254, 0.0001488095238095238095238095,
  0.002777777777777777777777778, 0.0003968253968253968253968254,
  0.00009920634920634920634920635, 0.001190476190476190476190476,
  0.0001488095238095238095238095, 0.0005952380952380952380952381,
  0.01666666666666666666666667, 0.002777777777777777777777778,
  0.0007936507936507936507936508, 0.0002976190476190476190476190,
  0.002777777777777777777777778, 0.0003968253968253968253968254,
  0.00009920634920634920634920635, 0.0007936507936507936507936508,
  0.00009920634920634920634920635, 0.0002976190476190476190476190,
  0.008333333333333333333333333, 0.001190476190476190476190476,
  0.0002976190476190476190476190, 0.001190476190476190476190476,
  0.0001488095238095238095238095, 0.0002976190476190476190476190,
  0.004761904761904761904761905, 0.0005952380952380952380952381,
  0.0005952380952380952380952381, 0.002976190476190476190476190
};


void TestExactness()
{
    DROPS::TetraBuilderCL tet( 0);
    DROPS::MultiGridCL mg( tet);
    TetraCL& s= *mg.GetAllTetraBegin();
//    s.DebugInfo( std::cout);
    std::cout.precision( 18);

    NewQuadT q;
    size_t c= 0;
    for (degz= 0; degz <= 5; ++degz) {
        for (degy= 0; degy + degz <= 5; ++degy) {
            for (degx= 0; degx + degy + degz <= 5; ++degx) {
                q.assign( s, f);
                std::cout << "degz: " << degz << "\tdegy: " << degy << "\tdegx: " << degx
                          << "\t\tI-Q_h: " << exactint[c++] - q.quad( 1.)
                          << "\tIntegral: " << q.quad( 1.) <<  "             ";
                for (size_t i= 0; i < q.size(); ++i)
                    std::cout << '\t' << q[i];
                std::cout << std::endl;
            }
        }
    }
}

void TestP2Mass()
{
    std::cout.precision( 18);
    for (size_t i= 0; i < 10; ++i) {
        for (size_t j= 0; j < 10; ++j) {
            NewQuadT q, q2;
            double m1= P2DiscCL::GetMass( i, j)*1.0;
            LocalP2CL<> phi;
            phi[i]= 1.0;
            q.assign( phi);
            double m2= q.quadP2( j, 1.);
            LocalP2CL<> phi2;
            phi2[j]= 1.0;
            q2.assign( phi2);
            q*= q2;
            double m3= q.quad( 1.);
            std::cout << "i: " << i << "\tj: " << j
                      << "\tm1: " << m1
                      << "\tm2: " << m2
                      << "\tm3: " << m3
                      << "\tm1-m2: " << m1 - m2
                      << "\tm1-m3: " << m1 - m3
                      << "\tm2-m3: " << m2 - m3
                      << '\n';
        }
    }
}


void TestTransform()
{
    std::cout << "\n\nTestTransform():\n";
    Uint rule = 0; // Choose a refinement rule
    double absdet=0.;
    double intval=0.;
    Quad5CL<> qf;
    TetraBuilderCL builder(rule);
    MultiGridCL mg( builder);

    DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), sit) {
        absdet= sit->GetVolume()*6.0;
        qf.assign( *sit, &hsq);
        intval+= qf.quad( absdet);
        }
    std::cout << std::setprecision(10) <<  "Integral: " << intval << '\n';

    SArrayCL<BaryCoordCL,4> M;
    for (Uint i=0; i<4; ++i)
        M[i][i]=2.;
    M[0][0]=1.; M[1][0]= M[2][0]= M[3][0]= -1.;

/*    M[1][1]=1.;
    M[1][2]=2.;
    M[2][2]=3.;
    M[3][1]=1.;
    M[3][2]=1.;
    M[3][3]=1.;*/
    BaryCoordCL* nodes;
    nodes= qf.TransformNodes( M);
    intval=0.;
    DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), sit) {
        absdet= sit->GetVolume()*6.0*8.0;

        qf.assign( *sit, &hsq, 0.0, nodes);
        intval+= qf.quad( absdet);
    }
    std::cout << "Integral: " << intval << '\n';
}


int main ()
{
  try {
    TestExactness();
    TestP2Mass();

    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                40, 40, 40);
    DROPS::MultiGridCL mg( brick);
    DROPS::IdxDescCL idx( P2_FE, Bnd);
    idx.CreateNumbering( 0, mg);
    DROPS::VecDescCL vd0( &idx);
    SetFun( vd0, mg, f);
    DROPS::VecDescCL vd1( &idx);
    SetFun( vd1, mg, g);
    DROPS::VecDescCL vd2( &idx);
    SetFun( vd2, mg, h);
    TimerCL time;
    time.Start();
    double q1= NewQuadrature( mg, vd0, vd1, vd2);
    time.Stop();
    std::cout << "integral: " << q1 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    DROPS::IdxDescCL vidx( vecP2_FE, VBnd);
    vidx.CreateNumbering( 0, mg);
    DROPS::VecDescCL vd3( &vidx);
    SetFun( vd3, mg, fv);
    NewVQuadT vq1, vq2;
    std::cout << "dot: "<< dot( vq1, vq2)[1] << std::endl;

    TestTransform();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
