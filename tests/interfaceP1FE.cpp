/// \file interfaceP1FE.cpp
/// \brief tests implementation of the numbering methods of interface finite elements
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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

#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/surfacetension.h"
#include <fstream>

using namespace DROPS;

void
TestSingleTetra()
{
    TetraBuilderCL tetra( 0); // unrefined reference tetra
    MultiGridCL mg( tetra);
    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    BndCondT bc[4]= { NoBC, NoBC, NoBC, NoBC };
    LsetBndDataCL::bnd_val_fun bfun[4]= { 0,0,0,0};
    LsetBndDataCL lsbnd( 4, bc, bfun);
    LevelsetP2CL lset( mg, lsbnd, sf);

    lset.idx.CreateNumbering( 0, mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Phi.Data= 1.0;

    IdxDescCL ifaceidx( P1IF_FE);
    std::cout << "Testing vertex numbering around no interface:" << std::endl;
    ifaceidx.CreateNumbering( 0, mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);

    std::cout << "Testing vertex numbering interface in 1 tetra:" << std::endl;
    lset.Phi.Data[0]= -1.0;
    ifaceidx.CreateNumbering( 0, mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);
}


double plane(const Point3DCL& p)
{
    return p[0] - 0.41;
}

double plane2(const Point3DCL& p)
{
    return p[0] - 0.4;
}

double sphere_1(const Point3DCL& p)
{
    Point3DCL c( 0.47);
    Point3DCL x= fabs( p - c);

    return x[0] + 1.11*x[1] + 1.111*x[2] - 0.31;
}

double sphere_2(const Point3DCL& p)
{
    Point3DCL c( 0.5);
    Point3DCL x= p - c;

    return x.norm() - 0.31;
}

double sphere_max(const Point3DCL& p)
{
    Point3DCL c( 0.47);
    Point3DCL x= fabs( p - c);

    return std::max( x[0], std::max( x[1], x[2])) - 0.31;
}

void
TestPlaneInCube()
{
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                10, 1, 1);
    MultiGridCL mg( brick);
    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    BndCondT bc[6]= { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
    LsetBndDataCL::bnd_val_fun bfun[6]= { 0,0,0,0,0,0};
    LsetBndDataCL lsbnd( 6, bc, bfun);
    LevelsetP2CL lset( mg, lsbnd, sf);

    lset.idx.CreateNumbering( 0, mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Phi.Data= 1.0;

    IdxDescCL ifaceidx( P1IF_FE);
    std::cout << "Testing vertex numbering around planar interface:" << std::endl;
    lset.Init( plane);
    ifaceidx.CreateNumbering( 0, mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);

    std::cout << "Testing vertex numbering around planar interface containing vertices:" << std::endl;
    lset.Init( plane2);
    ifaceidx.CreateNumbering( 0, mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);
}


void SetupInterfaceMassP1OnTriangle (const LocalP1CL<> p1[4],
    Quad5_2DCL<> q[4], MatrixBuilderCL& M, const IdxT Numb[4],
    const BaryCoordCL triangle[3], double det)
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);

    Quad5_2DCL<> m;
    double tmp;
    for (int i= 0; i < 4; ++i) {
        m= q[i]*q[i];
        M( Numb[i], Numb[i])+= m.quad( det);
        for(int j= 0; j < i; ++j) {
            m= (q[j]*q[i]);
            tmp= m.quad( det);
            M( Numb[i], Numb[j])+= tmp;
            M( Numb[j], Numb[i])+= tmp;
        }
    }
}

void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsbnd)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns();
    MatrixBuilderCL M( &matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();
    IdxT Numb[4];

    double det;
    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4];

    InterfaceTriangleCL triangle;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, it) {
    	triangle.Init( *it, ls, lsbnd);

        GetLocalNumbP1NoBnd( Numb, *it, *matM->RowIdx);
        for (int ch= 0; ch < 8; ++ch) {
            if (!triangle.ComputeForChild( ch)) // no patch for this child
                continue;

            det= triangle.GetAbsDet();
            SetupInterfaceMassP1OnTriangle( p1, q, M, Numb, &triangle.GetBary( 0), det);
            if (triangle.IsQuadrilateral()) {
                det*= triangle.GetAreaFrac();
                SetupInterfaceMassP1OnTriangle( p1, q, M, Numb, &triangle.GetBary( 1), det);
            }
        }
    }
    M.Build();
}

void SetupInterfaceMassP1OnTriangleNew (const LocalP2CL<> q[4],
    MatrixBuilderCL& M, const IdxT Numb[4],
    const InterfaceTriangleCL& triangle, Uint tri)
{
    LocalP2CL<> m;

    double tmp;
    for (int i= 0; i < 4; ++i) {
        m= q[i]*q[i];
        M( Numb[i], Numb[i])+= triangle.quad2D( m, tri);
        for(int j= 0; j < i; ++j) {
            m= (q[j]*q[i]);
            tmp= triangle.quad2D( m, tri);
            M( Numb[i], Numb[j])+= tmp;
            M( Numb[j], Numb[i])+= tmp;
        }
    }
}

void SetupInterfaceMassP1New (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsbnd)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns();
    MatrixBuilderCL M( &matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();
    IdxT Numb[4];

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    LocalP2CL<double> q[4];
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i]);

    InterfaceTriangleCL triangle;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, it) {
    	triangle.Init( *it, ls, lsbnd);

        GetLocalNumbP1NoBnd( Numb, *it, *matM->RowIdx);
        for (int ch= 0; ch < 8; ++ch) {
            if (!triangle.ComputeForChild( ch)) // no patch for this child
                continue;

            for (int tri=0; tri<triangle.GetNumTriangles(); ++tri)
                SetupInterfaceMassP1OnTriangleNew( q, M, Numb, triangle, tri);
        }
    }
    M.Build();
}


int main ()
{
  try {
    TestSingleTetra();
    TestPlaneInCube();

    std::cout << "SetupInterfaceMassP1:\n";
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                20, 20, 20);
    MultiGridCL mg( brick);
    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    BndCondT bc[6]= { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
    LsetBndDataCL::bnd_val_fun bfun[6]= { 0,0,0,0,0,0};
    LsetBndDataCL lsbnd( 6, bc, bfun);
    LevelsetP2CL lset( mg, lsbnd, sf);

    lset.idx.CreateNumbering( 0, mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( &sphere_2);

    IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.CreateNumbering( 0, mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    MatDescCL M( &ifaceidx, &ifaceidx);
    SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
    std::cout << "Writing matrix to m_iface.txt\n";
    std::ofstream fff( "m_iface.txt");
    fff.precision( 18);
    fff << M.Data << std::endl;
    fff.close();

    MatDescCL Mnew( &ifaceidx, &ifaceidx);
    SetupInterfaceMassP1New( mg, &Mnew, lset.Phi, lset.GetBndData());
    std::cout << "Writing new matrix to mnew_iface.txt\n";
    std::ofstream ggg( "mnew_iface.txt");
    ggg.precision( 18);
    ggg << Mnew.Data << std::endl;
    ggg.close();

    MatrixCL D;
    D.LinComb( 1., M.Data,  -1., Mnew.Data); // M - Mnew
    std::cout <<   "|| M-Mnew ||_sup  / || M ||_sup  = " << supnorm(D)/supnorm(M.Data)
              << "\n|| M-Mnew ||_frob / || M ||_frob = " << frobeniusnorm(D)/frobeniusnorm(M.Data) << std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
