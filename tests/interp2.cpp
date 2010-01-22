/// \file interp2.cpp
/// \brief tests implementation of RepareAfterRefineP2
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
#include "num/fe.h"
#include "misc/problem.h"

using namespace DROPS;

typedef double (*fun_ptr)(const SVectorCL<3>&);

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

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &Bnd, &mg);
    const Uint lvl= vd.GetLevel();
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord()));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5));
    }
}


int CheckResult(DROPS::P2EvalCL<double, BndCL,
                const DROPS::VecDescCL>& fun, fun_ptr f, OutputModeT om)
{
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const DROPS::Uint trilevel= fun.GetLevel();
    if (om!=SILENT) std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0.;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= fun.val( *sit) - f( sit->GetCoord());
        if ( std::abs(diff) > vmaxdiff) { ++ret; vmaxdiff= std::abs(diff); v= &*sit; }
    }
    if (om!=SILENT) std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin( trilevel),
         theend= mg.GetTriangEdgeEnd( trilevel); sit!=theend; ++sit) {
        diff = fun.val( *sit, .5) - f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5);
        if ( std::abs(diff) > emaxdiff) { ++ret; emaxdiff= std::abs(diff); e= &*sit; }
    }
    if (om!=SILENT) {
        std::cout << "maximale Differenz Vertices: " << vmaxdiff << " auf\n";
        if (v) v->DebugInfo( std::cout);
        std::cout << "maximale Differenz Edges: " << emaxdiff << " auf\n";
        if (e) e->DebugInfo( std::cout);
        std::cout << std::endl;
    }
    return ret;
}


DROPS::Uint Rule(DROPS::Uint r)
{
    return r < 64 ? r : 127;
}

// True, iff every bit in p is also set in q.
bool SubSuperPattern(DROPS::Uint p, DROPS::Uint q)
{
    return !((~q) & p);
}


// Checks every possible tetra-modification.
int TestReMark()
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair on single tetra-combinations:\n";
    int ttt, ret= 0;
    for (DROPS::Uint i= 0; i<=64; ++i) {
        for (DROPS::Uint j= 0; j<=64; ++j) {
            DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
//            std::cout << Rule( i) << "\t-->\t" << Rule( j) << " ";
            DROPS::TetraBuilderCL tet( Rule( i), DROPS::std_basis<3>( 1),
                                                 DROPS::std_basis<3>( 2),
                                                 DROPS::std_basis<3>( 3),
                                                 DROPS::Point3DCL( 1.0));
            DROPS::MultiGridCL mg( tet);
            DROPS::IdxDescCL i0(P2_FE, Bnd), i1(P2_FE, Bnd);
            i0.CreateNumbering( mg.GetLastLevel(), mg);
            DROPS::VecDescCL v0, v1;
            v0.SetIdx( &i0);
            SetFun( v0, mg, f);
//            SetFun( v0, mg, g2);
            tet.BogoReMark( mg, Rule( j));

            const Uint i1_Level= i0.TriangLevel() <= mg.GetLastLevel() ? i0.TriangLevel()
                                                                       : mg.GetLastLevel();
            i1.CreateNumbering( i1_Level, mg);
            v1.SetIdx( &i1);
            DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
            DROPS::RepairAfterRefineP2( fun0, v1);
            DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
            ttt= CheckResult( fun1, f, SILENT);
//            ttt= CheckResult( fun1, g2, NOISY);
            ret+= ttt;
            if (ttt != 0 && SubSuperPattern( Rule( i) & 63, Rule( j) & 63))
                std::cout << "Aerger: " << Rule( i) << "\t-->\t" << Rule( j) << " " << std::endl;
        }
    }
    return ret;
}


int TestRepairUniform()
{
    int ret= 0;
    DROPS::BrickBuilderCL brick( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                                 DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                                 1, 1, 1);
    DROPS::MultiGridCL mg(brick);
    DROPS::IdxDescCL i0( P2_FE, Bnd), i1( P2_FE, Bnd);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for uniform refinement with quadratic function:\n";
    for (DROPS::Uint i=0; i<5; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx( &i0);
        SetFun( v0, mg, f);
        MarkAll( mg);
        mg.Refine();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        v1.SetIdx( &i1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
        DROPS::RepairAfterRefineP2( fun0, v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY);
        i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
    }

    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for uniform coarsening with linear function:\n";
    for (DROPS::Uint i=0; i<5; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg, g);
        UnMarkAll( mg);
        mg.Refine();
        Uint i1_Level= i0.TriangLevel();
        if (mg.GetLastLevel() < i1_Level) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            --i1_Level;
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
        }
        i1.CreateNumbering( i1_Level, mg);
        v1.SetIdx( &i1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
        DROPS::RepairAfterRefineP2( fun0, v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, g, NOISY);
        if (mg.GetLastLevel() < i0.TriangLevel()) {
            Uint level= mg.GetLastLevel();
            DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllVertexBegin( level),
                                        mg.GetAllVertexEnd( level));
            DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllEdgeBegin( level),
                                        mg.GetAllEdgeEnd( level));
        }
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllVertexBegin( i1.TriangLevel()),
                                    mg.GetAllVertexEnd( i1.TriangLevel()));
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllEdgeBegin( i1.TriangLevel()),
                                    mg.GetAllEdgeEnd( i1.TriangLevel()));
    }
    return ret;
}

int TestRepair()
{
    int ret= 0;
    DROPS::BrickBuilderCL brick( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                                 DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                                 1, 1, 1);
    DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
    DROPS::MultiGridCL mg(brick);
    DROPS::IdxDescCL i0( P2_FE, Bnd), i1( P2_FE, Bnd);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for drop refinement with linear function:\n";
    for (DROPS::Uint i=0; i<8; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx( &i0);
        SetFun( v0, mg, g);
        MarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        v1.SetIdx( &i1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
        DROPS::RepairAfterRefineP2( fun0, v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, g, NOISY);
        i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
    }

    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for drop coarsening with linear function:\n";
    for (DROPS::Uint i=0; i<8; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg, g);
        UnMarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        Uint i1_Level= i0.TriangLevel();
        if (mg.GetLastLevel() < i1_Level) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            --i1_Level;
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
        }
        i1.CreateNumbering( i1_Level, mg);
        v1.SetIdx( &i1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
        DROPS::RepairAfterRefineP2( fun0, v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, g, NOISY);
        if (mg.GetLastLevel() < i0.TriangLevel())
            i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
    }
    return ret;
}


int TestInterpolateOld()
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting interpolation:\n";
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 4, 4, 4);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 3, 3, 3);
    DROPS::MultiGridCL mg(brick);
    MarkDrop(mg, 0);
//    MarkAll( mg);
    mg.Refine();
    MarkDrop(mg, 1);
    mg.Refine();

    DROPS::IdxDescCL i0( DROPS::P2_FE, Bnd), i1( DROPS::P2_FE, Bnd);
    i0.CreateNumbering( 0, mg);
    i1.CreateNumbering( 1, mg);

    VecDescBaseCL<VectorCL> v0, v1;
    v0.SetIdx( &i0);
    v1.SetIdx( &i1);
    SetFun( v0, mg, f);

    v1.Data.resize( v1.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL, const VecDescBaseCL<VectorCL> > fun0( &v0, &Bnd, &mg);
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun1( &v1, &Bnd, &mg);
    Interpolate( fun1, fun0);
//    DROPS::RepairAfterRefineP2( fun0, v1);
    std::cout << "Verts:" << std::endl;
    double diff;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(1),
         theend= mg.GetTriangVertexEnd(1); sit!=theend; ++sit) {
        diff= fun1.val(*sit) - f(sit->GetCoord());
        std::cout << diff << "\t";
        if (diff!=0.) return 1;
    }
    std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(1),
         theend= mg.GetTriangEdgeEnd(1); sit!=theend; ++sit) {
        diff = fun1.val( *sit, .5) - f( (sit->GetVertex(0)->GetCoord()+sit->GetVertex(1)->GetCoord())*.5);
        std::cout << diff << "\t";
        if (diff!=0.) return 1;
    }
    std::cout << std::endl;
    std::cout << std::endl << DROPS::SanityMGOutCL(mg) << std::endl;
    return 0;
}


int main ()
{
  try {
    int ret= TestRepairUniform();
    ret+= TestRepair();
    ret+= TestInterpolateOld();
    return ret + TestReMark();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
