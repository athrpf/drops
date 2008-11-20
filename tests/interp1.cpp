#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"

using namespace DROPS;

double f(const SVectorCL<3>& p)
//{ return 42; }
{ return p[0] +10.*p[1] +100.*p[2]+1; }

void MarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte( 0.5);

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It( mg.GetTriangTetraBegin(maxLevel)),
             ItEnd( mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It) {
        if ( (GetBaryCenter( *It)-Mitte).norm() <= std::max( 0.1, 1.5*std::pow( It->GetVolume(), 1.0/3.0)) )
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


void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg)
{
    vd.Data.resize(vd.RowIdx->NumUnknowns());
    P1EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun(&vd, &Bnd, &mg);
    const Uint lvl= vd.GetLevel();
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl), theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit)
    {
        if (!fun.GetBndData()->IsOnDirBnd( *sit)) fun.SetDoF(*sit, f(sit->GetCoord()));
    }
}


int CheckResult(DROPS::P1EvalCL<double, BndCL, DROPS::VecDescCL>& fun)
{
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const DROPS::Uint trilevel= fun.GetLevel();
    std::cout << "Verts:" << std::endl;
    double diff;
    for (MultiGridCL::const_TriangVertexIteratorCL sit= mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= fun.val( *sit) - f( sit->GetCoord());
        std::cout << diff << "\t";
        if (std::abs( diff) > 1e-9) return 1;
    }
    std::cout << std::endl;
    return 0;
}

DROPS::Uint Rule(DROPS::Uint r)
{
    return r < 64 ? r : 127;
}

// Checks every possible tetra-modification.
int TestReMark()
{
    int ttt, ret= 0;
    for (DROPS::Uint i= 0; i<=64; ++i) {
        for (DROPS::Uint j= 0; j<=64; ++j) {
//            std::cout << Rule( i) << "\t-->\t" << Rule( j) << " ";
            DROPS::TetraBuilderCL tet( Rule( i));
            DROPS::MultiGridCL mg( tet);
            DROPS::IdxDescCL i0( DROPS::P1_FE, Bnd), i1( DROPS::P1_FE, Bnd);
            i0.CreateNumbering( mg.GetLastLevel(), mg);
            DROPS::VecDescCL v0, v1;
            v0.SetIdx(&i0);
            SetFun(v0, mg);
            tet.BogoReMark( mg, Rule( j));

            i1.SetTriangLevel (i0.TriangLevel() <= mg.GetLastLevel() ? i0.TriangLevel()
                                                                : mg.GetLastLevel());
            i1.CreateNumbering( i1.TriangLevel(), mg);
            v1.SetIdx( &i1);
            DROPS::P1EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
            DROPS::RepairAfterRefineP1( fun0, v1);
            DROPS::P1EvalCL<double, BndCL, VecDescCL >       fun1( &v1, &Bnd, &mg);
            ttt= CheckResult( fun1);
            ret+= ttt;
            if (ttt != 0)
                std::cout << Rule( i) << "\t-->\t" << Rule( j) << " " << std::endl;
        }
    }
    return ret;
}


int TestRepair()
{
    int ret= 0;
    DROPS::BrickBuilderCL brick( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                                 DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                                 1, 1, 1);
    DROPS::MultiGridCL mg(brick);
    DROPS::IdxDescCL i0( DROPS::P1_FE, Bnd), i1( DROPS::P1_FE, Bnd);
    for (DROPS::Uint i=0; i<8; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg);
        MarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        v1.SetIdx( &i1);
        DROPS::P1EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
        DROPS::RepairAfterRefineP1( fun0, v1);
        DROPS::P1EvalCL<double, BndCL, VecDescCL >       fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1);
    }
    for (DROPS::Uint i=0; i<8; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg);
        UnMarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        if (mg.GetLastLevel() < i0.TriangLevel()) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            i1.SetTriangLevel( i0.TriangLevel()-1);
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
            i1.SetTriangLevel( i0.TriangLevel());
        }
        i1.CreateNumbering( i1.TriangLevel(), mg);
        v1.SetIdx( &i1);
        DROPS::P1EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
        DROPS::RepairAfterRefineP1( fun0, v1);
        DROPS::P1EvalCL<double, BndCL, VecDescCL >       fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1);
    }
    return ret;
}

int TestInterpolateOld()
{
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;
//    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);
    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 3, 3, 3);
    DROPS::MultiGridCL mg(brick);
    MarkDrop(mg, 0);
//    MarkAll( mg);
    mg.Refine();

    DROPS::IdxDescCL i0( DROPS::P1_FE, Bnd), i1( DROPS::P1_FE, Bnd);
    i0.CreateNumbering( 0, mg);
    i1.CreateNumbering( 1, mg);
    VecDescBaseCL<VectorCL> v0, v1;
    v0.SetIdx(&i0);
    v1.SetIdx(&i1);
    SetFun(v0, mg);

    DROPS::P1EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
    DROPS::P1EvalCL<double, BndCL, VecDescCL > fun1( &v1, &Bnd, &mg);
    Interpolate(fun1, fun0);
    return CheckResult( fun1);

}

int main()
{
  try {
    return TestRepair() + TestInterpolateOld() + TestReMark();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
