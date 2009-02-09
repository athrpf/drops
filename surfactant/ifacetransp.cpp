/// \file
/// \brief Discretization for PDEs on an interface.

#include "surfactant/ifacetransp.h"
#include "levelset/levelset.h"
#include "num/spmat.h"
#include <cstring>
#include <cmath>


namespace DROPS {

void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext)
{
    const Uint xidx( x.RowIdx->GetIdx()),
        xextidx( xext.RowIdx->GetIdx()),
        lvl( x.RowIdx->TriangLevel());
    xext.Data= 0.;

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            xext.Data[it->Unknowns( xextidx)]= x.Data[it->Unknowns( xidx)];
    }
}

void Restrict (const MultiGridCL& mg, const VecDescCL& xext, VecDescCL& x)
{
    const Uint xidx( x.RowIdx->GetIdx()),
        xextidx( xext.RowIdx->GetIdx()),
        lvl( x.RowIdx->TriangLevel());

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            x.Data[it->Unknowns( xidx)]= xext.Data[it->Unknowns( xextidx)];
    }
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
        if (Numb[i] == NoIdx) continue;
        m= q[i]*q[i];
        M( Numb[i], Numb[i])+= m.quad( det);
        for(int j= 0; j < i; ++j) {
            if (Numb[j] == NoIdx) continue;
            m= (q[j]*q[i]);
            tmp= m.quad( det);
            M( Numb[i], Numb[j])+= tmp;
            M( Numb[j], Numb[i])+= tmp;
        }
    }
}

void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls)
{
    const IdxT num_unks=  matM->RowIdx->NumUnknowns();
    MatrixBuilderCL M( &matM->Data, num_unks,  num_unks);

    const Uint lvl= matM->GetRowLevel();
    IdxT Numb[4];

    double det;
    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;

    InterfacePatchCL patch;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, it) {
        patch.Init( *it, ls);

        GetLocalNumbP1NoBnd( Numb, *it, *matM->RowIdx);
        for (int ch= 0; ch < 8; ++ch) {
            if (!patch.ComputeForChild( ch)) // no patch for this child
                continue;

            det= patch.GetFuncDet();
            SetupInterfaceMassP1OnTriangle( p1, q, M, Numb, &patch.GetBary( 0), det);
            if (patch.IsQuadrilateral()) {
                det*= patch.GetAreaFrac();
                SetupInterfaceMassP1OnTriangle( p1, q, M, Numb, &patch.GetBary( 1), det);
            }
        }
   }
    M.Build();
}

void SetupLBP1OnTriangle (InterfacePatchCL& patch, int tri, Point3DCL grad[4], double coup[4][4])
{
    Point3DCL surfgrad[4];
    for (int i= 0; i < 4; ++i)
        surfgrad[i]= patch.ApplyProj( grad[i]);
    for (int i= 0; i < 4; ++i)
        for (int j= 0; j <= i; ++j) {
            const double cLB= /*area of reference-triangle*/0.5*
                inner_prod( surfgrad[i], surfgrad[j])*patch.GetFuncDet( tri);
            coup[j][i]+= cLB;
            coup[i][j]= coup[j][i];
        }
}

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, double D)
{
    const IdxT num_unks= mat->RowIdx->NumUnknowns();
    MatrixBuilderCL M( &mat->Data, num_unks, num_unks);
    const Uint lvl = mat->GetRowLevel();

    IdxT num[4];

    std::cerr << "entering SetupLBP1: " << num_unks << " dof.\n";

    Point3DCL grad[4];

    double coup[4][4];
    double dummy;

    InterfacePatchCL patch;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *mat->RowIdx);
            P1DiscCL::GetGradients( grad, dummy, *it);
            std::memset( coup, 0, 4*4*sizeof( double));

            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri)
                    SetupLBP1OnTriangle( patch, tri, grad, coup);
            }

            for(int i= 0; i < 4; ++i) {// assemble row Numb[i]
                if (num[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (num[j] == NoIdx) continue;
                    M( num[i],   num[j])+= coup[j][i];
                }
            }
        }
    }
    M.Build();
    mat->Data*= D; // diffusion coefficient
    std::cerr << mat->Data.num_nonzeros() << " nonzeros in A_LB" << std::endl;
}

void SetupConvectionP1OnTriangle ( const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> qp1[4],
    const LocalP2CL<Point3DCL>& u, Point3DCL grad[4], double coup[4][4])
{
    for (int i= 0; i < 4; ++i)
        qp1[i].assign( p1[i], triangle);
    Quad5_2DCL<Point3DCL> qu( u, triangle);

    for (int i= 0; i < 4; ++i)
        for (int j= 0; j < 4; ++j) {
            const double c= Quad5_2DCL<>( dot( grad[j], qu)*qp1[i]).quad( det);
            coup[i][j]+= c;
        }
}

void SetupMassDivP1OnTriangle (const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> qp1[4],
    const LocalP2CL<Point3DCL>& u, LocalP1CL<Point3DCL> gradp2[10], const Point3DCL& n, double coup[4][4])
{
    for (int i= 0; i < 4; ++i)
        qp1[i].assign( p1[i], triangle);

    Quad5_2DCL<Point3DCL> qgradp2i;
    Quad5_2DCL<> qdivgamma_u;
    for (int i= 0; i < 10; ++i) {
        qgradp2i.assign( gradp2[i], triangle);
        qdivgamma_u+= dot(u[i], qgradp2i) - inner_prod( n, u[i])*dot( n, qgradp2i);
    } // Now qdivgamma_u contains the surface-divergence of u.

    for (int i= 0; i < 4; ++i)
        for (int j= 0; j < 4; ++j)
            coup[i][j]+= Quad5_2DCL<>( qdivgamma_u*qp1[i]*qp1[j]).quad( det);
}

void SetupInterfaceRhsP1OnTriangle (const LocalP1CL<> p1[4],
    Quad5_2DCL<> q[4],VectorCL& v, const IdxT Numb[4],
    const TetraCL& t, const BaryCoordCL triangle[3], double det,
    instat_scalar_fun_ptr f)
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);
    Quad5_2DCL<> qf( t, triangle, f), r;

    for (int i= 0; i < 4; ++i) {
        if (Numb[i] == NoIdx) continue;
        r= qf*q[i];
        v[Numb[i]]+= r.quad( det);
    }
}

void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls,instat_scalar_fun_ptr f)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    IdxT num[4];

    std::cerr << "entering SetupInterfaceRhsP1: " << num_unks << " dof.\n";

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;

    InterfacePatchCL patch;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri)
                    SetupInterfaceRhsP1OnTriangle( p1, q, v->Data, num,
                        *it, &patch.GetBary( tri), patch.GetFuncDet( tri), f);
            }
        }
    }
    std::cerr << " Rhs set up." << std::endl;
}


void
InterfaceP1RepairCL::post_refine ()
{
    VecDescCL loc_u;
    IdxDescCL loc_idx( P1_FE);

    loc_idx.CreateNumbering( fullp1idx_.TriangLevel(), mg_);
    loc_u.SetIdx( &loc_idx);
    DROPS::NoBndDataCL<> dummy;
    RepairAfterRefineP1( make_P1Eval( mg_, dummy, fullu_), loc_u);
    fullu_.Clear();
    fullp1idx_.swap( loc_idx);
    loc_idx.DeleteNumbering( mg_);
    fullu_.SetIdx( &fullp1idx_);
    fullu_.Data= loc_u.Data;
}

void
InterfaceP1RepairCL::pre_refine_sequence ()
{
    fullp1idx_.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    fullu_.SetIdx( &fullp1idx_);
    Extend( mg_, u_, fullu_);
}

void
InterfaceP1RepairCL::post_refine_sequence ()
{
    u_.RowIdx->DeleteNumbering( mg_);
    u_.RowIdx->CreateNumbering( fullp1idx_.TriangLevel(), mg_, &lset_vd_);
    u_.SetIdx( u_.RowIdx);

    Restrict( mg_, fullu_, u_);

    fullp1idx_.DeleteNumbering( mg_);
    fullu_.Clear();
}

void
Ensight6IfaceScalarCL::put (Ensight6OutCL& cf) const
{
    IdxDescCL p1idx;
    p1idx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    VecDescCL p1u( &p1idx);
    Extend( mg_, u_, p1u);
    BndDataCL<> bnd( 0);

    cf.putScalar( make_P1Eval( mg_, bnd, p1u), varName());

    p1idx.DeleteNumbering( mg_);
}

} // end of namespace DROPS
