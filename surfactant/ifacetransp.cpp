/// \file
/// \brief Discretization for PDEs on an interface.

#include "surfactant/ifacetransp.h"
#include "levelset/levelset.h"
#include "num/spmat.h"
#include <cstring>


namespace DROPS {

void CreateNumbOnInterfaceVertex (const Uint idx, IdxT& counter, Uint stride,
    const MultiGridCL::TriangTetraIteratorCL& begin,
    const MultiGridCL::TriangTetraIteratorCL& end,
    const VecDescCL& ls)
{
    if (stride == 0) return;

    InterfacePatchCL p;
    for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
        p.Init( *it, ls);
        if (!p.Intersects()) continue;

        const bool innercut( p.IntersectsInterior());
        for (Uint i= 0; i < NumVertsC; ++i) {
            UnknownHandleCL& u= const_cast<VertexCL*>( it->GetVertex( i))->Unknowns;
            if (innercut || p.GetSign( i) == 0) {
                u.Prepare( idx);
                if ( u( idx) == NoIdx) {
                    u( idx)= counter;
                    counter+= stride;
                }
            }
            else
                if (!u.Exist( idx)) {
                    u.Prepare( idx);
                    u( idx)= NoIdx;
                }
        }
    }
}

void CreateNumbOnInterface(Uint level, IdxDescCL& idx, MultiGridCL& mg,
    const VecDescCL& ls)
{
    // set up the index description
    idx.TriangLevel = level;
    idx.NumUnknowns = 0;

    const Uint idxnum= idx.GetIdx();
    // allocate space for indices; number unknowns in TriangLevel level
    if (idx.NumUnknownsVertex() != 0)
        CreateNumbOnInterfaceVertex( idxnum, idx.NumUnknowns, idx.NumUnknownsVertex(),
            mg.GetTriangTetraBegin( level), mg.GetTriangTetraEnd( level), ls);

    if (idx.NumUnknownsEdge() != 0 || idx.NumUnknownsFace() != 0 || idx.NumUnknownsTetra() != 0)
        throw DROPSErrCL( "CreateNumbOnInterface: Only vertex unknowns are implemented\n" );
}

void Extend (const MultiGridCL& mg, const VecDescCL& x, VecDescCL& xext)
{
    const Uint xidx( x.RowIdx->GetIdx()),
        xextidx( xext.RowIdx->GetIdx()),
        lvl( x.RowIdx->TriangLevel);
    xext.Data= 0.;

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( xidx) && it->Unknowns.Exist( xextidx))
            xext.Data[it->Unknowns( xextidx)]= x.Data[it->Unknowns( xidx)];
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
    const IdxT num_unks=  matM->RowIdx->NumUnknowns;
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

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls)
{
    const IdxT num_unks= mat->RowIdx->NumUnknowns;
    MatrixBuilderCL M( &mat->Data, num_unks, num_unks);
    const Uint lvl = mat->GetRowLevel();

    IdxT num[4];

    std::cerr << "entering SetupLBP1: " << num_unks << " dof.\n";

    Point3DCL grad[4];
    SMatrixCL<3,3> T;

    double coup[4][4];
    double det, dummy;

    InterfacePatchCL patch;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *mat->RowIdx);
            GetTrafoTr( T, det, *it);
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
    std::cerr << mat->Data.num_nonzeros() << " nonzeros in A_LB" << std::endl;
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
    const IdxT num_unks= v->RowIdx->NumUnknowns;
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

} // end of namespace DROPS
