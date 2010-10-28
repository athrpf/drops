/// \file
/// \brief Discretization for PDEs on an interface.
/// \author LNM RWTH Aachen: Joerg Grande

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

/// Implementation of Setup-Routines and class-methods.

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

void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    const IdxT num_unks=  matM->RowIdx->NumUnknowns();
    MatrixBuilderCL M( &matM->Data, num_unks,  num_unks);

    const Uint lvl= matM->GetRowLevel();
    IdxT Numb[4];

    double det;
    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;

    InterfaceTriangleCL triangle;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, it) {
        triangle.Init( *it, ls, lsetbnd);

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

void SetupLBP1OnTriangle (InterfaceTriangleCL& triangle, int tri, Point3DCL grad[4], double coup[4][4])
{
    Point3DCL surfgrad[4];
    for (int i= 0; i < 4; ++i)
        surfgrad[i]= triangle.ApplyProj( grad[i]);
    for (int i= 0; i < 4; ++i)
        for (int j= 0; j <= i; ++j) {
            const double cLB= /*area of reference-triangle*/0.5*
                inner_prod( surfgrad[i], surfgrad[j])*triangle.GetAbsDet( tri);
            coup[j][i]+= cLB;
            coup[i][j]= coup[j][i];
        }
}

void SetupLBP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double D)
{
    const IdxT num_rows= mat->RowIdx->NumUnknowns();
    const IdxT num_cols= mat->ColIdx->NumUnknowns();
    MatrixBuilderCL M( &mat->Data, num_rows, num_cols);
    const Uint lvl = mat->GetRowLevel();

    IdxT numr[4], numc[4];

    std::cout << "entering SetupLBP1: " << num_rows << " rows, " << num_cols << " cols. ";

    Point3DCL grad[4];

    double coup[4][4];
    double dummy;

    InterfaceTriangleCL triangle;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( numr, *it, *mat->RowIdx);
            GetLocalNumbP1NoBnd( numc, *it, *mat->ColIdx);
            P1DiscCL::GetGradients( grad, dummy, *it);
            std::memset( coup, 0, 4*4*sizeof( double));

            for (int ch= 0; ch < 8; ++ch) {
            	triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    SetupLBP1OnTriangle( triangle, tri, grad, coup);
            }

            for(int i= 0; i < 4; ++i) {// assemble row Numb[i]
                if (numr[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (numc[j] == NoIdx) continue;
                    M( numr[i],   numc[j])+= coup[j][i];
                }
            }
        }
    }
    M.Build();
    mat->Data*= D; // diffusion coefficient
    std::cout << mat->Data.num_nonzeros() << " nonzeros in A_LB" << std::endl;
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

void SetupMixedMassP1OnTriangle (const BaryCoordCL triangle[3], double det,
    const LocalP1CL<> p1[4], Quad5_2DCL<> q[4], double coup[4][4])
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);

    Quad5_2DCL<> m;
    for (int i= 0; i < 4; ++i) {
        m= q[i]*q[i];
        coup[i][i]+= m.quad( det);
        for(int j= 0; j < i; ++j) {
            m= q[j]*q[i];
            coup[i][j]+= m.quad( det);
            coup[j][i]+= m.quad( det);
        }
    }
}

void SetupMixedMassP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const BndDataCL<>& lsetbnd)
{
    const IdxT rows= mat->RowIdx->NumUnknowns(),
               cols= mat->ColIdx->NumUnknowns();
    MatrixBuilderCL m( &mat->Data, rows, cols);
    const Uint lvl= mat->GetRowLevel();
    IdxT rownum[4], colnum[4];

    std::cerr << "entering SetupMixedMassP1: " << rows << " rows, " << cols << " cols. ";

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> qp1[4];

    double coup[4][4];
    InterfaceTriangleCL triangle;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls, lsetbnd);
        if (!triangle.Intersects()) continue; // We are at the phase boundary.

        GetLocalNumbP1NoBnd( rownum, *it, *mat->RowIdx);
        GetLocalNumbP1NoBnd( colnum, *it, *mat->ColIdx);
        std::memset( coup, 0, 4*4*sizeof( double));

        for (int ch= 0; ch < 8; ++ch) {
        	triangle.ComputeForChild( ch);
            for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                SetupMixedMassP1OnTriangle ( &triangle.GetBary( tri), triangle.GetAbsDet( tri), p1, qp1, coup);
        }

        for(int i= 0; i < 4; ++i) {// assemble row Numb[i]
            if (rownum[i] == NoIdx) continue;
            for(int j= 0; j < 4; ++j) {
                if (colnum[j] == NoIdx) continue;
                m( rownum[i], colnum[j])+= coup[i][j];
            }
        }
    }
    m.Build();
    std::cerr << mat->Data.num_nonzeros() << " nonzeros in mixed mass-divergence matrix!" << std::endl;
}

void SetupInterfaceRhsP1 (const MultiGridCL& mg, VecDescCL* v,
    const VecDescCL& ls, const BndDataCL<>& lsetbnd, instat_scalar_fun_ptr f)
{
    const IdxT num_unks= v->RowIdx->NumUnknowns();
    const Uint lvl = v->GetLevel();

    IdxT num[4];

    std::cout << "entering SetupInterfaceRhsP1: " << num_unks << " dof... ";

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;

    InterfaceTriangleCL triangle;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        triangle.Init( *it, ls, lsetbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *v->RowIdx);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
                    SetupInterfaceRhsP1OnTriangle( p1, q, v->Data, num,
                        *it, &triangle.GetBary( tri), triangle.GetAbsDet( tri), f);
            }
        }
    }
    std::cout << " Rhs set up." << std::endl;
}

/// \todo This should be a generic function somewhere in num or misc.
void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, MultiGridCL& mg, double t)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( idx))
            ic.Data[it->Unknowns( idx)]= icf( it->GetCoord(), t);
    }
}

void SurfactantcGP1CL::Init (instat_scalar_fun_ptr icf)
{
    P1Init ( icf, ic, MG_, 0.);
}

void SurfactantcGP1CL::SetTimeStep (double dt, double theta)
{
    dt_= dt;
    if (theta >= 0. && theta <= 1.) theta_= theta;
}

void SurfactantcGP1CL::Update()
{
    // std::cout << "SurfactantcGP1CL::Update:\n";
    IdxDescCL* cidx= ic.RowIdx;

    M.Data.clear();
    M.SetIdx( cidx, cidx);
    DROPS::SetupInterfaceMassP1( MG_, &M, lset_vd_, lsetbnd_);
    // std::cout << "M is set up.\n";
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    DROPS::SetupLBP1( MG_, &A, lset_vd_, lsetbnd_, D_);
    // std::cout << "A is set up.\n";
    C.Data.clear();
    C.SetIdx( cidx, cidx);
    DROPS::SetupConvectionP1( MG_, &C, lset_vd_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, *v_));
    // std::cout << "C is set up.\n";
    Md.Data.clear();
    Md.SetIdx( cidx, cidx);
    DROPS::SetupMassDivP1( MG_, &Md, lset_vd_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, *v_));
    // std::cout << "Md is set up.\n";

    if (theta_ != 1.0) {
        M2.Data.clear();
        M2.SetIdx( cidx, cidx);
        DROPS::SetupInterfaceMassP1( MG_, &M2, oldls_, lsetbnd_);
        // std::cout << "M2 is set up.\n";
    }
    std::cout << "SurfactantP1CL::Update: Finished\n";
}

VectorCL SurfactantcGP1CL::InitStep ()
{

    // std::cout << "SurfactantcGP1CL::InitStep:\n";
    idx.CreateNumbering( oldidx_.TriangLevel(), MG_, &lset_vd_, &lsetbnd_); // InitOld deletes oldidx_ and swaps idx and oldidx_.
    std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;
    ic.SetIdx( &idx);

    MatDescCL m( &idx, &oldidx_);
    DROPS::SetupMixedMassP1( MG_, &m, lset_vd_, lsetbnd_);
    // std::cout << "mixed M on new interface is set up.\n";
    VectorCL rhs( theta_*(m.Data*oldic_));

    if (theta_ == 1.0) return rhs;

    m.Data.clear();
    DROPS::SetupMixedMassP1( MG_, &m, oldls_, lsetbnd_);
    // std::cout << "mixed M on old interface is set up.\n";
    rhs+= (1. - theta_)*(m.Data*oldic_);

    m.Data.clear();
    DROPS::SetupLBP1( MG_, &m, oldls_, lsetbnd_, D_);
    // std::cout << "mixed A on old interface is set up.\n";
    VectorCL rhs2( m.Data*oldic_);
    m.Data.clear();
    DROPS::SetupConvectionP1( MG_, &m, oldls_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, oldv_));
    // std::cout << "mixed C on old interface is set up.\n";
    rhs2+= m.Data*oldic_;
    m.Data.clear();
    DROPS::SetupMassDivP1( MG_, &m, oldls_, lsetbnd_, make_P2Eval( MG_, Bnd_v_, oldv_));
    // std::cout << "mixed Md on old interface is set up.\n";
    rhs2+= m.Data*oldic_;

    return VectorCL( rhs - ((1. - theta_)*dt_)*rhs2);
}

void SurfactantcGP1CL::DoStep (const VectorCL& rhs)
{
    Update();

    if (theta_ == 1.)
        L_.LinComb( theta_, M.Data, dt_*theta_, A.Data, dt_*theta_, Md.Data, dt_*theta_, C.Data);
    else {
        MatrixCL m;
        m.LinComb( theta_, M.Data, dt_*theta_, A.Data, dt_*theta_, Md.Data, dt_*theta_, C.Data);
        L_.LinComb( 1., m, 1. - theta_, M2.Data);
    }
    std::cout << "Before solve: res = " << norm( L_*ic.Data - rhs) << std::endl;
    gm_.Solve( L_, ic.Data, rhs);
    std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantcGP1CL::CommitStep ()
{
    return;
}

void SurfactantcGP1CL::DoStep (double new_t)
{
    VectorCL rhs( InitStep());
    ic.t= new_t;
    DoStep( rhs);
    CommitStep();
}

void SurfactantcGP1CL::InitOld ()
{
    if (oldidx_.NumUnknowns() > 0)
        oldidx_.DeleteNumbering( MG_);
    oldidx_.swap( idx);
    oldic_.resize( ic.Data.size());
    oldic_= ic.Data;
    oldls_.RowIdx= lset_vd_.RowIdx;
    oldls_.Data.resize( lset_vd_.Data.size());
    oldls_.Data= lset_vd_.Data;
    oldv_.SetIdx( v_->RowIdx);
    oldv_.Data= v_->Data;
    oldt_= ic.t;
}


void
InterfaceP1RepairCL::post_refine ()
{
    VecDescCL loc_u;
    IdxDescCL loc_idx( P1_FE);

    loc_idx.CreateNumbering( fullp1idx_.TriangLevel(), mg_);
//    loc_idx.CreateNumbering( mg_.GetLastLevel(), mg_);
    loc_u.SetIdx( &loc_idx);
    DROPS::NoBndDataCL<> dummy;

//    P1EvalCL<double, DROPS::NoBndDataCL<>, const VecDescCL> oldsol( &fullu_, &dummy, &mg_);
//    P1EvalCL<double, DROPS::NoBndDataCL<>,       VecDescCL>    sol( &loc_u , &dummy, &mg_);
//    Interpolate( sol, oldsol);
    RepairAfterRefineP1( make_P1Eval( mg_, dummy, fullu_), loc_u);

    fullu_.Clear( fullu_.t);
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
    u_.RowIdx->CreateNumbering( fullp1idx_.TriangLevel(), mg_, &lset_vd_, &lset_bnd_);
    u_.SetIdx( u_.RowIdx);

    Restrict( mg_, fullu_, u_);

    fullp1idx_.DeleteNumbering( mg_);
    fullu_.Clear( u_.t);
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

void
VTKIfaceScalarCL::put (VTKOutCL& cf) const
{
    IdxDescCL p1idx;
    p1idx.CreateNumbering( u_.RowIdx->TriangLevel(), mg_);
    VecDescCL p1u( &p1idx);
    Extend( mg_, u_, p1u);
    BndDataCL<> bnd( 0);

    cf.PutScalar( make_P1Eval( mg_, bnd, p1u), varName());

    p1idx.DeleteNumbering( mg_);
}

} // end of namespace DROPS
