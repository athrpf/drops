/// \file transport2phase.cpp
/// \brief Classes that constitute a 2-phase-transport-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt; SC RWTH Aachen:

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

#include "poisson/transport2phase.h"

namespace DROPS
{

void TransportP1CL::Init (instat_scalar_fun_ptr cneg, instat_scalar_fun_ptr cpos)
{
    const Uint lvl= ct.GetLevel(),
               idx= ct.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( idx)) {
            if (lset_.Phi.Data[it->Unknowns( lset_.Phi.RowIdx->GetIdx())] <= 0.)
                ct.Data[it->Unknowns( idx)]= H_*cneg( it->GetCoord(), 0.);
            else
                ct.Data[it->Unknowns( idx)]= cpos( it->GetCoord(), 0.);
        }
    }
}

void TransportP1CL::SetTimeStep (double dt, double theta)
{
    dt_= dt;
    if (theta >= 0. && theta <= 1.) theta_= theta;
}


void TransportP1CL::SetupLocalSystem (const TetraCL& t,
    double coupM[4][4], double coupA[4][4], double coupC[4][4], const double time,
    const LocalP2CL<> p2[4], const LocalP2CL<> pipj[4][4], const Quad5CL<> p[4]) const
{
    static LocalP2CL<> one( 1.);

    P2EvalCL<SVectorCL<3>, const VelBndDataT, VecDescCL> u( v_, &Bnd_v_, &MG_);

    double det;
    SMatrixCL<3,4> G;
    P1DiscCL::GetGradients( G, det, t);
    const double absdet= std::fabs( det);
    const SMatrixCL<4,4> GTG( GramMatrix( G));

    InterfaceTetraCL cut;
    cut.Init( t, lset_.Phi, lset_.GetBndData());
    if (!cut.Intersects()) {
        const Quad5CL<Point3DCL> u_loc( t, u, time);
        const double coeff_d = cut.GetSign( 0) == 1 ? D_[0] : D_[1];
        const double coeff_h = cut.GetSign( 0) == 1 ? 1. : 1./H_;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < i; ++j) {
                coupM[j][i]= coeff_h*P1DiscCL::GetMass( i, j)*absdet;
                coupM[i][j]= coupM[j][i];
                coupA[j][i]= coeff_h*coeff_d*GTG( i, j)/6.0*absdet;
                coupA[i][j]= coupA[j][i];
                coupC[j][i]= coeff_h*Quad5CL<>( dot( u_loc, Quad5CL<Point3DCL>( G.col( j)))*p[i]).quad( absdet);
                coupC[i][j]= coeff_h*Quad5CL<>( dot( u_loc, Quad5CL<Point3DCL>( G.col( i)))*p[j]).quad( absdet);
            }
            coupM[i][i]= coeff_h*P1DiscCL::GetMass( i, i)*absdet;
            coupA[i][i]= coeff_h*coeff_d*GTG( i, i)/6.0*absdet;
            coupC[i][i]= coeff_h*Quad5CL<>( dot( u_loc, Quad5CL<Point3DCL>( G.col( i)))*p[i]).quad( absdet);
        }
    }
    else {
        const LocalP2CL<Point3DCL> u_loc( t, u, time);
        LocalP2CL<> convection_interface;
        double iAp, iAn, iMp, iMn, iCp, iCn, integralp, integraln;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < 4; ++j) {
                integralp= integraln= 0.;
                iAp= iAn = iMp = iMn = iCp = iCn = 0.;
                convection_interface= dot( u_loc, LocalP2CL<Point3DCL>(G.col( j)))*p2[i];
                for (int ch= 0; ch < 8; ++ch) {
                    cut.ComputeCutForChild( ch);
                    cut.quadBothParts( integralp, integraln, pipj[i][j], absdet);
                    iMp+= integralp; iMn+= integraln;
                    cut.quadBothParts( integralp, integraln, one, absdet);
                    iAp+= integralp; iAn+= integraln;
                    cut.quadBothParts( integralp, integraln, convection_interface, absdet);
                    iCp+= integralp; iCn+= integraln;
                }
                coupM[j][i]= iMp + iMn/H_;
                coupA[j][i]= (iAp*D_[0] + iAn*D_[1]/H_)* GTG( j, i);
                coupC[j][i]= iCp + iCn/H_;
            }
        }
    }
}

void TransportP1CL::SetupInstatSystem (MatrixCL& matA, VecDescCL* cplA,
                        MatrixCL& matM, VecDescCL* cplM, MatrixCL& matC, VecDescCL* cplC,
                        IdxDescCL& RowIdx, const double time) const
{
    if (cplM != 0)
    {
        cplM->Data= 0.;
        cplA->Data= 0.;
        cplC->Data= 0.;
    }
    matM.clear();
    matA.clear();
    matC.clear();
    const IdxT num_unks=  RowIdx.NumUnknowns();
    MatrixBuilderCL A(&matA, num_unks,  num_unks),//diffusion
                    M(&matM, num_unks,  num_unks),//mass matrix
                    C(&matC, num_unks,  num_unks);// convection

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL n;
    double coupA[4][4], coupM[4][4], coupC[4][4];

    // The 16 products of the P1-shape-functions
    LocalP2CL<> p2[4], pipj[4][4], convection_interface;
    Quad5CL<> p[4];
    for(int i= 0; i < 4; ++i) {
        LocalP1CL<> p1;
        p1[i]= 1.;
        p[i].assign( p1);
        p2[i].assign( p1);
    }

    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
            pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
        }
        pipj[i][i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
                pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {
        SetupLocalSystem ( *sit, coupM, coupA, coupC, time, p2, pipj, p);
        n.assign( *sit, RowIdx, Bnd_);
        // write values into matrix
        for(int i= 0; i < 4; ++i)
            if (n.WithUnknowns( i))
                for(int j= 0; j < 4; ++j)
                    if (n.WithUnknowns( j)) {
                        M( n.num[i], n.num[j])+= coupM[j][i];
                        A( n.num[i], n.num[j])+= coupA[j][i];
                        C( n.num[i], n.num[j])+= coupC[j][i];
                    }
                    else if (cplM != 0) {
                        const double val= Bnd_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time);
                        cplM->Data[n.num[i]]-= coupM[j][i]*val;
                        cplA->Data[n.num[i]]-= coupA[j][i]*val;
                        cplC->Data[n.num[i]]-= coupC[j][i]*val;
                    }
    }
    A.Build();
    M.Build();
    C.Build();
}

void TransportP1CL::SetupInstatSystem (MLMatDescCL& matA, VecDescCL& cplA,
    MLMatDescCL& matM, VecDescCL& cplM, MLMatDescCL& matC, VecDescCL& cplC, const double time) const
{
    MLMatrixCL::iterator itA = matA.Data.begin();
    MLMatrixCL::iterator itM = matM.Data.begin();
    MLMatrixCL::iterator itC = matC.Data.begin();
    MLIdxDescCL::iterator it = matA.RowIdx->begin();
    for (size_t lvl=0; lvl < matA.Data.size(); ++lvl, ++itA, ++itM, ++itC, ++it)
        if (lvl != 0)
            SetupInstatSystem (*itA, &cplA, *itM, &cplM, *itC, &cplC, *it, time);
        else
            SetupInstatSystem (*itA, 0, *itM, 0, *itC, 0, *it, time);
}

void TransportP1CL::Update()
{
    MLIdxDescCL* cidx= &idx;

    c.SetIdx( cidx);
    ct2c();

    cplM.SetIdx( cidx);
    cplA.SetIdx( cidx);
    cplC.SetIdx( cidx);

    oldcplM.SetIdx( cidx);
    oldcplA.SetIdx( cidx);
    oldcplC.SetIdx( cidx);

    M.Data.clear();
    M.SetIdx( cidx, cidx);
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    C.Data.clear();
    C.SetIdx( cidx, cidx);

    SetupInstatSystem( A, oldcplA, M, oldcplM, C, oldcplC, c.t);
}

void TransportP1CL::InitStep (VectorCL& rhs)
{
    VectorCL rhs1( (1./dt_)*ct.Data);
    if (theta_ != 1.0) {
        VectorCL tmp( (1. - theta_)*(oldcplA.Data - A.Data*ct.Data + oldcplC.Data - C.Data*ct.Data));
        VectorCL rhs2( tmp.size());
        gm_.Solve( M.Data, rhs2, tmp);
        std::cout << "Inverting M_old: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
        rhs1+= rhs2;
    }
    SetupInstatSystem( A, cplA, M, cplM, C, cplC, c.t);
    // Todo: Add 1/dt_*M_new(c_{dirichlet,old}), if there are time-dependant Dirichlet-BC.
    rhs= M.Data*rhs1 + /*Todo: add together with the above (1./dt_)*cplM.Data + */ theta_*(cplA.Data + cplC.Data);

    MLMatrixCL L;
    L.LinComb( theta_, A.Data, theta_, C.Data);
    L_.clear();
    L_.LinComb( 1./dt_, M.Data, 1., L);
}

void TransportP1CL::DoStep (const VectorCL& rhs)
{
    std::cout << "Before solve: res = " << norm( L_*ct.Data - rhs) << std::endl;
    // std::cout << "ct:\n" << ct.Data << "\ncplC:\n" << cplC.Data << "\nC:\n" << C.Data << std::endl;
    // std::cout << "\ncplM:\n" << cplM.Data << "\nM:\n" << M.Data << std::endl;
    // std::cout << "\ncplA:\n" << cplA.Data << "\nA:\n" << A.Data << std::endl;
    gm_.Solve( L_, ct.Data, rhs);
    std::cout << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void TransportP1CL::CommitStep ()
{
    std::swap( cplM, oldcplM);
    std::swap( cplA, oldcplA);
    std::swap( cplC, oldcplC);
    ct2c();
}

void TransportP1CL::DoStep (double new_t)
{
    VectorCL rhs( c.Data.size());
    c.t= new_t;
    InitStep( rhs);
    DoStep( rhs);
    CommitStep();
}

void TransportP1CL::c2ct()
// For lset <= 0, multiply with H
{
    Uint cidx= c.RowIdx->GetIdx(),
         phiidx= lset_.Phi.RowIdx->GetIdx();
    ct.Data= c.Data;
    DROPS_FOR_TRIANG_VERTEX( MG_, c.RowIdx->TriangLevel(), sit)
        if (lset_.Phi.Data[sit->Unknowns( phiidx)] <= 0. && sit->Unknowns.Exist( cidx))
            ct.Data[sit->Unknowns( cidx)]*= H_;
}

void TransportP1CL::ct2c()
// For lset <= 0, divide by H
{
    Uint ctidx= ct.RowIdx->GetIdx(),
         phiidx= lset_.Phi.RowIdx->GetIdx();
    c.Data= ct.Data;
    DROPS_FOR_TRIANG_VERTEX( MG_, ct.RowIdx->TriangLevel(), sit)
        if (lset_.Phi.Data[sit->Unknowns( phiidx)] <= 0. && sit->Unknowns.Exist( ctidx))
            c.Data[sit->Unknowns( ctidx)]/= H_;
}

//*****************************************************************************
//                               TransportRepairCL
//*****************************************************************************
void
TransportRepairCL::post_refine ()
{
    VecDescCL loc_ct;
    IdxDescCL loc_cidx( P1_FE);
    VecDescCL& ct= c_.ct;
    match_fun match= mg_.GetBnd().GetMatchFun();

    loc_cidx.CreateNumbering( mg_.GetLastLevel(), mg_, c_.GetBndData(), match);
    loc_ct.SetIdx( &loc_cidx);
    RepairAfterRefineP1( c_.GetSolution( ct), loc_ct);

    ct.Clear( c_.c.t);
    ct.RowIdx->DeleteNumbering( mg_);
    c_.idx.GetFinest().swap( loc_cidx);
    ct.SetIdx( &c_.idx);
    ct.Data= loc_ct.Data;
}


} // end of namespace DROPS
