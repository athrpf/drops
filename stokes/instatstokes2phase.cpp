/// \file instatstokes2phase.cpp
/// \brief classes that constitute the 2-phase Stokes problem
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

#include "stokes/instatstokes2phase.h"
#include "num/accumulator.h"

namespace DROPS
{
// -----------------------------------------------------------------------------
//                        Routines for SetupSystem2
// -----------------------------------------------------------------------------


void SetupSystem2_P2P0( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData,
                        MatrixCL* B, VecDescCL* c, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
// P2 / P0 FEs for vel/pr
{
    SparseMatBuilderCL<double, SMatrixCL<1,3> >  mB( B, RowIdx->NumUnknowns(), ColIdx->NumUnknowns());
    if (c != 0) c->Clear( t);
    const Uint lvl= RowIdx->TriangLevel();
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;

    const Uint pidx= RowIdx->GetIdx();

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *ColIdx, BndData.Vel);
        const IdxT prNumbTetra= sit->Unknowns(pidx);

        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
            {
                tmp= Grad[vel].quad( absdet);
                mB( prNumbTetra, n.num[vel])  -=  SMatrixCL<1,3>(tmp);
            }
            else if (c != 0) // put coupling on rhs
            {
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                        : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                c->Data[ prNumbTetra]+= inner_prod( Grad[vel].quad( absdet), tmp);
            }
        }
    }
    mB.Build();
}


void SetupSystem2_P2P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, MatrixCL* B, VecDescCL* c,
                        IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
// P2 / P1 FEs (Taylor-Hood) for vel/pr
{
    SparseMatBuilderCL<double, SMatrixCL<1,3> > mB( B, RowIdx->NumUnknowns(), ColIdx->NumUnknowns());
    if (c != 0) c->Clear( t);
    const Uint lvl= RowIdx->TriangLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *ColIdx, BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  SMatrixCL<1,3>(tmp);
                }
            else if (c != 0)
            { // put coupling on rhs
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                        : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr, absdet), tmp);
            }
        }
    }
    mB.Build();
}

inline void ComputePgradV( LocalP2CL<Point3DCL>& PgradV, Uint pr, const Quad2CL<Point3DCL>& gradV)
{
    PgradV= Point3DCL();
    PgradV[pr]= gradV[pr];
    for (Uint vert=0; vert<4; ++vert)
        if (vert!=pr)
            PgradV[EdgeByVert(pr,vert)+4]= 0.25*(gradV[pr]+gradV[vert]);
}


void SetupSystem2_P2P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
// P2 / P1X FEs (X=extended) for vel/pr
{
    const ExtIdxDescCL& Xidx= RowIdx->GetXidx();
    SparseMatBuilderCL<double, SMatrixCL<1,3> > mB( B, RowIdx->NumUnknowns(), ColIdx->NumUnknowns());
    if (c != 0) c->Clear( t);
    const Uint lvl= RowIdx->TriangLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;
    InterfaceTetraCL cut;
    LocalP2CL<> loc_phi;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *ColIdx, BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  SMatrixCL<1,3>(tmp);
                }
            else if (c != 0) { // put coupling on rhs
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                        : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr, absdet), tmp);
            }
        }

        loc_phi.assign( *sit, lset.Phi, lset.GetBndData());
        cut.Init( *sit, loc_phi);
        if (!cut.Intersects()) continue; // extended basis functions have only support on tetra intersecting Gamma!
        for(int pr=0; pr<4; ++pr) {
            // compute the integrals
            // I = \int_{T_+} grad v_vel p_pr dx  -  C \int_{T}grad v_vel p_pr dx,
            // where C= (sign Phi_pr==1) \in {0,1} and T_+ = T \cap \Omega_2 (positive part)
            LocalP2CL<Point3DCL> PgradV;
            const IdxT xidx= Xidx[prNumb[pr]];
            if (xidx==NoIdx) continue;

            for(int vel=0; vel<10; ++vel)
            {
                ComputePgradV( PgradV, pr, Grad[vel]);

                Point3DCL integral;
                const bool is_pos= cut.GetSign(pr)==1;
                // for C=0 (<=> !is_pos) we have I = -\int_{T_-} grad v_vel p_pr dx
                for (int ch=0; ch<8; ++ch)
                {
                    cut.ComputeCutForChild(ch);
                    integral+= cut.quad( PgradV, absdet, !is_pos); // integrate on other part
                }

                // for C=0 we have I = -\int_{T_-} grad v_vel p_pr dx
                if (is_pos) integral= -integral;

                if (n.WithUnknowns( vel))
                    mB( xidx, n.num[vel])  -=  SMatrixCL<1,3>(integral);
                else if (c != 0)
                { // put coupling on rhs
                    typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                    bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                    tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                            : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                    c->Data[ xidx]+= inner_prod( integral, tmp);
                }
            }
        }
    }
    mB.Build();
}


void SetupSystem2_P2RP1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
// P2X / P1X FEs (X=extended) for vel/pr
{
    const ExtIdxDescCL& prXidx=    RowIdx->GetXidx();
    const ExtIdxDescCL& velXidx= ColIdx->GetXidx();
    SparseMatBuilderCL<double, SMatrixCL<1,3> >  mB( B, RowIdx->NumUnknowns(), ColIdx->NumUnknowns());
    if (c != 0) c->Clear( t);
    const Uint lvl= RowIdx->TriangLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], gradVx;
    LocalP2CL<> velR_p[4][8], velR_n[4][8];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;
    InterfaceTetraCL cut;
    LocalP2CL<> loc_phi;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *ColIdx, BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)

        // do setup for std FEM part
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  SMatrixCL<1,3>(tmp);
                }
            else if (c != 0) { // put coupling on rhs
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                        : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr, absdet), tmp);
            }
        }
        // now handle extended dofs
        loc_phi.assign( *sit, lset.Phi, lset.GetBndData());
        cut.Init( *sit, loc_phi);
        if (!cut.Intersects()) continue; // extended basis functions have only support on tetra intersecting Gamma!

        P2RidgeDiscCL::GetExtBasisOnChildren( velR_p, velR_n, loc_phi);
        for(int pr=0; pr<4; ++pr) {
            // compute the integrals
            // I = \int_{T_+} grad v_vel p_pr dx  -  C \int_{T}grad v_vel p_pr dx,
            // where C= (sign Phi_pr==1) \in {0,1} and T_+ = T \cap \Omega_2 (positive part)
            LocalP2CL<Point3DCL> PgradV;
            const IdxT xidx= prXidx[prNumb[pr]];

            for(int vel= (xidx!=NoIdx ? 0 : 10); vel<14; ++vel)
            {
                const IdxT stdvidx= vel<10 ? n.num[vel] : n.num[vel-10],
                        xvidx= (vel>=10 && stdvidx!=NoIdx) ? velXidx[stdvidx] : NoIdx,
                        vidx= vel<10 ? stdvidx : xvidx;
                if (vel>=10 && xvidx==NoIdx) continue; // no extended vel dof

                Point3DCL int_Px_gradV, int_P_gradV;
                const bool is_pos= cut.GetSign(pr)==1;

                if (vel<10) { // standard P2 FE
                    ComputePgradV( PgradV, pr, Grad[vel]);
                    for (int ch=0; ch<8; ++ch)
                    {
                        cut.ComputeCutForChild(ch);
                        int_Px_gradV+= cut.quad( PgradV, absdet, !is_pos); // integrate on other part
                    }
                    // for C=1 (<=> is_pos) we have I = -\int_{T_-} grad v_vel p_pr dx
                    if (is_pos) int_Px_gradV= -int_Px_gradV;
                } else { // ridge enrichment for vel
                    Point3DCL intNeg, intPos;
                    for (int ch=0; ch<8; ++ch)
                    {
                        cut.ComputeCutForChild(ch);
                        // integrate on pos. part
                        P2DiscCL::GetFuncGradient( gradVx, velR_p[vel-10][ch], Grad);
                        ComputePgradV( PgradV, pr, gradVx);
                        intPos+= cut.quad( PgradV, absdet, true);
                        // integrate on neg. part
                        P2DiscCL::GetFuncGradient( gradVx, velR_n[vel-10][ch], Grad);
                        ComputePgradV( PgradV, pr, gradVx);
                        intNeg+= cut.quad( PgradV, absdet, false);
                    }
                    int_P_gradV= intPos + intNeg;
                    // for C=1 (<=> is_pos) we have I = -\int_{T_-} grad v_vel p_pr dx
                    int_Px_gradV= is_pos ? -intNeg : intPos;
                }

                if (vidx!=NoIdx)
                {
                    if (xidx!=NoIdx) {
                        mB( xidx, vidx)  -=  SMatrixCL<1,3>(int_Px_gradV);
                    }
                    if (vel>=10) { // extended vel: write also entry for p_gradVx
                        const IdxT pidx= prNumb[pr];
                        mB( pidx, vidx)  -=  SMatrixCL<1,3>(int_P_gradV);
                    }
                }
                else if (c != 0 && vel<10)
                { // put coupling on rhs
                    typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                    bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                    tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                            : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                    c->Data[ xidx]+= inner_prod( int_Px_gradV, tmp);
                }
            }
        }
    }
    mB.Build();
}


void SetupSystem2_P2RP1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, MatrixCL* B, VecDescCL* c, const LevelsetP2CL& lset, IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
// P2X / P1 FEs (X=extended) for vel/pr
{
    const ExtIdxDescCL& velXidx= ColIdx->GetXidx();
    SparseMatBuilderCL<double, SMatrixCL<1,3> > mB( B, RowIdx->NumUnknowns(), ColIdx->NumUnknowns());
    if (c != 0) c->Clear( t);
    const Uint lvl= RowIdx->TriangLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], gradVx;
    LocalP2CL<> velR_p[4][8], velR_n[4][8];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;
    InterfaceTetraCL cut;
    LocalP2CL<> loc_phi;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *ColIdx, BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)

        // do setup for std FEM part
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  SMatrixCL<1,3>(tmp);
                }
            else if (c != 0) { // put coupling on rhs
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                        : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr, absdet), tmp);
            }
        }
        // now handle extended dofs
        loc_phi.assign( *sit, lset.Phi, lset.GetBndData());
        cut.Init( *sit, loc_phi);
        if (!cut.Intersects()) continue; // extended basis functions have only support on tetra intersecting Gamma!

        P2RidgeDiscCL::GetExtBasisOnChildren( velR_p, velR_n, loc_phi);
        for(int pr=0; pr<4; ++pr) {
            // compute the integrals
            // I = \int_{T_+} grad v_vel p_pr dx  -  C \int_{T}grad v_vel p_pr dx,
            // where C= (sign Phi_pr==1) \in {0,1} and T_+ = T \cap \Omega_2 (positive part)
            LocalP2CL<Point3DCL> PgradV;

            for(int xvel=0; xvel<4; ++xvel)
            {
                const IdxT stdvidx= n.num[xvel],
                    xvidx= stdvidx!=NoIdx ? velXidx[stdvidx] : NoIdx;
                if (xvidx==NoIdx) continue; // no extended vel dof

                Point3DCL int_P_gradV;

                Point3DCL intNeg, intPos;
                for (int ch=0; ch<8; ++ch)
                {
                    cut.ComputeCutForChild(ch);
                    // integrate on pos. part
                    P2DiscCL::GetFuncGradient( gradVx, velR_p[xvel][ch], Grad);
                    ComputePgradV( PgradV, pr, gradVx);
                    intPos+= cut.quad( PgradV, absdet, true);
                    // integrate on neg. part
                    P2DiscCL::GetFuncGradient( gradVx, velR_n[xvel][ch], Grad);
                    ComputePgradV( PgradV, pr, gradVx);
                    intNeg+= cut.quad( PgradV, absdet, false);
                }
                int_P_gradV= intPos + intNeg;

                // extended vel: write entry for p_gradVx
                const IdxT pidx= prNumb[pr];
                mB( pidx, xvidx)  -=  SMatrixCL<1,3>(int_P_gradV);
            }
        }
    }
    mB.Build();
}


void SetupSystem2_P2P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, MatrixCL* B, VecDescCL* c,
                         IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
// P2 / P1D FEs for vel/pr
{
    SparseMatBuilderCL<double, SMatrixCL<1,3> > mB( B, RowIdx->NumUnknowns(), ColIdx->NumUnknowns());
    if (c != 0) c->Clear( t);
    const Uint lvl= RowIdx->TriangLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
        send=MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *ColIdx, BndData.Vel);
        GetLocalNumbP1DNoBnd( prNumb, *sit, *RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1D( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  SMatrixCL<1,3>(tmp);
                }
            else if (c != 0)
            { // put coupling on rhs
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                        : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1D( pr, absdet), tmp);
            }
        }
    }
    mB.Build();
}


// -----------------------------------------------------------------------------
//                        Routines for SetupRhs2
// -----------------------------------------------------------------------------


void SetupRhs2_P2P0( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t)
{
    c->Clear( t);
    const Uint lvl= c->GetLevel();
    const Uint pidx= c->RowIdx->GetIdx();

    bool IsOnDirBnd[10];
    Quad2CL<Point3DCL> Grad_vel, GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);

  for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        // collect some bnd information about the edges and verts of the tetra
        // and save it in IsOnDirBnd
        for (int i=0; i<4; ++i)
            IsOnDirBnd[i]= BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) );
        for (int i=0; i<6; ++i)
            IsOnDirBnd[i+4]= BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) );
        const IdxT prNumbTetra= sit->Unknowns(pidx);

        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (IsOnDirBnd[vel]) { // put coupling on rhs
                P2DiscCL::GetGradient( Grad_vel, GradRef[vel], T);
                tmp= vel<4 ? BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                c->Data[ prNumbTetra]+= inner_prod( Grad_vel.quad( absdet), tmp);
            }
        }
    }
}


void SetupRhs2_P2P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t)
{
    c->Clear( t);

    const Uint lvl = c->GetLevel();
    const Uint pidx= c->RowIdx->GetIdx();

    IdxT prNumb[4];
    bool IsOnDirBnd[10];

    Quad2CL<Point3DCL> Grad_vel, GradRef[10];

    SMatrixCL<3,3> T;

    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        // collect some information about the edges and verts of the tetra
        // and save it in prNumb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            IsOnDirBnd[i]= BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) );
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for (int i=0; i<6; ++i)
            IsOnDirBnd[i+4]= BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) );

        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        // b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (IsOnDirBnd[vel]) { // put coupling on rhs
                P2DiscCL::GetGradient( Grad_vel, GradRef[vel], T);
                tmp= vel<4 ? BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad_vel.quadP1( pr, absdet), tmp);
            }
        }
    }
}


void SetupRhs2_P2P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL& lset, double t)
// P2 / P1X FEs (X=extended) for vel/pr
{
    c->Clear( t);
    const Uint lvl=  c->GetLevel();
    const Uint pidx= c->RowIdx->GetIdx();
    const ExtIdxDescCL& Xidx= c->RowIdx->GetXidx();
    IdxT prNumb[4];
    bool IsOnDirBnd[10];
    Quad2CL<Point3DCL> Grad_vel, GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;
    InterfaceTetraCL cut;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        // collect some information about the edges and verts of the tetra
        // and save it in prNumb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            IsOnDirBnd[i]= BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) );
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for (int i=0; i<6; ++i)
            IsOnDirBnd[i+4]= BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) );

        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);
        cut.Init( *sit, lset.Phi, lset.GetBndData());
        const bool nocut= !cut.Intersects();

        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (IsOnDirBnd[vel]) { // put coupling on rhs
                P2DiscCL::GetGradient( Grad_vel, GradRef[vel], T);
                tmp= vel<4 ? BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr) {
                    const Point3DCL int_grad_vel_pr= Grad_vel.quadP1( pr, absdet);
                    c->Data[ prNumb[pr]]+= inner_prod( int_grad_vel_pr, tmp);

                    if (nocut) continue;
                    // compute the integrals
                    // \int_{T_2} grad v_vel p_pr dx  -  C \int_{T}grad v_vel p_pr dx,
                    // where C = H(Phi(x_pr)) \in {0,1} and T_2 = T \cap \Omega_2
                    LocalP2CL<Point3DCL> gradv_p;
                    gradv_p[pr]= Grad_vel[pr];
                    for (int vert=0; vert<4; ++vert)
                        if (vert!=pr)
                            gradv_p[EdgeByVert(pr,vert)+4]= 0.25*(Grad_vel[pr]+Grad_vel[vert]);

                    Point3DCL integral= cut.GetSign(pr)==1 ? -int_grad_vel_pr : Point3DCL();

                    for (int ch=0; ch<8; ++ch)
                    {
                        cut.ComputeCutForChild(ch);
                        integral+= cut.quad( gradv_p, absdet, true); // integrate on positive part
                    }

                    c->Data[ Xidx[prNumb[pr]]]+= inner_prod( integral, tmp);
                }
            }
        }
    }
}


void SetupRhs2_P2P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL&, const StokesBndDataCL& BndData, VecDescCL* c, double t)
{
    c->Clear( t);

    const Uint lvl         = c->GetLevel();
    const Uint pidx        = c->RowIdx->GetIdx();

    IdxT prNumb[4];
    bool IsOnDirBnd[10];

    Quad2CL<Point3DCL> Grad_vel, GradRef[10];

    SMatrixCL<3,3> T;

    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        // collect some information about the edges and verts of the tetra
        // and save it in prNumb and IsOnDirBnd
        for (int i=0; i<4; ++i) {
            IsOnDirBnd[i]= BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) );
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for (int i=0; i<6; ++i)
            IsOnDirBnd[i+4]= BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) );

        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        // b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (IsOnDirBnd[vel]) { // put coupling on rhs
                P2DiscCL::GetGradient( Grad_vel, GradRef[vel], T);
                tmp= vel<4 ? BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad_vel.quadP1D( pr, absdet), tmp);
            }
        }
    }
}


// -----------------------------------------------------------------------------
//                        Routines for SetupPrMass
// -----------------------------------------------------------------------------


void SetupPrMass_P0(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset)
{
    const IdxT num_unks_pr=  RowIdx.NumUnknowns();
    MatrixBuilderCL M_pr(&matM, num_unks_pr,  num_unks_pr);

    const Uint lvl= RowIdx.TriangLevel();

    SmoothedJumpCL nu_invers( 1./Coeff.mu(0), 1./Coeff.mu(1), Coeff.mu);
    Quad2CL<double> nu_inv;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint ls_lvl = ls.GetLevel();
    LocalP2CL<> locallset;
    const Uint pidx= RowIdx.GetIdx();

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin(lvl),
         send= MG.GetTriangTetraEnd(lvl); sit != send; ++sit) {
        const double absdet= sit->GetVolume()*6.;
        if (ls_lvl != lvl) {
            locallset.assign( *sit, ls);
            nu_inv.assign( locallset);
        }
        else
            nu_inv.assign( *sit, ls);
        nu_inv.apply( nu_invers);

        const IdxT prNumbTetra= sit->Unknowns( pidx);
        M_pr( prNumbTetra, prNumbTetra)= nu_inv.quad( absdet);

    }
    M_pr.Build();
}


void SetupPrMass_P1(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset)
{
    const IdxT num_unks_pr=  RowIdx.NumUnknowns();
    MatrixBuilderCL M_pr(&matM, num_unks_pr,  num_unks_pr);

    const Uint lvl= RowIdx.TriangLevel();
    IdxT prNumb[4];
    double coup[4][4];
    const double nu_inv_p= 1./Coeff.mu( 1.0),
                 nu_inv_n= 1./Coeff.mu( -1.0);
    double integralp, integraln;
    InterfaceTetraCL cut;
    LocalP2CL<> pipj[4][4];
    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
            pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
        }
        pipj[i][i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
                pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }

    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        loc_phi.assign( *sit, lset.Phi, lset.GetBndData());
        cut.Init( *sit, loc_phi);
        const bool nocut= !cut.Intersects();
        GetLocalNumbP1NoBnd( prNumb, *sit, RowIdx);
        if (nocut) { // nu is constant in tetra
            const double nu_inv= cut.GetSign( 0) == 1 ? nu_inv_p : nu_inv_n;
            // write values into matrix
            for(int i=0; i<4; ++i)
                for(int j=0; j<4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= nu_inv*P1DiscCL::GetMass( i, j)*absdet;
        }
        else { // nu is discontinuous in tetra
            for(int i=0; i<4; ++i) {
                for(int j=0; j<=i; ++j) {
                    // compute the integrals
                    // \int_{T_i} p_i p_j dx,    where T_i = T \cap \Omega_i, i=1,2
                    integralp= integraln= 0.;
                    for (int ch= 0; ch < 8; ++ch) {
                        cut.ComputeCutForChild( ch);
                        integralp+= cut.quad( pipj[i][j], absdet, true);  // integrate on positive part
                        integraln+= cut.quad( pipj[i][j], absdet, false); // integrate on negative part
                    }
                    coup[j][i]= integralp*nu_inv_p + integraln*nu_inv_n;
                    coup[i][j]= coup[j][i];
                }
            }

            // write values into matrix
            for(int i=0; i<4; ++i)
                for(int j= 0; j < 4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= coup[i][j];
        }
    }
    M_pr.Build();
}


void SetupPrMass_P1X(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset)
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    const IdxT num_unks_pr=  RowIdx.NumUnknowns();
    MatrixBuilderCL M_pr(&matM, num_unks_pr,  num_unks_pr);

    const Uint lvl= RowIdx.TriangLevel();
    IdxT prNumb[4];
    double coup[4][4], coupT2[4][4];

    const double nu_inv_p= 1./Coeff.mu( 1.0),
                 nu_inv_n= 1./Coeff.mu( -1.0);
    double integralp, integraln;
    InterfaceTetraCL cut;
    bool sign[4];

    // The 16 products of the P1-shape-functions
    LocalP2CL<> pipj[4][4];
    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
            pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
        }
        pipj[i][i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
                pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }

    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        loc_phi.assign( *sit, lset.Phi, lset.GetBndData());
        cut.Init( *sit, loc_phi);
        const bool nocut= !cut.Intersects();
        GetLocalNumbP1NoBnd( prNumb, *sit, RowIdx);
        if (nocut) { // nu is constant in tetra
            const double nu_inv= cut.GetSign( 0) == 1 ? nu_inv_p : nu_inv_n;
            // write values into matrix
            for(int i=0; i<4; ++i)
                for(int j=0; j<4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= nu_inv*P1DiscCL::GetMass( i, j)*absdet;
        }
        else { // extended basis functions have only support on tetra intersecting Gamma!
            for(int i=0; i<4; ++i) {
                sign[i]= cut.GetSign(i) == 1;
                for(int j=0; j<=i; ++j) {
                    // compute the integrals
                    // \int_{T_2} p_i p_j dx,    where T_2 = T \cap \Omega_2
                    integralp= integraln= 0.;
                    for (int ch= 0; ch < 8; ++ch) {
                        cut.ComputeCutForChild( ch);
                        integralp+= cut.quad( pipj[i][j], absdet, true);  // integrate on positive part
                        integraln+= cut.quad( pipj[i][j], absdet, false); // integrate on negative part
                    }
                    coup[j][i]= integralp*nu_inv_p + integraln*nu_inv_n;
                    coup[i][j]= coup[j][i];
                    coupT2[j][i]= integralp*nu_inv_p;
                    coupT2[i][j]= coupT2[j][i];
                }
            }

            // write values into matrix
            for(int i=0; i<4; ++i)
            {
                const IdxT xidx_i= Xidx[prNumb[i]];
                for(int j= 0; j < 4; ++j) {
                    M_pr( prNumb[i], prNumb[j])+= coup[i][j];
                    // tetra intersects Gamma => Xidx defined for all DoFs
                    const IdxT xidx_j= Xidx[prNumb[j]];
                    if (xidx_j!=NoIdx)
                        M_pr( prNumb[i], xidx_j)+= coupT2[i][j] - sign[j]*coup[i][j];
                    if (xidx_i!=NoIdx)
                        M_pr( xidx_i, prNumb[j])+= coupT2[i][j] - sign[i]*coup[i][j];
                    if (xidx_i!=NoIdx && xidx_j!=NoIdx && sign[i]==sign[j])
                        M_pr( xidx_i, xidx_j)+= sign[i] ? coup[i][j] - coupT2[i][j] : coupT2[i][j];
                }
            }
        }
    }

    M_pr.Build();
}


void SetupPrMass_P1D(const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& matM, IdxDescCL& RowIdx, const LevelsetP2CL& lset)
{
    const IdxT num_unks_pr=  RowIdx.NumUnknowns();
    MatrixBuilderCL M_pr(&matM, num_unks_pr,  num_unks_pr);

    const Uint lvl= RowIdx.TriangLevel();
    IdxT prNumb[4];

    SmoothedJumpCL nu_invers( 1./Coeff.mu(0), 1./Coeff.mu(1), Coeff.mu);
    Quad2CL<double> nu_inv;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint ls_lvl = ls.GetLevel();
    LocalP2CL<> locallset;

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin(lvl),
        send=MG.GetTriangTetraEnd(lvl); sit != send; ++sit) {
        const double absdet= sit->GetVolume()*6.;
        if (ls_lvl != lvl) {
            locallset.assign( *sit, ls);
            nu_inv.assign( locallset);
        }
        else
            nu_inv.assign( *sit, ls);
        nu_inv.apply( nu_invers);

        GetLocalNumbP1DNoBnd( prNumb, *sit, RowIdx);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                M_pr( prNumb[i], prNumb[j])+= nu_inv.quadP1D( i, j, absdet);
    }
    M_pr.Build();
}


// -----------------------------------------------------------------------------
//                        Routines for SetupPrStiff
// -----------------------------------------------------------------------------


void SetupPrStiff_P1( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset)
{
    MatrixBuilderCL A( &A_pr, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    const Uint lvl= RowIdx.TriangLevel();
    const Uint idx= RowIdx.GetIdx();
    SMatrixCL<3,4> G;
    double coup[4][4];
    double det, absdet, IntRhoInv;
    IdxT UnknownIdx[4];

    const double rho_inv_p= 1./Coeff.rho(1.),
                 rho_inv_n= 1./Coeff.rho(-1.);
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> loc_phi, ones( 1.);
    InterfaceTetraCL cut;

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit)
    {
        P1DiscCL::GetGradients( G,det,*sit);
        absdet= std::fabs( det);
        loc_phi.assign( *sit, ls);
        cut.Init( *sit, loc_phi);
        if (!cut.Intersects()) {
            IntRhoInv= absdet/6*(cut.GetSign( 0) == 1 ? rho_inv_p : rho_inv_n);
        } else {
            double Vol_p=0, Vol_n=0;
            for (int ch= 0; ch < 8; ++ch) {
                cut.ComputeCutForChild( ch);
                Vol_p+= cut.quad( ones, absdet, true);  // integrate on positive part
                Vol_n+= cut.quad( ones, absdet, false); // integrate on negative part
            }
            IntRhoInv= Vol_p*rho_inv_p + Vol_n*rho_inv_n;
        }
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )*IntRhoInv;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetVertex( i)->Unknowns( idx);
        }
        for(int i=0; i<4; ++i)    // assemble row i
            for(int j=0; j<4;++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i];
    }
    A.Build();
}

/// \todo: As in SetupPrMass_P1X, replace the smoothed density-function with integration
///        over the inner and outer part.
void SetupPrStiff_P1X( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset)
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    MatrixBuilderCL A( &A_pr, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    const Uint lvl= RowIdx.TriangLevel();
    const Uint idx= RowIdx.GetIdx();
    SMatrixCL<3,4> G;
    double coup[4][4], coupT2[4][4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];
    bool sign[4];
    InterfaceTetraCL cut;

    const double rho_inv_p= 1./Coeff.rho(1.),
                 rho_inv_n= 1./Coeff.rho(-1.);
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> locallset, ones( 1.);

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit)
    {
        locallset.assign( *sit, ls);
        cut.Init( *sit, locallset);
        const bool nocut= !cut.Intersects();

        P1DiscCL::GetGradients( G,det,*sit);
        absdet= std::fabs( det);
        double IntRhoInv, IntRhoInv_p;


        if (nocut) {
            IntRhoInv_p= cut.GetSign( 0) == 1 ? absdet/6*rho_inv_p : 0;
            IntRhoInv= absdet/6*(cut.GetSign( 0) == 1 ? rho_inv_p : rho_inv_n);
        } else {
            double Vol_p=0, Vol_n=0;
            for (int ch= 0; ch < 8; ++ch) {
                cut.ComputeCutForChild( ch);
                Vol_p+= cut.quad( ones, absdet, true);  // integrate on positive part
                Vol_n+= cut.quad( ones, absdet, false); // integrate on negative part
            }
            IntRhoInv_p= Vol_p*rho_inv_p;
            IntRhoInv=   Vol_p*rho_inv_p + Vol_n*rho_inv_n;
        }

        // compute local matrices
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )*IntRhoInv;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetVertex( i)->Unknowns( idx);
            if (nocut) continue; // extended basis functions have only support on tetra intersecting Gamma!

            sign[i]= cut.GetSign(i)==1;
            for(int j=0; j<=i; ++j) {
                // compute the integrals
                // \int_{T_2} grad_i grad_j dx,    where T_2 = T \cap \Omega_2
                coupT2[j][i]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )*IntRhoInv_p;
                coupT2[i][j]= coupT2[j][i];
            }
        }

        // write values into matrix
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<4; ++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i];
            if (nocut) continue; // extended basis functions have only support on tetra intersecting Gamma!

            const IdxT xidx_i= Xidx[UnknownIdx[i]];
            for(int j=0; j<4; ++j) // write values for extended basis functions
            {
                const IdxT xidx_j= Xidx[UnknownIdx[j]];
                if (xidx_j!=NoIdx)
                    A( UnknownIdx[i], xidx_j)+= coupT2[i][j] - sign[j]*coup[i][j];
                if (xidx_i!=NoIdx)
                    A( xidx_i, UnknownIdx[j])+= coupT2[i][j] - sign[i]*coup[i][j];
                if (xidx_i!=NoIdx && xidx_j!=NoIdx && sign[i]==sign[j])
                    A( xidx_i, xidx_j)+= sign[i] ? coup[i][j] - coupT2[i][j] : coupT2[i][j];
            }
        }
    }
    A.Build();
}


void SetupPrStiff_P1D( const MultiGridCL& MG, const TwoPhaseFlowCoeffCL& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx, const LevelsetP2CL& lset)
{
    MatrixBuilderCL A( &A_pr,RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    const Uint lvl= RowIdx.TriangLevel();
    const Uint idx= RowIdx.GetIdx();
    SMatrixCL<3,4> G;
    double coup[4][4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];

    SmoothedJumpCL rho_invers( 1./Coeff.rho(0), 1./Coeff.rho(1), Coeff.rho);
    Quad2CL<double> rho_inv;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint ls_lvl = ls.GetLevel();
    LocalP2CL<> locallset;

    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit)
    {
        if (ls_lvl != lvl) {
            locallset.assign( *sit, ls);
            rho_inv.assign( locallset);
        }
        else
            rho_inv.assign( *sit, ls);
        rho_inv.apply( rho_invers);

        P1DDiscCL::GetGradients( G,det,*sit);
        absdet= std::fabs( det);
        const double IntRhoInv= rho_inv.quad( absdet);
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )*IntRhoInv;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetFace( i)->Unknowns( idx);
        }
        for(int i=0; i<4; ++i)    // assemble row i
            for(int j=0; j<4;++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i];
    }
    A.Build();
}

// =============================================================================
//                        InstatStokes2PhaseP2P1CL
// =============================================================================

// Create numbering
// ----------------

void InstatStokes2PhaseP2P1CL::CreateNumberingVel( Uint level, MLIdxDescCL* idx, match_fun match, const LevelsetP2CL* lsetp)
{
    idx->CreateNumbering( level, MG_, BndData_.Vel, match, lsetp ? &(lsetp->Phi) : 0, lsetp ? &(lsetp->GetBndData()) : 0);
}


void InstatStokes2PhaseP2P1CL::CreateNumberingPr ( Uint level, MLIdxDescCL* idx, match_fun match, const LevelsetP2CL* lsetp)
{
    idx->CreateNumbering( level, MG_, BndData_.Pr, match, lsetp ? &(lsetp->Phi) : 0, lsetp ? &(lsetp->GetBndData()) : 0);
}


void InstatStokes2PhaseP2P1CL::SmoothVel( VelVecDescCL* v, int num, double tau)
{
    const VectorCL diag= A.Data.GetFinest().GetDiag();

    for (int i=0; i<num; ++i)
        v->Data-= tau*((A.Data.GetFinest()*v->Data)/diag);
}


void InstatStokes2PhaseP2P1CL::SetupPrMass( MLMatDescCL* matM, const LevelsetP2CL& lset) const
/// Needed for preconditioning of the Schur complement. Uses natural
/// boundary conditions for the pressure unknowns.
{
    MLMatrixCL::iterator itM = matM->Data.begin();
    MLIdxDescCL::iterator itIdx = matM->RowIdx->begin();
    for (size_t lvl=0; lvl < matM->Data.size(); ++lvl, ++itM, ++itIdx)
    {
        switch (GetPrFE())
        {
        case P0_FE:
            SetupPrMass_P0( MG_, Coeff_, *itM, *itIdx, lset); break;
        case P1_FE:
            SetupPrMass_P1( MG_, Coeff_, *itM, *itIdx, lset); break;
        case P1X_FE:
            SetupPrMass_P1X( MG_, Coeff_, *itM, *itIdx, lset); break;
        case P1D_FE:
            SetupPrMass_P1D( MG_, Coeff_, *itM, *itIdx, lset); break;
        default:
            throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupPrMass not implemented for this FE type");
        }
    }
}


void InstatStokes2PhaseP2P1CL::SetupPrStiff( MLMatDescCL* A_pr, const LevelsetP2CL& lset) const
/// Needed for preconditioning of the Schur complement. Uses natural
/// boundary conditions for the pressure unknowns.
{
    MLMatrixCL::iterator itM = A_pr->Data.begin();
    MLIdxDescCL::iterator itRowIdx = A_pr->RowIdx->begin();
    MLIdxDescCL::iterator itColIdx = A_pr->ColIdx->begin();
    for (size_t lvl=0; lvl < A_pr->Data.size(); ++lvl, ++itM, ++itRowIdx, ++itColIdx)
    {
        switch (GetPrFE())
        {
        case P1_FE:
            SetupPrStiff_P1( MG_, Coeff_, *itM, *itRowIdx, *itColIdx, lset); break;
        case P1X_FE:
            SetupPrStiff_P1X( MG_, Coeff_, *itM, *itRowIdx, *itColIdx, lset); break;
        case P1D_FE:
            SetupPrStiff_P1D( MG_, Coeff_, *itM, *itRowIdx, *itColIdx, lset); break;
        default:
            throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupPrStiff not implemented for this FE type");
        }
    }
}


void InstatStokes2PhaseP2P1CL::InitVel(VelVecDescCL* vec, instat_vector_fun_ptr LsgVel, double t0) const
{
    VectorCL& lsgvel= vec->Data;
    vec->t = t0;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

    for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>( MG_).GetTriangVertexBegin( lvl),
         send= const_cast<const MultiGridCL&>( MG_).GetTriangVertexEnd( lvl);
         sit != send; ++sit) {
        if (!BndData_.Vel.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel(sit->GetCoord(), t0));
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL sit= const_cast<const MultiGridCL&>( MG_).GetTriangEdgeBegin( lvl),
         send= const_cast<const MultiGridCL&>( MG_).GetTriangEdgeEnd( lvl);
         sit != send; ++sit) {
        if (!BndData_.Vel.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t0));
    }
}


/// \brief Raw data for "system 1", both for one phase and two phases.
///
/// scalar-valued mass-matrix, scalar-valued mu-Laplacian, genuinely tensor-valued part of the deformation tensor and the integrals of \f$\rho\phi_i\f$ for the gravitation as load-vector
/// \todo: Precise description
struct LocalSystem1DataCL
{
    double         M [10][10];
    double         A [10][10];
    SMatrixCL<3,3> Ak[10][10];

    double rho_phi[10];
};

/// \brief Setup of the local "system 1" on a tetra in a single phase.
class LocalSystem1OnePhase_P2CL
{
  private:
    double mu_;
    double rho_;

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    Quad2CL<> Ones;

  public:
    LocalSystem1OnePhase_P2CL (double muarg= 0., double rhoarg= 0.)
        : mu_( muarg), rho_( rhoarg), Ones( 1.)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void   mu  (double new_mu)        { mu_= new_mu; }
    double mu  ()               const { return mu_; }
    void   rho (double new_rho)       { rho_= new_rho; }
    double rho ()               const { return rho_; }

    void setup (const SMatrixCL<3,3>& T, double absdet, LocalSystem1DataCL& loc);
};

void LocalSystem1OnePhase_P2CL::setup (const SMatrixCL<3,3>& T, double absdet, LocalSystem1DataCL& loc)
{
    P2DiscCL::GetGradients( Grad, GradRef, T);
    for (int i=0; i<10; ++i)
    {
        loc.rho_phi[i]= rho()*Ones.quadP2( i, absdet);
        for (int j=0; j<=i; ++j)
        {
            // M: As we are not at the phase-boundary this is exact.
            loc.M[j][i]= rho()*P2DiscCL::GetMass( j, i)*absdet;

            // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k) = \int mu *\nabla\phi_i \outerprod \nabla\phi_j
            loc.Ak[j][i]= mu()*Quad2CL< SMatrixCL<3,3> >( outer_product( Grad[i], Grad[j])).quad( absdet);
            // dot-product of the gradients
            loc.A[j][i]= trace( loc.Ak[j][i]);
            if (i!=j) { // The local matrices coupM, coupA, coupAk are symmetric.
                loc.M[i][j]= loc.M[j][i];
                loc.A[i][j]= loc.A[j][i];
                assign_transpose( loc.Ak[i][j], loc.Ak[j][i]);
            }
        }
    }
}

/// \brief Setup of the local "system 1" on a tetra intersected by the dividing surface.
class LocalSystem1TwoPhase_P2CL
{
  private:
    const double mu_p, mu_n;
    const double rho_p, rho_n;

    LocalP1CL<Point3DCL> GradRefLP1[10], GradLP1[10];
    LocalP2CL<> p2;

    BaryCoordCL nodes_q2[Quad2CL<>::NumNodesC];
    BaryCoordCL nodes_q5[Quad5DataCL::NumNodesC];
    double absdets[48]; //there exists maximally 8*6=48 SubTetras
    Quad5CL<> q[48][10]; //there exists maximally 8*6=48 SubTetras
    Quad2CL<Point3DCL> qA[48][10]; //there exists maximally 8*6=48 SubTetras

    double intpos, intneg;
    SMatrixCL<3,3> cAkp, cAkn;

  public:
    LocalSystem1TwoPhase_P2CL (double mup, double mun, double rhop, double rhon)
        : mu_p( mup), mu_n( mun), rho_p( rhop), rho_n( rhon)
    { P2DiscCL::GetGradientsOnRef( GradRefLP1); }

    double mu  (int sign) const { return sign > 0 ? mu_p  : mu_n; }
    double rho (int sign) const { return sign > 0 ? rho_p : rho_n; }

    void setup (const SMatrixCL<3,3>& T, double absdet, InterfaceTetraCL& tetra, LocalSystem1DataCL& loc);
};

void LocalSystem1TwoPhase_P2CL::setup (const SMatrixCL<3,3>& T, double absdet, InterfaceTetraCL& tetra, LocalSystem1DataCL& loc)
{
    std::memset( loc.rho_phi, 0, 10*sizeof( double));
    P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);

    tetra.ComputeSubTets();
    for (Uint k=0; k<tetra.GetNumTetra(); ++k) {
        absdets[k]= absdet*VolFrac( tetra.GetTetra( k));
        Quad2CL<>::TransformNodes( tetra.GetTetra( k), nodes_q2);
        Quad5CL<>::TransformNodes( tetra.GetTetra( k), nodes_q5);
        for (int i= 0; i < 10; ++i) {
            // For M
            p2[i]= 1.; p2[i==0 ? 9 : i - 1]= 0.;
            q[k][i].assign( p2, nodes_q5);
            // For A
            qA[k][i].assign( GradLP1[i], nodes_q2);
            // \int rho*phi_i
            loc.rho_phi[i]+= (k < tetra.GetNumNegTetra() ? rho_n : rho_p )*q[k][i].quad( absdets[k]);
        }
    }
    for (int i= 0; i < 10; ++i) {
        for (int j= 0; j <= i; ++j) {
            intneg= intpos = 0.;
            std::memset( &cAkn, 0, sizeof( SMatrixCL<3,3>));
            std::memset( &cAkp, 0, sizeof( SMatrixCL<3,3>));
            for (Uint k=0; k<tetra.GetNumTetra(); ++k) {
                (k<tetra.GetNumNegTetra() ? intneg : intpos)+= Quad5CL<>( q[k][i]*q[k][j]).quad( absdets[k]);
                // A: kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                (k < tetra.GetNumNegTetra() ? cAkn : cAkp)+= Quad2CL< SMatrixCL<3,3> >( outer_product( qA[k][i], qA[k][j])).quad( absdets[k]);
            }
            loc.M[j][i]= rho_p*intpos + rho_n*intneg;
            loc.Ak[j][i]= mu_p*cAkp + mu_n*cAkn;
            // dot-product of the gradients
            loc.A[j][i]= trace( loc.Ak[j][i]);
            if (i != j) { // The local stiffness matrices coupM, coupA, coupAk are symmetric.
                loc.M[i][j]= loc.M[j][i];
                loc.A[i][j]= loc.A[j][i];
                assign_transpose( loc.Ak[i][j], loc.Ak[j][i]);
            }
        }
    }
}

/// \brief Accumulator to set up the matrices A, M and, if requested the right-hand side b and cplM, cplA for two-phase flow.
class System1Accumulator_P2CL : public TetraAccumulatorCL
{
  private:
    const TwoPhaseFlowCoeffCL& Coeff;
    const StokesBndDataCL& BndData;
    const LevelsetP2CL& lset;
    double t;

    IdxDescCL& RowIdx;
    MatrixCL& A;
    MatrixCL& M;
    VecDescCL* cplA;
    VecDescCL* cplM;
    VecDescCL* b;

    SparseMatBuilderCL<double, SMatrixCL<3,3> >* mA_;
    SparseMatBuilderCL<double, SDiagMatrixCL<3> >* mM_;

    LocalSystem1OnePhase_P2CL local_onephase; ///< used on tetras in a single phase
    LocalSystem1TwoPhase_P2CL local_twophase; ///< used on intersected tetras
    LocalSystem1DataCL loc; ///< Contains the memory, in which the local operators are set up; former coupM, coupA, coupAk, rho_phi.

    LocalNumbP2CL n; ///< global numbering of the P2-unknowns

    SMatrixCL<3,3> T;
    double det, absdet;
    InterfaceTetraCL tetra;

    Quad2CL<Point3DCL> rhs;
    Point3DCL loc_b[10], dirichlet_val[10]; ///< Used to transfer boundary-values from local_setup() update_global_system().

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and, if required, the Dirichlet-values needed to eliminate the boundary-dof from the global system.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    System1Accumulator_P2CL (const TwoPhaseFlowCoeffCL& Coeff, const StokesBndDataCL& BndData_,
        const LevelsetP2CL& ls, IdxDescCL& RowIdx_, MatrixCL& A_, MatrixCL& M_,
        VecDescCL* b_, VecDescCL* cplA_, VecDescCL* cplM_, double t);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);
};

System1Accumulator_P2CL::System1Accumulator_P2CL (const TwoPhaseFlowCoeffCL& Coeff_, const StokesBndDataCL& BndData_,
    const LevelsetP2CL& lset_arg, IdxDescCL& RowIdx_, MatrixCL& A_, MatrixCL& M_,
    VecDescCL* b_, VecDescCL* cplA_, VecDescCL* cplM_, double t_)
    : Coeff( Coeff_), BndData( BndData_), lset( lset_arg), t( t_),
      RowIdx( RowIdx_), A( A_), M( M_), cplA( cplA_), cplM( cplM_), b( b_),
      local_twophase( Coeff.mu( 1.0), Coeff.mu( -1.0), Coeff.rho( 1.0), Coeff.rho( -1.0))
{}

void System1Accumulator_P2CL::begin_accumulation ()
{
    std::cout << "entering SetupSystem1_P2CL::begin_accumulation ()" << std::endl;
    const size_t num_unks_vel= RowIdx.NumUnknowns();
    mA_= new SparseMatBuilderCL<double, SMatrixCL<3,3> >( &A, num_unks_vel, num_unks_vel);
    mM_= new SparseMatBuilderCL<double, SDiagMatrixCL<3> >( &M, num_unks_vel, num_unks_vel);
    if (b != 0) {
        b->Clear( t);
        cplM->Clear( t);
        cplA->Clear( t);
    }
}

void System1Accumulator_P2CL::finalize_accumulation ()
{
    mA_->Build();
    delete mA_;
    mM_->Build();
    delete mM_;
#ifndef _PAR
    std::cout << A.num_nonzeros() << " nonzeros in A, "
              << M.num_nonzeros() << " nonzeros in M! " << std::endl;
#endif
    std::cout << "leaving SetupSystem1_P2CL::finalize_accumulation ()" << std::endl;
}

void System1Accumulator_P2CL::visit (const TetraCL& tet)
{
    local_setup( tet);
    update_global_system();
}

void System1Accumulator_P2CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    rhs.assign( tet, Coeff.volforce, t);
    n.assign( tet, RowIdx, BndData.Vel);

    tetra.Init( tet, lset.Phi, lset.GetBndData());
    if (tetra.Intersects())
        local_twophase.setup( T, absdet, tetra, loc);
    else {
        local_onephase.mu(  local_twophase.mu(  tetra.GetSign( 0)));
        local_onephase.rho( local_twophase.rho( tetra.GetSign( 0)));
        local_onephase.setup( T, absdet, loc);
    }
    add_transpose_kronecker_id( loc.Ak, loc.A);

    if (b != 0) {
        for (int i= 0; i < 10; ++i) {
            if (!n.WithUnknowns( i)) {
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[i]).GetBndFun();
                dirichlet_val[i]= i<4 ? bf( tet.GetVertex( i)->GetCoord(), t)
                    : bf( GetBaryCenter( *tet.GetEdge( i-4)), t);
            }
            else
                loc_b[i]= rhs.quadP2( i, absdet) + loc.rho_phi[i]*Coeff.g;
        }
    }
}

void System1Accumulator_P2CL::update_global_system ()
{
    SparseMatBuilderCL<double, SMatrixCL<3,3> >& mA= *mA_;
    SparseMatBuilderCL<double, SDiagMatrixCL<3> >& mM= *mM_;

    for(int i=0; i<10; ++i)    // assemble row Numb[i]
        if (n.WithUnknowns( i)) { // dof i is not on a Dirichlet boundary
            for(int j=0; j<10; ++j) {
                if (n.WithUnknowns( j)) { // dof j is not on a Dirichlet boundary
                    mA( n.num[i],   n.num[j]  )+= loc.Ak[i][j];
                    mM( n.num[i],   n.num[j]  )+= SDiagMatrixCL<3>( loc.M[j][i]);
                }
                else if (b != 0) { // right-hand side for eliminated Dirichlet-values
                    add_to_global_vector( cplA->Data, -loc.Ak[i][j]*dirichlet_val[j], n.num[i]);
                    add_to_global_vector( cplM->Data, -loc.M[j][i] *dirichlet_val[j], n.num[i]);
                }
            }
            if (b != 0) // assemble the right-hand side
                add_to_global_vector( b->Data, loc_b[i], n.num[i]);
       }
}

void SetupSystem1_P2( const MultiGridCL& MG_, const TwoPhaseFlowCoeffCL& Coeff_, const StokesBndDataCL& BndData_, MatrixCL& A, MatrixCL& M,
                      VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, IdxDescCL& RowIdx, double t)
/// Set up matrices A, M and rhs b (depending on phase bnd)
{
    // TimerCL time;
    // time.Start();
    System1Accumulator_P2CL accu( Coeff_, BndData_, lset, RowIdx, A, M, b, cplA, cplM, t);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accus( MG_.GetTriangTetraBegin( RowIdx.TriangLevel()), MG_.GetTriangTetraEnd( RowIdx.TriangLevel()));
    // time.Stop();
    // std::cout << "setup: " << time.GetTime() << std::endl;
}


void SetupSystem1_P2R( const MultiGridCL& MG_, const TwoPhaseFlowCoeffCL& Coeff_, const StokesBndDataCL& BndData_, MatrixCL& A, MatrixCL& M,
                         VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, IdxDescCL& RowIdx, double t)
/// Set up matrices A, M and rhs b (depending on phase bnd)
{
    const IdxT num_unks_vel= RowIdx.NumUnknowns();
    const ExtIdxDescCL xidx= RowIdx.GetXidx();

    MatrixBuilderCL mA( &A, num_unks_vel, num_unks_vel),
                    mM( &M, num_unks_vel, num_unks_vel);
    if (b != 0)
    {
        b->Clear( t);
        cplM->Clear( t);
        cplA->Clear( t);
    }

    const Uint lvl = RowIdx.TriangLevel();

    LocalNumbP2CL n;
    IdxT num[14];
#ifndef _PAR
    std::cout << "entering SetupSystem1: " << num_unks_vel << " vels. ";
#endif

    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    LocalP1CL<Point3DCL> GradRefLP1[10], GradLP1[10], gradxLP1;
    LocalP2CL<Point3DCL> GradLP2[10];
    Quad2CL<double> Ones( 1.), kreuzterm;
    const double mu_p= Coeff_.mu( 1.0),
                 mu_n= Coeff_.mu( -1.0),
                 rho_p= Coeff_.rho( 1.0),
                 rho_n= Coeff_.rho( -1.0);

    SMatrixCL<3,3> T;

    double coupA[14][14], coupAk[14][14][3][3],
           coupM[14][14], rho_phi[14];
    double det, absdet, cAp, cAn, cAkp, cAkn, intHat_p, intHat_n;
    Point3DCL tmp, intRhs[14];
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> aij_n, aij_p, akreuz_n[3][3], akreuz_p[3][3], phi_i, ones( 1.);

    P2DiscCL::GetGradientsOnRef( GradRef);
    P2DiscCL::GetGradientsOnRef( GradRefLP1);

    InterfaceTetraCL patch;
    BaryCoordCL* nodes;
    LocalP2CL<> p1[4], p2[10]; // 4 P1 and 10 P2 basis functions
    LocalP2CL<> Fabs_p, Fabs_n; // enrichment function on pos./neg. part (to be interpreted as isoP2 function)
    LocalP2CL<> p1abs_p[4][8], p1abs_n[4][8]; // extended basis functions on pos./neg. part, resp., for each of the 8 regular children
    double intpos, intneg;
    for (int k=0; k<10; ++k) {
        p2[k][k]=1.;
        if (k<4)
            p1[k][k]=1.;
        else { // set corresponding edge value of p1 hat functions of corresponding vertices
            p1[VertOfEdge(k-4,0)][k]= 0.5;
            p1[VertOfEdge(k-4,1)][k]= 0.5;
        }
    }
    Quad5CL<> q[10][48], qx_p[4][48], qx_n[4][48]; // quadrature for basis functions (there exist maximally 8*6=48 SubTetras)
    LocalP2CL<> loc_phi;

    for (MultiGridCL::const_TriangTetraIteratorCL sit = MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        rhs.assign( *sit, Coeff_.volforce, t);

        // collect some information about the edges and verts of the tetra
        // and save it n.
        n.assign( *sit, RowIdx, BndData_.Vel);
        loc_phi.assign( *sit, ls);
        patch.Init( *sit, loc_phi);
        const bool nocut= !patch.Intersects();
        if (nocut) {
            const double mu_const= patch.GetSign( 0) == 1 ? mu_p : mu_n;
            const double rho_const= patch.GetSign( 0) == 1 ? rho_p : rho_n;

            P2DiscCL::GetGradients( Grad, GradRef, T);
            // compute all couplings between HatFunctions on edges and verts
            for (int i=0; i<10; ++i)
            {
                rho_phi[i]= rho_const*Ones.quadP2( i, absdet);
                intRhs[i]= rhs.quadP2(i, absdet);
                for (int j=0; j<=i; ++j)
                {
                    // dot-product of the gradients
                    const double cA= mu_const*Quad2CL<>(dot( Grad[i], Grad[j])).quad( absdet);
                    coupA[i][j]= cA;
                    coupA[j][i]= cA;
                    // kreuzterm
                    for (int k= 0; k < 3; ++k)
                        for (int l= 0; l < 3; ++l) {
                            // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                            for (size_t m=0; m<kreuzterm.size();  ++m)
                                kreuzterm[m]= Grad[i][m][l] * Grad[j][m][k];

                            coupAk[i][j][k][l]= mu_const*kreuzterm.quad( absdet);
                            coupAk[j][i][l][k]= coupAk[i][j][k][l];
                        }
                    // As we are not at the phase-boundary this is exact:
                    const double cM= rho_const*P2DiscCL::GetMass( j, i)*absdet;
                    coupM[i][j]= cM;
                    coupM[j][i]= cM;
                }
            }
        }
        else { // We are at the phase boundary.
            // compute all couplings between HatFunctions on edges and verts
            std::memset( coupA, 0, 14*14*sizeof( double));
            std::memset( coupAk, 0, 14*14*3*3*sizeof( double));
            std::memset( rho_phi, 0, 14*sizeof( double));
            P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);
            for (int i=0; i<10; ++i)
                GradLP2[i].assign( GradLP1[i]);

            double Vol_p= 0;
            P2RidgeDiscCL::GetExtBasisOnChildren(p1abs_p, p1abs_n, loc_phi);
            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeCutForChild( ch);
                Vol_p+= patch.quad( ones, absdet, true);

                LocalP2CL<Point3DCL> gradx_n[4], gradx_p[4]; // gradients of extended basis functions
                for (int i=0; i<4; ++i) { // init gradients of extended basis functions
                    P2DiscCL::GetFuncGradient( gradxLP1, p1abs_p[i][ch], GradLP1);
                    gradx_p[i].assign(gradxLP1);
                    P2DiscCL::GetFuncGradient( gradxLP1, p1abs_n[i][ch], GradLP1);
                    gradx_n[i].assign(gradxLP1);
                }
                for (int i=0; i<10; ++i) {
                    patch.quadBothParts( intHat_p, intHat_n, p2[i], absdet);
                    rho_phi[i]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_i
                }
                for (int i=0; i<4; ++i) {
                    intHat_p= patch.quad( p1abs_p[i][ch], absdet, true);
                    intHat_n= patch.quad( p1abs_n[i][ch], absdet, false);
                    rho_phi[i+10]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_abs
                }
                // compute coupA, coupAk
                for (int i=0; i<14; ++i) {
                    LocalP2CL<Point3DCL> &gradi_n= i<10 ? GradLP2[i] : gradx_n[i-10],
                                         &gradi_p= i<10 ? GradLP2[i] : gradx_p[i-10];
                    for (int j=0; j<=i; ++j) {
                        LocalP2CL<Point3DCL> &gradj_n= j<10 ? GradLP2[j] : gradx_n[j-10],
                                             &gradj_p= j<10 ? GradLP2[j] : gradx_p[j-10];
                        aij_n= dot( gradi_n, gradj_n);
                        aij_p= dot( gradi_p, gradj_p);
                        for (int k= 0; k < 3; ++k)
                            for (int l= 0; l < 3; ++l) {
                                // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m= 0; m < /* #Components of akreuz[k][l] */ 10;  ++m) {
                                    akreuz_n[k][l][m]= gradi_n[m][l] * gradj_n[m][k];
                                    akreuz_p[k][l][m]= gradi_p[m][l] * gradj_p[m][k];
                                }
                            }

                        // integrate aij on positive and negative part
                        cAp= patch.quad( aij_p, absdet, true);
                        cAn= patch.quad( aij_n, absdet, false);
                        const double cA= cAp*mu_p + cAn*mu_n;
                        coupA[i][j]+= cA;
                        for (int k= 0; k < 3; ++k)
                            for (int l= 0; l < 3; ++l) {
                                // integrate akreuz on positive and negative part
                                cAkp= patch.quad( akreuz_p[k][l], absdet, true);
                                cAkn= patch.quad( akreuz_n[k][l], absdet, false);
                                const double cAk= cAkp*mu_p + cAkn*mu_n;
                                coupAk[i][j][k][l]+= cAk;
                            }
                    }
                }
            } // child loop

            // compute coupM
            patch.ComputeSubTets();
            for (Uint k=0; k<patch.GetNumTetra(); ++k)
            { // init quadrature objects for basis functions
                nodes = Quad5CL<>::TransformNodes(patch.GetTetra(k));
                const int ch= patch.GetChildIdx(k);
                for (Uint j=0; j<10; ++j)  // standard FE
                    q[j][k].assign(p2[j], nodes);
                for (Uint j=0; j<4; ++j) { // extended FE
                    qx_p[j][k].assign(p1abs_p[j][ch], nodes);
                    qx_n[j][k].assign(p1abs_n[j][ch], nodes);
                }
                delete[] nodes;
            }
            for (int i=0; i<14; ++i) {
                Quad5CL<> *qi_n= i<10 ? q[i] : qx_n[i-10],
                          *qi_p= i<10 ? q[i] : qx_p[i-10];
                for (int j=0; j<=i; ++j) {
                    // M
                    intpos = 0.;
                    intneg = 0.;
                    Quad5CL<> *qj_n= j<10 ? q[j] : qx_n[j-10],
                              *qj_p= j<10 ? q[j] : qx_p[j-10];
                    for (Uint k=0; k<patch.GetNumTetra(); k++)
                        if (k<patch.GetNumNegTetra())
                            intneg += Quad5CL<>(qi_n[k]*qj_n[k]).quad(absdet*VolFrac(patch.GetTetra(k)));
                        else
                            intpos += Quad5CL<>(qi_p[k]*qj_p[k]).quad(absdet*VolFrac(patch.GetTetra(k)));
                    coupM[i][j]= rho_p*intpos + rho_n*intneg;
                }
                if (b != 0) {
                    if (i<10)
                        intRhs[i]= rhs.quadP2(i, absdet);
                    else {
                        intRhs[i]= Point3DCL();
                        for (Uint k=0; k<patch.GetNumTetra(); k++) {
                            nodes= Quad5CL<>::TransformNodes(patch.GetTetra(k));
                            Quad5CL<Point3DCL> rhs5( *sit, Coeff_.volforce, t, nodes);

                            if (k<patch.GetNumNegTetra())
                                intRhs[i] += Quad5CL<Point3DCL>(qi_n[k]*rhs5).quad(absdet*VolFrac(patch.GetTetra(k)));
                            else
                                intRhs[i] += Quad5CL<Point3DCL>(qi_p[k]*rhs5).quad(absdet*VolFrac(patch.GetTetra(k)));
                            delete[] nodes;
                        }
                    }
                }
            }
        }

        for (int i=0; i<14; ++i) {
            // collect local numbering
            num[i]= i<10 ? (n.WithUnknowns(i) ? n.num[i] : NoIdx) // standard FE part
                         : (nocut || !n.WithUnknowns(i-10) ? NoIdx : xidx[n.num[i-10]]);      // extended FE part
            // copy computed entries, as local stiffness matrices coupA, coupAk, coupM are symmetric
            for (int j=0; j<i; ++j) {
                coupA[j][i]= coupA[i][j];
                coupM[j][i]= coupM[i][j];
                for (int k= 0; k < 3; ++k)
                    for (int l= 0; l < 3; ++l)
                        coupAk[j][i][l][k]= coupAk[i][j][k][l];
            }
        }

        for(int i=0; i<14; ++i) {   // assemble row Numb[i]
            const IdxT numi= num[i];
            if (numi != NoIdx)  // dof i exists
            {
                for(int j=0; j<14; ++j)
                {
                    if (num[j] != NoIdx) // dof j exists
                    {
                        const IdxT numj= num[j];
                        mA( numi,   numj  )+= coupA[j][i];
                        mA( numi+1, numj+1)+= coupA[j][i];
                        mA( numi+2, numj+2)+= coupA[j][i];
                        for (int k=0; k<3; ++k)
                            for (int l=0; l<3; ++l)
                                mA( numi+k, numj+l)+= coupAk[i][j][k][l];
                        mM( numi,   numj  )+= coupM[j][i];
                        mM( numi+1, numj+1)+= coupM[j][i];
                        mM( numi+2, numj+2)+= coupM[j][i];
                    }
                    else if (b != 0 && j<10) // put coupling on rhs
                    /// \todo Interpolation of boundary data w.r.t. extended dofs not clear
                    {
                        typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                        bnd_val_fun bf= BndData_.Vel.GetBndSeg( n.bndnum[j]).GetBndFun();
                        tmp= j<4 ? bf( sit->GetVertex( j)->GetCoord(), t)
                                : bf( GetBaryCenter( *sit->GetEdge( j-4)), t);
                        const double cA= coupA[j][i],
                                    cM= coupM[j][i];
                        for (int k=0; k<3; ++k)
                        {
                            cplA->Data[numi+k]-= cA*tmp[k];
                            for (int l=0; l<3; ++l)
                                cplA->Data[numi+k]-= coupAk[i][j][k][l]*tmp[l];
                        }
                        cplM->Data[numi  ]-= cM*tmp[0];
                        cplM->Data[numi+1]-= cM*tmp[1];
                        cplM->Data[numi+2]-= cM*tmp[2];
                    }
                }
                if (b != 0)
                {
                    tmp= intRhs[i] + rho_phi[i]*Coeff_.g;
                    b->Data[numi  ]+= tmp[0];
                    b->Data[numi+1]+= tmp[1];
                    b->Data[numi+2]+= tmp[2];
                }
            }
        }
    }

    mA.Build();
    mM.Build();
#ifndef _PAR
    std::cout << A.num_nonzeros() << " nonzeros in A, "
              << M.num_nonzeros() << " nonzeros in M! " << std::endl;
#endif
}


void InstatStokes2PhaseP2P1CL::SetupSystem1( MLMatDescCL* A, MLMatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const
{
    MLMatrixCL::iterator itA = A->Data.begin();
    MLMatrixCL::iterator itM = M->Data.begin();
    MLIdxDescCL::iterator it = A->RowIdx->begin();
    for (size_t lvl=0; lvl < A->Data.size(); ++lvl, ++itA, ++itM, ++it)
        if (it->GetFE()==vecP2_FE)
            SetupSystem1_P2 ( MG_, Coeff_, BndData_, *itA, *itM, lvl == A->Data.size()-1 ? b : 0, cplA, cplM, lset, *it, t);
        else if (it->GetFE()==vecP2R_FE)
            SetupSystem1_P2R( MG_, Coeff_, BndData_, *itA, *itM, lvl == A->Data.size()-1 ? b : 0, cplA, cplM, lset, *it, t);
        else
            throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem1 not implemented for this FE type");
}


void SetupRhs1_P2( const MultiGridCL& MG_, const TwoPhaseFlowCoeffCL& Coeff_, const StokesBndDataCL& BndData_, VecDescCL* b, const LevelsetP2CL& lset, double t)
{
    const Uint lvl = b->GetLevel();

    b->Clear( t);

    LocalNumbP2CL n;
    SMatrixCL<3,3> T;

    Quad2CL<Point3DCL> rhs;
    Quad2CL<double> Ones( 1.);
    LocalP2CL<> phi_i;

    const double rho_p= Coeff_.rho( 1.0),
                 rho_n= Coeff_.rho( -1.0);
    double rho_phi[10];
    double det, absdet, intHat_p, intHat_n;
    Point3DCL tmp;

    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    InterfaceTetraCL tetra;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        rhs.assign( *sit, Coeff_.volforce, t);

        // collect some information about the edges and verts of the tetra
        // and save it n.
        n.assign( *sit, *b->RowIdx, BndData_.Vel);
        tetra.Init( *sit, lset.Phi, lset.GetBndData());
        const bool nocut= !tetra.Intersects();
        if (nocut) {
            const double rho_const= tetra.GetSign( 0) == 1 ? rho_p : rho_n;

            // compute all couplings between HatFunctions on edges and verts
            for (int i=0; i<10; ++i)
            {
                rho_phi[i]= rho_const*Ones.quadP2( i, absdet);
            }
        }
        else { // We are at the phase boundary.
            // compute all couplings between HatFunctions on edges and verts
            std::memset( rho_phi, 0, 10*sizeof( double));
            for (int ch= 0; ch < 8; ++ch) {
                tetra.ComputeCutForChild( ch);
                for (int i=0; i<10; ++i) {
                    // init phi_i =  i-th P2 hat function
                    phi_i[i]= 1.; phi_i[i==0 ? 9 : i-1]= 0.;
                    tetra.quadBothParts( intHat_p, intHat_n, phi_i, absdet);
                    rho_phi[i]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_i
                }
            }
        }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (n.WithUnknowns( i))  // vert/edge i is not on a Dirichlet boundary
            {
                tmp= rhs.quadP2( i, absdet) + rho_phi[i]*Coeff_.g;
                b->Data[n.num[i]  ]+= tmp[0];
                b->Data[n.num[i]+1]+= tmp[1];
                b->Data[n.num[i]+2]+= tmp[2];
            }
    }
}


void SetupRhs1_P2R( const MultiGridCL& MG_, const TwoPhaseFlowCoeffCL& Coeff_, const StokesBndDataCL& BndData_, VecDescCL* b, const LevelsetP2CL& lset, double t)
/// \todo proper implementation missing, yet
{
    throw DROPSErrCL("SetupRhs1_P2R(...) is buggy, aborting.");
    const Uint lvl = b->GetLevel();

    b->Clear( t);

    LocalNumbP2CL n;
    const IdxDescCL& RowIdx= *b->RowIdx;
    const ExtIdxDescCL xidx= RowIdx.GetXidx();

    Quad2CL<Point3DCL> rhs;
    Quad2CL<double> Ones( 1.), kreuzterm;
    const double rho_p= Coeff_.rho( 1.0),
                 rho_n= Coeff_.rho( -1.0);

    SMatrixCL<3,3> T;

    double rho_phi[14];
    double det, absdet, intHat_p, intHat_n;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> phi_i;

    InterfaceTetraCL patch;
    BaryCoordCL* nodes;
    LocalP2CL<> p1[4], p2[10]; // 4 P1 and 10 P2 basis functions
    LocalP2CL<> Fabs_p, Fabs_n; // enrichment function on pos./neg. part (to be interpreted as isoP2 function)
    Quad5CL<> qx_p[4][48], qx_n[4][48]; // quadrature for basis functions (there exist maximally 8*6=48 SubTetras)
    LocalP2CL<> loc_phi;

    for (int k=0; k<10; ++k) {
        p2[k][k]=1.;
        if (k<4)
            p1[k][k]=1.;
        else { // set corresponding edge value of p1 hat functions of corresponding vertices
            p1[VertOfEdge(k-4,0)][k]= 0.5;
            p1[VertOfEdge(k-4,1)][k]= 0.5;
        }
    }

    for (MultiGridCL::const_TriangTetraIteratorCL sit = MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        rhs.assign( *sit, Coeff_.volforce, t);
        
        // collect some information about the edges and verts of the tetra
        // and save it n.
        n.assign( *sit, RowIdx, BndData_.Vel);

        loc_phi.assign( *sit, ls);
        patch.Init( *sit, loc_phi);
        const bool nocut= !patch.Intersects();
        if (nocut) {
            const double rho_const= patch.GetSign( 0) == 1 ? rho_p : rho_n;

            for(int i=0; i<10; ++i)    // init rho_phi
                rho_phi[i]= rho_const*Ones.quadP2( i, absdet);
        }
        else { // We are at the phase boundary.
            for (int i= 0; i < 10; ++i) {
                // init enrichment function (to be interpreted as isoP2 function)
                Fabs_p[i]= patch.GetSign(i)== 1 ? 0 : -2*loc_phi[i];
                Fabs_n[i]= patch.GetSign(i)==-1 ? 0 :  2*loc_phi[i];
            }

            LocalP2CL<> p1abs_p[4][8], p1abs_n[4][8]; // extended basis functions on pos./neg. part, resp., for each of the 8 regular children
            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeCutForChild( ch);
                LocalP2CL<> extFabs_p, extFabs_n; // extension of enrichment function from child to parent
                // extend P1 values on child (as Fabs has to be interpreted as isoP2 function) to whole parent
                ExtendP1onChild( Fabs_p, ch, extFabs_p);
                ExtendP1onChild( Fabs_n, ch, extFabs_n);
                for (int i=0; i<4; ++i) { // init extended basis functions and its gradients
                    p1abs_p[i][ch]= p1[i]*extFabs_p;
                    p1abs_n[i][ch]= p1[i]*extFabs_n;
                }
                for (int i=0; i<10; ++i) {
                    // init phi_i =  i-th P2 hat function
                    phi_i[i]= 1.; phi_i[i==0 ? 9 : i-1]= 0.;
                    patch.quadBothParts( intHat_p, intHat_n, phi_i, absdet);
                    rho_phi[i]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_i
                }
                for (int i=0; i<4; ++i) {
                    intHat_p= patch.quad( p1abs_p[i][ch], absdet, true);
                    intHat_n= patch.quad( p1abs_n[i][ch], absdet, false);
                    rho_phi[i+10]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_abs
                }
            }
            patch.ComputeSubTets();
            for (Uint k=0; k<patch.GetNumTetra(); ++k)
            { // init quadrature objects for basis functions
                nodes = Quad5CL<>::TransformNodes(patch.GetTetra(k));
                const int ch= patch.GetChildIdx(k);
                for (Uint j=0; j<4; ++j) { // extended FE
                    qx_p[j][k].assign(p1abs_p[j][ch], nodes);
                    qx_n[j][k].assign(p1abs_n[j][ch], nodes);
                }
                delete[] nodes;
            }
        }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (n.WithUnknowns( i))  // vert/edge i is not on a Dirichlet boundary
            {
                tmp= rhs.quadP2( i, absdet) + rho_phi[i]*Coeff_.g;
                b->Data[n.num[i]  ]+= tmp[0];
                b->Data[n.num[i]+1]+= tmp[1];
                b->Data[n.num[i]+2]+= tmp[2];

                if (i<4 && !nocut) // extended dof
                {
                    Point3DCL intRhs;
                    for (Uint k=0; k<patch.GetNumTetra(); k++) {
                        nodes= Quad5CL<Point3DCL>::TransformNodes(patch.GetTetra(k));
                        Quad5CL<Point3DCL> rhs5( *sit, Coeff_.volforce, t, nodes);
                        if (k<patch.GetNumNegTetra())
                            intRhs += Quad5CL<Point3DCL>(qx_n[i][k]*rhs5).quad(absdet*VolFrac(patch.GetTetra(k)));
                        else
                            intRhs += Quad5CL<Point3DCL>(qx_p[i][k]*rhs5).quad(absdet*VolFrac(patch.GetTetra(k)));
                        delete[] nodes;
                    }

                    tmp= intRhs + rho_phi[i+10]*Coeff_.g;
                    const IdxT xnum= xidx[n.num[i]];
                    b->Data[xnum  ]+= tmp[0];
                    b->Data[xnum+1]+= tmp[1];
                    b->Data[xnum+2]+= tmp[2];
                }
            }
    }
}


void InstatStokes2PhaseP2P1CL::SetupRhs1( VecDescCL* b, const LevelsetP2CL& lset, double t) const
/// Set up rhs b (depending on phase bnd)
{
    const FiniteElementT fe= b->RowIdx->GetFE();
    if (fe==vecP2_FE)
        SetupRhs1_P2 ( MG_, Coeff_, BndData_, b, lset, t);
    else if (fe==vecP2R_FE)
        SetupRhs1_P2R( MG_, Coeff_, BndData_, b, lset, t);
    else
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs1 not implemented for this FE type");
}


void SetupLB_P2( const MultiGridCL& MG_, const TwoPhaseFlowCoeffCL& Coeff_, const StokesBndDataCL& BndData_, MatrixCL& A, VelVecDescCL* cplA, const LevelsetP2CL& lset, IdxDescCL& RowIdx, double t)
/// Set up the Laplace-Beltrami-matrix
{
    const IdxT num_unks_vel= RowIdx.NumUnknowns();

    MatrixBuilderCL mA( &A, num_unks_vel, num_unks_vel);
    if (cplA != 0) cplA->Clear( t);

    const Uint lvl = RowIdx.TriangLevel();

    LocalNumbP2CL n;
#ifndef _PAR
    std::cout << "entering SetupLB: " << num_unks_vel << " vels. ";
#endif
    LocalP1CL<Point3DCL> GradRefLP1[10], GradLP1[10];
    LocalP2CL<Point3DCL> GradLP2[10];

    Quad5_2DCL<Point3DCL> surfGrad[10];
    Quad5_2DCL<> surfTension, LB;

    SMatrixCL<3,3> T;

    double coupA[10][10];
    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> aij;

    P2DiscCL::GetGradientsOnRef( GradRefLP1);

    InterfaceTriangleCL triangle;
    LocalP2CL<> loc_phi;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        loc_phi.assign( *sit, lset.Phi, lset.GetBndData());
        triangle.Init( *sit, loc_phi);
        if (triangle.Intersects()) { // We are at the phase boundary.
            n.assign( *sit, RowIdx, BndData_.Vel);
            GetTrafoTr( T, det, *sit);
            absdet= std::fabs( det);
            P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);
            std::memset( coupA, 0, 10*10*sizeof( double));
            for (int i= 0; i < 10; ++i)
                GradLP2[i].assign( GradLP1[i]);

            for (int ch= 0; ch < 8; ++ch) {
                triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++ tri) {
                    surfTension= Coeff_.SurfTens;
                    for (int i= 0; i < 10; ++i) {
                        surfGrad[i].assign( GradLP1[i], &triangle.GetBary( tri));
                        surfGrad[i].apply( triangle, &InterfaceTriangleCL::ApplyProj);
                    }
                    for (int i=0; i<10; ++i)
                        for (int j=0; j<=i; ++j) {
                            // Laplace-Beltrami... Stabilization
                            LB= surfTension*dot( surfGrad[i], surfGrad[j]);
                            const double cLB= LB.quad( triangle.GetAbsDet( tri));
                            coupA[j][i]+= cLB;
                            coupA[i][j]= coupA[j][i];
                        }
                }
            }

            for(int i=0; i<10; ++i)    // assemble row Numb[i]
                if (n.WithUnknowns( i))  // vert/edge i is not on a Dirichlet boundary
                {
                    for(int j=0; j<10; ++j)
                    {
                        if (n.WithUnknowns( j)) // vert/edge j is not on a Dirichlet boundary
                        {
                            mA( n.num[i],   n.num[j]  )+= coupA[j][i];
                            mA( n.num[i]+1, n.num[j]+1)+= coupA[j][i];
                            mA( n.num[i]+2, n.num[j]+2)+= coupA[j][i];
                        }
                        else if (cplA != 0)  // put coupling on rhs
                        {
                            typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                            bnd_val_fun bf= BndData_.Vel.GetBndSeg( n.bndnum[j]).GetBndFun();
                            tmp= j<4 ? bf( sit->GetVertex( j)->GetCoord(), t)
                                    : bf( GetBaryCenter( *sit->GetEdge( j-4)), t);
                            const double cA= coupA[j][i];
                            for (int k=0; k<3; ++k)
                                cplA->Data[n.num[i]+k]-= cA*tmp[k];
                        }
                    }
                }
        }
    }

    mA.Build();
#ifndef _PAR
    std::cout << A.num_nonzeros() << " nonzeros in A_LB" << std::endl;
#endif
}


void InstatStokes2PhaseP2P1CL::SetupLB (MLMatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const
{
    MLMatrixCL::iterator  itA = A->Data.begin();
    MLIdxDescCL::iterator it  = A->RowIdx->begin();
    for (size_t lvl=0; lvl < A->RowIdx->size(); ++lvl, ++itA, ++it)
        SetupLB_P2( MG_,  Coeff_, BndData_, *itA, lvl == A->Data.size()-1 ? cplA : 0, lset, *it, t);
}


void InstatStokes2PhaseP2P1CL::SetupSystem2( MLMatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const
// Set up matrix B and rhs c
{
    MLMatrixCL::iterator     itB   = B->Data.begin();
    MLIdxDescCL::iterator    itRow = B->RowIdx->begin();
    MLIdxDescCL::iterator    itCol = B->ColIdx->begin();
    if ( B->RowIdx->size() == 1 || B->ColIdx->size() == 1)
    { // setup B only on finest level, if row or column index has only 1 level
        itCol = B->ColIdx->GetFinestIter();
        itRow = B->RowIdx->GetFinestIter();
        itB   = B->Data.GetFinestIter();
    }
    for (; itB!=B->Data.end() && itRow!=B->RowIdx->end() && itCol!=B->ColIdx->end(); ++itB, ++itRow, ++itCol)
    {
#ifndef _PAR
        std::cout << "entering SetupSystem2: " << itRow->NumUnknowns() << " prs, " << itCol->NumUnknowns() << " vels. ";
#endif
        VecDescCL* rhsPtr= itB==B->Data.GetFinestIter() ? c : 0; // setup rhs only on finest level
        if (itCol->GetFE()==vecP2_FE)
            switch (GetPrFE())
            {
                case P0_FE:
                    SetupSystem2_P2P0 ( MG_, Coeff_, BndData_, &(*itB), rhsPtr, &(*itRow), &(*itCol), t); break;
                case P1_FE:
                    SetupSystem2_P2P1 ( MG_, Coeff_, BndData_, &(*itB), rhsPtr, &(*itRow), &(*itCol), t); break;
                case P1X_FE:
                    SetupSystem2_P2P1X( MG_, Coeff_, BndData_, &(*itB), rhsPtr, lset, &(*itRow), &(*itCol), t); break;
                case P1D_FE:
                    SetupSystem2_P2P1D( MG_, Coeff_, BndData_, &(*itB), rhsPtr, &(*itRow), &(*itCol), t); break;
                default:
                    throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2 not implemented for this FE type");
            }
        else if (itCol->GetFE()==vecP2R_FE)
            switch (GetPrFE())
            {
                case P1_FE:
                    SetupSystem2_P2RP1 ( MG_, Coeff_, BndData_, &(*itB), rhsPtr, lset, &(*itRow), &(*itCol), t); break;
                case P1X_FE:
                    SetupSystem2_P2RP1X( MG_, Coeff_, BndData_, &(*itB), rhsPtr, lset, &(*itRow), &(*itCol), t); break;
                default:
                    throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2 not implemented for this FE type");
            }
        else
            throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2 not implemented for this FE type");
#ifndef _PAR
        std::cout << itB->num_nonzeros() << " nonzeros in B!" << std::endl;
#endif
    }
}


void InstatStokes2PhaseP2P1CL::SetupRhs2( VecDescCL* c, const LevelsetP2CL& lset, double t) const
// Set up rhs c
{
    if (vel_idx.GetFinest().GetFE()==vecP2_FE)
        switch (GetPrFE())
        {
          case P0_FE:
            SetupRhs2_P2P0( MG_, Coeff_, BndData_, c, t); break;
          case P1_FE:
            SetupRhs2_P2P1( MG_, Coeff_, BndData_, c, t); break;
          case P1X_FE:
            SetupRhs2_P2P1X( MG_, Coeff_, BndData_, c, lset, t); break;
          case P1D_FE:
            SetupRhs2_P2P1D( MG_, Coeff_, BndData_, c, t); break;
          default:
            throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2 not implemented for this FE type");
        }
    else if (vel_idx.GetFinest().GetFE()==vecP2R_FE)
        switch (GetPrFE())
        {
//          case P1_FE:
//            SetupRhs2_P2RP1( MG_, Coeff_, BndData_, c, t); break;
//          case P1X_FE:
//            SetupRhs2_P2RP1X( MG_, Coeff_, BndData_, c, lset, t); break;
          default:
            throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2 not implemented for this FE type");
        }
    else
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2 not implemented for this FE type");
}


void InstatStokes2PhaseP2P1CL::SetupBdotv (VecDescCL* Bdotv, const VelVecDescCL* vel,
    const LevelsetP2CL& lset, double t) const
{
    Bdotv->Clear( t);
    const Uint lvl= Bdotv->GetLevel();
    IdxT prNumb[4];

    LocalP1CL<Point3DCL>  Grad[10], GradRef[10]; // Gradient of p2-hat-functions
    P2DiscCL::GetGradientsOnRef( GradRef);
    LocalP2CL<Point3DCL>  lp2Grad;
    LocalP2CL<Point3DCL> loc_u;
    LocalP2CL<>  divu;
    Quad5_2DCL<Point3DCL> qGrad;
    Quad5_2DCL<Point3DCL> n, qu;
    Quad5_2DCL<> qdivu, q1, q2;
    LocalP1CL<double> lp1[4]; // p1-hat-functions
    for (int i= 0; i < 4; ++i) lp1[i][i]= 1.;

    SMatrixCL<3,3> T;
    double det;
    InterfaceTriangleCL cut;
    const ExtIdxDescCL& p_xidx= Bdotv->RowIdx->GetXidx();

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {
        cut.Init( *sit, lset.Phi, lset.GetBndData());
        if (!cut.Intersects()) continue;

        GetLocalNumbP1NoBnd( prNumb, *sit, *Bdotv->RowIdx);
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        loc_u.assign( *sit, *vel, GetBndData().Vel);
        divu= 0.;
        for (int i= 0; i < 10; ++i) {
            lp2Grad.assign( Grad[i]);
            divu+= dot( LocalP2CL<Point3DCL>( loc_u[i]), lp2Grad);
        }
        for (int ch= 0; ch < 8; ++ch) {
            cut.ComputeForChild( ch);
            for (int t= 0; t < cut.GetNumTriangles(); ++t) {
                const BaryCoordCL* const p( &cut.GetBary( t));
                qu.assign( loc_u, p);
                qdivu.assign( divu, p);
                n= Point3DCL();
                for (int v= 0; v < 10; ++v) {
                    qGrad.assign( Grad[v], p);
                    n+= cut.GetPhi( v)*qGrad;
                }
                for (int i= 0; i < Quad5_2DDataCL::NumNodesC; ++i)
                    if (n[i].norm()>1e-8) n[i]/= n[i].norm();
                q1= dot( n, qu)*qdivu;
                for(int pr= 0; pr < 4; ++pr) {
                    const IdxT xidx( p_xidx[prNumb[pr]]);
                    if (xidx == NoIdx) continue;
                    q2.assign( lp1[pr], p);
                    q2*= q1;
                    // n is the outer normal of {lset <= 0}; we need the outer normal of supp(pr-hat-function)\cap \Gamma).
                    Bdotv->Data[xidx]-= (cut.GetSign( pr) > 0 ? -1. : 1.)*q2.quad( cut.GetAbsDet( t));
                }
            }
        }
    }
}


void InstatStokes2PhaseP2P1CL::SetIdx()
{
    MLIdxDescCL* vidx= &vel_idx;
    MLIdxDescCL* pidx= &pr_idx;

    b.SetIdx   ( vidx);
    c.SetIdx   ( pidx);

    A.SetIdx   ( vidx, vidx);
    B.SetIdx   ( pidx, vidx);
    prM.SetIdx ( pidx, pidx);
    prA.SetIdx ( pidx, pidx);
    M.SetIdx   ( vidx, vidx);
}


void InstatStokes2PhaseP2P1CL::SetNumVelLvl( size_t n)
{
#ifdef _PAR
    if (n>1)
        throw DROPSErrCL("Multilevel not implemented in parallel DROPS yet, sorry");
#endif
    match_fun match= MG_.GetBnd().GetMatchFun();
    const double bound = vel_idx.GetFinest().GetXidx().GetBound();
    vel_idx.resize( n, GetVelFE(), BndData_.Vel, match, bound);
    A.Data.resize   (vel_idx.size());
    M.Data.resize   (vel_idx.size());
}


void InstatStokes2PhaseP2P1CL::SetNumPrLvl( size_t n)
{
#ifdef _PAR
    if (n>1)
        throw DROPSErrCL("Multilevel not implemented in parallel DROPS yet, sorry");
#endif
    match_fun match= MG_.GetBnd().GetMatchFun();
    const double bound = pr_idx.GetFinest().GetXidx().GetBound();
    pr_idx.resize( n, GetPrFE(),  BndData_.Pr, match, bound);
    B.Data.resize   (pr_idx.size());
    prM.Data.resize (pr_idx.size());
    prA.Data.resize (pr_idx.size());
}


void InstatStokes2PhaseP2P1CL::GetPrOnPart( VecDescCL& p_part, const LevelsetP2CL& lset, bool posPart)
{
    const Uint lvl= p.RowIdx->TriangLevel(),
        idxnum= p.RowIdx->GetIdx();
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const MultiGridCL& mg= this->GetMG();
    const ExtIdxDescCL& Xidx= this->GetXidx();

    p_part.SetIdx( p.RowIdx);
    VectorCL& pp= p_part.Data;
    pp= p.Data;

    // add extended pressure
    for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
        end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (Xidx[nr]==NoIdx) continue;

        const bool is_pos= InterfacePatchCL::Sign( ls.val( *it))==1;
        if (posPart==is_pos) continue; // extended hat function ==0 on this part
        if (posPart)
            pp[nr]+= p.Data[Xidx[nr]];
        else
            pp[nr]-= p.Data[Xidx[nr]];
    }
}


double InstatStokes2PhaseP2P1CL::GetCFLTimeRestriction( LevelsetP2CL& lset)
{
    const Uint lvl= p.RowIdx->TriangLevel();
    const MultiGridCL& mg= this->GetMG();
    const_DiscVelSolCL vel= GetVelSolution();
    LocalP2CL<Point3DCL> velLoc;
    LocalNumbP2CL curvNumb;
    VecDescCL curv( &vel_idx);
    lset.AccumulateBndIntegral( curv);

    double convMax= -1, viscMax= -1., gravMax= -1, stMax= -1;
    const double rho_min= std::min( Coeff_.rho(-1.), Coeff_.rho(1.)),
            nu_max= std::max( Coeff_.mu(-1.)/Coeff_.rho(-1.), Coeff_.mu(1.)/Coeff_.rho(1.));

    for( MultiGridCL::const_TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl),
        end= mg.GetTriangTetraEnd(lvl); it != end; ++it)
    {
        velLoc.assign( *it, vel);
        curvNumb.assign( *it, *curv.RowIdx, BndData_.Vel);

        // compute average curvature
        double tauKappa= 0., visc= 0.,
            h_min=1e99;
        for (int i=0; i<10; ++i)
            if (curvNumb.WithUnknowns(i))
                tauKappa+= curv.Data[curvNumb.num[i]];

        tauKappa/= it->GetVolume();

        for ( TetraCL::const_EdgePIterator ed= it->GetEdgesBegin(), eend= it->GetEdgesEnd(); ed!=eend; ++ed)
        {
            const Point3DCL dir= (*ed)->GetVertex(0)->GetCoord() - (*ed)->GetVertex(1)->GetCoord();
            const double length= norm(dir);
            if (length < h_min) h_min= length;
            visc+= 1./length/length;

            const double grav= std::sqrt(std::abs( inner_prod( Coeff_.g, dir)/length/length));
            if (grav > gravMax) gravMax= grav;

            for (int i=0; i<10; ++i) {
                const double conv= std::abs( inner_prod( velLoc[i], dir)/length/length);
                if (conv > convMax) convMax= conv;
            }
        }

        visc*= nu_max;
        if (visc > viscMax) viscMax= visc;

        const double st= std::sqrt(std::abs(tauKappa)/rho_min/h_min/h_min);
        if (st > stMax) stMax= st;
    }

    const double dtMax= 2./(convMax + viscMax + std::sqrt( (convMax+viscMax)*(convMax+viscMax) + 4*gravMax*gravMax + 4*stMax*stMax));

    std::cout << "CFL factors: conv= " << convMax << "\tvisc = " << viscMax << "\tgrav = " << gravMax << "\tst = " << stMax
        << " \n\t dt < " << dtMax << std::endl;

    return dtMax;
}



P1XRepairCL::P1XRepairCL (MultiGridCL& mg, VecDescCL& p)
    : UsesXFEM_( p.RowIdx->IsExtended()), mg_( mg), idx_( P1_FE), ext_( &idx_), p_( p)
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    const ExtIdxDescCL& extidx= p_.RowIdx->GetXidx();
//     size_t extbegin( extidx.GetNumUnknownsStdFE());
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p.RowIdx->GetIdx());

    idx_.CreateNumbering( p.RowIdx->TriangLevel(), mg);
    // Save the extended unknown-values.
    ext_.SetIdx( &idx_);

    // Attach the extended index to the vertex, so that it survives grid modifications.
    // We assume that all vertices have p-unknowns (like e. g. the pressure).
    DROPS_FOR_TRIANG_VERTEX( mg, p.RowIdx->TriangLevel(), it) {
        if (!it->Unknowns.Exist( pidx)) continue;
        if ( (extunknown= extidx[it->Unknowns( pidx)]) != NoIdx ) {
            ext_.Data[ it->Unknowns( repairidx)]= p.Data[extunknown];
        }
        else{
            it->Unknowns( repairidx)= NoIdx;
        }
    }
}

void P1XRepairCL::operator() ()
{
    if (!UsesXFEM_) return; // Do nothing, if we do not use P1X-elements.

    // We assume that the caller has repaired p as a P1-FE-variable.
    // Thus we only repair the extended part of p.
    size_t ci= 0;
    size_t extunknown;
    Uint repairidx( idx_.GetIdx()),
         pidx( p_.RowIdx->GetIdx());
    const ExtIdxDescCL& extidx= p_.RowIdx->GetXidx();
    // We assume that all vertices in p's level hold a p1-value. Thus, it->Unknowns.Exist()
    // can be spared.
    DROPS_FOR_TRIANG_VERTEX( mg_, p_.RowIdx->TriangLevel(), it) {
        if (!it->Unknowns.Exist( pidx)) continue;
        if ( ((extunknown= extidx[it->Unknowns( pidx)]) != NoIdx) && it->Unknowns.Exist( repairidx) ) {
            p_.Data[extunknown]= ext_.Data[it->Unknowns( repairidx)];
            ++ci;
        }
    }
//     std::cout << "P1XRepairCL::(): #P1-unknowns: " << extidx.GetNumUnknownsStdFE()
//               << "\t#copied extended-dof: " << ci << '\n';
}


void SetupMassDiag_P1(const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const BndCondCL& bnd)
{
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL Numb;

    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        Numb.assign( *sit, RowIdx, bnd);
        for(int i=0; i<4; ++i)
            if (Numb.WithUnknowns( i))
                M[Numb.num[i]]+= P1DiscCL::GetMass( i, i)*absdet;
    }
}

void SetupMassDiag_P1X (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const VecDescCL& lset,
                        const BndDataCL<>& lsetbnd, const BndCondCL& bnd)
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL Numb;
    double coup[4], coupT2[4];

    double integralp;
    InterfaceTetraCL cut;
    bool sign[4];

    // The 4 squares of the P1-shape-functions
    LocalP2CL<> pi2[4];
    for(int i= 0; i < 4; ++i) {
        pi2[i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
            pi2[i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }

    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        loc_phi.assign( *sit, lset, lsetbnd);
        cut.Init( *sit, loc_phi);
        const bool nocut= !cut.Intersects();
        Numb.assign( *sit, RowIdx, bnd);
        if (nocut) {
            for(int i= 0; i < 4; ++i)
                if ( Numb.WithUnknowns( i))
                    M[Numb.num[i]]+= P1DiscCL::GetMass( i, i)*absdet;
        }
        else { // extended basis functions have only support on tetra intersecting Gamma!
            for(int i=0; i<4; ++i) {
                sign[i]= cut.GetSign(i) == 1;
                // compute the integrals
                // \int_{T_2} p_i^2 dx,    where T_2 = T \cap \Omega_2
                integralp= 0.;
                for (int ch= 0; ch < 8; ++ch) {
                    cut.ComputeCutForChild( ch);
                    integralp+= cut.quad( pi2[i], absdet, true);  // integrate on positive part
                }
                coup[i]= P1DiscCL::GetMass( i, i)*absdet;
                coupT2[i]= integralp;
            }

            // write values into matrix
            for(int i=0; i<4; ++i) {
                if (!Numb.WithUnknowns( i)) continue;
                M[Numb.num[i]]+= coup[i];
                const IdxT xidx_i= Xidx[Numb.num[i]];
                if (xidx_i!=NoIdx)
                    M[xidx_i]+= coupT2[i]*(1 - 2*sign[i]) + sign[i]*coup[i];
            }
        }
    }
}

void SetupMassDiag_vecP2(const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const BndCondCL& bnd)
{
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP2CL Numb;

    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        Numb.assign( *sit, RowIdx, bnd);
        for(int i=0; i<10; ++i)
            if (Numb.WithUnknowns( i)) {
            	const double contrib= P2DiscCL::GetMass( i, i)*absdet;
                M[Numb.num[i]  ]+= contrib;
                M[Numb.num[i]+1]+= contrib;
                M[Numb.num[i]+2]+= contrib;
            }
    }
}

void SetupMassDiag (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const BndCondCL& bnd, const VecDescCL* lsetp, const BndDataCL<>* lsetbnd)
{
    switch(RowIdx.GetFE())
    {
    case P1_FE:
        SetupMassDiag_P1( MG, M, RowIdx, bnd); break;
    case P1X_FE:
        SetupMassDiag_P1X( MG, M, RowIdx, *lsetp, *lsetbnd, bnd); break;
    case vecP2_FE:
        SetupMassDiag_vecP2( MG, M, RowIdx, bnd); break;
    default:
        throw DROPSErrCL("SetupMassDiag not implemented for this FE type");
    }
}



void SetupLumpedMass_P1(const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const BndCondCL& bnd)
{
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL Numb;

    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        Numb.assign( *sit, RowIdx, bnd);
        for(int i=0; i<4; ++i)
            if (Numb.WithUnknowns( i))
                M[Numb.num[i]]+= P1DiscCL::GetLumpedMass( i)*absdet;
    }
}

void SetupLumpedMass_P1X (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const VecDescCL& lset, const BndDataCL<>& lsetbnd, const BndCondCL& bnd)
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL Numb;
    double coup[4], coupT2[4];

    double integralp;
    InterfaceTetraCL cut;
    bool sign[4];

    // The 4 P1-shape-functions
    LocalP2CL<> pi[4];
    for(int i= 0; i < 4; ++i) {
        pi[i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
            pi[i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.5;
    }

    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        loc_phi.assign( *sit, lset, lsetbnd);
        cut.Init( *sit, loc_phi);
        const bool nocut= !cut.Intersects();
        Numb.assign( *sit, RowIdx, bnd);
        if (nocut) {
            for(int i= 0; i < 4; ++i)
                if ( Numb.WithUnknowns( i))
                    M[Numb.num[i]]+= P1DiscCL::GetMass( i, i)*absdet;
        }
        else { // extended basis functions have only support on tetra intersecting Gamma!
            for(int i=0; i<4; ++i) {
                sign[i]= cut.GetSign(i) == 1;
                // compute the integrals
                // \int_{T_2} p_i dx,    where T_2 = T \cap \Omega_2
                integralp= 0.;
                for (int ch= 0; ch < 8; ++ch) {
                    cut.ComputeCutForChild( ch);
                    integralp+= cut.quad( pi[i], absdet, true);  // integrate on positive part
                }
                coup[i]= P1DiscCL::GetLumpedMass( i)*absdet;
                coupT2[i]= integralp;
            }

            // write values into matrix
            for(int i=0; i<4; ++i) {
                if (!Numb.WithUnknowns( i)) continue;
                M[Numb.num[i]]+= coup[i];
                const IdxT xidx_i= Xidx[Numb.num[i]];
                if (xidx_i!=NoIdx)
                    M[xidx_i]+= coupT2[i]*(1 - 2*sign[i]) + sign[i]*coup[i];
            }
        }
    }
}

void SetupLumpedMass_vecP2(const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const BndCondCL& bnd)
{
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP2CL Numb;

    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        Numb.assign( *sit, RowIdx, bnd);
        for(int i=0; i<10; ++i)
            if (Numb.WithUnknowns( i)) {
                const double contrib= P2DiscCL::GetLumpedMass( i)*absdet;
                M[Numb.num[i]  ]+= contrib;
                M[Numb.num[i]+1]+= contrib;
                M[Numb.num[i]+2]+= contrib;
            }
    }
}

void SetupLumpedMass_vecP2R (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const VecDescCL& lset, const BndDataCL<>& lsetbnd, const BndCondCL& bnd)
{
    const ExtIdxDescCL& Xidx= RowIdx.GetXidx();
    M.resize( RowIdx.NumUnknowns());

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP2CL Numb;
    double contribExt[4];

    InterfaceTetraCL cut;
    LocalP2CL<> p1abs_p[4][8], p1abs_n[4][8]; // extended basis functions on pos./neg. part, resp., for each of the 8 regular children
    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        loc_phi.assign( *sit, lset, lsetbnd);
        cut.Init( *sit, loc_phi);
        Numb.assign( *sit, RowIdx, bnd);
        // write standard FE values into matrix
        for(int i=0; i<10; ++i)
            if (Numb.WithUnknowns( i)) {
                const double contrib= P2DiscCL::GetLumpedMass( i)*absdet;
                M[Numb.num[i]  ]+= contrib;
                M[Numb.num[i]+1]+= contrib;
                M[Numb.num[i]+2]+= contrib;
            }

        if (cut.Intersects()) { // extended basis functions have only support on tetra intersecting Gamma!
            P2RidgeDiscCL::GetExtBasisOnChildren(p1abs_p, p1abs_n, loc_phi);
            for(int i=0; i<4; ++i)
                contribExt[i]= 0;
            // compute integrals    int_T v_i^R dx
            for (int ch= 0; ch < 8; ++ch) {
                cut.ComputeCutForChild( ch);
                for(int i=0; i<4; ++i) {
                    contribExt[i]+= cut.quad( p1abs_p[i][ch], absdet, true);   // integrate on positive part
                    contribExt[i]+= cut.quad( p1abs_n[i][ch], absdet, false);  // integrate on negative part
                }
            }

            // write extended values into matrix
            for(int i=0; i<4; ++i) {
                if (!Numb.WithUnknowns( i)) continue;
                const IdxT xidx_i= Xidx[Numb.num[i]];
                if (xidx_i!=NoIdx) {
                    M[xidx_i  ]+= contribExt[i];
                    M[xidx_i+1]+= contribExt[i];
                    M[xidx_i+2]+= contribExt[i];
                }
            }
        }
    }
}

void SetupLumpedMass (const MultiGridCL& MG, VectorCL& M, const IdxDescCL& RowIdx, const BndCondCL& bnd, const VecDescCL* lsetp, const BndDataCL<>* lsetbnd)
{
    switch(RowIdx.GetFE())
    {
    case P1_FE:
        SetupLumpedMass_P1( MG, M, RowIdx, bnd); break;
    case P1X_FE:
        SetupLumpedMass_P1X( MG, M, RowIdx, *lsetp, *lsetbnd, bnd); break;
    case vecP2_FE:
        SetupLumpedMass_vecP2( MG, M, RowIdx, bnd); break;
    case vecP2R_FE:
        SetupLumpedMass_vecP2R( MG, M, RowIdx, *lsetp, *lsetbnd, bnd); break;
    default:
        throw DROPSErrCL("SetupLumpedMass not implemented for this FE type");
    }
}

} // end of namespace DROPS
