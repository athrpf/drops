/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#include "num/discretize.h"
#include <cstring>

namespace DROPS
{

// -----------------------------------------------------------------------------
//                        Routines for SetupSystem2
// -----------------------------------------------------------------------------

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupMatrix2( MatDescCL* B, MatDescCL* BT) const
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    MatrixBuilderCL mBT( &BT->Data, BT->RowIdx->NumUnknowns, BT->ColIdx->NumUnknowns);
    const Uint lvl= B->GetRowLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin( lvl),
         send= _MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *B->ColIdx, _BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *B->RowIdx);
// Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                tmp= Grad[vel].quadP1( pr, absdet);
                mB( prNumb[pr], n.num[vel])  -=  tmp[0];
                mB( prNumb[pr], n.num[vel]+1)-=  tmp[1];
                mB( prNumb[pr], n.num[vel]+2)-=  tmp[2];
                mBT( n.num[vel], prNumb[pr])  -=  tmp[0];
                mBT( n.num[vel]+1, prNumb[pr])-=  tmp[1];
                mBT( n.num[vel]+2, prNumb[pr])-=  tmp[2];                }
        }
    }
    mB.Build();
    mBT.Build();
}

template <class CoeffT>
void SetupSystem2_P2P0( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, double t)
// P2 / P0 FEs for vel/pr
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    c->Clear();
    const Uint lvl= B->GetRowLevel();
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;

    const Uint pidx= B->RowIdx->GetIdx();

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *B->ColIdx, BndData.Vel);
        const IdxT prNumbTetra= sit->Unknowns(pidx);

        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
        {
            tmp= Grad[vel].quad( absdet);
            mB( prNumbTetra, n.num[vel])  -=  tmp[0];
            mB( prNumbTetra, n.num[vel]+1)-=  tmp[1];
            mB( prNumbTetra, n.num[vel]+2)-=  tmp[2];
        }
        else { // put coupling on rhs
                typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                           : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                c->Data[ prNumbTetra]+= inner_prod( Grad[vel].quad( absdet), tmp);
            }
        }
    }
    mB.Build();
}

template <class CoeffT>
void SetupSystem2_P2P1( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, double t)
// P2 / P1 FEs (Taylor-Hood) for vel/pr
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    c->Clear();
    const Uint lvl= B->GetRowLevel();
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
        n.assign( *sit, *B->ColIdx, BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *B->RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  tmp[0];
                    mB( prNumb[pr], n.num[vel]+1)-=  tmp[1];
                    mB( prNumb[pr], n.num[vel]+2)-=  tmp[2];
                }
            else { // put coupling on rhs
                typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
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

template <class CoeffT>
void SetupSystem2_P2P1X( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, const ExtIdxDescCL& Xidx, double t)
// P2 / P1X FEs (X=extended) for vel/pr
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    c->Clear();
    const Uint lvl= B->GetRowLevel();
    IdxT prNumb[4];
    LocalNumbP2CL n;
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;
    InterfacePatchCL cut;

    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        n.assign( *sit, *B->ColIdx, BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *B->RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  tmp[0];
                    mB( prNumb[pr], n.num[vel]+1)-=  tmp[1];
                    mB( prNumb[pr], n.num[vel]+2)-=  tmp[2];
                }
            else { // put coupling on rhs
                typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[vel]).GetBndFun();
                tmp= vel<4 ? bf( sit->GetVertex( vel)->GetCoord(), t)
                           : bf( GetBaryCenter( *sit->GetEdge( vel-4)), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr, absdet), tmp);
            }
        }

        cut.Init( *sit, lset.Phi);
        if (!cut.Intersects()) continue; // extended basis functions have only support on tetra intersecting Gamma!
        for(int pr=0; pr<4; ++pr) {
            // compute the integrals
            // I = \int_{T_+} grad v_vel p_pr dx  -  C \int_{T}grad v_vel p_pr dx,
            // where C= (sign Phi_pr==1) \in {0,1} and T_+ = T \cap \Omega_2 (positive part)
            LocalP2CL<Point3DCL> gradv_p;
            const IdxT xidx= Xidx[prNumb[pr]];
            if (xidx==NoIdx) continue;

            for(int vel=0; vel<10; ++vel)
            {
                gradv_p[pr]= Grad[vel][pr];
                for (int vert=0; vert<4; ++vert)
                    if (vert!=pr)
                        gradv_p[EdgeByVert(pr,vert)+4]= 0.25*(Grad[vel][pr]+Grad[vel][vert]);

                Point3DCL integral;
                const bool is_pos= cut.GetSign(pr)==1;
                // for C=0 (<=> !is_pos) we have I = -\int_{T_-} grad v_vel p_pr dx
                for (int ch=0; ch<8; ++ch)
                {
                    cut.ComputeCutForChild(ch);
                    integral+= cut.quad( gradv_p, absdet, !is_pos); // integrate on other part
                }

                // for C=0 we have I = -\int_{T_-} grad v_vel p_pr dx
                if (is_pos) integral= -integral;

                if (n.WithUnknowns( vel))
                {
                    mB( xidx, n.num[vel])  -=  integral[0];
                    mB( xidx, n.num[vel]+1)-=  integral[1];
                    mB( xidx, n.num[vel]+2)-=  integral[2];
                }
                else { // put coupling on rhs
                    typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
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

template <class CoeffT>
void SetupSystem2_P2P1D( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, double t)
// P2 / P1D FEs for vel/pr
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    c->Clear();
    const Uint lvl= B->GetRowLevel();
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
        n.assign( *sit, *B->ColIdx, BndData.Vel);
        GetLocalNumbP1DNoBnd( prNumb, *sit, *B->RowIdx);
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel) {
            if (n.WithUnknowns( vel))
                for(int pr=0; pr<4; ++pr) {
                    tmp= Grad[vel].quadP1D( pr, absdet);
                    mB( prNumb[pr], n.num[vel])  -=  tmp[0];
                    mB( prNumb[pr], n.num[vel]+1)-=  tmp[1];
                    mB( prNumb[pr], n.num[vel]+2)-=  tmp[2];
                }
            else { // put coupling on rhs
                typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
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

template <class CoeffT>
void SetupRhs2_P2P0( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, double t)
{
    c->Clear();
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


template <class CoeffT>
void SetupRhs2_P2P1( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, double t)
{
    c->Clear();

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

template <class CoeffT>
void SetupRhs2_P2P1X( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL& lset, const ExtIdxDescCL& Xidx, double t)
// P2 / P1X FEs (X=extended) for vel/pr
{
    c->Clear();
    const Uint lvl=  c->GetLevel();
    const Uint pidx= c->RowIdx->GetIdx();
    IdxT prNumb[4];
    bool IsOnDirBnd[10];
    Quad2CL<Point3DCL> Grad_vel, GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    Point3DCL tmp;
    InterfacePatchCL cut;

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
        cut.Init( *sit, lset.Phi);
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

template <class CoeffT>
void SetupRhs2_P2P1D( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, double t)
{
    c->Clear();

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

template<class CoeffT>
void SetupPrMass_P0(const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* matM, const LevelsetP2CL& lset)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;
    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();

    SmoothedJumpCL nu_invers( 1./Coeff.mu(0), 1./Coeff.mu(1), Coeff.mu);
    Quad2CL<double> nu_inv;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint ls_lvl = ls.GetLevel();
    LocalP2CL<> locallset;
    const Uint pidx= matM->RowIdx->GetIdx();

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

template<class CoeffT>
void SetupPrMass_P1(const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* matM, const LevelsetP2CL& lset)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;
    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();
    IdxT prNumb[4];

    SmoothedJumpCL nu_invers( 1./Coeff.mu(0), 1./Coeff.mu(1), Coeff.mu);
    Quad2CL<double> nu_inv;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint ls_lvl = ls.GetLevel();
    LocalP2CL<> locallset;

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

        GetLocalNumbP1NoBnd( prNumb, *sit, *matM->RowIdx);
        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                M_pr( prNumb[i], prNumb[j])+= nu_inv.quadP1(i,j, absdet);
    }
    M_pr.Build();
}

template<class CoeffT>
void SetupPrMass_P1X(const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* matM, const LevelsetP2CL& lset, const ExtIdxDescCL& Xidx)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;
    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();
    IdxT prNumb[4];
    double coup[4][4], coupT2[4][4];

    const double nu_inv_p= 1./Coeff.mu( 1.0),
                 nu_inv_n= 1./Coeff.mu( -1.0);
    double integralp, integraln;
    InterfacePatchCL cut;
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

    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, sit) {
        const double absdet= sit->GetVolume()*6.;
        cut.Init( *sit, lset.Phi);
        const bool nocut= !cut.Intersects();
        GetLocalNumbP1NoBnd( prNumb, *sit, *matM->RowIdx);
        if (nocut) {
            const double nu_inv= cut.GetSign( 0) == 1 ? nu_inv_p : nu_inv_n;
            for(int i= 0; i < 4; ++i) {
                for(int j= 0; j < i; ++j) {
                    coup[j][i]= nu_inv*P1DiscCL::GetMass( i, j)*absdet;
                    coup[i][j]= coup[j][i];
                }
                coup[i][i]= nu_inv*P1DiscCL::GetMass( i, i)*absdet;
            }

            // write values into matrix
            for(int i=0; i<4; ++i)
                for(int j=0; j<4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= coup[i][j];
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
                        M_pr( xidx_i, xidx_j)+= coupT2[i][j]*(1-sign[i]-sign[j]) + sign[i]*sign[j]*coup[i][j];
                }
            }
        }
    }

    M_pr.Build();
}

template<class CoeffT>
void SetupPrMass_P1D(const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* matM, const LevelsetP2CL& lset)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;
    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();
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

        GetLocalNumbP1DNoBnd( prNumb, *sit, *matM->RowIdx);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                M_pr( prNumb[i], prNumb[j])+= nu_inv.quadP1D( i, j, absdet);
    }
    M_pr.Build();
}


// -----------------------------------------------------------------------------
//                        Routines for SetupPrStiff
// -----------------------------------------------------------------------------

template <class CoeffT>
void SetupPrStiff_P1( const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* A_pr, const LevelsetP2CL& lset)
{
    MatrixBuilderCL A( &A_pr->Data, A_pr->RowIdx->NumUnknowns, A_pr->ColIdx->NumUnknowns);
    const Uint lvl= A_pr->GetRowLevel();
    const Uint idx= A_pr->RowIdx->GetIdx();
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

        P1DiscCL::GetGradients( G,det,*sit);
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
            UnknownIdx[i]= sit->GetVertex( i)->Unknowns( idx);
        }
        for(int i=0; i<4; ++i)    // assemble row i
            for(int j=0; j<4;++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i];
    }
    A.Build();
}

// TODO: As in SetupPrMass_P1X, replace the smoothed density-function with integration
//       over the inner and outer part.
template <class CoeffT>
void SetupPrStiff_P1X( const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* A_pr, const LevelsetP2CL& lset, const ExtIdxDescCL& Xidx)
{
    MatrixBuilderCL A( &A_pr->Data, A_pr->RowIdx->NumUnknowns, A_pr->ColIdx->NumUnknowns);
    const Uint lvl= A_pr->GetRowLevel();
    const Uint idx= A_pr->RowIdx->GetIdx();
    SMatrixCL<3,4> G;
    double coup[4][4], coupT2[4][4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];
    bool sign[4];
    InterfacePatchCL cut;

    SmoothedJumpCL rho_invers( 1./Coeff.rho(0), 1./Coeff.rho(1), Coeff.rho);
    Quad2CL<double> rho_inv;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const Uint ls_lvl = ls.GetLevel();
    LocalP2CL<> locallset, rho_inv_2( 1./Coeff.rho(1));

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
        cut.Init( *sit, lset.Phi);
        const bool nocut= !cut.Intersects();

        P1DiscCL::GetGradients( G,det,*sit);
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
            UnknownIdx[i]= sit->GetVertex( i)->Unknowns( idx);
        }
        for(int i=0; i<4; ++i)    // assemble row i
            for(int j=0; j<4;++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i];

        if (nocut) continue; // extended basis functions have only support on tetra intersecting Gamma!

        for(int i=0; i<4; ++i) {
            sign[i]= cut.GetSign(i)==1;
            for(int j=0; j<=i; ++j) {
                // compute the integrals
                // \int_{T_2} grad_i grad_j dx,    where T_2 = T \cap \Omega_2
                double integral= 0;

                for (int ch=0; ch<8; ++ch)
                {
                    cut.ComputeCutForChild(ch);
                    integral+= cut.quad( rho_inv_2, absdet, true); // integrate on positive part
                }
                coupT2[j][i]= integral*( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) );
                coupT2[i][j]= coupT2[j][i];
            }
        }

        // write values into matrix
        for(int i=0; i<4; ++i)
        {
            const IdxT xidx_i= Xidx[UnknownIdx[i]];
            for(int j=0; j<4; ++j)
            { // tetra intersects Gamma => Xidx defined for all DoFs
                const IdxT xidx_j= Xidx[UnknownIdx[j]];
                if (xidx_j!=NoIdx)
                    A( UnknownIdx[i], xidx_j)+= coupT2[i][j] - sign[j]*coup[i][j];
                if (xidx_i!=NoIdx)
                    A( xidx_i, UnknownIdx[j])+= coupT2[i][j] - sign[i]*coup[i][j];
                if (xidx_i!=NoIdx && xidx_j!=NoIdx)
                    A( xidx_i, xidx_j)+=        coupT2[i][j]*(1-sign[i]-sign[j]) + sign[i]*sign[j]*coup[i][j];
            }
        }
    }
    A.Build();
}

template <class CoeffT>
void SetupPrStiff_P1D( const MultiGridCL& MG, const CoeffT& Coeff, MatDescCL* A_pr, const LevelsetP2CL& lset)
{
    MatrixBuilderCL A( &A_pr->Data, A_pr->RowIdx->NumUnknowns, A_pr->ColIdx->NumUnknowns);
    const Uint lvl= A_pr->GetRowLevel();
    const Uint idx= A_pr->RowIdx->GetIdx();
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


template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SmoothVel( VelVecDescCL* v, int num, double tau)
{
    const VectorCL diag= A.Data.GetDiag();

    for (int i=0; i<num; ++i)
        v->Data-= tau*((A.Data*v->Data)/diag);
}

template <class Coeff>
/// Needed for preconditioning of the Schur complement. Uses natural
/// boundary conditions for the pressure unknowns.
void InstatStokes2PhaseP2P1CL<Coeff>::SetupPrMass(MatDescCL* matM, const LevelsetP2CL& lset) const
{
    switch (prFE_)
    {
      case P0_FE:
        SetupPrMass_P0( _MG, _Coeff, matM, lset); break;
      case P1_FE:
        SetupPrMass_P1( _MG, _Coeff, matM, lset); break;
      case P1X_FE:
        SetupPrMass_P1X( _MG, _Coeff, matM, lset, Xidx_); break;
      case P1D_FE:
        SetupPrMass_P1D( _MG, _Coeff, matM, lset); break;
      default:
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupPrMass not implemented for this FE type");
    }
}


template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupPrStiff( MatDescCL* A_pr, const LevelsetP2CL& lset) const
/// Needed for preconditioning of the Schur complement. Uses natural
/// boundary conditions for the pressure unknowns.
{
    switch (prFE_)
    {
      case P1_FE:
        SetupPrStiff_P1( _MG, _Coeff, A_pr, lset); break;
      case P1X_FE:
        SetupPrStiff_P1X( _MG, _Coeff, A_pr, lset, Xidx_); break;
      case P1D_FE:
        SetupPrStiff_P1D( _MG, _Coeff, A_pr, lset); break;
      default:
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupPrStiff not implemented for this FE type");
    }
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::InitVel(VelVecDescCL* vec, instat_vector_fun_ptr LsgVel, double t0) const
{
    VectorCL& lsgvel= vec->Data;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

    for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>( _MG).GetTriangVertexBegin( lvl),
         send= const_cast<const MultiGridCL&>( _MG).GetTriangVertexEnd( lvl);
         sit != send; ++sit) {
        if (!_BndData.Vel.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel(sit->GetCoord(), t0));
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL sit= const_cast<const MultiGridCL&>( _MG).GetTriangEdgeBegin( lvl),
         send= const_cast<const MultiGridCL&>( _MG).GetTriangEdgeEnd( lvl);
         sit != send; ++sit) {
        if (!_BndData.Vel.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t0));
    }
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem1( MatDescCL* A, MatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, const LevelsetP2CL& lset, double t) const
// Set up matrices A, M and rhs b (depending on phase bnd)
{
    const IdxT num_unks_vel= A->RowIdx->NumUnknowns;

    MatrixBuilderCL mA( &A->Data, num_unks_vel, num_unks_vel),
                    mM( &M->Data, num_unks_vel, num_unks_vel);
    b->Clear();
    cplM->Clear();
    cplA->Clear();

    const Uint lvl = A->GetRowLevel();

    LocalNumbP2CL n;

    std::cerr << "entering SetupSystem1: " << num_unks_vel << " vels. ";

    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    LocalP1CL<Point3DCL> GradRefLP1[10], GradLP1[10];
    LocalP2CL<Point3DCL> GradLP2[10];
    Quad2CL<double> Ones( 1.), kreuzterm;
    const double mu_p= _Coeff.mu( 1.0),
                 mu_n= _Coeff.mu( -1.0),
                 rho_p= _Coeff.rho( 1.0),
                 rho_n= _Coeff.rho( -1.0);

    SMatrixCL<3,3> T;

    double coupA[10][10], coupAk[10][10][3][3],
           coupM[10][10], rho_phi[10];
    double det, absdet, cAp, cAn, cAkp, cAkn, intHat_p, intHat_n;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> aij, akreuz[3][3], phi_i, ones( 1.);

    P2DiscCL::GetGradientsOnRef( GradRef);
    P2DiscCL::GetGradientsOnRef( GradRefLP1);

    InterfacePatchCL patch;
    BaryCoordCL* nodes;
    LocalP2CL<>p2[10];
    double intpos, intneg;
    for (int k=0; k<10; ++k)
        p2[k][k]=1.;
    Quad5CL<> q[10][48]; //there exists maximally 8*6=48 SubTetras

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        rhs.assign( *sit, _Coeff.f, t);

        // collect some information about the edges and verts of the tetra
        // and save it n.
        n.assign( *sit, *A->RowIdx, _BndData.Vel);
        patch.Init( *sit, lset.Phi);
        const bool nocut= !patch.Intersects();
        if (nocut) {
            const double mu_const= patch.GetSign( 0) == 1 ? mu_p : mu_n;
            const double rho_const= patch.GetSign( 0) == 1 ? rho_p : rho_n;

            P2DiscCL::GetGradients( Grad, GradRef, T);
            // compute all couplings between HatFunctions on edges and verts
            for (int i=0; i<10; ++i)
            {
                rho_phi[i]= rho_const*Ones.quadP2( i, absdet);
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
            std::memset( coupA, 0, 10*10*sizeof( double));
            std::memset( coupAk, 0, 10*10*3*3*sizeof( double));
            std::memset( rho_phi, 0, 10*sizeof( double));
            P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);
            for (int i= 0; i < 10; ++i) {
                GradLP2[i].assign( GradLP1[i]);
            }

            double Vol_p= 0;
            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeCutForChild( ch);
                Vol_p+= patch.quad( ones, absdet, true);
                for (int i=0; i<10; ++i) {
                    // init phi_i =  i-th P2 hat function
                    phi_i[i]= 1.; phi_i[i==0 ? 9 : i-1]= 0.;
                    patch.quadBothParts( intHat_p, intHat_n, phi_i, absdet);
                    rho_phi[i]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_i
                }
            }
            patch.ComputeSubTets();
            for (Uint k=0; k<patch.GetNumTetra(); ++k)
            {
                nodes = Quad5CL<>::TransformNodes(patch.GetTetra(k));
                for (Uint j=0; j<10; ++j)
                    q[j][k].assign(p2[j], nodes);
                delete[] nodes;
            }
            for (int i=0; i<10; ++i)
                for (int j=0; j<=i; ++j) {
                    // M
                    intpos = 0.;
                    intneg = 0.;
                    for (Uint k=0; k<patch.GetNumTetra(); k++)
                        if (k<patch.GetNumNegTetra())
                            intneg += Quad5CL<>(q[i][k]*q[j][k]).quad(absdet*VolFrac(patch.GetTetra(k)));
                        else
                            intpos += Quad5CL<>(q[i][k]*q[j][k]).quad(absdet*VolFrac(patch.GetTetra(k)));
                    coupM[i][j]= rho_p*intpos + rho_n*intneg;
                    coupM[j][i]= rho_p*intpos + rho_n*intneg;

                    // A
                    aij= dot( GradLP2[i], GradLP2[j]);
                    for (int k= 0; k < 3; ++k)
                        for (int l= 0; l < 3; ++l) {
                            // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                            for (size_t m= 0; m < /* #Components of akreuz[k][l] */ 10;  ++m)
                                akreuz[k][l][m]= GradLP2[i][m][l] * GradLP2[j][m][k];
                        }

                    for (int ch= 0; ch < 8; ++ch) {
                        patch.ComputeCutForChild( ch);
                        patch.quadBothParts( cAp, cAn, aij, absdet);  // integrate on positive and negative part
                        const double cA= cAp*mu_p + cAn*mu_n;
                        coupA[j][i]+= cA;
                        for (int k= 0; k < 3; ++k)
                            for (int l= 0; l < 3; ++l) {
                                patch.quadBothParts( cAkp, cAkn, akreuz[k][l], absdet);  // integrate on positive and negative part
                                const double cAk= cAkp*mu_p + cAkn*mu_n;
                                coupAk[i][j][k][l]+= cAk;
                            }
                    }
                    if (i!=j)
                    { // copy computed entries, as local stiffness matrices coupA, coupAk are symmetric
                        coupA[i][j]= coupA[j][i];
                        for (int k= 0; k < 3; ++k)
                            for (int l= 0; l < 3; ++l)
                                coupAk[j][i][l][k]= coupAk[i][j][k][l];
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
                        for (int k=0; k<3; ++k)
                            for (int l=0; l<3; ++l)
                                mA( n.num[i]+k, n.num[j]+l)+= coupAk[i][j][k][l];
                        mM( n.num[i],   n.num[j]  )+= coupM[j][i];
                        mM( n.num[i]+1, n.num[j]+1)+= coupM[j][i];
                        mM( n.num[i]+2, n.num[j]+2)+= coupM[j][i];
                    }
                    else // put coupling on rhs
                    {
                        typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                        bnd_val_fun bf= _BndData.Vel.GetBndSeg( n.bndnum[j]).GetBndFun();
                        tmp= j<4 ? bf( sit->GetVertex( j)->GetCoord(), t)
                                 : bf( GetBaryCenter( *sit->GetEdge( j-4)), t);
                        const double cA= coupA[j][i],
                                     cM= coupM[j][i];
                        for (int k=0; k<3; ++k)
                        {
                            cplA->Data[n.num[i]+k]-= cA*tmp[k];
                            for (int l=0; l<3; ++l)
                                cplA->Data[n.num[i]+k]-= coupAk[i][j][k][l]*tmp[l];
                        }
                        cplM->Data[n.num[i]  ]-= cM*tmp[0];
                        cplM->Data[n.num[i]+1]-= cM*tmp[1];
                        cplM->Data[n.num[i]+2]-= cM*tmp[2];
                    }
                }
                tmp= rhs.quadP2( i, absdet) + rho_phi[i]*_Coeff.g;
                b->Data[n.num[i]  ]+= tmp[0];
                b->Data[n.num[i]+1]+= tmp[1];
                b->Data[n.num[i]+2]+= tmp[2];
            }
    }

    mA.Build();
    mM.Build();
    std::cerr << A->Data.num_nonzeros() << " nonzeros in A, "
              << M->Data.num_nonzeros() << " nonzeros in M! " << std::endl;
}


template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs1( VecDescCL* b, const LevelsetP2CL& lset, double t) const
// Set up rhs b (depending on phase bnd)
{
    const Uint lvl = b->GetLevel();

    b->Clear();

    LocalNumbP2CL n;
    SMatrixCL<3,3> T;

    Quad2CL<Point3DCL> rhs;
    Quad2CL<double> Ones( 1.);
    LocalP2CL<> phi_i;

    const double rho_p= _Coeff.rho( 1.0),
                 rho_n= _Coeff.rho( -1.0);
    double rho_phi[10];
    double det, absdet, intHat_p, intHat_n;
    Point3DCL tmp;

    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    InterfacePatchCL patch;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);

        rhs.assign( *sit, _Coeff.f, t);

        // collect some information about the edges and verts of the tetra
        // and save it n.
        n.assign( *sit, *b->RowIdx, _BndData.Vel);
        patch.Init( *sit, lset.Phi);
        const bool nocut= !patch.Intersects();
        if (nocut) {
            const double rho_const= patch.GetSign( 0) == 1 ? rho_p : rho_n;

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
                patch.ComputeCutForChild( ch);
                for (int i=0; i<10; ++i) {
                    // init phi_i =  i-th P2 hat function
                    phi_i[i]= 1.; phi_i[i==0 ? 9 : i-1]= 0.;
                    patch.quadBothParts( intHat_p, intHat_n, phi_i, absdet);
                    rho_phi[i]+= rho_p*intHat_p + rho_n*intHat_n; // \int rho*phi_i
                }
            }
        }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (n.WithUnknowns( i))  // vert/edge i is not on a Dirichlet boundary
            {
                tmp= rhs.quadP2( i, absdet) + rho_phi[i]*_Coeff.g;
                b->Data[n.num[i]  ]+= tmp[0];
                b->Data[n.num[i]+1]+= tmp[1];
                b->Data[n.num[i]+2]+= tmp[2];
            }
    }
}


template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupMatrices1( MatDescCL* A,
    MatDescCL* M, const LevelsetP2CL& lset, double t) const
// Set up matrices A, M (depending on phase bnd)
{
    const IdxT num_unks_vel= A->RowIdx->NumUnknowns;
    MatrixBuilderCL mA( &A->Data, num_unks_vel, num_unks_vel),
                    mM( &M->Data, num_unks_vel, num_unks_vel);
    const Uint lvl= A->GetRowLevel();
    LocalNumbP2CL locn;

    std::cerr << "entering SetupMatrices1: " << num_unks_vel << " vels. ";

    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    Quad2CL<double> rho, mu, Phi, kreuzterm;
    LocalP2CL<> ls_loc;

    SMatrixCL<3,3> T;

    double coupA[10][10], coupM[10][10];
    double det, absdet;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    P2DiscCL::GetGradientsOnRef( GradRef);

    InterfacePatchCL patch;
    BaryCoordCL* nodes;
    LocalP2CL<>p2[10];
    double intpos, intneg;
    for (int k=0; k<10; ++k)
        p2[k][k]=1.;
    Quad5CL<> q[10][48]; //there exist maximally 48 SubTetras
    const double rho_p= _Coeff.rho( 1.0),
                 rho_n= _Coeff.rho( -1.0);

    for (MultiGridCL::const_TriangTetraIteratorCL
            sit= const_cast<const MultiGridCL&>( _MG).GetTriangTetraBegin( lvl),
            send= const_cast<const MultiGridCL&>( _MG).GetTriangTetraEnd( lvl);
            sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        // Collect information about Numbering of unknowns and boundary conditions.
        locn.assign( *sit, *A->RowIdx, _BndData.Vel);
        ls_loc.assign( *sit, ls, t); // needed for restrictions
        Phi.assign( ls_loc);
        // rho = rho( Phi),    mu= mu( Phi)
        rho=   Phi;
        rho.apply( _Coeff.rho);
        mu= Phi;
        mu.apply( _Coeff.mu);

        // rhs = f + rho*g
        rhs.assign( *sit, _Coeff.f, t);
        rhs+= Quad2CL<Point3DCL>( _Coeff.g)*rho;

        patch.Init(*sit, ls_loc);
        patch.ComputeSubTets();
        for (Uint k=0; k<patch.GetNumTetra(); ++k)
        {
            nodes = Quad5CL<>::TransformNodes(patch.GetTetra(k));
            for (Uint j=0; j<10; ++j)
                q[j][k].assign(p2[j], nodes);
            delete[] nodes;
        }

        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const double cA= Quad2CL<>(dot( Grad[i], Grad[j]) * mu).quad( absdet);
                coupA[i][j]= cA;
                coupA[j][i]= cA;

                intpos = 0.;
                intneg = 0.;
                for (Uint k=0; k<patch.GetNumTetra(); k++)
                    if (patch.GetNumNegTetra()>k)
                        intneg += Quad5CL<>(q[i][k]*q[j][k]).quad(absdet*VolFrac(patch.GetTetra(k)));
                    else
                        intpos += Quad5CL<>(q[i][k]*q[j][k]).quad(absdet*VolFrac(patch.GetTetra(k)));
                coupM[i][j]= rho_p*intpos + rho_n*intneg;
                coupM[j][i]= rho_p*intpos + rho_n*intneg;
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (locn.WithUnknowns( i)) { // vert/edge i is not on a Dirichlet boundary
                for(int j=0; j<10; ++j) {
                    if (locn.WithUnknowns( j)) { // vert/edge j is not on a Dirichlet boundary
                        mA( locn.num[i],   locn.num[j]  )+= coupA[j][i];
                        mA( locn.num[i]+1, locn.num[j]+1)+= coupA[j][i];
                        mA( locn.num[i]+2, locn.num[j]+2)+= coupA[j][i];
                        for (int k=0; k<3; ++k)
                            for (int l=0; l<3; ++l) {
                                // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm[m]= Grad[i][m][l] * Grad[j][m][k] * mu[m];
                                mA( locn.num[i]+k, locn.num[j]+l)+= kreuzterm.quad( absdet);
                            }
                        mM( locn.num[i],   locn.num[j]  )+= coupM[j][i];
                        mM( locn.num[i]+1, locn.num[j]+1)+= coupM[j][i];
                        mM( locn.num[i]+2, locn.num[j]+2)+= coupM[j][i];
                    }
                }
            }
    }

    mA.Build();
    mM.Build();
    std::cerr << A->Data.num_nonzeros() << " nonzeros in A, "
              << M->Data.num_nonzeros() << " nonzeros in M! " << std::endl;
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupMatricesMG (MGDataCL* matMG, const LevelsetP2CL& lset, double dt, double theta) const
// Set up the MG-hierarchy
{
    for(MGDataCL::iterator it= matMG->begin(); it!=matMG->end(); ++it) {
        MGLevelDataCL& tmp= *it;
        MatDescCL A;
        A.SetIdx( &tmp.Idx, &tmp.Idx);
        tmp.A.SetIdx( &tmp.Idx, &tmp.Idx);
        tmp.Mvel.SetIdx( &tmp.Idx, &tmp.Idx);
        std::cerr << "Create StiffMatrix for "
                  << (&tmp.Idx)->NumUnknowns << " unknowns." << std::endl;
        if(&tmp != &matMG->back()) {
            SetupMatrices1( &A, &tmp.Mvel, lset, t);
            tmp.A.Data.LinComb( 1./dt, tmp.Mvel.Data, theta, A.Data);
            if (it != matMG->begin())
                tmp.Mvel.Data.clear();
        }
        tmp.ABlock = &tmp.A.Data;
        if (matMG->StokesMG()) {
            tmp.B.SetIdx( &tmp.IdxPr , &tmp.Idx );
            tmp.BT.SetIdx( &tmp.Idx, &tmp.IdxPr );
            tmp.Mpr.SetIdx( &tmp.IdxPr, &tmp.IdxPr);
            std::cerr << "Create StokesMatrices2 for " << (&tmp.IdxPr)->NumUnknowns << " unknown" << std::endl;
            if (&tmp != &matMG->back()) {
                SetupMatrix2( &tmp.B, &tmp.BT);
                SetupPrMass(&tmp.Mpr, lset);
            }
            else
                transpose (tmp.B.Data, tmp.BT.Data);
        }
    }
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupLB (MatDescCL* A, VecDescCL* cplA, const LevelsetP2CL& lset, double t) const
// Set up the Laplace-Beltrami-matrix
{
    const IdxT num_unks_vel= A->RowIdx->NumUnknowns;

    MatrixBuilderCL mA( &A->Data, num_unks_vel, num_unks_vel);
    cplA->Clear();

    const Uint lvl = A->GetRowLevel();

    LocalNumbP2CL n;

    std::cerr << "entering SetupLB: " << num_unks_vel << " vels. ";

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

    InterfacePatchCL patch;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        patch.Init( *sit, lset.Phi);
        if (patch.Intersects()) { // We are at the phase boundary.
            n.assign( *sit, *A->RowIdx, _BndData.Vel);
            GetTrafoTr( T, det, *sit);
            absdet= std::fabs( det);
            P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);
            std::memset( coupA, 0, 10*10*sizeof( double));
            for (int i= 0; i < 10; ++i)
                GradLP2[i].assign( GradLP1[i]);

            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++ tri) {
                    surfTension= _Coeff.SurfTens;
                    for (int i= 0; i < 10; ++i) {
                        surfGrad[i].assign( GradLP1[i], &patch.GetBary( tri));
                        surfGrad[i].apply( patch, &InterfacePatchCL::ApplyProj);
                    }
                    for (int i=0; i<10; ++i)
                        for (int j=0; j<=i; ++j) {
                            // Laplace-Beltrami... Stabilization
                            LB= surfTension*dot( surfGrad[i], surfGrad[j]);
                            const double cLB= LB.quad( patch.GetFuncDet( tri));
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
                        else // put coupling on rhs
                        {
                            typedef typename StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                            bnd_val_fun bf= _BndData.Vel.GetBndSeg( n.bndnum[j]).GetBndFun();
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
    std::cerr << A->Data.num_nonzeros() << " nonzeros in A_LB" << std::endl;
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2( MatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const
// Set up matrix B and rhs c
{
    std::cerr << "entering SetupSystem2: " << B->RowIdx->NumUnknowns << " prs. ";
    switch (prFE_)
    {
      case P0_FE:
        SetupSystem2_P2P0( _MG, _Coeff, _BndData, B, c, t); break;
      case P1_FE:
        SetupSystem2_P2P1( _MG, _Coeff, _BndData, B, c, t); break;
      case P1X_FE:
        SetupSystem2_P2P1X( _MG, _Coeff, _BndData, B, c, lset, Xidx_, t); break;
      case P1D_FE:
        SetupSystem2_P2P1D( _MG, _Coeff, _BndData, B, c, t); break;
      default:
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2 not implemented for this FE type");
    }
    std::cerr << B->Data.num_nonzeros() << " nonzeros in B!" << std::endl;
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2( VecDescCL* c, const LevelsetP2CL& lset, double t) const
// Set up rhs c
{
    switch (prFE_)
    {
      case P0_FE:
        SetupRhs2_P2P0( _MG, _Coeff, _BndData, c, t); break;
      case P1_FE:
        SetupRhs2_P2P1( _MG, _Coeff, _BndData, c, t); break;
      case P1X_FE:
        SetupRhs2_P2P1X( _MG, _Coeff, _BndData, c, lset, Xidx_, t); break;
      case P1D_FE:
        SetupRhs2_P2P1D( _MG, _Coeff, _BndData, c, t); break;
      default:
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2 not implemented for this FE type");
    }
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupBdotv (VecDescCL* Bdotv, const VelVecDescCL* vel,
    const ExtIdxDescCL& v_idx, const LevelsetP2CL& lset, double t) const
{
    Bdotv->Clear();
    const Uint lvl= Bdotv->GetLevel();
    IdxT prNumb[4];
    LocalNumbP2CL num;

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
    InterfacePatchCL cut;

    DROPS_FOR_TRIANG_TETRA( _MG, lvl, sit) {
        cut.Init( *sit, lset.Phi);
        if (!cut.Intersects()) continue;

        num.assign( *sit, *v_idx.Idx, _BndData.Vel);
        GetLocalNumbP1NoBnd( prNumb, *sit, *Bdotv->RowIdx);
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        loc_u.assign( *sit, *vel, GetBndData().Vel, t);
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
                for (int i= 0; i < Quad5_2DCL<>::NumNodesC; ++i)
                    if (n[i].norm()>1e-8) n[i]/= n[i].norm();
                q1= dot( n, qu)*qdivu;
                for(int pr= 0; pr < 4; ++pr) {
                    const IdxT xidx( v_idx[prNumb[pr]]);
                    if (xidx == NoIdx) continue;
                    q2.assign( lp1[pr], p);
                    q2*= q1;
                    // n is the outer normal of {lset <= 0}; we need the outer normal of supp(pr-hat-function)\cap \Gamma).
                    Bdotv->Data[xidx]-= (cut.GetSign( pr) > 0 ? -1. : 1.)*q2.quad( cut.GetFuncDet( t));
                }
            }
        }
    }
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetIdx()
{
    IdxDescCL* vidx= &vel_idx;
    IdxDescCL* pidx= &pr_idx;
    b.SetIdx  ( vidx);
    c.SetIdx  ( pidx);
    A.SetIdx  ( vidx, vidx);
    B.SetIdx  ( pidx, vidx);
    prM.SetIdx( pidx, pidx);
    prA.SetIdx( pidx, pidx);
    M.SetIdx  ( vidx, vidx);
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::GetPrOnPart( VecDescCL& p_part, const LevelsetP2CL& lset, bool posPart)
{
    const Uint lvl= p.RowIdx->TriangLevel,
        idxnum= p.RowIdx->GetIdx();
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const MultiGridCL& mg= this->GetMG();

    p_part.SetIdx( p.RowIdx);
    VectorCL& pp= p_part.Data;
    pp= p.Data;

    // add extended pressure
    for( MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
        end= mg.GetTriangVertexEnd(lvl); it != end; ++it)
    {
        const IdxT nr= it->Unknowns(idxnum);
        if (Xidx_[nr]==NoIdx) continue;

        const bool is_pos= InterfacePatchCL::Sign( ls.val( *it))==1;
        if (posPart==is_pos) continue; // extended hat function ==0 on this part
        if (posPart)
            pp[nr]+= p.Data[Xidx_[nr]];
        else
            pp[nr]-= p.Data[Xidx_[nr]];
    }
}

//*****************************************************************************
//                               VelocityRepairCL
//*****************************************************************************
template<class StokesT>
  inline void
  VelocityRepairCL<StokesT>::post_refine ()
{
    VelVecDescCL loc_v;
    IdxDescCL    loc_vidx( vecP2_FE);
    VelVecDescCL& v= stokes_.v;
    Uint LastLevel= stokes_.GetMG().GetLastLevel();
    match_fun match= stokes_.GetMG().GetBnd().GetMatchFun();

    stokes_.CreateNumberingVel( LastLevel, &loc_vidx, match);
    if (LastLevel != v.RowIdx->TriangLevel) {
        std::cout << "LastLevel: " << LastLevel
                  << " old v->TriangLevel: " << v.RowIdx->TriangLevel << std::endl;
        throw DROPSErrCL( "VelocityRepairCL::post_refine: Sorry, not yet implemented.");
    }
    loc_v.SetIdx( &loc_vidx);
    RepairAfterRefineP2( stokes_.GetVelSolution( v), loc_v);
    v.Clear();
    stokes_.DeleteNumbering( v.RowIdx);

    stokes_.vel_idx.swap( loc_vidx);
    v.SetIdx( &stokes_.vel_idx);
    v.Data= loc_v.Data;
}

//*****************************************************************************
//                               PressureRepairCL
//*****************************************************************************
template<class StokesT>
  inline void
  PressureRepairCL<StokesT>::post_refine ()
{
    VecDescCL loc_p;
    IdxDescCL loc_pidx( P1_FE);
    VecDescCL& p= stokes_.p;
    match_fun match= stokes_.GetMG().GetBnd().GetMatchFun();

    stokes_.CreateNumberingPr( stokes_.GetMG().GetLastLevel(), &loc_pidx, match);
    loc_p.SetIdx( &loc_pidx);
    RepairAfterRefineP1( stokes_.GetPrSolution( p), loc_p);
    p.Clear();
    stokes_.DeleteNumbering( p.RowIdx);
    stokes_.pr_idx.swap( loc_pidx);
    p.SetIdx( &stokes_.pr_idx);
    p.Data= loc_p.Data;
}

template<class StokesT>
  inline void
  PressureRepairCL<StokesT>::pre_refine_sequence ()
{
    p1xrepair_= std::auto_ptr<P1XRepairCL>( new P1XRepairCL( stokes_.UsesXFEM(), stokes_.GetMG(), stokes_.p, stokes_.GetXidx()));
}

template<class StokesT>
  inline void
  PressureRepairCL<StokesT>::post_refine_sequence ()
{
    (*p1xrepair_)( ls_);
    p1xrepair_.reset();
}

} // end of namespace DROPS
