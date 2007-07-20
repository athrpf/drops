/// \file
/// \brief classes that constitute the 2-phase Stokes problem

#include "num/discretize.h"

namespace DROPS
{

// -----------------------------------------------------------------------------
//                        Routines for SetupSystem2
// -----------------------------------------------------------------------------

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
    IdxT prNumb[4], numVerts= std::distance( MG.GetTriangVertexBegin(lvl), MG.GetTriangVertexEnd(lvl));
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
    Quad2CL<double> rho, Phi, kreuzterm;
    const double mu_p= _Coeff.mu( 1.0),
                 mu_n= _Coeff.mu( -1.0),
                 rho_p= _Coeff.rho( 1.0),
                 rho_n= _Coeff.rho( -1.0);

    SMatrixCL<3,3> T;

    double coupA[10][10], coupAk[10][10][3][3],
           coupM[10][10];
    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    LocalP2CL<> aij, akreuz[3][3];

    P2DiscCL::GetGradientsOnRef( GradRef);
    P2DiscCL::GetGradientsOnRef( GradRefLP1);

    InterfacePatchCL patch;

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
            // rhs = f + rho*g
            rhs+= Quad2CL<Point3DCL>( rho_const*_Coeff.g);

            P2DiscCL::GetGradients( Grad, GradRef, T);
            // compute all couplings between HatFunctions on edges and verts
            for (int i=0; i<10; ++i)
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
        else { // We are at the phase boundary.
            //  rhs = f + rho*g    with rho = rho( Phi)
            // TODO: We should integrate on the children, but we want to use the special quad-rule .quadP2 below.
            Phi.assign( *sit, ls, t);
            rho= Phi;
            rho.apply( _Coeff.rho);
            rhs+= Quad2CL<Point3DCL>( _Coeff.g)*rho;

            // compute all couplings between HatFunctions on edges and verts
            std::memset( coupA, 0, 10*10*sizeof( double));
            std::memset( coupAk, 0, 10*10*3*3*sizeof( double));
            P2DiscCL::GetGradients( GradLP1, GradRefLP1, T);
            for (int i= 0; i < 10; ++i) {
                GradLP2[i].assign( GradLP1[i]);
            }

            for (int i=0; i<10; ++i)
                for (int j=0; j<=i; ++j) {
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
                        const double cAp= patch.quad( aij, absdet, true);       // integrate on positive part
                        const double cAn= P1DiscCL::Quad( aij)*absdet/8. - cAp; // integrate on negative part
                        coupA[j][i]+= cAp*mu_p + cAn*mu_n;
                        coupA[i][j]= coupA[j][i];
                        for (int k= 0; k < 3; ++k)
                            for (int l= 0; l < 3; ++l) {
                                const double cAkp= patch.quad( akreuz[k][l], absdet, true);  // integrate on positive part
                                const double cAkn= P1DiscCL::Quad( akreuz[k][l])*absdet/8. - cAkp; // integrate on negative part
                                coupAk[i][j][k][l]+= cAkp*mu_p + cAkn*mu_n;
                                coupAk[j][i][l][k]= coupAk[i][j][k][l];
                            }
                    }
                }

            for (int i=0; i<10; ++i) // FIXME: This belongs in the preceding loop.
                for (int j=0; j<=i; ++j) {
                    // M
                    // FIXME: Ad hoc integration
                    const double cM= rho[4]*P2DiscCL::GetMass( j, i)*absdet;
                    coupM[i][j]= cM;
                    coupM[j][i]= cM;
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
                                mA( n.num[i]+k, n.num[j]+l)+= coupAk[j][i][l][k];

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
                                cplA->Data[n.num[i]+k]-= coupAk[j][i][l][k]*tmp[l];
                        }
                        cplM->Data[n.num[i]  ]-= cM*tmp[0];
                        cplM->Data[n.num[i]+1]-= cM*tmp[1];
                        cplM->Data[n.num[i]+2]-= cM*tmp[2];
                    }
                }
                tmp= rhs.quadP2( i, absdet);
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

        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const double cA= Quad2CL<>(dot( Grad[i], Grad[j]) * mu).quad( absdet);
                coupA[i][j]= cA;
                coupA[j][i]= cA;

// FIXME:                const double cM= rho.quadP2(i,j, absdet);
                const double cM= rho[4]*P2DiscCL::GetMass( j, i)*absdet;
                coupM[i][j]= cM;
                coupM[j][i]= cM;
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

} // end of namespace DROPS
