/// \file 
/// \brief classes that constitute the 2-phase Stokes problem

#include "num/discretize.h"

namespace DROPS
{

// -----------------------------------------------------------------------------
//                        Routines for SetupSystem2
// -----------------------------------------------------------------------------
  
template <class CoeffT>
void SetupSystem2_P2P0( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, const LevelsetP2CL&, double t)
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
        absdet= fabs( det);
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
void SetupSystem2_P2P1( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, const LevelsetP2CL&, double t)
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
        absdet= fabs( det);
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
void SetupSystem2_P2P1D( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, MatDescCL* B, VecDescCL* c, const LevelsetP2CL&, double t)
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
        absdet= fabs( det);
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
void SetupRhs2_P2P0( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL&, double t)
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
        absdet= fabs( det);

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
void SetupRhs2_P2P1( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL&, double t)
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
        absdet= fabs( det);
    
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
void SetupRhs2_P2P1D( const MultiGridCL& MG, const CoeffT&, const StokesBndDataCL& BndData, VecDescCL* c, const LevelsetP2CL&, double t)
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
        absdet= fabs( det);
    
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
        absdet= fabs( det);
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
        absdet= fabs( det);
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
    Quad2CL<double> rho, mu, Phi, kreuzterm;
        
    SMatrixCL<3,3> T;
    
    double coupA[10][10], coupM[10][10];
    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
    
        rhs.assign( *sit, _Coeff.f, t);
        Phi.assign( *sit, ls, t);
        // collect some information about the edges and verts of the tetra
        // and save it n.
        n.assign( *sit, *A->RowIdx, _BndData.Vel);

        // rho = rho( Phi),    mu= mu( Phi)
        rho= Phi;     rho.apply( _Coeff.rho);
        mu=  Phi;     mu.apply( _Coeff.mu);

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( _Coeff.g)*rho;
        
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const Quad2CL<double> dotGrad( dot( Grad[i], Grad[j]) * mu);
                coupA[i][j]= coupA[j][i]= dotGrad.quad( absdet);
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j, absdet);
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
                            {
                                // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm[m]= Grad[i][m][l] * Grad[j][m][k] * mu[m];

                                mA( n.num[i]+k, n.num[j]+l)+= kreuzterm.quad( absdet);
                            }
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
                            {
                                // kreuzterm = \int mu * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm[m]= Grad[i][m][l] * Grad[j][m][k] * mu[m];
                                cplA->Data[n.num[i]+k]-= kreuzterm.quad( absdet)*tmp[l];
                            }
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
        absdet= fabs( det);
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
                const Quad2CL<double> dotGrad( dot( Grad[i], Grad[j]) * mu);
                coupA[i][j]= coupA[j][i]= dotGrad.quad( absdet);
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j, absdet);
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
void InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2( MatDescCL* B, VecDescCL* c, const LevelsetP2CL& lset, double t) const
// Set up matrix B and rhs c
{
    std::cerr << "entering SetupSystem2: " << B->RowIdx->NumUnknowns << " prs. ";
    switch (prFE_)
    {
      case P0_FE:
        SetupSystem2_P2P0( _MG, _Coeff, _BndData, B, c, lset, t); break;
      case P1_FE: 
        SetupSystem2_P2P1( _MG, _Coeff, _BndData, B, c, lset, t); break;
      case P1D_FE: 
        SetupSystem2_P2P1D( _MG, _Coeff, _BndData, B, c, lset, t); break;
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
        SetupRhs2_P2P0( _MG, _Coeff, _BndData, c, lset, t); break;
      case P1_FE: 
        SetupRhs2_P2P1( _MG, _Coeff, _BndData, c, lset, t); break;
      case P1D_FE: 
        SetupRhs2_P2P1D( _MG, _Coeff, _BndData, c, lset, t); break;
      default:
        throw DROPSErrCL("InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2 not implemented for this FE type");
    }
}

} // end of namespace DROPS
