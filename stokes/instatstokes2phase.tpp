//**************************************************************************
// File:    instatstokes2phase.tpp                                         *
// Content: classes that constitute the 2-phase stokes-problem             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "num/discretize.h"

namespace DROPS
{

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::DeleteNumberingVel(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetAllVertexBegin(level), _MG.GetAllVertexEnd(level) );
    DeleteNumbOnSimplex( idxnum, _MG.GetAllEdgeBegin(level), _MG.GetAllEdgeEnd(level) );
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::DeleteNumberingPr(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetAllVertexBegin(level), _MG.GetAllVertexEnd(level) );
}


template<class _Cont, class DataT>
  inline DataT
      P2(const _Cont& dof, const DataT& , double v1, double v2, double v3)
{
    return dof[0] * FE_P2CL::H0( v1, v2, v3)
        + dof[1] * FE_P2CL::H1( v1, v2, v3) + dof[2] * FE_P2CL::H2( v1, v2, v3)
        + dof[3] * FE_P2CL::H3( v1, v2, v3) + dof[4] * FE_P2CL::H4( v1, v2, v3)
        + dof[5] * FE_P2CL::H5( v1, v2, v3) + dof[6] * FE_P2CL::H6( v1, v2, v3)
        + dof[7] * FE_P2CL::H7( v1, v2, v3) + dof[8] * FE_P2CL::H8( v1, v2, v3)
        + dof[9] * FE_P2CL::H9( v1, v2, v3);
}


template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupPrMass(MatDescCL* matM, const LevelsetP2CL& lset, double nu1, double nu2) const
// Sets up the mass matrix for the pressure
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;

    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl    = matM->GetRowLevel();
    const Uint pidx   = matM->RowIdx->GetIdx();

    IdxT prNumb[4];
    SmoothedJumpCL nu_invers( 1./nu1, 1./nu2, _Coeff.rho);
    Quad2CL<double> nu_inv;
    LevelsetP2CL::DiscSolCL ls= lset.GetSolution();
    double lsarray[10];
    
    // compute all couplings between HatFunctions on verts:
    // I( i, j) = int ( psi_i*psi_j, T_ref) * absdet
//    const double coupl_ii= 1./60.,  // = 1/120 + 2/15 * 1/16,
//                 coupl_ij= 1./120.; // =         2/15 * 1/16;

    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6.;
        if (ls.GetLevel()!=lvl)
            RestrictP2( *sit, ls, lsarray);
        else
            ls.GetDoF( *sit, lsarray);
        
        for(int i=0; i<4; ++i)
        {
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
            nu_inv.val[i]= lsarray[i];
        }
        nu_inv.val[4]= P2( lsarray, 1.0/*type-dummy*/, 0.25, 0.25, 0.25); 
        nu_inv.apply( nu_invers);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= //(i!=j ? coupl_ij : coupl_ii) * absdet;
                                                  nu_inv.quadP1(i,j, absdet);
    }
    M_pr.Build();
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupPrStiff( MatDescCL* A_pr) const
// Assumes, that indices for A_pr are set up. We know, there are only natural
// boundary conditions.
{
    MatrixBuilderCL A( &A_pr->Data, A_pr->RowIdx->NumUnknowns, A_pr->ColIdx->NumUnknowns);
    const Uint lvl= A_pr->GetRowLevel();
    const Uint idx= A_pr->RowIdx->GetIdx();
    SMatrixCL<3,4> G;
    double coup[4][4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];

    for (MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const DROPS::MultiGridCL&>( _MG).GetTriangTetraBegin( lvl),
         send= const_cast<const DROPS::MultiGridCL&>( _MG).GetTriangTetraEnd( lvl);
         sit != send; ++sit) 
    {
        P1DiscCL::GetGradients( G,det,*sit);
        absdet= fabs( det);
        for(int i=0; i<4; ++i) 
        {
            for(int j=0; j<=i; ++j) 
            {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )/6.0*absdet;
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

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::InitVel(VelVecDescCL* vec, vector_instat_fun_ptr LsgVel, double t0) const
{
    VectorCL& lsgvel= vec->Data;
    Uint lvl        = vec->GetLevel(),
         vidx       = vec->RowIdx->GetIdx();
    
    SVectorCL<3> tmp;
    
    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
            tmp= LsgVel(sit->GetCoord(), t0);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= tmp[i];
        }
    }
    
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangEdgeBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
            tmp= LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t0);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= tmp[i];
        }
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
    
    const Uint lvl         = A->GetRowLevel(),
               vidx        = A->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupSystem1: " << num_unks_vel << " vels. ";

    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    Quad2CL<double> rho, mu_Re, Phi, kreuzterm;
        
    SMatrixCL<3,3> T;
    
    double coupA[10][10], coupM[10][10];
    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::DiscSolCL ls= lset.GetSolution();

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
    
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            rhs.val[i]= _Coeff.f( sit->GetVertex(i)->GetCoord(), t);
            Phi.val[i]= ls.val( *sit->GetVertex(i));
        }
        for (int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        rhs.val[4]= _Coeff.f( GetBaryCenter( *sit), t);
        Phi.val[4]= ls.val( *sit, 0.25, 0.25, 0.25);

        // rho = rho( Phi),    mu_Re= mu( Phi)/Re
        rho=   Phi;     rho.apply( _Coeff.rho);
        mu_Re= Phi;     mu_Re.apply( _Coeff.mu);     mu_Re*= 1./_Coeff.Re;

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( _Coeff.g)*rho;
        
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const Quad2CL<double> dotGrad= dot( Grad[i], Grad[j]) * mu_Re;
                coupA[i][j]= coupA[j][i]= dotGrad.quad( absdet);
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j, absdet);
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        mA( Numb[i],   Numb[j]  )+= coupA[j][i];
                        mA( Numb[i]+1, Numb[j]+1)+= coupA[j][i];
                        mA( Numb[i]+2, Numb[j]+2)+= coupA[j][i];
                        for (int k=0; k<3; ++k)
                            for (int l=0; l<3; ++l)
                            {
                                // kreuzterm = \int mu/Re * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm.val[m]= Grad[i].val[m][l] * Grad[j].val[m][k] * mu_Re.val[m];

                                mA( Numb[i]+k, Numb[j]+l)+= kreuzterm.quad( absdet);
                            }
                        mM( Numb[i],   Numb[j]  )+= coupM[j][i];
                        mM( Numb[i]+1, Numb[j]+1)+= coupM[j][i];
                        mM( Numb[i]+2, Numb[j]+2)+= coupM[j][i];
                    }
                    else // put coupling on rhs
                    {
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j), t)
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4), t);
                        const double cA= coupA[j][i],
                                     cM= coupM[j][i];
                        for (int k=0; k<3; ++k)
                        {
                            cplA->Data[Numb[i]+k]-= cA*tmp[k];
                            
                            for (int l=0; l<3; ++l)
                            {
                                // kreuzterm = \int mu/Re * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm.val[m]= Grad[i].val[m][l] * Grad[j].val[m][k] * mu_Re.val[m];
                                cplA->Data[Numb[i]+k]-= kreuzterm.quad( absdet)*tmp[l];
                            }
                        }
                        cplM->Data[Numb[i]  ]-= cM*tmp[0];
                        cplM->Data[Numb[i]+1]-= cM*tmp[1];
                        cplM->Data[Numb[i]+2]-= cM*tmp[2];
                    }
                }
                tmp= rhs.quadP2( i, absdet);
                b->Data[Numb[i]  ]+= tmp[0];
                b->Data[Numb[i]+1]+= tmp[1];
                b->Data[Numb[i]+2]+= tmp[2];

                if ( i<4 ? _BndData.Vel.IsOnNeuBnd(*sit->GetVertex(i))
                         : _BndData.Vel.IsOnNeuBnd(*sit->GetEdge(i-4)) ) // vert/edge i is on Neumann boundary
                {
                    Uint face;
                    const int num_faces= i<4 ? 3 : 2;
                    for (int f=0; f < num_faces; ++f)
                    {
                        face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);
                        if ( sit->IsBndSeg(face))
                        {   // TODO: FIXME: Hier muss doch eigentlich eine 2D-Integrationsformel fuer P2-Elemente stehen, oder?
/*
                            tmp= Quad2D(*sit, face, i, _BndData.Vel.GetSegData(sit->GetBndIdx(face)).GetBndFun(), t);
                            a[Numb[i]]+=          tmp[0];
                            a[Numb[i]+stride]+=   tmp[1];
                            a[Numb[i]+2*stride]+= tmp[2];
*/                            
                        }
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
void InstatStokes2PhaseP2P1CL<Coeff>::SetupMatrices1( MatDescCL* A,
    MatDescCL* M, const LevelsetP2CL& lset, double t) const
// Set up matrices A, M (depending on phase bnd)
{
    const IdxT num_unks_vel= A->RowIdx->NumUnknowns;

    MatrixBuilderCL mA( &A->Data, num_unks_vel, num_unks_vel),
                    mM( &M->Data, num_unks_vel, num_unks_vel);
    
    const Uint lvl= A->GetRowLevel(),
               vidx= A->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupMatrices1: " << num_unks_vel << " vels. ";

    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    Quad2CL<double> rho, mu_Re, Phi, kreuzterm;
        
    SMatrixCL<3,3> T;
    
    double coupA[10][10], coupM[10][10];
    double lsarray[10];
    double det, absdet;
    LevelsetP2CL::DiscSolCL ls= lset.GetSolution();
    
    for (MultiGridCL::const_TriangTetraIteratorCL 
            sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl),
            send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
            sit != send; ++sit) {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        P2DiscCL::GetGradientsOnRef( GradRef);
        if (ls.GetLevel()!=lvl)
            RestrictP2( *sit, ls, lsarray);
        else
            ls.GetDoF( *sit, lsarray);
    
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            rhs.val[i]= _Coeff.f( sit->GetVertex(i)->GetCoord(), t);
            Phi.val[i]= lsarray[i];
        }
        for (int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        rhs.val[4]= _Coeff.f( GetBaryCenter( *sit), t);
        Phi.val[4]= P2( lsarray, 1.0/*type-dummy*/, 0.25, 0.25, 0.25);

        // rho = rho( Phi),    mu_Re= mu( Phi)/Re
        rho=   Phi;     rho.apply( _Coeff.rho);
        mu_Re= Phi;     mu_Re.apply( _Coeff.mu);     mu_Re*= 1./_Coeff.Re;

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( _Coeff.g)*rho;
        
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const Quad2CL<double> dotGrad= dot( Grad[i], Grad[j]) * mu_Re;
                coupA[i][j]= coupA[j][i]= dotGrad.quad( absdet);
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j, absdet);
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        mA( Numb[i],   Numb[j]  )+= coupA[j][i];
                        mA( Numb[i]+1, Numb[j]+1)+= coupA[j][i];
                        mA( Numb[i]+2, Numb[j]+2)+= coupA[j][i];
                        for (int k=0; k<3; ++k)
                            for (int l=0; l<3; ++l)
                            {
                                // kreuzterm = \int mu/Re * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm.val[m]= Grad[i].val[m][l] * Grad[j].val[m][k] * mu_Re.val[m];

                                mA( Numb[i]+k, Numb[j]+l)+= kreuzterm.quad( absdet);
                            }
                        mM( Numb[i],   Numb[j]  )+= coupM[j][i];
                        mM( Numb[i]+1, Numb[j]+1)+= coupM[j][i];
                        mM( Numb[i]+2, Numb[j]+2)+= coupM[j][i];
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
void InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2( MatDescCL* B, VecDescCL* c, double t) const
// Set up matrix B and rhs c (independent of phase bnd, but c is time dependent)
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    c->Clear();
    
    const Uint lvl         = B->GetRowLevel();
    const Uint pidx        = B->RowIdx->GetIdx(),
               vidx        = B->ColIdx->GetIdx();

    IdxT Numb[10], prNumb[4];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupSystem2: " << B->RowIdx->NumUnknowns << " prs. ";

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
        
    SMatrixCL<3,3> T;
    
    double det, absdet;
    Point3DCL tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
    
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for (int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for (int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel)
        {
            if (!IsOnDirBnd[vel])
                for(int pr=0; pr<4; ++pr)
                {
                    tmp= Grad[vel].quadP1( pr, absdet);
                    mB( prNumb[pr], Numb[vel])  -=  tmp[0];
                    mB( prNumb[pr], Numb[vel]+1)-=  tmp[1];
                    mB( prNumb[pr], Numb[vel]+2)-=  tmp[2];
                } 
            else // put coupling on rhs
            {
                tmp= vel<4 ? _BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : _BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr, absdet), tmp);
            }
        }
    }

    mB.Build();
    std::cerr << B->Data.num_nonzeros() << " nonzeros in B!" << std::endl;
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupRhs2( VecDescCL* c, double t) const
// Set up rhs c (time dependent)
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
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        // collect some information about the edges and verts of the tetra
        // and save it in prNumb and IsOnDirBnd
        for (int i=0; i<4; ++i)
        {
            IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) );
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for (int i=0; i<6; ++i)
            IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) );

        GetTrafoTr( T, det, *sit);
        absdet= fabs( det);
    
        // b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel)
        {
            if (IsOnDirBnd[vel])
            { // put coupling on rhs
                P2DiscCL::GetGradient( Grad_vel, GradRef[vel], T);
                tmp= vel<4 ? _BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : _BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad_vel.quadP1( pr, absdet), tmp);
            }
        }
    }
}


} // end of namespace DROPS
