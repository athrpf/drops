//**************************************************************************
// File:    instatstokes2phase.tpp                                         *
// Content: classes that constitute the 2-phase stokes-problem             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "num/discretize.h"

namespace DROPS
{

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::CreateNumberingPr(Uint level, IdxDescCL* idx)
// used for numbering of the Unknowns depending on the index IdxDesc[idxnum].
// sets up the description of the index idxnum in IdxDesc[idxnum],
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, thatIdxDesc[idxnr].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx->TriangLevel = level;
    idx->NumUnknowns = 0;

    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL

    // allocate space for indices; number unknowns in TriangLevel level; "true" generates unknowns on
    // the boundary
    CreateNumbOnVertex( idxnum, idx->NumUnknowns, idx->NumUnknownsVertex,
                        _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level),
                        _BndData.Pr );
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::CreateNumberingVel(Uint level, IdxDescCL* idx)
// used for numbering of the Unknowns depending on the index IdxDesc[idxnum].
// sets up the description of the index idxnum in IdxDesc[idxnum],
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, thatIdxDesc[idxnr].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx->TriangLevel = level;
    idx->NumUnknowns = 0;

    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL

    // allocate space for indices; number unknowns in TriangLevel level
    CreateNumbOnVertex( idxnum, idx->NumUnknowns, idx->NumUnknownsVertex,
                        _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level),
                        _BndData.Vel );
    CreateNumbOnEdge( idxnum, idx->NumUnknowns, idx->NumUnknownsEdge,
                      _MG.GetTriangEdgeBegin(level), _MG.GetTriangEdgeEnd(level),
                      _BndData.Vel );
}


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

template <class Coeff>
double H_eps( double s)
{
    if (s <= -Coeff::eps)
        return 0;
    if (s >= Coeff::eps)
        return 1;
    // -eps < s < eps
    s/= Coeff::eps;
    const double s2= s*s, s3= s2*s;
    return 0.5 + 1.40625*s - 1.5625*s3 + 0.65625*s2*s3;
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::SetupPrMass(MatDescCL* matM) const
// Sets up the mass matrix for the pressure
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;

    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl    = matM->RowIdx->TriangLevel;
    const Uint pidx   = matM->RowIdx->GetIdx();

    IdxT prNumb[4];

    // compute all couplings between HatFunctions on verts:
    // I( i, j) = int ( psi_i*psi_j, T_ref) * absdet
    const double coupl_ii= 1./60.,  // = 1/120 + 2/15 * 1/16,
                 coupl_ij= 1./120.; // =         2/15 * 1/16;

    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6.;
        
        for(int i=0; i<4; ++i)
        {
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= (i!=j ? coupl_ij : coupl_ii) * absdet;
                    //mu_Re.quadP1(i,j) * absdet;
    }
    M_pr.Build();
}

template <class Coeff>
void InstatStokes2PhaseP2P1CL<Coeff>::InitVel(VelVecDescCL* vec, vector_instat_fun_ptr LsgVel, double t0) const
{
    VectorCL& lsgvel= vec->Data;
    Uint lvl        = vec->RowIdx->TriangLevel,
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
void InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem1( MatDescCL* A, MatDescCL* M, VecDescCL* b, VecDescCL* cplM, const LevelsetP2CL<_self>& lset, double t) const
// Set up matrices A, M and rhs b (depending on phase bnd)
{
    const IdxT num_unks_vel= A->RowIdx->NumUnknowns;

    MatrixBuilderCL mA( &A->Data, num_unks_vel, num_unks_vel),
                    mM( &M->Data, num_unks_vel, num_unks_vel);
    b->Clear();
    cplM->Clear();
    
    const Uint lvl         = A->RowIdx->TriangLevel,
               vidx        = A->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupSystem1: " << num_unks_vel << " vels." << std::endl;                            

    Quad2CL<Point3DCL> Grad[10], GradRef[10], rhs;
    Quad2CL<double> rho, mu_Re, H, kreuzterm;
        
    SMatrixCL<3,3> T;
    
    double coupA[10][10], coupM[10][10];
    double det, absdet;
    Point3DCL tmp;
    typename LevelsetP2CL<_self>::DiscSolCL ls= lset.GetSolution();

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
            H.val[i]= ls.val( *sit->GetVertex(i));
        }
        for (int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        rhs.val[4]= _Coeff.f( GetBaryCenter( *sit), t);
        H.val[4]= ls.val( *sit, 0.25, 0.25, 0.25);
        H.apply( H_eps<Coeff>);

        // rho = rho1 + (rho2-rho1)*H
        rho= H * (_Coeff.rho2 - _Coeff.rho1);
        rho+= _Coeff.rho1;

        // mu_Re = (mu1 + (mu2-mu1)*H) / Re
        mu_Re= H * ((_Coeff.mu2 - _Coeff.mu1) / _Coeff.Re);
        mu_Re+= _Coeff.mu1 / _Coeff.Re;

        // rhs = f + rho*g
        rhs+= Quad2CL<Point3DCL>( _Coeff.g)*rho;
        
        // compute all couplings between HatFunctions on edges and verts
        for (int i=0; i<10; ++i)
            for (int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const Quad2CL<double> dotGrad= dot( Grad[i], Grad[j]) * mu_Re;
                coupA[i][j]= coupA[j][i]= dotGrad.quad()*absdet;
                coupM[i][j]= coupM[j][i]= rho.quadP2(i,j)*absdet;
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

                                mA( Numb[i]+k, Numb[j]+l)+= kreuzterm.quad()*absdet;
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
                            b->Data[Numb[i]+k]-= cA*tmp[k];
                            
                            for (int l=0; l<3; ++l)
                            {
                                // kreuzterm = \int mu/Re * (dphi_i / dx_l) * (dphi_j / dx_k)
                                for (size_t m=0; m<kreuzterm.size();  ++m)
                                    kreuzterm.val[m]= Grad[i].val[m][l] * Grad[j].val[m][k] * mu_Re.val[m];
                                b->Data[Numb[i]+k]-= kreuzterm.quad()*absdet*tmp[l];
                            }
                        }
                        cplM->Data[Numb[i]  ]-= cM*tmp[0];
                        cplM->Data[Numb[i]+1]-= cM*tmp[1];
                        cplM->Data[Numb[i]+2]-= cM*tmp[2];
                    }
                }
                tmp= rhs.quadP2( i)*absdet;
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
void InstatStokes2PhaseP2P1CL<Coeff>::SetupSystem2( MatDescCL* B, VecDescCL* c, double t) const
// Set up matrix B and rhs c (independent of phase bnd, but c is time dependent)
{
    MatrixBuilderCL mB( &B->Data, B->RowIdx->NumUnknowns, B->ColIdx->NumUnknowns);
    c->Clear();
    
    const Uint lvl         = B->RowIdx->TriangLevel;
    const Uint pidx        = B->RowIdx->GetIdx(),
               vidx        = B->ColIdx->GetIdx();

    IdxT Numb[10], prNumb[4];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupSystem2: " << B->RowIdx->NumUnknowns << " prs." << std::endl;                            

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
                    tmp= Grad[vel].quadP1( pr)*absdet;
                    mB( prNumb[pr], Numb[vel])  -=  tmp[0];
                    mB( prNumb[pr], Numb[vel]+1)-=  tmp[1];
                    mB( prNumb[pr], Numb[vel]+2)-=  tmp[2];
                } 
            else // put coupling on rhs
            {
                tmp= vel<4 ? _BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : _BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr)
                    c->Data[ prNumb[pr]]+= inner_prod( Grad[vel].quadP1( pr), tmp)*absdet;
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
    
    const Uint lvl         = c->RowIdx->TriangLevel;
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
                    c->Data[ prNumb[pr]]+= inner_prod( Grad_vel.quadP1( pr), tmp)*absdet;
            }
        }
    }
}


} // end of namespace DROPS
