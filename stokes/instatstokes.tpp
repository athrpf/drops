//**************************************************************************
// File:    instatstokes.tpp                                               *
// Content: classes that constitute the stokes-problem                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Sep, 13 2001                                           *
//**************************************************************************

#include "num/discretize.h"

namespace DROPS
{

template <class Coeff>
void InstatStokesP2P1CL<Coeff>::CreateNumberingPr(Uint level, IdxDescCL* idx)
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
void InstatStokesP2P1CL<Coeff>::CreateNumberingVel(Uint level, IdxDescCL* idx)
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
void InstatStokesP2P1CL<Coeff>::DeleteNumberingVel(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetAllVertexBegin(level), _MG.GetAllVertexEnd(level) );
    DeleteNumbOnSimplex( idxnum, _MG.GetAllEdgeBegin(level), _MG.GetAllEdgeEnd(level) );
}

template <class Coeff>
void InstatStokesP2P1CL<Coeff>::DeleteNumberingPr(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetAllVertexBegin(level), _MG.GetAllVertexEnd(level) );
}




template <class Coeff>
void InstatStokesP2P1CL<Coeff>::GetDiscError(vector_instat_fun_ptr LsgVel, scalar_instat_fun_ptr LsgPr, double t) const
{
    Uint lvl= A.RowIdx->TriangLevel,
        vidx= A.RowIdx->GetIdx(),
        pidx= B.RowIdx->GetIdx();
    VectorCL lsgvel(A.RowIdx->NumUnknowns);
    VectorCL lsgpr( B.RowIdx->NumUnknowns);
    SVectorCL<3> tmp;

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
            tmp= LsgVel(sit->GetCoord(), t);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= tmp[i];
        }
    }
    
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangEdgeBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
            tmp= LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= tmp[i];
        }
    }
    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
        lsgpr[sit->Unknowns(pidx)]= LsgPr(sit->GetCoord(), t);

    std::cerr << "discretization error to check the system (x,y = continuous solution): "<<std::endl;
    VectorCL res= A.Data*lsgvel + transp_mul(B.Data,lsgpr)-b.Data; 
    std::cerr <<"|| Ax + BTy - f || = "<< res.norm()<<", max "<<res.supnorm()<<std::endl;
    VectorCL resB= B.Data*lsgvel - c.Data; 
    std::cerr <<"|| Bx - g || = "<< resB.norm()<<", max "<<resB.supnorm()<<std::endl;
}



inline SVectorCL<3> Quad( const TetraCL& tetra, vector_instat_fun_ptr coeff, int i, double t)
{
    SVectorCL<3> f[5];
    
    if (i<4) // hat function on vert
    {
        f[0]= coeff( tetra.GetVertex(i)->GetCoord(), t);
        for (int k=0, l=1; k<4; ++k)
            if (k!=i) f[l++]= coeff( tetra.GetVertex(k)->GetCoord(), t);
        f[4]= coeff( GetBaryCenter(tetra), t);
        return f[0]/504. - (f[1] + f[2] + f[3])/1260. - f[4]/126.;
    }
    else  // hat function on edge
    {
        const double ve= 4./945.,  // coeff for verts of edge
                     vn= -1./756.,  // coeff for other verts
                     vs= 26./945.;   // coeff for barycenter
        double a[4];
        a[VertOfEdge(i-4,0)]= a[VertOfEdge(i-4,1)]= ve;
        a[VertOfEdge(OppEdge(i-4),0)]= a[VertOfEdge(OppEdge(i-4),1)]= vn;

        SVectorCL<3> sum= vs * coeff( GetBaryCenter(tetra), t);
        for(int k=0; k<4; ++k)
            sum+= a[k] * coeff( tetra.GetVertex(k)->GetCoord(), t);

        return sum;
    }
}

inline InstatStokesVelBndDataCL::bnd_type Quad2D( const TetraCL& tetra, Uint face, Uint vert, 
                                                  InstatStokesVelBndDataCL::bnd_val_fun bfun, double t)
// Integrate bnd_val * phi_vert over face
{
    const VertexCL* v[3];
    
    v[0]= tetra.GetVertex(vert);
    for (int i=0, k=1; i<3; ++i)
    {
        if (VertOfFace(face,i)!=vert)
            v[k++]= tetra.GetVertex( VertOfFace(face,i) );
    }
    const InstatStokesVelBndDataCL::bnd_type f0= bfun(v[0]->GetCoord(), t);
    const InstatStokesVelBndDataCL::bnd_type f1= bfun(v[1]->GetCoord(), t) +  bfun( v[2]->GetCoord(), t);
    const InstatStokesVelBndDataCL::bnd_type f2= bfun((v[0]->GetCoord() + v[1]->GetCoord() + v[2]->GetCoord())/3.0, t);    //Barycenter of Face
    const double absdet= FuncDet2D(v[1]->GetCoord() - v[0]->GetCoord(), v[2]->GetCoord() - v[0]->GetCoord());
    return (11./240.*f0 + 1./240.*f1 + 9./80.*f2) * absdet;
}


template <class Coeff>
void InstatStokesP2P1CL<Coeff>::SetupSystem( MatDescCL* matA, VelVecDescCL* vecA, 
                                             MatDescCL* matB, VelVecDescCL* vecB, double t) const
// Sets up the stiffness matrices and right hand sides
{
    vecA->Clear();
    vecB->Clear();

    const IdxT num_unks_vel= matA->RowIdx->NumUnknowns;
    const IdxT num_unks_pr=  matB->RowIdx->NumUnknowns;

    MatrixBuilderCL A(&matA->Data, num_unks_vel, num_unks_vel), 
                    B(&matB->Data, num_unks_pr,  num_unks_vel);
    VelVecDescCL& b   = *vecA;
    VelVecDescCL& c   = *vecB;
    const Uint lvl    = matA->RowIdx->TriangLevel;
    const Uint vidx   = matA->RowIdx->GetIdx(),
               pidx   = matB->RowIdx->GetIdx();

    IdxT Numb[10], prNumb[4];
    bool IsOnDirBnd[10];
    
    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns
                            
    std::cerr << "entering SetupSystem: " <<num_unks_vel<<" vels, "<<num_unks_pr<<" prs"<< std::endl;                            

    // fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup[10][10];
    double det, absdet;
    SVectorCL<3> tmp;

    GetGradientsOnRef(GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }

        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= _Coeff.nu * QuadGrad( Grad, i, j)*absdet;
                coup[i][j]+= Quad(*sit, &_Coeff.q, i, j)*absdet;
                coup[j][i]= coup[i][j];
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        A(Numb[i],          Numb[j])+=          coup[j][i]; 
                        A(Numb[i]+stride,   Numb[j]+stride)+=   coup[j][i]; 
                        A(Numb[i]+2*stride, Numb[j]+2*stride)+= coup[j][i]; 
                    }
                    else // coupling with vert/edge j on right-hand-side
                    {
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j), t)
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4), t);
                        b.Data[Numb[i]]-=          coup[j][i] * tmp[0];
                        b.Data[Numb[i]+stride]-=   coup[j][i] * tmp[1];
                        b.Data[Numb[i]+2*stride]-= coup[j][i] * tmp[2];
                    }
                }
                tmp= Quad(*sit, &_Coeff.f, i, t)*absdet;
                b.Data[Numb[i]]+=          tmp[0];
                b.Data[Numb[i]+stride]+=   tmp[1];
                b.Data[Numb[i]+2*stride]+= tmp[2];

                if ( i<4 ? _BndData.Vel.IsOnNeuBnd(*sit->GetVertex(i))
                         : _BndData.Vel.IsOnNeuBnd(*sit->GetEdge(i-4)) ) // vert/edge i is on Neumann boundary
                {
                    Uint face;
                    for (int f=0; f < 3; ++f)
                    {
                        face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);
                        if ( sit->IsBndSeg(face))
                        {
                            tmp= Quad2D(*sit, face, i, _BndData.Vel.GetSegData(sit->GetBndIdx(face)).GetBndFun(), t);
                            b.Data[Numb[i]]+=          tmp[0];
                            b.Data[Numb[i]+stride]+=   tmp[1];
                            b.Data[Numb[i]+2*stride]+= tmp[2];
                        }
                    }
                }
            }
            
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel)
        {
            if (!IsOnDirBnd[vel])
                for(int pr=0; pr<4; ++pr)
                {
                    // numeric integration is exact: psi_i * div( phi_j) is of degree 2 !
                    // hint: psi_i( barycenter ) = 0.25
                    B(prNumb[pr],Numb[vel])-=           (Grad[vel](0,pr) / 120. + Grad[vel](0,4)*0.25 * 2./15. )*absdet;
                    B(prNumb[pr],Numb[vel]+stride)-=    (Grad[vel](1,pr) / 120. + Grad[vel](1,4)*0.25 * 2./15. )*absdet;
                    B(prNumb[pr],Numb[vel]+2*stride)-=  (Grad[vel](2,pr) / 120. + Grad[vel](2,4)*0.25 * 2./15. )*absdet;
                }
            else // put coupling on rhs
            {
                tmp= vel<4 ? _BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), t)
                           : _BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), t);
                for(int pr=0; pr<4; ++pr)
                {
                    // numeric integration is exact: psi_i * div( phi_j) is of degree 2 !
                    c.Data[prNumb[pr]]+= (Grad[vel](0,pr) / 120. + Grad[vel](0,4)*0.25 * 2./15. )*absdet*tmp[0]
                                         +(Grad[vel](1,pr) / 120. + Grad[vel](1,4)*0.25 * 2./15. )*absdet*tmp[1]
                                         +(Grad[vel](2,pr) / 120. + Grad[vel](2,4)*0.25 * 2./15. )*absdet*tmp[2];
                }
            }
        }
    }
    std::cerr << "done: value part fill" << std::endl;

    A.Build();
    B.Build();
    std::cerr << matA->Data.num_nonzeros() << " nonzeros in A, "
              << matB->Data.num_nonzeros() << " nonzeros in B! " << std::endl;
}


template <class Coeff>
void InstatStokesP2P1CL<Coeff>::SetupPrMass(MatDescCL* matM) const
// Sets up the mass matrix for the pressure
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns;

    MatrixBuilderCL M_pr(&matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl    = matM->RowIdx->TriangLevel;
    const Uint pidx   = matM->RowIdx->GetIdx();

    IdxT prNumb[4];

    // compute all couplings between HatFunctions on verts:
    // I( i, j) = int ( psi_i*psi_j, T_ref) * absdet
    const double coupl_ii= 1./120. + 2./15./16.,
                 coupl_ij=           2./15./16.;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6.;
        
        for(int i=0; i<4; ++i)
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                    M_pr( prNumb[i], prNumb[j])+= (i==j ? coupl_ii : coupl_ij) * absdet;
    }
    M_pr.Build();
}



inline double OneFct( const Point3DCL&) { return 1.;}


template <class Coeff>
void InstatStokesP2P1CL<Coeff>::SetupInstatSystem(MatDescCL* matA, MatDescCL* matB, MatDescCL* matI) const
// Sets up the stiffness matrices and right hand sides
{
    const IdxT num_unks_vel= matA->RowIdx->NumUnknowns;
    const IdxT num_unks_pr=  matB->RowIdx->NumUnknowns;

    MatrixBuilderCL A(&matA->Data, num_unks_vel, num_unks_vel), 
                    B(&matB->Data, num_unks_pr,  num_unks_vel),
                    I(&matI->Data, num_unks_vel, num_unks_vel);

    const Uint lvl    = matA->RowIdx->TriangLevel;
    const Uint vidx   = matA->RowIdx->GetIdx(),
               pidx   = matB->RowIdx->GetIdx();

    IdxT Numb[10], prNumb[4];
    bool IsOnDirBnd[10];
    
    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns
                            
    std::cerr << "entering SetupSystem: " <<num_unks_vel<<" vels, "<<num_unks_pr<<" prs"<< std::endl;                            

    // fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup[10][10];
    double det, absdet;

    GetGradientsOnRef(GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }

        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= _Coeff.nu * QuadGrad( Grad, i, j)*absdet;
//                coup[i][j]+= Quad(*sit, &_Coeff.q, i, j)*absdet;
                coup[j][i]= coup[i][j];
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        A(Numb[i],          Numb[j])+=          coup[j][i]; 
                        A(Numb[i]+stride,   Numb[j]+stride)+=   coup[j][i]; 
                        A(Numb[i]+2*stride, Numb[j]+2*stride)+= coup[j][i];
                        I(Numb[i],          Numb[j])
                        = I(Numb[i]+stride,   Numb[j]+stride)
                        = I(Numb[i]+2*stride, Numb[j]+2*stride)+= Quad( *sit, &OneFct, i, j)*absdet;
                    }
                }

            }
            
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel)
        {
            if (!IsOnDirBnd[vel])
                for(int pr=0; pr<4; ++pr)
                {
                    // numeric integration is exact: psi_i * div( phi_j) is of degree 2 !
                    // hint: psi_i( barycenter ) = 0.25
                    B(prNumb[pr],Numb[vel])-=           (Grad[vel](0,pr) / 120. + Grad[vel](0,4)*0.25 * 2./15. )*absdet;
                    B(prNumb[pr],Numb[vel]+stride)-=    (Grad[vel](1,pr) / 120. + Grad[vel](1,4)*0.25 * 2./15. )*absdet;
                    B(prNumb[pr],Numb[vel]+2*stride)-=  (Grad[vel](2,pr) / 120. + Grad[vel](2,4)*0.25 * 2./15. )*absdet;
                } // else put coupling on rhs
        }
    }
    std::cerr << "done: value part fill" << std::endl;

    A.Build();
    B.Build();
    I.Build();
    std::cerr << matA->Data.num_nonzeros() << " nonzeros in A, "
              << matB->Data.num_nonzeros() << " nonzeros in B, "
              << matI->Data.num_nonzeros() << " nonzeros in I! " << std::endl;
}

template <class Coeff>
void InstatStokesP2P1CL<Coeff>::SetupInstatRhs( VelVecDescCL* vecA, VelVecDescCL* vecB, 
                                                VelVecDescCL* vecI, double tA, 
                                                VelVecDescCL* vecf, double tf) const
// Sets up the couplings with hat fcts on bnd (for matrices A, I) and discretizes the PDE coeff f(t)
{
    vecA->Clear();
    vecB->Clear();
    vecf->Clear();
    vecI->Clear();

    VectorCL& a    = vecA->Data;
    VectorCL& c    = vecB->Data;
    VectorCL& f    = vecf->Data;
    VectorCL& id   = vecI->Data;
    const Uint lvl = vecA->RowIdx->TriangLevel;
    const Uint vidx= vecA->RowIdx->GetIdx(),
               pidx= vecB->RowIdx->GetIdx();

    IdxT Numb[10], prNumb[4];
    bool IsOnDirBnd[10];
    
    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns
                            
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup[10][10];
    double det, absdet;
    SVectorCL<3> tmp;

    GetGradientsOnRef(GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }

        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= _Coeff.nu * QuadGrad( Grad, i, j)*absdet;
//                coup[i][j]+= Quad(*sit, &_Coeff.q, i, j)*absdet;
                coup[j][i]= coup[i][j];
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (IsOnDirBnd[j]) // vert/edge j is on a Dirichlet boundary
                    { // coupling with vert/edge j on right-hand-side
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j), tA)
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4), tA);
                        a[Numb[i]]-=          coup[j][i] * tmp[0];
                        a[Numb[i]+stride]-=   coup[j][i] * tmp[1];
                        a[Numb[i]+2*stride]-= coup[j][i] * tmp[2];

                        const double val= Quad( *sit, &OneFct, i, j)*absdet;
                        id[Numb[i]]-=          val*tmp[0];
                        id[Numb[i]+stride]-=   val*tmp[1];
                        id[Numb[i]+2*stride]-= val*tmp[2];
                    }
                }
                tmp= Quad(*sit, &_Coeff.f, i, tf)*absdet;
                f[Numb[i]]+=          tmp[0];
                f[Numb[i]+stride]+=   tmp[1];
                f[Numb[i]+2*stride]+= tmp[2];

                if ( i<4 ? _BndData.Vel.IsOnNeuBnd(*sit->GetVertex(i))
                         : _BndData.Vel.IsOnNeuBnd(*sit->GetEdge(i-4)) ) // vert/edge i is on Neumann boundary
                {
                    Uint face;
                    for (int f=0; f < 3; ++f)
                    {
                        face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);
                        if ( sit->IsBndSeg(face))
                        {
                            tmp= Quad2D(*sit, face, i, _BndData.Vel.GetSegData(sit->GetBndIdx(face)).GetBndFun(), tA);
                            a[Numb[i]]+=          tmp[0];
                            a[Numb[i]+stride]+=   tmp[1];
                            a[Numb[i]+2*stride]+= tmp[2];
                        }
                    }
                }
            }
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel)
        {
            if (IsOnDirBnd[vel])
            {    // put coupling on rhs
                tmp= vel<4 ? _BndData.Vel.GetDirBndValue( *sit->GetVertex(vel), tA)
                           : _BndData.Vel.GetDirBndValue( *sit->GetEdge(vel-4), tA);
                for(int pr=0; pr<4; ++pr)
                {
                    // numeric integration is exact: psi_i * div( phi_j) is of degree 2 !
                    c[prNumb[pr]]+=( (Grad[vel](0,pr) / 120. + Grad[vel](0,4)*0.25 * 2./15. )*tmp[0]
                                    +(Grad[vel](1,pr) / 120. + Grad[vel](1,4)*0.25 * 2./15. )*tmp[1]
                                    +(Grad[vel](2,pr) / 120. + Grad[vel](2,4)*0.25 * 2./15. )*tmp[2] )*absdet;
                }
            }
        }
    }
}


template <class Coeff>
void InstatStokesP2P1CL<Coeff>::InitVel(VelVecDescCL* vec, vector_instat_fun_ptr LsgVel, double t0) const
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
void InstatStokesP2P1CL<Coeff>::CheckSolution(const VelVecDescCL* lsgvel, const VecDescCL* lsgpr, 
                                 vector_instat_fun_ptr LsgVel, scalar_instat_fun_ptr LsgPr, double t) const
{
    double diff, maxdiff=0, norm2= 0;
    Uint lvl=lsgvel->RowIdx->TriangLevel,
         vidx=lsgvel->RowIdx->GetIdx();
    
    VectorCL res1= A.Data*lsgvel->Data + transp_mul( B.Data, lsgpr->Data ) - b.Data;
    VectorCL res2= B.Data*lsgvel->Data - c.Data;

    std::cerr << "\nChecken der Loesung...\n";    
    std::cerr << "|| Ax + BTy - F || = " << res1.norm() << ", max. " << res1.supnorm() << std::endl;
    std::cerr << "||       Bx - G || = " << res2.norm() << ", max. " << res2.supnorm() << std::endl<<std::endl;
    
    DiscVelSolCL vel(lsgvel, &_BndData.Vel, &_MG, t);
    double L1_div= 0, L2_div= 0;
    SMatrixCL<3,3> T;
    double det, absdet;

    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double div[5]= {0., 0., 0., 0., 0.};   // Divergenz in den Verts und im BaryCenter
        GetTrafoTr(T,det,*sit);
        absdet= fabs(det);
        for(Uint i= 0; i<10; ++i)
        {
            SVectorCL<3> value= i<4 ? vel.val(*sit->GetVertex(i))
                                    : vel.val(*sit->GetEdge(i-4));
            div[0]+= inner_prod(T*FE_P2CL::DHRef(i,0,0,0),value);
            div[1]+= inner_prod(T*FE_P2CL::DHRef(i,1,0,0),value);
            div[2]+= inner_prod(T*FE_P2CL::DHRef(i,0,1,0),value);
            div[3]+= inner_prod(T*FE_P2CL::DHRef(i,0,0,1),value);
            div[4]+= inner_prod(T*FE_P2CL::DHRef(i,0.25,0.25,0.25),value);
        }
        L1_div+= ( (fabs(div[0])+fabs(div[1])+fabs(div[2])+fabs(div[3]))/120 + fabs(div[4])*2./15. ) * absdet;
        L2_div+= ( (div[0]*div[0]+div[1]*div[1]+div[2]*div[2]+div[3]*div[3])/120 + div[4]*div[4]*2./15. ) * absdet;
    }
    L2_div= ::sqrt(L2_div);
    std::cerr << "|| div x ||_L1 = " << L1_div << std::endl;
    std::cerr << "|| div x ||_L2 = " << L2_div << std::endl << std::endl;

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
           {
               diff= fabs( LsgVel(sit->GetCoord(), t)[i] - lsgvel->Data[sit->Unknowns(vidx)+i]);
               norm2+= diff*diff;
               if (diff>maxdiff)
               {
                   maxdiff= diff;
               }
           }
        }
    }
    
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangEdgeBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
           {
               diff= fabs( LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t)[i] - lsgvel->Data[sit->Unknowns(vidx)+i]);
               norm2+= diff*diff;
               if (diff>maxdiff)
               {
                   maxdiff= diff;
               }
           }
        }
    }
    norm2= ::sqrt(norm2 / lsgvel->Data.size());

    Point3DCL L1_vel(0.0), L2_vel(0.0);
    for(MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send= const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
        sit != send; ++sit)
    {
	Point3DCL sum(0.0), diff, Diff[5];
	for(int i=0; i<4; ++i)
	{
	    Diff[i]= diff= LsgVel(sit->GetVertex(i)->GetCoord(), t) - vel.val(*sit->GetVertex(i));
	    diff[0]= fabs(diff[0]); diff[1]= fabs(diff[1]); diff[2]= fabs(diff[2]);
	    sum+= diff;
	}
	sum/= 120;
	Diff[4]= diff= LsgVel(GetBaryCenter(*sit), t) - vel.val(*sit, 0.25, 0.25, 0.25);
	diff[0]= fabs(diff[0]); diff[1]= fabs(diff[1]); diff[2]= fabs(diff[2]);
	sum+= diff*2./15.;
	sum*= sit->GetVolume()*6;
	L1_vel+= sum;

	for(int i=0; i<10; ++i)
	{
	    sum= Quad(Diff, i)*sit->GetVolume()*6;
	    diff= i<4 ? Diff[i] : LsgVel( (sit->GetEdge(i-4)->GetVertex(0)->GetCoord() 
                                         + sit->GetEdge(i-4)->GetVertex(1)->GetCoord() )/2, t) 
                                - vel.val(*sit->GetEdge(i-4));
	    sum[0]*= diff[0]; sum[1]*= diff[1]; sum[2]*= diff[2];
	    L2_vel+= sum;
	}
    }
    L2_vel= sqrt(L2_vel);
    std::cerr << "Geschwindigkeit: Abweichung von der tatsaechlichen Loesung:\n"
              << "w-2-Norm= " << norm2 << std::endl
              << " L2-Norm= (" << L2_vel[0]<<", "<<L2_vel[1]<<", "<<L2_vel[2]<<")" << std::endl
              << " L1-Norm= (" << L1_vel[0]<<", "<<L1_vel[1]<<", "<<L1_vel[2]<<")" << std::endl
              << "max-Norm= " << maxdiff << std::endl;

    norm2= 0; maxdiff= 0; double mindiff= 1000;

    // Compute the pressure-coefficient in direction of 1/sqrt(meas(Omega)), which eliminates
    // the allowed offset of the pressure by setting it to 0.
    double L1_pr= 0, L2_pr= 0, MW_pr= 0, vol= 0;
    DiscPrSolCL pr(lsgpr, &_BndData.Pr, &_MG);
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double sum= 0;
        for(int i=0; i<4; ++i)
            sum+= pr.val(*sit->GetVertex(i)) - LsgPr(sit->GetVertex(i)->GetCoord(), t);
        sum/= 120;
        sum+= 2./15.* (pr.val(*sit, .25, .25, .25) - LsgPr(GetBaryCenter(*sit), t));
        MW_pr+= sum * sit->GetVolume()*6.;
        vol+= sit->GetVolume();
    }
    const double c_pr= MW_pr / vol;
    std::cerr << "\nconstant pressure offset is " << c_pr<<", volume of cube is " << vol<<std::endl;;

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        diff= fabs( c_pr + LsgPr(sit->GetCoord(), t) - pr.val(*sit));
        norm2+= diff*diff;
        if (diff>maxdiff)
            maxdiff= diff;
        if (diff<mindiff)
            mindiff= diff;
    }
    norm2= ::sqrt(norm2 / lsgpr->Data.size() );

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double sum= 0, sum1= 0;
        for(int i=0; i<4; ++i)
        {
            diff= c_pr + LsgPr(sit->GetVertex(i)->GetCoord(), t) - pr.val(*sit->GetVertex(i));
            sum+= diff*diff; sum1+= fabs(diff);
        }
        sum/= 120;   sum1/= 120;
        diff= c_pr + LsgPr(GetBaryCenter(*sit), t) - pr.val(*sit, .25, .25, .25);
        sum+= 2./15.*diff*diff;   sum1+= 2./15.*fabs(diff);
        L2_pr+= sum * sit->GetVolume()*6.;
        L1_pr+= sum1 * sit->GetVolume()*6.;
    }
    L2_pr= sqrt( L2_pr);


    std::cerr << "Druck: Abweichung von der tatsaechlichen Loesung:\n"
              << "w-2-Norm= " << norm2 << std::endl
              << " L2-Norm= " << L2_pr << std::endl
              << " L1-Norm= " << L1_pr << std::endl
              << "Differenz liegt zwischen " << mindiff << " und " << maxdiff << std::endl;
}

} // end of namespace DROPS
