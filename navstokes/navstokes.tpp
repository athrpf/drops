//**************************************************************************
// File:    navstokes.tpp                                                  *
// Content: classes that constitute the navier-stokes-problem              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

namespace DROPS
{

template <class Coeff>
void NavierStokesP2P1CL<Coeff>::GetDiscError(vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr)
{
    Uint lvl= A.RowIdx->TriangLevel,
        vidx= A.RowIdx->GetIdx(),
        pidx= B.RowIdx->GetIdx();
    VecDescCL veldesc( A.RowIdx);
    VectorCL& lsgvel= veldesc.Data;
    VectorCL  lsgpr( B.RowIdx->NumUnknowns);

    for (MultiGridCL::TriangVertexIteratorCL sit=_MG.GetTriangVertexBegin(lvl), send=_MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
               lsgvel[sit->Unknowns(vidx)+i]= LsgVel(sit->GetCoord())[i];
        }
    }
    
    for (MultiGridCL::TriangEdgeIteratorCL sit=_MG.GetTriangEdgeBegin(lvl), send=_MG.GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2.)[i];
        }
    }
    for (MultiGridCL::TriangVertexIteratorCL sit=_MG.GetTriangVertexBegin(lvl), send=_MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
        lsgpr[sit->Unknowns(pidx)]= LsgPr(sit->GetCoord());
    VecDescCL rhsN( A.ColIdx);
    N.SetIdx( A.RowIdx, A.ColIdx);
    N.Data.clear();
    SetupNonlinear( &N, &veldesc, &rhsN);
    std::cerr << "discretization error to check the system (x,y = continuos solution): "<<std::endl;
    VectorCL res= A.Data*lsgvel + N.Data*lsgvel + transp_mul(B.Data,lsgpr) - b.Data - rhsN.Data; 
    std::cerr <<"|| Ax + Nx + BTy - f || = "<< res.norm()<<", max "<<res.supnorm()<<std::endl;
    VectorCL resB= B.Data*lsgvel - c.Data; 
    std::cerr <<"|| Bx - g || = "<< resB.norm()<<", max "<<resB.supnorm()<<std::endl;
//    std::cerr << res << std::endl;
/*    VectorCL resC= N.Data*lsgvel - rhsN.Data;
    std::cerr <<"|| Nx - rhsN || = "<< resC.norm()<<", max "<<resC.supnorm()<<std::endl;
*/}


template <class Coeff>
void NavierStokesP2P1CL<Coeff>::CheckSolution(const VelVecDescCL* lsgvel, const VecDescCL* lsgpr, 
                                 vector_fun_ptr LsgVel,      scalar_fun_ptr LsgPr)
{
    double diff, maxdiff=0, norm2= 0;
    Uint lvl=lsgvel->RowIdx->TriangLevel,
         vidx=lsgvel->RowIdx->GetIdx();
    
    VecDescCL rhsN( lsgvel->RowIdx);
    N.Data.clear();
    SetupNonlinear( &N, lsgvel, &rhsN);
    VectorCL res1= A.Data*lsgvel->Data + N.Data*lsgvel->Data + transp_mul( B.Data, lsgpr->Data)
                 - b.Data - rhsN.Data;
    VectorCL res2= B.Data*lsgvel->Data - c.Data;

    std::cerr << "\nChecken der Loesung...\n";    
    std::cerr << "|| Ax + Nx + BTy - f || = " << res1.norm() << ", max. " << res1.supnorm() << std::endl;
    std::cerr << "||      Bx       - g || = " << res2.norm() << ", max. " << res2.supnorm() << std::endl<<std::endl;
    
    P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL>  vel(lsgvel, &_BndData.Vel, &_MG);
    double L1_div= 0, L2_div= 0;
    SMatrixCL<3,3> T;
    double det, absdet;

    
    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
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

    for (MultiGridCL::TriangVertexIteratorCL sit=_MG.GetTriangVertexBegin(lvl), send=_MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
           {
               diff= fabs( LsgVel(sit->GetCoord())[i] - lsgvel->Data[sit->Unknowns(vidx)+i] );
               norm2+= diff*diff;
               if (diff>maxdiff)
               {
                   maxdiff= diff;
               }
           }
        }
    }
    
    for (MultiGridCL::TriangEdgeIteratorCL sit=_MG.GetTriangEdgeBegin(lvl), send=_MG.GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
           {
               diff= fabs( LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2.)[i] - lsgvel->Data[sit->Unknowns(vidx)+i] );
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
    for(MultiGridCL::TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(lvl), send= _MG.GetTriangTetraEnd(lvl);
        sit != send; ++sit)
    {
	Point3DCL sum(0.0), diff, Diff[5];
	for(int i=0; i<4; ++i)
	{
	    Diff[i]= diff= LsgVel(sit->GetVertex(i)->GetCoord()) - vel.val(*sit->GetVertex(i));
	    diff[0]= fabs(diff[0]); diff[1]= fabs(diff[1]); diff[2]= fabs(diff[2]);
	    sum+= diff;
	}
	sum/= 120;
	Diff[4]= diff= LsgVel(GetBaryCenter(*sit)) - vel.val(*sit, 0.25, 0.25, 0.25);
	diff[0]= fabs(diff[0]); diff[1]= fabs(diff[1]); diff[2]= fabs(diff[2]);
	sum+= diff*2./15.;
	sum*= sit->GetVolume()*6;
	L1_vel+= sum;

	for(int i=0; i<10; ++i)
	{
	    sum= Quad(Diff, i)*sit->GetVolume()*6;
	    diff= i<4 ? Diff[i] : LsgVel( (sit->GetEdge(i-4)->GetVertex(0)->GetCoord() + sit->GetEdge(i-4)->GetVertex(1)->GetCoord() )/2) - vel.val(*sit->GetEdge(i-4));
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
    P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr(lsgpr, &_BndData.Pr, &_MG);
    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double sum= 0;
        for(int i=0; i<4; ++i)
            sum+= pr.val(*sit->GetVertex(i)) - LsgPr(sit->GetVertex(i)->GetCoord());
        sum/= 120;
        sum+= 2./15.* (pr.val(*sit, .25, .25, .25) - LsgPr(GetBaryCenter(*sit)));
        MW_pr+= sum * sit->GetVolume()*6.;
        vol+= sit->GetVolume();
    }
    const double c_pr= MW_pr / vol;
    std::cerr << "\nconstant pressure offset is " << c_pr<<", volume of cube is " << vol<<std::endl;;

    for (MultiGridCL::TriangVertexIteratorCL sit=_MG.GetTriangVertexBegin(lvl), send=_MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        diff= fabs( c_pr + LsgPr(sit->GetCoord()) - pr.val(*sit));
        norm2+= diff*diff;
        if (diff>maxdiff)
            maxdiff= diff;
        if (diff<mindiff)
            mindiff= diff;
    }
    norm2= ::sqrt(norm2 / lsgpr->Data.size() );

    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double sum= 0, sum1= 0;
        for(int i=0; i<4; ++i)
        {
            diff= c_pr + LsgPr(sit->GetVertex(i)->GetCoord()) - pr.val(*sit->GetVertex(i));
            sum+= diff*diff; sum1+= fabs(diff);
        }
        sum/= 120;   sum1/= 120;
        diff= c_pr + LsgPr(GetBaryCenter(*sit)) - pr.val(*sit, .25, .25, .25);
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


/*****************************************************************************************************
* formulas for   n u m e r i c   i n t e g r a t i o n   on the reference tetrahedron
*****************************************************************************************************/

inline double Quad(double f[5], int i)
{
    if (i<4) // hat function on vert
    {
        std::swap( f[0], f[i]);
//        return f[0]*2./81. + (f[1] + f[2] + f[3])*71./3240. - f[4]*8./81.;
        return f[0]/504. - (f[1] + f[2] + f[3])/1260. - f[4]/126.;
    }
    else  // hat function on edge
    {
//        const double ve= 19./810.,  // coeff for verts of edge
//                     vn= 29./1620.,  // coeff for other verts
//                     vs= -4./81.;   // coeff for barycenter
        const double ve= 4./945.,  // coeff for verts of edge
                     vn= -1./756.,  // coeff for other verts
                     vs= 26./945.;   // coeff for barycenter
        double a[4];
        a[VertOfEdge(i-4,0)]= a[VertOfEdge(i-4,1)]= ve;
        a[VertOfEdge(OppEdge(i-4),0)]= a[VertOfEdge(OppEdge(i-4),1)]= vn;

        double sum= vs * f[4];
        for(int k=0; k<4; ++k)
            sum+= a[k] * f[k];

        return sum;
    }
}

template <class Coeff>
void NavierStokesP2P1CL<Coeff>::SetupNonlinear( MatDescCL* matN, const VelVecDescCL* velvec, VelVecDescCL* vecb) const
// Sets up the stiffness matrices and right hand sides
{
    vecb->Clear();
    
    const IdxT num_unks_vel= matN->RowIdx->NumUnknowns;
    MatrixBuilderCL N( &matN->Data, num_unks_vel, num_unks_vel);

    typename _base::DiscVelSolCL u( velvec, &_BndData.Vel, &_MG);
    VectorCL& b= vecb->Data;
    const Uint lvl    = matN->RowIdx->TriangLevel;
    const Uint vidx   = matN->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns
                            
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup;
    double det, absdet;
    SVectorCL<3> tmp;
    double func[5];

    GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl),
                                                 send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
        for(int i=0; i<6; ++i)
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);

        for(int j=0; j<10; ++j)
        {
            // N(u)_ij = int( phi_i *( u1 phi_j_x + u2 phi_j_y + u3 phi_j_z ) )
            //                       \______________ func __________________/
            
            for( Uint pt=0; pt<5; ++pt) // eval func at 5 points:
            {
                tmp= pt<4 ? u.val( *sit->GetVertex(pt))      // value of u in vert pt
                          : u.val( *sit, 0.25, 0.25, 0.25);  // value of u in barycenter
                func[pt]= Grad[j](0,pt) * tmp[0]
                        + Grad[j](1,pt) * tmp[1]
                        + Grad[j](2,pt) * tmp[2];
            }
            for(int i=0; i<10; ++i)    // assemble row Numb[i]
                if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
                {
                    coup= Quad( func, i) * absdet;
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        N(Numb[i],          Numb[j])+=          coup; 
                        N(Numb[i]+stride,   Numb[j]+stride)+=   coup; 
                        N(Numb[i]+2*stride, Numb[j]+2*stride)+= coup; 
                    }
                    else // coupling with vert/edge j on right-hand-side
                    {
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j))
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4));
                        b[Numb[i]]-=          coup * tmp[0];
                        b[Numb[i]+stride]-=   coup * tmp[1];
                        b[Numb[i]+2*stride]-= coup * tmp[2];
                    }
                }
        }
    }
    N.Build();
}


/*
template <class Coeff>
void NavierStokesP2P1CL<Coeff>::SetupNonlinearRhs( const VelVecDescCL* velvec, VelVecDescCL* vecb) const
// Sets up the right hand sides corresponding to the sgtiffness matrix for N
{
    vecb->Clear();
    
    typename _base::DiscVelSolCL u( velvec, &_BndData.Vel, &_MG);
    VectorCL& b= vecb->Data;
    const Uint lvl    = velvec->RowIdx->TriangLevel;
    const Uint vidx   = velvec->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns
                            
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup;
    double det, absdet;
    SVectorCL<3> tmp;
    double func[5];

    GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl),
                                                 send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
            if(!(IsOnDirBnd[i]= _BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
        for(int i=0; i<6; ++i)
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);

        for(int j=0; j<10; ++j)
        {
            // N(u)_ij = int( phi_i *( u1 phi_j_x + u2 phi_j_y + u3 phi_j_z ) )
            //                       \______________ func __________________/
            
            for( Uint pt=0; pt<5; ++pt) // eval func at 5 points:
            {
                tmp= pt<4 ? u.val( *sit->GetVertex(pt))      // value of u in vert pt
                          : u.val( *sit, 0.25, 0.25, 0.25);  // value of u in barycenter
                func[pt]= Grad[j](0,pt) * tmp[0]
                        + Grad[j](1,pt) * tmp[1]
                        + Grad[j](2,pt) * tmp[2];
            }
            for(int i=0; i<10; ++i)    // assemble row Numb[i]
                if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
                {
                    coup= Quad( func, i) * absdet;
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                    }
                    else // coupling with vert/edge j on right-hand-side
                    {
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j))
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4));
                        b[Numb[i]]-=          coup * tmp[0];
                        b[Numb[i]+stride]-=   coup * tmp[1];
                        b[Numb[i]+2*stride]-= coup * tmp[2];
                    }
                }
        }
    }
}
*/

} // end of namespace DROPS
