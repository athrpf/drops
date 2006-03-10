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
  void
  NavierStokesP2P1CL<Coeff>::CheckSolution(
    const VelVecDescCL* lsgvel, const VecDescCL* lsgpr, 
    instat_vector_fun_ptr LsgVel, instat_scalar_fun_ptr LsgPr,
    double t)
{
    double diff, maxdiff=0, norm2= 0;
    Uint lvl=lsgvel->GetLevel(),
         vidx=lsgvel->RowIdx->GetIdx();
    
    VecDescCL rhsN( lsgvel->RowIdx);
    MatDescCL myN( lsgvel->RowIdx, lsgvel->RowIdx);
    SetupNonlinear( &myN, lsgvel, &rhsN, t, t);
    VectorCL res1( A.Data*lsgvel->Data + myN.Data*lsgvel->Data + transp_mul( B.Data, lsgpr->Data)
                 - b.Data - rhsN.Data);
    VectorCL res2( B.Data*lsgvel->Data - c.Data);

    // XXX For instationary equations, this ist wrong, because we ignore the
    // mass-matrix-parts from the time-discretization. Also, InstatNavStokesThetaSchemeCL
    // swaps some right-hand-sides (b) with local vectors, so that we use the wrong one
    // in every second timestep.
    std::cerr << "\nChecken der Loesung des LGS:\n";    
    std::cerr << "|| Ax + Nx + BTy - f || = " << norm( res1) << ", max. " << supnorm( res1) << std::endl;
    std::cerr << "||      Bx       - g || = " << norm( res2) << ", max. " << supnorm( res2) << std::endl<<std::endl;
    
    typename _base::const_DiscVelSolCL vel(lsgvel, &_BndData.Vel, &_MG, t);
    double L1_div= 0, L2_div= 0;
    SMatrixCL<3,3> T;
    double det, absdet;

    
    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double div[5]= {0., 0., 0., 0., 0.};   // Divergenz in den Verts und im BaryCenter
        GetTrafoTr(T,det,*sit);
        absdet= std::fabs(det);
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
        L1_div+= ( (std::fabs(div[0])+std::fabs(div[1])+std::fabs(div[2])+std::fabs(div[3]))/120 + std::fabs(div[4])*2./15. ) * absdet;
        L2_div+= ( (div[0]*div[0]+div[1]*div[1]+div[2]*div[2]+div[3]*div[3])/120 + div[4]*div[4]*2./15. ) * absdet;
    }
    L2_div= std::sqrt(L2_div);
    std::cerr << "|| div x ||_L1 = " << L1_div << std::endl;
    std::cerr << "|| div x ||_L2 = " << L2_div << std::endl << std::endl;

    for (MultiGridCL::TriangVertexIteratorCL sit=_MG.GetTriangVertexBegin(lvl), send=_MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!_BndData.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
           {
               diff= std::fabs( LsgVel(sit->GetCoord(), t)[i] - lsgvel->Data[sit->Unknowns(vidx)+i] );
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
               diff= std::fabs( LsgVel( GetBaryCenter( *sit), t)[i] - lsgvel->Data[sit->Unknowns(vidx)+i] );
               norm2+= diff*diff;
               if (diff>maxdiff)
               {
                   maxdiff= diff;
               }
           }
        }
    }
    norm2= std::sqrt(norm2 / lsgvel->Data.size());

    Point3DCL L1_vel(0.0), L2_vel(0.0);
    for(MultiGridCL::TriangTetraIteratorCL sit= _MG.GetTriangTetraBegin(lvl), send= _MG.GetTriangTetraEnd(lvl);
        sit != send; ++sit)
    {
	Point3DCL sum(0.0), diff, Diff[10];
        const double absdet= sit->GetVolume()*6;
	for(int i=0; i<4; ++i)
	{
	    Diff[i]= diff= LsgVel(sit->GetVertex(i)->GetCoord(), t) - vel.val(*sit->GetVertex(i));
	    sum+= fabs( diff);
	}
	sum/= 120;
	diff= LsgVel(GetBaryCenter(*sit), t) - vel.val(*sit, 0.25, 0.25, 0.25);
	sum+= fabs( diff)*2./15.;
	sum*= absdet;
	L1_vel+= sum;

	for(int i=4; i<10; ++i) // differences on edges
	    Diff[i]= LsgVel( GetBaryCenter(*sit->GetEdge(i-4)), t) - vel.val(*sit->GetEdge(i-4));

	for(int i=0; i<10; ++i)
	{
	    sum= P2DiscCL::Quad(Diff, i)*absdet;
	    sum*= Diff[i];
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
    typename _base::const_DiscPrSolCL pr(lsgpr, &_BndData.Pr, &_MG);
    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
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

    VertexCL* maxvert= 0;
    for (MultiGridCL::TriangVertexIteratorCL sit=_MG.GetTriangVertexBegin(lvl), send=_MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        diff= std::fabs( c_pr + LsgPr(sit->GetCoord(), t) - pr.val(*sit));
        norm2+= diff*diff;
        if (diff>maxdiff)
	{
            maxdiff= diff;
	    maxvert= &*sit;
	}
        if (diff<mindiff)
            mindiff= diff;
    }
    norm2= std::sqrt(norm2 / lsgpr->Data.size() );
    std::cout << "Maximaler Druckfehler: ";
    maxvert->DebugInfo(std::cout);
    std::cout<<std::endl;

    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double sum= 0, sum1= 0;
        for(int i=0; i<4; ++i)
        {
            diff= c_pr + LsgPr(sit->GetVertex(i)->GetCoord(), t) - pr.val(*sit->GetVertex(i));
            sum+= diff*diff; sum1+= std::fabs(diff);
        }
        sum/= 120;   sum1/= 120;
        diff= c_pr + LsgPr(GetBaryCenter(*sit), t) - pr.val(*sit, .25, .25, .25);
        sum+= 2./15.*diff*diff;   sum1+= 2./15.*std::fabs(diff);
        L2_pr+= sum * sit->GetVolume()*6.;
        L1_pr+= sum1 * sit->GetVolume()*6.;
    }
    L2_pr= std::sqrt( L2_pr);


    std::cerr << "Druck: Abweichung von der tatsaechlichen Loesung:\n"
              << "w-2-Norm= " << norm2 << std::endl
              << " L2-Norm= " << L2_pr << std::endl
              << " L1-Norm= " << L1_pr << std::endl
              << "Differenz liegt zwischen " << mindiff << " und " << maxdiff << std::endl;
}


template <class Coeff>
  void
  NavierStokesP2P1CL<Coeff>::SetupNonlinear(MatDescCL* matN, const VelVecDescCL* velvec,
      VelVecDescCL* vecb, double t, double t2) const
// Sets up the approximation of the nonlinear term and the corresponding right hand side.
{
    vecb->Clear();
    
    const IdxT num_unks_vel= matN->RowIdx->NumUnknowns;
    MatrixBuilderCL N( &matN->Data, num_unks_vel, num_unks_vel);

    typename _base::const_DiscVelSolCL u( velvec, &_BndData.Vel, &_MG, t);
    VectorCL& b= vecb->Data;
    const Uint lvl    = matN->GetRowLevel();
    const Uint vidx   = matN->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns
                            
    Quad2CL<Point3DCL> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    Quad2CL<Point3DCL> u_loc;
    SMatrixCL<3,3> T;
    double det, absdet;
    SVectorCL<3> tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl),
                                                 send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        P2DiscCL::GetGradients(Grad, GradRef, T);
        absdet= std::fabs(det);
        u_loc.assign( *sit, u, t);
        
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
            
            for(int i=0; i<10; ++i)    // assemble row Numb[i]
                if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
                {
                    const double N_ij= Quad2CL<>(dot( u_loc, Grad[j])).quadP2( i, absdet);
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        N(Numb[i],          Numb[j])+=          N_ij; 
                        N(Numb[i]+stride,   Numb[j]+stride)+=   N_ij; 
                        N(Numb[i]+2*stride, Numb[j]+2*stride)+= N_ij; 
                    }
                    else // coupling with vert/edge j on right-hand-side
                    {
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j), t2)
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4), t2);
                        b[Numb[i]]-=          N_ij * tmp[0];
                        b[Numb[i]+stride]-=   N_ij * tmp[1];
                        b[Numb[i]+2*stride]-= N_ij * tmp[2];
                    }
                }
        }
    }
    N.Build();
}

} // end of namespace DROPS
