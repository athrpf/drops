/// \file
/// \brief classes that constitute the Navier-Stokes problem
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

/// History: begin - March, 16 2001

namespace DROPS
{

template <class Coeff>
  void
  NavierStokesP2P1CL<Coeff>::CheckSolution(
    const VelVecDescCL* lsgvel, MLIdxDescCL* idx, const VecDescCL* lsgpr,
    instat_vector_fun_ptr LsgVel, instat_scalar_fun_ptr LsgPr)
{
    const double t =  lsgvel->t;
    double diff, maxdiff=0, norm2= 0;
    Uint lvl=lsgvel->GetLevel(),
         vidx=lsgvel->RowIdx->GetIdx();
#ifdef _PAR
    ExchangeCL ExVel= GetEx(base_::velocity);
    ExchangeCL ExPr = GetEx(base_::pressure);
#endif

    VecDescCL rhsN( lsgvel->RowIdx);
    MLMatDescCL myN( idx, idx);
    SetupNonlinear( &myN, lsgvel, &rhsN, t);
    VectorCL res1( A.Data*lsgvel->Data + myN.Data*lsgvel->Data + transp_mul( B.Data, lsgpr->Data)
                 - b.Data - rhsN.Data);
    VectorCL res2( B.Data*lsgvel->Data - c.Data);

    // XXX For instationary equations, this ist wrong, because we ignore the
    // mass-matrix-parts from the time-discretization. Also, InstatNavStokesThetaSchemeCL
    // swaps some right-hand-sides (b) with local vectors, so that we use the wrong one
    // in every second timestep.
    double norm_res1, norm_res2,
           sup_res1 = supnorm(res1),
           sup_res2 = supnorm(res2);
#ifndef _PAR
    norm_res1 = norm(res1);
    norm_res2 = norm(res2);
#else
    norm_res1 =  ExVel.Norm(res1, false);
    norm_res2 =  ExPr.Norm(res2, false);
    sup_res1  =  GlobalMax(sup_res1);
    sup_res2  =  GlobalMax(sup_res2);
#endif

    IF_MASTER
      std::cout << "\nChecken der Loesung des LGS:\n"
                << "|| Ax + Nx + BTy - f || = " << norm_res1 << ", max. " << sup_res1 << std::endl
                << "||      Bx       - g || = " << norm_res2 << ", max. " << sup_res2 << std::endl<<std::endl;

    typename base_::const_DiscVelSolCL vel(lsgvel, &BndData_.Vel, &MG_);
    double L1_div= 0, L2_div= 0;
    SMatrixCL<3,3> T;
    double det, absdet;


    for (MultiGridCL::TriangTetraIteratorCL sit=MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
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
#ifdef _PAR
    L1_div = GlobalSum(L1_div);
    L2_div = GlobalSum(L2_div);
#endif

    L2_div= std::sqrt(L2_div);

    IF_MASTER
      std::cout << "|| div x ||_L1 = " << L1_div << '\n'
                << "|| div x ||_L2 = " << L2_div << '\n' << std::endl;

    for (MultiGridCL::TriangVertexIteratorCL sit=MG_.GetTriangVertexBegin(lvl), send=MG_.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!BndData_.Vel.IsOnDirBnd(*sit))
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
    for (MultiGridCL::TriangEdgeIteratorCL sit=MG_.GetTriangEdgeBegin(lvl), send=MG_.GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!BndData_.Vel.IsOnDirBnd(*sit))
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
    IdxT vel_size = lsgvel->Data.size();
#ifdef _PAR
    vel_size = GlobalSum(vel_size);
    norm2    = GlobalSum(norm2);
#endif

    norm2= std::sqrt(norm2 / vel_size);

    Point3DCL L1_vel(0.0), L2_vel(0.0);
    for(MultiGridCL::TriangTetraIteratorCL sit= MG_.GetTriangTetraBegin(lvl), send= MG_.GetTriangTetraEnd(lvl);
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
#ifdef _PAR
    for (size_t i=0; i<3; ++i)
        L2_vel[i]= GlobalSum(L2_vel[i]);
    maxdiff = GlobalMax(maxdiff);
#endif
    L2_vel= sqrt(L2_vel);

    IF_MASTER
      std::cout << "Geschwindigkeit: Abweichung von der tatsaechlichen Loesung:\n"
                << "w-2-Norm= " << norm2 << std::endl
                << " L2-Norm= (" << L2_vel[0]<<", "<<L2_vel[1]<<", "<<L2_vel[2]<<")" << std::endl
                << " L1-Norm= (" << L1_vel[0]<<", "<<L1_vel[1]<<", "<<L1_vel[2]<<")" << std::endl
                << "max-Norm= " << maxdiff << std::endl;

    norm2= 0; maxdiff= 0; double mindiff= 1000;

    // Compute the pressure-coefficient in direction of 1/sqrt(meas(Omega)), which eliminates
    // the allowed offset of the pressure by setting it to 0.
    double L1_pr= 0, L2_pr= 0, MW_pr= 0, vol= 0;
    typename base_::const_DiscPrSolCL pr(lsgpr, &BndData_.Pr, &MG_);
    for (MultiGridCL::TriangTetraIteratorCL sit=MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
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
#ifdef _PAR
    vol  = GlobalSum(vol);
    MW_pr= GlobalSum(MW_pr);
#endif
    const double c_pr= MW_pr / vol;
    IF_MASTER
      std::cout << "\nconstant pressure offset is " << c_pr<<", volume of cube is " << vol<<std::endl;;

    VertexCL* maxvert= 0;
    for (MultiGridCL::TriangVertexIteratorCL sit=MG_.GetTriangVertexBegin(lvl), send=MG_.GetTriangVertexEnd(lvl);
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
    IdxT pr_size = lsgpr->Data.size();
#ifdef _PAR
    pr_size= GlobalSum(pr_size);
    norm2  = GlobalSum(norm2);
    mindiff= GlobalMin(mindiff);
    maxdiff= GlobalMax(maxdiff);
#endif
    norm2= std::sqrt(norm2 / pr_size );

#ifndef _PAR
    std::cout << "Maximaler Druckfehler: ";
    maxvert->DebugInfo(std::cout);
    std::cout<<std::endl;
#endif

    for (MultiGridCL::TriangTetraIteratorCL sit=MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
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
#ifdef _PAR
    L1_pr = GlobalSum(L1_pr);
    L2_pr = GlobalSum(L2_pr);
#endif
    L2_pr= std::sqrt( L2_pr);

    IF_MASTER
      std::cout << "Druck: Abweichung von der tatsaechlichen Loesung:\n"
                << "w-2-Norm= " << norm2 << std::endl
                << " L2-Norm= " << L2_pr << std::endl
                << " L1-Norm= " << L1_pr << std::endl
                << "Differenz liegt zwischen " << mindiff << " und " << maxdiff << std::endl;
}


template <class Coeff>
  void
  NavierStokesP2P1CL<Coeff>::SetupNonlinear_P2( MatrixCL& matN, const VelVecDescCL* velvec,
                                                VelVecDescCL* vecb, IdxDescCL& RowIdx, double t2) const
// Sets up the approximation of the nonlinear term and the corresponding right hand side.
{
    if (vecb != 0) vecb->Clear( t2);

    const IdxT num_unks_vel= RowIdx.NumUnknowns();
    MatrixBuilderCL N( &matN, num_unks_vel, num_unks_vel);

    typename base_::const_DiscVelSolCL u( velvec, &BndData_.Vel, &MG_);
    VectorCL& b= vecb->Data;
    const Uint lvl    = RowIdx.TriangLevel();
    LocalNumbP2CL n;

    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns

    Quad5CL<Point3DCL> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    Quad5CL<Point3DCL> u_loc;
    SMatrixCL<3,3> T;
    double det, absdet;
    SVectorCL<3> tmp;

    P2DiscCL::GetGradientsOnRef( GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>( MG_).GetTriangTetraBegin(lvl),
                                                 send=const_cast<const MultiGridCL&>( MG_).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        P2DiscCL::GetGradients(Grad, GradRef, T);
        absdet= std::fabs(det);
        u_loc.assign( *sit, u);

        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        n.assign( *sit, RowIdx, BndData_.Vel);

        for (int i= 0; i < 10; ++i) // assemble row n.num[i]
            if (n.WithUnknowns( i)) // vert/edge i is not on a Dirichlet boundary
                for (int j= 0; j < 10; ++j)
                { // N(u)_ij = int( phi_i *( u1 phi_j_x + u2 phi_j_y + u3 phi_j_z ) )
                    const double N_ij= Quad5CL<>(dot( u_loc, Grad[j])).quadP2( i, absdet);
                    if (n.WithUnknowns( j)) // vert/edge j is not on a Dirichlet boundary
                    {
                        N(n.num[i],          n.num[j])+=          N_ij;
                        N(n.num[i]+stride,   n.num[j]+stride)+=   N_ij;
                        N(n.num[i]+2*stride, n.num[j]+2*stride)+= N_ij;
                    }
                    else // coupling with vert/edge j on right-hand-side
                        if (vecb != 0)
                        {
                            tmp= j<4 ? BndData_.Vel.GetDirBndValue(*sit->GetVertex(j), t2)
                                    : BndData_.Vel.GetDirBndValue(*sit->GetEdge(j-4), t2);
                            b[n.num[i]]-=          N_ij * tmp[0];
                            b[n.num[i]+stride]-=   N_ij * tmp[1];
                            b[n.num[i]+2*stride]-= N_ij * tmp[2];
                        }
                }
    }
    N.Build();
}

template <class Coeff>
  void
  NavierStokesP2P1CL<Coeff>::SetupNonlinear(MLMatDescCL* matN, const VelVecDescCL* velvec,
      VelVecDescCL* vecb, double t2) const
{
    MLIdxDescCL::iterator itRow= matN->RowIdx->begin();
    MLMatrixCL::iterator itN= matN->Data.begin();
    for (size_t lvl=0; lvl < matN->Data.size(); ++lvl, ++itN, ++itRow)
    {
        if (lvl != matN->Data.size()-1)
            SetupNonlinear_P2( *itN, 0, 0, *itRow, t2);
        else
            SetupNonlinear_P2( *itN, velvec, vecb, *itRow, t2);
    }
}

} // end of namespace DROPS
