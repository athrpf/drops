//**************************************************************************
// File:    instatnavstokes2phase.tpp                                      *
// Content: classes that constitute the 2-phase navier-stokes-problem      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "num/discretize.h"

namespace DROPS
{

template<class Coeff>
void InstatNavierStokes2PhaseP2P1CL<Coeff>::SetupNonlinear 
    ( MatDescCL* N, const VelVecDescCL* vel, VelVecDescCL* cplN, 
      const LevelsetP2CL& lset, double t) const
{
    const IdxT num_unks_vel= N->RowIdx->NumUnknowns;

    MatrixBuilderCL mN( &N->Data, num_unks_vel, num_unks_vel);
    cplN->Clear();
    
    const Uint lvl         = N->GetRowLevel(),
               vidx        = N->RowIdx->GetIdx();

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupNonlinear: " << num_unks_vel << " vels. ";

    Quad2CL<Point3DCL> Grad[10], GradRef[10], u_loc, u_rho;
    Quad2CL<double> rho, Phi;
        
    SMatrixCL<3,3> T;
    
    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::DiscSolCL ls= lset.GetSolution();
    DiscVelSolCL u( vel, &GetBndData().Vel, &GetMG(), t);

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
            Phi.val[i]= ls.val( *sit->GetVertex(i));
            u_loc.val[i]= u.val( *sit->GetVertex(i));
        }
        for (int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }
        u_loc.val[4]= u.val( *sit, 0.25, 0.25, 0.25);
        Phi.val[4]= ls.val( *sit, 0.25, 0.25, 0.25);

        // rho = rho( Phi)
        rho= Phi;    rho.apply( _Coeff.rho);

        u_rho= u_loc*rho;
        
        for(int i=0; i<10; ++i)  // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {   // N(u)_ij = int( phi_i * rho*u * grad phi_j )
                    const double N_ij= Quad2CL<double>(dot( u_rho, Grad[j])).quadP2( i, absdet);
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        mN( Numb[i],   Numb[j]  )+= N_ij;
                        mN( Numb[i]+1, Numb[j]+1)+= N_ij;
                        mN( Numb[i]+2, Numb[j]+2)+= N_ij;
                    }
                    else // put coupling on rhs
                    {
                        tmp= j<4 ? _BndData.Vel.GetDirBndValue(*sit->GetVertex(j), t)
                                 : _BndData.Vel.GetDirBndValue(*sit->GetEdge(j-4), t);
                        cplN->Data[Numb[i]  ]-= N_ij*tmp[0];
                        cplN->Data[Numb[i]+1]-= N_ij*tmp[1];
                        cplN->Data[Numb[i]+2]-= N_ij*tmp[2];
                    }
                }

            }
    }

    mN.Build();
    std::cerr << N->Data.num_nonzeros() << " nonzeros in N!" << std::endl;
}

} // end of namespace DROPS
