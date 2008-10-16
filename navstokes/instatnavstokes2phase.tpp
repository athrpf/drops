/// \file
/// \brief classes that constitute the 2-phase Navier-Stokes problem

#include "num/discretize.h"

namespace DROPS
{

template<class Coeff>
void InstatNavierStokes2PhaseP2P1CL<Coeff>::SetupNonlinear
    ( MatDescCL* N, const VelVecDescCL* vel, VelVecDescCL* cplN,
      const LevelsetP2CL& lset, double t) const
/// Couplings with dirichlet BCs are accumulated in cplN,
/// so call cplN->Clear() before if only couplings are needed.
{
    const IdxT num_unks_vel= N->RowIdx->NumUnknowns;

    MatrixBuilderCL mN( &N->Data, num_unks_vel, num_unks_vel);
    cplN->Clear();

    const Uint lvl= N->GetRowLevel();
    LocalNumbP2CL n;

    std::cerr << "entering SetupNonlinear: " << num_unks_vel << " vels. ";

    Quad5CL<Point3DCL> Grad[10], GradRef[10], u_loc, u_rho;
    Quad5CL<double> rho, Phi;

    SMatrixCL<3,3> T;

    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const_DiscVelSolCL u( vel, &GetBndData().Vel, &GetMG(), t);

    P2DiscCL::GetGradientsOnRef( GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        Phi.assign( *sit, ls, t);
        u_loc.assign( *sit, u, t);
        n.assign( *sit, *N->RowIdx, _BndData.Vel);

        // rho = rho( Phi)
        rho= Phi;
        rho.apply( _Coeff.rho);

        u_rho= u_loc*rho;

        for(int i= 0; i < 10; ++i)  // assemble row n.num[i]
            if (n.WithUnknowns( i))  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j= 0; j < 10; ++j)
                {   // N(u)_ij = int( phi_i * rho*u * grad phi_j )
                    const double N_ij= Quad5CL<double>( dot( u_rho, Grad[j])).quadP2( i, absdet);
                    if (n.WithUnknowns( j)) // vert/edge j is not on a Dirichlet boundary
                    {
                        mN( n.num[i],   n.num[j]  )+= N_ij;
                        mN( n.num[i]+1, n.num[j]+1)+= N_ij;
                        mN( n.num[i]+2, n.num[j]+2)+= N_ij;
                    }
                    else // put coupling on rhs
                    {
                        tmp= j < 4 ? _BndData.Vel.GetDirBndValue( *sit->GetVertex( j), t)
                                   : _BndData.Vel.GetDirBndValue( *sit->GetEdge( j - 4), t);
                        cplN->Data[n.num[i]    ]-= N_ij*tmp[0];
                        cplN->Data[n.num[i] + 1]-= N_ij*tmp[1];
                        cplN->Data[n.num[i] + 2]-= N_ij*tmp[2];
                    }
                }

            }
    }

    mN.Build();
    std::cerr << N->Data.num_nonzeros() << " nonzeros in N!" << std::endl;
}

template<class Coeff>
void InstatNavierStokes2PhaseP2P1CL<Coeff>::SetupN
    ( MatDescCL* N, const VelVecDescCL* vel, const LevelsetP2CL& lset, double t) const
{
    const IdxT num_unks_vel= N->RowIdx->NumUnknowns;

    MatrixBuilderCL mN( &N->Data, num_unks_vel, num_unks_vel);

    const Uint lvl= N->GetRowLevel();
    LocalNumbP2CL n;

    std::cerr << "entering SetupNonlinear: " << num_unks_vel << " vels. ";

    Quad5CL<Point3DCL> Grad[10], GradRef[10], u_loc, u_rho;
    Quad5CL<double> rho, Phi;

    SMatrixCL<3,3> T;

    double det, absdet;
    Point3DCL tmp;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const_DiscVelSolCL u( vel, &GetBndData().Vel, &GetMG(), t);

    P2DiscCL::GetGradientsOnRef( GradRef);
    LocalP2CL<> p2Phi;
    LocalP2CL<Point3DCL> p2u;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        p2Phi.assign(*sit, ls, t);
        Phi.assign( p2Phi);
        p2u.assign(*sit, u, t);
        u_loc.assign(p2u);

        //u_loc.assign( *sit, u, t);
        n.assign( *sit, *N->RowIdx, _BndData.Vel);

        // rho = rho( Phi)
        rho= Phi;
        rho.apply( _Coeff.rho);

        u_rho= u_loc*rho;

        for(int i= 0; i < 10; ++i)  // assemble row n.num[i]
            if (n.WithUnknowns( i))  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j= 0; j < 10; ++j)
                {   // N(u)_ij = int( phi_i * rho*u * grad phi_j )
                    const double N_ij= Quad5CL<double>( dot( u_rho, Grad[j])).quadP2( i, absdet);
                    if (n.WithUnknowns( j)) // vert/edge j is not on a Dirichlet boundary
                    {
                        mN( n.num[i],   n.num[j]  )+= N_ij;
                        mN( n.num[i]+1, n.num[j]+1)+= N_ij;
                        mN( n.num[i]+2, n.num[j]+2)+= N_ij;
                    }
                }

            }
    }

    mN.Build();
    std::cerr << N->Data.num_nonzeros() << " nonzeros in N!" << std::endl;
}

template <class Coeff>
void InstatNavierStokes2PhaseP2P1CL<Coeff>::SetupNMatricesMG (MGDataCL* matMG, const VelVecDescCL& vel, double alpha) const
// Set up the MG-hierarchy
{
    for(MGDataCL::iterator it= matMG->begin(); it!=matMG->end(); ++it) {
        MGLevelDataCL& tmp= *it;
        MatDescCL N;
        if(&tmp != &matMG->back()) {
            N.SetIdx( &tmp.Idx, &tmp.Idx);
            tmp.AN.SetIdx( &tmp.Idx, &tmp.Idx);
            std::cerr << "Create Nonlinear for "
                    << (&tmp.Idx)->NumUnknowns << " unknowns." << std::endl;
            SetupN( &N, &vel, *ls_, t);
            tmp.AN.Data.LinComb( 1., tmp.A.Data, alpha, N.Data);
        }
        tmp.ABlock = &tmp.AN.Data;
    }
}


} // end of namespace DROPS
