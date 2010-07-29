/// \file params.h
/// \brief test reparametrization of the levelset function
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#include <fstream>
#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "levelset/surfacetension.h"
#include "out/ensightOut.h"

double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha*alpha);
}

double DistFct( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL midpt( 0.5);
    midpt[1]= 0.25;
    return (midpt-p).norm()-0.2;
}

double Phi0( const DROPS::Point3DCL& p)
{
    double dist= DistFct( p);
    return SmoothedSign( dist, 0.1);
}

double Phi1( const DROPS::Point3DCL& p)
{
    return 6*DistFct( p);
}

double Phi2( const DROPS::Point3DCL& p)
{
    const double depth= 0.8-p[1],
                 alpha= 0.1;

    if (depth>0)
        return SmoothedSign( DistFct( p), alpha);
    else
        return SmoothedSign( depth, alpha);
}

double sigma (const DROPS::Point3DCL&, double) { return 0.; }

namespace DROPS
{  // for strategy

double func_abs( double x) { return std::abs(x); }

void SetupReparamSystem( LevelsetP2CL& lset, MatrixCL& M_, MatrixCL& R_, const VectorCL& Psi, VectorCL& b, double diff)
// M, R, b describe the following terms used for reparametrization:
// b_i  = ( S(Phi0),           v_i     )
// M_ij = ( v_j,               v_i     )
// R_ij = ( w(Psi) grad v_j,   v_i     )
//      + (|S(Phi0)| grad v_j, grad v_i) * diff
// where v_i, v_j denote the ansatz functions
// and w(Psi) = sign(Phi0) * grad Psi / |grad Psi| the scaled gradient of Psi
{
    const IdxT num_unks= lset.Phi.RowIdx->NumUnknowns();
    const Uint lvl= lset.Phi.GetLevel();

    SparseMatBuilderCL<double> R(&R_, num_unks, num_unks);
    SparseMatBuilderCL<double> M(&M_, num_unks, num_unks);
    b.resize( 0);
    b.resize( num_unks);
    std::cout << "entering SetupReparamSystem: " << num_unks << " levelset unknowns. ";

    Quad2CL<double>    Sign_Phi;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], w_loc;
    SMatrixCL<3,3>     T;
    P2DiscCL::GetGradientsOnRef( GradRef);

    IdxT         Numb[10];
    SVectorCL<3> grad_Psi[4];
    LevelsetP2CL::const_DiscSolCL    phi= lset.GetSolution();
    double det, absdet;
    const double alpha= 0.1;  // for smoothing of signum fct

    const MultiGridCL& MG = lset.GetMG();
    for (MultiGridCL::const_TriangTetraIteratorCL sit=(MG).GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        GetLocalNumbP2NoBnd( Numb, *sit, *lset.Phi.RowIdx);

        // init Sign_Phi, w_loc
        for (int i=0; i<4; ++i)
        {
            Sign_Phi[i]= SmoothedSign( phi.val( *sit->GetVertex(i)), alpha);
            grad_Psi[i]= Point3DCL();  // init with zero
            for (int l=0; l<10; ++l)
                grad_Psi[i]+= Psi[ Numb[l]] * Grad[l][i];
            w_loc[i]= (Sign_Phi[i]/grad_Psi[i].norm() )*grad_Psi[i];
        }
        // values in barycenter
        Point3DCL gr= 0.25*(grad_Psi[0]+grad_Psi[1]+grad_Psi[2]+grad_Psi[3]);
        Sign_Phi[4]= SmoothedSign( phi.val( *sit, 0.25, 0.25, 0.25), alpha);
        w_loc[4]=    Sign_Phi[4]*gr/gr.norm();

        Quad2CL<> w_Grad[10];
        for (int i=0; i<10; ++i)
            w_Grad[i]= dot(w_loc, Grad[i]);

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            // b_i  = ( S(Phi0),         v_i + SD * w(Psi) grad v_i )
            b[ Numb[i]]+= Sign_Phi.quadP2( i, absdet);
            for(int j=0; j<10; ++j)
            {
                // M_ij = ( v_j, v_i )
                M( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j)*absdet;
                // R_ij = ( w(Psi) grad v_j, v_i + SD * w(Psi) grad v_i )
                R( Numb[i], Numb[j])+= w_Grad[j].quadP2(i, absdet)
                    + diff*Quad2CL<>(dot( Grad[j]*Sign_Phi.apply( func_abs), Grad[i])).quad( absdet);
            }
        }
    }
    M.Build();
    R.Build();
    std::cout << M_.num_nonzeros() << " nonzeros in M, "
              << R_.num_nonzeros() << " nonzeros in R!" << std::endl;
}

void Reparam( LevelsetP2CL& lset, Uint steps, double dt, double theta, double diff)
// Reparametrization of the levelset function Phi
{
    VectorCL Psi= lset.Phi.Data, b;
    MatrixCL L, R, M;

    for (Uint i=0; i<steps; ++i)
    {
        SetupReparamSystem( lset, M, R, Psi, b, diff);
        L.LinComb( 1./dt, M, theta, R);

        b+= (1./dt)*(M*Psi) - (1.-theta) * (R*Psi);
        typedef GMResSolverCL<SSORPcCL> LsetSolverT;
        SSORPcCL ssorpc;
        GMResSolverCL<SSORPcCL> gm( ssorpc, 100, 1000, 1e-7);
        gm.Solve( L, Psi, b);
        std::cout << "Reparam: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
    }

    lset.Phi.Data= Psi;
}


template<class ProblemT>
void Strategy( ProblemT& prob, BndDataCL<>& lsetbnd, double dt, int num_steps, double diff, int bsp)
{
    MultiGridCL& mg= prob.GetMG();
    // Levelset: sigma= 0, theta= 0.5, SD= 0
    SurfaceTensionCL sf( sigma, 0);
    LevelsetP2CL lset( mg, lsetbnd, sf, 0.5);

    IdxDescCL& lidx= lset.idx;
    MLIdxDescCL& vidx= prob.vel_idx;
    VecDescCL& vel=  prob.v;

    prob.CreateNumberingVel( mg.GetLastLevel(), &vidx);
    lset.CreateNumbering(    mg.GetLastLevel(), &lidx);
    vel.SetIdx( &vidx);
    lset.Phi.SetIdx( &lidx);
    switch (bsp)
    {
        case 0:  lset.Init( Phi0); break;
        case 1:  lset.Init( Phi1); break;
        case 2:  lset.Init( Phi2); break;
        default: lset.Init( DistFct);
    }
    lset.SetupSystem( prob.GetVelSolution(), dt);

    // Initialize Ensight6 output
    std::string ensf( "ensight/rep");
    Ensight6OutCL ensight( "rep.case", num_steps + 2);
    ensight.Register( make_Ensight6Geom  ( mg, mg.GetLastLevel(), "Cube",     ensf + ".geo"));
    ensight.Register( make_Ensight6Scalar( lset.GetSolution(),    "Levelset", ensf + ".scl", true));
    ensight.Write( 0.);

    TimerCL time;
    time.Start();
    lset.Reparam();
    time.Stop();
    std::cout << time.GetTime() << " sec for Fast Marching\n";
    ensight.Write( dt/2);

    for (int i=1; i<=num_steps; ++i)
    {
        Reparam( lset, 1, dt, 0.5, diff);
        ensight.Write( i*dt);
    }
}

class DummyStokesCoeffCL {};

} // end of namespace DROPS

int main( int argc, char **argv)
{
  try{
    double dt= 0.01, diff= 1e-4;
    int bsp= 0, num_steps= 100;

    if (argc>1)
        diff= std::atof( argv[1]);
    if (argc>2)
        dt= std::atof( argv[2]);
    if (argc>3)
        num_steps= std::atoi( argv[3]);

    std::cout << num_steps << " steps of lenght dt = " << dt
              << ", diff = " << diff << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    typedef DROPS::StokesP2P1CL<DROPS::DummyStokesCoeffCL> StokesOnBrickCL;
    typedef StokesOnBrickCL                                                  MyStokesCL;

    int num= 0;
    std::cout << "# Unterteilungen: "; std::cin >> num;
    std::cout << "Beispielnr.: "; std::cin >> bsp;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, num, num, num);
    const bool IsNeumann[6]=
        {true, true, true, true, true, true};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel };

    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    MyStokesCL prob(brick, DROPS::DummyStokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));

    Strategy( prob, lsbnd, dt, num_steps, diff, bsp);

    return 0;
  } catch(DROPS::DROPSErrCL err) { err.handle(); }
}
