/// \file instatnavstokes2phase.cpp
/// \brief classes that constitute the 2-phase Navier-Stokes problem
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt, Christoph Lehrenfeld; SC RWTH Aachen: Oliver Fortmeier

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

#include "navstokes/instatnavstokes2phase.h"

namespace DROPS
{

/// \brief Local element matrix (this object is more or less only introduced for consistency reason)
struct LocalNonlConvDataCL
{
    double         C [10][10];
};

/// \brief Setup of the local (nonlinear) convection term on a tetra in a single phase.
class LocalNonlConvSystemOnePhase_P2CL
{
  private:
    double rho_;

    Quad5CL<Point3DCL> Grad[10], GradRef[10];
    Quad5CL<Point3DCL> vel_;
    const SVectorCL<Quad2DataCL::NumNodesC> Ones;

  public:
    LocalNonlConvSystemOnePhase_P2CL ()
        : Ones( 1.)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void   velocity  (const LocalP2CL<Point3DCL> & velp2)          { vel_.assign(velp2); }
    const Quad5CL<Point3DCL> & velocity  () const { return vel_;        }

    void   rho (double new_rho)                   { rho_= new_rho;      }
    double rho () const                           { return rho_;        }

    void setup (const SMatrixCL<3,3>& T, double absdet, LocalNonlConvDataCL& loc);
};

void LocalNonlConvSystemOnePhase_P2CL::setup (const SMatrixCL<3,3>& T, double absdet, LocalNonlConvDataCL& loc)
{
    P2DiscCL::GetGradients( Grad, GradRef, T);
    for (Uint i= 0; i < 10; ++i) {
        for (Uint j= 0; j < 10; ++j) {
            loc.C[i][j]= rho() * Quad5CL<>( dot(velocity(),Grad[j])).quadP2(i,absdet);
        }
    }
}


/// \brief Setup of the local (nonl.) convection term on a tetra with two phases
///  using smoothed versions of disc. density 
class LocalNonlConvSystemSmoothedJumps_P2CL
{
  private:
    const SmoothedJumpCL & rho_;
    
    Quad5CL<Point3DCL> Grad[10], GradRef[10];
    Quad5CL<Point3DCL> vel_;
    Quad5CL<double> lset_;
    Quad5CL<double> qrho_;
    const SVectorCL<Quad2DataCL::NumNodesC> Ones;

  public:
    LocalNonlConvSystemSmoothedJumps_P2CL (const SmoothedJumpCL & rho_arg)
        : rho_( rho_arg), Ones( 1.)
    { 
        P2DiscCL::GetGradientsOnRef( GradRef); 

    }

    void   velocity  (const LocalP2CL<Point3DCL> & velp2) { vel_.assign(velp2); }
    const Quad5CL<Point3DCL> & velocity  () const   { return vel_;        }

    void   levelset  (const LocalP2CL<double> & lsetp2) { 
        lset_.assign(lsetp2); 
        qrho_ = lset_;
        qrho_.apply( rho_);
    }
    const Quad5CL<double> & levelset () const     { return lset_;  }
    const Quad5CL<double> & smoothed_rho () const { return qrho_; }

    void setup (const SMatrixCL<3,3>& T, double absdet, LocalNonlConvDataCL& loc);
};

void LocalNonlConvSystemSmoothedJumps_P2CL::setup (const SMatrixCL<3,3>& T, double absdet, LocalNonlConvDataCL& loc)
{
    P2DiscCL::GetGradients( Grad, GradRef, T);
    for (Uint i= 0; i < 10; ++i) {
        for (Uint j= 0; j < 10; ++j) {
            loc.C[i][j]= Quad5CL<>( smoothed_rho() * dot(velocity(),Grad[j])) .quadP2(i,absdet);
        }
    }
}



/// \brief Setup of the local (nonl.) convection term on a tetra with two phases
///  using integration rules which are exact w.r.t. to the _approximated_ interface.
class LocalNonlConvSystemTwoPhase_P2CL
{
  private:
    const PrincipalLatticeCL& lat;

    const double rho_p, rho_n;

    LocalP1CL<Point3DCL> GradRef[10], Grad[10];
    LocalP2CL<> p2;

    std::valarray<double> ls_loc; //level set values in partition int. points
    TetraPartitionCL partition;   //the partitioning
    QuadDomainCL q5dom;

    std::valarray<double> qshape[10];
    GridFunctionCL<Point3DCL> qdshape[10];

    double intpos, intneg;
    SMatrixCL<3,3> cAkp, cAkn;


  public:
    LocalNonlConvSystemTwoPhase_P2CL (double rhop, double rhon)
        : lat( PrincipalLatticeCL::instance( 2)), rho_p( rhop), rho_n( rhon), ls_loc( 10)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    double rho (int sign) const                   { return sign > 0 ? rho_p : rho_n; }

    void setup (const SMatrixCL<3,3>& T, double absdet, const LocalP2CL<Point3DCL> & velp2, const LocalP2CL<>& ls, LocalNonlConvDataCL& loc);
};

void LocalNonlConvSystemTwoPhase_P2CL::setup (const SMatrixCL<3,3>& T, double absdet, const LocalP2CL<Point3DCL> & velp2, const LocalP2CL<>& ls, LocalNonlConvDataCL& loc)
{
    P2DiscCL::GetGradients( Grad, GradRef, T);

    evaluate_on_vertexes( ls, lat, Addr( ls_loc));
    partition.make_partition<SortedVertexPolicyCL, MergeCutPolicyCL>( lat, ls_loc);
    make_CompositeQuad5Domain( q5dom, partition);
    GridFunctionCL<Point3DCL> velocity;
    resize_and_evaluate_on_vertexes( velp2, q5dom, velocity);
    for (int i= 0; i < 10; ++i) {
        p2[i]= 1.; p2[i==0 ? 9 : i - 1]= 0.;
        resize_and_evaluate_on_vertexes( p2,      q5dom,  qshape[i]); // shape 
        resize_and_evaluate_on_vertexes( Grad[i], q5dom, qdshape[i]); // gradient shape
    }
    for (int i= 0; i < 10; ++i) {
        for (int j= 0; j < 10; ++j) {
            quad( qshape[i]*dot(velocity,qdshape[j]), absdet, q5dom, intneg, intpos);
            loc.C[i][j]= rho_p*intpos + rho_n*intneg;
        }
    }
}

/// \brief Accumulator to set up the matrices N and, if requested the right-hand side cplN for two-phase flow.
class NonlConvSystemAccumulator_P2CL : public TetraAccumulatorCL
{
  private:
    const bool smoothed; /// use smoothed coefficients jumps instead of really jumps (=true coincides with the old version!)
    const TwoPhaseFlowCoeffCL& Coeff;
    const StokesBndDataCL& BndData;
    const MultiGridCL& MG;
    const VelVecDescCL & vel;
    const LevelsetP2CL& lset;
   
    double t;

    IdxDescCL& RowIdx;
    MatrixCL& N;
    VecDescCL* cplN;

    SparseMatBuilderCL<double, SDiagMatrixCL<3> >* mN_;

    LocalNonlConvSystemOnePhase_P2CL local_onephase; ///< used on tetras in a single phase
    LocalNonlConvSystemTwoPhase_P2CL local_twophase; ///< used on intersected tetras
    LocalNonlConvSystemSmoothedJumps_P2CL local_smoothed_twophase; ///< used on intersected tetras
    LocalNonlConvDataCL loc; ///< Contains the memory, in which the local operators are set up; former coupM, coupA, coupAk, rho_phi.

    LocalNumbP2CL n; ///< global numbering of the P2-unknowns

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;
    LocalP2CL<Point3DCL> vel_loc;

    Point3DCL dirichlet_val[10]; ///< Used to transfer boundary-values from local_setup() update_global_system().

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and, if required, the Dirichlet-values needed to eliminate the boundary-dof from the global system.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    NonlConvSystemAccumulator_P2CL (const TwoPhaseFlowCoeffCL& Coeff, const MultiGridCL& MG_, const StokesBndDataCL& BndData_, 
                                    const VelVecDescCL& vel_, const LevelsetP2CL& ls, IdxDescCL& RowIdx_, 
                                    MatrixCL& N_, VecDescCL* cplN_, double t, bool smoothed=false);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    TetraAccumulatorCL* clone (int /*tid*/) { return new NonlConvSystemAccumulator_P2CL ( *this); };
};

NonlConvSystemAccumulator_P2CL::NonlConvSystemAccumulator_P2CL (const TwoPhaseFlowCoeffCL& Coeff_, const MultiGridCL& MG_, const StokesBndDataCL& BndData_,
                                                                const VelVecDescCL& vel_, const LevelsetP2CL& lset_arg, IdxDescCL& RowIdx_, 
                                                                MatrixCL& N_, VecDescCL* cplN_, double t_, bool smoothed_)
    : smoothed(smoothed_), Coeff( Coeff_), BndData( BndData_), MG(MG_),
      vel(vel_), lset( lset_arg), t( t_),
      RowIdx( RowIdx_), N( N_), cplN( cplN_), 
      local_twophase( Coeff.rho( 1.0), Coeff.rho( -1.0)),
      local_smoothed_twophase( Coeff.rho)
{}

void NonlConvSystemAccumulator_P2CL::begin_accumulation ()
{
    std::cout << "entering NonlConvSystemP2CL";
    if (smoothed)
        std::cout << " [smoothed]";
    const size_t num_unks_vel= RowIdx.NumUnknowns();
    mN_= new SparseMatBuilderCL<double, SDiagMatrixCL<3> >( &N, num_unks_vel, num_unks_vel);
    if (cplN != 0) {
        cplN->Clear( t);
    }
}

void NonlConvSystemAccumulator_P2CL::finalize_accumulation ()
{
    mN_->Build();
    delete mN_;
#ifndef _PAR
    std::cout << ": " << N.num_nonzeros() << " nonzeros in N, ";
#endif
    std::cout << '\n';
}

void NonlConvSystemAccumulator_P2CL::visit (const TetraCL& tet)
{
    local_setup( tet);
    update_global_system();
}

void NonlConvSystemAccumulator_P2CL::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    n.assign( tet, RowIdx, BndData.Vel);

    ls_loc.assign( tet, lset.Phi, lset.GetBndData());
    vel_loc.assign( tet, vel, BndData.Vel);
    if (equal_signs( ls_loc)) {
        local_onephase.velocity( vel_loc);
        local_onephase.rho( local_twophase.rho( sign( ls_loc[0])));
        local_onephase.setup( T, absdet, loc);
    }
    else {
        if (!smoothed)
            local_twophase.setup( T, absdet, vel_loc, ls_loc, loc);
        else {
            local_smoothed_twophase.velocity( vel_loc);
            local_smoothed_twophase.levelset( ls_loc);
            local_smoothed_twophase.setup( T, absdet, loc);
        }
    }
    
    if (cplN != 0) {
        for (int i= 0; i < 10; ++i) {
            if (!n.WithUnknowns( i)) {
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[i]).GetBndFun();
                dirichlet_val[i]= i<4 ? bf( tet.GetVertex( i)->GetCoord(), t)
                    : bf( GetBaryCenter( *tet.GetEdge( i-4)), t);
            }
        }
    }
}

void NonlConvSystemAccumulator_P2CL::update_global_system ()
{
    SparseMatBuilderCL<double, SDiagMatrixCL<3> >& mN= *mN_;    
    
    for(int i= 0; i < 10; ++i)    // assemble row Numb[i]
        if (n.WithUnknowns( i)) { // dof i is not on a Dirichlet boundary
            for(int j= 0; j < 10; ++j) {
                if (n.WithUnknowns( j)) { // dof j is not on a Dirichlet boundary
                    mN( n.num[i], n.num[j])+= SDiagMatrixCL<3>( loc.C[i][j]);
                }
                else if (cplN != 0) { // right-hand side for eliminated Dirichlet-values
                    add_to_global_vector( cplN->Data, -loc.C[i][j] *dirichlet_val[j], n.num[i]);
                }
            }
        }
}

void InstatNavierStokes2PhaseP2P1CL::SetupNonlinear_P2(MatrixCL& N, const VelVecDescCL* vel, VelVecDescCL* cplN, const LevelsetP2CL& lset, IdxDescCL& RowIdx, double t) const
/// Set up matrix N
{
    // TimerCL time;
    // time.Start();
    NonlConvSystemAccumulator_P2CL accu( Coeff_, MG_, BndData_, *vel, lset, RowIdx, N, cplN, t);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    if ( omp_get_max_threads() > 1)
        accus( MG_.GetColorClasses( RowIdx.TriangLevel()));
    else
        accus( MG_.GetTriangTetraBegin( RowIdx.TriangLevel()), MG_.GetTriangTetraEnd( RowIdx.TriangLevel()));
    // time.Stop();
    // std::cout << "setup: " << time.GetTime() << std::endl;
}


void InstatNavierStokes2PhaseP2P1CL::SetupNonlinear
    ( MLMatDescCL* N, const VelVecDescCL* vel, VelVecDescCL* cplN,
      const LevelsetP2CL& lset, double t) const
/// Couplings with dirichlet BCs are accumulated in cplN,
/// so call cplN->Clear() before if only couplings are needed.
{
    MLMatrixCL::iterator  itN = N->Data.begin();
    MLIdxDescCL::iterator it  = N->RowIdx->begin();
    for (size_t lvl=0; lvl < N->Data.size(); ++lvl, ++itN, ++it)
        SetupNonlinear_P2( *itN, vel, lvl == N->Data.size()-1 ? cplN : 0, lset,*it, t);
}

} // end of namespace DROPS
