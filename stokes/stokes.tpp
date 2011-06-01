/// \file stokes.tpp
/// \brief classes that constitute the stokes-problem
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

#include "num/discretize.h"
#include "misc/problem.h"
#include "num/accumulator.h"
#include "num/quadrature.h"
#include "num/lattice-eval.h"
#include <vector>
#include <numeric>

#ifdef _PAR
#  include "parallel/parallel.h"
#endif

namespace DROPS
{
/**************************************************************************************************
* member functions to create and delete numbering
**************************************************************************************************/
template <class Coeff>
void StokesP2P1CL<Coeff>::CreateNumberingVel( Uint level, MLIdxDescCL* idx, match_fun match)
{
    idx->CreateNumbering( level, MG_, BndData_.Vel, match);
}


template <class Coeff>
void StokesP2P1CL<Coeff>::CreateNumberingPr ( Uint level, MLIdxDescCL* idx, match_fun match)
{
    idx->CreateNumbering( level, MG_, BndData_.Pr, match);
}


/**************************************************************************************************
* member functions to handle with index descriptions
**************************************************************************************************/

template <class Coeff>
  void
  StokesP2P1CL<Coeff>::GetDiscError(instat_vector_fun_ptr LsgVel,
      instat_scalar_fun_ptr LsgPr, double t) const
{
    const Uint lvl= A.GetRowLevel(),
              vidx= A.RowIdx->GetIdx(),
              pidx= B.RowIdx->GetIdx();
    VectorCL lsgvel(A.RowIdx->NumUnknowns());
    VectorCL lsgpr( B.RowIdx->NumUnknowns());
    SVectorCL<3> tmp;

    for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>( MG_).GetTriangVertexBegin( lvl),
         send= const_cast<const MultiGridCL&>( MG_).GetTriangVertexEnd( lvl);
         sit != send; ++sit) {
        if (!BndData_.Vel.IsOnDirBnd( *sit)) {
            tmp= LsgVel(sit->GetCoord(), t);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns( vidx)+i]= tmp[i];
        }
    }

    for (MultiGridCL::const_TriangEdgeIteratorCL sit= const_cast<const MultiGridCL&>( MG_).GetTriangEdgeBegin( lvl),
         send= const_cast<const MultiGridCL&>( MG_).GetTriangEdgeEnd( lvl);
         sit != send; ++sit) {
        if (!BndData_.Vel.IsOnDirBnd( *sit)) {
            tmp= LsgVel( GetBaryCenter( *sit), t);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns( vidx)+i]= tmp[i];
        }
    }
    for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>( MG_).GetTriangVertexBegin( lvl),
         send= const_cast<const MultiGridCL&>( MG_).GetTriangVertexEnd( lvl);
         sit != send; ++sit)
        lsgpr[sit->Unknowns( pidx)]= LsgPr( sit->GetCoord(), t);

    std::cout << "discretization error to check the system (x,y = continuous solution): " << std::endl;
    VectorCL res( A.Data*lsgvel + transp_mul(B.Data, lsgpr) - b.Data);
    std::cout << "|| Ax + BTy - f || = "<< norm( res)<< ", max "<< supnorm( res) << std::endl;
    VectorCL resB( B.Data*lsgvel - c.Data);
    std::cout << "|| Bx - g || = " <<  norm( resB) << ", max " << supnorm( resB) << std::endl;
}



/*****************************************************************************************************
* formulas for   n u m e r i c   i n t e g r a t i o n   on the reference tetrahedron
*****************************************************************************************************/

// ToDo: Kann man auf diese haessliche Formel verzichten??
inline double Quad( const TetraCL& s, instat_scalar_fun_ptr f, int i, int j,  double t= 0.0)
// cubatur formula for int f(x)*phi_i*phi_j dx, exact up to degree 1
{
    double a[5];
    if (i>j) std::swap(i,j);
    switch(i*10+j)
    {
      case  0: a[0]= 1./1260.; a[1]= a[2]= a[3]= 0.; a[4]= 1./630.; break;
      case  1: a[0]= a[1]= -1./8505.; a[2]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case  2: a[0]= a[2]= -1./8505.; a[1]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case  3: a[0]= a[3]= -1./8505.; a[1]= a[2]= 11./136080.; a[4]= 4./8505; break;
      case  4: a[0]= 1./2520.; a[1]= -1./2520.; a[2]= a[3]= 0.; a[4]= -1./630.; break;
      case  5: a[0]= 1./2520.; a[2]= -1./2520.; a[1]= a[3]= 0.; a[4]= -1./630.; break;
      case  6: a[0]= a[3]= 1./8505.; a[1]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case  7: a[0]= 1./2520.; a[3]= -1./2520.; a[1]= a[2]= 0.; a[4]= -1./630.; break;
      case  8: a[0]= a[2]= 1./8505.; a[1]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case  9: a[0]= a[1]= 1./8505.; a[2]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 11: a[1]= 1./1260.; a[0]= a[2]= a[3]= 0.; a[4]= 1./630.; break;
      case 12: a[1]= a[2]= -1./8505.; a[0]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case 13: a[1]= a[3]= -1./8505.; a[0]= a[2]= 11./136080.; a[4]= 4./8505; break;
      case 14: a[1]= 1./2520.; a[0]= -1./2520.; a[2]= a[3]= 0.; a[4]= -1./630.; break;
      case 15: a[1]= a[3]= 1./8505.; a[0]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 16: a[1]= 1./2520.; a[2]= -1./2520.; a[0]= a[3]= 0.; a[4]= -1./630.; break;
      case 17: a[1]= a[2]= 1./8505.; a[0]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 18: a[1]= 1./2520.; a[3]= -1./2520.; a[0]= a[2]= 0.; a[4]= -1./630.; break;
      case 19: a[0]= a[1]= 1./8505.; a[2]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 22: a[2]= 1./1260.; a[0]= a[1]= a[3]= 0.; a[4]= 1./630.; break;
      case 23: a[2]= a[3]= -1./8505.; a[0]= a[1]= 11./136080.; a[4]= 4./8505; break;
      case 24: a[2]= a[3]= 1./8505.; a[0]= a[1]= -19./68040.; a[4]= -1./486.; break;
      case 25: a[2]= 1./2520.; a[0]= -1./2520.; a[1]= a[3]= 0.; a[4]= -1./630.; break;
      case 26: a[2]= 1./2520.; a[1]= -1./2520.; a[0]= a[3]= 0.; a[4]= -1./630.; break;
      case 27: a[2]= a[1]= 1./8505.; a[0]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 28: a[2]= a[0]= 1./8505.; a[1]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 29: a[2]= 1./2520.; a[3]= -1./2520.; a[0]= a[1]= 0.; a[4]= -1./630.; break;
      case 33: a[3]= 1./1260.; a[0]= a[1]= a[2]= 0.; a[4]= 1./630.; break;
      case 34: a[3]= a[2]= 1./8505.; a[0]= a[1]= -19./68040.; a[4]= -1./486.; break;
      case 35: a[3]= a[1]= 1./8505.; a[0]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 36: a[3]= a[0]= 1./8505.; a[1]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 37: a[3]= 1./2520.; a[0]= -1./2520.; a[1]= a[2]= 0.; a[4]= -1./630.; break;
      case 38: a[3]= 1./2520.; a[1]= -1./2520.; a[0]= a[2]= 0.; a[4]= -1./630.; break;
      case 39: a[3]= 1./2520.; a[2]= -1./2520.; a[0]= a[1]= 0.; a[4]= -1./630.; break;
      case 44: a[0]= a[1]= 37./17010.; a[2]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 45: a[0]= 1./972.; a[1]= a[2]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 46: a[1]= 1./972.; a[0]= a[2]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 47: a[0]= 1./972.; a[1]= a[3]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 48: a[1]= 1./972.; a[0]= a[3]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 49: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 55: a[0]= a[2]= 37./17010.; a[1]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 56: a[2]= 1./972.; a[0]= a[1]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 57: a[0]= 1./972.; a[2]= a[3]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 58: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 59: a[2]= 1./972.; a[0]= a[3]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 66: a[1]= a[2]= 37./17010.; a[0]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 67: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 68: a[1]= 1./972.; a[2]= a[3]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 69: a[2]= 1./972.; a[1]= a[3]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 77: a[0]= a[3]= 37./17010.; a[1]= a[2]= -17./17010.; a[4]= 88./8505.; break;
      case 78: a[3]= 1./972.; a[0]= a[1]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 79: a[3]= 1./972.; a[0]= a[2]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 88: a[1]= a[3]= 37./17010.; a[0]= a[2]= -17./17010.; a[4]= 88./8505.; break;
      case 89: a[3]= 1./972.; a[1]= a[2]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 99: a[2]= a[3]= 37./17010.; a[0]= a[1]= -17./17010.; a[4]= 88./8505.; break;
      default: throw DROPSErrCL("Quad(i,j): no such shape function");
    }
    double sum= a[4]*f(GetBaryCenter(s),t);
    for(Uint i=0; i<4; ++i)
        sum+= a[i]*f(s.GetVertex(i)->GetCoord(),t);
    return sum;
}

/*****************************************************************************************************
*                                   setup routines
*****************************************************************************************************/

/// \brief Raw data for "system 1", both for one phase and two phases.
///
/// scalar-valued mass-matrix, scalar-valued mu-Laplacian, genuinely tensor-valued part of the deformation tensor and the integrals of \f$\rho\phi_i\f$ for the gravitation as load-vector
/// \todo: Precise description
struct LocalStokesSystem1DataCL
{
    double         M [10][10];
    double         A [10][10];
};

/// \brief Setup of the local "system 1" on a tetra in a single phase.
class LocalStokesSystem1OnePhase_P2CL
{
  private:
    double mu_;
    double rho_;

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    const SVectorCL<Quad2DataCL::NumNodesC> Ones;

  public:
    LocalStokesSystem1OnePhase_P2CL (double muarg= 0., double rhoarg= 0.)
        : mu_( muarg), rho_( rhoarg), Ones( 1.)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void   mu  (double new_mu)        { mu_= new_mu; }
    double mu  ()               const { return mu_; }
    void   rho (double new_rho)       { rho_= new_rho; }
    double rho ()               const { return rho_; }

    void setup (const SMatrixCL<3,3>& T, double absdet, LocalStokesSystem1DataCL& loc)
    {
        P2DiscCL::GetGradients( Grad, GradRef, T);
        for (Uint i= 0; i < 10; ++i) {
            for (Uint j= 0; j <= i; ++j) {
                // M: As we are not at the phase-boundary this is exact.
                loc.M[j][i]= rho()*P2DiscCL::GetMass( j, i)*absdet;
                loc.A[j][i]= mu()*quad( dot( Grad[i], Grad[j]), absdet, make_Quad2Data());

                if (i != j) { // The local matrices coupM, coupA, coupAk are symmetric.
                    loc.M[i][j]= loc.M[j][i];
                    loc.A[i][j]= loc.A[j][i];
                }
            }
        }
    }

};


/// \brief Accumulator to set up the matrices A, M and, if requested the right-hand side b and cplM, cplA for Stokes flow.
template< class CoeffT>
class StokesSystem1Accumulator_P2CL : public TetraAccumulatorCL
{
  private:
    const CoeffT& Coeff;
    const StokesBndDataCL& BndData;
    double t;

    IdxDescCL& RowIdx;
    MatrixCL& A;
    MatrixCL& M;
    VecDescCL* cplA;
    VecDescCL* cplM;
    VecDescCL* b;

    SparseMatBuilderCL<double, SDiagMatrixCL<3> >* mA_;
    SparseMatBuilderCL<double, SDiagMatrixCL<3> >* mM_;

    LocalStokesSystem1OnePhase_P2CL local_onephase; ///< used on tetras in a single phase
    LocalStokesSystem1DataCL loc; ///< Contains the memory, in which the local operators are set up; former coupM, coupA, coupAk, rho_phi.

    LocalNumbP2CL n; ///< global numbering of the P2-unknowns

    SMatrixCL<3,3> T;
    double det, absdet;
    LocalP2CL<> ls_loc;

    Quad2CL<Point3DCL> rhs;
    Point3DCL loc_b[10], dirichlet_val[10]; ///< Used to transfer boundary-values from local_setup() update_global_system().

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and, if required, the Dirichlet-values needed to eliminate the boundary-dof from the global system.
    void local_setup (const TetraCL& tet);
    ///\brief Update the global system.
    void update_global_system ();

  public:
    StokesSystem1Accumulator_P2CL (const CoeffT& Coeff, const StokesBndDataCL& BndData_,
        IdxDescCL& RowIdx_, MatrixCL& A_, MatrixCL& M_,
        VecDescCL* b_, VecDescCL* cplA_, VecDescCL* cplM_, double t);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);
    
    TetraAccumulatorCL* clone (){ return new StokesSystem1Accumulator_P2CL ( *this); };
};

template< class CoeffT>
StokesSystem1Accumulator_P2CL<CoeffT>::StokesSystem1Accumulator_P2CL (const CoeffT& Coeff_, const StokesBndDataCL& BndData_,
    IdxDescCL& RowIdx_, MatrixCL& A_, MatrixCL& M_,
    VelVecDescCL* b_, VelVecDescCL* cplA_, VelVecDescCL* cplM_, double t_)
    : Coeff( Coeff_), BndData( BndData_), t( t_),
      RowIdx( RowIdx_), A( A_), M( M_), cplA( cplA_), cplM( cplM_), b( b_)
{}

template< class CoeffT>
void StokesSystem1Accumulator_P2CL<CoeffT>::begin_accumulation ()
{
    const size_t num_unks_vel= RowIdx.NumUnknowns();
    mA_= new SparseMatBuilderCL<double, SDiagMatrixCL<3> >( &A, num_unks_vel, num_unks_vel);
    mM_= new SparseMatBuilderCL<double, SDiagMatrixCL<3> >( &M, num_unks_vel, num_unks_vel);
    if (b != 0) {
        b->Clear( t);
        cplM->Clear( t);
        cplA->Clear( t);
    }
}

template< class CoeffT>
void StokesSystem1Accumulator_P2CL<CoeffT>::finalize_accumulation ()
{
    mA_->Build();
    delete mA_;
    mM_->Build();
    delete mM_;
}

template< class CoeffT>
void StokesSystem1Accumulator_P2CL<CoeffT>::visit (const TetraCL& tet)
{
    local_setup( tet);
    update_global_system();
}

template< class CoeffT>
void StokesSystem1Accumulator_P2CL<CoeffT>::local_setup (const TetraCL& tet)
{
    GetTrafoTr( T, det, tet);
    absdet= std::fabs( det);

    rhs.assign( tet, Coeff.f, t);
    n.assign( tet, RowIdx, BndData.Vel);

    local_onephase.mu(  Coeff.nu);
    local_onephase.rho( 1.0);
    local_onephase.setup( T, absdet, loc);

    if (b != 0) {
        for (int i= 0; i < 10; ++i) {
            if (!n.WithUnknowns( i)) {
                typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
                bnd_val_fun bf= BndData.Vel.GetBndSeg( n.bndnum[i]).GetBndFun();
                dirichlet_val[i]= i<4 ? bf( tet.GetVertex( i)->GetCoord(), t)
                    : bf( GetBaryCenter( *tet.GetEdge( i-4)), t);
            }
            else
                loc_b[i]= rhs.quadP2( i, absdet);
        }
    }
}

template< class CoeffT>
void StokesSystem1Accumulator_P2CL<CoeffT>::update_global_system ()
{
    SparseMatBuilderCL<double, SDiagMatrixCL<3> >& mA= *mA_;
    SparseMatBuilderCL<double, SDiagMatrixCL<3> >& mM= *mM_;

    for(int i= 0; i < 10; ++i)    // assemble row Numb[i]
        if (n.WithUnknowns( i)) { // dof i is not on a Dirichlet boundary
            for(int j= 0; j < 10; ++j) {
                if (n.WithUnknowns( j)) { // dof j is not on a Dirichlet boundary
                    mA( n.num[i], n.num[j])+= SDiagMatrixCL<3>( loc.A[j][i]);
                    mM( n.num[i], n.num[j])+= SDiagMatrixCL<3>( loc.M[j][i]);
                }
                else if (b != 0) { // right-hand side for eliminated Dirichlet-values
                    add_to_global_vector( cplA->Data, -loc.A[j][i]*dirichlet_val[j], n.num[i]);
                    add_to_global_vector( cplM->Data, -loc.M[j][i]*dirichlet_val[j], n.num[i]);
                }
            }
            if (b != 0) // assemble the right-hand side
                add_to_global_vector( b->Data, loc_b[i], n.num[i]);
       }
}

template< class CoeffT>
void SetupSystem1_P2( const MultiGridCL& MG_, const CoeffT& Coeff_, const StokesBndDataCL& BndData_, MatrixCL& A, MatrixCL& M,
                      VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, IdxDescCL& RowIdx, double t)
/// Set up matrices A, M and rhs b (depending on phase bnd)
{
    StokesSystem1Accumulator_P2CL<CoeffT> accu( Coeff_, BndData_, RowIdx, A, M, b, cplA, cplM, t);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accus( MG_.GetTriangTetraBegin( RowIdx.TriangLevel()), MG_.GetTriangTetraEnd( RowIdx.TriangLevel()));
}


template< class CoeffT>
void StokesP2P1CL<CoeffT>::SetupSystem1( MLMatDescCL* A, MLMatDescCL* M, VecDescCL* b, VecDescCL* cplA, VecDescCL* cplM, double t) const
{
    MLMatrixCL::iterator itA = A->Data.begin();
    MLMatrixCL::iterator itM = M->Data.begin();
    MLIdxDescCL::iterator it = A->RowIdx->begin();
    for (size_t lvl=0; lvl < A->Data.size(); ++lvl, ++itA, ++itM, ++it) {
#ifndef _PAR
        std::cout << "entering SetupSystem1: " << it->NumUnknowns() << " unknowns ";
#endif
        SetupSystem1_P2 ( MG_, Coeff_, BndData_, *itA, *itM, lvl == A->Data.size()-1 ? b : 0, cplA, cplM, *it, t);
#ifndef _PAR
        std::cout << itA->num_nonzeros() << " nonzeros in A, "
                  << itM->num_nonzeros() << " nonzeros in M! " << std::endl;
#endif
    }
}


/// \brief Shared data for "system 2" between P1 and P1X.
/// All members are setup by System2Accumulator_P2P1CL::visit.
struct LocalSystem2_sharedDataCL
{
    IdxT          prNumb[4];  ///< global numbering of the P1-unknowns
    LocalNumbP2CL n;          ///< global numbering of the P2-unknowns

    Point3DCL dirichlet_val[10]; ///< Dirichlet values, filled in only on the Dirichlet-boundary.

    SparseMatBuilderCL<double, SMatrixCL<1,3> >* mB;
    VecDescCL*                                   c;

    SMatrixCL<3,3> T;
    double         absdet;
};


/// \brief Accumulator to set up the matrix B and, if requested the right-hand side C for two-phase flow.
template<class CoeffT>
class System2Accumulator_P2P1CL : public TetraAccumulatorCL
{
  private:
    const CoeffT& coeff;
    const StokesBndDataCL& BndData;
    const double t;

    const IdxDescCL& RowIdx;
    const IdxDescCL& ColIdx;
    MatrixCL& B;

    LocalSystem2_sharedDataCL loc;

    Quad2CL<Point3DCL> GradRef[10],
                       Grad[10];
    SMatrixCL<1,3>     locB[10][4];

    ///\brief Computes the mapping from local to global data "n", the local matrices in loc and, if required, the Dirichlet-values needed to eliminate the boundary-dof from the global system.
    void local_setup ();
    ///\brief Update the global system.
    void update_global_system ();

  public:
    System2Accumulator_P2P1CL (const CoeffT& coeff_arg, const StokesBndDataCL& BndData_arg,
        const IdxDescCL& RowIdx_arg, const IdxDescCL& ColIdx_arg,
        MatrixCL& B_arg, VecDescCL* c_arg, double t_arg);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    void visit (const TetraCL& sit);

    LocalSystem2_sharedDataCL& GetLocalData() { return loc; }
    
    TetraAccumulatorCL* clone (){ return new System2Accumulator_P2P1CL ( *this); };
};

template< class CoeffT>
System2Accumulator_P2P1CL<CoeffT>::System2Accumulator_P2P1CL ( const CoeffT& coeff_arg, const StokesBndDataCL& BndData_arg,
    const IdxDescCL& RowIdx_arg, const IdxDescCL& ColIdx_arg,
    MatrixCL& B_arg, VecDescCL* c_arg, double t_arg)
    : coeff( coeff_arg), BndData( BndData_arg), t( t_arg), RowIdx( RowIdx_arg), ColIdx( ColIdx_arg), B( B_arg)
{
    loc.c = c_arg;
    P2DiscCL::GetGradientsOnRef( GradRef);
}

template< class CoeffT>
void System2Accumulator_P2P1CL<CoeffT>::begin_accumulation ()
{
    loc.mB = new SparseMatBuilderCL<double, SMatrixCL<1,3> > ( &B, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    if (loc.c != 0) loc.c->Clear( t);
}

template< class CoeffT>
void System2Accumulator_P2P1CL<CoeffT>::finalize_accumulation ()
{
    loc.mB->Build();
    delete loc.mB;
}

template< class CoeffT>
void System2Accumulator_P2P1CL<CoeffT>::visit (const TetraCL& tet)
{
    double det;
    GetTrafoTr( loc.T, det, tet);
    P2DiscCL::GetGradients( Grad, GradRef, loc.T);
    loc.absdet= std::fabs( det);
    loc.n.assign( tet, ColIdx, BndData.Vel);
    GetLocalNumbP1NoBnd( loc.prNumb, tet, RowIdx);

    if (loc.c != 0) {
        typedef StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_val_fun;
        for (int i= 0; i < 10; ++i)
            if (!loc.n.WithUnknowns( i)) {
                bnd_val_fun bf= BndData.Vel.GetBndSeg( loc.n.bndnum[i]).GetBndFun();
                loc.dirichlet_val[i]= i<4 ? bf( tet.GetVertex( i)->GetCoord(), t)
                    : bf( GetBaryCenter( *tet.GetEdge( i-4)), t);
            }
    }
    local_setup();
    update_global_system();
}

template< class CoeffT>
void System2Accumulator_P2P1CL<CoeffT>::local_setup ()
{
    // b(i,j) =  -\int psi_i * div( phi_j)
    for(int vel=0; vel<10; ++vel) {
        for(int pr=0; pr<4; ++pr)
            locB[vel][pr]= SMatrixCL<1,3>( quad( Grad[vel], loc.absdet, Quad2Data_Mul_P1_CL(), pr));
    }
}

template< class CoeffT>
void System2Accumulator_P2P1CL<CoeffT>::update_global_system ()
{
    SparseMatBuilderCL<double, SMatrixCL<1,3> >& mB= *loc.mB;

    for(int vel=0; vel<10; ++vel) {
        if (loc.n.WithUnknowns( vel))
            for(int pr=0; pr<4; ++pr)
                mB( loc.prNumb[pr], loc.n.num[vel])-= locB[vel][pr];
        else if (loc.c != 0) { // put coupling on rhs
            for(int pr=0; pr<4; ++pr)
                loc.c->Data[loc.prNumb[pr]]+= inner_prod( locB[vel][pr], loc.dirichlet_val[vel]); // operator* returns SVectorCL<1>.
        }
    }
}

template <class CoeffT>
void SetupSystem2_P2P1( const MultiGridCL& MG, const CoeffT& coeff, const StokesBndDataCL& BndData, MatrixCL* B, VecDescCL* c,
        IdxDescCL* RowIdx, IdxDescCL* ColIdx, double t)
/// Set up matrices B and rhs c
{
    System2Accumulator_P2P1CL<CoeffT> accu( coeff, BndData, *RowIdx, *ColIdx, *B, c, t);
    TetraAccumulatorTupleCL accus;
    accus.push_back( &accu);
    accus( MG.GetTriangTetraBegin( RowIdx->TriangLevel()), MG.GetTriangTetraEnd( RowIdx->TriangLevel()));
}

template <class CoeffT>
void StokesP2P1CL<CoeffT>::SetupSystem2( MLMatDescCL* B, VecDescCL* c, double t) const
// Set up matrix B and rhs c
{
    MLMatrixCL::iterator     itB   = B->Data.begin();
    MLIdxDescCL::iterator    itRow = B->RowIdx->begin();
    MLIdxDescCL::iterator    itCol = B->ColIdx->begin();
    if ( B->RowIdx->size() == 1 || B->ColIdx->size() == 1)
    { // setup B only on finest level, if row or column index has only 1 level
        itCol = B->ColIdx->GetFinestIter();
        itRow = B->RowIdx->GetFinestIter();
        itB   = B->Data.GetFinestIter();
    }
    for (; itB!=B->Data.end() && itRow!=B->RowIdx->end() && itCol!=B->ColIdx->end(); ++itB, ++itRow, ++itCol)
    {
#ifndef _PAR
        std::cout << "entering SetupSystem2: " << itRow->NumUnknowns() << " prs, " << itCol->NumUnknowns() << " vels. ";
#endif
        VecDescCL* rhsPtr= itB==B->Data.GetFinestIter() ? c : 0; // setup rhs only on finest level
        SetupSystem2_P2P1 ( MG_, Coeff_, BndData_, &(*itB), rhsPtr, &(*itRow), &(*itCol), t);
#ifndef _PAR
        std::cout << itB->num_nonzeros() << " nonzeros in B!" << std::endl;
#endif
    }
}


template <class CoeffT>
void SetupPrStiff_P1_Nolst( const MultiGridCL& MG, const CoeffT& Coeff, MatrixCL& A_pr, IdxDescCL& RowIdx, IdxDescCL& ColIdx)
{
    MatrixBuilderCL A( &A_pr, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    const Uint lvl= RowIdx.TriangLevel();
    const Uint idx= RowIdx.GetIdx();
    SMatrixCL<3,4> G;
    double coup[4][4];
    double det, absdet, IntRhoInv;
    IdxT UnknownIdx[4];

    const double rho_inv= 1./Coeff.rho;


    for (MultiGridCL::const_TriangTetraIteratorCL sit= MG.GetTriangTetraBegin( lvl),
         send= MG.GetTriangTetraEnd( lvl); sit != send; ++sit)
    {
        P1DiscCL::GetGradients( G,det,*sit);
        absdet= std::fabs( det);

        IntRhoInv= absdet/6*rho_inv;

        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )*IntRhoInv;
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
void StokesP2P1CL<Coeff>::SetupPrStiff( MLMatDescCL* A_pr) const
/// Needed for preconditioning of the Schur complement. Uses natural
/// boundary conditions for the pressure unknowns.
{
    MLMatrixCL::iterator itM = A_pr->Data.begin();
    MLIdxDescCL::iterator itRowIdx = A_pr->RowIdx->begin();
    MLIdxDescCL::iterator itColIdx = A_pr->ColIdx->begin();
    for (size_t lvl=0; lvl < A_pr->Data.size(); ++lvl, ++itM, ++itRowIdx, ++itColIdx)
        SetupPrStiff_P1_Nolst( MG_, Coeff_, *itM, *itRowIdx, *itColIdx);
}

//====================================================================


template <class CoeffT>
  void
  SetupPrMass_P2P1( const MultiGridCL& MG, const CoeffT&, MatrixCL& matM, IdxDescCL& RowIdx)
/// Sets up the mass matrix for the pressure
{
    const IdxT num_unks_pr=  RowIdx.NumUnknowns();

    MatrixBuilderCL M_pr(&matM, num_unks_pr,  num_unks_pr);

    const Uint lvl    = RowIdx.TriangLevel();
    const Uint pidx   = RowIdx.GetIdx();

    IdxT prNumb[4];

    // compute all couplings between HatFunctions on verts:
    // M( i, j) = int ( psi_i*psi_j, T_ref) * absdet

    for (MultiGridCL::const_TriangTetraIteratorCL sit=MG.GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6.;

        for(int i=0; i<4; ++i)
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                M_pr( prNumb[i], prNumb[j])+= P1DiscCL::GetMass(i,j) * absdet;
    }
    M_pr.Build();
}

template <class Coeff>
void StokesP2P1CL<Coeff>::SetupPrMass(MLMatDescCL* matM) const
{
    MLMatrixCL::iterator  itM   = matM->Data.begin();
    MLIdxDescCL::iterator itRow = matM->RowIdx->begin();
    for ( size_t lvl=0; lvl < matM->Data.size(); ++lvl, ++itM, ++itRow)
        SetupPrMass_P2P1( MG_, Coeff_, *itM, *itRow);
}


template <class Coeff>
void StokesP2P1CL<Coeff>::SetupInstatRhs( VelVecDescCL* vecA, VelVecDescCL* vecB,
                                                VelVecDescCL* vecI, double tA,
                                                VelVecDescCL* vecf, double tf) const
/// Sets up the couplings with hat fcts on bnd (for matrices A, I) and discretizes the PDE coeff f(t)
{
    vecA->Clear( tA);
    vecB->Clear( tA);
    vecf->Clear( tf);
    vecI->Clear( tA);

#ifndef _PAR
    __UNUSED__ const IdxT allnum_unks_vel= vecA->RowIdx->NumUnknowns();
#else
    __UNUSED__ const IdxT allnum_unks_vel= ProcCL::GlobalSum(vecA->RowIdx->NumUnknowns());
#endif
    Comment("entering StokesP2P1CL::SetupInstatSystem: "<<allnum_unks_vel<< " velocity unknowns.\n", DebugDiscretizeC);


    VectorCL& a    = vecA->Data;
    VectorCL& c    = vecB->Data;
    VectorCL& f    = vecf->Data;
    VectorCL& id   = vecI->Data;
    const Uint lvl = vecA->GetLevel();
    const Uint vidx= vecA->RowIdx->GetIdx(),
               pidx= vecB->RowIdx->GetIdx();

    IdxT Numb[10], prNumb[4];
    bool IsOnDirBnd[10];

    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns

    Quad2CL<Point3DCL> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    Quad2CL<Point3DCL> rhs;
    SMatrixCL<3,3> T;
    double coup[10][10], coupMass[10][10];
    double det, absdet;
    SVectorCL<3> tmp;

    P2DiscCL::GetGradientsOnRef(GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        P2DiscCL::GetGradients(Grad, GradRef, T);
        absdet= std::fabs(det);
        rhs.assign( *sit, Coeff::f, tf);

        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= BndData_.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= BndData_.Vel.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(vidx);
        }

        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                const double c= Coeff_.nu * Quad2CL<>( dot( Grad[i], Grad[j])).quad(absdet);
//                c+= Quad(*sit, &Coeff::q, i, j)*absdet;
                coup[i][j]= c;
                coup[j][i]= c;

                const double cM= P2DiscCL::GetMass( j, i)*absdet;
                coupMass[i][j]= cM;
                coupMass[j][i]= cM;
            }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (IsOnDirBnd[j]) // vert/edge j is on a Dirichlet boundary
                    { // coupling with vert/edge j on right-hand-side
                        tmp= j<4 ? BndData_.Vel.GetDirBndValue(*sit->GetVertex(j), tA)
                                 : BndData_.Vel.GetDirBndValue(*sit->GetEdge(j-4), tA);
                        a[Numb[i]]-=          coup[j][i] * tmp[0];
                        a[Numb[i]+stride]-=   coup[j][i] * tmp[1];
                        a[Numb[i]+2*stride]-= coup[j][i] * tmp[2];

                        const double val= coupMass[i][j];
                        id[Numb[i]]-=          val*tmp[0];
                        id[Numb[i]+stride]-=   val*tmp[1];
                        id[Numb[i]+2*stride]-= val*tmp[2];
                    }
                }
                tmp= rhs.quadP2( i, absdet);//P2DiscCL::Quad(*sit, Coeff::f, i, tf)*absdet;
                f[Numb[i]]+=          tmp[0];
                f[Numb[i]+stride]+=   tmp[1];
                f[Numb[i]+2*stride]+= tmp[2];

                if ( i<4 ? BndData_.Vel.IsOnNatBnd(*sit->GetVertex(i))
                         : BndData_.Vel.IsOnNatBnd(*sit->GetEdge(i-4)) ) // vert/edge i is on natural boundary
                {
                    Uint face;
                    for (int f=0; f < (i<4?3:2); ++f)
                    {
                        face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);
                        if ( sit->IsBndSeg(face))
                        {   // TODO: FIXME: Hier muss doch eigentlich eine 2D-Integrationsformel fuer P2-Elemente stehen, oder?
/*                            tmp= Quad2D(*sit, face, i, BndData_.Vel.GetSegData(sit->GetBndIdx(face)).GetBndFun(), tA);
                            a[Numb[i]]+=          tmp[0];
                            a[Numb[i]+stride]+=   tmp[1];
                            a[Numb[i]+2*stride]+= tmp[2];
*/                        }
                    }
                }
            }
        // Setup B:   b(i,j) =  -\int psi_i * div( phi_j)
        for(int vel=0; vel<10; ++vel)
        {
            if (IsOnDirBnd[vel])
            {    // put coupling on rhs
                const Point3DCL bndval= vel<4 ? BndData_.Vel.GetDirBndValue( *sit->GetVertex(vel), tA)
                                              : BndData_.Vel.GetDirBndValue( *sit->GetEdge(vel-4), tA);
                for(int pr=0; pr<4; ++pr)
                {
                    // numeric integration is exact: psi_i * div( phi_j) is of degree 2 !
                    tmp= Grad[vel].quadP1( pr, absdet);
                    c[prNumb[pr]]+= inner_prod( tmp, bndval);
                }
            }
        }
    }
}


template <class Coeff>
void StokesP2P1CL<Coeff>::InitVel(VelVecDescCL* vec, instat_vector_fun_ptr LsgVel, double t0) const
{
    VectorCL& lsgvel= vec->Data;
    Uint lvl        = vec->GetLevel(),
         vidx       = vec->RowIdx->GetIdx();

    SVectorCL<3> tmp;

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!BndData_.Vel.IsOnDirBnd(*sit))
        {
            tmp= LsgVel(sit->GetCoord(), t0);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= tmp[i];
        }
    }

    for (MultiGridCL::const_TriangEdgeIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangEdgeBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangEdgeEnd(lvl);
         sit != send; ++sit)
    {
        if (!BndData_.Vel.IsOnDirBnd(*sit))
        {
            tmp= LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t0);
            for(int i=0; i<3; ++i)
                lsgvel[sit->Unknowns(vidx)+i]= tmp[i];
        }
    }
}


// CheckSolution
template <class Coeff>
void StokesP2P1CL<Coeff>::CheckSolution(const VelVecDescCL* lsgvel, const VecDescCL* lsgpr,
    instat_vector_fun_ptr LsgVel, instat_matrix_fun_ptr DLsgVel, instat_scalar_fun_ptr LsgPr, bool is_stat) const
{
    double t = lsgpr->t;
#ifdef _PAR
    const ExchangeCL& exV = vel_idx.GetEx();
    const ExchangeCL& exP = pr_idx.GetEx();
#endif
    Uint lvl=lsgvel->GetLevel();

    if (is_stat)
    {
        VectorCL res1( A.Data*lsgvel->Data + transp_mul( B.Data, lsgpr->Data ) - b.Data);
        VectorCL res2( B.Data*lsgvel->Data - c.Data);

        #ifndef _PAR
            const double norm_res1      = norm(res1);
            const double norm_res2      = norm(res2);
            const double norm_sup_res_1 = supnorm(res1);
            const double norm_sup_res_2 = supnorm(res2);
        #else
            VectorCL res1_acc(res1);
            VectorCL res2_acc(res2);
            const double norm_res1      = std::sqrt(exV.ParDotAcc(res1_acc,res1));
            const double norm_res2      = std::sqrt(exP.ParDotAcc(res2_acc,res2));
            const double norm_sup_res_1 = ProcCL::GlobalMax(supnorm(res1_acc));
            const double norm_sup_res_2 = ProcCL::GlobalMax(supnorm(res2_acc));
        #endif

        IF_MASTER
            std::cout << "\nCheck the solution..."
                      << "\n|| Ax + BTy - F || = " << norm_res1 << ", max. " << norm_sup_res_1
                      << "\n||       Bx - G || = " << norm_res2 << ", max. " << norm_sup_res_2
                      << '\n' << std::endl;
    }

    double L2_div= 0;
    SMatrixCL<3,3> T;
    double det;

    // Compute the pressure-coefficient in direction of 1/std::sqrt(meas(Omega)), which eliminates
    // the allowed offset of the pressure by setting it to 0.
    double MW_pr= 0, vol= 0;
    Quad5CL<> q5_pr, q5_pr_exact;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl),
         send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl); sit != send; ++sit)
    {
        const double volT= sit->GetVolume();
        LocalP1CL<> loc_pr( *sit, make_P1Eval(MG_,BndData_.Pr,*lsgpr));
        q5_pr.assign(loc_pr);
        q5_pr_exact.assign(*sit,LsgPr,t);
        MW_pr+= Quad5CL<> (q5_pr-q5_pr_exact).quad(6.*volT);
        vol+= volT;
    }
#ifdef _PAR
    vol   = ProcCL::GlobalSum(vol);
    MW_pr = ProcCL::GlobalSum(MW_pr);
#endif
    const double c_pr= MW_pr/vol;
    IF_MASTER
        std::cout << "\n  mean value of exact pressure - discrete pressure is " << c_pr << ", volume of cube is " << vol << std::endl;

    // Some norms of velocities: u_h - u
    // number of nodes for Quad5CL rule is 15 (see discretize.h Quad5_DataCL NumNodesC =15)
    double Frob_Dvel(0.0), L2_vel(0.0);
    double L2_pr(0.0), H1_vel(0.0);

    Quad5CL<Point3DCL> q5_vel, q5_vel_exact;
    Quad5CL<Point3DCL> Grad[10], GradRef[10];
    Quad5CL<SMatrixCL<3,3> > q5_dvel, q5_dvel_exact;
    P2DiscCL::GetGradientsOnRef( GradRef);
    for (MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl),
        send= const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl); sit != send; ++sit)
    {
         GetTrafoTr(T,det,*sit);
         const double absdet= std::fabs(det);
         LocalP2CL<Point3DCL> loc_vel(*sit, make_P2Eval(MG_,BndData_.Vel,*lsgvel));
         LocalP1CL<> loc_pr(*sit, make_P1Eval(MG_,BndData_.Pr,*lsgpr));
         q5_pr.assign(loc_pr);
         q5_pr_exact.assign(*sit,LsgPr,t);
         q5_vel.assign(loc_vel);
         q5_vel_exact.assign(*sit,LsgVel,t);
         if( DLsgVel != NULL)
             q5_dvel_exact.assign(*sit,DLsgVel,t);
         Quad5CL<Point3DCL> q5_vel_diff( q5_vel-q5_vel_exact);
         L2_vel += Quad5CL<> (dot(q5_vel_diff,q5_vel_diff)).quad(absdet);

         L2_pr  += Quad5CL<> (std::pow(q5_pr-q5_pr_exact-c_pr,2)).quad(absdet);

	 P2DiscCL::GetGradients( Grad, GradRef, T);
	 q5_dvel= SMatrixCL<3,3>();
	 for (int i=0; i<10; i++)
	 {
             q5_dvel += outer_product(loc_vel[i],Grad[i]);
         }
         if( DLsgVel != NULL)
	 {
	     Quad5CL< SMatrixCL<3,3> > q5_dvel_diff( q5_dvel-q5_dvel_exact);
	     Frob_Dvel += Quad5CL<> (frobenius_norm_sq(q5_dvel_diff)).quad(absdet);
	     L2_div += Quad5CL<> (trace(q5_dvel)).quad(absdet);
	 }
         
     }
     L2_pr = std::sqrt(L2_pr);
     if( DLsgVel != NULL)
     {
         H1_vel = std::sqrt(L2_vel+Frob_Dvel); //L2_vel is squared here, Frob_Dvel is squared.
         Frob_Dvel = std::sqrt(Frob_Dvel);     //Frob_Dvel is the true value.
	 L2_div = std::sqrt(std::fabs(L2_div));
     }
     L2_vel = std::sqrt(L2_vel);               //L2_vel is the true value.  


     if( DLsgVel != NULL)
         std::cout << "|| (u_h, p_h) - (u, p) ||_X = " << H1_vel
                   << ", || u_h - u ||_L2 = " <<  L2_vel << ", || Du_h - Du ||_L2 = " << Frob_Dvel
                   << ", || p_h - p ||_L2 = " << L2_pr
                   << "\n|| div x ||_L2 = " << L2_div << '\n' << std::endl;
     else
         std::cout << ", || u_h - u ||_L2 = " <<  L2_vel 
                   << ", || p_h - p ||_L2 = " << L2_pr << '\n' << std::endl;
}


#ifndef _PAR
template <class Coeff>
  double
  StokesP2P1CL<Coeff>::ResidualErrEstimator(const TetraCL& s, const const_DiscPrSolCL& pr,
      const const_DiscVelSolCL& vel, double t)
{
    double det;
    SMatrixCL<3,3> M;
    GetTrafoTr(M, det, s);
    const double absdet= std::fabs(det);
    const double vol= absdet/6.;
    Point3DCL cc; // unused; circumcenter of T dummy
    double hT; // radius of circumcircle of T
    circumcircle(s, cc, hT);
    // P_0(f) := (f, 1_T)_T * 1_T = int(f, T)/|T|
    const SVectorCL<3> P0f= Quad3PosWeightsCL::Quad(s, Coeff::f, t)*6.;

    const Uint numpts= Quad3PosWeightsCL::GetNumPoints();
    double* vals= new double[numpts];
    double err_sq= 0.0;

    // Now eta_T^2:
    // Put velocity degrees of freedom in veldof
    std::vector< SVectorCL<3> > veldof( 10);
    vel.GetDoF(s, veldof);
    // Put pressure degrees of freedom in prdof
    std::vector< double > prdof(NumVertsC);
    pr.GetDoF(s, prdof);

    // || div(u_h) ||_L2(T) squared; the integrand is linear, thus we only need a quadrature-formula,
    // which is exact up to degree 2 (squares in L2-norm...); optimize this later.
    for (Uint i=0; i<numpts; ++i)
    {
        double tmp= 0.;
        for (Uint j=0; j<10; ++j)
        {
            tmp+= inner_prod( veldof[j], M*FE_P2CL::DHRef(j, Quad3PosWeightsCL::GetPoints()[i][0], Quad3PosWeightsCL::GetPoints()[i][1], Quad3PosWeightsCL::GetPoints()[i][2]) );
        }
        vals[i]= tmp*tmp;
    }
    err_sq+= Quad3PosWeightsCL::Quad(vals)*absdet;
    delete[] vals;

    // hT^2*int((-laplace(u) + grad(p) - P0f)^2, T) -- the sign of grad(p) is due to the implemented version of stokes eq.
    // the integrand is a constant for P2P1-discretisation...
    SVectorCL<3> tmp= -P0f;
    tmp+= M*(  prdof[0]*FE_P1CL::DH0Ref() + prdof[1]*FE_P1CL::DH1Ref()
             + prdof[2]*FE_P1CL::DH2Ref() + prdof[3]*FE_P1CL::DH3Ref() );
    for(Uint i=0; i<10; ++i)
    {
        // TODO: 1.0 stands for Stokes:GetCoeff().nu!! How do I obtain this here?
        tmp-= 1.0*veldof[i]*FE_P2CL::Laplace(i, M);
    }
    err_sq+= 4.*hT*hT*inner_prod(tmp, tmp)*vol;


    // Sum over all hF*int( [nu*nE*grad(u) + p*nE]_jumpoverface^2, face) not on the boundary
    const Uint lvl= vel.GetLevel();
    const Uint numpts2= FaceQuad2CL::GetNumPoints();
    vals= new double[numpts2];
    for (Uint f=0; f<NumFacesC; ++f)
    {
        const FaceCL& face= *s.GetFace(f);
        if ( !face.IsOnBoundary() )
        {
            const TetraCL& neigh= *face.GetNeighInTriang(&s, lvl);
            const Uint f_n= face.GetFaceNumInTetra(&neigh);
            // Put velocity degrees of freedom of neigh in veldof_n
            std::vector< SVectorCL<3> > veldof_n( 10);
            vel.GetDoF(neigh, veldof_n);
            SMatrixCL<3,3> M_n;
            double ndet, dir;
            SVectorCL<3> n; // normal of the face; is the same on t and neigh,
                            // if the triangulation is consistently numbered.
            const double absdet2D= s.GetNormal(f, n, dir);
            GetTrafoTr(M_n, ndet, neigh);
            for (Uint pt=0; pt<numpts2; ++pt)
            {
                SVectorCL<3> tmp_me(0.);
                SVectorCL<3> tmp_n(0.);
                const SVectorCL<3> pos= FaceToTetraCoord(s, f, FaceQuad2CL::GetPoints()[pt]);
                const SVectorCL<3> pos_n= FaceToTetraCoord(neigh, f_n, FaceQuad2CL::GetPoints()[pt]);
                for (Uint i=0; i<10; ++i)
                {
                    const double gr= inner_prod( n, M*FE_P2CL::DHRef(i, pos[0], pos[1], pos[2]) );
                    const double ngr= inner_prod( n, M_n*FE_P2CL::DHRef(i, pos_n[0], pos_n[1], pos_n[2]) );
                    // TODO: 1.0 for nu; see above
                    tmp_me+= 1.0*veldof[i]*gr;
                    tmp_n+= 1.0*veldof_n[i]*ngr;
                }
                for (Uint i=0; i<4; ++i)
                {
                    tmp_me+= pr.val(s, pos[0], pos[1], pos[2])*n;
                    tmp_n+= pr.val(neigh, pos_n[0], pos_n[1], pos_n[2])*n;
                }
                vals[pt]= (tmp_me-tmp_n).norm_sq();
            }
            Point3DCL ccF; // Dummy
            double rF; // Radius of the Face
            circumcircle(s, f, ccF, rF);
            err_sq+= 2.*rF*FaceQuad2CL::Quad(vals)*absdet2D;
        }
    }
    delete[] vals;
    return err_sq;
}
#endif // end of ifndef _PAR

template <class Coeff>
void StokesP2P1CL<Coeff>::SetNumVelLvl( size_t n)
{
#ifdef _PAR
    if (n>1)
        throw DROPSErrCL("StokesP2P1CL::SetNumVelLvl: Multilevel not implemented in parallel DROPS, yet, sorry");
#endif
    match_fun match= MG_.GetBnd().GetMatchFun();
    vel_idx.resize( n, vecP2_FE, BndData_.Vel, match);
    A.Data.resize ( vel_idx.size());
    M.Data.resize ( vel_idx.size());
}

template <class Coeff>
void StokesP2P1CL<Coeff>::SetNumPrLvl( size_t n)
{
#ifdef _PAR
    if (n>1)
        throw DROPSErrCL("StokesP2P1CL::SetNumVelLvl: Multilevel not implemented in parallel DROPS, yet, sorry");
#endif
    match_fun match= MG_.GetBnd().GetMatchFun();
    pr_idx.resize( n, P1_FE,  BndData_.Pr, match);
    B.Data.resize( pr_idx.size());
    prM.Data.resize( pr_idx.size());
    prA.Data.resize( pr_idx.size());
}

template <class Coeff>
void StokesP2P1CL<Coeff>::SetIdx()
{
    MLIdxDescCL* vidx= &vel_idx;
    MLIdxDescCL* pidx= &pr_idx;

    b.SetIdx   ( vidx);
    c.SetIdx   ( pidx);

    A.SetIdx   ( vidx, vidx);
    B.SetIdx   ( pidx, vidx);
    prM.SetIdx ( pidx, pidx);
    prA.SetIdx ( pidx, pidx);
    M.SetIdx   ( vidx, vidx);
}

//*********************************************************************
//**************P1Bubble-P1 - discretisation***************************
//*********************************************************************

#ifndef _PAR
//
template <class Coeff>
void StokesP1BubbleP1CL<Coeff>::GetDiscError(instat_vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const
{
    Uint lvl= A.GetRowLevel(),
        vidx= A.RowIdx->GetIdx(),
        pidx= B.RowIdx->GetIdx();
    VectorCL lsgvel(A.RowIdx->NumUnknowns());
    VectorCL lsgpr( B.RowIdx->NumUnknowns());

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!BndData_.Vel.IsOnDirBnd(*sit))
        {
           for(int i=0; i<3; ++i)
               lsgvel[sit->Unknowns(vidx)+i]= LsgVel(sit->GetCoord(), 0.)[i];
        }
    }

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {  // Tetras are never on the boundary.
       for(int i=0; i<3; ++i)
       {
            lsgvel[sit->Unknowns(vidx)+i]= LsgVel( GetBaryCenter(*sit), 0. )[i];
            // The coefficient of the bubble-function is not the value of the solution
            // in the barycenter, as the linear shape-functions contribute to the value
            // there.
            lsgvel[sit->Unknowns(vidx)+i]-= .25*LsgVel( sit->GetVertex(0)->GetCoord(), 0. )[i];
            lsgvel[sit->Unknowns(vidx)+i]-= .25*LsgVel( sit->GetVertex(1)->GetCoord(), 0. )[i];
            lsgvel[sit->Unknowns(vidx)+i]-= .25*LsgVel( sit->GetVertex(2)->GetCoord(), 0. )[i];
            lsgvel[sit->Unknowns(vidx)+i]-= .25*LsgVel( sit->GetVertex(3)->GetCoord(), 0. )[i];
       }
    }
    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
        lsgpr[sit->Unknowns(pidx)]= LsgPr(sit->GetCoord());

    std::cout << "discretization error to check the system (x,y = continuous solution): "<<std::endl;
    VectorCL res( A.Data*lsgvel + transp_mul(B.Data,lsgpr)-b.Data);
    std::cout <<"|| Ax + BTy - f || = "<< norm( res)<<", max "<< supnorm( res) << std::endl;
    VectorCL resB( B.Data*lsgvel - c.Data);
    std::cout <<"|| Bx - g || = "<< norm( resB)<<", max "<< supnorm( resB) << std::endl;
}


/**************************************************************************************************
* member functions to handle with index descriptions
**************************************************************************************************/

//
inline double QuadGradP1Bubble(const SMatrixCL<3,3>& T, Uint i, Uint j)
// integral of the product of the gradients of shape function i and j
// This is exact, for i==j==4 see P1Bubble.txt - Maple is your friend :-)
// You still have to multiply the result of this call with std::fabs(det(T))!
{
    if (i<4 && j<4)
    { // both shape-functions are linear
        return inner_prod( T*FE_P1BubbleCL::DHRef(i), T*FE_P1BubbleCL::DHRef(j) )/6.;
    }
    else if (i==4 && j==4)
    { // the bubble-function in both args
        return 4096./1825.*( T(0,0)*T(0,0) + (T(0,0) + T(0,1))*(T(0,1) + T(0,2)) + T(0,2)*T(0,2)
                            +T(1,0)*T(1,0) + (T(1,0) + T(1,1))*(T(1,1) + T(1,2)) + T(1,2)*T(1,2)
                            +T(2,0)*T(2,0) + (T(2,0) + T(2,1))*(T(2,1) + T(2,2)) + T(2,2)*T(2,2) );
    }
    // linear function with bubble function; the integral of the gradients over the tetra vanishes
    // (proof via partial integration on tetra.)
    return 0.;
}

//
inline double GetB(Uint pr, Uint vel, const SMatrixCL<3,3>& T, Uint num)
{
    if (vel<4)
    { // both shape-functions are linear; evaluation of P1-functions in
      // the barycenter always yields .25;
        const SVectorCL<3> gradient( FE_P1CL::DHRef(vel) );
        return .25/6.*(T(num,0)*gradient[0] + T(num,1)*gradient[1] + T(num,2)*gradient[2]);
    }
    // i==4: we have the bubble-function for the velocity
    return 16./315.*(pr==0 ? T(num,0) + T(num,1) + T(num,2) : -T(num,pr-1) );
}

template <class CoeffT>
void SetupSystem_P1BubbleP1( const MultiGridCL& MG, const CoeffT& Coeff, const StokesBndDataCL& BndData, MatrixCL& matA,
                             VelVecDescCL* vecA, MatrixCL& matB, VelVecDescCL* vecB, IdxDescCL& RowIdxA, IdxDescCL& RowIdxB)
{
// Sets up the stiffness matrices and right hand sides
    if (vecA != 0)
    {
        vecA->Clear( 0.0);
        vecB->Clear( 0.0);
    }

    const IdxT num_unks_vel= RowIdxA.NumUnknowns();
    const IdxT num_unks_pr=  RowIdxB.NumUnknowns();

    MatrixBuilderCL A(&matA, num_unks_vel, num_unks_vel),
                    B(&matB, num_unks_pr,  num_unks_vel);
    VelVecDescCL& b   = *vecA;
    VelVecDescCL& c   = *vecB;
    const Uint lvl    = RowIdxA.TriangLevel();
    const Uint vidx   = RowIdxA.GetIdx(),
               pidx   = RowIdxB.GetIdx();

    IdxT Numb[5], prNumb[4];
    bool IsOnDirBnd[5];

    const IdxT stride= 1;   // stride between unknowns on same simplex, which
                            // depends on numbering of the unknowns

    std::cout << "entering SetupSystem: " <<num_unks_vel<<" vels, "<<num_unks_pr<<" prs"<< std::endl;

    // fill value part of matrices
    SMatrixCL<3,3> T;
    double coup[5][5];
    double det, absdet;
    SVectorCL<3> tmp;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=MG.GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        absdet= std::fabs(det);

        // collect some information about the verts and the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= BndData.Vel.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(vidx);
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);
        }
        Numb[4]= sit->Unknowns(vidx);
        IsOnDirBnd[4]= false; // the bubble-function is never on a boundary

        // compute all couplings between HatFunctions on tetra and verts
        for(int i=0; i<5; ++i)
            for(int j=0; j<=i; ++j)
            {
                // negative dot-product of the gradients
                coup[i][j]= Coeff.nu * QuadGradP1Bubble(T, i, j)*absdet;
//                coup[i][j]+= P1BubbleDiscCL::Quad(*sit, &Coeff::q, i, j)*absdet;
                coup[j][i]= coup[i][j];
            }

        for(int i=0; i<5; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<5; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        A(Numb[i],          Numb[j])+=          coup[j][i];
                        A(Numb[i]+stride,   Numb[j]+stride)+=   coup[j][i];
                        A(Numb[i]+2*stride, Numb[j]+2*stride)+= coup[j][i];
                    }
                    else // coupling with vert/edge j on right-hand-side; j is always <4!
                        if (vecA != 0)
                        {
                            tmp= BndData.Vel.GetDirBndValue(*sit->GetVertex(j));
                            b.Data[Numb[i]]-=          coup[j][i] * tmp[0];
                            b.Data[Numb[i]+stride]-=   coup[j][i] * tmp[1];
                            b.Data[Numb[i]+2*stride]-= coup[j][i] * tmp[2];
                        }
                }
                if (vecA != 0)
                {
                    tmp= P1BubbleDiscCL::Quad(*sit, &CoeffT::f, i)*absdet;
                    b.Data[Numb[i]]+=          tmp[0];
                    b.Data[Numb[i]+stride]+=   tmp[1];
                    b.Data[Numb[i]+2*stride]+= tmp[2];

                    if ( i<4 && BndData.Vel.IsOnNatBnd(*sit->GetVertex(i)) ) // vert i is on natural boundary
                    // i==4 is the bubble-function, which is never on any boundary.
                    {
                        Uint face;
                        for (int f=0; f < 3; ++f)
                        {
                            face= FaceOfVert(i,f);
                            if ( sit->IsBndSeg(face))
                            {
                                tmp= P1DiscCL::Quad2D(*sit, face, BndData.Vel.GetBndSeg(sit->GetBndIdx(face)).GetBndFun(), i);
                                b.Data[Numb[i]]+=          tmp[0];
                                b.Data[Numb[i]+stride]+=   tmp[1];
                                b.Data[Numb[i]+2*stride]+= tmp[2];
                            }
                        }
                    }
                }
            }

        // Setup B:   B(i,j) = int( psi_i * div( phi_j) )
        for(int vel=0; vel<5; ++vel)
        {
            if (!IsOnDirBnd[vel])
                for(int pr=0; pr<4; ++pr)
                {
                    B(prNumb[pr],Numb[vel])+=           GetB(pr, vel, T, 0)*absdet;
                    B(prNumb[pr],Numb[vel]+stride)+=    GetB(pr, vel, T, 1)*absdet;
                    B(prNumb[pr],Numb[vel]+2*stride)+=  GetB(pr, vel, T, 2)*absdet;
                }
            else // put coupling on rhs
                if (vecB != 0)
                {
                    tmp= BndData.Vel.GetDirBndValue( *sit->GetVertex(vel));
                    for(int pr=0; pr<4; ++pr)
                    {
                        // numeric integration is exact: psi_i * div( phi_j) is of degree 2 !
                        c.Data[prNumb[pr]]-= GetB(pr, vel, T, 0)*absdet*tmp[0]
                                            +GetB(pr, vel, T, 1)*absdet*tmp[1]
                                            +GetB(pr, vel, T, 2)*absdet*tmp[2];
                    }
                }
        }
    }
    std::cout << "done: value part fill" << std::endl;

    A.Build();
    B.Build();
    std::cout << matA.num_nonzeros() << " nonzeros in A, "
              << matB.num_nonzeros() << " nonzeros in B! " << std::endl;
}


template <class Coeff>
void StokesP1BubbleP1CL<Coeff>::SetupSystem(MLMatDescCL* matA, VelVecDescCL* vecA, MLMatDescCL* matB, VelVecDescCL* vecB) const
{
    MLMatrixCL::iterator  itA    = matA->Data.begin();
    MLMatrixCL::iterator  itB    = matB->Data.begin();
    MLIdxDescCL::iterator itRowA = matA->RowIdx->begin();
    MLIdxDescCL::iterator itRowB = matB->RowIdx->begin();
    for ( size_t lvl=0; lvl < matA->Data.size(); ++lvl, ++itA, ++itB, ++itRowA, ++itRowB)
    {
        if (lvl != matA->Data.size()-1)
            SetupSystem_P1BubbleP1( MG_, Coeff_, BndData_, *itA, 0, *itB, 0, *itRowA, *itRowB);
        else
            SetupSystem_P1BubbleP1( MG_, Coeff_, BndData_, *itA, vecA, *itB, vecB, *itRowA, *itRowB);
    }
}

template <class CoeffT>
void SetupPrMass_P1BubbleP1( const MultiGridCL& MG, const CoeffT&, MatrixCL& matM, IdxDescCL& RowIdx)
// Sets up the pressure mass matrix
{
    const IdxT num_unks_pr=  RowIdx.NumUnknowns();

    MatrixBuilderCL M(&matM, num_unks_pr,  num_unks_pr);

    const Uint lvl    = RowIdx.TriangLevel();
    const Uint pidx   = RowIdx.GetIdx();

    IdxT prNumb[4];

    // compute all couplings between HatFunctions on verts:
    // M( i, j) = int ( psi_i*psi_j, T_ref) * absdet

    for (MultiGridCL::const_TriangTetraIteratorCL sit=MG.GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double absdet= sit->GetVolume()*6;

        for(int i=0; i<4; ++i)
            prNumb[i]= sit->GetVertex(i)->Unknowns(pidx);

        for(int i=0; i<4; ++i)    // assemble row prNumb[i]
            for(int j=0; j<4; ++j)
                    M( prNumb[i], prNumb[j])+= P1DiscCL::GetMass( i, j) * absdet;
    }
    M.Build();
}

template <class Coeff>
void StokesP1BubbleP1CL<Coeff>::SetupPrMass(MLMatDescCL* matM) const
{
    MLMatrixCL::iterator  itM   = matM->Data.begin();
    MLIdxDescCL::iterator itRow = matM->RowIdx->begin();
    for ( size_t lvl=0; lvl < matM->Data.size(); ++lvl, ++itM, ++itRow)
        SetupPrMass_P1BubbleP1( MG_, Coeff_, *itM, *itRow);
}

template <class Coeff>
void StokesP1BubbleP1CL<Coeff>::SetNumVelLvl( size_t n)
{
    match_fun match= MG_.GetBnd().GetMatchFun();
    vel_idx.resize( vecP1Bubble_FE, BndData_.Vel, match, n);
    A.Data.resize ( vel_idx.size());
}

template <class Coeff>
void StokesP1BubbleP1CL<Coeff>::SetNumPrLvl( size_t n)
{
    match_fun match= MG_.GetBnd().GetMatchFun();
    pr_idx.resize( P1_FE,  BndData_.Pr, match, n);
    B.Data.resize( pr_idx.size());
}

template <class Coeff>
  double
  StokesP1BubbleP1CL<Coeff>::ResidualErrEstimator(const TetraCL& t,
      const const_DiscPrSolCL& pr, const const_DiscVelSolCL& vel, double)
{
    Uint lvl= vel.GetLevel();
    double err_sq= 0.0;
    double det, absdet, vol;
    SMatrixCL<3,3> M;

    GetTrafoTr(M, det, t);
    absdet= std::fabs(det);
    vol= absdet/6.;

    // || div(u_h,l) ||_L2 squared
    for (Uint i=0; i<4; ++i)
        for (Uint j=0; j<4; ++j)
        {
             err_sq+= inner_prod( vel.val(*t.GetVertex(i)), M*FE_P1BubbleCL::DHRef(i) )
                     *inner_prod( vel.val(*t.GetVertex(j)), M*FE_P1BubbleCL::DHRef(j) )/6.*absdet;
        }

    // || P_0(f) - grad(p_h) ||_L2 squared * Vol(T)
    const SVectorCL<3> P0f= Quad3PosWeightsCL::Quad(t, &Coeff::f)*absdet/std::sqrt(vol); // == *absdet/std::sqrt(vol)....
    const SVectorCL<3> gp= pr.val(t.GetVertex(0)->GetCoord())*M*FE_P1CL::DH0Ref()
                          +pr.val(t.GetVertex(1)->GetCoord())*M*FE_P1CL::DH1Ref()
                          +pr.val(t.GetVertex(2)->GetCoord())*M*FE_P1CL::DH2Ref()
                          +pr.val(t.GetVertex(3)->GetCoord())*M*FE_P1CL::DH3Ref();
    err_sq+= (P0f - gp).norm_sq()*vol*vol;

// TODO: non-null dirichlet-boundary values, natural boundary
    // .5*sum_{E \in \partial T \cap \Omega} |E|*|| Sprung von u_h,l ||_L2,E^2
    SMatrixCL<3,4> vertval;
    for (Uint i=0; i<NumVertsC; ++i)
    {
        const SVectorCL<3> v= vel.val(*t.GetVertex(i));
        vertval(0,i)= v[0];
        vertval(1,i)= v[1];
        vertval(2,i)= v[2];
    }
    for (Uint f=0; f<NumFacesC; ++f)
    {
        const FaceCL& face= *t.GetFace(f);
        if ( !face.IsOnBoundary() )
        {
            const TetraCL& neigh= *face.GetNeighInTriang(&t, lvl);
            SMatrixCL<3,3> nM;
            double ndet, dir;
            SVectorCL<3> n;
            t.GetNormal(f, n, dir);
            GetTrafoTr(nM, ndet, neigh);
            const double absdet2D= FuncDet2D( face.GetVertex(1)->GetCoord() - face.GetVertex(0)->GetCoord(),
                                              face.GetVertex(2)->GetCoord() - face.GetVertex(0)->GetCoord() );
            SMatrixCL<3,4> nvertval;
            for (Uint i=0; i<NumVertsC; ++i)
            {
                const SVectorCL<3> v= vel.val(*neigh.GetVertex(i));
                nvertval(0,i)= v[0];
                nvertval(1,i)= v[1];
                nvertval(2,i)= v[2];
            }
            SVectorCL<3> grad_me(0.0);
            SVectorCL<3> grad_n(0.0);
            for (Uint i=0; i<NumVertsC; ++i)
            {
                const double gr= inner_prod( n, M*FE_P1CL::DHRef(i) );
                const double ngr= inner_prod( n, nM*FE_P1CL::DHRef(i) );
                grad_me[0]+= vertval(0,i)*gr;
                grad_me[1]+= vertval(1,i)*gr;
                grad_me[2]+= vertval(2,i)*gr;
                grad_n[0]+= nvertval(0,i)*ngr;
                grad_n[1]+= nvertval(1,i)*ngr;
                grad_n[2]+= nvertval(2,i)*ngr;
            }
            err_sq+= (grad_n - grad_me).norm_sq()/8.*absdet2D*absdet2D;
        }
    }
//    std::cout << err_sq << '\n';
    return err_sq;
}

template <class _TetraEst, class _ProblemCL>
void StokesDoerflerMarkCL<_TetraEst, _ProblemCL>::Init(const const_DiscPrSolCL& pr, const const_DiscVelSolCL& vel)
{
    MultiGridCL& mg= _Problem.GetMG();
    const Uint lvl= vel.GetLevel();

    double tmp;
    _InitGlobErr= 0.;
    for (MultiGridCL::TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        tmp= _Estimator(*sit, pr, vel, 0.0);
        _InitGlobErr+= tmp;
    }
    _InitGlobErr= std::sqrt(_InitGlobErr);
    _ActGlobErr= _InitGlobErr;
    if (_outp)
        *_outp << "ErrorEstimator is initialized now; initial error: " << _InitGlobErr
               << " We try to reduce by a factor of " << _RelReduction
               << " by marking the tetrahedrons with largest error-estimates, until the error marked"
               << " is at least " << _Threshold << " of the actual global error." << std::endl;
}

template <class _TetraEst, class _ProblemCL>
bool StokesDoerflerMarkCL<_TetraEst, _ProblemCL>::Estimate(const const_DiscPrSolCL& pr, const const_DiscVelSolCL& vel)
{
    Err_ContCL err_est;
    const MultiGridCL& mg= _Problem.GetMG();
    const Uint lvl= vel.GetLevel();
//    const VecDescCL& lsg= *vel.GetSolution();

    _NumLastMarkedForRef= 0;
    _NumLastMarkedForDel= 0;
    for (MultiGridCL::const_TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double localerr= _Estimator(*sit, pr, vel, 0.0);
        err_est.push_back( std::make_pair(&*sit, localerr) );
    }
    const double globalerr_sq= std::accumulate(err_est.begin(), err_est.end(), 0.0, AccErrCL() );
    const double globalerr= std::sqrt(globalerr_sq);
    const double ref_threshold_sq= globalerr_sq*_Threshold*_Threshold;
    if (globalerr>=_InitGlobErr*_RelReduction && _DoMark)
    {
        std::sort( err_est.begin(), err_est.end(), Err_Pair_GTCL() );
        double akt_ref_err_sq= 0;
        for (Err_ContCL::iterator it= err_est.begin(), theend= err_est.end();
             it != theend && akt_ref_err_sq < ref_threshold_sq; ++it)
            if ( it->first->IsUnrefined() )
            {
                it->first->SetRegRefMark();
                akt_ref_err_sq+= it->second;
                ++_NumLastMarkedForRef;
            }
    }
    if (_outp)
        *_outp << "Estimated global |v|1 + ||pr||0-error: " << globalerr << ". This is a reduction of " << globalerr/_ActGlobErr
               << " compared to the last estimation. Marked " << _NumLastMarkedForRef
               << ", unmarked " << _NumLastMarkedForDel << " out of " << err_est.size() << " tetrahedrons, "
               << "which account for " << _Threshold << " of the global error."
               << std::endl;
    _ActGlobErr= globalerr;
    return _NumLastMarkedForRef || _NumLastMarkedForDel;
}


template <class Coeff>
void StokesP1BubbleP1CL<Coeff>::CheckSolution(const VelVecDescCL* lsgvel, const VecDescCL* lsgpr,
    instat_vector_fun_ptr LsgVel, scalar_fun_ptr LsgPr) const
{
    double diff, maxdiff=0, norm2= 0;
    Uint lvl=lsgvel->GetLevel(),
         vidx=lsgvel->RowIdx->GetIdx();

    VectorCL res1( A.Data*lsgvel->Data + transp_mul( B.Data, lsgpr->Data ) - b.Data);
    VectorCL res2( B.Data*lsgvel->Data - c.Data);

    std::cout << "\nChecken der Loesung...\n";
    std::cout << "|| Ax + BTy - F || = " << norm( res1) << ", max. " << supnorm( res1) << std::endl;
    std::cout << "||       Bx - G || = " << norm( res2) << ", max. " << supnorm( res2) << std::endl << std::endl;

    const_DiscVelSolCL vel(lsgvel, &BndData_.Vel, &MG_);
    double L1_div= 0, L2_div= 0;
    SMatrixCL<3,3> T;
    double det, absdet;

// We only use the linear part of the velocity-solution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double* pvals= new double[Quad3PosWeightsCL::GetNumPoints()];
    double* pvals_sq= new double[Quad3PosWeightsCL::GetNumPoints()];
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        absdet= std::fabs(det);

        for (Uint i=0; i<Quad3PosWeightsCL::GetNumPoints(); ++i)
        {
            pvals[i]= std::fabs( inner_prod(T*FE_P1BubbleCL::DH0Ref(), vel.val(*sit->GetVertex(0)) )
                           +inner_prod(T*FE_P1BubbleCL::DH1Ref(), vel.val(*sit->GetVertex(1)) )
                           +inner_prod(T*FE_P1BubbleCL::DH2Ref(), vel.val(*sit->GetVertex(2)) )
                           +inner_prod(T*FE_P1BubbleCL::DH3Ref(), vel.val(*sit->GetVertex(3)) ) );
            pvals_sq[i]= pvals[i]*pvals[i];
        }
        L1_div+= Quad3PosWeightsCL::Quad(pvals)*absdet;
        L2_div+= Quad3PosWeightsCL::Quad(pvals_sq)*absdet;
    }
    L2_div= std::sqrt(L2_div);
    std::cout << "|| div x ||_L1 = " << L1_div << std::endl;
    std::cout << "|| div x ||_L2 = " << L2_div << std::endl << std::endl;
    delete[] pvals_sq; delete[] pvals;

    Uint countverts=0;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (!BndData_.Vel.IsOnDirBnd(*sit))
        {
           ++countverts;
           for(int i=0; i<3; ++i)
           {
               diff= std::fabs( LsgVel(sit->GetCoord(), 0.)[i] - lsgvel->Data[sit->Unknowns(vidx)+i] );
               norm2+= diff*diff;
               if (diff>maxdiff)
               {
                   maxdiff= diff;
               }
           }
        }
    }
    norm2= std::sqrt(norm2/countverts);

    Point3DCL L1_vel(0.0), L2_vel(0.0);
    SVectorCL<3>* vvals= new SVectorCL<3>[Quad3PosWeightsCL::GetNumPoints()];
    SVectorCL<3>* vvals_sq= new SVectorCL<3>[Quad3PosWeightsCL::GetNumPoints()];
    for(MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send= const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
        sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        absdet= std::fabs(det);
//        Point3DCL sum(0.0), diff, Diff[5];
        for (Uint i=0; i<Quad3PosWeightsCL::GetNumPoints(); ++i)
        {
            const Point3DCL& pt= Quad3PosWeightsCL::GetPoints()[i];
            vvals[i]= fabs(LsgVel(GetWorldCoord(*sit, pt), 0.) - vel.lin_val(*sit, pt[0], pt[1], pt[2]));
            vvals_sq[i]= vvals[i]*vvals[i];
        }
        L1_vel+= Quad3PosWeightsCL::Quad(vvals)*absdet;
        L2_vel+= Quad3PosWeightsCL::Quad(vvals_sq)*absdet;
    }
    L2_vel= sqrt(L2_vel);
    delete[] vvals_sq; delete[] vvals;
    std::cout << "Geschwindigkeit: Abweichung von der tatsaechlichen Loesung:\n"
              << "w-2-Norm= " << norm2 << std::endl
              << " L2-Norm= (" << L2_vel[0]<<", "<<L2_vel[1]<<", "<<L2_vel[2]<<")" << std::endl
              << " L1-Norm= (" << L1_vel[0]<<", "<<L1_vel[1]<<", "<<L1_vel[2]<<")" << std::endl
              << "max-Norm= " << maxdiff << std::endl;

    norm2= 0; maxdiff= 0; double mindiff= 1000;

    // Compute the pressure-coefficient in direction of 1/std::sqrt(meas(Omega)), which eliminates
    // the allowed offset of the pressure by setting it to 0.
    double L1_pr= 0, L2_pr= 0, MW_pr= 0, vol= 0;
    P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr(lsgpr, &BndData_.Pr, &MG_);
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
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
    std::cout << "\nconstant pressure offset is " << c_pr<<", volume of cube is " << vol<<std::endl;;

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        diff= std::fabs( c_pr + LsgPr(sit->GetCoord()) - pr.val(*sit));
        norm2+= diff*diff;
        if (diff>maxdiff)
            maxdiff= diff;
        if (diff<mindiff)
            mindiff= diff;
    }
    norm2= std::sqrt(norm2 / lsgpr->Data.size() );

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double sum= 0, sum1= 0;
        for(int i=0; i<4; ++i)
        {
            diff= c_pr + LsgPr(sit->GetVertex(i)->GetCoord()) - pr.val(*sit->GetVertex(i));
            sum+= diff*diff; sum1+= std::fabs(diff);
        }
        sum/= 120;   sum1/= 120;
        diff= c_pr + LsgPr(GetBaryCenter(*sit)) - pr.val(*sit, .25, .25, .25);
        sum+= 2./15.*diff*diff;   sum1+= 2./15.*std::fabs(diff);
        L2_pr+= sum * sit->GetVolume()*6.;
        L1_pr+= sum1 * sit->GetVolume()*6.;
    }
    L2_pr= std::sqrt( L2_pr);


    std::cout << "Druck: Abweichung von der tatsaechlichen Loesung:\n"
              << "w-2-Norm= " << norm2 << std::endl
              << " L2-Norm= " << L2_pr << std::endl
              << " L1-Norm= " << L1_pr << std::endl
              << "Differenz liegt zwischen " << mindiff << " und " << maxdiff << std::endl;
}
#endif // end of ifndef _PAR
} // end of namespace DROPS
