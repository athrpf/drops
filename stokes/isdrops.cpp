/// \file isdrops.cpp
/// \brief instationary stokes problem
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "geom/multigrid.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "geom/builder.h"
#include "num/stokessolver.h"
#include "stokes/stokes.h"
#include "stokes/integrTime.h"
#include <fstream>
#include <sstream>


struct InstatStokesCL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])*t*t;
        ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])*t*t;
        ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*t;
        return ret/3.;
    }

    static double LsgPr(const DROPS::Point3DCL& p, double t)
    {
        return std::cos(p[0])*std::sin(p[1])*std::sin(p[2])*t*t
               -(std::sin( 1.) -2.*std::sin( 1.)*std::cos( 1.) + std::sin( 1.)*std::pow( std::cos( 1.), 2))*t*t; // (...)==0.1778213062
    }

    // du/dt + q*u - nu*laplace u + Dp = f
    //                          -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double t)
        {
            DROPS::SVectorCL<3> ret;
            ret[0]= 2./3.*t*std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
            ret[1]= -2./3.*t*std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
            ret[2]= std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*(4./3. + 3.*t);
            return ret;
        }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };
    static StokesCoeffCL Coeff;
};

InstatStokesCL::StokesCoeffCL InstatStokesCL::Coeff;


namespace DROPS
{
template <class StokesT, class SolverT>
class MyInstatStokesThetaSchemeCL
//**************************************************************************
//  for solving the instationary Stokes equation of type StokesT with a
//  1-step-theta-scheme: theta=1   -> impl. Euler (BDF1, backward Euler)
//                       theta=1/2 -> Crank-Nicholson (Trapezregel)
//
//  Inner stationary Stokes-type problems are solved with a SolverT-solver.
//  The matrices A, B, M and the rhs b, c of the Stokes class have to be set
//  properly! After construction, SetTimeStep has to be called once. Then
//  every DoStep performs one step in time. Changing time steps require
//  further calls to SetTimeStep.
//**************************************************************************
{
  private:
    StokesT& _Stokes;
    SolverT& _solver;

    VelVecDescCL *_b, *_old_b;        // rhs + couplings with poisson matrix A
    VelVecDescCL *_cplM, *_old_cplM;  // couplings with mass matrix M
    VectorCL      _rhs;
    MLMatrixCL    _mat;               // (1./dt)*M + theta*A

    double _theta, _dt;

  public:
    MyInstatStokesThetaSchemeCL( StokesT& Stokes, SolverT& solver, double theta= 0.5)
        : _Stokes( Stokes), _solver( solver), _b( &Stokes.b), _old_b( new VelVecDescCL),
          _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), _rhs( Stokes.b.RowIdx->NumUnknowns),
          _theta( theta)
    {
        _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
        _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    }

    ~MyInstatStokesThetaSchemeCL()
    {
        if (_old_b == &_Stokes.b)
            delete _b;
        else
            delete _old_b;
        delete _cplM; delete _old_cplM;
    }

    double GetTheta()    const { return _theta; }
    double GetTime()     const { return _Stokes.t; }
    double GetTimeStep() const { return _dt; }

    void SetTimeStep( double dt)
    {
        _dt= dt;
        _mat.LinComb( 1./dt, _Stokes.M.Data, _theta, _Stokes.A.Data);
    }

    void DoStep( VectorCL& v, VectorCL& p);
};

template <class StokesT, class SolverT>
void MyInstatStokesThetaSchemeCL<StokesT,SolverT>::DoStep( VectorCL& v, VectorCL& p)
{
    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * v;
    _rhs*= (_theta-1.);
    _rhs+= (1./_dt)*(_Stokes.M.Data*v + _cplM->Data - _old_cplM->Data)
         +  _theta*_b->Data + (1.-_theta)*_old_b->Data;

    _solver.Solve( _mat, _Stokes.B.Data, v, p, _rhs, _Stokes.c.Data);

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
}


template <typename PressurePreT, typename PoissonSolver2T>
class MyUzawaSolver2CL : public SolverBaseCL
{
  private:
    PressurePreT&    pr_pre_;
    PoissonSolver2T& poissonSolver2_;
    MatrixCL& M_;
    double    tau_;

  public:
    MyUzawaSolver2CL (PressurePreT& pre, PoissonSolver2T& solver2,
                      MatrixCL& M, int maxiter, double tol, double tau= 1.)
        : SolverBaseCL( maxiter, tol), pr_pre_( pre), poissonSolver2_( solver2),
          M_( M), tau_( tau) {}

    double GetTau()            const { return tau_; }
    void   SetTau( double tau)       { tau_= tau; }

    void Solve( const MLMatrixCL& A, const MLMatrixCL& B, VectorCL& v, VectorCL& p,
                const VectorCL& b, const VectorCL& c);
};

template <class PoissonSolverT, class PoissonSolver2T>
void MyUzawaSolver2CL<PoissonSolverT, PoissonSolver2T>::Solve(
    const MLMatrixCL& A, const MLMatrixCL& B,
    VectorCL& v, VectorCL& p, const VectorCL& b, const VectorCL& c)
{
    VectorCL v_corr( v.size()),
             p_corr( p.size()),
             res1( v.size()),
             res2( p.size());
    double tol= _tol;
    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( _iter=0; _iter<_maxiter; ++_iter) {
        z_xpay(res2, B*v, -1.0, c);
        res2_norm= norm_sq( res2);
        pr_pre_.Apply( M_, p_corr, res2);
//        p+= _tau * p_corr;
        axpy(tau_, p_corr, p);
//        res1= A*v + transp_mul(B,p) - b;
        z_xpaypby2( res1, A*v, 1.0, transp_mul( B, p), -1.0, b);
        res1_norm= norm_sq( res1);
        if (res1_norm + res2_norm < tol) {
            _res= std::sqrt( res1_norm + res2_norm);
            return;
        }
        if( (_iter%output)==0)
            std::cout << "step " << _iter << ": norm of 1st eq= " << std::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << std::sqrt( res2_norm) << std::endl;

        poissonSolver2_.Apply( A, v_corr, res1);
//        poissonSolver2_.SetTol( std::sqrt( res1_norm)/20.0);
//        poissonSolver2_.Solve( A, v_corr, res1);
//        std::cout << "velocity: iterations: " << poissonSolver2_.GetIter()
//                  << "\tresidual: " << poissonSolver2_.GetResid() << std::endl;
        v-= v_corr;
    }
    _res= std::sqrt( res1_norm + res2_norm );
}

void
ZeroMean(DROPS::P1EvalCL< double,
                          const DROPS::StokesPrBndDataCL,
                          DROPS::VecDescCL>& f)
{
    const DROPS::Uint lvl= f.GetLevel();
    DROPS::MultiGridCL& mg= const_cast<DROPS::MultiGridCL&>( f.GetMG());
    double MV= 0., vol= 0., sum;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin( lvl),
         send= mg.GetTriangTetraEnd( lvl); sit != send; ++sit) {
        sum= 0.;
        for(int i=0; i<4; ++i)
            sum+= f.val( *sit->GetVertex( i));
        sum/= 120;
        sum+= 2./15.*f.val( *sit, .25, .25, .25);
        MV+= sum * sit->GetVolume()*6.;
        vol+= sit->GetVolume();
    }
    const double c= MV/vol;
    std::cout << "\nconstant pressure offset: " << c << ", volume of domain: " << vol
              << std::endl;
    for (DROPS::MultiGridCL::TriangVertexIteratorCL sit= mg.GetTriangVertexBegin( lvl),
         send= mg.GetTriangVertexEnd( lvl); sit != send; ++sit) {
        f.SetDoF( *sit, f.val( *sit) - c);
    }
}

typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > APcT;
typedef BlockPreCL<APcT, ISMGPreCL> PcT;
typedef PMResSolverCL<PLanczosONBCL<DROPS::VectorCL,PcT> > SolverT;
class PMinresSP_FullMG_CL : public BlockMatrixSolverCL<SolverT>
{
  private:
    SolverT solver_;
    PcT pre_;
    SSORsmoothCL smoother_;
    SSORPcCL ssor_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc;
    APcT Apc_;
    ISMGPreCL Spc_;
    PLanczosONBCL<DROPS::VectorCL, PcT> q_;

  public:
    PMinresSP_FullMG_CL( DROPS::MLMatrixCL& MGApr,
                         DROPS::MLMatrixCL& Mpr, double kA, double kM,
                         int iter_vel, int iter_prA, int iter_prM, int maxiter, double tol)
        :BlockMatrixSolverCL<SolverT> (solver_), solver_(q_, maxiter, tol),
         pre_( Apc_, Spc_), smoother_(1.0), coarsesolver_ ( ssor_, 500, 1e-14),
         mgc( smoother_, coarsesolver_, iter_vel, 1e-14, false), Apc_( mgc),
         Spc_( MGApr, Mpr, kA, kM, iter_prA, iter_prM), q_( pre_)
    {}
    MLMatrixCL* GetPVel() { return mgc.GetProlongation(); }
    MLMatrixCL* GetPPr()  { return Spc_.GetProlongation(); }
};


} // end of namespace DROPS


// Assumes, that indices for A_pr are set up. We know, there are only natural
// boundary conditions.
void
SetupPoissonPressure( DROPS::MultiGridCL& mg, DROPS::MatrixCL& A_pr, DROPS::IdxDescCL& RowIdx, DROPS::IdxDescCL& ColIdx)
{
    DROPS::MatrixBuilderCL A( &A_pr, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    const DROPS::Uint lvl= RowIdx.TriangLevel();
    const DROPS::Uint idx= RowIdx.GetIdx();
    DROPS::SMatrixCL<3,4> G;
    double coup[4][4];
    double det;
    double absdet;
    DROPS::IdxT UnknownIdx[4];

    for (DROPS::MultiGridCL::const_TriangTetraIteratorCL sit= const_cast<const DROPS::MultiGridCL&>( mg).GetTriangTetraBegin( lvl),
         send= const_cast<const DROPS::MultiGridCL&>( mg).GetTriangTetraEnd( lvl);
         sit != send; ++sit) {
        DROPS::P1DiscCL::GetGradients( G,det,*sit);
        absdet= std::fabs( det);
        for(int i=0; i<4; ++i) {
            for(int j=0; j<=i; ++j) {
                // dot-product of the gradients
                coup[i][j]= ( G( 0, i)*G( 0, j) + G( 1, i)*G( 1, j) + G( 2, i)*G( 2, j) )/6.0*absdet;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetVertex( i)->Unknowns( idx);
        }
        for(int i=0; i<4; ++i)    // assemble row i
            for(int j=0; j<4;++j)
                A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i];
    }
    A.Build();
    std::cout << A_pr.num_nonzeros() << " nonzeros in A_pr.\n";
}

// We know, there are only natural boundary conditions.
template<class Coeff>
void
SetupPoissonPressureMG(DROPS::StokesP2P1CL<Coeff>& stokes, DROPS::MLMatDescCL& MGData)
{
    DROPS::MultiGridCL& mg= stokes.GetMG();
    DROPS::MLIdxDescCL::iterator itRow= MGData.RowIdx->begin();
    DROPS::MLIdxDescCL::iterator itCol= MGData.ColIdx->begin();
    for (DROPS::MLMatrixCL::iterator itA = MGData.Data.begin(); itA != MGData.Data.end(); ++itA)
    {
        SetupPoissonPressure( mg, *itA, *itRow, *itCol);
        ++itRow; ++itCol;
    }
    std::cout << "Check MG-Data..." << std::endl;
    std::cout << "                begin     " << MGData.RowIdx->GetCoarsest().NumUnknowns() << std::endl;
    std::cout << "                end       " << MGData.RowIdx->GetFinest().NumUnknowns() << std::endl;
    //CheckMGData( MGData.begin(), MGData.end());
}

namespace DROPS
{

class PSchur2_PCG_Pr_CL: public PSchurSolver2CL<PCG_SsorCL,
                                                PCGSolverCL<ISPreCL> >
{
  private:
    SSORPcCL             ssor_;
    PCG_SsorCL           PCGsolver_;
    PCGSolverCL<ISPreCL> PCGsolver2_;

  public:
    PSchur2_PCG_Pr_CL(ISPreCL& Spc, int outer_iter, double outer_tol,
                                    int inner_iter, double inner_tol)
        : PSchurSolver2CL<PCG_SsorCL, PCGSolverCL<ISPreCL> >(
              PCGsolver_, PCGsolver2_, outer_iter, outer_tol
          ),
          PCGsolver_( ssor_, inner_iter, inner_tol),
          PCGsolver2_( Spc, outer_iter, outer_tol)
        {}
};

class PSchur2_PCG_Pr_MG_CL: public PSchurSolver2CL<PCG_SsorCL,
                                                   PCGSolverCL<ISMGPreCL> >
{
  private:
    SSORPcCL               ssor_;
    PCG_SsorCL             PCGsolver_;
    PCGSolverCL<ISMGPreCL> PCGsolver2_;

  public:
    PSchur2_PCG_Pr_MG_CL(ISMGPreCL& Spc, int outer_iter, double outer_tol,
                                         int inner_iter, double inner_tol)
        : PSchurSolver2CL<PCG_SsorCL, PCGSolverCL<ISMGPreCL> >(
              PCGsolver_, PCGsolver2_, outer_iter, outer_tol
          ),
          PCGsolver_( ssor_, inner_iter, inner_tol),
          PCGsolver2_( Spc, outer_iter, outer_tol)
        {}
};

class PSchur2_Full_MG_CL: public PSchurSolver2CL<MGSolverCL<SSORsmoothCL, PCG_SsorCL>,
                                                 PCGSolverCL<ISMGPreCL> >
{
  private:
    SSORsmoothCL smoother_;
    SSORPcCL     ssor_;
    PCG_SsorCL   coarsesolver_;
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> solver_;
    PCGSolverCL<ISMGPreCL> solver2_;

  public:
    PSchur2_Full_MG_CL( ISMGPreCL& Spc,
                        int outer_iter, double outer_tol,
                        int inner_iter, double inner_tol)
        : PSchurSolver2CL<MGSolverCL<SSORsmoothCL, PCG_SsorCL>, PCGSolverCL<ISMGPreCL> >(
              solver_, solver2_, outer_iter, outer_tol),
          smoother_(1.0), coarsesolver_(ssor_, 500, inner_tol),
          solver_( smoother_, coarsesolver_, inner_iter, inner_tol),
          solver2_( Spc, outer_iter, outer_tol)
        {}
    MLMatrixCL* GetPVel() { return solver_.GetProlongation(); }
};

} // end of namespace DROPS


typedef InstatStokesCL MyPdeCL;

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&, double);

int
CheckVel(DROPS::P2EvalCL< DROPS::SVectorCL<3>,
                          const DROPS::StokesVelBndDataCL,
                          DROPS::VelVecDescCL>& fun,
         fun_ptr f)
{
    using namespace DROPS;
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const double t= fun.GetTime();
    const DROPS::Uint trilevel= fun.GetLevel();
    std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0.;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= (fun.val( *sit) - f( sit->GetCoord(), t)).norm();
        if ( std::abs( diff) > vmaxdiff) { ++ret; vmaxdiff= std::abs( diff); v= &*sit; }
    }
    std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin( trilevel),
         theend= mg.GetTriangEdgeEnd( trilevel); sit!=theend; ++sit) {
        diff = (fun.val( *sit, .5) - f( (sit->GetVertex( 0)->GetCoord()
                                        +sit->GetVertex( 1)->GetCoord())*0.5, t)).norm();
        if (std::abs( diff) > emaxdiff) { ++ret; emaxdiff= std::abs( diff); e= &*sit; }
    }
    {
        std::cout << "maximale Differenz Vertices: " << vmaxdiff << " auf\n";
        if (v) v->DebugInfo( std::cout);
        std::cout << "maximale Differenz Edges: " << emaxdiff << " auf\n";
        if (e) e->DebugInfo( std::cout);
        std::cout << std::endl;
    }
    return ret;
}


const double radiusorbit= 0.3; // Radius of the drops' orbit.
const double radiusdrop= 0.15; // Initial radius of the drop.

// positive outside the drop, negative inside the drop.
double
SignedDistToInterface(const DROPS::Point3DCL& p, double t)
{
   DROPS::Point3DCL c;
   c[0]= 0.5 + radiusorbit*std::cos( 2.*M_PI*t);
   c[1]= 0.5 + radiusorbit*std::sin( 2.*M_PI*t);
   c[2]= 0.5;
   return (p-c).norm() - radiusdrop;
}

typedef double (*signed_dist_fun)(const DROPS::Point3DCL& p, double t);


bool
ModifyGridStep(DROPS::MultiGridCL& mg,
               const signed_dist_fun Dist,
               const double width,         // Thickness of refined shell on each side of the interface
               const DROPS::Uint c_level,  // Outside the shell, use this level
               const DROPS::Uint f_level,  // Inside the shell, use this level
               const double t)             // Time of evaluation
// One step of grid change; returns true if modifications were necessary,
// false, if nothing changed.
{
    using namespace DROPS;
    bool shell_not_ready= false;
        for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
             theend= mg.GetTriangTetraEnd(); it!=theend; ++it) {
            double d= 1.;
            for (Uint j=0; j<4; ++j) {
                d= std::min( d, std::abs( Dist( it->GetVertex( j)->GetCoord(), t)));
            }
            const Uint l= it->GetLevel();
            if (d<=width) { // In the shell; level should be f_level.
                if (l < f_level) {
                    shell_not_ready= true;
                    it->SetRegRefMark();
                }
                else
                    if (l > f_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
            else { // Outside the shell; level should be c_level;
                if (l < c_level) { it->SetRegRefMark(); }
                else
                    if (l> c_level) { it->SetRemoveMark(); }
                    else {} // nothing
            }
        }
        mg.Refine();
    return shell_not_ready;
}

template<class Coeff>
void
UpdateTriangulation(DROPS::StokesP2P1CL<Coeff>& NS,
                    const signed_dist_fun Dist,
                    const double t,
                    const double width,         // Thickness of refined shell on eache side of the interface
                    const DROPS::Uint c_level,  // Outside the shell, use this level
                    const DROPS::Uint f_level,  // Inside the shell, use this level
                    DROPS::VelVecDescCL* v1,
                    DROPS::VecDescCL* p1)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;
    Assert( c_level<=f_level, "UpdateTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    MultiGridCL& mg= NS.GetMG();
    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= v1->RowIdx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx1= p1->RowIdx;
    IdxDescCL* pidx2= &loc_pidx;
    VelVecDescCL  loc_v;
    VecDescCL     loc_p;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p2= &loc_p;
    vidx2->SetFE( vecP2_FE);
    pidx2->SetFE( P1_FE);
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    const StokesBndDataCL& BndData= NS.GetBndData();
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i) {
        shell_not_ready= ModifyGridStep( mg, Dist, width, c_level, f_level, t);
        // Repair velocity
        std::swap( v2, v1);
        std::swap( vidx2, vidx1);
        match_fun match= mg.GetBnd().GetMatchFun();
        vidx1->CreateNumbering( mg.GetLastLevel(), mg, NS.GetBndData().Vel, match);
        if ( mg.GetLastLevel() != vidx2->TriangLevel()) {
            std::cout << "LastLevel: " << mg.GetLastLevel()
                      << " vidx2->TriangLevel: " << vidx2->TriangLevel() << std::endl;
            throw DROPSErrCL( "Strategy: Sorry, not yet implemented.");
        }
        v1->SetIdx( vidx1);
        P2EvalCL< SVectorCL<3>, const StokesVelBndDataCL,
                  const VelVecDescCL> funv2( v2, &BndData.Vel, &mg, t);
        RepairAfterRefineP2( funv2, *v1);
        v2->Clear();
        vidx2->DeleteNumbering( mg);
//P2EvalCL< SVectorCL<3>, const StokesVelBndDataCL,
//          VelVecDescCL> funv1( v1, &BndData.Vel, &mg, t);
//CheckVel( funv1, &MyPdeCL::LsgVel);
        // Repair pressure
        std::swap( p2, p1);
        std::swap( pidx2, pidx1);
        pidx1->CreateNumbering( mg.GetLastLevel(), mg, NS.GetBndData().Pr, match);
        p1->SetIdx( pidx1);
        typename StokesCL::const_DiscPrSolCL oldfunpr( p2, &BndData.Pr, &mg);
        RepairAfterRefineP1( oldfunpr, *p1);
        p2->Clear();
        pidx2->DeleteNumbering( mg);
    }
    // We want the solution to be where v1, p1 point to.
    if (v1 == &loc_v) {
        NS.vel_idx.GetFinest().swap( loc_vidx);
        NS.pr_idx.GetFinest().swap( loc_pidx);
        NS.v.SetIdx( &NS.vel_idx);
        NS.p.SetIdx( &NS.pr_idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
    }
    time.Stop();
    std::cout << "UpdateTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

void
MakeInitialTriangulation(DROPS::MultiGridCL& mg,
                         const signed_dist_fun Dist,
                         const double width,         // Thickness of refined shell on eache side of the interface
                         const DROPS::Uint c_level,  // Outside the shell, use this level
                         const DROPS::Uint f_level)  // Inside the shell, use this level
{
    using namespace DROPS;
    Assert( c_level<=f_level, "MakeInitialTriangulation: Levels are cheesy.\n", ~0);
    TimerCL time;

    time.Reset();
    time.Start();
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i)
        shell_not_ready=  ModifyGridStep( mg, Dist, width, c_level, f_level, 0.);
    time.Stop();
    std::cout << "MakeTriangulation: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg.GetLastLevel() << '\n';
    mg.SizeInfo( std::cout);
}

template<class Coeff>
void
SetMatVecIndices(DROPS::StokesP2P1CL<Coeff>& NS,
                 DROPS::MLIdxDescCL* const vidx,
                 DROPS::MLIdxDescCL* const pidx)
{
    std::cout << "#Druck-Unbekannte: " << pidx->NumUnknowns() << std::endl;
    std::cout << "#Geschwindigkeitsunbekannte: " << vidx->NumUnknowns ()<< std::endl;
    NS.b.SetIdx( vidx);
    NS.c.SetIdx( pidx);
    NS.A.SetIdx( vidx, vidx);
    NS.B.SetIdx( pidx, vidx);
    NS.M.SetIdx( vidx, vidx);
}

template<class Coeff>
void
ResetSystem(DROPS::StokesP2P1CL<Coeff>& NS)
{
    NS.A.Reset(); NS.B.Reset();
    NS.M.Reset();
    NS.b.Reset(); NS.c.Reset();
}


template<class Coeff>
void
StrategyMRes(DROPS::StokesP2P1CL<Coeff>& NS,
             int stokes_maxiter, double stokes_tol,
             double theta,
             DROPS::Uint num_timestep,
             double kA, double kM,
             double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
    MLMatDescCL  ML_pr, MG_pr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;

//    typedef PMinresSP_Id_CL StatsolverCL;
    typedef PMinresSP_FullMG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    typedef InstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
//    typedef MyInstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
    InstatsolverCL* instatsolver= 0;

    MakeInitialTriangulation( mg, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.SetNumVelLvl( mg.GetNumLevel());
    NS.SetNumPrLvl ( mg.GetNumLevel());
    M_pr.Data.resize( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    // Initialize Ensight6 output
    std::string ensf( "insa");
    Ensight6OutCL ensight( ensf + ".case", num_timestep + 1);
    ensight.Register( make_Ensight6Geom      ( mg, mg.GetLastLevel(),   "insa geometry", ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( NS.GetPrSolution(),  "p",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( NS.GetVelSolution(), "v",      ensf + ".vel", true));
    ensight.Write( 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt) {
        std::cout << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
//        if (timestep%(num_timestep/10) == 0) { // modify the grid
        if (timestep == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
//                M_pr.Reset();
                ResetSystem( NS);
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            MG_pr.SetIdx( pidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cout << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);

            SetupPoissonPressureMG( NS, MG_pr);

//            statsolver= new StatsolverCL( stokes_maxiter, stokes_tol);
            statsolver= new PMinresSP_FullMG_CL( MG_pr.Data, M_pr.Data, kA, kM,
                                                 1, 1, 1, stokes_maxiter, stokes_tol);
            SetupP2ProlongationMatrix( mg, *(statsolver->GetPVel()), vidx1, vidx1);
            SetupP1ProlongationMatrix( mg, *(statsolver->GetPPr()), pidx1, pidx1);
            instatsolver= new InstatsolverCL( NS, *statsolver, theta);
        }
        instatsolver->SetTimeStep( dt);
        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,
             DROPS::VecDescCL> pr( p1, &NS.GetBndData().Pr, &mg);
        ZeroMean( pr);
        std::cout << "Before timestep." << std::endl;
        instatsolver->DoStep( v1->Data, p1->Data);
//        instatsolver->DoStep( v1->Data, p1->Data, timestep==0);
        std::cout << "After timestep." << std::endl;
        std::cout << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, t+dt);
        ensight.Write( t+dt);
    }
    delete instatsolver; instatsolver= 0;
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
//    M_pr.Reset();
}

template<class Coeff>
void
StrategyUzawa(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         double theta,
         DROPS::Uint num_timestep,
         double kA, double kM,
         double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
    MLMatDescCL  MG_pr;

    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;
    SSORsmoothCL smoother(1.0);
    SSORPcCL ssor;
    PCG_SsorCL   coarsesolver(ssor, 500, 1e-14);
    MGSolverCL<SSORsmoothCL, PCG_SsorCL> mgc ( smoother, coarsesolver, 1, -1., false);
    typedef SolverAsPreCL<MGSolverCL<SSORsmoothCL, PCG_SsorCL> > MGPCT;
    MGPCT MGPC (mgc);
//    typedef MyUzawaSolver2CL<ISPreCL, PCG_SsorCL> StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPCT> > StatsolverCL;
//    typedef MyUzawaSolver2CL<ISPreCL, MGPCT> StatsolverCL;
    typedef MyUzawaSolver2CL<ISMGPreCL, MGPCT> StatsolverCL;
    StatsolverCL* statsolver= 0;
    typedef InstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
    InstatsolverCL* instatsolver= 0;
    ISMGPreCL* ispcp= 0;
//    ISPreCL* ispcp= 0;
    MakeInitialTriangulation( mg, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    // Initialize Ensight6 output
    std::string ensf( "insa");
    Ensight6OutCL ensight( ensf + ".case", num_timestep + 1);
    ensight.Register( make_Ensight6Geom      ( mg, mg.GetLastLevel(),   "insa geometry", ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( NS.GetPrSolution(),  "p",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( NS.GetVelSolution(), "v",      ensf + ".vel", true));
    ensight.Write( 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt) {
        std::cout << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
//        if (timestep%(num_timestep/10) == 0) { // modify the grid
        if (timestep == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                delete ispcp; ispcp= 0;
                ResetSystem( NS);
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            SetupP2ProlongationMatrix( mg, *mgc.GetProlongation(), vidx1, vidx1);
            M_pr.SetIdx( pidx1, pidx1);
            MG_pr.SetIdx( pidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cout << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);
//            A_pr.SetIdx( pidx1, pidx1);
//            SetupPoissonPressure( mg, A_pr);
//            ispcp= new ISPreCL( A_pr.Data, M_pr.Data, kA, kM, 1.0);
            SetupPoissonPressureMG( NS, MG_pr);
            ispcp= new ISMGPreCL( MG_pr.Data, M_pr.Data, kA, kM, 1);
            SetupP1ProlongationMatrix( mg, *ispcp->GetProlongation(), pidx1, pidx1);
//            PCGSolverCL<ISPreCL> sol1( ispc, stokes_maxiter, stokes_tol);
//            PCG_SsorCL sol2( SSORPcCL( 1.0), stokes_maxiter, stokes_tol);
//            PCGSolverCL<MGPCT> sol2( MGPC, stokes_maxiter, stokes_tol);
//            statsolver= new MyUzawaSolver2CL<ISPreCL, PCG_SsorCL>(
//                               ispc,
//                                sol2,
//                                M_pr.Data, stokes_maxiter, stokes_tol);
//            statsolver= new MyUzawaSolver2CL<ISPreCL, PCGSolverCL<MGPCT> >(
//                                ispc,
//                                sol2,
//                                M_pr.Data, stokes_maxiter, stokes_tol);
//            statsolver= new MyUzawaSolver2CL<ISPreCL,MGPCT>(
//                                *ispcp,
//                                MGPC,
//                                M_pr.Data, stokes_maxiter, stokes_tol);
            statsolver= new MyUzawaSolver2CL<ISMGPreCL, MGPCT>(
                                *ispcp,
                                MGPC,
                                M_pr.Data.GetFinest(), stokes_maxiter, stokes_tol, 1.0);
            instatsolver= new InstatsolverCL( NS, *statsolver, theta);
        }
        instatsolver->SetTimeStep( dt);
        std::cout << "Before timestep." << std::endl;
        instatsolver->DoStep( v1->Data, p1->Data);
        std::cout << "After timestep." << std::endl;
        std::cout << "StatSolver: iterations: " << statsolver->GetIter() << '\n';
        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,
             DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
        ZeroMean( pr);
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, t+dt);
        ensight.Write( t+dt);
    }
    delete instatsolver;
    delete statsolver;
    delete ispcp;
    ResetSystem( NS);
}

template<class Coeff>
void
Strategy(DROPS::StokesP2P1CL<Coeff>& NS,
         int stokes_maxiter, double stokes_tol,
         int poi_maxiter, double poi_tol,
         double theta,
         DROPS::Uint num_timestep,
         double kA, double kM,
         double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
{
    using namespace DROPS;
    typedef StokesP2P1CL<Coeff> StokesCL;

    MultiGridCL& mg= NS.GetMG();
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MLMatDescCL  M_pr;
    MLMatDescCL  MG_pr;
    vidx1->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    TimerCL time;
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;

//    typedef PSchur_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_CL StatsolverCL;
//    typedef PSchur2_PCG_Pr_MG_CL StatsolverCL;
    typedef PSchur2_Full_MG_CL StatsolverCL;
    StatsolverCL* statsolver= 0;
    typedef InstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
//    typedef MyInstatStokesThetaSchemeCL<StokesCL, StatsolverCL> InstatsolverCL;
    InstatsolverCL* instatsolver= 0;
    MakeInitialTriangulation( mg, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.SetNumVelLvl  ( mg.GetNumLevel());
    NS.SetNumPrLvl   ( mg.GetNumLevel());
    M_pr.Data.resize ( mg.GetNumLevel());
    MG_pr.Data.resize( mg.GetNumLevel());
    NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    // Initialize Ensight6 output
    std::string ensf( "insa");
    Ensight6OutCL ensight( ensf + ".case", num_timestep + 1);
    ensight.Register( make_Ensight6Geom      ( mg, mg.GetLastLevel(),   "insa geometry", ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( NS.GetPrSolution(),  "p",      ensf + ".pr",  true));
    ensight.Register( make_Ensight6Vector    ( NS.GetVelSolution(), "v",      ensf + ".vel", true));
    ensight.Write( 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt) {
        std::cout << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
//        if (timestep%(num_timestep/10) == 0) { // modify the grid
        if (timestep == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;

                ResetSystem( NS);
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            M_pr.SetIdx( pidx1, pidx1);
            MG_pr.SetIdx( pidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cout << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);
//            A_pr.SetIdx( pidx1, pidx1);
//            SetupPoissonPressure( mg, A_pr);
//            ISPreCL ispc( A_pr.Data, M_pr.Data, kA, kM, 1.0);
            SetupPoissonPressureMG( NS, MG_pr);
            ISMGPreCL ispc( MG_pr.Data, M_pr.Data, kA, kM, 1);
//            statsolver= new PSchur_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                           poi_maxiter, poi_tol);
//            statsolver= new PSchur2_PCG_CL( M_pr.Data, stokes_maxiter, stokes_tol,
//                                            poi_maxiter, poi_tol);
//            statsolver= new PSchur2_PCG_Pr_CL( ispc, stokes_maxiter, stokes_tol,
//                                               poi_maxiter, poi_tol);
//            statsolver= new PSchur2_PCG_Pr_MG_CL( ispc, stokes_maxiter, stokes_tol,
//                                                  poi_maxiter, poi_tol);
            statsolver= new PSchur2_Full_MG_CL( ispc, stokes_maxiter, stokes_tol,
                                                poi_maxiter, poi_tol);
            SetupP2ProlongationMatrix( mg, *(statsolver->GetPVel()), vidx1, vidx1);
            SetupP1ProlongationMatrix( mg, *ispc.GetProlongation(), pidx1, pidx1);
            instatsolver= new InstatsolverCL( NS, *statsolver, theta);
        }
        instatsolver->SetTimeStep( dt);
        std::cout << "Before timestep." << std::endl;
        instatsolver->DoStep( v1->Data, p1->Data);
//        instatsolver->DoStep( v1->Data, p1->Data, timestep==0);
        std::cout << "After timestep." << std::endl;
        DROPS::P1EvalCL<double, const DROPS::StokesPrBndDataCL,
             DROPS::VecDescCL> pr( &NS.p, &NS.GetBndData().Pr, &mg);
        ZeroMean( pr);
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr, t+dt);
        ensight.Write( t+dt);
    }
    delete instatsolver; instatsolver= 0;
    delete statsolver; statsolver= 0;
    ResetSystem( NS);
}


int main (int argc, char** argv)
{
  try
  {
    if (argc!=13) {
        std::cout <<
"Usage (isadrops): <stokes_maxiter> <stokes_tol> <poi_maxiter> <poi_tol>\n"
"    <theta> <num_timestep> <kA> <kM> <shell_width> <c_level> <f_level>\n"
"    <method>" << std::endl;
        return 1;
    }
    // No C-IO here
    std::ios::sync_with_stdio(false);
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                2,2,2);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
          &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };

    int stokes_maxiter= std::atoi( argv[1]);
    double stokes_tol= std::atof( argv[2]);
    int poi_maxiter= std::atoi( argv[3]);
    double poi_tol= std::atof( argv[4]);
    double theta= std::atof( argv[5]);
    int num_timestep= std::atoi( argv[6]);
    double kA= std::atof( argv[7]);
    double kM= std::atof( argv[8]);
    double shell_width= std::atof( argv[9]);
    int c_level= std::atoi( argv[10]);
    int f_level= std::atoi( argv[11]);
    int method= std::atoi( argv[12]);
    std::cout << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cout << "stokes_tol: " << stokes_tol << ", ";
    std::cout << "poi_maxiter: " << poi_maxiter << ", ";
    std::cout << "poi_tol: " << poi_tol << ", ";
    std::cout << "theta: " << theta << ", ";
    std::cout << "num_timestep: " << num_timestep <<  ", ";
    std::cout << "kA: " << kA <<  ", ";
    std::cout << "kM: " << kM <<  ", ";
    std::cout << "shell_width: " << shell_width <<  ", ";
    std::cout << "c_level: " << c_level << ", ";
    std::cout << "f_level: " << f_level << ", ";
    std::cout << "method: " << method << std::endl;

    typedef DROPS::StokesP2P1CL<MyPdeCL::StokesCoeffCL> NSOnBrickCL;
    typedef NSOnBrickCL MyStokesCL;
    MyStokesCL prob( brick, MyPdeCL::StokesCoeffCL(),
                     DROPS::StokesBndDataCL( 6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();

    switch (method) {
      case 0:
        StrategyUzawa( prob, stokes_maxiter, stokes_tol,
                       theta, num_timestep, kA, kM, shell_width, c_level, f_level);
        break;
      case 1:
        Strategy( prob, stokes_maxiter, stokes_tol, poi_maxiter, poi_tol,
                  theta, num_timestep, kA, kM, shell_width, c_level, f_level);
        break;
      case 2:
        StrategyMRes( prob, stokes_maxiter, stokes_tol,
                      theta, num_timestep, kA, kM, shell_width, c_level, f_level);
        break;
      default:
        std::cout << "Unknown method.\n";
        break;
    }
    std::cout << "hallo" << std::endl;
    std::cout << DROPS::SanityMGOutCL( mg) << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
