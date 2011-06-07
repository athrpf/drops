/// \file
/// \brief Solve a stationary Navier-Stokes-problem
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
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include "num/nssolver.h"
#include "navstokes/navstokes.h"
#include <fstream>

typedef double    (*instat_scalar_fun_ptr)(const DROPS::Point3DCL&, double);
typedef DROPS::Point3DCL (*instat_vector_fun_ptr)(const DROPS::Point3DCL&, double);

static double Reaction(const DROPS::Point3DCL&, double =0)
{
    return 0.0;
}
   
DROPS::Point3DCL Source( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.0);
    ret[2]= 3*p[2];
    return ret; 
}


static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]= p[0];
    ret[1]= p[1];
    ret[2]= -2*p[2];
    return ret;
}
    
static double LsgPr(const DROPS::Point3DCL& p, double)
{
    return -(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2;
}

// q*u - nu*laplace u + (u*D)u + Dp = f
//                           -div u = 0
class StokesCoeffCL
{
  public:
    //reaction
    static instat_scalar_fun_ptr q;
    //source term
    static instat_vector_fun_ptr f;
        
    const double nu;

    StokesCoeffCL(): nu(1.0)
    {
       q = &Reaction,
       f = &Source;
    }

};

instat_scalar_fun_ptr StokesCoeffCL::q;
instat_vector_fun_ptr StokesCoeffCL::f;


typedef DROPS::StokesP2P1CL<StokesCoeffCL>
        StokesOnBrickCL;
typedef StokesOnBrickCL MyStokesCL;


namespace DROPS // for Strategy
{

using ::MyStokesCL;

class Uzawa_PCG_CL : public UzawaSolverCL<PCG_SsorCL>
{
  private:
    SSORPcCL   _ssor;
    PCG_SsorCL _PCGsolver;
  public:
    Uzawa_PCG_CL( MatrixCL& M, int outer_iter, double outer_tol, int inner_iter, double inner_tol, double tau= 1., double omega=1.)
        : UzawaSolverCL<PCG_SsorCL>( _PCGsolver, M, outer_iter, outer_tol, tau),
          _ssor( omega), _PCGsolver(_ssor, inner_iter, inner_tol)
        {}
};

template <class NavStokesT>
class AFPDeCo_Uzawa_PCG_CL: public AdaptFixedPtDefectCorrCL<NavStokesT>
{
  private:
    Uzawa_PCG_CL _uzawaSolver;

  public:
    AFPDeCo_Uzawa_PCG_CL( NavStokesT& NS, MatrixCL& M, int fp_maxiter, double fp_tol, int stokes_maxiter,
                          int poiss_maxiter, double poiss_tol, double reduction= 0.1)
        : AdaptFixedPtDefectCorrCL<NavStokesT>( NS, _uzawaSolver, fp_maxiter, fp_tol, reduction),
          _uzawaSolver( M, stokes_maxiter, fp_tol, poiss_maxiter, poiss_tol) // outer_tol will be set by the AFPDeCo-solver!
        {}
};

template<class Coeff>
void StrategyNavSt(NavierStokesP2P1CL<Coeff>& NS, int maxStep, double fp_tol, int fp_maxiter,
                                                 double uzawa_red, double poi_tol, int poi_maxiter)
// flow control
{
    typedef NavierStokesP2P1CL<Coeff> NavStokesCL;
    MultiGridCL& MG= NS.GetMG();

    MLIdxDescCL  loc_vidx, loc_pidx;
    MLIdxDescCL* vidx1= &NS.vel_idx;
    MLIdxDescCL* pidx1= &NS.pr_idx;
    MLIdxDescCL* vidx2= &loc_vidx;
    MLIdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &NS.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &NS.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &NS.b;
    VelVecDescCL* c= &NS.c;

    MLMatDescCL* A= &NS.A;
    MLMatDescCL* B= &NS.B;
    MLMatDescCL* N= &NS.N;
    MLMatDescCL* M= &NS.M;
    VelVecDescCL  cplM( vidx1);
    int step= 0;

    vidx1->SetFE( vecP2_FE);
    vidx2->SetFE( vecP2_FE);
    pidx1->SetFE( P1_FE);
    pidx2->SetFE( P1_FE);

    TimerCL time;
    do
    {
        MG.Refine();
        NS.CreateNumberingVel( MG.GetLastLevel(), vidx1);
        NS.CreateNumberingPr ( MG.GetLastLevel(), pidx1);
        std::cout << "altes und neues TriangLevel: " << vidx2->TriangLevel() << ", "
                  << vidx1->TriangLevel() << std::endl;
        MG.SizeInfo(std::cout);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        cplM.SetIdx(vidx1);
        std::cout << "Anzahl der Druck-Unbekannten: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cout << "Anzahl der Geschwindigkeitsunbekannten: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;
        if (v2->RowIdx)
        {
            const StokesBndDataCL& BndData= NS.GetBndData();
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
            Interpolate(pr1, pr2);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        M->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        N->SetIdx(vidx1, vidx1);
        time.Reset();
        time.Start();

        NS.SetupSystem1( A, M, b, b, &cplM, NS.v.t);
        NS.SetupSystem2( B, c, NS.v.t);
        time.Stop();
        std::cout << time.GetTime() << " seconds for setting up all systems!" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cout << " A*x took " << time.GetTime() << " seconds!" << std::endl;
        time.Reset();
        time.Start();
        transp_mul( A->Data, v1->Data);
        time.Stop();
        std::cout << "AT*x took " << time.GetTime() << " seconds!" << std::endl;

//        { // write system in files for MatLab
//            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), Ndat("Nmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
//            Adat << A->Data;   Bdat << B->Data; Ndat << N->Data;   bdat << b->Data;    cdat << c->Data;
//        }

        // adaptive fixedpoint defect correction
        //---------------------------------------
        time.Reset();

        NS.prM.SetIdx( pidx1, pidx1);
        NS.SetupPrMass( &NS.prM);
        AFPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver(NS, NS.prM.Data.GetFinest(), fp_maxiter, fp_tol,
                                                     2000, poi_maxiter, poi_tol, uzawa_red);
//        FPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver(NS, NS.prM.Data, fp_maxiter, fp_tol,
//                                                  2000, poi_maxiter, poi_tol, uzawa_red);
//        AFPDeCo_Schur_PCG_CL<NavStokesCL> statsolver(NS, fp_maxiter, fp_tol,
//                                                   2000, poi_maxiter, poi_tol, uzawa_red);
//        FPDeCo_Schur_PCG_CL<NavStokesCL> statsolver(NS, fp_maxiter, fp_tol,
//                                                  2000, poi_maxiter, poi_tol, uzawa_red);
        VelVecDescCL old_v1( vidx1);
        old_v1.Data= v1->Data;
        VelVecDescCL rhsN( vidx1);
        NS.SetupNonlinear( N, v1, &rhsN);
        statsolver.Solve( A->Data, B->Data, *v1, p1->Data, b->Data, rhsN, c->Data, 1.);

/*
        VectorCL d( vidx1->NumUnknowns), e( pidx1->NumUnknowns),
                 w( vidx1->NumUnknowns), q( pidx1->NumUnknowns);
        VelVecDescCL rhsN( vidx1), v_omw( vidx1);
        MLMatDescCL M;
        M.SetIdx( pidx1, pidx1);
        NS.SetupMass( &M);
        double omega= 1, res; // initial value (no damping)
        Uzawa_IPCG_CL uzawaSolver(M.Data, 500, -1., poi_maxiter, poi_tol, 1.);
        for(int fp_step=0; fp_step<fp_maxiter; ++fp_step)
        {
            NS.SetupNonlinear( N, v1, &rhsN);
            MatrixCL AN;
            AN.LinComb( 1., A->Data, 1., N->Data);

            // calculate defect:
            d= AN*v1->Data + transp_mul( B->Data, p1->Data) - b->Data - rhsN.Data;
//            e= B->Data*v1->Data                             - c->Data;
            z_xpay(e, B->Data*v1->Data, -1.0, c->Data);

            std::cout << "fp_step: " << fp_step << ", res = " << (res= std::sqrt(d*d + e*e) ) << std::endl;
            if (res < fp_tol )
                break;

            // solve correction:
            double uzawa_tol= res/uzawa_red;
            if (uzawa_tol < fp_tol) uzawa_tol= fp_tol;
            uzawaSolver.SetTol(uzawa_tol);
            uzawaSolver.Init_A_Pc(AN); // only for Uzawa_IPCG_CL.
            uzawaSolver.Solve(AN, B->Data, w, q, d, e);
            std::cout << "iteration stopped after step " << uzawaSolver.GetIter()
                      << " with res = " << uzawaSolver.GetResid() << std::endl;

            // calculate adaption:
            N->Data.clear();
//            v_omw.Data= v1->Data - omega*w;
            z_xpay(v_omw.Data, v1->Data, -omega, w);
            NS.SetupNonlinear( N, &v_omw, &rhsN);

//            d= A->Data*w + N->Data*w + transp_mul( B->Data, q);
            z_xpaypby2(d, A->Data*w, 1.0, N->Data*w, 1.0, transp_mul(B->Data, q) );
            e= B->Data*w;
            omega= d*(A->Data*v1->Data) + d*(N->Data*v1->Data) + d*transp_mul( B->Data, p1->Data) + e*(B->Data*v1->Data)
                 - d*b->Data - d*rhsN.Data - e*c->Data;
            omega/= d*d + e*e;
            std::cout << "omega = " << omega << std::endl;

            // update solution:
//            v1->Data-= omega*w;
            axpy(-omega, w, v1->Data);
//            p1->Data-= omega*q;
            axpy(-omega, q, p1->Data);
        }
*/

        time.Stop();
        std::cout << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
        NS.CheckSolution(v1, vidx1, p1, &LsgVel, &LsgPr);
        MarkAll(MG);

        A->Reset();
        M->Reset();
        B->Reset();
        b->Reset();
        c->Reset();
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cout << std::endl;
    }
    while (++step<maxStep);
    // we want the solution to be in NS.v, NS.pr
    if (v2 == &loc_v)
    {
        NS.vel_idx.swap( loc_vidx);
        NS.pr_idx.swap( loc_pidx);
        NS.v.SetIdx(&NS.vel_idx);
        NS.p.SetIdx(&NS.pr_idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=7)
    {
        std::cout << "Usage (navstokes): <fp_tol> <poi_tol> <fp_maxiter> <poi_maxiter> <uzawa_red> <num_refinement>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 3, 3, 3);
    const bool IsNeumann[6]=
        {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel };

    DROPS::RBColorMapperCL colormap;

        double fp_tol= std::atof(argv[1]);
        double poi_tol= std::atof(argv[2]);
        int fp_maxiter= std::atoi(argv[3]);
        int poi_maxiter= std::atoi(argv[4]);
        double uzawa_red= std::atof(argv[5]);
        int num_ref= std::atoi(argv[6]);
        std::cout << "fp_tol: " << fp_tol<< ", ";
        std::cout << "poi_tol: " << poi_tol << ", ";
        std::cout << "fp_maxiter: " << fp_maxiter << ", ";
        std::cout << "poi_maxiter: " << poi_maxiter << ", ";
        std::cout << "uzawa_red: " << uzawa_red << ", ";
        std::cout << "num_ref: " << num_ref << std::endl;

        typedef DROPS::NavierStokesP2P1CL<StokesCoeffCL>
                NSOnBrickCL;
        typedef NSOnBrickCL MyNavierStokesCL;

        MyNavierStokesCL prob(brick, StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
        DROPS::MultiGridCL& mg = prob.GetMG();

        StrategyNavSt(prob, num_ref, fp_tol, fp_maxiter, uzawa_red, poi_tol, poi_maxiter);

        std::cout << "hallo" << std::endl;
        std::cout << DROPS::SanityMGOutCL(mg) << std::endl;
        std::ofstream fil("navstokespr.off");
        double min= prob.p.Data.min(),
               max= prob.p.Data.max();
        fil << DROPS::GeomSolOutCL<MyNavierStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;
        fil.close();

        DROPS::MLIdxDescCL tecIdx( DROPS::P1_FE);
        prob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);

        std::ofstream v2d("navstokestec2D.dat");
        DROPS::TecPlot2DSolOutCL< MyNavierStokesCL::const_DiscVelSolCL, MyNavierStokesCL::const_DiscPrSolCL>
            tecplot2d( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx.GetFinest(), -1, 1, 0.5); // cutplane is y=0.5
        v2d << tecplot2d;
        v2d.close();
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
