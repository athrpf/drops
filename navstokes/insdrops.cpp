#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/instatstokes.h"
#include "num/nssolver.h"
#include "navstokes/instatnavstokes.h"
#include "navstokes/integrTime.h"
#include <fstream>


struct InstatNSCL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret;
	ret[0]= (2.*t - 1.)*p[0];
	ret[1]= (2.*t - 1.)*p[1];
	ret[2]= (2. - 4.*t)*p[2];
	return ret;	
    }
    
    // int_{x=0..1, y=0..1,z=0..1} p(x,y,z,t) dx dy dz = 0 for all t.
    static double LsgPr(const DROPS::Point3DCL& p, double t)
    {
        return (0.5 - t)*p.norm_sq() + t - 0.5;
    }

    // Jacobi-matrix of exact solution (only in the spatial variables)
    static inline DROPS::SMatrixCL<3, 3> DxLsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SMatrixCL<3, 3> ret(0.0);
        ret(0,0)= 2.*t - 1.;
        ret(1,1)= 2.*t - 1.;
        ret(2,2)= 2. - 4.*t;
        return ret;
    }

    // Time-derivative of exact solution
    static inline DROPS::SVectorCL<3> DtLsgVel(const DROPS::Point3DCL& p, double t)
    {
        DROPS::SVectorCL<3> ret(0.0);
        ret[0]= 2.*p[0];
        ret[1]= 2.*p[1];
        ret[2]= -4.*p[2];
        return ret;
    }
    // u_t + q*u - nu*laplace u + (u*D)u + Dp = f
    //                                 -div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&, double) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p, double t)
        {
            DROPS::SVectorCL<3> ret;
	    ret[0]= (4.*t*t - 6.*t + 4.)*p[0];
	    ret[1]= (4.*t*t - 6.*t + 4.)*p[1];
	    ret[2]= (16.*t*t -18.*t + 1.)*p[2];
	    return ret; 
	}
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };
    
    static StokesCoeffCL Coeff;
};

InstatNSCL::StokesCoeffCL InstatNSCL::Coeff;

typedef InstatNSCL MyPdeCL;

namespace DROPS // for Strategy
{

using ::MyPdeCL;


template<class Coeff>
void Strategy(InstatNavierStokesP2P1CL<Coeff>& NS, int num_ref, double fp_tol, int fp_maxiter, 
              double deco_red, int stokes_maxiter, double poi_tol, int poi_maxiter,
	      double theta, double dt)
// flow control
{
    typedef InstatNavierStokesP2P1CL<Coeff> NavStokesCL;
    
    MultiGridCL& MG= NS.GetMG();

    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx2= &loc_pidx;

    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &NS.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &NS.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &NS.b;
    VelVecDescCL* cplN= &NS.cplN;
    VelVecDescCL* cplM= &NS.cplM;
    VelVecDescCL* c= &NS.c;

    MatDescCL* A= &NS.A;
    MatDescCL* B= &NS.B;
    MatDescCL* N= &NS.N;
    MatDescCL* M= &NS.M;
    int refstep= 0;

    vidx1->Set( 3, 3, 0, 0);
    vidx2->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    pidx2->Set( 1, 0, 0, 0);
    
    TimerCL time;
    do
    {
        MG.Refine();
        NS.CreateNumberingVel(MG.GetLastLevel(), vidx1);    
        NS.CreateNumberingPr(MG.GetLastLevel(), pidx1);    
        std::cerr << "altes und neues TriangLevel: " << (refstep!=0 ?
	vidx2->TriangLevel : -1) << ", " << vidx1->TriangLevel << std::endl;
        MG.SizeInfo(std::cerr);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        cplN->SetIdx(vidx1);
        cplM->SetIdx(vidx1);
        std::cerr << "#Druck-Unbekannte: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cerr << "#Geschwindigkeitsunbekannte: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;
        if (v2->RowIdx)
        {
            const InstatStokesBndDataCL& BndData= NS.GetBndData();
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesBndDataCL::PrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
            Interpolate(pr1, pr2);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        N->SetIdx(vidx1, vidx1);
        NS.InitVel(v1, &MyPdeCL::LsgVel);
        time.Reset();
        time.Start();
        NS.SetupInstatSystem(A, B, M);
        time.Stop();
        std::cerr << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cerr << " A*x: " << time.GetTime() << " seconds" << std::endl;
        time.Reset();

        MatDescCL M_pr;
        M_pr.SetIdx( pidx1, pidx1);
        NS.SetupPrMass( &M_pr);  
//        AFPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                     stokes_maxiter, poi_maxiter, poi_tol, deco_red);
        FPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
                                                    stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        FPDeCo_Uzawa_CG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                   stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        FPDeCo_Uzawa_SGSPCG_CL<NavStokesCL> statsolver(NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                   stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        AFPDeCo_Schur_PCG_CL<NavStokesCL> statsolver(NS, fp_maxiter, fp_tol,
//                                                   stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//        FPDeCo_Schur_PCG_CL<NavStokesCL> statsolver(NS, fp_maxiter, fp_tol,
//                                                    stokes_maxiter, poi_maxiter, poi_tol, deco_red);
        // If the saddlepoint-problem is solved via an Uzawa-method, the mass-matrix alone is
        // not an appropriate preconditioner for the Schur-Complement-Matrix. M has to be scaled
        // by 1/(theta*dt).
        statsolver.GetStokesSolver().SetTau(theta*dt);
        time.Reset();
        time.Start();
        NS.SetupNonlinear(N, v1, cplN, 0., 0.);
	time.Stop();
        std::cerr << "SetupNonlinear: " << time.GetTime() << " seconds" << std::endl;
        NS.SetupInstatRhs( b, c, cplM, 0., b, 0.);
//        InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Schur_PCG_CL<NavStokesCL> >
//	    instatsolver(NS, statsolver, v1, theta);
        InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Uzawa_PCG_CL<NavStokesCL> >
	    instatsolver(NS, statsolver, v1, theta);
        std::cerr << "After constructor." << std::endl;
	double t= 0.;
	NS.t= 0;
	for (int timestep=0; t<1.; ++timestep, t+= dt, NS.t+= dt)
	{
//if (timestep==25) return;
	    std::cerr << "------------------------------------------------------------------"
	              << std::endl << "t: " << t << std::endl;
            NS.SetTime(t+dt); // We have to set the new time!
	    // GetDiscError calls SetupNonlinear for exact solution, thus we call it
	    // prior to the actual call for the time-evolution.
//            NS.GetDiscError(&MyPdeCL::LsgVel, &MyPdeCL::DtLsgVel, &MyPdeCL::LsgPr, t);
            if (timestep==0) // Check the initial solution, at least velocities.
	        NS.CheckSolution(v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DtLsgVel, &MyPdeCL::LsgPr, t);
	    instatsolver.SetTimeStep(dt);
            std::cerr << "Before timestep." << std::endl;
	    instatsolver.DoStep(*v1, p1->Data);
            std::cerr << "After timestep." << std::endl;
            NS.CheckSolution(v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DtLsgVel, &MyPdeCL::LsgPr, t+dt);
	}
	MarkAll(MG);
        
        A->Reset();
        B->Reset();
	M->Reset();
	N->Reset();
//	M_pr.Reset();
        b->Reset();
        c->Reset();
	cplN->Reset();
	cplM->Reset();
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cerr << std::endl;
    }
    while (refstep++ < num_ref);
    // we want the solution to be in _v
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
    if (argc!=10)
    {
        std::cerr << "Usage (insdrops): <num_refinement> <fp_tol> <fp_maxiter> "
	          << "<deco_red> <stokes_maxiter> <poi_tol> <poi_maxiter> "
                  << "<theta> <dt>" << std::endl;
        return 1;
    }
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
				DROPS::std_basis<3>(2),
				DROPS::std_basis<3>(3),
				4,4,4);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::InstatStokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
	  &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };
    DROPS::RBColorMapperCL colormap;

    int num_ref= atoi(argv[1]);
    double fp_tol= atof(argv[2]);
    int fp_maxiter= atoi(argv[3]);
    double deco_red= atof(argv[4]);
    int stokes_maxiter= atoi(argv[5]);
    double poi_tol= atof(argv[6]);
    int poi_maxiter= atoi(argv[7]);
    double theta= atof(argv[8]);
    double dt= atof(argv[9]);
    std::cerr << "num_ref: " << num_ref << ", ";
    std::cerr << "fp_tol: " << fp_tol<< ", ";
    std::cerr << "fp_maxiter: " << fp_maxiter << ", ";
    std::cerr << "deco_red: " << deco_red << ", ";
    std::cerr << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cerr << "poi_tol: " << poi_tol << ", ";
    std::cerr << "poi_maxiter: " << poi_maxiter << ", ";
    std::cerr << "theta: " << theta << ", ";
    std::cerr << "dt: " << dt << std::endl;

    typedef DROPS::InstatNavierStokesP2P1CL<MyPdeCL::StokesCoeffCL> 
    	    NSOnBrickCL;
    typedef NSOnBrickCL MyNavierStokesCL;
    MyNavierStokesCL prob(brick, MyPdeCL::StokesCoeffCL(),
                          DROPS::InstatStokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    
    Strategy(prob, num_ref, fp_tol, fp_maxiter, deco_red, stokes_maxiter, poi_tol, poi_maxiter,
             theta, dt);

    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("navstokespr.off");
    double min= prob.p.Data.min(),
    	   max= prob.p.Data.max();
    fil << DROPS::GeomSolOutCL<MyNavierStokesCL::DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;
    fil.close();

    DROPS::IdxDescCL tecIdx;
    tecIdx.Set( 1, 0, 0, 0);
    prob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);    
    std::ofstream v2d("navstokestec2D.dat");
    DROPS::TecPlot2DSolOutCL< MyNavierStokesCL::DiscVelSolCL, MyNavierStokesCL::DiscPrSolCL>
    	tecplot2d( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx, -1, 2, 0.5); // cutplane is z=0.5
    v2d << tecplot2d;
    v2d.close();
    v2d.open("navstokestec2D2.dat");
    DROPS::TecPlot2DSolOutCL< MyNavierStokesCL::DiscVelSolCL, MyNavierStokesCL::DiscPrSolCL>
    	tecplot2d2( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx, -1, 1, 0.5); // cutplane is z=0.5
    v2d << tecplot2d2;
    v2d.close();
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
