//**************************************************************************
// File:    risingBubbleAdap.cpp                                           *
// Content: gravity driven flow of a rising bubble, grid adaptivity        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/instatstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include <fstream>

double      delta_t= 0.05;
DROPS::Uint num_steps= 5;
const int   FPsteps= -1;

// rho*du/dt - mu/Re*laplace u + Dp = f + rho*g - okn
//                          -div u = 0
//                               u = u0, t=t0

// Tropfendaten:
DROPS::Point3DCL Mitte(0.5);
double           Radius= 0.25;

// Glaettungszone fuer Dichte-/Viskositaetssprung
const double sm_eps= 0.05;


class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    DROPS::SmoothedJumpCL rho, mu;
    const double Re, We;
    DROPS::Point3DCL g;

    ZeroFlowCL() 
      : rho( DROPS::JumpCL( 1, 10), DROPS::H_sm, sm_eps),
         mu( DROPS::JumpCL( 1, 2), DROPS::H_sm, sm_eps),
        Re(1), We(1) 
    { g[2]= -9.81; }
};

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)   { return DROPS::SVectorCL<3>(0.); }

double DistanceFct( const DROPS::Point3DCL& p)
{
    const DROPS::Point3DCL d= Mitte-p;
    return d.norm()-Radius;
}


namespace DROPS // for Strategy
{

class AdapTriangCL
{
  private:
    MultiGridCL& mg_;
    double width_;
    Uint c_level_, f_level_;
    bool modified_;
    
    template <class DistFctT>
    double GetValue( DistFctT& dist, const VertexCL& v)
    {
        return dist.val( v);
    }
    
    template <class DistFctT>
    double GetValue( DistFctT& dist, const TetraCL& t)
    {
        return dist.val( t, 0.25, 0.25, 0.25);
    }
    
    double GetValue( scalar_fun_ptr dist, const VertexCL& v)
    {
        return dist( v.GetCoord() );
    }

    double GetValue( scalar_fun_ptr dist, const TetraCL& t)
    {
        return dist( GetBaryCenter( t) );
    }

  public:
    AdapTriangCL( MultiGridCL& mg, double width, Uint c_level, Uint f_level)
      : mg_(mg), width_(width), c_level_(c_level), f_level_(f_level), modified_(false) 
      { Assert( 0<=c_level && c_level<=f_level, "AdapTriangCL: Levels are cheesy.\n", ~0); }
    
    template <class DistFctT>
    void MakeInitialTriang( DistFctT& Dist)
    {
        TimerCL time;

        time.Reset();
        time.Start();
        const Uint min_ref_num= f_level_ - c_level_;
        Uint i;
        for (i=0; i<2*min_ref_num; ++i)
            ModifyGridStep( Dist);
        time.Stop();
        std::cout << "MakeTriang: " << i
                  << " refinements in " << time.GetTime() << " seconds\n"
                  << "last level: " << mg_.GetLastLevel() << '\n';
        mg_.SizeInfo( std::cout);
    }
    
    template <class DistFctT>
    bool ModifyGridStep( DistFctT& Dist)
    // One step of grid change; returns true if modifications were necessary,
    // false, if nothing changed.
    {
        bool modified= false;
        for (MultiGridCL::TriangTetraIteratorCL it= mg_.GetTriangTetraBegin(),
             end= mg_.GetTriangTetraEnd(); it!=end; ++it) 
        {
            double d= 1.;
            for (Uint j=0; j<4; ++j) 
                d= std::min( d, std::abs( GetValue( Dist, *it->GetVertex( j)) ));
            d= std::min( d, std::abs( GetValue( Dist, *it)));
            const Uint l= it->GetLevel();
	    // In the shell:      level should be f_level_.
            // Outside the shell: level should be c_level_.
            const Uint soll_level= d<=width_ ? f_level_ : c_level_;
            
            if (l !=  soll_level)
	    { // tetra will be marked for refinement/remove
	        modified= true;
                if (l < soll_level) 
                    it->SetRegRefMark();
                else // l > soll_level 
                    it->SetRemoveMark();
            }
        }
        if (modified) 
	    mg_.Refine();
        return modified;
    }

    template <class StokesT>
    void UpdateTriang( StokesT& NS, LevelsetP2CL& lset)
    {
        TimerCL time;

        time.Reset();
        time.Start();
        VelVecDescCL  loc_v;
        VecDescCL     loc_p;
	VecDescCL     loc_l;
        VelVecDescCL *v1= &NS.v, 
                     *v2= &loc_v;
        VecDescCL    *p1= &NS.p,
                     *p2= &loc_p,
		     *l1= &lset.Phi,
		     *l2= &loc_l;
        IdxDescCL  loc_vidx( 3, 3), loc_pidx( 1), loc_lidx( 1, 1);
        IdxDescCL  *vidx1= v1->RowIdx,
                   *vidx2= &loc_vidx,
		   *pidx1= p1->RowIdx,
		   *pidx2= &loc_pidx,
		   *lidx1= l1->RowIdx,
		   *lidx2= &loc_lidx;
        modified_= false;
        const Uint min_ref_num= f_level_ - c_level_;
        const InstatStokesBndDataCL& BndData= NS.GetBndData();
        Uint i, LastLevel= mg_.GetLastLevel();
        
        for (i=0; i<2*min_ref_num; ++i)
        {            
	    LevelsetP2CL::DiscSolCL sol(l1,&lset.GetBndData(),&mg_);
            if (!ModifyGridStep(sol))
                break;
            LastLevel= mg_.GetLastLevel();
            modified_= true;
	    
            // Repair velocity
            std::swap( v2, v1);
            std::swap( vidx2, vidx1);
            NS.CreateNumberingVel( LastLevel, vidx1);
            if ( LastLevel != vidx2->TriangLevel) 
            {
                std::cout << "LastLevel: " << LastLevel
                          << " vidx2->TriangLevel: " << vidx2->TriangLevel << std::endl;
                throw DROPSErrCL( "AdapTriangCL::UpdateTriang: Sorry, not yet implemented.");
            }
            v1->SetIdx( vidx1);
            typename StokesT::DiscVelSolCL funvel( v2, &BndData.Vel, &mg_, NS.t);
            RepairAfterRefineP2( funvel, *v1);
            v2->Clear();
            NS.DeleteNumberingVel( vidx2);

            // Repair pressure
            std::swap( p2, p1);
            std::swap( pidx2, pidx1);
            NS.CreateNumberingPr( LastLevel, pidx1);
            p1->SetIdx( pidx1);
            typename StokesT::DiscPrSolCL funpr( p2, &BndData.Pr, &mg_);
            RepairAfterRefineP1( funpr, *p1);
            p2->Clear();
            NS.DeleteNumberingPr( pidx2);

            // Repair levelset
            std::swap( l2, l1);
            std::swap( lidx2, lidx1);
            lset.CreateNumbering( LastLevel, lidx1);
            l1->SetIdx( lidx1);
            LevelsetP2CL::DiscSolCL funlset( l2, &lset.GetBndData(), &mg_);
            RepairAfterRefineP2( funlset, *l1);
            l2->Clear();
            lset.DeleteNumbering( lidx2);
        }
        // We want the solution to be in NS.v, NS.pr, lset.Phi
        if (v1 == &loc_v) 
        {
            NS.vel_idx.swap( loc_vidx);
            NS.pr_idx.swap( loc_pidx);
	    lset.idx.swap( loc_lidx);
            NS.v.SetIdx( &NS.vel_idx);
            NS.p.SetIdx( &NS.pr_idx);
	    lset.Phi.SetIdx( &lset.idx);

            NS.v.Data= loc_v.Data;
            NS.p.Data= loc_p.Data;
	    lset.Phi.Data= loc_l.Data;
        }
        time.Stop();
        std::cout << "UpdateTriang: " << i
                  << " refinements/interpolations in " << time.GetTime() << " seconds\n"
                  << "last level: " << LastLevel << '\n';
        mg_.SizeInfo( std::cout);
    }

    bool WasModified() const { return modified_; }
};

template<class Coeff>
void Strategy( InstatStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adap, double inner_iter_tol, double sigma)
// flow control
{
    typedef InstatStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    LevelsetP2CL lset( MG, sigma, 0.5, 0.1); // Crank-Nicholson, SD=0.1

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;
    MatDescCL  prM;
    TimerCL time;
    
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);    
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx);    
    lset.CreateNumbering(      MG.GetLastLevel(), lidx);

    lset.Phi.SetIdx( lidx);
    Stokes.p.SetIdx( pidx);
    Stokes.v.SetIdx( vidx);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";
    prM.SetIdx( pidx, pidx);
    Stokes.SetupPrMass( &prM);
    
    Stokes.InitVel( &Stokes.v, Null);
    lset.Init( DistanceFct);
    
    time.Reset();

    double outer_tol;
    std::cerr << "tol = "; std::cin >> outer_tol;

    lset.GetSolver().SetTol( 1e-14);
    lset.GetSolver().SetMaxIter( 50000);

    IdxDescCL ens_idx( 1, 1);
    lset.CreateNumbering( MG.GetLastLevel(), &ens_idx);
    EnsightP2SolOutCL ensight( MG, &ens_idx);
    
    const char datgeo[]= "ensight/risebubbleadap.geo", 
               datpr[] = "ensight/risebubbleadap.pr",
               datvec[]= "ensight/risebubbleadap.vec",
               datscl[]= "ensight/risebubbleadap.scl";
    ensight.CaseBegin( "risebubbleadap.case", num_steps+1);
    ensight.DescribeGeom( "zero flow", datgeo,  true);
    ensight.DescribeScalar( "Levelset", datscl, true); 
    ensight.DescribeScalar( "Pressure", datpr,  true); 
    ensight.DescribeVector( "Velocity", datvec, true); 
    ensight.putGeom( datgeo, 0);
    ensight.putVector( datvec, Stokes.GetVelSolution(), 0);
    ensight.putScalar( datpr,  Stokes.GetPrSolution(), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);

    PSchur_GSPCG_CL schurSolver( prM.Data, 200, outer_tol, 200, inner_iter_tol);

    CouplLevelsetStokes2PhaseCL<StokesProblemT, PSchur_GSPCG_CL> 
        cpl( Stokes, lset, schurSolver);
    cpl.SetTimeStep( delta_t);

    for (Uint step= 1; step<=num_steps; ++step)
    {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        cpl.DoStep( FPsteps);
/*        if ((step%5)==0) 
        {
            ensight.putGeom( datgeo, (step-0.01)*delta_t);
            ensight.putScalar( datpr, Stokes.GetPrSolution(), (step-0.01)*delta_t);
            ensight.putVector( datvec, Stokes.GetVelSolution(), (step-0.01)*delta_t);
            ensight.putScalar( datscl, lset.GetSolution(), (step-0.01)*delta_t);
            lset.ReparamSaveIF();
        }
*/        ensight.putGeom( datgeo, step*delta_t);
        ensight.putScalar( datpr, Stokes.GetPrSolution(), step*delta_t);
        ensight.putVector( datvec, Stokes.GetVelSolution(), step*delta_t);
        ensight.putScalar( datscl, lset.GetSolution(), step*delta_t);
	if (step<num_steps) // omit in last step
	{
//            LevelsetP2CL::DiscSolCL sol= lset.GetSolution();
            lset.DeleteNumbering( &ens_idx);
            adap.UpdateTriang( Stokes, lset);
            lset.CreateNumbering( MG.GetLastLevel(), &ens_idx);
/*            
        ensight.putGeom( datgeo, (step+0.01)*delta_t);
        ensight.putScalar( datpr, Stokes.GetPrSolution(), (step+0.01)*delta_t);
        ensight.putVector( datvec, Stokes.GetVelSolution(), (step+0.01)*delta_t);
        ensight.putScalar( datscl, lset.GetSolution(), (step+0.01)*delta_t);
*/
            if (adap.WasModified() )
            {
                cpl.Update();
                // don't forget to update the pr mass matrix for the schur compl. preconditioner!!
                prM.SetIdx( pidx, pidx);
                Stokes.SetupPrMass( &prM);
            }
	}
    }
    ensight.CaseEnd();
    
    std::cerr << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc<4)
    {
        std::cerr << "You have to specify at least three parameters:\n\t" 
                  << argv[0] << " <inner_iter_tol> <num_subdiv> <surf.tension> [<dt> <num_steps>]" << std::endl;
        return 1;
    }
    double inner_iter_tol= atof(argv[1]);
    int sub_div= atoi(argv[2]);
    double sigma= atof(argv[3]);
    if (argc>4) delta_t= atof(argv[4]);
    if (argc>5) num_steps= atoi(argv[5]);

    std::cerr << "inner iter tol:  " << inner_iter_tol << std::endl;
    std::cerr << "sub divisions:   " << sub_div << std::endl;
    std::cerr << "surface tension: " << sigma << std::endl;
    std::cerr << num_steps << " time steps of size " << delta_t << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= 1.; e3[2]= 2.;

    typedef DROPS::InstatStokes2PhaseP2P1CL<ZeroFlowCL> 
            StokesOnBrickCL;
    typedef StokesOnBrickCL MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, sub_div, sub_div, DROPS::Uint(sub_div*e3[2]));

    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::InstatStokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Null }; 
        
    StokesOnBrickCL prob(brick, ZeroFlowCL(), DROPS::InstatStokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::AdapTriangCL adap( mg, 0.2, 0, 3);

    adap.MakeInitialTriang( DistanceFct);
    Strategy( prob, adap, inner_iter_tol, sigma);
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
