#include "geom/multigrid.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "geom/builder.h"
#include "stokes/instatstokes.h"
#include "num/nssolver.h"
#include "navstokes/instatnavstokes.h"
#include "navstokes/integrTime.h"
#include <fstream>
#include <sstream>


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
    static inline DROPS::SMatrixCL<3, 3> DxLsgVel(const DROPS::Point3DCL&, double t)
    {
        DROPS::SMatrixCL<3, 3> ret(0.0);
        ret(0,0)= 2.*t - 1.;
        ret(1,1)= 2.*t - 1.;
        ret(2,2)= 2. - 4.*t;
        return ret;
    }

    // Time-derivative of exact solution
    static inline DROPS::SVectorCL<3> DtLsgVel(const DROPS::Point3DCL& p, double)
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


// Interface to EnsightP2SolOutCL for time-dependent geometry,
// velocity and pressure
class EnsightWriterCL
{
  private:
    std::string casefile_;
    std::string geomfile_;
    std::string prfile_;
    std::string velfile_;
    DROPS::MultiGridCL& MG_;
    DROPS::IdxDescCL ensightidx_;
    DROPS::EnsightP2SolOutCL ensight_;
    bool have_idx_;

  public:
    EnsightWriterCL(DROPS::MultiGridCL& MG, DROPS::Uint num_timestep,
                    std::string casefile, std::string geomfile,
                    std::string prfile, std::string velfile);
    ~EnsightWriterCL();

    // Call after a grid-change before writing data. The constructor calls this, too.
    void CreateNumbering(int level= -1);
    // Destroy the index before modifying the grid.
    void DeleteNumbering();

    // To write geometry, pressure and velocity at time t.
    template<class InstatNSCL>
    void
    WriteAtTime(const InstatNSCL& NS, const double t);
};

EnsightWriterCL::EnsightWriterCL(DROPS::MultiGridCL& MG, DROPS::Uint num_timestep,
                                 std::string casefile, std::string geomfile,
                                 std::string prfile, std::string velfile)
    :casefile_( casefile), geomfile_( geomfile), prfile_( prfile), velfile_( velfile),
     MG_( MG), ensight_( MG, &ensightidx_), have_idx_( false)
{
    ensightidx_.Set( 1,1,0,0);
    this->CreateNumbering();
    have_idx_= true;
    ensight_.CaseBegin( casefile_.c_str(), num_timestep);
    ensight_.DescribeGeom( "insa_geometry", geomfile_, true);
    ensight_.DescribeScalar( "p", prfile_, true);
    ensight_.DescribeVector( "v", velfile_, true);
}

EnsightWriterCL::~EnsightWriterCL()
{
    if (!have_idx_)
        std::cerr << "EnsightWriter::~EnsightWriterCL: Error; no index found.\n";
    ensight_.CaseEnd();
    this->DeleteNumbering();
}

template<class InstatNSCL>
void
EnsightWriterCL::WriteAtTime(const InstatNSCL& NS, const double t)
{
    if (!have_idx_)
        throw DROPS::DROPSErrCL( "EnsightWriter::WriteAtTime: Call CreateNumbering first.");
    ensight_.putGeom( geomfile_, t);
    typename InstatNSCL::DiscPrSolCL ensightp( &NS.p, &NS.GetBndData().Pr, &MG_);
    ensight_.putScalar( prfile_, ensightp, t);
    typename InstatNSCL::DiscVelSolCL ensightv( &NS.v, &NS.GetBndData().Vel, &MG_, t);
    ensight_.putVector( velfile_, ensightv, t);
}

void
EnsightWriterCL::CreateNumbering( int level)
{
    if (have_idx_)
        throw DROPS::DROPSErrCL( "EnsightWriter::CreateIndex: Already done.");
    DROPS::NoBndDataCL<> ensightbnd;
    ensightidx_.TriangLevel= level < 0 ? MG_.GetLastLevel(): level;
    DROPS::CreateNumbOnVertex( ensightidx_.GetIdx(), ensightidx_.NumUnknowns, 1,
                               MG_.GetTriangVertexBegin( ensightidx_.TriangLevel),
                               MG_.GetTriangVertexEnd( ensightidx_.TriangLevel),
                               ensightbnd);
    DROPS::CreateNumbOnEdge( ensightidx_.GetIdx(), ensightidx_.NumUnknowns, 1,
                             MG_.GetTriangEdgeBegin( ensightidx_.TriangLevel),
                             MG_.GetTriangEdgeEnd( ensightidx_.TriangLevel),
                             ensightbnd);
    have_idx_= true;
}

void
EnsightWriterCL::DeleteNumbering()
{
    if (!have_idx_)
        throw DROPS::DROPSErrCL( "EnsightWriter::WriteAtTime: Call CreateNumbering first.");
    DROPS::DeleteNumbOnSimplex( ensightidx_.GetIdx(),
                                MG_.GetAllVertexBegin( ensightidx_.TriangLevel),
                                MG_.GetAllVertexEnd( ensightidx_.TriangLevel));
    DROPS::DeleteNumbOnSimplex( ensightidx_.GetIdx(),
                                MG_.GetAllEdgeBegin( ensightidx_.TriangLevel),
                                MG_.GetAllEdgeEnd( ensightidx_.TriangLevel));
    ensightidx_.NumUnknowns= 0;
    have_idx_= false;
}


InstatNSCL::StokesCoeffCL InstatNSCL::Coeff;

typedef InstatNSCL MyPdeCL;

typedef DROPS::SVectorCL<3> (*fun_ptr)(const DROPS::SVectorCL<3>&, double);

int
CheckVel(DROPS::InstatP2EvalCL< DROPS::SVectorCL<3>,
                                const DROPS::InstatStokesVelBndDataCL,
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

void
SetVel(DROPS::InstatP2EvalCL< DROPS::SVectorCL<3>,
                              const DROPS::InstatStokesVelBndDataCL,
                              DROPS::VelVecDescCL>& fun,
       double t)
{
    const DROPS::Uint lvl= fun.GetLevel();
    DROPS::MultiGridCL& mg= const_cast<DROPS::MultiGridCL&>( fun.GetMG());
    for (DROPS::MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        if (!fun.GetBndData()->IsOnDirBnd( *sit))
            fun.SetDoF( *sit, InstatNSCL::LsgVel( sit->GetCoord(), t));
    }
    for (DROPS::MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        if (!fun.GetBndData()->IsOnDirBnd( *sit))
            fun.SetDoF( *sit, InstatNSCL::LsgVel( (sit->GetVertex( 0)->GetCoord()
                                                   + sit->GetVertex( 1)->GetCoord())*0.5, t));
    }
}

void
SetPr(DROPS::P1EvalCL< double,
                       const DROPS::InstatStokesPrBndDataCL,
                       DROPS::VecDescCL>& fun,
      double t)
{
    const DROPS::Uint lvl= fun.GetLevel();
    DROPS::MultiGridCL& mg= const_cast<DROPS::MultiGridCL&>( fun.GetMG());
    for (DROPS::MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, InstatNSCL::LsgPr( sit->GetCoord(), t));
    }
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
                    else; // nothing
            }
            else { // Outside the shell; level should be c_level;
                if (l < c_level) { it->SetRegRefMark(); }
                else
                    if (l> c_level) { it->SetRemoveMark(); }
                    else; // nothing
            }
        }
        mg.Refine();
    return shell_not_ready;
}

template<class Coeff>
void
UpdateTriangulation(DROPS::InstatNavierStokesP2P1CL<Coeff>& NS,
                    const signed_dist_fun Dist,
                    const double t,
                    const double width,         // Thickness of refined shell on eache side of the interface
                    const DROPS::Uint c_level,  // Outside the shell, use this level
                    const DROPS::Uint f_level,  // Inside the shell, use this level
                    DROPS::VelVecDescCL* v1,
                    DROPS::VecDescCL* p1)
{
    using namespace DROPS;
    Assert( 0<=c_level && c_level<=f_level, "UpdateTriangulation: Levels are cheesy.\n", ~0);
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
    vidx2->Set( 3, 3, 0, 0);
    pidx2->Set( 1, 0, 0, 0);
    bool shell_not_ready= true;
    const Uint min_ref_num= f_level - c_level;
    const InstatStokesBndDataCL& BndData= NS.GetBndData();
    Uint i;
    for(i=0; shell_not_ready || i<min_ref_num; ++i) {
        shell_not_ready= ModifyGridStep( mg, Dist, width, c_level, f_level, t);
        // Repair velocity
        std::swap( v2, v1);
        std::swap( vidx2, vidx1);
        NS.CreateNumberingVel( mg.GetLastLevel(), vidx1);
        if ( mg.GetLastLevel() != vidx2->TriangLevel) {
            std::cout << "LastLevel: " << mg.GetLastLevel()
                      << " vidx2->TriangLevel: " << vidx2->TriangLevel << std::endl;
            throw DROPSErrCL( "Strategy: Sorry, not yet implemented.");
        }
        v1->SetIdx( vidx1);
        InstatP2EvalCL< SVectorCL<3>, const InstatStokesVelBndDataCL, 
                        const VelVecDescCL> funv2( v2, &BndData.Vel, &mg, t);
        RepairAfterRefineP2( funv2, *v1);
        v2->Clear();
        NS.DeleteNumberingVel( vidx2);
//InstatP2EvalCL< SVectorCL<3>, const InstatStokesVelBndDataCL, 
//                VelVecDescCL> funv1( v1, &BndData.Vel, &mg, t);
//CheckVel( funv1, &MyPdeCL::LsgVel);
        // Repair pressure
        std::swap( p2, p1);
        std::swap( pidx2, pidx1);
        NS.CreateNumberingPr( mg.GetLastLevel(), pidx1);
        p1->SetIdx( pidx1);
        P1EvalCL<double, const InstatStokesPrBndDataCL,
                 const VecDescCL> funpr( p2, &BndData.Pr, &mg);
        RepairAfterRefineP1( funpr, *p1);
        p2->Clear();
        NS.DeleteNumberingPr( pidx2);
    }
    // We want the solution to be where v1, p1 point to.
    if (v1 == &loc_v) {
        NS.vel_idx.swap( loc_vidx);
        NS.pr_idx.swap( loc_pidx);
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
    Assert( 0<=c_level && c_level<=f_level, "MakeInitialTriangulation: Levels are cheesy.\n", ~0);
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
SetMatVecIndices(DROPS::InstatNavierStokesP2P1CL<Coeff>& NS,
                 DROPS::IdxDescCL* const vidx,
                 DROPS::IdxDescCL* const pidx)
{
    std::cout << "#Druck-Unbekannte: " << pidx->NumUnknowns << std::endl;
    std::cout << "#Geschwindigkeitsunbekannte: " << vidx->NumUnknowns << std::endl;
    NS.b.SetIdx( vidx);
    NS.c.SetIdx( pidx);
    NS.cplM.SetIdx( vidx);
    NS.cplN.SetIdx( vidx);
    NS.A.SetIdx( vidx, vidx);
    NS.B.SetIdx( pidx, vidx);
    NS.M.SetIdx( vidx, vidx);
    NS.N.SetIdx( vidx, vidx);
}

template<class Coeff>
void
ResetSystem(DROPS::InstatNavierStokesP2P1CL<Coeff>& NS)
{
    NS.A.Reset(); NS.B.Reset();
    NS.M.Reset(); NS.N.Reset();
    NS.b.Reset(); NS.c.Reset();
    NS.cplM.Reset(); NS.cplN.Reset();
}


template<class Coeff>
void
Strategy(DROPS::InstatNavierStokesP2P1CL<Coeff>& NS,
         double fp_tol, int fp_maxiter, 
         double deco_red, int stokes_maxiter,
         double poi_tol, int poi_maxiter,
         double theta,
         DROPS::Uint num_timestep,
         double shell_width, DROPS::Uint c_level, DROPS::Uint f_level)
{
    using namespace DROPS;
    typedef InstatNavierStokesP2P1CL<Coeff> NavStokesCL;
    
    MultiGridCL& MG= NS.GetMG();
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    VelVecDescCL* v1= &NS.v;
    VecDescCL*    p1= &NS.p;
    MatDescCL  M_pr;
    vidx1->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    TimerCL time;
    double t= 0.;
    const double dt= 1./num_timestep;
    NS.t= 0;
    Uint timestep= 0;

    FPDeCo_Uzawa_PCG_CL<NavStokesCL>* statsolver= 0; 
//    FPDeCo_Uzawa_CG_CL<NavStokesCL>* statsolver= 0;
//    FPDeCo_Uzawa_SGSPCG_CL<NavStokesCL>* statsolver= 0;
//    AFPDeCo_Schur_PCG_CL<NavStokesCL>* statsolver= 0;
//    FPDeCo_Schur_PCG_CL<NavStokesCL>* statsolver= 0;
//    InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Schur_PCG_CL<NavStokesCL> >* instatsolver= 0;
    InstatNavStokesThetaSchemeCL<NavStokesCL, FPDeCo_Uzawa_PCG_CL<NavStokesCL> >* instatsolver= 0;

    MakeInitialTriangulation( MG, &SignedDistToInterface, shell_width, c_level, f_level);
    NS.CreateNumberingVel( MG.GetLastLevel(), vidx1);    
    v1->SetIdx( vidx1);
    NS.InitVel( v1, &MyPdeCL::LsgVel);
    NS.CreateNumberingPr( MG.GetLastLevel(), pidx1);
    p1->SetIdx( pidx1);

    EnsightWriterCL ensightout( MG, num_timestep, "insa.case", "insa_geo", "insa_pr", "insa_vel");
    ensightout.WriteAtTime( NS, 0.);

    for (; timestep<num_timestep; ++timestep, t+= dt, NS.t+= dt) {
        std::cerr << "----------------------------------------------------------------------------"
                  << std::endl << "t: " << t << std::endl;
        if (timestep%(num_timestep/10) == 0) { // modify the grid
            if (timestep>0) { // perform cleanup, which is not neccessary for t==0.
                delete statsolver; statsolver= 0;
                delete instatsolver; instatsolver= 0;
                M_pr.Reset();
                ResetSystem( NS);
                ensightout.DeleteNumbering();
                UpdateTriangulation( NS, &SignedDistToInterface, t,
                                     shell_width, c_level, f_level, v1, p1);
                ensightout.CreateNumbering();
            }
            SetMatVecIndices( NS, vidx1, pidx1);
            time.Reset(); time.Start();
            NS.SetupInstatSystem( &NS.A, &NS.B, &NS.M);
            time.Stop();
            std::cerr << "SetupInstatSystem: " << time.GetTime() << " seconds" << std::endl;
            time.Reset();
            M_pr.SetIdx( pidx1, pidx1);
            NS.SetupPrMass( &M_pr);  
//            AFPDeCo_Uzawa_PCG_CL<NavStokesCL> statsolver( NS, M_pr.Data, fp_maxiter, fp_tol,
//                                                          stokes_maxiter, poi_maxiter, poi_tol, deco_red);
            statsolver= new FPDeCo_Uzawa_PCG_CL<NavStokesCL> ( NS, M_pr.Data, fp_maxiter, fp_tol,
                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new FPDeCo_Uzawa_CG_CL<NavStokesCL>( NS, M_pr.Data, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new FPDeCo_Uzawa_SGSPCG_CL<NavStokesCL>( NS, M_pr.Data, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new AFPDeCo_Schur_PCG_CL<NavStokesCL>( NS, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
//            statsolver= new FPDeCo_Schur_PCG_CL<NavStokesCL>( NS, fp_maxiter, fp_tol,
//                           stokes_maxiter, poi_maxiter, poi_tol, deco_red);
            // If the saddlepoint-problem is solved via an Uzawa-method, the mass-matrix alone is
            // not an appropriate preconditioner for the Schur-Complement-Matrix. M has to be scaled
            // by 1/(theta*dt).
            statsolver->GetStokesSolver().SetTau( theta*dt);
            time.Reset(); time.Start();
            NS.SetupNonlinear( &NS.N, v1, &NS.cplN, t, t);
            time.Stop();
            std::cerr << "SetupNonlinear: " << time.GetTime() << " seconds" << std::endl;
            NS.SetupInstatRhs( &NS.b, &NS.c, &NS.cplM, t, &NS.b, t);
//            instatsolver= new InstatNavStokesThetaSchemeCL<NavStokesCL,
//                              FPDeCo_Schur_PCG_CL<NavStokesCL> >(NS, *statsolver, theta);
            instatsolver= new InstatNavStokesThetaSchemeCL<NavStokesCL,
                             FPDeCo_Uzawa_PCG_CL<NavStokesCL> >( NS, *statsolver, theta, t);
            if (timestep == 0) // check initial velocities
                NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DtLsgVel, &MyPdeCL::LsgPr, t);
        }
        NS.SetTime( t+dt); // We have to set the new time!
        instatsolver->SetTimeStep( dt);
        std::cerr << "Before timestep." << std::endl;
        instatsolver->DoStep( *v1, p1->Data);
        std::cerr << "After timestep." << std::endl;
        NS.CheckSolution( v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::DtLsgVel, &MyPdeCL::LsgPr, t+dt);
        ensightout.WriteAtTime( NS, t+dt);
    }
    delete statsolver; statsolver= 0;
    delete instatsolver; instatsolver= 0;
    ResetSystem( NS);
    M_pr.Reset();
}


int main (int argc, char** argv)
{
  try
  {
    if (argc!=12)
    {
        std::cerr << "Usage (insadrops):  <fp_tol> <fp_maxiter> "
	          << "<deco_red> <stokes_maxiter> <poi_tol> <poi_maxiter> "
                  << "<theta> <num_timestep> <shell_width> <c_level> <f_level>" << std::endl;
        return 1;
    }
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
				DROPS::std_basis<3>(2),
				DROPS::std_basis<3>(3),
				2,2,2);
    const bool IsNeumann[6]= {false, false, false, false, false, false};
    const DROPS::InstatStokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel,
	  &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel };
    DROPS::RBColorMapperCL colormap;

    double fp_tol= atof(argv[1]);
    int fp_maxiter= atoi(argv[2]);
    double deco_red= atof(argv[3]);
    int stokes_maxiter= atoi(argv[4]);
    double poi_tol= atof(argv[5]);
    int poi_maxiter= atoi(argv[6]);
    double theta= atof(argv[7]);
    int num_timestep= atoi(argv[8]);
    double shell_width= atof(argv[9]);
    int c_level= atoi(argv[10]);
    int f_level= atoi(argv[11]);

    std::cerr << "fp_tol: " << fp_tol<< ", ";
    std::cerr << "fp_maxiter: " << fp_maxiter << ", ";
    std::cerr << "deco_red: " << deco_red << ", ";
    std::cerr << "stokes_maxiter: " << stokes_maxiter << ", ";
    std::cerr << "poi_tol: " << poi_tol << ", ";
    std::cerr << "poi_maxiter: " << poi_maxiter << ", ";
    std::cerr << "theta: " << theta << ", ";
    std::cerr << "num_timestep: " << num_timestep <<  ", ";
    std::cerr << "shell_width: " << shell_width <<  ", ";
    std::cerr << "c_level: " << c_level << ", ";
    std::cerr << "f_level: " << f_level << std::endl;

    typedef DROPS::InstatNavierStokesP2P1CL<MyPdeCL::StokesCoeffCL> 
    	    NSOnBrickCL;
    typedef NSOnBrickCL MyNavierStokesCL;
    MyNavierStokesCL prob(brick, MyPdeCL::StokesCoeffCL(),
                          DROPS::InstatStokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    
    Strategy(prob, fp_tol, fp_maxiter, deco_red, stokes_maxiter, poi_tol, poi_maxiter,
             theta, num_timestep, shell_width, c_level, f_level);

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
