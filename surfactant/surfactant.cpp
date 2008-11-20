#include "surfactant/ifacetransp.h"
#include "surfactant/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "out/ensightOut.h"

#include <fstream>

using namespace DROPS;

DROPS::ParamSurfactantCL C;
std::string filename,
            datgeo,
            datscl,
            datvel,
            datsurf,
            datsol,
            datlset;

DROPS::Point3DCL u_func (const DROPS::Point3DCL&, double)
{
    return C.Velocity;
}

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc[6]= {
    DROPS::DirBC, DROPS::DirBC,
    DROPS::DirBC, DROPS::DirBC,
    DROPS::DirBC, DROPS::DirBC
};

bnd_val_fun bf[6]= {
    &u_func, &u_func, &u_func, &u_func, &u_func, &u_func
};

double sphere_2 (const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL x( p - C.Mitte);

    return x.norm() - C.Radius[0];
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (C.Mitte + t*u_func(p, t)));

    return x.norm() - C.Radius[0];
}


// TestCase == 0: Sphere around 0, radius 1, v == 0
// A right hand side from C.J. Heine...
const double a( -13./8.*std::sqrt( 35./M_PI));
double rhs0 (const DROPS::Point3DCL& p, double)
{
    return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
// ...and the corresponding solution (extended)
double sol0 (const DROPS::Point3DCL& p, double)
{
    return p.norm_sq()/(12. + p.norm_sq())*rhs0( p, 0.);
}

double sol0t (const DROPS::Point3DCL& p, double t)
{
    const Point3DCL q( p - (C.Mitte + t*u_func(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

    return q.norm_sq()/(12. + q.norm_sq())*val;
}

typedef DROPS::P1EvalCL<double, const DROPS::NoBndDataCL<>, DROPS::VecDescCL> DiscP1FunT;

template<class DiscP1FunType>
double L2_error (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls,
    const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfacePatchCL patch;
    DROPS::Quad5_2DCL<> qsol, qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &patch.GetBary( tri), extsol, t);
                    qdiscsol.assign(  *it, &patch.GetBary( tri), discsol);
                    d+= DROPS::Quad5_2DCL<>( std::pow( qdiscsol - qsol, 2)).quad( patch.GetFuncDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls,
    DROPS::instat_scalar_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfacePatchCL patch;
    DROPS::Quad5_2DCL<> qsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &patch.GetBary( tri), extsol, t);
                    d+= DROPS::Quad5_2DCL<>( qsol*qsol).quad( patch.GetFuncDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::scalar_fun_ptr d)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord());

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= ls.Data[it->Unknowns( idx)]=
            0.5*(ls.Data[it->GetVertex( 0)->Unknowns( idx)] + ls.Data[it->GetVertex( 1)->Unknowns( idx)]);
}

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
}

void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, instat_vector_fun_ptr LsgVel, double t)
{
    VectorCL& lsgvel= vec->Data;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

   DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel(sit->GetCoord(), t));
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t));
    }
}


class TransportP1FunctionCL
/// Solve D/Dt u = 0, u(t^0) given, with SDFEM
{
  public:
    typedef BndDataCL<>    BndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    const IdxDescCL&    p1idx;

  private:
    MultiGridCL&        MG_;
    double              SD_,       ///< streamline diffusion
                        theta_,
                        dt_;
    MatrixCL            L_;
    BndDataT            Bnd_;
    GSPcCL              pc_;
    GMResSolverCL<GSPcCL>  gm_;

  public:
    MatrixCL            E_old, E_, H_old, H_;

    template<class DiscVelSolT>
    TransportP1FunctionCL (MultiGridCL& mg, const DiscVelSolT& v_old, const IdxDescCL& thep1idx, double dt, double theta= 0.5, double SD= 0., int iter= 1000, double tol= 1e-7)
    : p1idx( thep1idx), MG_( mg), SD_( SD),
        theta_( theta), dt_( dt), Bnd_( BndDataT( mg.GetBnd().GetNumBndSeg())),
        gm_( pc_, 500, iter, tol, true)
    {
        if (theta_ != 1.) SetupSystem( v_old, E_old, H_old);
    }

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    template<class DiscVelSolT>
    void SetupSystem (const DiscVelSolT&, MatrixCL& E, MatrixCL& H);
    /// perform one time step
    template <class DiscVelSolT>
    void DoStep (VectorCL& u, const DiscVelSolT& /* new velocity*/);

};

template<class DiscVelSolT>
void TransportP1FunctionCL::SetupSystem( const DiscVelSolT& vel, MatrixCL& E, MatrixCL& H)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= p1idx.NumUnknowns();
    const Uint lvl= p1idx.TriangLevel();

    SparseMatBuilderCL<double> bE(&E, num_unks, num_unks),
                               bH(&H, num_unks, num_unks);
    IdxT Numb[4];

    std::cerr << "entering TransportP1Function::SetupSystem: " << num_unks << "  unknowns. ";
    std::cerr << "SD_: " << SD_ << " dt_: " << dt_ << " theta_ :" << theta_ << "\n";

    // fill value part of matrices
    Quad5CL<Point3DCL>  u_loc;
    Point3DCL Grad[4];
    Quad5CL<> u_Grad[4], // fuer u grad v_i
              p1[4];
    double det, absdet, h_T;

    LocalP1CL<> p1dummy;
    for (int i= 0; i < 4; ++i) {
        p1dummy[i]= 1.0;
        p1[i].assign( p1dummy);
        p1dummy[i]= 0.0;
    }

    DROPS_FOR_TRIANG_CONST_TETRA( const_cast<const MultiGridCL&>( MG_), lvl, sit) {
        P1DiscCL::GetGradients( Grad, det, *sit);
        absdet= std::fabs( det);
        h_T= std::pow( absdet, 1./3.);

        // save information about the edges and verts of the tetra in Numb
        GetLocalNumbP1NoBnd( Numb, *sit, p1idx);

        // save velocities inside tetra for quadrature in u_loc
        u_loc.assign( *sit, vel);
        for(int i=0; i<4; ++i)
            u_Grad[i]= dot( Grad[i], u_loc);

        /// \todo fixed limit for maxV (maxV_limit), any better idea?
        double maxV = 1.; // scaling of SD parameter
        //const double limit= 1e-3;
//         for(int i= 0; i < Quad5CL<>::NumNodesC; ++i)
//            maxV = std::max( maxV, u_loc[i].norm());
        // const double maxV= 1.; // no scaling
//        const double sd_fact= SD_ * h_T; // / std::max( maxV, limit / h_T);
        for(int i= 0; i < 4; ++i)    // assemble row Numb[i]
            for(int j= 0; j < 4; ++j) {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
               bE( Numb[i], Numb[j])+= P1DiscCL::GetMass(i,j) * absdet
                                       + Quad5CL<>( u_Grad[i]*p1[j]).quad( absdet)*SD_/maxV*h_T;

               // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
               bH( Numb[i], Numb[j])+= Quad5CL<>( u_Grad[j]*p1[i]).quad( absdet)
                                       + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * SD_/maxV*h_T;
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
/*                bE( Numb[i], Numb[j])+= P1DiscCL::GetMass(i,j) * absdet * (1. + sd_fact/std::max(dt_, 1.2*h_T/maxV))
                                        + Quad5CL<>( u_Grad[i]*p1[j]).quad( absdet) * sd_fact * theta_;

                // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                bH( Numb[i], Numb[j])+= Quad5CL<>( u_Grad[j]*p1[i]).quad( absdet) * (1. + sd_fact/std::max(dt_, 1.2*h_T/maxV))
                                      + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * sd_fact * theta_;*/
            }
    }
    bE.Build();
    bH.Build();
    std::cerr << E.num_nonzeros() << " nonzeros in E, "
              << H.num_nonzeros() << " nonzeros in H! " << std::endl;
}

template <class DiscVelSolT>
void TransportP1FunctionCL::DoStep (VectorCL& u, const DiscVelSolT& vel)
{
    SetupSystem( vel, E_, H_);
    L_.clear();
    const int N= 1;
    L_.LinComb( 1./(dt_/N), E_, theta_, H_);

    for (int i= 0; i < N ; ++i) {
	    VectorCL rhs( (1./(dt_/N))*u);
	    if (theta_ != 1.) {
		GMResSolverCL<GSPcCL> gm( gm_);
		VectorCL tmp( rhs.size());
		gm.Solve( E_old, tmp, VectorCL( H_old*u));
		std::cerr << "TransportP1FunctionCL::DoStep rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
		rhs-= (1. - theta_)*tmp;
	    }
	    gm_.Solve( L_, u, VectorCL(E_*rhs));
            if (N != 1) std::cerr << "substep " << i << ": iter: " << gm_.GetIter() << " ";
    }
    std::cerr << "TransportP1FunctionCL::DoStep: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

typedef BndDataCL<Point3DCL> VelBndDataT;

/// \brief P1-discretization and solution of the transport equation on the interface
class SurfactantP1CL
{
  public:
//    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef NoBndDataCL<>                                     BndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;

    IdxDescCL idx, full_idx;
    VecDescCL ic; ///< concentration on the interface
    MatDescCL A,  ///< diffusion matrix
              M,  ///< mass matrix
              C,  ///< convection matrix
              Md; ///< mass matrix with interface-divergence of velocity

  private:
    MatrixCL      L_;              ///< sum of matrices
    MultiGridCL&  MG_;
    double        D_,              ///< diffusion coefficient
                  theta_, dt_, t_; ///< time scheme parameter, time step and time

    BndDataT            Bnd_;
    const VelBndDataT&  Bnd_v_;  ///< Boundary condition for the velocity
    VecDescCL*          v_;      ///< velocity at current time step
    LevelsetP2CL&       lset_;   ///< levelset at current time step
    TransportP1FunctionCL* fulltransport_;

    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;
    double omit_bound_;

  public:
    SurfactantP1CL (MultiGridCL& mg, const VelBndDataT& Bnd_v,
        double theta, double D, VecDescCL* v, LevelsetP2CL& lset,
        double t, double dt, int iter= 1000, double tol= 1e-7, double omit_bound= -1.)
    : idx( P1_FE), full_idx( P1_FE), MG_( mg), D_( D), theta_( theta), dt_( dt), t_( t),
        Bnd_v_( Bnd_v), v_( v), lset_( lset), fulltransport_( 0), gm_( pc_, 100, iter, tol, true),
        omit_bound_( omit_bound)
    {}

    const MultiGridCL& GetMG() const { return MG_; }
    GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }

     /// initialize the interface concentration
    void Init (instat_scalar_fun_ptr);

    /// \remarks call SetupSystem \em before calling SetTimeStep!
    void SetTimeStep( double dt, double theta=-1);

    /// perform one time step
    void DoStep (double new_t);

    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &ic, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myic) const
        { return const_DiscSolCL( &Myic, &Bnd_, &MG_); }
    ///@}

    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the coupling navstokes-levelset. They should not be called by a common user.
    /// Use DoStep() instead.
    ///@{
    VectorCL InitStep ();
    void DoStep (const VectorCL&);
    void CommitStep ();
    void Update ();
    ///@}
};

void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, MultiGridCL& mg, double t)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( idx))
            ic.Data[it->Unknowns( idx)]= icf( it->GetCoord(), t);
    }
}

void SurfactantP1CL::Init (instat_scalar_fun_ptr icf)
{
    P1Init ( icf, ic, MG_, 0.);
}

void SurfactantP1CL::SetTimeStep (double dt, double theta)
{
    dt_= dt;
    if (theta >= 0. && theta <= 1.) theta_= theta;
}

void SurfactantP1CL::Update()
{
    IdxDescCL* cidx= ic.RowIdx;

    M.Data.clear();
    M.SetIdx( cidx, cidx);
    DROPS::SetupInterfaceMassP1( MG_, &M, lset_.Phi);
    std::cerr << "M is set up.\n";
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    DROPS::SetupLBP1( MG_, &A, lset_.Phi, D_);
    std::cerr << "A is set up.\n";
    //C.Data.clear();
    //C.SetIdx( cidx, cidx);
    //DROPS::SetupConvectionP1( MG_, &C, lset_.Phi, make_P2Eval( MG_, Bnd_v_, *v_, t_));
    //std::cerr << "C is set up.\n";
    Md.Data.clear();
    Md.SetIdx( cidx, cidx);
    DROPS::SetupMassDivP1( MG_, &Md, lset_.Phi, make_P2Eval( MG_, Bnd_v_, *v_, t_));
    std::cerr << "Md is set up.\n";
    std::cerr << "SurfactantP1CL::Update: Finished\n";
}

VectorCL SurfactantP1CL::InitStep ()
{
    full_idx.CreateNumbering( idx.TriangLevel(), MG_);
    std::cout << "full NumUnknowns: " << full_idx.NumUnknowns() << std::endl;

    fulltransport_= new TransportP1FunctionCL( MG_, make_P2Eval( MG_, Bnd_v_, *v_, t_), full_idx, dt_,  /*theta=*/theta_, /*SD=*/ ::C.surf_SD, /*iter=*/ 2000, /*tol=*/ 0.1*gm_.GetTol());

    VecDescCL rhs( &idx);
    rhs.Data= (1./dt_)*ic.Data;
    if (theta_ != 1.) {
        GMResSolverCL<GSPcCL> gm( gm_);
        VectorCL tmp( rhs.Data.size());
        gm.Solve( M.Data, tmp, VectorCL( A.Data*ic.Data + Md.Data*ic.Data));
        std::cerr << "SurfactantP1CL::InitStep: rhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
        rhs.Data-= (1. - theta_)*tmp;
    }
    // rhs.Data*= std::sqrt( M.Data.GetDiag()); // scaling
    DROPS::VecDescCL rhsext( &full_idx);
    DROPS::Extend( MG_, rhs, rhsext);
    return rhsext.Data;
}

void SurfactantP1CL::DoStep (const VectorCL& rhsext)
{
    idx.DeleteNumbering( MG_);
    CreateNumbOnInterface( idx.TriangLevel(), idx, MG_, lset_.Phi, omit_bound_);
    std::cout << "new NumUnknowns: " << idx.NumUnknowns() << std::endl;

    VecDescCL transp_rhs( &idx),
              transp_rhsext( &full_idx);
    transp_rhsext.Data= rhsext;
    fulltransport_->DoStep( transp_rhsext.Data, make_P2Eval( MG_, Bnd_v_, *v_, t_)); // XXX new velocity and old velocity in constrt.
    Restrict( MG_, transp_rhsext, transp_rhs);

    ic.SetIdx( &idx);
    Update();
    // transp_rhs.Data/= std::sqrt( M.Data.GetDiag()); // scaling
    L_.LinComb( 1./dt_, M.Data, theta_, A.Data, theta_, Md.Data);
    const VectorCL therhs( M.Data*transp_rhs.Data);
    std::cerr << "Before solve: res = " << norm( L_*ic.Data - therhs) << std::endl;
    gm_.Solve( L_, ic.Data, therhs);
    std::cerr << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void SurfactantP1CL::CommitStep ()
{
    full_idx.DeleteNumbering( MG_);
    delete fulltransport_;
}

void SurfactantP1CL::DoStep (double new_t) // XXX inkonsistenz bei velocity
{
    VectorCL rhs( InitStep());
    t_= new_t;
    DoStep( rhs);
    CommitStep();
}


void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset,
    DROPS::EnsightP2SolOutCL& ensight)
{
    using namespace DROPS;

    LSInit( mg, lset.Phi, &sphere_2move, 0.);

    ensight.DescribeVector( "Velocity", datvel, true);
    ensight.DescribeScalar( "TrueSol", datsol, true);

    BndDataCL<> lsetbnd2( 6);
    DROPS::LevelsetP2CL lset2( mg, lsetbnd2, 0, 0, C.lset_theta, C.lset_SD); // Only for output
    lset2.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset2.Phi.SetIdx( &lset2.idx);
    LSInit( mg, lset2.Phi, &sphere_2move, 0.);
    ensight.DescribeScalar( "Lset", datlset, true);

    const double Vol= lset.GetVolume();
    std::cerr << "droplet volume: " << Vol << std::endl;

    BndDataCL<Point3DCL> Bnd_v( 6, bc, bf);
    IdxDescCL vidx( vecP2_FE);
    vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( mg, &v, Bnd_v, u_func, 0.);

    lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, 0.));
    lset2.SetTimeStep( C.dt);

    SurfactantP1CL timedisc( mg, Bnd_v, C.theta_surf, C.muI, &v, lset, /* t */ 0., C.dt, C.surf_iter, C.surf_tol);

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    InterfaceP1RepairCL ic_repair( mg, lset.Phi, timedisc.ic);
    adap.push_back( &ic_repair);
    LevelsetRepairCL lset2repair( lset2);
    adap.push_back( &lset2repair);

    // Init Interface-Sol
    DROPS::CreateNumbOnInterface( mg.GetLastLevel(), timedisc.idx, mg, lset.Phi, C.surf_omit_bound);
    std::cout << "NumUnknowns: " << timedisc.idx.NumUnknowns() << std::endl;
    timedisc.ic.SetIdx( &timedisc.idx);
    timedisc.Init( &sol0);
    timedisc.Update();

    ensight.putGeom( datgeo, 0);
    ensight.putScalar( datscl,  lset.GetSolution(), 0);
    ensight.putVector( datvel, make_P2Eval( mg, Bnd_v, v, 0.), 0);
    DROPS::IdxDescCL ifacefullidx( P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);
    DROPS::VecDescCL icext( &ifacefullidx);
    DROPS::Extend( mg, timedisc.ic, icext);
    DROPS::NoBndDataCL<> nobnd;
    ensight.putScalar( datsurf, make_P1Eval( mg, nobnd, icext, 0.), 0);
    P1Init ( sol0, icext, mg, 0.);
    ensight.putScalar( datsol, make_P1Eval( mg, nobnd, icext, 0.), 0);
    ensight.putScalar( datlset,  lset2.GetSolution(), 0);

    // timedisc.SetTimeStep( C.dt, C.theta_surf);
//    std::cerr << "L_2-error: " << L2_error( mg, lset.Phi, timedisc.GetSolution(), &sol0t, 0.)
//              << " norm of true solution: " << L2_norm( mg, lset.Phi, &sol0t, 0.)
//              << std::endl;

    for (int step= 1; step <= C.num_steps; ++step) {
        std::cerr << "======================================================== step " << step << ":\n";

        LSInit( mg, lset.Phi, &sphere_2move, step*C.dt);
        timedisc.DoStep( step*C.dt);

        //lset2.DoStep();
//        VectorCL rhs( lset2.Phi.Data.size());
//        lset2.ComputeRhs( rhs);
//        lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, step*C.dt));
//        lset2.SetTimeStep( C.dt);
//        lset2.DoStep( rhs);

        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        if (C.VolCorr) {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        }
        //if (C.RepFreq && step%C.RepFreq==0) { // reparam levelset function
            // lset.ReparamFastMarching( C.RepMethod);
        if (C.ref_freq != 0 && step%C.ref_freq == 0) {
            if (C.ref_freq != 0) {
                adap.UpdateTriang( lset);
                vidx.DeleteNumbering( mg);
                vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
                v.SetIdx( &vidx);
                InitVel( mg, &v, Bnd_v, u_func, step*C.dt);
                LSInit( mg, lset.Phi, &sphere_2move, step*C.dt);
                timedisc.Update();
                ifacefullidx.DeleteNumbering( mg);
                ifacefullidx.CreateNumbering( mg.GetLastLevel(),  mg);
                icext.SetIdx( &ifacefullidx);

                lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, step*C.dt));
                lset2.SetTimeStep( C.dt);
            }
            std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (C.VolCorr) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
        }
        ensight.putGeom( datgeo, step*C.dt);
        ensight.putScalar( datscl,  lset.GetSolution(), step*C.dt);
        Extend( mg, timedisc.ic, icext);
        ensight.putScalar( datsurf, timedisc.GetSolution( icext), step*C.dt);
        P1Init ( sol0t, icext, mg, step*C.dt);
        ensight.putScalar( datsol, timedisc.GetSolution( icext), step*C.dt);
        ensight.putVector( datvel, make_P2Eval( mg, Bnd_v, v, step*C.dt), step*C.dt);
        ensight.putScalar( datlset,  lset2.GetSolution(), step*C.dt);
//        std::cerr << "L_2-error: " << L2_error( mg, lset.Phi, timedisc.GetSolution(), &sol0t, step*C.dt)
//                  << " norm of true solution: " << L2_norm( mg, lset.Phi, &sol0t, step*C.dt)
//                  << std::endl;
    }
    std::cerr << std::endl;
}


int main (int argc, char* argv[])
{
  try {
    std::ifstream param;
    if (argc!=2) {
        std::cerr << "Using default parameter file: surfactant.param\n";
        param.open( "surfactant.param");
    }
    else
        param.open( argv[1]);
    if (!param) {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    // ensight files
    filename= C.EnsDir + "/" + C.EnsCase;
    datgeo=   filename+".geo";
    datscl=   filename+".scl";
    datvel=   filename+".vel";
    datsurf=  filename+".sur";
    datsol =  filename+".sol";
    datlset=  filename+".lset";

    std::cerr << "Setting up interface-PDE:\n";
    DROPS::BrickBuilderCL brick( DROPS::MakePoint3D( -2., -2., -2.),
                                 4.*DROPS::std_basis<3>( 1),
                                 4.*DROPS::std_basis<3>( 2),
                                 4.*DROPS::std_basis<3>( 3),
                                 C.cdiv, C.cdiv, C.cdiv);
    DROPS::MultiGridCL mg( brick);
    DROPS::AdapTriangCL adap( mg, C.ref_width, 0, C.ref_flevel);
    adap.MakeInitialTriang( sphere_2);

    DROPS::LevelsetP2CL lset( mg);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    LinearLSInit( mg, lset.Phi, &sphere_2);
//    lset.Init( &sphere_2);

    DROPS::EnsightP2SolOutCL ensight( mg, &lset.idx);
    ensight.CaseBegin( std::string( C.EnsCase+".case").c_str(), C.num_steps + 1);
    ensight.DescribeGeom( "Geometrie", datgeo, true);
    ensight.DescribeScalar( "Levelset", datscl, true);
//    ensight.DescribeVector( "Velocity",      datvel, true);
    ensight.DescribeScalar( "InterfaceSol", datsurf,  true);

    if (C.TestCase > 0) { // Time dependent tests in Strategy
        Strategy( mg, adap, lset, ensight);
        return 0;
    }

    DROPS::IdxDescCL ifaceidx( P1_FE);
    DROPS::CreateNumbOnInterface( mg.GetLastLevel(), ifaceidx, mg, lset.Phi, C.surf_omit_bound);
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi);
    std::cerr << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi, C.muI);
    DROPS::MatrixCL L;
    L.LinComb( 1.0, A.Data, 1.0, M.Data);
    DROPS::VecDescCL b( &ifaceidx);
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, rhs0);

    DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    DROPS::WriteToFile( b.Data, "rhs_iface.txt", "rhs");

    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, C.surf_iter, C.surf_tol, true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data);
    std::cerr << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    DROPS::WriteToFile( x.Data, "x_iface.txt", "solution");

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);
    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);
    DROPS::NoBndDataCL<> nobnd;
    DiscP1FunT surfsol( &xext, &nobnd, &mg, 0.);
    ensight.putGeom( datgeo, 0);
    ensight.putScalar( datscl,  lset.GetSolution(), 0);
    ensight.putScalar( datsurf, surfsol, 0);
    ensight.Commit();

    double L2_err( L2_error( mg, lset.Phi, surfsol, &sol0));
    std::cerr << "L_2-error: " << L2_err << std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
