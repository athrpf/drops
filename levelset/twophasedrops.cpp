/// \file
/// \brief flow in measurement cell or brick
/// \author Sven Gross, Joerg Grande, Patrick Esser, IGPM

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/coupling.h"
#include "levelset/params.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "poisson/transport2phase.h"
#include "num/stokessolverfactory.h"
#include <fstream>
#include <sstream>


DROPS::ParamMesszelleNsCL C;
// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0


//brickflow.cpp + brick_transp.cpp + brick_ns_adap.cpp
DROPS::SVectorCL<3> InflowBrick( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = p[0]*(2*C.r_inlet-p[0]) / (C.r_inlet*C.r_inlet),
                 z = p[2]*(2*C.r_inlet-p[2]) / (C.r_inlet*C.r_inlet);

    ret[1]= x * z * C.Anstroem * (1-C.inflow_ampl*std::cos(2*M_PI*C.inflow_freq*t));
    return ret;
}

//mzelle_ns_adap.cpp + mzelle_instat.cpp
DROPS::SVectorCL<3> InflowCell( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double s2= C.r_inlet*C.r_inlet,
                 r2= p.norm_sq() - p[C.flow_dir]*p[C.flow_dir];
    ret[C.flow_dir]= -(r2-s2)/s2*C.Anstroem;
    return ret;
}

typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;
const DROPS::BndCondT c_bc[6]= {
    DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC,
    DROPS::OutflowBC, DROPS::OutflowBC, DROPS::OutflowBC
};
const c_bnd_val_fun c_bfun[6]= {0, 0, 0, 0, 0, 0};

double Initialcneg (const DROPS::Point3DCL& , double)
{
    return C.transp_cNeg;
}

double Initialcpos (const DROPS::Point3DCL& , double)
{
    return C.transp_cNeg;
}

namespace DROPS // for Strategy
{

template<class Coeff>
void WriteMatrices (InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, int i)
{
    std::string pfad( "matrices/");
    std::ostringstream n;
    n << pfad << "A" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    std::ofstream s( n.str().c_str());
    s << std::setprecision( 15);
    s << Stokes.A.Data;
    s.close();

    n.str( "");
    n << pfad << "B" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.B.Data;
    s.close();

    n.str( "");
    n << pfad << "M" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.M.Data;
    s.close();

    n.str( "");
    n << pfad << "prA" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.prA.Data;
    s.close();

    n.str( "");
    n << pfad << "prM" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.prM.Data;
    s.close();

    n.str( "");
    n << pfad << "N" << std::setfill( '0') << std::setw( 4) << i << ".txt";
    s.open( n.str().c_str());
    s << Stokes.N.Data;
    s.close();
}

template< class StokesProblemT>
TimeDisc2PhaseCL<StokesProblemT>* CreateTimeDisc(StokesProblemT& Stokes, LevelsetP2CL& lset,
    NSSolverBaseCL<StokesProblemT>* solver, ParamMesszelleNsCL& C, bool usematMG, MGDataCL* matMG)
{
    switch (C.scheme)
    {
        case 1 : 
            return (new LinThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_stab, usematMG, matMG));
        break;
        case 2 :
            return (new RecThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_proj, C.cpl_stab, usematMG, matMG));
        break;
        case 3 :
            return (new ThetaScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.theta, C.nonlinear, C.cpl_proj, C.cpl_stab, usematMG, matMG));
        break;
        case 4 :
            return (new OperatorSplitting2PhaseCL<StokesProblemT, StokesSolverBaseCL>
                        (Stokes, lset, solver->GetStokesSolver(), C.inner_iter, C.inner_tol, C.nonlinear));
        break;
        case 5 :
            return (new CrankNicolsonScheme2PhaseCL<StokesProblemT, NSSolverBaseCL<StokesProblemT> >
                        (Stokes, lset, *solver, C.nonlinear, C.cpl_proj, C.cpl_stab, usematMG, matMG));
        break;
        default : throw DROPSErrCL("Unknown TimeDiscMethod");
    }
}

class FunAsP2EvalCL
{
  private:
    instat_scalar_fun_ptr f_;
  public:
    FunAsP2EvalCL( instat_scalar_fun_ptr f): f_(f)
    {}
    double val(const VertexCL& v) const {return f_(v.GetCoord(), 0.0);}
    double val(const EdgeCL& e, double) const {return f_(GetBaryCenter(e), 0.0);}
};


class EnsightWriterCL
{
  private:
    MultiGridCL& MG_;
    EnsightP2SolOutCL ensight_;
    std::string datgeo_,
                datpr_,
                datvec_,
                datscl_,
                datprx_,
                datxidx_,
                datsf_,
                datc_,
                datct_;
  public:
    EnsightWriterCL (MultiGridCL& MG, const IdxDescCL* idx, const ParamMesszelleNsCL& C);
    ~EnsightWriterCL();

    // To write at time t.
    template<class InstatNSCL>
    void
    WriteAtTime (const InstatNSCL& Stokes, const LevelsetP2CL& lset,
        instat_scalar_fun_ptr sigmap, const TransportP1CL& c, const double t);
};

EnsightWriterCL::EnsightWriterCL (MultiGridCL& MG, const IdxDescCL* idx, const ParamMesszelleNsCL& C)
    : MG_( MG), ensight_( MG, idx)
{
    if (C.EnsCase == "none") return; // no ensight output
    const std::string filename= C.EnsDir + "/" + C.EnsCase;
    datgeo_  = filename+".geo";
    datpr_   = filename+".pr" ;
    datvec_  = filename+".vel";
    datscl_  = filename+".scl";
    datprx_  = filename+".prx";
    datxidx_ = filename+".xidx";
    datsf_   = filename+".sf";
    datc_    = filename+".c";
    datct_   = filename+".ct";
    ensight_.CaseBegin( std::string( C.EnsCase+".case").c_str(), C.num_steps + 1);
    ensight_.DescribeGeom  ( "Messzelle",     datgeo_, true);
    ensight_.DescribeScalar( "Levelset",      datscl_, true);
    ensight_.DescribeScalar( "Pressure",      datpr_,  true);
    ensight_.DescribeVector( "Velocity",      datvec_, true);
    ensight_.DescribeScalar( "Surfaceforce",  datsf_,  true);
    if (C.transp_do)
    {
        ensight_.DescribeScalar( "Concentration", datc_,   true);
        ensight_.DescribeScalar( "TransConc",     datct_, true);
    }
}

EnsightWriterCL::~EnsightWriterCL()
{
    if (C.EnsCase == "none") return;
    ensight_.CaseEnd();
}

template<class InstatNSCL>
void
EnsightWriterCL::WriteAtTime (const InstatNSCL& Stokes, const LevelsetP2CL& lset,
    instat_scalar_fun_ptr sigmap, const TransportP1CL& c, const double t)
{
    if (C.EnsCase == "none") return;
    ensight_.putGeom( datgeo_, t);
    ensight_.putVector( datvec_, Stokes.GetVelSolution(), t);
    ensight_.putScalar( datpr_,  Stokes.GetPrSolution(), t);
    ensight_.putScalar( datscl_, lset.GetSolution(), t);
    FunAsP2EvalCL sf( sigmap);
    ensight_.putScalar( datsf_,  sf, t);
    if (C.transp_do)
    {
        ensight_.putScalar( datc_,  c.GetSolution(), t);
        ensight_.putScalar( datct_, c.GetSolution( c.ct), t);
    }
    ensight_.Commit();
    if (Stokes.UsesXFEM()) {
        std::string datprxnow( datprx_);
        ensight_.AppendTimecode( datprxnow);
        std::ofstream fff( datprxnow.c_str());
        fff.precision( 16);
        size_t num_prx= Stokes.pr_idx.NumUnknowns - Stokes.GetXidx().GetNumUnknownsP1();
        out( fff, VectorCL( Stokes.p.Data[std::slice( Stokes.GetXidx().GetNumUnknownsP1(), num_prx, 1)]));

        std::string datxidxnow( datxidx_);
        ensight_.AppendTimecode( datxidxnow);
        std::ofstream fxidx( datxidxnow.c_str());
        for (Uint i=0; i< Stokes.GetXidx().Xidx.size(); ++i)
            fxidx << Stokes.GetXidx().Xidx[i] << '\n';
    }
}

template<class Coeff>
void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, AdapTriangCL& adap)
// flow control
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    MultiGridCL& MG= Stokes.GetMG();
    sigma= Stokes.GetCoeff().SurfTens;
    eps= C.st_jumpWidth;    lambda= C.st_relPos;    sigma_dirt_fac= C.st_red;
    instat_scalar_fun_ptr sigmap  = 0;
    instat_vector_fun_ptr gsigmap = 0;
    if (C.st_var)
    {
        sigmap  = &sigma_step;
        gsigmap = &gsigma_step;
    }
    else
    {
        sigmap  = &sigmaf;
        gsigmap = &gsigma;
    }
    LevelsetP2CL lset( MG, sigmap, gsigmap, C.lset_theta, C.lset_SD,
        -1, C.lset_iter, C.lset_tol, C.CurvDiff);

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset);
    adap.push_back( &prrepair);

    IdxDescCL* lidx= &lset.idx;
    IdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    if (C.st_var)
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);

    EnsightWriterCL writer( MG, lidx, C);

    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, ZeroVel);
    switch (C.IniCond)
    {
      case  1: //flow without droplet
          lset.Init( &One);
      break;
      case -1: // read from file
      {
        ReadEnsightP2SolCL reader( MG);
        reader.ReadScalar( C.IniData+".scl", lset.Phi, lset.GetBndData());
        reader.ReadVector( C.IniData+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset, /*NumberingChanged*/ false);
        Stokes.p.SetIdx( pidx);
        // reads the P1-part of the pressure
        reader.ReadScalar( C.IniData+".pr",  Stokes.p, Stokes.GetBndData().Pr);
        // reads the P1X-part of the pressure
        if (Stokes.UsesXFEM()) {
            std::ifstream fff  ( (C.IniData+".prx").c_str());
            std::ifstream fxidx( (C.IniData+".xidx").c_str());
            if (fff && fxidx) {
                size_t NumP1XUnknowns;
                fff >> NumP1XUnknowns;
                if (NumP1XUnknowns != (pidx->NumUnknowns - Stokes.GetXidx().GetNumUnknownsP1()))
                    throw (DROPSErrCL("error while reading P1X unknowns"));
                for (Uint i=Stokes.GetXidx().GetNumUnknownsP1(); i < pidx->NumUnknowns; ++i)
                    fff >> Stokes.p.Data[i];
                for (Uint k=0; k<Stokes.GetXidx().GetNumUnknownsP1(); ++k)
                    fxidx >> Stokes.GetXidx().Xidx[k];
            }
        }
      } break;
      default:
        lset.Init( EllipsoidCL::DistanceFct);
    }
    MG.SizeInfo( std::cerr);
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";

    cBndDataCL Bnd_c( 6, c_bc, c_bfun);
    double D[2] = {C.transp_cPos, C.transp_cNeg};
    TransportP1CL c( MG, Bnd_c, Stokes.GetBndData().Vel, C.transp_theta, D, C.transp_H, &Stokes.v, lset,
        /*t*/ 0., C.dt, C.transp_iter, C.transp_tol);
    TransportRepairCL transprepair(c, MG);
    if (C.transp_do)
    {
        adap.push_back(&transprepair);
        IdxDescCL* cidx= &c.idx;
        c.CreateNumbering( MG.GetLastLevel(), cidx);
        c.ct.SetIdx( cidx);
        if (C.IniCond != -1)
            c.Init( &Initialcneg, &Initialcpos);
        else
        {
            ReadEnsightP2SolCL reader( MG);
            reader.ReadScalar( C.IniData+".ct", c.ct, c.GetBndData());
        }
        c.Update();
        std::cerr << c.c.Data.size() << " concentration unknowns,\n";
    }

    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    writer.WriteAtTime( Stokes, lset, sigmap, c, 0.);

    // Stokes-Solver
    StokesSolverFactoryCL<StokesProblemT, ParamMesszelleCL> stokessolverfactory(Stokes, C);
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver();

    // Navier-Stokes-Solver
    NSSolverBaseCL<StokesProblemT>* navstokessolver = 0;
    if (C.nonlinear==0.0)
        navstokessolver = new NSSolverBaseCL<StokesProblemT>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<StokesProblemT>(Stokes, *stokessolver, C.ns_iter, C.ns_tol, C.ns_red);

    // Time discretisation + coupling
    MGDataCL& velMG = stokessolverfactory.GetVelMG();
    bool mgused = stokessolverfactory.MGUsed();
    TimeDisc2PhaseCL<StokesProblemT>* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, C, mgused, &velMG);
    timedisc->SetTimeStep( C.dt);

    stokessolverfactory.SetMatrixA( timedisc->GetUpperLeftBlock());
    bool second = false;
    std::ofstream infofile((C.EnsCase+".info").c_str());
    double lsetmaxGradPhi, lsetminGradPhi;
    std::ofstream out("lset.info");
    for (int step= 1; step<=C.num_steps; ++step)
    {
        std::cerr << "======================================================== step " << step << ":\n";

        IFInfo.Update( lset, Stokes.GetVelSolution());
        IFInfo.Write(infofile);
        timedisc->DoStep( C.cpl_iter);
        if (C.transp_do) c.DoStep( step*C.dt);

//        WriteMatrices( Stokes, step);
        std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
        out << lsetmaxGradPhi << '\t' << lsetminGradPhi << std::endl;

        if (C.RepFreq && step%C.RepFreq==0) { // reparam levelset function
            if (C.VolCorr) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
            lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
            if (lsetmaxGradPhi > C.MaxGrad || lsetminGradPhi < C.MinGrad) {
                std::cerr << "before reparametrisation: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n'; 
                lset.ReparamFastMarching( C.RepMethod);
                lset.GetMaxMinGradPhi( lsetmaxGradPhi, lsetminGradPhi);
                std::cerr << "after  reparametrisation: minGradPhi " << lsetminGradPhi << "\tmaxGradPhi " << lsetmaxGradPhi << '\n'; 
            }
        }
        if (step%C.ref_freq == 0) { // grid modification
            if (C.VolCorr) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cerr << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
            adap.UpdateTriang( lset);
            if (adap.WasModified()) {
                timedisc->Update();
                if (C.transp_do) c.Update();
            }
            if (C.serialization_file != "none" && C.ref_freq != 0) {
                std::stringstream filename;
                filename << C.serialization_file;
                if (second) filename << "0";
                second = !second;
                MGSerializationCL ser( MG, filename.str().c_str());
                ser.WriteMG();
            }
        }
        if (C.VolCorr) {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cerr << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cerr << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        }
        writer.WriteAtTime( Stokes, lset, sigmap, c, step*C.dt);
    }
    std::cerr << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
}

} // end of namespace DROPS

void CreateGeom (DROPS::MultiGridCL* &mgp, DROPS::StokesBndDataCL* &bnddata)
{
    if (C.GeomType == 0) {
        std::ifstream meshfile( C.meshfile.c_str());
        if (!meshfile)
            throw DROPS::DROPSErrCL ("error while opening mesh file\n");

        DROPS::ReadMeshBuilderCL builder (meshfile);
        if (C.deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( builder);
        else {
            DROPS::FileBuilderCL filebuilder( C.deserialization_file, &builder);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        const DROPS::BoundaryCL& bnd= mgp->GetBnd();
        const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

        DROPS::BndCondT* bc = new DROPS::BndCondT[num_bnd];
        DROPS::StokesVelBndDataCL::bnd_val_fun* bnd_fun = new  DROPS::StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= builder.GetBC( i))==DROPS::DirBC ? &InflowCell : &DROPS::ZeroVel;
            std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
        }
        bnddata = new DROPS::StokesBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
    }
    if (C.GeomType == 1) {
        int nx, ny, nz;
        double dx, dy, dz;
        std::string mesh( C.meshfile), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
        if (!brick_info)
            DROPS::DROPSErrCL("error while reading geometry information: " + mesh);
        C.r_inlet= dx/2;
        DROPS::Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;
        DROPS::BrickBuilderCL builder ( orig, px, py, pz, nx, ny, nz);
        if (C.deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( builder);
        else {
            DROPS::FileBuilderCL filebuilder( C.deserialization_file, &builder);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        DROPS::BndCondT bc[6]= { DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC };
        DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[6]=
            { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel };
        switch (C.bnd_type)
        {
            case 1 : // hom. Dirichlet bnd-data
              break;
            case 2 : // inhom. Dirichlet bnd-data for in-/outflow
            {
                bc[2]= bc[3]= DROPS::DirBC;
                bfun[2]= bfun[3]= &InflowBrick;
            } break;
            case 3 : // tube/canal
            {
                bc[3]= DROPS::DirBC;
                bc[2]= DROPS::NatBC; //Rohr
                //bc[2]= bc[4]= bc[5]= DROPS::NatBC;          //Kanal
                bfun[2]= &DROPS::ZeroVel;
                //bfun[2]=bfun[4]=bfun[5]= &DROPS::ZeroVel;   //Kanal
                bfun[2]= 
                bfun[3]= &InflowBrick;
            } break;
            default: throw DROPS::DROPSErrCL("Unknown boundary data type");
        }
        bnddata = new DROPS::StokesBndDataCL(6, bc, bfun);
    }
}

int main (int argc, char** argv)
{
  try
  {
    std::ifstream param;
    if (argc!=2)
    {
        std::cerr << "Using default parameter file: risingdroplet.param\n";
        param.open( "risingdroplet.param");
    }
    else
        param.open( argv[1]);
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    DROPS::MultiGridCL* mg= 0;
    DROPS::StokesBndDataCL* bnddata= 0;

    CreateGeom(mg, bnddata);
    EllipsoidCL::Init( C.Mitte, C.Radius);

    DROPS::AdapTriangCL adap( *mg, C.ref_width, 0, C.ref_flevel);

    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (C.deserialization_file == "none")
        adap.MakeInitialTriang( EllipsoidCL::DistanceFct);

    std::cerr << DROPS::SanityMGOutCL(*mg) << std::endl;
    MyStokesCL prob(*mg, ZeroFlowCL(C), *bnddata, C.XFEMStab<0 ? DROPS::P1_FE : DROPS::P1X_FE, C.XFEMStab);

    Strategy( prob, adap);    // do all the stuff

    delete mg;
    delete bnddata;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
