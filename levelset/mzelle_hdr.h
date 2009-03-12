//**************************************************************************
// File:    mzelle_hdr.h                                                   *
// Content: header for all main programs simulating the measuring cell     *
// Author:  Sven Gross, Joerg Grande, IGPM RWTH Aachen                     *
//**************************************************************************

#ifndef DROPS_MZELLE_HDR_H
#define DROPS_MZELLE_HDR_H

#include "geom/multigrid.h"
#include "levelset/params.h"
#include "levelset/levelset.h"
#include "num/discretize.h"
#include "stokes/instatstokes2phase.h"

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

/// \brief Parameter class describing a zero flow
class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleNsCL& C)
      : rho( DROPS::JumpCL( C.rhoD, C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muD,  C.muF),   DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( const DROPS::ParamMesszelleNsCL& C)
      : rho( DROPS::JumpCL( 1., C.rhoF/C.rhoD ), DROPS::H_sm, C.sm_eps),
        mu ( DROPS::JumpCL( 1., C.muF/C.muD),    DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma/C.rhoD), g( C.g)    {}
};

class EllipsoidCL
{
  private:
    static DROPS::Point3DCL Mitte_;
    static DROPS::Point3DCL Radius_;

  public:
    EllipsoidCL( const DROPS::Point3DCL& Mitte, const DROPS::Point3DCL& Radius)
    { Init( Mitte, Radius); }
    static void Init( const DROPS::Point3DCL& Mitte, const DROPS::Point3DCL& Radius)
    { Mitte_= Mitte;    Radius_= Radius; }
    static double DistanceFct( const DROPS::Point3DCL& p)
    {
        DROPS::Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() { return 4./3.*M_PI*Radius_[0]*Radius_[1]*Radius_[2]; }
};

DROPS::Point3DCL EllipsoidCL::Mitte_;
DROPS::Point3DCL EllipsoidCL::Radius_;

// collision setting (rising butanol droplet in water)
//  RadDrop1 =  1.50e-3  1.500e-3  1.50e-3
//  PosDrop1 =  6.00e-3  3.000e-3  6.00e-3
//  RadDrop2 =  0.75e-3  0.750e-3  0.75e-3
//  PosDrop2 =  6.00e-3  5.625e-3  6.00e-3
//  MeshFile =  12e-3x30e-3x12e-3@4x10x4

class TwoEllipsoidCL
{
  private:
    static DROPS::Point3DCL Mitte1_,  Mitte2_;
    static DROPS::Point3DCL Radius1_, Radius2_;

  public:
    TwoEllipsoidCL( const DROPS::Point3DCL& Mitte1, const DROPS::Point3DCL& Radius1,
                    const DROPS::Point3DCL& Mitte2, const DROPS::Point3DCL& Radius2)
    { Init( Mitte1, Radius1, Mitte2, Radius2); }
    static void Init( const DROPS::Point3DCL& Mitte1, const DROPS::Point3DCL& Radius1,
                      const DROPS::Point3DCL& Mitte2, const DROPS::Point3DCL& Radius2)
    { Mitte1_= Mitte1;    Radius1_= Radius1;
      Mitte2_= Mitte2;    Radius2_= Radius2;}
    static double DistanceFct( const DROPS::Point3DCL& p)
    {
        DROPS::Point3DCL d1= p - Mitte1_;
        const double avgRad1= cbrt(Radius1_[0]*Radius1_[1]*Radius1_[2]);
        d1/= Radius1_;
        DROPS::Point3DCL d2= p - Mitte2_;
        const double avgRad2= cbrt(Radius2_[0]*Radius2_[1]*Radius2_[2]);
        d2/= Radius2_;
        return std::min(std::abs( avgRad1)*d1.norm() - avgRad1, std::abs( avgRad2)*d2.norm() - avgRad2);
    }
    static double GetVolume() { return 4./3.*M_PI*Radius1_[0]*Radius1_[1]*Radius1_[2]
        + 4./3.*M_PI*Radius2_[0]*Radius2_[1]*Radius2_[2]; }
};

DROPS::Point3DCL TwoEllipsoidCL::Mitte1_;
DROPS::Point3DCL TwoEllipsoidCL::Radius1_;
DROPS::Point3DCL TwoEllipsoidCL::Mitte2_;
DROPS::Point3DCL TwoEllipsoidCL::Radius2_;

class InterfaceInfoCL
{
  public:
    DROPS::Point3DCL bary, vel, min, max;
    double maxGrad, Vol, h_min, h_max;

    template<class DiscVelSolT>
    void Update (const DROPS::LevelsetP2CL& ls, const DiscVelSolT& u) {
        ls.GetInfo( maxGrad, Vol, bary, vel, u, min, max);
        std::pair<double, double> h= h_interface( ls.GetMG().GetTriangEdgeBegin( ls.Phi.RowIdx->TriangLevel()), ls.GetMG().GetTriangEdgeEnd( ls.Phi.RowIdx->TriangLevel()), ls.Phi);
        h_min= h.first; h_max= h.second;
    }
    void WriteHeader(std::ofstream& file) {
        file << "# time maxGradPhi volume bary_drop min_drop max_drop vel_drop h_min h_max" << std::endl;
    }
    void Write (double time, std::ofstream& file) {
        file << time << " " << maxGrad << " " << Vol << " " << bary << " " << min << " " << max << " " << vel << " " << h_min << " " << h_max << std::endl;
    }
} IFInfo;

double eps=5e-4, // Sprungbreite
    lambda=1.5, // Position des Sprungs zwischen Oberkante (lambda=0) und Schwerpunkt (lambda=1)
    sigma_dirt_fac= 0.8; // gesenkte OFspannung durch Verunreinigungen im unteren Teil des Tropfens
double sigma;

double sm_step(const DROPS::Point3DCL& p)
{
    double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
        y= p[1] - y_mid;
    if (y > eps) return sigma;
    if (y < -eps) return sigma_dirt_fac*sigma;
    const double z=y/eps*M_PI/2.;
    return sigma_dirt_fac*sigma + (sigma - sigma_dirt_fac*sigma) * (std::sin(z)+1)/2;
}

DROPS::Point3DCL grad_sm_step (const DROPS::Point3DCL& p)
{
    double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
        y= p[1] - y_mid;
    DROPS::Point3DCL ret;
    if (y > eps) return ret;
    if (y < -eps) return ret;
    const double z=y/eps*M_PI/2.;
    ret[1]= (sigma - sigma_dirt_fac*sigma) * std::cos(z)/eps*M_PI/4;
    return ret;
}

double lin(const DROPS::Point3DCL& p)
{
    const double y_top= IFInfo.max[1],
                 y_bot= IFInfo.bary[1],
                 y_slope= sigma*(1 - sigma_dirt_fac)/(y_top - y_bot);
    return sigma + (p[1] - y_top)*y_slope;
}

DROPS::Point3DCL grad_lin (const DROPS::Point3DCL&)
{
    const double y_top= IFInfo.max[1],
                 y_bot= IFInfo.bary[1],
                 y_slope= sigma*(1 - sigma_dirt_fac)/(y_top - y_bot);
    DROPS::Point3DCL ret;
    ret[1]= y_slope;
    return ret;
}

double sigmaf (const DROPS::Point3DCL&, double) { return sigma; }
DROPS::Point3DCL gsigma (const DROPS::Point3DCL&, double) { return DROPS::Point3DCL(); }

double sigma_step(const DROPS::Point3DCL& p, double) { return sm_step( p); }
DROPS::Point3DCL gsigma_step (const DROPS::Point3DCL& p, double) { return grad_sm_step( p); }

namespace DROPS{
/// \brief Timedepending zero function as a helper function
SVectorCL<3> Null( const Point3DCL&, double) { return DROPS::SVectorCL<3>(0.); }
double       One ( const Point3DCL&)         { return 1.; }
}

/// \brief Create geometry of a Mzelle of a brick
void CreateGeom (DROPS::MultiGridCL* &mgp, DROPS::StokesBndDataCL* &bnddata,
                 DROPS::instat_vector_fun_ptr inflow,
                 const std::string& meshfile_name,
                 int GeomType, int bnd_type,
                 const std::string& deserialization_file, double& r_inlet)
{
#ifdef _PAR
    DROPS::ParMultiGridCL::InstancePtr();
#endif
    if (GeomType == 0) {
        std::ifstream meshfile( meshfile_name.c_str());
        if (!meshfile)
            throw DROPS::DROPSErrCL ("error while opening mesh file\n");

        DROPS::ReadMeshBuilderCL *mgb= 0;       // builder of the multigrid

        // read geometry information from a file and create the multigrid
        IF_MASTER
            mgb = new DROPS::ReadMeshBuilderCL( meshfile );
        IF_NOT_MASTER
            mgb = new DROPS::EmptyReadMeshBuilderCL( meshfile );
        // Create the multigrid
        if (deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( *mgb);
        else {
            DROPS::FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        const DROPS::BoundaryCL& bnd= mgp->GetBnd();
        const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

        DROPS::BndCondT* bc = new DROPS::BndCondT[num_bnd];
        DROPS::StokesVelBndDataCL::bnd_val_fun* bnd_fun = new DROPS::StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= mgb->GetBC( i))==DROPS::DirBC ? inflow : &DROPS::ZeroVel;
            std::cerr << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cerr);
        }
        bnddata = new DROPS::StokesBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
        delete   mgb;
    }
    if (GeomType == 1) {
        int nx, ny, nz;
        double dx, dy, dz;
        std::string mesh( meshfile_name), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
        if (!brick_info)
            throw DROPS::DROPSErrCL("error while reading geometry information: " + mesh);
        r_inlet= dx/2;
        DROPS::Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;

        DROPS::BrickBuilderCL *mgb = 0;
        IF_MASTER
            mgb = new DROPS::BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
        IF_NOT_MASTER
            mgb = new DROPS::EmptyBrickBuilderCL(orig, px, py, pz);

        if (deserialization_file == "none")
            mgp= new DROPS::MultiGridCL( *mgb);
        else {
            DROPS::FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new DROPS::MultiGridCL( filebuilder);
        }
        DROPS::BndCondT bc[6]= { DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC, DROPS::Dir0BC };
        DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[6]=
            { &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel, &DROPS::ZeroVel };
        switch (bnd_type)
        {
            case 1 : // hom. Dirichlet bnd-data
              break;
            case 2 : // inhom. Dirichlet bnd-data for in-/outflow
            {
                bc[2]= bc[3]= DROPS::DirBC;
                bfun[2]= bfun[3]= inflow;
            } break;
            case 3 : // tube/canal
            {
                bc[3]= DROPS::DirBC;
                bc[2]= DROPS::NatBC; //Rohr
                //bc[2]= bc[4]= bc[5]= DROPS::NatBC;          //Kanal
                bfun[2]= &DROPS::ZeroVel;
                //bfun[2]=bfun[4]=bfun[5]= &DROPS::ZeroVel;   //Kanal
                bfun[3]= inflow;
            } break;
            default: throw DROPS::DROPSErrCL("Unknown boundary data type");
        }
        bnddata = new DROPS::StokesBndDataCL(6, bc, bfun);
        delete mgb;
    }
}

/// \brief Display a detailed list of unknowns
template <typename StokesT, typename LevelsetT>
  void DisplayUnks(const StokesT& Stokes, const LevelsetT& levelset, __UNUSED__ const DROPS::MultiGridCL& MG)
/** This functions write information about unknowns on the display. These
    informations are for the level-set-, pressure- and velocity-DOF:
    - global DOF
    - accumulated DOF
    - max and min DOF on a single processor (and the ratio)
    - max and min number of distributed DOF on a processor (and the ratio to the remaining DOF)
*/
{
#ifndef _PAR
    std::cerr << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cerr << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cerr << levelset.Phi.Data.size() << " levelset unknowns.\n";
#else
    using namespace DROPS;
    const MLIdxDescCL* vidx = &Stokes.vel_idx,
                     * pidx = &Stokes.pr_idx;
    const IdxDescCL*   lidx = &levelset.idx;
    const ExchangeCL& ExV = Stokes.vel_idx.GetEx(),
                    & ExP = Stokes.pr_idx.GetEx(),
                    & ExL = levelset.idx.GetEx();

    // local number on unknowns
    Ulint Psize      = pidx->NumUnknowns();
    Ulint Vsize      = vidx->NumUnknowns();
    Ulint Lsize      = lidx->NumUnknowns();

    // global number of unknowns
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);
    Ulint GLsize     = lidx->GetGlobalNumUnknowns(MG);

    // accumulated size of unknwons
    Ulint Psize_acc = ProcCL::GlobalSum(Psize);
    Ulint Vsize_acc = ProcCL::GlobalSum(Vsize);
    Ulint Lsize_acc = ProcCL::GlobalSum(Lsize);

    // maximal and minimal number of unknowns
    Ulint P_min= ProcCL::GlobalMin(Psize); Ulint P_max= ProcCL::GlobalMax(Psize);
    Ulint V_min= ProcCL::GlobalMin(Vsize); Ulint V_max= ProcCL::GlobalMax(Vsize);
    Ulint L_min= ProcCL::GlobalMin(Lsize); Ulint L_max= ProcCL::GlobalMax(Lsize);

    // ratios between maximal number of unknowns/proc and minimal number
    double P_ratio   = (double)P_max/(double)P_min;
    double V_ratio   = (double)V_max/(double)V_min;
    double L_ratio   = (double)L_max/(double)L_min;

    // number on boundaries
    Ulint P_accmax= ProcCL::GlobalMax(ExP.AccDistIndex.size()), P_accmin= ProcCL::GlobalMin(ExP.AccDistIndex.size());
    Ulint V_accmax= ProcCL::GlobalMax(ExV.AccDistIndex.size()), V_accmin= ProcCL::GlobalMin(ExV.AccDistIndex.size());
    Ulint L_accmax= ProcCL::GlobalMax(ExL.AccDistIndex.size()), L_accmin= ProcCL::GlobalMin(ExL.AccDistIndex.size());

    // ratio of these unknowns
    double P_accratio= (double)P_accmax / (double)P_accmin;
    double V_accratio= (double)V_accmax / (double)V_accmin;
    double L_accratio= (double)L_accmax / (double)L_accmin;

    // output on screen
    if (ProcCL::IamMaster()){
        std::cerr << "  + Number of DOF\n        "
                  << std::setw(10)<<"global"<<std::setw(10)<<"accum"<<std::setw(10)
                  << "max"<<std::setw(10)<<"min"<<std::setw(10)<<"ratio"<<"  |  "
                  << std::setw(10)<<"max_acc" <<std::setw(10)<<"min_acc"<<std::setw(10)<<"ratio_acc"<<std::endl;

        std::cerr << "    "<<"pr  "
                  << std::setw(10)<<GPsize<<std::setw(10)<<Psize_acc<<std::setw(10)<<P_max
                  << std::setw(10)<<P_min<< std::setw(10)<<P_ratio<<"  |  "
                  << std::setw(10)<<P_accmax<<std::setw(10)<<P_accmin<<std::setw(10)<<P_accratio<<std::endl;

        std::cerr << "    "<<"vel "
                  << std::setw(10)<<GVsize<<std::setw(10)<<Vsize_acc<<std::setw(10)<<V_max
                  << std::setw(10)<<V_min<< std::setw(10)<<V_ratio<<"  |  "
                  << std::setw(10)<<V_accmax<<std::setw(10)<<V_accmin<<std::setw(10)<<V_accratio<<std::endl;

        std::cerr << "    "<<"scl "
                  << std::setw(10)<<GLsize<<std::setw(10)<<Lsize_acc<<std::setw(10)<<L_max
                  << std::setw(10)<<L_min<< std::setw(10)<<L_ratio<<"  |  "
                  << std::setw(10)<<L_accmax<<std::setw(10)<<L_accmin<<std::setw(10)<<L_accratio<<std::endl;

        std::cerr << std::endl;
    }
#endif
}

void DisplayDetailedGeom(DROPS::MultiGridCL& mg)
{
#ifndef _PAR
    mg.SizeInfo( std::cerr);
#else
    const DROPS::Uint level=mg.GetLastLevel();
    DROPS::Uint *numTetrasAllProc=0;
    DROPS::Uint *numFacesAllProc=0;
    DROPS::Uint *numDistFaceAllProc=0;
    if (DROPS::ProcCL::IamMaster()){
        numTetrasAllProc  = new DROPS::Uint[DROPS::ProcCL::Size()];
        numFacesAllProc   = new DROPS::Uint[DROPS::ProcCL::Size()];
        numDistFaceAllProc= new DROPS::Uint[DROPS::ProcCL::Size()];
    }
    // Gather information about distribution on master processor
    DROPS::ProcCL::Gather(mg.GetNumTriangTetra(level),      numTetrasAllProc,   DROPS::ProcCL::Master());
    DROPS::ProcCL::Gather(mg.GetNumTriangFace(level),       numFacesAllProc,    DROPS::ProcCL::Master());
    DROPS::ProcCL::Gather(mg.GetNumDistributedFaces(level), numDistFaceAllProc, DROPS::ProcCL::Master());

    // Display information
    if (DROPS::ProcCL::IamMaster()){
        double ratioTetra       =  (double)*std::max_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size())
                                  /(double)*std::min_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size());
        DROPS::Uint allTetra    =  std::accumulate(numTetrasAllProc, numTetrasAllProc+DROPS::ProcCL::Size(), 0),
                    allFace     =  std::accumulate(numFacesAllProc, numFacesAllProc+DROPS::ProcCL::Size(), 0),
                    allDistFace =  std::accumulate(numDistFaceAllProc, numDistFaceAllProc+DROPS::ProcCL::Size(), 0);
        double      *ratioDistFace=new double[DROPS::ProcCL::Size()];

        // global information
        std::cerr << "Detailed information about the parallel multigrid:\n"
                  << "#(master tetras on finest level):    "<<allTetra<<'\n'
                  << "#(all Faces on finest level):        "<<allFace<<'\n'
                  << "#(distributed Faces on fines level): "<<allDistFace<<'\n';

        // local information for all processors
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            ratioDistFace[i]= ((double)numDistFaceAllProc[i]/(double)numFacesAllProc[i]*100.);

        double maxRatio= *std::max_element(ratioDistFace, ratioDistFace+DROPS::ProcCL::Size());
        std::cerr << "Ratio between max/min Tetra: "<<ratioTetra
                  <<" max Ratio DistFace/AllFace: "<<maxRatio<<std::endl;

        std::cerr << std::setw(6)  <<  "Proc"
                  << std::setw(8)  << "#Tetra"
                  << std::setw(8)  << "#Faces"
                  << std::setw(12) << "#DistFaces"
                  << std::setw(12) << "%DistFaces"
                  << '\n';
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            std::cerr << std::setw(6)  << i
                      << std::setw(8)  << numTetrasAllProc[i]
                      << std::setw(8)  << numFacesAllProc[i]
                      << std::setw(12) << numDistFaceAllProc[i]
                      << std::setw(12) << ratioDistFace[i] << std::endl;

        // free memory
        if (numTetrasAllProc)   delete[] numTetrasAllProc;
        if (numFacesAllProc)    delete[] numFacesAllProc;
        if (numDistFaceAllProc) delete[] numDistFaceAllProc;
        if (ratioDistFace)      delete[] ratioDistFace;
    }
#endif
}

namespace DROPS{
/// \brief Class for serializing a two-phase flow problem, i.e., storing
///    the multigrid and the numerical data
/// \todo  Storing of transport data!
template <typename StokesT>
class TwoPhaseStoreCL
{
  private:
    MultiGridCL&        mg_;
    const StokesT&      Stokes_;
    const LevelsetP2CL& lset_;
    std::string         path_;
    Uint                numRecoverySteps_;
    Uint                recoveryStep_;
    IdxDescCL           p1idx_;             ///< used to split P1X into two P1 functions
    VecDescCL           pneg_,              ///< one part of P1X
                        ppos_;              ///< other part of P1X


    /// \brief Write a numerical data, stored in v, in a file, named filename
    template <typename VecDescT>
    void WriteFEToFile(VecDescT& v, std::string filename)
    {
#ifdef _PAR
        ProcCL::AppendProcNum( filename);
#endif
        std::ofstream file( filename.c_str());
        if (!file) throw DROPSErrCL("TwoPhaseStoreCL::WriteFEToFile: Cannot open file from writing");
        v.Write( mg_, file);
    }

  public:
      /// \brief Construct a class for storing a two-phase flow problem in files
      /** This class generates multiple files, all with prefix path, for storing
       *  the geometric as well as the numerical data.
       *  \param recoverySteps number of backup steps before overwriting files
       *  */
    TwoPhaseStoreCL(MultiGridCL& mg, const StokesT& Stokes, const LevelsetP2CL& lset,
                    const std::string& path, Uint recoverySteps=2)
      : mg_(mg), Stokes_(Stokes), lset_(lset), path_(path), numRecoverySteps_(recoverySteps),
        recoveryStep_(0), p1idx_( P1_FE), pneg_( &p1idx_), ppos_( &p1idx_) {}

    /// \brief Write all information in a file
    void Write()
    {
        // Create filename
        std::stringstream filename;
        filename << path_ << ((recoveryStep_++)%numRecoverySteps_);

        // write multigrid
        MGSerializationCL ser( mg_, filename.str().c_str());
        ser.WriteMG();

#ifdef _PAR
        // write numerical data
        // Since serial DROPS can read ensight files, this is omitted for serial DROPS.
        WriteFEToFile(Stokes_.v, filename.str() + "velocity");
        WriteFEToFile(lset_.Phi, filename.str() + "levelset");
        if ( Stokes_.UsesXFEM()){
            if ( p1idx_.NumUnknowns()!=0)
                p1idx_.DeleteNumbering( mg_);
            p1idx_.CreateNumbering( Stokes_.p.RowIdx->TriangLevel(), mg_); // Create a P1 describer for splitting P1X function
            P1XtoP1 ( *Stokes_.p.RowIdx, Stokes_.p.Data, p1idx_, ppos_.Data, pneg_.Data, lset_.Phi, mg_); // also resizes data of ppos_ and pneg_
            WriteFEToFile(pneg_, filename.str() + "pressureNeg");
            WriteFEToFile(ppos_, filename.str() + "pressurePos");
        }
        else{
            WriteFEToFile(Stokes_.p, filename.str() + "pressure");
        }
#endif
    }
};

/// Read a serialized finite element function from a file
/// \param VecDescT VecDescBaseCL corresponding to the finite element function
/// \pre CreateNumering of v.RowIdx must have been called
template <typename VecDescT>
void ReadFEFromFile(VecDescT& v, const MultiGridCL& mg, std::string filename)
{
    std::cout << "Read FE "<<filename<<std::endl;
#ifdef _PAR
    ProcCL::AppendProcNum( filename);
#endif
    std::ifstream file( filename.c_str());
    if (!file) throw DROPSErrCL("ReadFEFromFile: Cannot open file");
    v.Read( mg, file);
}
}   // end of namespace DROPS

#endif
