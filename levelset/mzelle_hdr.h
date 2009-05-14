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
#include "poisson/transport2phase.h"

namespace DROPS
{

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0

/// \brief Parameter class describing a zero flow
class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static Point3DCL f(const Point3DCL&, double)
        { Point3DCL ret(0.0); return ret; }
    const SmoothedJumpCL rho, mu;
    const double SurfTens;
    const Point3DCL g;

    ZeroFlowCL( const ParamMesszelleNsCL& C)
      : rho( JumpCL( C.rhoD, C.rhoF ), H_sm, C.sm_eps),
        mu(  JumpCL( C.muD,  C.muF),   H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static Point3DCL f(const Point3DCL&, double)
        { Point3DCL ret(0.0); return ret; }
    const SmoothedJumpCL rho, mu;
    const double SurfTens;
    const Point3DCL g;

    DimLessCoeffCL( const ParamMesszelleNsCL& C)
      : rho( JumpCL( 1., C.rhoF/C.rhoD ), H_sm, C.sm_eps),
        mu ( JumpCL( 1., C.muF/C.muD),    H_sm, C.sm_eps),
        SurfTens( C.sigma/C.rhoD), g( C.g)    {}
};

class EllipsoidCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;

  public:
    EllipsoidCL( const Point3DCL& Mitte, const Point3DCL& Radius)
    { Init( Mitte, Radius); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& Radius)
    { Mitte_= Mitte;    Radius_= Radius; }
    static double DistanceFct( const Point3DCL& p)
    {
        Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() { return 4./3.*M_PI*Radius_[0]*Radius_[1]*Radius_[2]; }
};

Point3DCL EllipsoidCL::Mitte_;
Point3DCL EllipsoidCL::Radius_;

// collision setting (rising butanol droplet in water)
//  RadDrop1 =  1.50e-3  1.500e-3  1.50e-3
//  PosDrop1 =  6.00e-3  3.000e-3  6.00e-3
//  RadDrop2 =  0.75e-3  0.750e-3  0.75e-3
//  PosDrop2 =  6.00e-3  5.625e-3  6.00e-3
//  MeshFile =  12e-3x30e-3x12e-3@4x10x4

class TwoEllipsoidCL
{
  private:
    static Point3DCL Mitte1_,  Mitte2_;
    static Point3DCL Radius1_, Radius2_;

  public:
    TwoEllipsoidCL( const Point3DCL& Mitte1, const Point3DCL& Radius1,
                    const Point3DCL& Mitte2, const Point3DCL& Radius2)
    { Init( Mitte1, Radius1, Mitte2, Radius2); }
    static void Init( const Point3DCL& Mitte1, const Point3DCL& Radius1,
                      const Point3DCL& Mitte2, const Point3DCL& Radius2)
    { Mitte1_= Mitte1;    Radius1_= Radius1;
      Mitte2_= Mitte2;    Radius2_= Radius2;}
    static double DistanceFct( const Point3DCL& p)
    {
        Point3DCL d1= p - Mitte1_;
        const double avgRad1= cbrt(Radius1_[0]*Radius1_[1]*Radius1_[2]);
        d1/= Radius1_;
        Point3DCL d2= p - Mitte2_;
        const double avgRad2= cbrt(Radius2_[0]*Radius2_[1]*Radius2_[2]);
        d2/= Radius2_;
        return std::min(std::abs( avgRad1)*d1.norm() - avgRad1, std::abs( avgRad2)*d2.norm() - avgRad2);
    }
    static double GetVolume() { return 4./3.*M_PI*Radius1_[0]*Radius1_[1]*Radius1_[2]
        + 4./3.*M_PI*Radius2_[0]*Radius2_[1]*Radius2_[2]; }
};

Point3DCL TwoEllipsoidCL::Mitte1_;
Point3DCL TwoEllipsoidCL::Radius1_;
Point3DCL TwoEllipsoidCL::Mitte2_;
Point3DCL TwoEllipsoidCL::Radius2_;

class InterfaceInfoCL
{
  private:
    std::ofstream* file_;    ///< write information, to this file

  public:
    Point3DCL bary, vel, min, max;
    double maxGrad, Vol, h_min, h_max;


    template<class DiscVelSolT>
    void Update (const LevelsetP2CL& ls, const DiscVelSolT& u) {
        ls.GetInfo( maxGrad, Vol, bary, vel, u, min, max);
        std::pair<double, double> h= h_interface( ls.GetMG().GetTriangEdgeBegin( ls.Phi.RowIdx->TriangLevel()), ls.GetMG().GetTriangEdgeEnd( ls.Phi.RowIdx->TriangLevel()), ls.Phi);
        h_min= h.first; h_max= h.second;
    }
    void WriteHeader() {
        if (file_)
          (*file_) << "# time maxGradPhi volume bary_drop min_drop max_drop vel_drop h_min h_max" << std::endl;
    }
    void Write (double time) {
        if (file_)
          (*file_) << time << " " << maxGrad << " " << Vol << " " << bary << " " << min << " " << max << " " << vel << " " << h_min << " " << h_max << std::endl;
    }
    /// \brief Set file for writing
    void Init(std::ofstream* file) { file_= file; }
} IFInfo;

double eps=5e-4, // halbe Sprungbreite
    lambda=1.5, // Position des Sprungs zwischen Oberkante (lambda=0) und Schwerpunkt (lambda=1)
    sigma_dirt_fac= 0.8; // gesenkte OFspannung durch Verunreinigungen im unteren Teil des Tropfens
double sigma;

double sm_step(const Point3DCL& p)
{
    double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
        y= p[1] - y_mid;
    if (y > eps) return sigma;
    if (y < -eps) return sigma_dirt_fac*sigma;
    const double z=y/eps*M_PI/2.;
    return sigma_dirt_fac*sigma + (sigma - sigma_dirt_fac*sigma) * (std::sin(z)+1)/2;
}

Point3DCL grad_sm_step (const Point3DCL& p)
{
    double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
        y= p[1] - y_mid;
    Point3DCL ret;
    if (y > eps) return ret;
    if (y < -eps) return ret;
    const double z=y/eps*M_PI/2.;
    ret[1]= (sigma - sigma_dirt_fac*sigma) * std::cos(z)/eps*M_PI/4;
    return ret;
}

double lin(const Point3DCL& p)
{
    const double y_top= IFInfo.max[1],
                 y_bot= IFInfo.bary[1],
                 y_slope= sigma*(1 - sigma_dirt_fac)/(y_top - y_bot);
    return sigma + (p[1] - y_top)*y_slope;
}

Point3DCL grad_lin (const Point3DCL&)
{
    const double y_top= IFInfo.max[1],
                 y_bot= IFInfo.bary[1],
                 y_slope= sigma*(1 - sigma_dirt_fac)/(y_top - y_bot);
    Point3DCL ret;
    ret[1]= y_slope;
    return ret;
}

double sigmaf (const Point3DCL&, double) { return sigma; }
Point3DCL gsigma (const Point3DCL&, double) { return Point3DCL(); }

double sigma_step(const Point3DCL& p, double) { return sm_step( p); }
Point3DCL gsigma_step (const Point3DCL& p, double) { return grad_sm_step( p); }

/// \brief Timedepending zero function as a helper function
SVectorCL<3> Null( const Point3DCL&, double) { return SVectorCL<3>(0.); }
double       One ( const Point3DCL&)         { return 1.; }

/// \brief Create geometry of a Mzelle of a brick
void CreateGeom (MultiGridCL* &mgp, StokesBndDataCL* &bnddata,
                 instat_vector_fun_ptr inflow,
                 const std::string& meshfile_name,
                 int GeomType, int bnd_type,
                 const std::string& deserialization_file, double& r_inlet)
{
#ifdef _PAR
    ParMultiGridCL::InstancePtr();
#endif
    if (GeomType == 0) {
        std::ifstream meshfile( meshfile_name.c_str());
        if (!meshfile)
            throw DROPSErrCL ("error while opening mesh file\n");

        ReadMeshBuilderCL *mgb= 0;       // builder of the multigrid

        // read geometry information from a file and create the multigrid
        IF_MASTER
            mgb = new ReadMeshBuilderCL( meshfile );
        IF_NOT_MASTER
            mgb = new EmptyReadMeshBuilderCL( meshfile );
        // Create the multigrid
        if (deserialization_file == "none")
            mgp= new MultiGridCL( *mgb);
        else {
            FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new MultiGridCL( filebuilder);
        }
        const BoundaryCL& bnd= mgp->GetBnd();
        const BndIdxT num_bnd= bnd.GetNumBndSeg();

        BndCondT* bc = new BndCondT[num_bnd];
        StokesVelBndDataCL::bnd_val_fun* bnd_fun = new StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= mgb->GetBC( i))==DirBC ? inflow : &ZeroVel;
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
        }
        bnddata = new StokesBndDataCL(num_bnd, bc, bnd_fun);
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
            throw DROPSErrCL("error while reading geometry information: " + mesh);
        r_inlet= dx/2;
        Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;

        BrickBuilderCL *mgb = 0;
        IF_MASTER
            mgb = new BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
        IF_NOT_MASTER
            mgb = new EmptyBrickBuilderCL(orig, px, py, pz);

        if (deserialization_file == "none")
            mgp= new MultiGridCL( *mgb);
        else {
            FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new MultiGridCL( filebuilder);
        }
        BndCondT bc[6]= { Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC };
        StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[6]=
            { &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel };
        switch (bnd_type)
        {
            case 1 : // hom. Dirichlet bnd-data
              break;
            case 2 : // inhom. Dirichlet bnd-data for in-/outflow
            {
                bc[2]= bc[3]= DirBC;
                bfun[2]= bfun[3]= inflow;
            } break;
            case 3 : // tube/canal
            {
                bc[3]= DirBC;
                bc[2]= NatBC; //Rohr
                //bc[2]= bc[4]= bc[5]= NatBC;          //Kanal
                bfun[2]= &ZeroVel;
                //bfun[2]=bfun[4]=bfun[5]= &ZeroVel;   //Kanal
                bfun[3]= inflow;
            } break;
            default: throw DROPSErrCL("Unknown boundary data type");
        }
        bnddata = new StokesBndDataCL(6, bc, bfun);
        delete mgb;
    }
}

/// \brief Display a detailed list of unknowns
template <typename StokesT, typename LevelsetT>
  void DisplayUnks(const StokesT& Stokes, const LevelsetT& levelset, __UNUSED__ const MultiGridCL& MG)
/** This functions write information about unknowns on the display. These
    informations are for the level-set-, pressure- and velocity-DOF:
    - global DOF
    - accumulated DOF
    - max and min DOF on a single processor (and the ratio)
    - max and min number of distributed DOF on a processor (and the ratio to the remaining DOF)
*/
{
#ifndef _PAR
    std::cout << Stokes.p.Data.size() << " pressure unknowns,\n";
    std::cout << Stokes.v.Data.size() << " velocity unknowns,\n";
    std::cout << levelset.Phi.Data.size() << " levelset unknowns.\n";
#else
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
        std::cout << "  + Number of DOF\n        "
                  << std::setw(10)<<"global"<<std::setw(10)<<"accum"<<std::setw(10)
                  << "max"<<std::setw(10)<<"min"<<std::setw(10)<<"ratio"<<"  |  "
                  << std::setw(10)<<"max_acc" <<std::setw(10)<<"min_acc"<<std::setw(10)<<"ratio_acc"<<std::endl;

        std::cout << "    "<<"pr  "
                  << std::setw(10)<<GPsize<<std::setw(10)<<Psize_acc<<std::setw(10)<<P_max
                  << std::setw(10)<<P_min<< std::setw(10)<<P_ratio<<"  |  "
                  << std::setw(10)<<P_accmax<<std::setw(10)<<P_accmin<<std::setw(10)<<P_accratio<<std::endl;

        std::cout << "    "<<"vel "
                  << std::setw(10)<<GVsize<<std::setw(10)<<Vsize_acc<<std::setw(10)<<V_max
                  << std::setw(10)<<V_min<< std::setw(10)<<V_ratio<<"  |  "
                  << std::setw(10)<<V_accmax<<std::setw(10)<<V_accmin<<std::setw(10)<<V_accratio<<std::endl;

        std::cout << "    "<<"scl "
                  << std::setw(10)<<GLsize<<std::setw(10)<<Lsize_acc<<std::setw(10)<<L_max
                  << std::setw(10)<<L_min<< std::setw(10)<<L_ratio<<"  |  "
                  << std::setw(10)<<L_accmax<<std::setw(10)<<L_accmin<<std::setw(10)<<L_accratio<<std::endl;

        std::cout << std::endl;
    }
#endif
}

void DisplayDetailedGeom(MultiGridCL& mg)
{
#ifndef _PAR
    mg.SizeInfo( std::cout);
#else
    const Uint level=mg.GetLastLevel();
    Uint *numTetrasAllProc=0;
    Uint *numFacesAllProc=0;
    Uint *numDistFaceAllProc=0;
    if (ProcCL::IamMaster()){
        numTetrasAllProc  = new Uint[ProcCL::Size()];
        numFacesAllProc   = new Uint[ProcCL::Size()];
        numDistFaceAllProc= new Uint[ProcCL::Size()];
    }
    // Gather information about distribution on master processor
    ProcCL::Gather(mg.GetNumTriangTetra(level),      numTetrasAllProc,   ProcCL::Master());
    ProcCL::Gather(mg.GetNumTriangFace(level),       numFacesAllProc,    ProcCL::Master());
    ProcCL::Gather(mg.GetNumDistributedFaces(level), numDistFaceAllProc, ProcCL::Master());

    // Display information
    if (ProcCL::IamMaster()){
        double ratioTetra       =  (double)*std::max_element(numTetrasAllProc,   numTetrasAllProc+ProcCL::Size())
                                  /(double)*std::min_element(numTetrasAllProc,   numTetrasAllProc+ProcCL::Size());
        Uint allTetra    =  std::accumulate(numTetrasAllProc, numTetrasAllProc+ProcCL::Size(), 0),
                    allFace     =  std::accumulate(numFacesAllProc, numFacesAllProc+ProcCL::Size(), 0),
                    allDistFace =  std::accumulate(numDistFaceAllProc, numDistFaceAllProc+ProcCL::Size(), 0);
        double      *ratioDistFace=new double[ProcCL::Size()];

        // global information
        std::cout << "Detailed information about the parallel multigrid:\n"
                  << "#(master tetras on finest level):    "<<allTetra<<'\n'
                  << "#(all Faces on finest level):        "<<allFace<<'\n'
                  << "#(distributed Faces on fines level): "<<allDistFace<<'\n';

        // local information for all processors
        for (int i=0; i<ProcCL::Size(); ++i)
            ratioDistFace[i]= ((double)numDistFaceAllProc[i]/(double)numFacesAllProc[i]*100.);

        double maxRatio= *std::max_element(ratioDistFace, ratioDistFace+ProcCL::Size());
        std::cout << "Ratio between max/min Tetra: "<<ratioTetra
                  <<" max Ratio DistFace/AllFace: "<<maxRatio<<std::endl;

        std::cout << std::setw(6)  <<  "Proc"
                  << std::setw(8)  << "#Tetra"
                  << std::setw(8)  << "#Faces"
                  << std::setw(12) << "#DistFaces"
                  << std::setw(12) << "%DistFaces"
                  << '\n';
        for (int i=0; i<ProcCL::Size(); ++i)
            std::cout << std::setw(6)  << i
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

/// \brief Write finite element function, stored in \a v, in a file, named \a filename
void WriteFEToFile( const VecDescCL& v, MultiGridCL& mg, std::string filename, bool binary=false, const VecDescCL* lsetp=0)
{
    if (!v.RowIdx->IsExtended()) {
#ifdef _PAR
        ProcCL::AppendProcNum( filename);
#endif
        std::ofstream file( filename.c_str());
        if (!file) throw DROPSErrCL("WriteFEToFile: Cannot open file "+filename+" for writing");
        v.Write( file, binary);
    }
    else { // extended FE
        IdxDescCL p1( P1_FE);
        p1.CreateNumbering( v.RowIdx->TriangLevel(), mg, *v.RowIdx);
        VecDescCL vpos(&p1), vneg(&p1);
        P1XtoP1 ( *v.RowIdx, v.Data, p1, vpos.Data, vneg.Data, *lsetp, mg);
        WriteFEToFile(vneg, mg, filename + "Neg");
        WriteFEToFile(vpos, mg, filename + "Pos");
        p1.DeleteNumbering(mg);
    }
}

/// Read a serialized finite element function from a file
/// \pre CreateNumbering of v.RowIdx must have been called before
void ReadFEFromFile( VecDescCL& v, MultiGridCL& mg, std::string filename, bool binary=false, const VecDescCL* lsetp=0)
{
    if (!v.RowIdx->IsExtended()) {

        std::cout << "Read FE "<<filename<<std::endl;
#ifdef _PAR
        ProcCL::AppendProcNum( filename);
#endif
        std::ifstream file( filename.c_str());
        if (!file) throw DROPSErrCL("ReadFEFromFile: Cannot open file "+filename);
        v.Read( file, binary);
    }
    else { // extended FE
        IdxDescCL p1( P1_FE);
        p1.CreateNumbering( v.RowIdx->TriangLevel(), mg, *v.RowIdx);
        VecDescCL vpos(&p1), vneg(&p1);
        ReadFEFromFile(vneg, mg, filename + "Neg");
        ReadFEFromFile(vpos, mg, filename + "Pos");
        P1toP1X ( *v.RowIdx, v.Data, p1, vpos.Data, vneg.Data, *lsetp, mg);
        p1.DeleteNumbering(mg);
    }
}

/// \brief Class for serializing a two-phase flow problem, i.e., storing
///    the multigrid and the numerical data
/// \todo  Storing of transport data!
template <typename StokesT>
class TwoPhaseStoreCL
{
  private:
    MultiGridCL&         mg_;
    const StokesT&       Stokes_;
    const LevelsetP2CL&  lset_;
    const TransportP1CL* transp_;
    std::string          path_;
    Uint                 numRecoverySteps_;
    Uint                 recoveryStep_;

    /// \brief Write time info
    void WriteTime( std::string filename)
    {
        std::ofstream file( filename.c_str());
        if (!file) throw DROPSErrCL("TwoPhaseStoreCL::WriteTime: Cannot open file "+filename+" for writing");
        file << Stokes_.t << "\n";
        file.close();
    }

  public:
      /// \brief Construct a class for storing a two-phase flow problem in files
      /** This class generates multiple files, all with prefix path, for storing
       *  the geometric as well as the numerical data.
       *  \param recoverySteps number of backup steps before overwriting files
       *  */
    TwoPhaseStoreCL(MultiGridCL& mg, const StokesT& Stokes, const LevelsetP2CL& lset, const TransportP1CL* transp,
                    const std::string& path, Uint recoverySteps=2)
      : mg_(mg), Stokes_(Stokes), lset_(lset), transp_(transp), path_(path), numRecoverySteps_(recoverySteps),
        recoveryStep_(0) {}

    /// \brief Write all information in a file
    void Write()
    {
        // Create filename
        std::stringstream filename;
        filename << path_ << ((recoveryStep_++)%numRecoverySteps_);
        // first master writes time info
        IF_MASTER
            WriteTime( filename.str() + "time");

        // write multigrid
        MGSerializationCL ser( mg_, filename.str());
        ser.WriteMG();

        // write numerical data
        WriteFEToFile(Stokes_.v, mg_, filename.str() + "velocity");
        WriteFEToFile(lset_.Phi, mg_, filename.str() + "levelset");
        WriteFEToFile(Stokes_.p, mg_, filename.str() + "pressure", false, &lset_.Phi); // pass also level set, as p may be extended
        if (transp_) WriteFEToFile(transp_->ct, mg_, filename.str() + "concentrationTransf");
    }
};

}   // end of namespace DROPS

#endif
