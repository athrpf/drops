//**************************************************************************
// File:    two_phase_hdr.h                                                *
// Content: header for all parallel main programs simulating a             *
//          two-phase flow                                                 *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   03. September 2008                                             *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file two_phase_hdr.h
/// \brief Header for all parallel main programs simulating a two-phase flow

#ifndef DROPS_TWO_PHASE_HDR_H
#define DROPS_TWO_PHASE_HDR_H

#include "geom/multigrid.h"
#include "misc/utils.h"
#include "num/discretize.h"
#include "stokes/stokes.h"
#include "partests/params.h"
#include "levelset/levelset.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

/// \brief Display a detailed list of unknowns
template <typename StokesT, typename LevelsetT>
  void DisplayUnks(const StokesT& Stokes, const LevelsetT& levelset, const DROPS::MultiGridCL& MG)
/** This functions write information about unknowns on the display. These
    informations are for the level-set-, pressure- and velocity-DOF:
    - global DOF
    - accumulated DOF
    - max and min DOF on a single processor (and the ratio)
    - max and min number of distributed DOF on a processor (and the ratio to the remaining DOF)
*/
{
    using namespace DROPS;
    const MLIdxDescCL* vidx = &Stokes.vel_idx,
                     * pidx = &Stokes.pr_idx;
    const IdxDescCL*   lidx = &levelset.idx;
    const ExchangeCL& ExV = Stokes.GetEx(Stokes.velocity),
                    & ExP = Stokes.GetEx(Stokes.pressure),
                    & ExL = levelset.GetEx();

    // local number on unknowns
    Ulint Psize      = pidx->NumUnknowns();
    Ulint Vsize      = vidx->NumUnknowns();
    Ulint Lsize      = lidx->NumUnknowns();

    // global number of unknowns
    Ulint GPsize     = pidx->GetGlobalNumUnknowns(MG);
    Ulint GVsize     = vidx->GetGlobalNumUnknowns(MG);
    Ulint GLsize     = lidx->GetGlobalNumUnknowns(MG);

    // accumulated size of unknwons
    Ulint Psize_acc = GlobalSum(Psize);
    Ulint Vsize_acc = GlobalSum(Vsize);
    Ulint Lsize_acc = GlobalSum(Lsize);

    // maximal and minimal number of unknowns
    Ulint P_min= GlobalMin(Psize); Ulint P_max= GlobalMax(Psize);
    Ulint V_min= GlobalMin(Vsize); Ulint V_max= GlobalMax(Vsize);
    Ulint L_min= GlobalMin(Lsize); Ulint L_max= GlobalMax(Lsize);

    // ratios between maximal number of unknowns/proc and minimal number
    double P_ratio   = (double)P_max/(double)P_min;
    double V_ratio   = (double)V_max/(double)V_min;
    double L_ratio   = (double)L_max/(double)L_min;

    // number on boundaries
    Ulint P_accmax=GlobalMax(ExP.AccDistIndex.size()), P_accmin=GlobalMin(ExP.AccDistIndex.size());
    Ulint V_accmax=GlobalMax(ExV.AccDistIndex.size()), V_accmin=GlobalMin(ExV.AccDistIndex.size());
    Ulint L_accmax=GlobalMax(ExL.AccDistIndex.size()), L_accmin=GlobalMin(ExL.AccDistIndex.size());

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
}

namespace DROPS{
/// \brief Parameter class describing a zero flow
class ZeroFlowCL
{
// \Omega_1 = droplet,    \Omega_2 = bulk phase
  public:
    static Point3DCL f(const Point3DCL&, double)
      { Point3DCL ret(0.0); return ret; }
    const SmoothedJumpCL rho, mu;
    const double SurfTens;
    const Point3DCL g;

    ZeroFlowCL( const ParamMesszelleCL& C)
      : rho( JumpCL( C.rhoD, C.rhoF ), H_sm, C.sm_eps),
        mu(  JumpCL( C.muD,  C.muF),   H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};
}

namespace DROPS{
/// \brief Information about a droplet
class InterfaceInfoCL
{
  private:
    std::ofstream *infofile_;   ///< Pointer to a file, where to write out information

  public:
    InterfaceInfoCL() : infofile_(0) {}

    void Init(std::ofstream *info){
        infofile_=info;
    }

    Point3DCL bary;             ///< barycenter of a droplet
    Point3DCL min;              ///< bottom of a droplet
    Point3DCL max;              ///< top of a droplet
    Point3DCL vel;              ///< (accumulated) velocity of the droplet
    double maxGrad;             ///< maximal 2-norm of the gradient of the level set function
    double Vol;                 ///< volume of a droplet

    /// \brief Update
    template <typename LsetCL, typename VelDesc>
    void Update(const LsetCL& ls, const VelDesc& u){
        ls.GetInfo( maxGrad, Vol, bary, vel, u, min, max);
    }

    /// \brief Write information in a file
    void Write(double time){
        if (!infofile_) return;
        (*infofile_) << time << '\t' << maxGrad << '\t' << Vol << '\t' << bary << '\t' << vel
                     << '\t' << min << '\t' << max << std::endl;
    }
} IFInfo;
}

namespace DROPS{
/// \brief Class for describing surface tension forces
/** This class needs access to an InterfaceInfoCL named as IFInfo.*/
class SurfaceTensionCL
{
  public:
    static double eps;              ///< depth of the jump
    static double lambda;           ///< position of the jump (top: lambda=0, barycenter lambda=1)
    static double sigma;            ///< surface tension at the top of the drop
    static double sigma_dirt_fac;   ///< factor of surface tension at the bottom of the drop

    /// \brief Constant surface tension
    static double sigmaf (const Point3DCL&, double) {
        return sigma;
    }

    /// \brief Gradient of constant surface tension
    static Point3DCL gsigma (const Point3DCL&, double) {
        return Point3DCL();
    }

    /// \brief Surface tension modeled by a step
    static double sigma_step(const Point3DCL& p, double)
    {
        double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
               y    = p[1] - y_mid;
        if (y > eps) return sigma;
        if (y < -eps) return sigma_dirt_fac*sigma;
        const double z=y/eps*M_PI/2.;
        return sigma_dirt_fac*sigma + (sigma - sigma_dirt_fac*sigma) * (std::sin(z)+1)/2;
    }

    /// \brief Gradient of surface tension modeled by a step
    static Point3DCL gsigma_step (const Point3DCL& p, double)
    {
        double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
               y    = p[1] - y_mid;
        Point3DCL ret;
        if (y > eps) return ret;
        if (y < -eps) return ret;
        const double z=y/eps*M_PI/2.;
        ret[1]= (sigma - sigma_dirt_fac*sigma) * std::cos(z)/eps*M_PI/4;
        return ret;
    }
};
}

/// \brief Class describing an ellipsoid
class EllipsoidCL
{
  private:
    static DROPS::Point3DCL center_;        // center of the ellipsoid
    static DROPS::Point3DCL radius_;        // radi of the ellipsoid

  public:
    EllipsoidCL( const DROPS::Point3DCL& center, const DROPS::Point3DCL& radi){
        Init( center, radi);
    }

    /// \brief Initialization of an ellipsoid
    static void Init( const DROPS::Point3DCL& center, const DROPS::Point3DCL& radi){
        center_= center;
        radius_= radi;
    }

    /// \brief Get (signed) distance to the surface of an ellipsoid
    static double DistanceFct( const DROPS::Point3DCL& p)
    {
        DROPS::Point3DCL d= p - center_;
        const double avgRad= cbrt(radius_[0]*radius_[1]*radius_[2]);
        d/= radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }

    /// \brief Get Volume of an ellipsoid
    static double GetVolume() {
        return 4./3.*M_PI*radius_[0]*radius_[1]*radius_[2];
    }
};

DROPS::Point3DCL EllipsoidCL::center_;
DROPS::Point3DCL EllipsoidCL::radius_;

namespace DROPS{
/// \brief Timedepending zero function as a helper function
SVectorCL<3> Null( const Point3DCL&, double) { return DROPS::SVectorCL<3>(0.); }
double       One ( const Point3DCL&)         { return 1.; }
}

/*******************************************************************************
*    C R E A T E   G E O M E T R I E S                                         *
*******************************************************************************/


enum GeometryType{
    newMZelle,             // Messzelle-geometry
    fileMZelle,            // Messzelle-geometry from file
    newBrick,              // brick-geometry
    fileBrick              // brick-geometry from file
};

namespace DROPS{

/// \brief Create geometry of a Mzelle of a brick
void CreateGeom( MultiGridCL* &mg, ParMultiGridCL* &pmg,
                 LoadBalHandlerCL* &lb, StokesBndDataCL* &bnddata,
                 instat_vector_fun_ptr inflow,
                 const std::string& meshfile_name,
                 const Uint refineStrategy,
                 const GeometryType& geomType,
                 double &newInlet
                )
{
    if (ProcCL::IamMaster())
        std::cerr << "Creating initial geometry for "<<ProcCL::Size()<<" processors" << std::endl;

    // parallel multigrid, that is able to handle 3 dof
    pmg = new ParMultiGridCL;

    StokesVelBndDataCL::bnd_val_fun *bnd_fun = 0;       // boundary information for velocity
    BndCondT* bc= 0;
    BndIdxT num_bnd= 0;

    // Create a Mzelle-shaped geometry
    if (geomType==newMZelle || geomType==fileMZelle){
        if (geomType==fileMZelle)
            throw DROPSErrCL("Sorry, read serialized input for mzelle not yet implemented");

        // open file of geometry information
        std::ifstream meshfile( meshfile_name.c_str());
        if (!meshfile)
            throw DROPSErrCL("Error opening mesh file");
        ReadMeshBuilderCL *mgb= 0;       // builder of the multigrid

        // read geometry information from a file and create the multigrid
        if (DROPS::ProcCL::IamMaster())
            mgb = new ReadMeshBuilderCL( meshfile );
        else
            mgb = new EmptyReadMeshBuilderCL( meshfile );
        // Create the multigrid
        mg = new MultiGridCL( *mgb );

        // Get boundary conditions
        num_bnd= mg->GetBnd().GetNumBndSeg();
        bc = new BndCondT[num_bnd];
        for (BndIdxT i=0; i<num_bnd; ++i)
            bc[i]= mgb->GetBC( i);

        // Create boundary information for velocity (pressure and level-set do have natural Neumann boundary conditions)
        //const BoundaryCL& bnd    = mg->GetBnd();
        bnd_fun = new StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (DROPS::BndIdxT i=0; i<num_bnd; ++i)
            bnd_fun[i]= (bc[i]= mgb->GetBC( i))==DirBC ? inflow : &Null;

        // Free memory
        if (mgb) delete mgb;
    }

    // Create a brick-shaped geometry
    if (geomType==newBrick || geomType==fileBrick){
        int nx, ny, nz;
        double dx, dy, dz;
        std::string mesh( meshfile_name), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
        if (!brick_info)
            throw DROPSErrCL("error while reading geometry information");
        newInlet= dx/2;
        DROPS::Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;
        MGBuilderCL    *mgb=0;
        BrickBuilderCL *builder = 0;
        if (DROPS::ProcCL::IamMaster()){
            if (geomType==newBrick){
                mgb = new DROPS::BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
            }
            else if (geomType==fileBrick){   // read geometry stored in a file
                builder = new DROPS::BrickBuilderCL(orig, px, py, pz, nx, ny, nz);
                mgb = new DROPS::FileBuilderCL(meshfile_name, builder);
            }
        }
        else{
            mgb = new DROPS::EmptyBrickBuilderCL(orig, px, py, pz);
        }
        // Create the multigrid
        mg = new MultiGridCL( *mgb );

        // Set boundary conditions
        num_bnd = 6;
        bc = new BndCondT[num_bnd];
        bc[0]=WallBC; bc[1]=WallBC;       // wall at x=0, x=x_1
        bc[2]=Nat0BC; bc[3]=DirBC;        // outflow at the bottom, inflow from the top
        bc[4]=WallBC; bc[5]=WallBC;       // wall at z=0, z=z_1

        // Create boundary information for velocity (pressure and level-set do have natural Neumann boundary conditions)
        bnd_fun = new StokesVelBndDataCL::bnd_val_fun[num_bnd];
        bnd_fun[0]= &Null; bnd_fun[1]= &Null;
        bnd_fun[2]= &Null; bnd_fun[3]= inflow;
        bnd_fun[4]= &Null; bnd_fun[5]= &Null;

        // Free memory
        if (mgb)     delete mgb;
        if (builder) delete builder;
    }

    // Create the multigrid and tell parallel multigrid about the geometry
    pmg->AttachTo(*mg);

    // Create a load balancer and do initial distribution of the geometry
    lb = new LoadBalHandlerCL(*mg);
    lb->DoInitDistribution(ProcCL::Master());
    switch (refineStrategy){
        case 0 : lb->SetStrategy(NoMig);     break;
        case 1 : lb->SetStrategy(Adaptive);  break;
        case 2 : lb->SetStrategy(Recursive); break;
    }

    // Create boundary date
    bnddata = new DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun);

    if (bc)      delete[] bc;
    if (bnd_fun) delete[] bnd_fun;
}
}

namespace DROPS{

/// \brief Create all numberings for stokes DOF and level-set DOF and assign them
template<typename StokesT, typename LevelsetT>
void CreateIdxAndAssignIdx(StokesT& Stokes, LevelsetT& lset, const MultiGridCL& mg)
{
    // Create numbering
    Stokes.CreateNumberingVel( mg.GetLastLevel(), &Stokes.vel_idx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), &Stokes.pr_idx);
    lset.CreateNumbering     ( mg.GetLastLevel(), &lset.idx);

    // Tell matrices and vectors about the numbering
    Stokes.v.SetIdx(   &Stokes.vel_idx);
    Stokes.p.SetIdx(   &Stokes.pr_idx);
    Stokes.b.SetIdx(   &Stokes.vel_idx);
    Stokes.c.SetIdx(   &Stokes.pr_idx);
    Stokes.A.SetIdx(   &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.B.SetIdx(   &Stokes.pr_idx,  &Stokes.vel_idx);
    Stokes.M.SetIdx(   &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.N.SetIdx(   &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.prM.SetIdx( &Stokes.pr_idx,  &Stokes.pr_idx);
    Stokes.prA.SetIdx( &Stokes.pr_idx,  &Stokes.pr_idx);
    lset.Phi.SetIdx(   &lset.idx);
}
}


/*******************************************************************************
*    G L O B A L   A N D   S T A T I C   V A R I A B L E S                     *
*******************************************************************************/
namespace DROPS{

const char line[] ="------------------------------------------------------------";

// Init of static members
double SurfaceTensionCL::eps           = 5e-4;
double SurfaceTensionCL::lambda        = 1.5;
double SurfaceTensionCL::sigma         = 0.0;
double SurfaceTensionCL::sigma_dirt_fac= 0.8;
}

#endif
