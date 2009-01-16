//**************************************************************************
// File:    TestBrickflowPar.cpp                                           *
// Content: parallel solver for a simple 2-phase problem                   *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestBrickflowPar.cpp
/// \brief Testing parallel solver for a simple 2-phase problem

// include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include "parallel/exchange.h"
#include "parallel/logger.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "num/parsolver.h"
#include "num/parprecond.h"
#include "num/stokessolver.h"
#include "num/parstokessolver.h"
#include "stokes/integrTime.h"
#include "levelset/adaptriang.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "parallel/parmgserialization.h"

 // include problem class
#include "navstokes/instatnavstokes2phase.h"
#include "levelset/coupling.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#endif


enum TimePart{
    T_create_pmg,
    T_disc_init,
    T_create_ex,
    T_solve_NS_init
};

DROPS::TimeStoreCL Times(4);

void SetDescriber()
{
    Times.SetDescriber(T_create_pmg,    "Initialization of parallel multigrid");
    Times.SetDescriber(T_disc_init,     "Discretizing Stokes/Curv for initial velocities");
    Times.SetDescriber(T_create_ex,     "Create ExchangeCL");
    Times.SetDescriber(T_solve_NS_init, "Solving NavStokes for initial velocities");
    Times.SetCounterDescriber("Moved MultiNodes");
}

DROPS::ParamParBrickFlowCL C;

const char line[] ="------------------------------------------------------------";

/****************************************************************************
* P R I N T  M U L T I G R I D  I N  A S C I I  F I L E                     *
****************************************************************************/


bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream check(dat);
    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);
    check.close();
    return DROPS::Check(pmg_sane && mg_sane);
}

/// \brief Display number of unknowns
void DisplayNumUnknownsSparse(const DROPS::MultiGridCL& MG, const DROPS::VecDescCL& x, const std::string& name)
/// accumulated unknowns: accumulation of the unknowns of all processors. Some unknowns
///   are counted multiple due to the overlapping of edges and vertices
/// global unknowns: all unknowns are just once counted
{
    const DROPS::Ulint acc_num_unk = DROPS::GlobalSum(x.Data.size()),
                       glo_num_unk = x.RowIdx->GetGlobalNumUnknowns(MG);
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    if (DROPS::ProcCL::IamMaster())
        std::cerr << "  + Number of "<<name<<" DOF with index "<<idx<<" (accumulated/global):  "
                  <<acc_num_unk<< "/" <<glo_num_unk<< std::endl;
}

/// \brief Display a detailed list of unknowns
template <typename StokesT, typename LevelsetT>
  void DisplayUnks(const StokesT& Stokes, const LevelsetT& levelset, const DROPS::MultiGridCL& MG)
/** This functions write information about unknowns on the display. These
    informations are for the level-set-, pressure- and velocity-DOF:
    <ul>
     <li> global DOF
     <li> accumulated DOF
     <li> max and min DOF on a single processor (and the ratio)
     <li> max and min number of distributed DOF on a processor (and the ratio to the remaining DOF)
    </ul>
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

void DisplayDetailedGeom(DROPS::MultiGridCL& mg)
{
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
    DROPS::Gather(mg.GetNumTriangTetra(level),      numTetrasAllProc,   DROPS::ProcCL::Master());
    DROPS::Gather(mg.GetNumTriangFace(level),       numFacesAllProc,    DROPS::ProcCL::Master());
    DROPS::Gather(mg.GetNumDistributedFaces(level), numDistFaceAllProc, DROPS::ProcCL::Master());

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
}

class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
      { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( C.rhoD, C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muD,  C.muF),   DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

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

    /// \brief Surface tension modeled by a step
    static double sm_step(const Point3DCL& p, double){
        double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
               y    = p[1] - y_mid;
        if (y > eps) return sigma;
        if (y < -eps) return sigma_dirt_fac*sigma;
        const double z=y/eps*M_PI/2.;
        return sigma_dirt_fac*sigma + (sigma - sigma_dirt_fac*sigma) * (std::sin(z)+1)/2;
    }

    /// \brief Gradient of surface tension modeled by a step
    static Point3DCL grad_sm_step (const Point3DCL& p, double)
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

// Init of static members
double SurfaceTensionCL::eps           = 5e-4;
double SurfaceTensionCL::lambda        = 1.5;
double SurfaceTensionCL::sigma         = 0.0;
double SurfaceTensionCL::sigma_dirt_fac= 0.8;
}

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

DROPS::SVectorCL<3> Inflow( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = (p[0]/C.r_inlet)*(p[0]/C.r_inlet)-1,
                 z = (p[2]/C.r_inlet)*(p[2]/C.r_inlet)-1;

//     Inflow with waves
//     const double freq= 50, ampl= 0.1;
//     ret[1]= x * z * C.Anstroem * (1-ampl*std::cos(2*M_PI*freq*t));  // Rohr

    // constant inflow
    ret[1]= x * z * C.Anstroem;

    return ret;
}

// distance function of a droplet droplet
double DistanceFct1( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL d= C.Mitte-p;
    const double avgRad = cbrt(C.Radius[0]*C.Radius[1]* C.Radius[2]);
    d/= C.Radius/avgRad;
    return d.norm()-avgRad;
}

// middle
double DistanceFct( const DROPS::Point3DCL& p)
{
    return p[1]-0.015;
}

namespace DROPS{

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

template<class Coeff>
  void InitProblemWithDrop(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    ParTimerCL time;
    double duration;
    // droplet-information

    ExchangeCL& ExV = Stokes.GetEx(Stokes.velocity);
    ExchangeCL& ExP = Stokes.GetEx(Stokes.pressure);
//     ExchangeCL& ExL = lset.GetEx();
    MultiGridCL& mg=pmg.GetMG();

    switch (C.IniCond)
    {
      // zero flow
      case 0:
      {
        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        DisplayUnks(Stokes, lset, mg);

        // Initial velocity is zero
        Stokes.InitVel( &Stokes.v, Null);
        // Initial levelset function is distance to the drop
        lset.Init( DistanceFct1);

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);
      } break;
      // stationary flow with/without drop
      case 1: case 2:
      {
        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        // Initial velocity is zero
        Stokes.InitVel( &Stokes.v, Null);
        // Initial levelset function is distance to the drop
        lset.Init( DistanceFct1);

        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        // Setting up solvers
//         typedef ParDummyPcCL<ExchangeCL> SPcT;
//         SPcT ispc(ExP);
        typedef ISBBTPreCL<ExchangeCL, ExchangeCL> SPcT;
        SPcT ispc( Stokes.B.Data.GetFinestPtr(), Stokes.prM.Data.GetFinestPtr(), Stokes.M.Data.GetFinestPtr(),
                   ExP, ExV, /*kA*/1.0/C.dt, /*kM_*/C.theta, 1e-4, 1e-4);
        typedef ParJac0CL<ExchangeCL>  APcPcT;
        APcPcT Apcpc(ExV);
        typedef ParPCGSolverCL<APcPcT,ExchangeCL> ASolverT;
        ASolverT Asolver( 500, 0.02, ExV, Apcpc, /*relative=*/ true, /*accur*/ true);
        typedef SolverAsPreCL<ASolverT> APcT;
        APcT Apc( Asolver/*, &std::cerr*/);
        typedef ParInexactUzawaCL<APcT, SPcT, APC_SYM, ExchangeCL, ExchangeCL> OseenSolverT;
        OseenSolverT schurSolver( Apc, ispc, ExV, ExP, C.outer_iter, C.outer_tol, 0.1);

        VelVecDescCL curv(&Stokes.vel_idx),
                     cplN(&Stokes.vel_idx);

        time.Reset();
        Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &curv, lset, Stokes.t);
        Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.t);
        Stokes.SetupPrMass( &Stokes.prM, lset);
        curv.Clear();
        curv.Data=0.;
        lset.AccumulateBndIntegral( curv);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
        {
            std::cerr << "- Discretizing Stokes/Curv for initialization "<<duration<<" sec.\n";
            // Log measured duration.
            DROPS_LOGGER_SETVALUE("StokesCurv",duration);
        }

        //Solve initial problem
        double theta= C.theta;
        time.Reset();
        int step=0;
        int iters=0;
        do
        {
            Stokes.SetupNonlinear( &Stokes.N, &Stokes.v, &cplN, lset, Stokes.t);
            cplN.Data-= (1-theta) * (Stokes.N.Data * Stokes.v.Data);
            schurSolver.Solve( Stokes.A.Data, Stokes.B.Data,
                               Stokes.v.Data, Stokes.p.Data, Stokes.b.Data, Stokes.c.Data);
            if (ProcCL::IamMaster())
                std::cerr << "- Solving lin. Stokes ("<<step<<"): iter "<<schurSolver.GetIter()
                            <<", resid "<<schurSolver.GetResid()<<std::endl;
            ++step; iters+= schurSolver.GetIter();
        } while (schurSolver.GetIter() > 0);
        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster())
        {
            std::cerr << "- Solving Stokes for initialization took "<<duration<<" sec, "
                        << "steps "<<(step-1)<<", iter "<<iters<<", resid "<<schurSolver.GetResid()<<'\n';
            DROPS_LOGGER_SETVALUE("SolStokesTime",duration);
            DROPS_LOGGER_SETVALUE("SolStokesIter",iters);
        }
      }break;

      case 3:
      {
        // Create and number (on master) indices
        MLIdxDescCL read_vel_idx( vecP2_FE), read_pr_idx( P1_FE);
        IdxDescCL read_lset_idx( P2_FE);

        if (ProcCL::IamMaster()){
            if (C.checkMG){
                std::ofstream sanity_file("ser_sanity.txt");
                if (mg.IsSane(sanity_file))
                    std::cerr << "Multigrid is sane on master\n";
                else
                    std::cerr << "Multigrid is not sane on master\n ==> Check right dimension of the computational domain!";
            }
            Stokes.CreateNumberingVel( mg.GetLastLevel(), &read_vel_idx,  false);
            Stokes.CreateNumberingPr ( mg.GetLastLevel(), &read_pr_idx,   false);
            lset.CreateNumbering     ( mg.GetLastLevel(), &read_lset_idx, false);
        }

        // Allocate memory for unknowns
        VecDescCL read_vel_vec(&read_vel_idx), read_pr_vec(&read_pr_idx), read_lset_vec(&read_lset_idx);

        // Read DOF
        if (ProcCL::IamMaster()){
            ReadDOF(mg, &read_vel_vec,  std::string(C.ser_dir+"velocity"));
            ReadDOF(mg, &read_pr_vec,   std::string(C.ser_dir+"pressure"));
            ReadDOF(mg, &read_lset_vec, std::string(C.ser_dir+"level-set"));
        }

        // Distribute MG
        lb.DoMigration();
        if (C.checkMG && !Check( CheckParMultiGrid(pmg)) )
            throw DROPSErrCL("MultiGrid is incorrect!");

        // Create indices
        CreateIdxAndAssignIdx(Stokes, lset, mg);
        /// \todo Hier muessen noch die Level der Indices gesetzt werden. Warum?
//         read_vel_idx.TriangLevel = Stokes.vel_idx.TriangLevel;
//         read_pr_idx.TriangLevel  = Stokes.pr_idx.TriangLevel;
//         read_lset_idx.TriangLevel= lset.idx.TriangLevel;

        pmg.HandleNewIdx(&read_vel_idx,  &Stokes.v);
        pmg.HandleNewIdx(&read_pr_idx,   &Stokes.p);
        pmg.HandleNewIdx(&read_lset_idx, &lset.Phi);
        pmg.DelAllUnkRecv();
        pmg.DeleteRecvBuffer();

        Stokes.DeleteNumbering(&read_vel_idx);
        Stokes.DeleteNumbering( &read_pr_idx);
        lset.DeleteNumbering(     &read_lset_idx);

        if (ProcCL::IamMaster())
                std::cerr << "- Initial Conditions successfull read\n";
      } break;
    }

}

template<typename Coeff>
  void SolveCoupledNS(InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, LevelsetP2CL& lset, AdapTriangCL& adap)
{
    ParTimerCL time;
    double duration;
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;

    // writer for vtk-format
    typedef TwoPhaseVTKCL<StokesProblemT, LevelsetP2CL> VTKWriterT;
    VTKWriterT vtkwriter(adap.GetMG(), Stokes, lset, (C.vtk ? C.num_steps/C.vtk+1 : 0),
                         std::string(C.vtkDir + "/" + C.vtkName), C.vtkBinary);
    // writer for ensight format
    typedef Ensight2PhaseOutCL<StokesProblemT, LevelsetP2CL> EnsightWriterT;
    EnsightWriterT ensightwriter( adap.GetMG(), lset.Phi.RowIdx, Stokes, lset, C.ensDir, C.ensCase, C.geomName, /*adaptive=*/true,
                                  (C.ensight? C.num_steps/C.ensight+1 : 0), C.binary, C.masterOut);

    // Serialization
    typedef TwoPhaseSerializationCL<StokesProblemT, LevelsetP2CL> SerializationT;
    SerializationT serializer(adap.GetMG(), C.ser_dir, Stokes, lset, C.overwrite, C.num_steps/(C.overwrite+1));

    if (C.vtk)
        vtkwriter.write();
    if (C.ensight)
        ensightwriter.write();

    const double Vol= 4./3.*M_PI*C.Radius[0]*C.Radius[1]*C.Radius[2];
    double relVol = lset.GetVolume()/Vol;

    ExchangeCL& ExV = Stokes.GetEx(Stokes.velocity);
    ExchangeCL& ExP = Stokes.GetEx(Stokes.pressure);
    ExchangeCL& ExL = lset.GetEx();

    // linear solvers
    typedef ISBBTPreCL<ExchangeCL, ExchangeCL> SPcT;
    SPcT ispc( Stokes.B.Data.GetFinestPtr(), Stokes.prM.Data.GetFinestPtr(), Stokes.M.Data.GetFinestPtr(),
               ExP, ExV, /*kA*/1.0/C.dt, /*kM_*/C.theta, 1e-4, 1e-4);
//     typedef ParDummyPcCL<ExchangeCL> SPcT;
//     SPcT ispc(ExP);

//     typedef ParDummyPcCL<ExchangeCL> SPcPcT;
//     SPcPcT SPcPcA(ExP), SPcPcM(ExP);
//     typedef ParPCGSolverCL<SPcPcT,ExchangeCL> SsolverT;
//     SsolverT PrASolver(50, 0.02, ExP, SPcPcA, true);
//     SsolverT PrMSolver(50, 0.02, ExP, SPcPcM, true);
//     typedef ISNonlinearPreCL<SsolverT, SsolverT>  SPcT;
//     SPcT ispc(PrASolver, PrMSolver, Stokes.prA.Data, Stokes.prM.Data);

    typedef ParJac0CL<ExchangeCL>  APcPcT;
    APcPcT Apcpc(ExV);
    typedef ParPreGMResSolverCL<APcPcT, ExchangeCL>    ASolverT;        // GMRes-based APcT
    ASolverT Asolver( /*restart*/    100, /*iter*/ 500,  /*tol*/   2e-5, ExV, Apcpc,
                       /*relative=*/ true, /*acc*/ true, /*modGS*/ false, RightPreconditioning, /*mod4par*/true);

//     ASolverT Asolver( C.navstokes_pc_restart, C.navstokes_pc_iter, C.navstokes_pc_reltol, ExV, Apcpc,
//                         /*relative=*/ true, /*acc*/ true, /*modGS*/false, RightPreconditioning, /*mod4par*/false);

    typedef SolverAsPreCL<ASolverT> APcT;
    APcT Apc( Asolver/*, &std::cerr*/);

    // stokes solver
    typedef ParInexactUzawaCL<APcT, SPcT, APC_OTHER, ExchangeCL, ExchangeCL> OseenSolverT;
    OseenSolverT oseensolver( Apc, ispc, ExV, ExP, C.outer_iter, C.outer_tol, 0.1, 500, &std::cerr);

    // Navstokes solver
    typedef AdaptFixedPtDefectCorrCL<StokesProblemT> NSSolverT;
    NSSolverT nssolver( Stokes, oseensolver, C.ns_iter, C.ns_tol, C.ns_red);

    // coupling levelset NavStokes
    time.Reset();
    // Coupling Navier-Stokes with Levelset
//     typedef CouplLsNsFracStep2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
    typedef LinThetaScheme2PhaseCL <StokesProblemT, NSSolverT> CouplingT;
    CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear, true);
//     typedef RecThetaScheme2PhaseCL<StokesProblemT, NSSolverT> CouplingT;
//     CouplingT cpl( Stokes, lset, nssolver, C.theta, C.nonlinear, false);


    time.Stop();
    duration=time.GetMaxTime();
    if (ProcCL::IamMaster())
        std::cerr << "- Updating discretization took "<<duration<<" sec.\n";

    // Set time step and create matrices
    cpl.SetTimeStep( C.dt);

    for (int step= 1; step<=C.num_steps; ++step)
    {
        ParTimerCL step_time;
        step_time.Reset();
        if (ProcCL::IamMaster())
            std::cerr << "=================================================================================== Schritt " << step << ":\n"
                      << " Idx for vel  "<<Stokes.v.RowIdx->GetIdx()
                      << "\n Idx for pr   "<<Stokes.p.RowIdx->GetIdx()
                      << "\n Idx for lset "<<lset.Phi.RowIdx->GetIdx()<<std::endl;
        time.Reset();

        if (C.ref_freq && step%C.ref_freq==0)
        {
            if (ProcCL::IamMaster())
                std::cerr << "==> Adaptive Refinement of MultiGrid"<<std::endl;

            adap.UpdateTriang( lset);

            if (C.checkMG && !Check( CheckParMultiGrid(adap.GetPMG())) )
                throw DROPSErrCL("MultiGrid is incorrect!");

            if (adap.WasModified() )
            {
                cpl.Update();
                // don't forget to update the pr mass/stiff matrix for the schur compl. preconditioner!!
				// / \todo (of) Ueberfluessig!
                Stokes.prM.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);
                Stokes.SetupPrMass( &Stokes.prM, lset);
                Stokes.prA.SetIdx( &Stokes.pr_idx, &Stokes.pr_idx);
                Stokes.SetupPrStiff( &Stokes.prA, lset);
            }
        }

        if (C.printNumUnk)
            DisplayUnks(Stokes, lset, adap.GetMG());

        if (ProcCL::IamMaster())
            std::cerr << "==> Solving coupled Levelset-Navier-Stokes problem ....\n";

        cpl.DoStep( C.cpl_iter);

        time.Stop(); duration=time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cerr << "- Solving coupled Levelset-Navier-Stokes problem took "<<duration<<" sec."<<std::endl;
            // Store measured values
            DROPS_LOGGER_SETVALUE("NavStokesCoupledLevelset", duration);
        }

        // Write out solution
        if (C.ensight && step%C.ensight==0)
            ensightwriter.write();
        if (C.vtk && step%C.vtk==0)
            vtkwriter.write();

        // Write out droplet information
        IFInfo.Update(lset, Stokes.GetVelSolution());
        IFInfo.Write(Stokes.t);

        // Reparametrization of levelset function
        if (C.RepFreq && step%C.RepFreq==0)
        {
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Reparametrization\n"
                          << "- rel. Volume: " << relVol << std::endl;

            time.Reset();
            lset.ReparamFastMarching( ExL, C.RepMethod, C.RepMethod==3);
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cerr << "- Reparametrization took "<<duration<<" sec."<<std::endl;
                DROPS_LOGGER_SETVALUE("Reparametrization",duration);
            }
        }
        if (ProcCL::IamMaster())
            std::cerr << "- rel. Volume: " << relVol << std::endl;

        // Correction of the volume
        if (C.VolCorr)
        {
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Adjust volume ...\n";
            time.Reset();
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            lset.Phi.Data+= dphi;
            time.Stop(); duration=time.GetMaxTime();
            relVol = lset.GetVolume()/Vol;
            if (ProcCL::IamMaster()){
                std::cerr << "- Volume correction "<<dphi<<", new rel. Volume is " <<relVol<< " took "<<duration<<" sec"<<std::endl;
                DROPS_LOGGER_SETVALUE("VolumeCorrection",duration);
            }
        }

        // Serialization
        if (C.serialization && step%C.serialization==0)
        {
            if (ProcCL::IamMaster())
                std::cerr << "\n==> Serialize data ...\n";
            time.Reset();
            serializer.Serialize(step);
            time.Stop(); duration=time.GetMaxTime();
            if (ProcCL::IamMaster())
                std::cerr << "- Serialization took " <<duration<< " sec"<<std::endl;
        }

        step_time.Stop();
        duration=step_time.GetMaxTime();
        if (ProcCL::IamMaster()){
            std::cerr <<"========> Step "<<step<<" took "<<duration<<" sec."<<std::endl;
            DROPS_LOGGER_SETVALUE("TotalStep",duration);
            DROPS_LOGGER_NEXTSTEP();
        }
    }
}


template<class Coeff>
  void Strategy( InstatNavierStokes2PhaseP2P1CL<Coeff>& Stokes, ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    typedef InstatNavierStokes2PhaseP2P1CL<Coeff> StokesProblemT;
    MultiGridCL& mg= pmg.GetMG();
    ParTimerCL time;

    // Set parameter of the surface tension
    SurfaceTensionCL::eps           = C.st_jumpWidth;
    SurfaceTensionCL::lambda        = C.st_relPos;
    SurfaceTensionCL::sigma         = Stokes.GetCoeff().SurfTens;
    SurfaceTensionCL::sigma_dirt_fac= C.st_red;

//    instat_vector_fun_ptr gsigma= &(SurfaceTensionCL::grad_sm_step);
    LevelsetP2CL lset( mg, &SurfaceTensionCL::sm_step, &SurfaceTensionCL::grad_sm_step,
                       C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff, C.NarrowBand);

    AdapTriangCL adapt( pmg, lb, C.ref_width, 0, C.ref_flevel);

    if (C.IniCond!=3){
        adapt.MakeInitialTriang(::DistanceFct1);
        if (C.checkMG && !Check( CheckParMultiGrid(pmg)) )
            throw DROPSErrCL("MultiGrid is incorrect!");
    }

    DisplayDetailedGeom(mg);

    LevelsetRepairCL lsetrepair( lset, pmg);
    adapt.push_back( &lsetrepair);
    VelocityRepairCL<StokesProblemT> velrepair( Stokes, pmg);
    adapt.push_back( &velrepair);
    PressureRepairCL<StokesProblemT> prrepair( Stokes, lset, pmg);
    adapt.push_back( &prrepair);

    /// \todo (of) Testen von SF_ImprovedLBVar!!!
    lset.SetSurfaceForce( SF_ImprovedLB);

    std::ofstream *infofile=0;
    if (ProcCL::IamMaster())
        infofile = new std::ofstream( string(C.ensCase + ".info").c_str());

    IFInfo.Init(infofile);

    //Setup initial problem
    if (ProcCL::IamMaster())
        std::cerr << "=================================================================================== Init:\n"
                  << "==> Initialize Problem\n";
    InitProblemWithDrop(Stokes, lset, pmg, lb);

    SolveCoupledNS(Stokes, lset, adapt);

    if (ProcCL::IamMaster()){
        infofile->close();
        delete infofile;
    }
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  DROPS::ProcCL Proc(&argc, &argv);
  try
  {
    if (argc!=2)
    {
        IF_MASTER
          std::cerr << "You have to specify one parameter:\n\t"
                    << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        IF_MASTER
          std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    IF_MASTER
      std::cerr << C << std::endl;

    DROPS::ParTimerCL alltime;
    SetDescriber();
    DROPS::ParMultiGridCL pmg;

    typedef ZeroFlowCL                                    CoeffT;
    typedef DROPS::InstatNavierStokes2PhaseP2P1CL<CoeffT> MyStokesCL;

    int nx, ny, nz;
    double dx, dy, dz;
    std::string mesh( C.meshfile), delim("x@");
    size_t idx;
    while ((idx= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
    if (!brick_info || dx!=dz)
    {
        std::cerr << "error while reading geometry information: " << mesh << "\n";
        return 1;
    }

    C.r_inlet= dx/2;
    DROPS::Point3DCL orig, px, py, pz;
    px[0]= dx; py[1]= dy; pz[2]= dz;


    DROPS::ParTimerCL time;
    DROPS::MGBuilderCL    *mgb     = 0;
    DROPS::BrickBuilderCL *builder = 0;
    if (DROPS::ProcCL::IamMaster()){
        if (C.IniCond!=3){
            mgb = new DROPS::BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
        }
        else{   // read geometry out of a file
            builder = new DROPS::BrickBuilderCL(orig, px, py, pz, nx, ny, nz);
            mgb = new DROPS::FileBuilderCL(C.ser_dir, builder);
        }
    }
    else{
        mgb = new DROPS::EmptyBrickBuilderCL(orig, px, py, pz);
    }

    const bool bc[6]=
      {false, false, true, false, false, false};    // Rohr
//    {false, false, true, false, true, true};      // Kanal
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
      { &Null, &Null, &Null, &Inflow, &Null, &Null};

    DROPS::MultiGridCL mg(*mgb);
    pmg.AttachTo(mg);

    DROPS::LoadBalHandlerCL lb(mg);
    if (C.IniCond!=3)
        lb.DoInitDistribution(DROPS::ProcCL::Master());

    switch (C.refineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
    }


    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();


    MyStokesCL prob(mg, ZeroFlowCL(C), DROPS::StokesBndDataCL( num_bnd, bc, bnd_fun), DROPS::P1_FE, 0);

    Strategy( prob, pmg, lb);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cout);

    // Writeout logging-results
    if (DROPS::ProcCL::IamMaster()){
        // DROPS_LOGGER_WRITEOUT("TestBrickflowPar", (DROPS::ProcCL::Size()==1));

        // build specific name for logfile.
        char refLevel[32];
        sprintf(refLevel,"%d",C.ref_flevel);
        std::string refStr(&refLevel[0]);
        std::string basename("./TestBrickflowPar_Level_");
       // std::string basename("/home/th179319/DropsTimeMeas/results/TestBrickflowPar_Rad3_Level_");
        // std::string basename("/home/th179319/DropsTimeMeas/results/TestBrickflowPar_Level_");
        basename.append(refLevel);
        basename.append("_Proc");
        DROPS_LOGGER_WRITEOUT(basename.c_str(), true);
    }

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

