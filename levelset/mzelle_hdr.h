/// \file mzelle_hdr.h
/// \brief header for all main programs simulating the measuring cell
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_MZELLE_HDR_H
#define DROPS_MZELLE_HDR_H

#include "geom/multigrid.h"
#include "levelset/params.h"
#include "levelset/levelset.h"
#include "num/discretize.h"
#include "stokes/instatstokes2phase.h"
#include "poisson/transport2phase.h"
#include "poisson/poisson.h"
#include "geom/geomselect.h"
#include "num/nssolver.h"
#include "levelset/coupling.h"

namespace DROPS
{

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
    static Point3DCL& GetCenter() { return Mitte_; }
    static Point3DCL& GetRadius() { return Radius_; }
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
    double maxGrad, Vol, h_min, h_max, surfArea, sphericity;


    template<class DiscVelSolT>
    void Update (const LevelsetP2CL& ls, const DiscVelSolT& u) {
        ls.GetInfo( maxGrad, Vol, bary, vel, u, min, max, surfArea);
        std::pair<double, double> h= h_interface( ls.GetMG().GetTriangEdgeBegin( ls.Phi.RowIdx->TriangLevel()), ls.GetMG().GetTriangEdgeEnd( ls.Phi.RowIdx->TriangLevel()), ls.Phi);
        h_min= h.first; h_max= h.second;
        // sphericity is the ratio of surface area of a sphere of same volume and surface area of the approximative interface
        sphericity= std::pow(6*Vol, 2./3.)*std::pow(M_PI, 1./3.)/surfArea;
    }
    void WriteHeader() {
        if (file_)
          (*file_) << "# time maxGradPhi volume bary_drop min_drop max_drop vel_drop h_min h_max surfArea sphericity" << std::endl;
    }
    void Write (double time) {
        if (file_)
          (*file_) << time << " " << maxGrad << " " << Vol << " " << bary << " " << min << " " << max << " " << vel << " " << h_min << " " << h_max << " " << surfArea << " " << sphericity << std::endl;
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


/// \brief factory for the time discretization schemes
template<class LevelSetSolverT>
TimeDisc2PhaseCL* CreateTimeDisc( InstatNavierStokes2PhaseP2P1CL& Stokes, LevelsetP2CL& lset,
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* stokessolver, LevelSetSolverT* lsetsolver, ParamMesszelleNsCL& C, LevelsetModifyCL& lsetmod)
{
    if (C.tm_NumSteps == 0) return 0;
    switch (C.tm_Scheme)
    {
        case 1 :
            return (new LinThetaScheme2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.stk_Theta, C.lvs_Theta, C.ns_Nonlinear, C.cpl_Stab));
        break;
        case 3 :
            std::cout << "[WARNING] use of ThetaScheme2PhaseCL is deprecated using RecThetaScheme2PhaseCL instead\n";
        case 2 :
            return (new RecThetaScheme2PhaseCL<LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.stk_Theta, C.lvs_Theta, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab));
        break;
        case 4 :
            return (new OperatorSplitting2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, stokessolver->GetStokesSolver(), *lsetsolver, lsetmod, C.stk_InnerIter, C.stk_InnerTol, C.ns_Nonlinear));
        break;
        case 6 :
            return (new SpaceTimeDiscTheta2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.stk_Theta, C.lvs_Theta, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab, false));
        break;
        case 7 :
            return (new SpaceTimeDiscTheta2PhaseCL< LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.stk_Theta, C.lvs_Theta, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab, true));
        break;
        case 8 :
            return (new EulerBackwardScheme2PhaseCL<LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab));
        break;
        case 9 :
            return (new CrankNicolsonScheme2PhaseCL<RecThetaScheme2PhaseCL, LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab));
        break;
        case 10 :
            return (new CrankNicolsonScheme2PhaseCL<SpaceTimeDiscTheta2PhaseCL, LevelSetSolverT>
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab));
        break;
        case 11 :
            return (new FracStepScheme2PhaseCL<RecThetaScheme2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab, -1));
        break;
        case 12 :
            return (new FracStepScheme2PhaseCL<SpaceTimeDiscTheta2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab, -1));
        break;
        case 13 :
            return (new Frac2StepScheme2PhaseCL<RecThetaScheme2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab, -1));
        break;
        case 14 :
            return (new Frac2StepScheme2PhaseCL<SpaceTimeDiscTheta2PhaseCL, LevelSetSolverT >
                        (Stokes, lset, *stokessolver, *lsetsolver, lsetmod, C.cpl_Tol, C.ns_Nonlinear, C.cpl_Projection, C.cpl_Stab, -1));
        break;
        default : throw DROPSErrCL("Unknown TimeDiscMethod");
    }
}

template <class StokesT>
void SolveStatProblem( StokesT& Stokes, const LevelsetP2CL& lset,
                       NSSolverBaseCL<StokesT >& solver)
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;
    time.Reset();
    VelVecDescCL cplM, cplN;
    VecDescCL curv;
    cplM.SetIdx( &Stokes.vel_idx);
    cplN.SetIdx( &Stokes.vel_idx);
    curv.SetIdx( &Stokes.vel_idx);
    Stokes.SetIdx();
    Stokes.SetLevelSet( lset);
    lset.AccumulateBndIntegral( curv);
    Stokes.SetupSystem1( &Stokes.A, &Stokes.M, &Stokes.b, &Stokes.b, &cplM, lset, Stokes.v.t);
    Stokes.SetupPrStiff( &Stokes.prA, lset);
    Stokes.SetupPrMass ( &Stokes.prM, lset);
    Stokes.SetupSystem2( &Stokes.B, &Stokes.c, lset, Stokes.v.t);
    time.Stop();
    duration = time.GetTime();
    std::cout << "Discretizing took "<< duration << " sec.\n";
    time.Reset();
    Stokes.b.Data += curv.Data;
    solver.Solve( Stokes.A.Data, Stokes.B.Data, Stokes.v, Stokes.p.Data, Stokes.b.Data, cplN, Stokes.c.Data, 1.0);
    time.Stop();
    duration = time.GetTime();
    std::cout << "Solving (Navier-)Stokes took "<<  duration << " sec.\n";
    std::cout << "iter: " << solver.GetIter() << "\tresid: " << solver.GetResid() << std::endl;
}

void SetInitialLevelsetConditions( LevelsetP2CL& lset, MultiGridCL& MG, const ParamMesszelleNsCL& C)
{
    switch (C.dmc_InitialCond)
    {
#ifndef _PAR
      case -10: // read from ensight-file [deprecated]
      {
        std::cout << "[DEPRECATED] read from ensight-file [DEPRECATED]\n";
        ReadEnsightP2SolCL reader( MG);
        reader.ReadScalar( C.dmc_InitialFile+".scl", lset.Phi, lset.GetBndData());
      } break;
#endif
      case -1: // read from file
      {
        ReadFEFromFile( lset.Phi, MG, C.dmc_InitialFile+"levelset", C.rst_Binary);
      } break;
      case 0: case 1:// zero initial condition // stationary flow
          lset.Init( EllipsoidCL::DistanceFct);
        break;
      case  2: //flow without droplet
          lset.Init( &One);
      break;
      default : throw DROPSErrCL("Unknown initial condition");
    }
}

template <typename StokesT>
void SetInitialConditions(StokesT& Stokes, const LevelsetP2CL& lset, MultiGridCL& MG, const ParamMesszelleNsCL& C)
{
    MLIdxDescCL* pidx= &Stokes.pr_idx;
    switch (C.dmc_InitialCond)
    {
#ifndef _PAR
      case -10: // read from ensight-file [deprecated]
      {
        std::cout << "[DEPRECATED] read from ensight-file [DEPRECATED]\n";
        ReadEnsightP2SolCL reader( MG);
        reader.ReadVector( C.dmc_InitialFile+".vel", Stokes.v, Stokes.GetBndData().Vel);
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        if (Stokes.UsesXFEM()) {
            VecDescCL pneg( pidx), ppos( pidx);
            reader.ReadScalar( C.dmc_InitialFile+".prNeg", pneg, Stokes.GetBndData().Pr);
            reader.ReadScalar( C.dmc_InitialFile+".prPos", ppos, Stokes.GetBndData().Pr);
            P1toP1X ( pidx->GetFinest(), Stokes.p.Data, pidx->GetFinest(), ppos.Data, pneg.Data, lset.Phi, MG);
        }
        else
            reader.ReadScalar( C.dmc_InitialFile+".pr", Stokes.p, Stokes.GetBndData().Pr);
      } break;
#endif
      case -1: // read from file
      {
        ReadFEFromFile( Stokes.v, MG, C.dmc_InitialFile+"velocity", C.rst_Binary);
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
        ReadFEFromFile( Stokes.p, MG, C.dmc_InitialFile+"pressure", C.rst_Binary, &lset.Phi); // pass also level set, as p may be extended
      } break;
      case 0: // zero initial condition
          Stokes.UpdateXNumbering( pidx, lset);
          Stokes.p.SetIdx( pidx);
        break;
      case 1: // stationary flow
      {
        Stokes.UpdateXNumbering( pidx, lset);
        Stokes.p.SetIdx( pidx);
#ifdef _PAR
        ParJac0CL jacpc( Stokes.vel_idx.GetFinest());
        typedef ParPCGSolverCL<ParJac0CL> PCGSolverT;
        typedef SolverAsPreCL<PCGSolverT> PCGPcT;
        PCGSolverT PCGSolver(200, 1e-2, Stokes.vel_idx.GetFinest(), jacpc, /*rel*/ true, /*acc*/ true);
        PCGPcT     apc(PCGSolver);
        ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), Stokes.vel_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
        ParInexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM> inexactuzawasolver( apc, bbtispc, Stokes.vel_idx.GetFinest(), Stokes.pr_idx.GetFinest(),
                                                                           C.stk_OuterIter, C.stk_OuterTol, 0.6, 50, &std::cout);
#else
        SSORPcCL ssorpc;
        PCG_SsorCL PCGsolver( ssorpc, 200, 1e-2, true);
        typedef SolverAsPreCL<PCG_SsorCL> PCGPcT;
        PCGPcT apc( PCGsolver);
        ISBBTPreCL bbtispc( &Stokes.B.Data.GetFinest(), &Stokes.prM.Data.GetFinest(), &Stokes.M.Data.GetFinest(), Stokes.pr_idx.GetFinest(), 0.0, 1.0, 1e-4, 1e-4);
        InexactUzawaCL<PCGPcT, ISBBTPreCL, APC_SYM> inexactuzawasolver( apc, bbtispc, C.stk_OuterIter, C.stk_OuterTol, 0.6, 50);
#endif
        NSSolverBaseCL<StokesT> stokessolver( Stokes, inexactuzawasolver);
        SolveStatProblem( Stokes, lset, stokessolver);
      } break;
      case  2: //flow without droplet
          Stokes.UpdateXNumbering( pidx, lset);
          Stokes.p.SetIdx( pidx);
      break;
      default : throw DROPSErrCL("Unknown initial condition");
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
    bool                 binary_;

    /// \brief Write time info
    void WriteTime( std::string filename)
    {
        std::ofstream file( filename.c_str());
        if (!file) throw DROPSErrCL("TwoPhaseStoreCL::WriteTime: Cannot open file "+filename+" for writing");
        file << Stokes_.v.t << "\n";
        file.close();
    }

  public:
      /// \brief Construct a class for storing a two-phase flow problem in files
      /** This class generates multiple files, all with prefix path, for storing
       *  the geometric as well as the numerical data.
       *  \param recoverySteps number of backup steps before overwriting files
       *  */
    TwoPhaseStoreCL(MultiGridCL& mg, const StokesT& Stokes, const LevelsetP2CL& lset, const TransportP1CL* transp,
                    const std::string& path, Uint recoverySteps=2, bool binary= false)
      : mg_(mg), Stokes_(Stokes), lset_(lset), transp_(transp), path_(path), numRecoverySteps_(recoverySteps),
        recoveryStep_(0), binary_( binary) {}

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
        WriteFEToFile(Stokes_.v, mg_, filename.str() + "velocity", binary_);
        WriteFEToFile(lset_.Phi, mg_, filename.str() + "levelset", binary_);
        WriteFEToFile(Stokes_.p, mg_, filename.str() + "pressure", binary_, &lset_.Phi); // pass also level set, as p may be extended
        if (transp_) WriteFEToFile(transp_->ct, mg_, filename.str() + "concentrationTransf", binary_);
    }
};

}   // end of namespace DROPS

#endif

