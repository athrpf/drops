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
      : rho( JumpCL( C.mat_DensDrop, C.mat_DensFluid ), H_sm, C.mat_SmoothZone),
        mu(  JumpCL( C.mat_ViscDrop,  C.mat_ViscFluid),   H_sm, C.mat_SmoothZone),
        SurfTens( C.sft_SurfTension), g( C.exp_Gravity)    {}
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
      : rho( JumpCL( 1., C.mat_DensFluid/C.mat_DensDrop ), H_sm, C.mat_SmoothZone),
        mu ( JumpCL( 1., C.mat_ViscFluid/C.mat_ViscDrop),    H_sm, C.mat_SmoothZone),
        SurfTens( C.sft_SurfTension/C.mat_DensDrop), g( C.exp_Gravity)    {}
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

