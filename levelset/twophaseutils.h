/// \file twophaseutils.h
/// \brief utilities for twophasedrops.cpp (debug, measuring,...)
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef TWOPHASEUTILS_H_
#define TWOPHASEUTILS_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "misc/problem.h"
#include "levelset/levelset.h"
#ifdef _PAR
#include "parallel/exchange.h"
#endif

namespace DROPS {

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

/// \brief send geometry informations to std::cout
void DisplayDetailedGeom(MultiGridCL& mg);

///\brief Writes all Navier-Stokes matricies to file
template<class StokesT>
void WriteMatrices (StokesT& Stokes, int i)
{
    std::string path( "matrices/");
    std::ostringstream suffix;
    suffix << std::setfill( '0') << std::setw( 4) << i << ".txt";
    WriteToFile( Stokes.A.Data.GetFinest(),   path + "A"   + suffix.str(), "A");
    WriteToFile( Stokes.B.Data.GetFinest(),   path + "B"   + suffix.str(), "B");
    WriteToFile( Stokes.M.Data.GetFinest(),   path + "M"   + suffix.str(), "M");
    WriteToFile( Stokes.prA.Data.GetFinest(), path + "prA" + suffix.str(), "prA");
    WriteToFile( Stokes.prM.Data.GetFinest(), path + "prM" + suffix.str(), "prM");
    WriteToFile( Stokes.N.Data.GetFinest(),   path + "N"   + suffix.str(), "N");

    WriteToFile( Stokes.v.Data, path + "v" + suffix.str(), "v");
    WriteToFile( Stokes.p.Data, path + "p" + suffix.str(), "p");
}

/// For a two-level MG-solver: P2P1 -- P2P1X; canonical prolongations
void MakeP1P1XProlongation (size_t NumUnknownsVel, size_t NumUnknownsPr, size_t NumUnknownsPrP1,
    MatrixCL& PVel, MatrixCL& PPr);

/// For given exact solution u and discrete P2 solution uh on positive and negative subdomains, computes L2 and H1 errors || u - u_h ||
template<class ValT, class DiscSolT>
void ComputeErrorsP2(  ValT (*uPos)(const Point3DCL&, double), ValT (*uNeg)(const Point3DCL&, double), const DiscSolT& uhPos, const DiscSolT& uhNeg, ValT& L2, ValT& H1, const LevelsetP2CL& lset, double t);

void ComputeErrorsP2R( const instat_vector_fun_ptr uPos, const instat_vector_fun_ptr uNeg, const VecDescCL& uh, const BndDataCL<Point3DCL>& bnd, Point3DCL& L2, Point3DCL& H1, const LevelsetP2CL& lset, double t=0.);

// ================ inline/template definitions ==================



inline void AccumulateH1( double& H1, const LocalP2CL<>& diff, const LocalP1CL<Point3DCL> Grad[10], double absdet)
{
    LocalP1CL<Point3DCL> gdiff;
    LocalP2CL<Point3DCL> gdiffLP2;
    P2DiscCL::GetFuncGradient( gdiff, diff, Grad);
    gdiffLP2.assign( gdiff);
    H1+= Quad5CL<>(LocalP2CL<>(dot(gdiffLP2, gdiffLP2))).quad( absdet);
}

inline void AccumulateH1( Point3DCL& H1, const LocalP2CL<Point3DCL>& diff, const LocalP1CL<Point3DCL> Grad[10], double absdet)
{
    LocalP2CL<> diff_comp; // components of diff
    LocalP1CL<Point3DCL> gdiff;
    LocalP2CL<Point3DCL> gdiffLP2;
    for (int k=0; k<3; ++k) {
        ExtractComponent( diff, diff_comp, k);
        P2DiscCL::GetFuncGradient( gdiff, diff_comp, Grad);
        gdiffLP2.assign( gdiff);
        H1[k]+= Quad5CL<>(LocalP2CL<>(dot(gdiffLP2, gdiffLP2))).quad( absdet);
    }
}

inline void AccumulateH1( double& H1, const LocalP2CL<>& diff_p, const LocalP2CL<>& diff_n, InterfaceTetraCL& patch, const LocalP1CL<Point3DCL> Grad[10], double absdet)
{
    LocalP1CL<Point3DCL> gdiff_p, gdiff_n;
    LocalP2CL<Point3DCL> gdiffLP2;

    P2DiscCL::GetFuncGradient( gdiff_p, diff_p, Grad);
    P2DiscCL::GetFuncGradient( gdiff_n, diff_n, Grad);
    for (int ch= 0; ch < 8; ++ch) {
        patch.ComputeCutForChild( ch);
        gdiffLP2.assign( gdiff_p);
        H1+= patch.quad( LocalP2CL<>(dot( gdiffLP2, gdiffLP2)), absdet, true);
        gdiffLP2.assign( gdiff_n);
        H1+= patch.quad( LocalP2CL<>(dot( gdiffLP2, gdiffLP2)), absdet, false);
    }
}

inline void AccumulateH1( Point3DCL& H1, const LocalP2CL<Point3DCL>& diff_p, const LocalP2CL<Point3DCL>& diff_n, InterfaceTetraCL& patch, const LocalP1CL<Point3DCL> Grad[10], double absdet)
{
    LocalP2CL<> diff_comp; // components of diff
    LocalP1CL<Point3DCL> gdiff_p[3], gdiff_n[3];
    LocalP2CL<Point3DCL> gdiffLP2;
    for (int k=0; k<3; ++k) {
        ExtractComponent( diff_p, diff_comp, k);
        P2DiscCL::GetFuncGradient( gdiff_p[k], diff_comp, Grad);
        ExtractComponent( diff_n, diff_comp, k);
        P2DiscCL::GetFuncGradient( gdiff_n[k], diff_comp, Grad);
    }
    for (int ch= 0; ch < 8; ++ch) {
        patch.ComputeCutForChild( ch);
        for (int k=0; k<3; ++k) {
            gdiffLP2.assign( gdiff_p[k]);
            H1[k]+= patch.quad( LocalP2CL<>(dot( gdiffLP2, gdiffLP2)), absdet, true);
            gdiffLP2.assign( gdiff_n[k]);
            H1[k]+= patch.quad( LocalP2CL<>(dot( gdiffLP2, gdiffLP2)), absdet, false);
        }
    }
}

inline void AccumulateH1( Point3DCL& H1, const LocalP2CL<Point3DCL>& diff_p, const LocalP2CL<Point3DCL>& diff_n, const LocalP2CL<Point3DCL> ext_p[8], const LocalP2CL<Point3DCL> ext_n[8],
                          InterfaceTetraCL& patch, const LocalP1CL<Point3DCL> Grad[10], double absdet)
{
    LocalP2CL<> comp; // component
    LocalP1CL<Point3DCL> gdiff_p[3], gdiff_n[3], gext;
    LocalP2CL<Point3DCL> gdiffLP2;
    for (int k=0; k<3; ++k) {
        ExtractComponent( diff_p, comp, k);
        P2DiscCL::GetFuncGradient( gdiff_p[k], comp, Grad);
        ExtractComponent( diff_n, comp, k);
        P2DiscCL::GetFuncGradient( gdiff_n[k], comp, Grad);
    }
    for (int ch= 0; ch < 8; ++ch) {
        patch.ComputeCutForChild( ch);
        for (int k=0; k<3; ++k) {
            ExtractComponent( ext_p[ch], comp, k);
            P2DiscCL::GetFuncGradient( gext, comp, Grad);
            gdiffLP2.assign( LocalP1CL<Point3DCL>(gdiff_p[k] - gext));
            H1[k]+= patch.quad( LocalP2CL<>(dot( gdiffLP2, gdiffLP2)), absdet, true);

            ExtractComponent( ext_n[ch], comp, k);
            P2DiscCL::GetFuncGradient( gext, comp, Grad);
            gdiffLP2.assign(  LocalP1CL<Point3DCL>(gdiff_n[k] - gext));
            H1[k]+= patch.quad( LocalP2CL<>(dot( gdiffLP2, gdiffLP2)), absdet, false);
        }
    }
}

template<class ValT, class DiscSolT>
void ComputeErrorsP2(  ValT (*uPos)(const Point3DCL&, double), ValT (*uNeg)(const Point3DCL&, double), const DiscSolT& uhPos, const DiscSolT& uhNeg, ValT& L2, ValT& H1, const LevelsetP2CL& lset, double t)
{
    const int lvl= uhPos.GetLevel();
    const MultiGridCL& MG= lset.GetMG();
    LocalP2CL<ValT> u_p, u_n, uh_p, uh_n, diff, diff_p, diff_n;
    Quad5CL<ValT> u5_p, u5_n, diff5;
    LocalP1CL<Point3DCL> Grad[10], GradRef[10];
    LocalP2CL<> loc_phi;
    BaryCoordCL *nodes;
    double det, absdet;
    SMatrixCL<3,3> T;
    InterfaceTetraCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();

    L2= H1= ValT();

    for (MultiGridCL::const_TriangTetraIteratorCL sit = MG.GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= std::fabs( det);
        P2DiscCL::GetGradients( Grad, GradRef, T);

        loc_phi.assign( *sit, ls, t);
        patch.Init( *sit, loc_phi);
        const bool nocut= !patch.Intersects();
        if (nocut) {
            if (patch.GetSign( 0) == 1) { // pos. part
                u_p.assign(  *sit, uPos, t);
                u5_p.assign(  *sit, uPos, t);
                uh_p.assign( *sit, uhPos, t);
                diff= u_p - uh_p;
                diff5= u5_p - Quad5CL<ValT>(uh_p);
            } else { // neg. part
                u_n.assign(  *sit, uNeg, t);
                u5_n.assign( *sit, uNeg, t);
                uh_n.assign( *sit, uhNeg, t);
                diff= u_n - uh_n;
                diff5= u5_n - Quad5CL<ValT>(uh_n);
            }
            L2+= Quad5CL<ValT>(diff5*diff5).quad( absdet);

            AccumulateH1( H1, diff, Grad, absdet);
        }
        else { // We are at the phase boundary.
            u_p.assign( *sit, uPos, t);
            u_n.assign( *sit, uNeg, t);
            uh_p.assign( *sit, uhPos, t);
            uh_n.assign( *sit, uhNeg, t);
            diff_p= u_p - uh_p;
            diff_n= u_n - uh_n;
            // compute L2 norm
            patch.ComputeSubTets();
            for (Uint k=0; k<patch.GetNumTetra(); ++k)
            { // init quadrature objects for integrals on sub tetras
                nodes = Quad5CL<ValT>::TransformNodes(patch.GetTetra(k));
                Quad5CL<ValT> diff5cut( k<patch.GetNumNegTetra() ? diff_n : diff_p, nodes);
                L2+= Quad5CL<ValT>( diff5cut*diff5cut).quad(absdet*VolFrac(patch.GetTetra(k)));
                delete[] nodes;
            }

            AccumulateH1( H1, diff_p, diff_n, patch, Grad, absdet);
        }
    }
    H1+= L2;
    L2= sqrt( L2);
    H1= sqrt( H1);
}

} //end of namespace DROPS

#endif /* TWOPHASEUTILS_H_ */
