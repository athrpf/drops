/// \file twophaseutils.cpp
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

#include "levelset/twophaseutils.h"
#include "num/discretize.h"

namespace DROPS {

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

/// For a two-level MG-solver: P2P1 -- P2P1X; canonical prolongations
void MakeP1P1XProlongation (size_t NumUnknownsVel, size_t NumUnknownsPr, size_t NumUnknownsPrP1,
    MatrixCL& PVel, MatrixCL& PPr)
{
    // finest level
    //P2-Prolongation (Id)
    PVel= MatrixCL( std::valarray<double>(  1.0, NumUnknownsVel));
    //P1-P1X-Prolongation
    VectorCL diag( 0., NumUnknownsPr);
    diag[std::slice(0, NumUnknownsPrP1, 1)]= 1.;
    PPr= MatrixCL( diag);
}

void ComputeErrorsP2R( const instat_vector_fun_ptr uPos, const instat_vector_fun_ptr uNeg, const VecDescCL& u_h, const BndDataCL<Point3DCL>& bnd, Point3DCL& L2, Point3DCL& H1, const LevelsetP2CL& lset, double t)
{
    const IdxDescCL& idx= *u_h.RowIdx;
    const ExtIdxDescCL& extIdx= idx.GetXidx();
    const int lvl= idx.TriangLevel();
    const MultiGridCL& MG= lset.GetMG();

    LocalP2CL<Point3DCL> u_p, u_n, uh, diff, diff_p, diff_n;
    Quad5CL<Point3DCL> u5_p, u5_n, diff5;
    LocalP1CL<Point3DCL> gdiff, gdiff_p, gdiff_n;
    LocalP2CL<Point3DCL> gdiffLP2;
    LocalP1CL<Point3DCL> Grad[10], GradRef[10];
    LocalP2CL<> loc_phi;
    LocalP2CL<> p1r_p[4][8], p1r_n[4][8];    // extended basis functions per child
    LocalP2CL<Point3DCL> ext_p[8], ext_n[8]; // extended part per child
    BaryCoordCL *nodes;
    double det, absdet;
    Point3DCL xval;
    SMatrixCL<3,3> T;
    LocalNumbP1CL numb; // std P1 numbering
    IdxT xnr[4];        // extended numbering
    InterfaceTetraCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();

    L2= H1= Point3DCL();

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
                uh.assign( *sit, u_h, bnd);
                diff= u_p - uh;
                diff5= u5_p - Quad5CL<Point3DCL>(uh);
            } else { // neg. part
                u_n.assign(  *sit, uNeg, t);
                u5_n.assign( *sit, uNeg, t);
                uh.assign( *sit, u_h, bnd);
                diff= u_n - uh;
                diff5= u5_n - Quad5CL<Point3DCL>(uh);
            }
            L2+= Quad5CL<Point3DCL>(diff5*diff5).quad( absdet);

            AccumulateH1( H1, diff, Grad, absdet);
        }
        else { // We are at the phase boundary.
            u_p.assign( *sit, uPos, t);
            u_n.assign( *sit, uNeg, t);
            uh.assign( *sit, u_h, bnd);
            diff_p= u_p - uh;
            diff_n= u_n - uh;
            // compute extended part
            P2RidgeDiscCL::GetExtBasisOnChildren( p1r_p, p1r_n, loc_phi);
            numb.assign( *sit, idx, idx.GetBndInfo());
            for (int v=0; v<4; ++v) {
                xnr[v]= numb.num[v]!=NoIdx ? extIdx[numb.num[v]] : NoIdx;
            }
            for (int ch=0; ch<8; ++ch) {
                ext_p[ch]= Point3DCL();
                ext_n[ch]= Point3DCL();
                for (int v=0; v<4; ++v)
                    if (xnr[v] != NoIdx) {
                        for (int k=0; k<3; ++k)
                            xval[k]= u_h.Data[xnr[v]+k];
                        ext_p[ch]+= xval * p1r_p[v][ch];
                        ext_n[ch]+= xval * p1r_n[v][ch];
                    }
            }
            // compute L2 norm
            patch.ComputeSubTets();
            for (Uint k=0; k<patch.GetNumTetra(); ++k)
            { // init quadrature objects for integrals on sub tetras
                nodes = Quad5CL<Point3DCL>::TransformNodes(patch.GetTetra(k));
                const int ch= patch.GetChildIdx(k);
                Quad5CL<Point3DCL> diff5cut( k<patch.GetNumNegTetra() ? LocalP2CL<Point3DCL>(diff_n - ext_n[ch])
                                                                      : LocalP2CL<Point3DCL>(diff_p - ext_p[ch]), nodes);
                L2+= Quad5CL<Point3DCL>( diff5cut*diff5cut).quad(absdet*VolFrac(patch.GetTetra(k)));
                delete[] nodes;
            }

            AccumulateH1( H1, diff_p, diff_n, ext_p, ext_n, patch, Grad, absdet);
        }
    }
    H1+= L2;
    L2= sqrt( L2);
    H1= sqrt( H1);
}


}// end of namespace DROPS
