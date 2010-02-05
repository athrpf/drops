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


}// end of namespace DROPS
