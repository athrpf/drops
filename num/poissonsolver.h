//**************************************************************************
// File:    poissonsolver.h                                                *
// Content: solvers for the Poisson problem                                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - Nov, 20 2001                                           *
//**************************************************************************

#ifndef _POISSONSOLVER_H_
#define _POISSONSOLVER_H_

#include "solver.h"

namespace DROPS
{

//==================================================
//        derived classes for easier use
//==================================================

typedef PCGSolverCL<VectorCL, double, SsorPcCL<VectorCL, double> > 
        PCG_T;

typedef PCGSolverCL<VectorCL, double, ImprovedSsorPcCL<VectorCL, double> > 
        IPCG_T;

typedef PCGSolverCL<VectorCL, double, SGSPcCL<VectorCL, double> > 
        GSPCG_T;


}    // end of namespace DROPS

#endif
