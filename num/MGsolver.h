//**************************************************************************
// File:    MGpoisson.h                                                    *
// Content: classes that constitute the poisson-problem with MG-solver     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 16 2001                                         *
//**************************************************************************

#ifndef DROPS_MGSOLVER_H
#define DROPS_MGSOLVER_H

#include "misc/problem.h"
#include <list>


namespace DROPS
{



struct MGLevelDataCL  // data for one triang level
{
    IdxDescCL Idx;    // index description
    MatDescCL A, P;   // stiffness matrix / prolongation
};

typedef std::list<MGLevelDataCL> MGDataCL;
typedef MGDataCL::iterator       MGDataIterCL;
typedef MGDataCL::const_iterator const_MGDataIterCL;

template<class SmootherCL, class DirectSolverCL>
void MGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b, 
          const SmootherCL& Smoother, Uint smoothSteps, 
          DirectSolverCL& Solver, int numLevel, int numUnknDirect);
// Multigrid method, V-cycle, beginning from level 'fine' 
// numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
//      after 'numLevel' visited levels or if #Unknowns <= 'numUnknDirect'
// If one of the parameters is -1, it will be neglected. 
// If the coarsest level 'begin' has been reached, the direct solver is used too.
// NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)
    
void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end);
// checks  A_coarse= PT * A_fine * P on each level

//===================================
// definition of template functions
//===================================

template<class SmootherCL, class DirectSolverCL>
void MGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b, 
          const SmootherCL& Smoother, Uint smoothSteps, 
          DirectSolverCL& Solver, int numLevel, int numUnknDirect)
// Multigrid method, V-cycle. If numLevel==0 or #Unknowns <= numUnknDirect, the direct solver Solver is used.
// If one of the parameters is -1, it will be neglected. If MGData.begin() has been reached, the direct solver is used too.
{
    const_MGDataIterCL coarse= fine;
    --coarse;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) ) 
       || fine==begin)
    { // use direct solver
        Solver.Solve( fine->A.Data, x, b);
        return;
    }
    VectorCL d(coarse->Idx.NumUnknowns), e(coarse->Idx.NumUnknowns);
    // presmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( fine->A.Data, x, b);
    // restriction of defect
    d= transp_mul( fine->P.Data, b - fine->A.Data*x );
    // calculate coarse grid correction
    MGM( begin, coarse, e, d, Smoother, smoothSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    // add coarse grid correction
    x+= fine->P.Data * e;
    // postsmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( fine->A.Data, x, b);
}


} // end of namespace DROPS

#endif
