//**************************************************************************
// File:    instatstokes.cpp                                               *
// Content: classes that constitute the stokes-problem                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Sep, 13 2001                                           *
//**************************************************************************

#include "stokes/instatstokes.h"
#include "num/discretize.h"

namespace DROPS
{


//==== SchurComplMatrixCL ====

VectorCL operator* (const SchurComplNoPcMatrixCL& M, const VectorCL& v)
{
    double tol= M._tol;
    int maxiter= 1000;
    VectorCL x( M._matA.num_cols());

    CG(M._matA, x, transp_mul(M._matB, v), maxiter, tol);
    if (maxiter > 990)
        Comment(     "VectorCL operator* (const SchurComplNoPcMatrixCL& M, const VectorCL& v): "
                  << "Needed more than 990 iterations! tol: " << tol << std::endl,
                  DebugNumericC);
//    std::cerr << "Inner iteration took " << maxiter << " steps, residuum is " << tol << std::endl;
    return M._matB*x;
}    


} // end of namespace DROPS
