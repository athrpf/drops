//**************************************************************************
// File:    stokessolver.cpp                                               *
// Content: solvers for the Stokes problem                                 *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 15 2004                                           *
//**************************************************************************


#include "num/stokessolver.h"

namespace DROPS
{

void
InexactUzawa_CL::Solve( const MatrixCL& A, const MatrixCL& B, VectorCL& v, VectorCL& p,
    const VectorCL& b, const VectorCL& c)
{
    _res=  _tol;
    _iter= _maxiter;
    InexactUzawa( A, B, v, p, b, c, Apc_, Spc_, _iter, _res);
}


} // end of namespace DROPS
