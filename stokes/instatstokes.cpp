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

/*
StokesBndDataCL::bnd_type
bnd_val_e2e3(const Point2DCL& p)
{
    SVectorCL<3> ret;
    Point2DCL q= 2*M_PI*p;
    ret[0]= 0.; ret[1]= -cos(q[0])*sin(q[1]); ret[2]= 2.*sin(q[0])*cos(q[1]);
    return ret;
}

StokesBndDataCL::bnd_type
bnd_val_e1e3(const Point2DCL& p)
{
    SVectorCL<3> ret;
    Point2DCL q= 2*M_PI*p;
    ret[0]= 0.; ret[1]= -cos(q[0])*sin(q[1]); ret[2]= 0.;
    return ret;
}

StokesBndDataCL::bnd_type
bnd_val_e1e2(const Point2DCL& p)
{
    SVectorCL<3> ret;
    Point2DCL q= 2*M_PI*p;
    ret[0]= 0.; ret[1]= 0.; ret[2]= 2.*cos(q[0])*sin(q[1]);
    return ret;
}
*/

    

//==== SchurComplMatrixCL ====

VectorCL operator* (const SchurComplNoPcMatrixCL& M, const VectorCL& v)
{
    double tol= M._tol;
    int maxiter= 1000;
    VectorCL x( M._matA.num_cols());

    CG(M._matA, x, transp_mul(M._matB, v), maxiter, tol);
//    std::cerr << "Inner iteration took " << maxiter << " steps, residuum is " << tol << std::endl;
    return M._matB*x;
}    


} // end of namespace DROPS
