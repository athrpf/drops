//**************************************************************************
// File:    stokes.cpp                                                     *
// Content: classes that constitute the stokes-problem                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

#ifndef _STOKES_CPP_
#define _STOKES_CPP_

#include "stokes/stokes.h"
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
/*
double GradPrEstimator(const TetraCL& t, const DiscPrSolCL& p, const DiscVelSolCL&)
{
    double diff[3];
    double p0= p.val( t->GetVertex(0));
    for (Uint i=0; i<3; ++i)
        diff[i]= fabs( p.val( t->GetVertex(i+1)) -  p0);
        
    return std::max( diff) * t->GetVolume();    
}
*/

    
void Uzawa(const MatrixCL& A, const MatrixCL& B, const MatrixCL& M, VectorCL& x, VectorCL& y, const VectorCL& f, const VectorCL& g, 
           double tau, int& max_iter, double& tol, Uint inner_iter, double inner_iter_tol)
{
    VectorCL x_corr(x.size()),
             y_corr(y.size()),
             res1(x.size()),
             res2(y.size());
    SsorPcCL<VectorCL, double> pc(1);

    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm, res2_norm;
    for( int step=0; step<max_iter; ++step)
    {
        int inner_max_iter= inner_iter;
        double inner_tol= inner_iter_tol;

        PCG( M, y_corr, res2= B*x - g, pc, inner_max_iter, inner_tol);
        
        y+= tau * y_corr;
        res1= A*x + transp_mul(B,y) - f;

        res1_norm= res1.norm2();
        res2_norm= res2.norm2();
        if (res1_norm + res2_norm < tol)
        {
            tol= ::sqrt( res1_norm + res2_norm );
            max_iter= step;
            return;
        }   
        if( (step%output)==0 )
            std::cerr << "step " << step << ": norm of 1st eq= " << ::sqrt( res1_norm)
                      << ", norm of 2nd eq= " << ::sqrt( res2_norm) << std::endl;

        inner_max_iter= inner_iter;
        inner_tol= inner_iter_tol;  
        
        PCG( A, x_corr, res1, pc, inner_max_iter, inner_tol);
        x-= x_corr;
    }
    tol= ::sqrt( res1_norm + res2_norm );
}
/*
void Uzawa2(const MatrixCL& A, const MatrixCL& B, VectorCL& x, VectorCL& y, const VectorCL& f, const VectorCL& g, 
            int& max_iter, double& tol, Uint inner_iter, double inner_iter_tol)
{
    VectorCL x_corr(x.size()),
             y_corr(y.size());
    SsorPcCL<VectorCL, double> pc(1);
    SchurComplMatrixCL H (A, B, inner_iter_tol, 1);

    tol*= tol;
    Uint output= 100;//max_iter/20;  // nur 20 Ausgaben pro Lauf
    
    for(int step=0; step<max_iter; ++step)
    {
        int inner_max_iter= inner_iter;
        double inner_tol= 1e-20;
        PCG( A, x_corr, f - (A*x + transp_mul(B,y)), pc, inner_max_iter, inner_tol);
        x+= x_corr;

        inner_max_iter= inner_iter;
        inner_tol     = 1e-20;
        CG( H, y_corr, B*x - g, inner_max_iter, inner_tol);
        y+= y_corr;
    
        double x_corr_norm= x_corr.norm2(),
               y_corr_norm= y_corr.norm2();
        if( (step%output)==0 )
            std::cerr << "step: " << step << ": norm_sq of x_corr = " << x_corr_norm
                      << ", norm_sq of y_corr = " << y_corr_norm << std::endl;
        if (x_corr_norm + y_corr_norm < tol)
        {
            tol= ::sqrt( x_corr_norm + y_corr_norm );
            max_iter= step;
            return;
        }   
    }
    tol= ::sqrt( x_corr.norm2() + y_corr.norm2() );
}

void Uzawa3(const MatrixCL& A, const MatrixCL& B, VectorCL& x, VectorCL& y, const VectorCL& f, const VectorCL& g, 
            int& max_iter, double& tol, Uint inner_iter, double inner_iter_tol)
{
    VectorCL r(x.size()),
             s(y.size());
    VectorCL v(x.size()),
             c(y.size()), d(y.size());
    SsorPcCL<VectorCL, double> pc(1);
    SchurComplMatrixCL H (A, B, inner_iter_tol, 1);

    tol*= tol;
    Uint output= 100; // nach jeder 100. Iteration eine Ausgabe
         //max_iter/20;  // nur 20 Ausgaben pro Lauf
    
    for(int step=0; step<max_iter; ++step)
    {
        r= f - A*x - transp_mul(B,y);
        s= B*x - g;
        
        int inner_max_iter= inner_iter;
        double inner_tol= 1e-20;
        PCG( A, v, r, pc, inner_max_iter, inner_tol);

        inner_max_iter= inner_iter;
        inner_tol     = 1e-20;
        CG( H, d, s + B*v, inner_max_iter, inner_tol);
    
        inner_max_iter= inner_iter;
        inner_tol= 1e-20;
        PCG( A, v, r - transp_mul(B,d), pc, inner_max_iter, inner_tol);
        
        x+= v;
        y+= d;

        double x_corr_norm= v.norm2(),
               y_corr_norm= d.norm2();
        if( (step%output)==0 )
            std::cerr << "step: " << step << ": norm_sq of x_corr = " << x_corr_norm
                      << ", norm_sq of y_corr = " << y_corr_norm << std::endl;
        if (x_corr_norm + y_corr_norm < tol)
        {
            tol= ::sqrt( x_corr_norm + y_corr_norm );
            max_iter= step;
            return;
        }   
    }
    tol= ::sqrt( v.norm2() + d.norm2() );
}
*/


//==== SchurComplMatrixCL ====

VectorCL operator* (const SchurComplMatrixCL& M, const VectorCL& v)
{
    double tol= M._tol;
    int maxiter= 100;
    VectorCL x( M._matA.num_cols());
    
    PCG(M._matA, x, transp_mul(M._matB, v), M._pc, maxiter, tol);
//    std::cerr << "Inner iteration took " << maxiter << " steps, residuum is " << tol << std::endl;
    return M._matB*x;
}    

} // end of namespace DROPS

#endif



