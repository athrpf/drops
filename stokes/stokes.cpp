//**************************************************************************
// File:    stokes.cpp                                                     *
// Content: classes that constitute the stokes-problem                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

#include "stokes/stokes.h"
#include "num/discretize.h"

namespace DROPS
{

    
void Uzawa(const MatrixCL& A, const MatrixCL& B, const MatrixCL& M, VectorCL& x, VectorCL& y, const VectorCL& f, const VectorCL& g, 
           double tau, int& max_iter, double& tol, Uint inner_iter, double inner_iter_tol)
{
    VectorCL x_corr(x.size()),
             y_corr(y.size()),
             res1(x.size()),
             res2(y.size());
    SSORPcCL pc(1);

    tol*= tol;
    Uint output= 50;//max_iter/20;  // nur 20 Ausgaben pro Lauf

    double res1_norm= 0., res2_norm= 0.;
    for( int step=0; step<max_iter; ++step)
    {
        int inner_max_iter= inner_iter;
        double inner_tol= inner_iter_tol;

//        PCG( M, y_corr, res2= B*x - g, pc, inner_max_iter, inner_tol);
        z_xpay(res2, B*x, -1.0, g);
        PCG( M, y_corr, res2, pc, inner_max_iter, inner_tol);

//        y+= tau * y_corr;
        axpy(tau, y_corr, y);
//        res1= A*x + transp_mul(B,y) - f;
        z_xpaypby2(res1, A*x, 1.0, transp_mul(B,y), -1.0, f);

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

} // end of namespace DROPS
