//**************************************************************************
// File:    MGpoisson.cpp                                                  *
// Content: classes that constitute the poisson-problem with MG-solver     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 16 2001                                         *
//**************************************************************************

#include "num/solver.h"
#include "num/MGsolver.h"

namespace DROPS
{

void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end)
{
    Uint lvl= 1;
    for(const_MGDataIterCL coarse= begin, fine=++begin; fine!=end; ++coarse, ++fine, ++lvl)
    {
//        if (lvl==1) std::cerr << fine->A.Data << std::endl;
        const Uint nu= coarse->Idx.NumUnknowns;
        VectorCL ei(nu), Ai(nu);

        std::cerr << nu << " unknowns on level " << lvl << std::endl;
        
        
        if (nu>700)
        {
            std::cerr << "Check skipped: too many unknowns, too much time..." << std::endl;
            continue;
        }
        
        
        for(Uint i=0; i<nu; ++i)
        {
            if (i!=0) ei[i-1]=0; else ei[nu-1]=0; ei[i]=1;
            Ai= transp_mul( fine->P.Data, fine->A.Data * (fine->P.Data * ei) );
            for(Uint j=0; j<nu; ++j)
            {
                if (fabs(Ai[j] - coarse->A.Data(i,j))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": A_H(" << i << ", " << j <<")= "
                              << coarse->A.Data(i,j) << " != " << Ai[j] << std::endl;
                }
            }
        }    
    }
}

// NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//template <typename Vec> 
void MG( const MGDataCL& MGData, VectorCL& x, const VectorCL& b, 
         int maxiter, double tol )
{
    const_MGDataIterCL finest= --MGData.end();

    Uint   sm    =  2; // how many smoothing steps?
    int    lvl   = -1; // how many levels? (-1=all)
    double omega = 1.; // relaxation parameter for smoother
    
    Uint nit;
    double resid, old_resid;
    
//    JORsmoothCL smoother(omega);  // Jacobi
//    GSsmoothCL smoother(omega);  // Gauss-Seidel
//    SGSsmoothCL smoother(omega);  // symmetric Gauss-Seidel
    SORsmoothCL smoother(omega);  // Gauss-Seidel with over-relaxation
//    SSORsmoothCL smoother(omega);  // symmetric Gauss-Seidel with over-relaxation
    CGSolverCL  solver(200, tol); //CG-Verfahren
    for(Uint k=0; k<sm; ++k)
    {
        std::cerr << "x.size = " << x.size() <<std::endl;
        resid= (b - finest->A.Data * x).norm();
        std::cerr << "initial residuum: " << resid <<std::endl;
	nit = 0;
        do
        {
            MGM( MGData.begin(), finest, x, b, smoother, sm, solver, lvl, -1);
            old_resid= resid;
            resid= (b - finest->A.Data * x).norm();
	    nit = nit+1;
            std::cerr << "iteration: " << nit 
	              << "\tresiduum: " << resid 
	              << "\tred. " << resid/old_resid << std::endl;
        } while ( resid > tol);
    }
}

} // end of namespace DROPS
