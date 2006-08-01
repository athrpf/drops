//**************************************************************************
// File:    MGpoisson.cpp                                                  *
// Content: classes that constitute the poisson-problem with MG-solver     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 16 2001                                         *
//**************************************************************************

#include "num/MGsolver.h"

namespace DROPS
{

void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end)
{
    Uint lvl= 0;
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
                if (std::fabs(Ai[j] - coarse->A.Data(i,j))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": A_H(" << i << ", " << j <<")= "
                              << coarse->A.Data(i,j) << " != " << Ai[j] << std::endl;
                }
            }
        }
    }
}


//template <typename Vec>
void MG(const MGDataCL& MGData, VectorCL& x, const VectorCL& b,
        int& maxiter, double& tol, const bool residerr)
{
    const_MGDataIterCL finest= --MGData.end();
    Uint   sm   =  1; // how many smoothing steps?
    int    lvl  = -1; // how many levels? (-1=all)
    double omega= 1.; // relaxation parameter for smoother
    double resid= -1;
    double old_resid;
    VectorCL tmp;
    if (residerr == true) {
        resid= norm( b - finest->A.Data * x);
        std::cerr << "initial residual: " << resid << '\n';
    }
    else
        tmp.resize( x.size());

//    JORsmoothCL smoother( omega); // Jacobi
//    GSsmoothCL smoother( omega); // Gauss-Seidel
//    SGSsmoothCL smoother( omega); // symmetric Gauss-Seidel
//    SORsmoothCL smoother( omega); // Gauss-Seidel with over-relaxation
    SSORsmoothCL smoother( omega); // symmetric Gauss-Seidel with over-relaxation
//    CGSolverCL  solver( 200, tol); //CG-Verfahren
    SSORPcCL directpc; PCG_SsorCL solver( directpc, 200, tol);
    int it;
    for (it= 0; it<maxiter; ++it) {
        if (residerr == true) {
            if (resid <= tol) break;
        }
        else tmp= x;
        MGM( MGData.begin(), finest, x, b, smoother, sm, solver, lvl, -1);
        if (residerr == true) {
            old_resid= resid;
            resid= norm( b - finest->A.Data * x);
//            std::cerr << "iteration: " << it  << "\tresidual: " << resid;
//            std::cerr << "\treduction: " << resid/old_resid;
//            std::cerr << '\n';
        }
        else if ((resid= norm( tmp - x)) <= tol) break;
    }
    maxiter= it;
    tol= resid;
}

} // end of namespace DROPS
