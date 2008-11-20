/// \file
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Maxim Larin, Patrick Esser, IGPM

#include "num/MGsolver.h"

namespace DROPS
{

void CheckMGData( const MLMatrixCL& A, const MLMatrixCL& P)
{
    Uint lvl= 0;
    MLMatrixCL::const_iterator coarse = A.begin(), fine = ++A.begin();
    MLMatrixCL::const_iterator coarseP = P.begin(), fineP = ++P.begin();
    for(; fine!=A.end(); ++coarse, ++fine, ++lvl)
    {
//        if (lvl==1) std::cerr << fine->A.Data << std::endl;
        const Uint nu= coarse->num_cols();
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
            Ai= transp_mul( *fineP, *fine * (*fineP * ei) );
            for(Uint j=0; j<nu; ++j)
            {
                if (std::fabs(Ai[j] - (*coarse)(i,j))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": A_H(" << i << ", " << j <<")= "
                              << (*coarse)(i,j) << " != " << Ai[j] << std::endl;
                }
            }
        }
    }
}

} // end of namespace DROPS
