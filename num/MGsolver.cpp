//**************************************************************************
// File:    MGpoisson.cpp                                                  *
// Content: classes that constitute the poisson-problem with MG-solver     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - April, 16 2001                                         *
//**************************************************************************

#ifndef _MGSOLVER_CPP_
#define _MGSOLVER_CPP_

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


} // end of namespace DROPS

#endif

