/// \file MGsolver.cpp
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Maxim Larin, Volker Reichelt; SC RWTH Aachen:

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "num/MGsolver.h"

namespace DROPS
{

void CheckMGData( const MLMatrixCL& A, const MLMatrixCL& P)
{
    Uint lvl= 0;
    MLMatrixCL::const_iterator coarse = A.begin(), fine = ++A.begin();
    MLMatrixCL::const_iterator fineP = ++P.begin();
    for(; fine!=A.end(); ++coarse, ++fine, ++fineP, ++lvl)
    {
//        if (lvl==1) std::cout << fine->A.Data << std::endl;
        const Uint nu= coarse->num_cols();
        VectorCL ei(nu), Ai(nu);

        std::cout << nu << " unknowns on level " << lvl << std::endl;

        if (nu>700)
        {
            std::cout << "Check skipped: too many unknowns, too much time..." << std::endl;
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
                    std::cout << "Auf coarse-Level " << lvl << ": A_H(" << i << ", " << j <<")= "
                              << (*coarse)(i,j) << " != " << Ai[j] << std::endl;
                }
            }
        }
    }
}

} // end of namespace DROPS
