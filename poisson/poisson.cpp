/// \file poisson.cpp
/// \brief classes that constitute the poisson-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "num/discretize.h"
#include "poisson/poisson.h"

namespace DROPS
{

//===================================================
//
//                   Error estimators
//
//===================================================

double SimpleGradEstimator (const TetraCL& t, const VecDescCL& lsg, const PoissonBndDataCL& Bnd)
{
    IdxT UnknownIdx[4];
    const Uint Idx= lsg.RowIdx->GetIdx();
    double val1, val2, diff, maxdiff=0;

    for(int i=0; i<4; ++i)
        UnknownIdx[i]= t.GetVertex(i)->Unknowns.Exist() ? t.GetVertex(i)->Unknowns(Idx) : -1ul;
    for(int i=1; i<4; ++i)
        for(int j=0; j<i;++j)
        {
            val1= (t.GetVertex(i)->Unknowns.Exist()) ? lsg.Data[UnknownIdx[i]]
                 : Bnd.GetDirBndValue(*t.GetVertex(i));
            val2= (t.GetVertex(j)->Unknowns.Exist()) ? lsg.Data[UnknownIdx[j]]
                 : Bnd.GetDirBndValue(*t.GetVertex(j));
            diff= std::fabs(val1-val2);
            if (diff>maxdiff)
                maxdiff= diff;
        }
    return maxdiff;
}


} // end of namespace DROPS
