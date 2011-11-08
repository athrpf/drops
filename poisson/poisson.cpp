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

instat_scalar_fun_ptr StripTimeCL::_func= NULL;
double                StripTimeCL::_t= 0;

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

double Py_product(MultiGridCL& mg, IdxDescCL& Idx, MatDescCL& A, MatDescCL& M,
instat_scalar_fun_ptr f1, instat_scalar_fun_ptr f2, double t, bool H1)
{
    VecDescCL vf1, vf2;
    vf1.SetIdx(&Idx);
    vf2.SetIdx(&Idx);

    double ret;
    const Uint lvl = vf1.GetLevel(),
               idx = vf1.RowIdx->GetIdx();
    IdxT UnknownIdx[4];
    for (MultiGridCL::const_TriangTetraIteratorCL
            sit=const_cast<const MultiGridCL&>(mg).GetTriangTetraBegin(lvl),
            send=const_cast<const MultiGridCL&>(mg).GetTriangTetraEnd(lvl);
            sit != send; ++sit)
    {
        for(int i=0; i<4; ++i)
        {
            UnknownIdx[i]= sit->GetVertex(i)->Unknowns.Exist(idx) ? sit->GetVertex(i)->Unknowns(idx)
                    : NoIdx;
        }
        for(int i=0; i<4; ++i)
        {
            if (UnknownIdx[i] != NoIdx)  //not on a dirichlet boundary condition
            {
                vf1.Data[UnknownIdx[i]] = f1(sit->GetVertex(i)->GetCoord(), t);
                vf2.Data[UnknownIdx[i]] = f2(sit->GetVertex(i)->GetCoord(), t);
            }
        }
    }
    if(H1)
    {
        ret = dot(M.Data*vf1.Data, vf2.Data) + dot(A.Data* vf1.Data, vf2.Data);
    }
    else
        ret = dot(M.Data*vf1.Data, vf2.Data);

    return ret;
}


} // end of namespace DROPS
