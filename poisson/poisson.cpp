//**************************************************************************
// File:    poisson.cpp                                                    *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

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
