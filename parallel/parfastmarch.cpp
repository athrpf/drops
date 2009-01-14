// **************************************************************************
// File:    parfastmarch.h                                                  *
// Content: Classes for performing a parallel fastmarching algorithm        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen     *
//          Oliver Fortmeier, RZ RWTH Aachen                                *
// Version: 0.1                                                             *
// Date:                                                                    *
// Begin:   August  21th, 2006                                              *
// **************************************************************************
/// \author Oliver Fortmeier
/// \file parfastmarch.cpp

#include "parallel/parfastmarch.h"

namespace DROPS
{

//------------------------
// R E P R  T E T R A  C L
//------------------------

bool ReprTetraCL::HasGhost() const
{
    for (Uint i=0; i<4; ++i)
        if (IsGhost(i))
            return true;
    return false;
}

bool operator==(const ReprTetraCL& tetra1, const ReprTetraCL& tetra2)
{
    std::vector<IdxT> DoF1, DoF2;
    for (Uint i=0; i<4; ++i){
        DoF1.push_back(tetra1.data_[i].second);
        DoF2.push_back(tetra2.data_[i].second);
    }
    std::sort(DoF1.begin(), DoF1.end());
    std::sort(DoF2.begin(), DoF2.end());

    for (std::vector<IdxT>::iterator d1(DoF1.begin()), d2(DoF2.begin()), end(DoF1.end()); d1!=end; ++d1,++d2)
        if (*d1!=*d2)
            return false;
    return true;
}

//--------------------------
// F M  T R A N S F E R  C L
//--------------------------


} // end of namespace DROPS
