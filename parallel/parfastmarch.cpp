/// \file parfastmarch.cpp
/// \brief parallel version of fast marching method (line by line parallelization)
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier

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

/// This method is, from a performance point of view, not good. So, this
/// file may become deprecated.

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
