/// \file unknowns.cpp
/// \brief Implementation of the mapping from simplices to indices to
///    linear-algebra data-structures.)
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#include "num/unknowns.h"
#include <algorithm>


namespace DROPS
{


UnknownIdxCL::UnknownIdxCL( const UnknownIdxCL& orig)
    : _Idx( orig._Idx) {}


UnknownIdxCL& UnknownIdxCL::operator=( const UnknownIdxCL& rhs)
{
    if(&rhs == this) return *this;
    _Idx= rhs._Idx;
    return *this;
}


} // end of namespace DROPS
