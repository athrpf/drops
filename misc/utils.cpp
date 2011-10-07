/// \file utils.cpp
/// \brief Useful stuff that fits nowhere else.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Martin Horsky, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/utils.h"
#include <iostream>

#ifdef DROPS_WIN
#include <direct.h>
double cbrt(double arg)
{
    return pow(arg, 1.0/3);
}
#endif

namespace DROPS
{

#ifndef _PAR
// Definition of parallel version, see parallel/parallel.cpp
std::ostream&
DROPSErrCL::what(std::ostream& out) const
{
    out << _ErrMesg << std::endl;
    return out;
}

void
DROPSErrCL::handle() const
{
    what(std::cerr );
    std::cerr.flush();
    std::abort();
}
#endif

PermutationT
invert_permutation (const PermutationT& p)
{
    PermutationT pi( p.size());
    for (size_t i= 0; i < p.size(); ++i)
        pi[p[i]]= i;

    return pi;
}

PermutationT
compose_permutations (const PermutationT& p, const PermutationT& q)
{
    if (q.empty())
        return p;
    if (p.empty())
        return q;

    PermutationT ret( p.size());
    for (size_t i= 0; i < p.size(); ++i)
        ret[i]= p[q[i]];

    return ret;
}


int CreateDirectory(std::string path)
/** Used to create directories
    \param path path of new directory
    \return returns error code of mkdir
*/
{
#ifdef DROPS_WIN
    throw DROPSErrCL("Creating directories is not implemented for Windows OS");
    return _mkdir( path.c_str());
#else
    return mkdir(path.c_str(), 0777);
#endif
}

int DeleteFile(std::string file)
/** Used to delete files
    \param file name of the file
    \return returns error code of remove
*/
{
#ifdef DROPS_WIN
    throw DROPSErrCL("Deleting Files is not implemented for Windows OS");
    return 0;
#else
    return remove(file.c_str());
#endif
}

void reverseByteOrder(int size,char field[])
{
    std::vector<char> temp(size);
    for(int i=0;i<size;++i)
        temp[i]= field[i];

    for(int i=0;i<size;++i)
        field[(size-1)-i] = temp[i];
}


} // end of namespace DROPS
