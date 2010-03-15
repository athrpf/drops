/// \file extendP1onChild.cpp
/// \brief tests implementation of ExtendP1onChild
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

using namespace DROPS;

int Test( int seed)
{
//    std::cout << "running test " << seed << "...\n";
    int err= 0;

    LocalP1CL<> p1func;
    srand( seed);
    for (int i=0; i<4; ++i)
        p1func[i]= rand();

    LocalP2CL<> P1asP2( p1func), extension, diff;

    for (int ch=0; ch<8; ++ch) {
        ExtendP1onChild( P1asP2, ch, extension);
        diff= P1asP2 - extension; // should be zero
        if (diff.max() > 1e-18) {
            std::cout << "maximum difference = " << diff.max() << " when extending on child " << ch << "!\n";
            ++err;
        }
    }
    return err;
}

int main()
{
  try {
    int err= 0;
    for (int i=0; i<10000; ++i)
        err+= Test(i);
    return err;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
