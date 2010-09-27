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

int Test1( int seed)
{ // check that a P1 function is correctly extended
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
        if (norm(diff) > 1e-18) {
            std::cout << "test 1: maximum difference = " << diff.max() << " when extending on child " << ch << "!\n";
            ++err;
        }
    }
    return err;
}

int Test2( int seed)
{ // check that the extension is indeed a P1 function
//    std::cout << "running test " << seed << "...\n";
    int err= 0;
    RefRuleCL data= GetRefRule(RegRefRuleC);
    LocalP2CL<> p2func;
    srand( seed);
    for (int i=0; i<10; ++i)
        p2func[i]= rand();

    LocalP2CL<> extension, diff;

    for (int ch=0; ch<8; ++ch) {
        ExtendP1onChild( p2func, ch, extension);
        // check that values are the same on child's vertices
        const ChildDataCL child= GetChildData( data.Children[ch]);
        for (int i=0; i<4; ++i) {
            const int v= child.Vertices[i];
            if (std::abs(p2func[v] - extension[v]) > 1e-18) {
                std::cout << "test 2: values are not equal for child " << ch << " and dof " << v << ", rel. difference = " << std::abs(p2func[v] - extension[v])/p2func[v] << "!\n";
                ++err;
            }
        }
        // check linearity of extension
        LocalP1CL<> p1func;
        for (int i=0; i<4; ++i)
            p1func[i]= extension[i];
        LocalP2CL<> P1asP2( p1func);
        diff= P1asP2 - extension; // should be zero
        if (norm(diff) > 1e-18) {
            std::cout << "test 2: maximum difference = " << diff.max() << " when extending on child " << ch << "!\n";
            ++err;
        }
    }
    return err;
}

int main()
{
  try {
    int err1= 0, err2= 0;
    for (int i=0; i<10000; ++i) {
        err1+= Test1(i);
        err2+= Test2(i);
    }
    return err1+err2;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
