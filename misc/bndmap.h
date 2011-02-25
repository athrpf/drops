/// \file bndmap.h
/// \brief global map for boundary functions
/// \author LNM RWTH Aachen: Martin Horsky; SC RWTH Aachen:
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

#ifndef BNDMAP_H
#define BNDMAP_H

#include <map>
#include <string>
#include "num/discretize.h"
#include <sstream>

namespace DROPS
{
template<class T>
class SingletonMapCL : public std::map<std::string, T>
{
  private:
    SingletonMapCL() {}                         // von aussen keine Instanzen erzeugbar
    SingletonMapCL(const SingletonMapCL&) : std::map<std::string, T>() { }  // nicht kopierbar
    ~SingletonMapCL() {}
  public:
    static SingletonMapCL& getInstance();
    T operator[](std::string s);
};

typedef SingletonMapCL<DROPS::instat_vector_fun_ptr> InVecMap;
typedef SingletonMapCL<DROPS::instat_scalar_fun_ptr> InScaMap;
typedef SingletonMapCL<DROPS::scalar_fun_ptr> ScaMap;
typedef SingletonMapCL<DROPS::match_fun> MatchMap;

Point3DCL TestFunction(const Point3DCL& , double);

class RegisterVectorFunction
{
  public:
    RegisterVectorFunction(std::string, instat_vector_fun_ptr);
};

class RegisterScalarFunction
{
  public:
    RegisterScalarFunction(std::string, instat_scalar_fun_ptr);
    RegisterScalarFunction(std::string, scalar_fun_ptr);
};

class RegisterMatchingFunction
{
  public:
    RegisterMatchingFunction(std::string, match_fun);
};

} //end of namespace DROPS

#endif
