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

#include<map>
#include<string>
#include "num/discretize.h"
namespace DROPS
{
template<class T>
class SingletonBaseCL : public std::map<std::string, T>
{
  private:
     SingletonBaseCL() {}                         // von außen keine Instanzen erzeugbar
     SingletonBaseCL(const SingletonBaseCL&) {}  // nicht kopierbar
     ~SingletonBaseCL() {}
  public:
     static SingletonBaseCL& getInstance();
};

typedef SingletonBaseCL<DROPS::instat_vector_fun_ptr> InVecMap;
typedef SingletonBaseCL<DROPS::instat_scalar_fun_ptr> InScaMap;

Point3DCL TestFunction(const Point3DCL& , double);

class RegisterVelFunction
{
	public:
	RegisterVelFunction(std::string, instat_vector_fun_ptr);
};

class RegisterScalarFunction
{
	public:
	RegisterScalarFunction(std::string, instat_scalar_fun_ptr);
};

} //end of namespace DROPS
