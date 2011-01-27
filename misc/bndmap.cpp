/// \file bndmap.cpp
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
#include "misc/bndmap.h"

namespace DROPS
{

template<class T>
  SingletonBaseCL<T>& SingletonBaseCL<T>::getInstance()
  {
     static SingletonBaseCL instance;
     return instance;
  }

  RegisterVelFunction::RegisterVelFunction(std::string name, instat_vector_fun_ptr fptr){
    InVecMap::getInstance().insert(std::make_pair(name,fptr));
  }

  RegisterScalarFunction::RegisterScalarFunction(std::string name, instat_scalar_fun_ptr fptr){
	InScaMap::getInstance().insert(std::make_pair(name,fptr));
  }

} //end of namespace DROPS
