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
    SingletonMapCL<T>& SingletonMapCL<T>::getInstance()
    {
       static SingletonMapCL instance;
       return instance;
    }

    template<class T>
    T& SingletonMapCL<T>::operator[](std::string s){
        if (this->find(s) == this->end()){
            std::ostringstream os;
            os << "function with the name \"" << s << "\" not found in container!";
            throw DROPSErrCL(os.str());
        }
        return this->find(s)->second;
    }

    RegisterVectorFunction::RegisterVectorFunction(std::string name, instat_vector_fun_ptr fptr){
        InVecMap::getInstance().insert(std::make_pair(name,fptr));
    }

    RegisterScalarFunction::RegisterScalarFunction(std::string name, instat_scalar_fun_ptr fptr){
        InScaMap::getInstance().insert(std::make_pair(name,fptr));
    }

    RegisterScalarFunction::RegisterScalarFunction(std::string name, scalar_fun_ptr fptr){
        ScaMap::getInstance().insert(std::make_pair(name,fptr));
    }
    
    RegisterMatchingFunction::RegisterMatchingFunction(std::string name, match_fun fptr){
        MatchMap::getInstance().insert(std::make_pair(name,fptr));
    }
    
    RegisterMatrixFunction::RegisterMatrixFunction(std::string name, instat_matrix_fun_ptr fptr){
        InMatMap::getInstance().insert(std::make_pair(name,fptr));
    }



    template class SingletonMapCL<DROPS::instat_scalar_fun_ptr>;
    template class SingletonMapCL<DROPS::instat_vector_fun_ptr>;
    template class SingletonMapCL<DROPS::scalar_fun_ptr>;
    template class SingletonMapCL<DROPS::instat_matrix_fun_ptr>;
    template class SingletonMapCL<DROPS::match_fun>;

} //end of namespace DROPS
