/// \file   kd_tree_utils.h
/// \brief  Utils for implementing the kd-tree
/// \author Oliver Fortmeier, fortmeier@sc.rwth-aachen.de

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


#ifndef KD_TREE_UTILS_H
#define KD_TREE_UTILS_H

// include all stl headers used in this project
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <queue>
#include <string>
#include <valarray>
#include <vector>


/// \name OpenMP helper functions
//@{
#ifdef _OPENMP
   #include <omp.h>
    /// \brief get the current wall time
    inline double get_time() { return omp_get_wtime(); }
    /// \brief get the number of available threads
    inline int get_num_threads()
    {
        return omp_get_max_threads();
    }
    /// \brief get the id of the calling thread
    inline int get_thread_id()
    {
        return omp_get_thread_num();
    }
#ifdef _MSC_VER
#include <windows.h>
    /// \brief bind a thread to a core
    inline void set_thread_affinity () 
    {
      #pragma omp parallel default(shared)
      {
        DWORD_PTR mask = (1 << omp_get_thread_num());
        SetThreadAffinityMask( GetCurrentThread(), mask );
      }
    }
#undef max
#undef min
#endif
#else
   inline double get_time() { return 0; }
   inline int get_num_threads() { return 1; }
   inline int get_thread_id() { return 0; }
   inline void set_thread_affinity () {}
#endif
//@}

namespace DROPS{
namespace KDTree{
    typedef unsigned short int usint;       /// type of small ints like dimension
}
}

#endif
