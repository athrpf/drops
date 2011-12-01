/// \file   doc.h
/// \brief  Documentation of kd-trees
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


#ifndef DOC_H
#define DOC_H

/** \mainpage Oliver Fortmeier's implementation of a kd-tree that is parallelized by OpenMP
   
    This kd-tree represents points in R^K where K is the dimension of the points.

    The main classes are the following ones.
    <ul>
      <li> 
        The tree is given by the class \ref TreeCL. This class stores the kd-tree.
      </li>
      <li>
        The tree is built by the builder pattern by the class \ref TreeBuilderCL.
      </li>
      <li>
        Various search operations like search m nearest neighbors, neighbors in a given
        sphere, etc., are implemented by the strategy pattern. The search operations
        are given by derived classes of \ref BaseSearchCL.
      </li>
    </ul>
    
    In this project, some metrics are predefined, see metric.h. These metrics include
    the following metrics which are choosen by a template parameter \a metric.
    <ul>
      <li>
        Euclidian metric, \a metric=2, \f$ \| a-b \|_2 = \sqrt{\sum^K_{i=1} (a_i-b_i)^2} \f$
      </li>
      <li>
        one metric, \a metric=1, \f$ \| a-b \|_1 = \sum^K_{i=1} |a_i-b_i| \f$ (not implemented so far)
      </li>
      <li>
        sup metric, \a metric=0, \f$ \| a-b \|_\infty = \max_{i=1,\dots,K} |a_i-b_i| \f$ 
      </li>
      <li>
        Markowski metric, \a metric=p, \f$ \| a-b \|_p = \sqrt[p]{\sum^K_{i=1} |a_i-b_i|^p} \f$ (not implemented so far)
      </li>
    </ul>

    However, the user is free to implement own metrics. Therefore, create a class, say MyMetric, the has two functions
    distance and intersects. For an example, we refer to \ref EuclideanMetricCL.

    For a tutorial, we refer to the functions \ref small_tutorial, \ref large_tutorial, and \ref equal_points_tutorial.

    <p> This software has been tested on various hardware architectures using the following compilers.
    <ul> 
      <li> Visual Studio 10.0 </li>
      <li> intel version 11.1 - 12.0 </li>
      <li> gcc version 4.2 - 4.6 </li>
      <li> pgi 11.1 </li>
      <li> oracle studio 12.2 </li>
    </ul>
*/

namespace DROPS{
/// \brief Namespace including all functionality of the kd-tree.
namespace KDTree{
    /// \brief Internal helpers for implementing the kd-tree
    /** Internal helper classes and functions for the kd-tree.*/
    namespace internal{}
}
}
#endif
