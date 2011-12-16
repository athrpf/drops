/// \file   bucket.h
/// \brief  Bucket of a node of a kd-tree
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


#ifndef BUCKET_H
#define BUCKET_H

#include "misc/kd-tree/kd_tree_utils.h"

namespace DROPS{
namespace KDTree{

    /// \brief no index
    const size_t NoIdx= size_t(-1);

    namespace internal{

        /// \brief A bucket that containing a BucketSize indices
        /** The indices of the data is stored. This is not only specified 
            by [first,last) in order to make it more general (i.e., inserting points 
            should be possible)
            \tparam BucketSize   size of the buckets
        */
        template <int BucketSize>
        class BucketCL
        {
        private:  // ------- member variables ------- 
            size_t p_vals[BucketSize];                                      ///< the indices in this bucket

        public:  // ------- member functions ------- 
            BucketCL() {}                                                   ///< constructor
            ~BucketCL() {}                                                  ///< destructor
            const size_t& operator[] (int i) const { return p_vals[i]; }    ///< access the elements in the bucket
            size_t& operator[] (int i) { return p_vals[i]; }                ///< set a value in the bucket
            static size_t memory() { return BucketSize*sizeof(size_t); }    ///< memory for storing the bucket
        };
    }
}
}
#endif
