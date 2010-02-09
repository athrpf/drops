/// \file addeddata.h
/// \brief manage the transfer of numerical data within a load balancing step
/// \author LNM RWTH Aachen: Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

/// \todo (of) Addeddata classes are inefficient, because to every transfered numerical
///   data the type of the numerical data is submitted too. Maybe it would make more
///   sence, to use XferAddDataX!

#include "parallel/distributeddatatypes.h"
#include "misc/utils.h"
#include "misc/container.h"

#ifndef DROPS_ADDEDDATA_H
#define DROPS_ADDEDDATA_H

namespace DROPS
{
/****************************************************************************
* A D D E D  D A T A  C L A S S                                             *
****************************************************************************/
/// \brief Manages the transfer of numerical data
/** This class stores unknowns of SkalarType */
class AddedScalCL
{
  private:
    static TypeT _dddT;
    static void Declare();          // implemented in "parallel/parmultigrid.cpp"
    static void Define();

  public:
    Uint   idxVecDesc_;             ///< index, in which VecDescCL the to be transfered unknowns stands in the ParMultiGridCL of the sender
    double data_;                   ///< the skalar unknown

    AddedScalCL();                  // std-constructor
    AddedScalCL(Uint idx, double val) : idxVecDesc_(idx), data_(val) {}

    Uint     GetIdx()  const {return idxVecDesc_;}
    double   GetData() const {return data_;}
    static TypeT GetType()       {return _dddT;}

    void SetIdx(Uint i) {idxVecDesc_=i;}
    void SetData(double d) {data_=d;}

    friend class ParMultiGridCL;
};

class AddedVecCL
{
  private:
    static TypeT _dddT;
    static void Declare();      // implemented in "parallel/parmultigrid.cpp"
    static void Define();

  public:
    Uint      idxVecDesc_;      ///< index, in which VecDescCL the to be transfered unknowns stands in the ParMultiGridCL of the sender
    Point3DCL data_;            ///< the vectoriel unknown

    AddedVecCL() {}
    AddedVecCL(Uint idx, Point3DCL val) : idxVecDesc_(idx), data_(val) {}

    Uint       GetIdx()  const {return idxVecDesc_;}
    Point3DCL& GetData()       {return data_;}
    static TypeT   GetType()       {return _dddT;}

    void SetIdx(Uint i)              {idxVecDesc_=i;}
    void SetData(const Point3DCL& d) {data_=d;}


    friend class ParMultiGridCL;
};

} // end of namespace DROPS

#endif
