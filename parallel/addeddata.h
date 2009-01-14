//**************************************************************************
// File:    addeddata.h                                                    *
// Content: Class that manages the transfer of numerical data              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, RZ RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   Januar, 03rd, 2006                                             *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file addeddata.h
/// \brief Classes to manage the transfer of numerical data within a load balancing step
/// \todo (of) Addeddata classes are inefficient, because to every transfered numerical data the type of the numerical data is submitted too.
///  Maybe it would make more sence, to use XferAddDataX!

#include <ddd.h>
#include "misc/utils.h"
#include "misc/container.h"

#ifndef _ADDEDDATA_H_
#define _ADDEDDATA_H_

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
    static DDD_TYPE _dddT;
    static void Declare();          // implemented in "parallel/parmultigrid.cpp"
    static void Define();

  public:
    Uint   idxVecDesc_;             ///< index, in which VecDescCL the to be transfered unknowns stands in the ParMultiGridCL of the sender
    double data_;                   ///< the skalar unknown

    AddedScalCL();                  // std-constructor
    AddedScalCL(Uint idx, double val) : idxVecDesc_(idx), data_(val) {}

    Uint     GetIdx()  const {return idxVecDesc_;}
    double   GetData() const {return data_;}
    static DDD_TYPE GetType()       {return _dddT;}

    void SetIdx(Uint i) {idxVecDesc_=i;}
    void SetData(double d) {data_=d;}

    friend class ParMultiGridCL;
};

class AddedVecCL
{
  private:
    static DDD_TYPE _dddT;
    static void Declare();      // implemented in "parallel/parmultigrid.cpp"
    static void Define();

  public:
    Uint      idxVecDesc_;      ///< index, in which VecDescCL the to be transfered unknowns stands in the ParMultiGridCL of the sender
    Point3DCL data_;            ///< the vectoriel unknown

    AddedVecCL() {}
    AddedVecCL(Uint idx, Point3DCL val) : idxVecDesc_(idx), data_(val) {}

    Uint       GetIdx()  const {return idxVecDesc_;}
    Point3DCL& GetData()       {return data_;}
    static DDD_TYPE   GetType()       {return _dddT;}

    void SetIdx(Uint i)              {idxVecDesc_=i;}
    void SetData(const Point3DCL& d) {data_=d;}


    friend class ParMultiGridCL;
};

} // end of namespace DROPS

#endif
