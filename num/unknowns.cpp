//**************************************************************************
// File:    unknowns.cpp                                                   *
// Content:                                                                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt IGPM RWTH Aachen     *
// Version: 0.1                                                            *
// History: begin - March, 14 2001                                         *
//                                                                         *
//**************************************************************************

#include "num/unknowns.h"
#include <algorithm>


namespace DROPS
{


UnknownIdxCL::UnknownIdxCL( const UnknownIdxCL& orig)
    : _Idx( orig._Idx) {}


UnknownIdxCL& UnknownIdxCL::operator=( const UnknownIdxCL& rhs)
{
    if(&rhs == this) return *this;
    _Idx= rhs._Idx;
    return *this;
}


} // end of namespace DROPS
