//**************************************************************************
// File:    unknowns.cpp                                                   *
// Content:                                                                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt IGPM RWTH Aachen     *
// Version: 0.1                                                            *
// History: begin - March, 14 2001                                         *
//                                                                         *
//**************************************************************************

#ifndef _UNKNOWNS_CPP_
#define _UNKNOWNS_CPP_

#include "num/unknowns.h"
#include <algorithm>


namespace DROPS
{


UnknownIdxCL::UnknownIdxCL(const UnknownIdxCL& orig)
    : _Idx( orig._Idx.size() ), _Len(orig._Len)
{
    for (Uint i=0; i<_Idx.size(); ++i)
    {
        _Idx[i]= new IdxT[_Len[i]];
        std::copy(orig._Idx[i], orig._Idx[i]+_Len[i], _Idx[i]);
    }
}


UnknownIdxCL& UnknownIdxCL::operator=(const UnknownIdxCL& rhs)
{
    if(&rhs == this) return *this;

    for (Uint i=0; i<_Idx.size(); ++i) Dealloc(i);
    _Len= rhs._Len;
    _Idx.resize( rhs._Idx.size() );
    for (Uint i=0; i<_Idx.size(); ++i)
    {
        _Idx[i]= new IdxT[_Len[i]];
        std::copy(rhs._Idx[i], rhs._Idx[i]+_Len[i], _Idx[i]);
    }
    return *this;
}


} // end of namespace DROPS


#endif
