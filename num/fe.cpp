//**************************************************************************
// File:    fe.cpp                                                         *
// Content: description of various finite-element functions                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - May, 2 2001                                            *
//**************************************************************************

#ifndef _FE_CPP_
#define _FE_CPP_

#include "num/fe.h"

namespace DROPS
{

const double FE_P1CL::_gradient[4][3]=
    { {-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} };

} // end of namespace DROPS

#endif
