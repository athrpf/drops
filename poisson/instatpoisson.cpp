//**************************************************************************
// File:    instatpoisson.cpp                                              *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - November, 11 2002                                      *
//**************************************************************************

#ifndef _INSTAT_POISSON_CPP_
#define _INSTAT_POISSON_CPP_

#include "poisson/instatpoisson.h"
#include "num/discretize.h"

namespace DROPS
{

  scalar_instat_fun_ptr StripTimeCL::_func= NULL;
  double                StripTimeCL::_t= 0;

} // end of namespace DROPS

#endif



