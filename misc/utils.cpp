/// \file
/// \brief Useful stuff that fits nowhere else.


#include "misc/utils.h"
#include <iostream>

namespace DROPS
{

std::ostream&
DROPSErrCL::what(std::ostream& out) const
{
    out << _ErrMesg << std::endl;
    return out;
}

void
DROPSErrCL::handle() const
{
    what(std::cerr);
    std::cerr.flush();
    abort();
}

double TimerCL::_gtime= 0;
clock_t TimerCL::_gt_begin= 0, 
	TimerCL::_gt_end=0;


} // end of namespace DROPS
