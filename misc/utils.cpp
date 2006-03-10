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
    std::abort();
}

} // end of namespace DROPS
