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


PermutationT
invert_permutation (const PermutationT& p)
{
    PermutationT pi( p.size());
    for (size_t i= 0; i < p.size(); ++i)
        pi[p[i]]= i;

    return pi;
}

} // end of namespace DROPS
