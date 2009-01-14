/// \file
/// \brief Useful stuff that fits nowhere else.


#include "misc/utils.h"
#include <iostream>

namespace DROPS
{

#ifndef _PAR
// Definition of parallel version, see parallel/parallel.cpp
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
#endif

PermutationT
invert_permutation (const PermutationT& p)
{
    PermutationT pi( p.size());
    for (size_t i= 0; i < p.size(); ++i)
        pi[p[i]]= i;

    return pi;
}

int CreateDirectory(std::string path)
/** Used to create directories
    \param path path of new directory
    \return returns error code of mkdir
*/
{
    return mkdir(path.c_str(), 0777);
}

int DeleteFile(std::string file)
/** Used to delete files
    \param file name of the file
    \return returns error code of remove
*/
{
    return remove(file.c_str());
}

} // end of namespace DROPS
