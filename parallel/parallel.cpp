/// \file parallel.cpp
/// \brief interface to MPI and helper functions to easily use communication
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "parallel/parallel.h"
#include "parallel/pardistributeddata.h"
#include <limits>
#include "misc/utils.h"

namespace DROPS
{
/***************************************************************************
*   P R O C - C L A S S                                                    *
***************************************************************************/
Uint    ProcCL::my_rank_=0;
Uint    ProcCL::size_   =0;             // if _size==0, then this proc has not created a ProcCL
int     ProcCL::procDigits_=0;
ProcCL* ProcCL::instance_=0;            // only one instance of ProcCL may exist (Singleton-Pattern)
MuteStdOstreamCL* ProcCL::mute_=0;

#ifdef _MPICXX_INTERFACE
    const ProcCL::CommunicatorT& ProcCL::Communicator_ = MPI::COMM_WORLD;
    const ProcCL::DatatypeT      ProcCL::NullDataType  = MPI::Datatype(MPI_DATATYPE_NULL);

    const ProcCL::DatatypeT& ProcCL::MPI_TT<int>::dtype    = MPI::INT;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<Uint>::dtype   = MPI::UNSIGNED;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<Ulint>::dtype  = MPI::UNSIGNED_LONG;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<Usint>::dtype  = MPI::UNSIGNED_SHORT;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<double>::dtype = MPI::DOUBLE;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<char>::dtype   = MPI::CHAR;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<byte>::dtype   = MPI::CHAR;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<float>::dtype  = MPI::FLOAT;
#ifdef DROPS_WIN
    const ProcCL::DatatypeT& ProcCL::MPI_TT<size_t>::dtype = MPI::UNSIGNED;
#endif
#else
    const ProcCL::CommunicatorT& ProcCL::Communicator_ = MPI_COMM_WORLD;
    const ProcCL::DatatypeT      ProcCL::NullDataType  = MPI_DATATYPE_NULL;

    const ProcCL::DatatypeT& ProcCL::MPI_TT<int>::dtype    = MPI_INT;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<Uint>::dtype   = MPI_UNSIGNED;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<Ulint>::dtype  = MPI_UNSIGNED_LONG;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<Usint>::dtype  = MPI_UNSIGNED_SHORT;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<double>::dtype = MPI_DOUBLE;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<char>::dtype   = MPI_CHAR;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<byte>::dtype   = MPI_CHAR;
    const ProcCL::DatatypeT& ProcCL::MPI_TT<float>::dtype  = MPI_FLOAT;
#ifdef DROPS_WIN
# ifdef WIN64
    const ProcCL::DatatypeT& ProcCL::MPI_TT<size_t>::dtype = MPI_UNSIGNED_LONG;
# endif
#endif
#endif

ProcCL::ProcCL(int* argc, char*** argv)
{
    Assert(size_==0, DROPSErrCL("ProcCL instanciated multiple times"), DebugParallelC);

    DynamicDataInterfaceCL::Init(argc, argv);               // DDD Initialisieren und die Informationen beziehen
    my_rank_ = DynamicDataInterfaceCL::InfoMe();
    size_    = DynamicDataInterfaceCL::InfoProcs();
    procDigits_= 1;
    int procs  = Size();
    while( procs>9){
        ++procDigits_; procs/=10;
    }
    mute_    = new MuteStdOstreamCL();
    MuteStdOstreams();
}

ProcCL::~ProcCL()
{
	DynamicDataInterfaceCL::Exit();             // Logoff from DDD
    size_=0;                // Now, this class can be initialized again...
    RecoverStdOstreams();
    delete mute_;
}

void ProcCL::Prompt(int me)
{
    char in;
    if (MyRank()==me)
        std::cin >> in;
    Barrier();
}

void ProcCL::AppendProcNum( std::string& str)
{
    char format[]= ".%0Xi", postfix[8];
    format[3]= '0' + char(procDigits_);
    std::sprintf( postfix, format, ProcCL::MyRank());
    str+= postfix;
}

#ifdef _PAR
std::ostream&
DROPSErrCL::what(std::ostream& out) const
{
    out << "["<<ProcCL::MyRank()<<"] ";
    out << _ErrMesg << std::endl;
    return out;
}

void
DROPSErrCL::handle() const
{
    ProcCL::RecoverStdOstreams();
    what(std::cout);
    std::cout.flush();
    ProcCL::Abort(-1);
    std::abort();
}
#endif

} // namespace DROPS


