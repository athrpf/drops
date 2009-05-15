/***************************************************************************
*  File:    parallel.cpp                                                   *
*  Content: Interface for parallel support                                 *
*           ProcCL - Management of the procs                               *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           10.11.2005                                             *
*  last modified:   10.11.2005                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file parallel.cpp

#include "parallel/parallel.h"
#include <ddd.h>
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
    const ProcCL::DatatypeT& ProcCL::MPI_TT<size_t>::dtype = MPI_UNSIGNED;
#endif
#endif

ProcCL::ProcCL(int* argc, char*** argv)
{
    Assert(size_==0, DROPSErrCL("ProcCL instanciated multiple times"), DebugParallelC);

    DDD_Init(argc, argv);               // DDD Initialisieren und die Informationen beziehen
    my_rank_ = DDD_InfoMe();
    size_    = DDD_InfoProcs();
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
    DDD_Exit();             // Logoff from DDD
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


