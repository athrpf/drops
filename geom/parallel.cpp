//**************************************************************************
// File:     parallel.cpp                                                  *
// Content:  Management of multiple processors via MPI                     *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.1                                                           *
// History:  begin - August, 2 2000                                        *
//                                                                         *
// Remarks:                                                                *
//**************************************************************************

#include "geom/parallel.h"


namespace DROPS
{

MultiProcCL MultiProc;


//**************************************************************************
// Class:    MultiProcCL                                                   *
// Purpose:  Does MPI Initialization                                       *
//           and holds data for parallel grid management                   *
// Remarks:  Init() needs the commandline of the programme for MPI         *
//           initialization. It must be called prior to any MPI-routine    *
//           Finalize() must be called once before programme termination.  *
//**************************************************************************

void MultiProcCL::Init ( int &/*pargc*/, char **&/*pargv*/ )
{
    if (_IsValid) throw DROPSErrCL("MultiProcCL::Init: Already initialized!");
    else _IsValid=true;
    try
    {
// TODO:
//    MPI::Init(pargc,pargv);
//    NumberProcs=MPI::COMM_WORLD.Get_size();
//    MyProc=MPI::COMM_WORLD.Get_rank();
//
//    VertexCopies=new ListCL<VertexCL*[_NumberProcs];
// TODO: etc.
        _NumberProcs = 1;
        _MyProc = 0;
    }
// TODO: MPI-Fehlerklasse bzw. Fehler von new abfangen
    catch (DROPSErrCL& err) { err.handle(); }
}


void MultiProcCL::Finalize ()
{
    if (!_IsValid) throw DROPSErrCL("MultiProcCL::Finalize: Already finalized!");
// TODO:
//    delete [] VertexCopies; etc.
//    MPI::Finalize(...);
    _IsValid=false;
}

} // end of namespace DROPS
