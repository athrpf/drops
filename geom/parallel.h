//**************************************************************************
// File:     parallel.h                                                    *
// Content:  Management of multiple processors via MPI                     *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.1                                                           *
// History:  begin - August, 2 2000                                        *
//                                                                         *
// Remarks:                                                                *
//**************************************************************************

#ifndef _PARALLEL_H_
#define _PARALLEL_H_


#include "misc/utils.h"


namespace DROPS
{


//**************************************************************************
// Class:    MultiProcCL                                                   *
// Purpose:  Does MPI Initialization                                       *
//           and holds data for parallel grid management                   *
// Remarks:  Init() needs the commandline of the programme for MPI         *
//           initialization. It must be called prior to any MPI-routine    *
//           Finalize() must be called once before programme termination.  *
//**************************************************************************

class MultiProcCL
{
private:
    Uint _NumberProcs;
    Uint _MyProc;
    bool _IsValid;

//    ListCL<VertexCL*> *VertexCopies;
//    ListCL<EdgeCL*>   *EdgeCopies;  // TODO: Wirklich notwendig?
//    ListCL<TetraCL*>  *TetraCopies;
    void Check       () { if (!_IsValid) throw DROPSErrCL("MultiProcCL not yet initialized!"); }

public:
    MultiProcCL() : _NumberProcs(1), _MyProc(0), _IsValid(false) {}
    MultiProcCL(const MultiProcCL&);  // wird nicht implementiert!
    ~MultiProcCL() {}

    void Init        ( int &, char **& );
    void Finalize    ();
    Uint NumberProcs () { return _NumberProcs; }
    Uint GetProc     () { return _MyProc; }
    bool IamRoot     () { return !_MyProc; }
};

extern MultiProcCL MultiProc;


//**************************************************************************
// Class:    IdCL                                                          *
// Purpose:  provides a unique identifier for an object.                   *
// Remarks:  We use the template argument to specify the class whose       *
//           objects will carry an Id.                                     *
//**************************************************************************

template <class type>
class IdCL
{
private:
    static Ulint _Counter;

    Usint _Proc;
    Ulint _Identity;

public:
    IdCL () : _Proc(MultiProc.GetProc()), _Identity(_Counter++) {}
    IdCL (Usint Proc, Ulint Identity) : _Proc(Proc), _Identity(Identity) {}
    // Default Copy-ctor

    Ulint GetCounter () const { return _Counter; }
    Usint GetProc    () const { return _Proc; }
    Ulint GetIdent   () const { return _Identity; }

    bool operator == (const IdCL<type>& Id) const
        { return Id._Identity == _Identity && Id._Proc == _Proc; }
    bool operator != (const IdCL<type>& Id) const { return !(*this==Id); }
    bool operator <  (const IdCL<type>& Id) const
        { return _Identity!=Id._Identity ? _Identity<Id._Identity : _Proc<Id._Proc; }
};


// Initialize the counter only once!
template <class type> Ulint IdCL<type>::_Counter = 0;

} //end of namespace DROPS

#endif
