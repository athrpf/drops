//**************************************************************************
// File:    unknowns.h                                                     *
// Content:                                                                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt IGPM RWTH Aachen     *
// Version: 0.1                                                            *
// History: begin - March, 14 2001                                         *
//                                                                         *
//**************************************************************************


#ifndef DROPS_UNKNOWNS_H
#define DROPS_UNKNOWNS_H

#include "misc/utils.h"
#include <limits>
#include <vector>


namespace DROPS
{


typedef Ulint IdxT;

// Value of all unset indices; if we ever have so many unknowns, there will of course
// be a problem -- but probably there would be still one more index and the wraparound
// would kill us anyways.
const IdxT NoIdx= std::numeric_limits<IdxT>::max();


class UnknownIdxCL
// implementation-detail of UnknownHandleCL
{
  private:
    std::vector<IdxT> _Idx;

  public:
    UnknownIdxCL( Uint numsys) : _Idx( numsys, NoIdx) {}
    UnknownIdxCL() {}
    UnknownIdxCL(const UnknownIdxCL&);
    ~UnknownIdxCL() {}
    UnknownIdxCL& operator=( const UnknownIdxCL&);

    IdxT& GetIdx( Uint sysnum)
    { 
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Idx[sysnum];
    }
    IdxT  GetIdx( Uint sysnum) const
    { 
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Idx[sysnum];
    }
    Uint GetNumSystems()            const { return _Idx.size(); }

    void resize( Uint size, const IdxT defaultIdx= NoIdx)   
    { 
        _Idx.resize(size, defaultIdx);
    }
    
    void push_back(IdxT idx= NoIdx) { _Idx.push_back( idx); }
};


class UnknownHandleCL
// Maps a simplex and a "sysnum" (system number) on an index for
// accessing numerical data.
// Every simplex has a public member Unknowns of type UnknownHandleCL, which
// behaves as a container of indices (for numerical data) that can be accessed
// via a system number.
// This class is only a handle for an UnknownIdxCL.
{
  private:
    UnknownIdxCL* _unk;
  
  public:
    UnknownHandleCL() : _unk(0) {}
    UnknownHandleCL( const UnknownHandleCL& orig)
    {
        _unk= orig._unk ? new UnknownIdxCL( *orig._unk)
                        : 0;
    }
        
    UnknownHandleCL& operator=( const UnknownHandleCL& rhs)
    {
        if (this==&rhs) return *this;
        delete _unk;
        _unk= rhs._unk ? new UnknownIdxCL( *rhs._unk)
                       : 0;
        return *this;
    }
        
    ~UnknownHandleCL() { delete _unk; }
    
    void Init(Uint numsys= 0)   
    { 
        Assert( _unk==0, DROPSErrCL("UnknownHandleCL: Init was called twice"), DebugUnknownsC);
        _unk= new UnknownIdxCL(numsys);
    }

    void Destroy() { delete _unk; _unk= 0; }
    
    // True, iff this instance has already acquired an UnknownIdxCL-object.
    bool Exist()             const { return _unk; }
    // True, iff the system sysnum exists and has a valid index-entry.
    bool Exist( Uint sysnum) const { return sysnum<_unk->GetNumSystems() && _unk->GetIdx(sysnum) != NoIdx; }

    // Effectively deletes the index belonging to system sysnum. 
    void Invalidate( Uint sysnum) { _unk->GetIdx(sysnum)= NoIdx; }

    UnknownIdxCL* Get() const { return _unk; }
    // Retrieves the index for sysnum for writing.
    IdxT&        operator() ( Uint i)       { return _unk->GetIdx(i); }
    // Retrieves the index for sysnum for reading.
    IdxT         operator() ( Uint i) const { return _unk->GetIdx(i); }

    // Allocates memory for a system with number sysnum.  Afterwards, an index
    // can be stored for sysnum.
    // The initial index is set to NoIdx. Thus, .Exist( sysnum)==false.
    void Prepare( Uint sysnum)
    {
        if (!_unk) _unk= new UnknownIdxCL( sysnum+1);
        else if ( !(sysnum < _unk->GetNumSystems()) )
            _unk->resize( sysnum+1, NoIdx); 
    }
};


} // end of namespace DROPS


#endif
