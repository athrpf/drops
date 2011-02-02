/// \file unknowns.h
/// \brief Implementation of the mapping from simplices to indices to
///    linear-algebra data-structures.
/// \author LNM RWTH Aachen: Sven Gross, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

#ifndef DROPS_UNKNOWNS_H
#define DROPS_UNKNOWNS_H

#include "misc/utils.h"
#include <limits>
#include <vector>


namespace DROPS
{

/// The type of indices used to access the components of vectors,
/// matrices, etc.
typedef Ulint IdxT;

/// \brief Value of all unset indices.
///
/// \note If we ever have as many unknowns, there will of course be a
/// problem -- but probably there would be still one more index and the
/// wraparound would kill us anyways.
const IdxT NoIdx= std::numeric_limits<IdxT>::max();


/// Implementation-detail of UnknownHandleCL.
class UnknownIdxCL
{
  private:
    std::vector<IdxT> _Idx;
#ifdef _PAR
    // This flag array is used for remembering if an unknowns has just been received
    // or if the unknown has been exist before the refinement and migration
    // algorithm has been performed. (sorry for the missleading name giving)
    mutable std::vector<bool> UnkRecieved_;
#endif

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

#ifdef _PAR
    void SetUnkRecv(IdxT i) const
    {
        if (UnkRecieved_.size()<=i)
            UnkRecieved_.resize(i+1,false);
        UnkRecieved_[i]=true;
    }

    bool GetUnkRecv(IdxT i) const
    {
        if (UnkRecieved_.size()<=i)
            return false;
        else
            return UnkRecieved_[i];
    }

    void ResetUnkRecv()
    {
        for (Uint i=0; i<UnkRecieved_.size(); ++i)
            ResetUnkRecv(i);
        UnkRecieved_.resize(0);
    }

    void ResetUnkRecv( IdxT i )
    {
        if (UnkRecieved_.size()<=i)
            return;
        UnkRecieved_[i]=false;
    }

    bool HasUnkRecv() const
    {
        for (Uint i=0; i<UnkRecieved_.size(); ++i)
            if (UnkRecieved_[i])
                return true;
        return false;
    }
#endif
};


/// \brief Maps a simplex and a "sysnum" (system number) on an index for
///     accessing numerical data.
///
/// Every simplex has a public member Unknowns of type UnknownHandleCL,
/// which behaves as a container of indices (for numerical data) that
/// can be accessed via a system number.  This class is only a handle
/// for an UnknownIdxCL.
class UnknownHandleCL
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

    /// True, iff this instance has already acquired an UnknownIdxCL-object.
    bool Exist()             const { return _unk; }
    /// True, iff the system sysnum exists and has a valid index-entry.
    bool Exist( Uint sysnum) const { return _unk && sysnum<_unk->GetNumSystems() && _unk->GetIdx(sysnum) != NoIdx; }

    /// Effectively deletes the index belonging to system sysnum.
    void Invalidate( Uint sysnum) { _unk->GetIdx(sysnum)= NoIdx; }

    UnknownIdxCL* Get() const { return _unk; }
    /// Retrieves the index for sysnum for writing.
    IdxT&        operator() ( Uint i)       { return _unk->GetIdx(i); }
    /// Retrieves the index for sysnum for reading.
    IdxT         operator() ( Uint i) const { return _unk->GetIdx(i); }

    /// Allocates memory for a system with number sysnum.  Afterwards, an index
    /// can be stored for sysnum.
    /// The initial index is set to NoIdx. Thus, .Exist( sysnum)==false.
    void Prepare( Uint sysnum)
    {
        if (!_unk) _unk= new UnknownIdxCL( sysnum+1);
        else if ( !(sysnum < _unk->GetNumSystems()) )
            _unk->resize( sysnum+1, NoIdx);
    }

#ifdef _PAR
    /// Remember if an unknown of an index is just recieved
    void SetUnkRecieved( IdxT i ) const
    {
        Assert(_unk!=0, DROPSErrCL("UnknownHandleCL: Cannot set UnkRecieved before this class is init"), DebugUnknownsC | DebugParallelC);
        _unk->SetUnkRecv(i);
    }
    /// Get information if the unknown of the index is recieved
    bool UnkRecieved( IdxT i ) const
    {
        return _unk && _unk->GetUnkRecv(i);
    }
    /// Forget about the recieved information about all unknowns
    void ResetUnkRecieved() const
    {
        if (_unk!=0)
            _unk->ResetUnkRecv();
    }
    /// Forget about the recieved information about one index
    void ResetUnkRecieved( IdxT i) const
    {
        if (_unk!=0)
            _unk->ResetUnkRecv(i);
    }
    /// For Debugging Purpose: Check if there is an UnkRecv-Flag
    bool HasUnkRecieved() const
    {
        if (_unk==0)
            return false;
        else
            return _unk->HasUnkRecv();
    }
#endif
};

} // end of namespace DROPS

#endif
