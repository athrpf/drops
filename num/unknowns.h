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
#include <vector>


namespace DROPS
{


typedef Ulint IdxT;

class UnknownIdxCL
{
  private:
    std::vector<IdxT*> _Idx;
    std::vector<Uint>  _Len;

  public:
    UnknownIdxCL(Uint numsys)
        : _Idx(numsys), _Len(numsys) {}
    UnknownIdxCL() {}
    UnknownIdxCL(const UnknownIdxCL&);
    ~UnknownIdxCL()
    {
        for( std::vector<IdxT*>::iterator vit=_Idx.begin(); vit!=_Idx.end(); ++vit)
        {
            delete[] *vit;
        }
    }
    UnknownIdxCL& operator=(const UnknownIdxCL&);

    IdxT*& GetIdx(Uint sysnum)
    { 
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Idx[sysnum];
    }
    IdxT*  GetIdx(Uint sysnum) const
    { 
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Idx[sysnum];
    }
    Uint GetNumSystems()            const { return _Idx.size(); }
    Uint GetSystemSize(Uint sysnum) const
    {
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        return _Len[sysnum];
    }
    // Alloc geht davon aus, dass die sysnum schon existiert und fuer diese
    // sysnum noch kein Speicher alloziert worden ist...
    void Alloc(Uint sysnum, Uint size)
    {
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range"), DebugUnknownsC);
        if( !(_Idx[sysnum]= new IdxT[size]) )
            throw DROPSErrCL("Bad Alloc in UnknownIdxCL::Alloc");
        _Len[sysnum]= size;
    }

    void Dealloc(Uint sysnum)
    {
        Assert( sysnum<GetNumSystems(), DROPSErrCL("UnknownIdxCL: Sysnum out of range."), DebugUnknownsC);
        delete[] _Idx[sysnum];
        _Idx[sysnum]= 0;
        _Len[sysnum]= 0;
    }

    void Dealloc(const std::vector<Uint>& systems)
    {
        for (Uint sysnum=0; sysnum<systems.size(); ++sysnum) Dealloc(sysnum);
    }

    void resize(Uint size)   
    { 
        for(Uint sysnum= size; sysnum<_Idx.size(); ++sysnum)
            delete[] _Idx[sysnum];
        _Idx.resize(size); _Len.resize(size); 
    }
    
    void push_back() { _Idx.push_back(0); _Len.push_back(0); }
};


class UnknownHandleCL
{
  private:
    UnknownIdxCL* _unk;
  
  public:
    UnknownHandleCL() : _unk(0) {}
    UnknownHandleCL(const UnknownHandleCL& orig)
    {
        _unk= orig._unk ? new UnknownIdxCL(*orig._unk)
                        : 0;
    }
        
    UnknownHandleCL& operator= (const UnknownHandleCL& rhs)
    {
        if (this==&rhs) return *this;
        delete _unk;
        if (rhs._unk)
            _unk= new UnknownIdxCL(*rhs._unk);
        else
            _unk= 0;
        return *this;
    }
        
    ~UnknownHandleCL() { delete _unk; }
    
    void Init(Uint numsys= 0)   
    { 
        Assert( _unk==0, DROPSErrCL("UnknownHandleCL: Init was called twice"), DebugUnknownsC);
        _unk= new UnknownIdxCL(numsys);
    }

    void Destroy() { delete _unk; _unk= 0; }
    
    bool Exist() const            { return _unk; }
    bool Exist(Uint sysnum) const { return sysnum<_unk->GetNumSystems() && _unk->GetIdx(sysnum); }
    UnknownIdxCL* Get() const { return _unk; }
    IdxT*&        operator() (Uint i)       { return _unk->GetIdx(i); }
    IdxT*         operator() (Uint i) const { return _unk->GetIdx(i); }
    
    void Prepare(Uint sysnum, Uint NumUnknown)
    {
        if (!_unk) _unk= new UnknownIdxCL(sysnum+1);
        else if ( !(sysnum < _unk->GetNumSystems()) )
            _unk->resize(sysnum+1); 
        if ( !Exist(sysnum) )
            _unk->Alloc(sysnum, NumUnknown);
    }
};


} // end of namespace DROPS


#endif
