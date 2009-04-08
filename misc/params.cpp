/// \file
/// \brief read parameters from file.

#include "misc/params.h"

namespace DROPS
{

// =====================================================
//                    ReadParamsCL
// =====================================================

void ReadParamsCL::SetInfo( string& bez, char typ)
{
    bez= group_ + bez;
    if (!IsKnown(bez))
    {
        info_[bez].first= typ;
        info_[bez].second= false;
    }
    else
        throw DROPSErrCL( "ReadParamsCL: Parameter "+bez+" already registered");
}

void ReadParamsCL::RegInt( int& ref, string bez)
{
    SetInfo( bez, 'i');
    iparam_[bez]= &ref;
}

void ReadParamsCL::RegDouble( double& ref, string bez)
{
    SetInfo( bez, 'd');
    dparam_[bez]= &ref;
}

void ReadParamsCL::RegCoord( Point3DCL& ref, string bez)
{
    SetInfo( bez, 'c');
    cparam_[bez]= &ref;
}

void ReadParamsCL::RegString( string& ref, string bez)
{
    SetInfo( bez, 's');
    sparam_[bez]= &ref;
}

void ReadParamsCL::BeginGroup( const string& s)
{
    group_+= s;
    group_+= ":";
}

void ReadParamsCL::EndGroup()
{
    if (group_.empty())
        throw DROPSErrCL("ReadParamsCL::EndGroup: missing BeginGroup!");
    int pos= group_.rfind( ':', group_.size()-2);
    group_.erase( pos+1);
}


inline void SkipComment( std::istream& s)
{
  if (!s.eof()) s >> std::ws;
  while (!s.eof() && s.peek()=='#')
  {
    while (!s.eof() && s.get()!='\n') {}
    if (!s.eof()) s >> std::ws;
  }
}

void ReadParamsCL::ReadEntry( std::istream& is)
{
    SkipComment( is);
    while (!is.eof())
    {
        switch(is.peek())
        {
          case '#': SkipComment( is); break;
          case '=': is.get(); SkipComment( is);
                    name_= group_ + name_;
                    ReadData( is);
                    name_.clear();
                    return;
          case '{': is.get(); if (!name_.empty())
                               { BeginGroup( name_); name_.clear(); }
                               else throw DROPSErrCL("ReadParamsCL: group name missing before '{'");
                               break;
          case '}': is.get(); if (name_.empty()) EndGroup();
                               else throw DROPSErrCL("ReadParamsCL: end of group "+group_+" before defining parameter "+name_);
                               break;
          default: name_+= is.get(); break;
        }
        if (!is.eof()) is >> std::ws;
    }
}

void ReadParamsCL::ReadData( std::istream& is)
{
    if (!IsKnown( name_))
    {
        std::cout << "Skipping unknown parameter " << name_ << std::endl;
        is.ignore( 256, '\n');
        return;
    }
    int i= 0; double d= 0; Point3DCL p(0.); string s;
    char typ= info_[name_].first;
    switch (typ)
    {
        case 'i': is >> i; break;
        case 'd': is >> d; break;
        case 'c': is >> p[0] >> p[1] >> p[2]; break;
        case 's': is >> s;
    }
    if (!is)
    {
        throw DROPSErrCL( "ReadParamsCL: reading of data failed for parameter "+name_);
    }
    switch (typ)
    {
        case 'i': *iparam_[name_]= i; break;
        case 'd': *dparam_[name_]= d; break;
        case 'c': *cparam_[name_]= p; break;
        case 's': *sparam_[name_]= s; break;
    }
    info_[name_].second= true;
}

void ReadParamsCL::Clear()
{
    info_.clear(); iparam_.clear(); dparam_.clear(); cparam_.clear(); sparam_.clear();
}

void ReadParamsCL::ReadParams( std::istream& is)
{
    if (!group_.empty())
        throw DROPSErrCL("ReadParamCL: group "+group_+" not terminated properly");

    if (!is)
        throw DROPSErrCL("ReadParamsCL: file error");

    name_.clear(); group_.clear();
    // reset info: all parameters uninitialized
    for (InfoT::iterator it= info_.begin(), end= info_.end(); it!=end; ++it)
        it->second.second= false;

    while (!is.eof())
    {
        ReadEntry( is);
        if (!name_.empty())
            throw DROPSErrCL("ReadParamsCL: error while reading parameter "+name_);
    }
    PrintWarning();
    if (!group_.empty())
        throw DROPSErrCL("ReadParamCL: group "+group_+" not terminated properly");
}

/// \note This routine can be used to generate a standard parameter file
///
void ReadParamsCL::WriteParams( std::ostream& os) const
{
    os << "#=============================================\n"
       << "#    DROPS parameter file\n"
       << "#=============================================\n\n";
    for (InfoT::const_iterator it= info_.begin(), end= info_.end();
        it!=end; ++it)
    {
        string bez= it->first;
        char typ= it->second.first;
        os << bez << "\t=\t";
        switch(typ)
        {
          case 'i': os << *(iparam_.find(bez)->second); break;
          case 'd': os << *(dparam_.find(bez)->second); break;
          case 'c': os << *(cparam_.find(bez)->second); break;
          case 's': os << *(sparam_.find(bez)->second); break;
        }
        os << '\n';
    }
}

void ReadParamsCL::PrintWarning() const
{
    bool AllOk= true;
    for (InfoT::const_iterator it= info_.begin(), end= info_.end(); it!=end; ++it)
        if(!it->second.second) // parameter uninitialized
        {
            std::cout << "WARNING: Parameter " << it->first << " uninitialized!\n";
            AllOk= false;
        }
    if (!AllOk) throw DROPSErrCL("ReadParamsCL: Parameters above are missing in parameter file!");
}

} // end of namespace DROPS



