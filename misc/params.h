//**************************************************************************
// File:    params.h                                                       *
// Content: read parameters from file                                      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_MISC_PARAMS_H
#define DROPS_MISC_PARAMS_H

#include "misc/container.h"
#include <string>
#include <map>
#include <fstream>

namespace DROPS
{

using std::string;

//   ReadParamsCL: Parser for parameter files. 
//   Usage:
//   * register all parameters via RegXXX
//   * read parameters from file via ReadParams
class ReadParamsCL
{
  private:
    typedef std::map<string,std::pair<char,bool> > InfoT; // (type,initialized)
  
    string                      group_, name_;
    InfoT                       info_;    

    std::map<string,int*>       iparam_;
    std::map<string,double*>    dparam_;
    std::map<string,Point3DCL*> cparam_;
    std::map<string,string*>    sparam_;
    
    bool IsKnown( const string& s) { return info_.find( s) != info_.end(); }
    void ReadEntry( std::istream&);
    void ReadData ( std::istream&);
    void SetInfo( string&, char);
    void PrintWarning() const;
    
  public:
    ReadParamsCL() {}
    
    // register parameters
    void RegInt   ( int&,       string);
    void RegDouble( double&,    string);
    void RegCoord ( Point3DCL&, string);
    void RegString( string&,    string);
    
    // build groups (may be used recursively)
    void BeginGroup( const string& group);
    void EndGroup();
    
    // IO
    void ReadParams( std::istream&);
    void WriteParams( std::ostream&) const;
    
    //cleanup
    void Clear(); // deallocate memory
};

class ParamBaseCL
{
  protected:
    ReadParamsCL rp_;
    
  public:
    void Clear() { rp_.Clear(); }		// cleanup

    friend std::istream& operator >> ( std::istream& is, ParamBaseCL& P)       { P.rp_.ReadParams(  is); return is; }
    friend std::ostream& operator << ( std::ostream& os, const ParamBaseCL& P) { P.rp_.WriteParams( os); return os; }
};

} // end of namespace DROPS

#endif
    
    
    
    
