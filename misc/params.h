/// \file
/// \brief read parameters from file.

#ifndef DROPS_MISC_PARAMS_H
#define DROPS_MISC_PARAMS_H

#include "misc/container.h"
#include <string>
#include <map>
#include <fstream>

namespace DROPS
{

using std::string;

///   \brief Parser for parameter files used by ParamBaseCL. 
///
///   Usage:
///   - register all parameters via register functions RegXXX
///   - read parameters from file via ReadParams
///   \todo Describe syntax of a parameter file
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
    
    /// \name Registration
    //@{
    /// Register parameters under a certain name.
    void RegInt   ( int&,       string);
    void RegDouble( double&,    string);
    void RegCoord ( Point3DCL&, string);
    void RegString( string&,    string);
    //@}
    
    /// \name Groups 
    ///
    /// Parameters can be arranged in groups by enclosing Begin/EndGroup-calls 
    /// during the registration. Recursive calls lead to nested groups.
    //@{
    void BeginGroup( const string& group);
    void EndGroup();
    //@}
    
    /// \name Input/Output
    //@{
    void ReadParams( std::istream&);
    void WriteParams( std::ostream&) const;
    //@}
    
    /// Cleanup: deallocate memory 
    void Clear();
};

/// \brief Base class for parameter classes
///
/// All problem dependent parameter classes should be derivated from this class.
/// For an example see e.g. ParamMesszelleCL in levelset/params.h
class ParamBaseCL
{
  protected:
    /// Does all the work: registrating parameters, parsing parameter file
    ReadParamsCL rp_;
    
  public:
    /// Cleanup: deallocate memory
    void Clear() { rp_.Clear(); }

    friend std::istream& operator >> ( std::istream& is, ParamBaseCL& P)       { P.rp_.ReadParams(  is); return is; }
    friend std::ostream& operator << ( std::ostream& os, const ParamBaseCL& P) { P.rp_.WriteParams( os); return os; }
};

} // end of namespace DROPS

#endif
    
    
    
    
