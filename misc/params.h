/// \file params.h
/// \brief read parameters from file.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_MISC_PARAMS_H
#define DROPS_MISC_PARAMS_H

#include "misc/container.h"
#include <string>
#include <map>
#include <fstream>

namespace DROPS
{

///   \brief Parser for parameter files used by ParamBaseCL.
///
///   Usage:
///   - register all parameters via register functions RegXXX
///   - read parameters from file via ReadParams
///   \todo Describe syntax of a parameter file
class ReadParamsCL
{
  private:
    typedef std::map<std::string,std::pair<char,bool> > InfoT; // (type,initialized)

    std::string                        group_, name_;
    InfoT                              info_;

    std::map<std::string,int*>         iparam_;
    std::map<std::string,double*>      dparam_;
    std::map<std::string,Point3DCL*>   cparam_;
    std::map<std::string,std::string*> sparam_;
    
    std::map<std::string,int>       iparam_def_;
    std::map<std::string,double>    dparam_def_;
    std::map<std::string,Point3DCL> cparam_def_;
    std::map<std::string,std::string>    sparam_def_;

    bool IsKnown( const std::string& s) { return info_.find( s) != info_.end(); }
    void ReadEntry( std::istream&);
    void ReadData ( std::istream&);
    
    void SetInfo( std::string&, char);
    
    void UseDefaults();
    
    void PrintWarning() const;

  public:
    ReadParamsCL() {}

    /// \name Registration
    //@{
    /// Register parameters under a certain name.
    void RegInt   ( int&,         std::string);
    void RegDouble( double&,      std::string);
    void RegCoord ( Point3DCL&,   std::string);
    void RegString( std::string&, std::string);
    /// Register parameters under a certain name and allow for 
    /// a default initialization 
    void RegInt   ( int&,         std::string,         int defaultvalue);
    void RegDouble( double&,      std::string,      double defaultvalue);
    void RegCoord ( Point3DCL&,   std::string,   Point3DCL defaultvalue);
    void RegString( std::string&, std::string, std::string defaultvalue);
    //@}

    /// \name Groups
    ///
    /// Parameters can be arranged in groups by enclosing Begin/EndGroup-calls
    /// during the registration. Recursive calls lead to nested groups.
    //@{
    void BeginGroup( const std::string& group);
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




