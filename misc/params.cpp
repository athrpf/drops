/// \file params.cpp
/// \brief read parameters from file.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt, Thorolf Schulte ; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/params.h"

namespace DROPS
{
  // =====================================================
  //                    ParamCL
  // =====================================================

  ParamCL::ParamCL() {}

  void ParamCL::open(const std::string path)
  {
    boost::property_tree::read_json(path, this->pt);
  }

  std::istream &operator>>(std::istream& stream, ParamCL& P)
  {
    boost::property_tree::read_json(stream, P.pt);

    return stream;
  }

  std::ostream &operator<<(std::ostream& stream, ParamCL& P)
  {
    return P.print(stream);
  }



  //explicit implementation of get for Point3DCL due to lacking support of this by PropertyTree
  template<>
  Point3DCL ParamCL::get<Point3DCL>(const std::string & pathInPT) const
  {
     Point3DCL point;

     using boost::property_tree::ptree;
     int i=0;

     ptree curPT = this->pt.get_child(pathInPT);
     for (ptree::const_iterator it = curPT.begin(); it != curPT.end(); ++it) {
         point[i++] = it->second.get_value<double>();
     }

     return point;
  }

  std::ostream& ParamCL::print(std::ostream& s)
  {
    using boost::property_tree::ptree;

    for (ptree::const_iterator it = this->pt.begin(); it != this->pt.end(); ++it) {
        s << it->first << ": " << it->second.get_value<std::string>() << "\n";
        print(it->second, std::string("\t"), s);
    }
    return s;
  }

  void ParamCL::print(boost::property_tree::ptree child, std::string level, std::ostream& s)
  {

    using boost::property_tree::ptree;

    for (ptree::const_iterator it = child.begin(); it != child.end(); ++it) {
        s << level << it->first << ": " << it->second.get_value<std::string>() << "\n";
        print(it->second, level+"\t", s);
    }
  }


  // =====================================================
  //                    ReadParamsCL
  // =====================================================

  void ReadParamsCL::SetInfo( std::string& bez, char typ)
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

  void ReadParamsCL::RegInt( int& ref, std::string bez)
  {
      SetInfo( bez, 'i');
      iparam_[bez]= &ref;
  }

  void ReadParamsCL::RegDouble( double& ref, std::string bez)
  {
      SetInfo( bez, 'd');
      dparam_[bez]= &ref;
  }

  void ReadParamsCL::RegCoord( Point3DCL& ref, std::string bez)
  {
      SetInfo( bez, 'c');
      cparam_[bez]= &ref;
  }

  void ReadParamsCL::RegString( std::string& ref, std::string bez)
  {
      SetInfo( bez, 's');
      sparam_[bez]= &ref;
  }

  void ReadParamsCL::RegInt( int& ref, std::string bez, int defaultvalue)
  {
      SetInfo( bez, 'i');
      iparam_[bez]= &ref;
      iparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::RegDouble( double& ref, std::string bez, double defaultvalue)
  {
      SetInfo( bez, 'd');
      dparam_[bez]= &ref;
      dparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::RegCoord( Point3DCL& ref, std::string bez, Point3DCL defaultvalue)
  {
      SetInfo( bez, 'c');
      cparam_[bez]= &ref;
      cparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::RegString( std::string& ref, std::string bez, std::string defaultvalue)
  {
      SetInfo( bez, 's');
      sparam_[bez]= &ref;
      sparam_def_[bez]= defaultvalue;
  }

  void ReadParamsCL::BeginGroup( const std::string& s)
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
      int i= 0; double d= 0; Point3DCL p(0.); std::string s;
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
      UseDefaults();
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
          std::string bez= it->first;
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

  void ReadParamsCL::UseDefaults()
  {
      for (InfoT::const_iterator it= info_.begin(), end= info_.end(); it!=end; ++it)
          if(!it->second.second) // parameter uninitialized
          {
              switch(it->second.first){
                  case 'i':
                  {
                      std::map<std::string,int>::iterator f = iparam_def_.find(it->first);
                      if ( f != iparam_def_.end()){
                       *iparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  case 'd':
                  {
                      std::map<std::string,double>::iterator f = dparam_def_.find(it->first);
                      if ( f != dparam_def_.end()){
                       *dparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  case 'c':
                  {
                      std::map<std::string,Point3DCL>::iterator f = cparam_def_.find(it->first);
                      if ( f != cparam_def_.end()){
                       *cparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  case 's':
                  {
                      std::map<std::string,std::string>::iterator f = sparam_def_.find(it->first);
                      if ( f != sparam_def_.end()){
                       *sparam_[f->first]= f->second;
                       info_[f->first].second= true;
                       std::cout << "WARNING: Parameter " << f->first << " initialized by default value " << f->second << "!\n";
                      }
                      break;
                  }
                  default:
                      break;

              }
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


