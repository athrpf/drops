/// \file logger.cpp
/// \brief log certain values like e.g. timespan measured during discretization and integer values etc.
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Timo Henrich

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

#include "parallel/logger.h"
#include <stdio.h>
#include <stdlib.h>

namespace DROPS
{

// Define used map
std::map<std::string,std::vector<double> > LoggerCL::storedVars;
Uint LoggerCL::stepCnt=1;

void LoggerCL::setValue(const std::string& varName, const double& varValue)
{
    if(LoggerCL::storedVars[varName].size()==LoggerCL::stepCnt)
    {
        LoggerCL::storedVars[varName][LoggerCL::stepCnt]=varValue;
    }
    else
    {
        LoggerCL::storedVars[varName].push_back(varValue);
    }
}

void LoggerCL::addValue(const std::string& varName, const double& varValue)
{
    storedVars[varName][(LoggerCL::stepCnt-1)]+= varValue;
}

double LoggerCL::getValue(const std::string& varName)
{
    return LoggerCL::storedVars[varName][(LoggerCL::stepCnt-1)];
}

void LoggerCL::writeOut(const std::string& fileName,bool printHeader)
{
    // Open file for writting
    char procNumStr[32];
    sprintf(procNumStr,"%d",ProcCL::Size());
    std::string procStr(&procNumStr[0]);
    procStr.append(".logger");
    std::ofstream outFile((fileName+procStr).c_str());

    if(!outFile.is_open())
    {
       std::cout << " Error: Can not open file '"<<fileName<<"' for writting!"<<std::endl;
       return;
    }


    // Write out comments - if wished
    // Generiere Kopfzeile
    if(printHeader)
    {
        outFile<<"%Procs ";
        for(std::map<std::string,std::vector<double> >::iterator it = storedVars.begin() ; it != storedVars.end() ; ++it)
        {
            if((*it).second.size()==1)
                outFile<< (*it).first <<" ";
            else
                for(Uint i=1; i <= (*it).second.size() ; ++i)
                {
                 outFile<< (*it).first<<i<<" ";
                }
        }
        outFile<<std::endl;
    }


    // Write out all saved data

    // Schreibe Anzahl der Prozessoren aus
    outFile<<ProcCL::Size()<<" ";

    for(std::map<std::string,std::vector<double> >::iterator it = LoggerCL::storedVars.begin() ; it != LoggerCL::storedVars.end() ; it++)
    {
        // Write out all data-sets for this label
        for(std::vector<double>::iterator ix = (*it).second.begin(); ix != (*it).second.end(); ++ix)
        {
            outFile<< (*ix) <<" ";
        }
    }

    outFile<<std::endl;

    outFile.close();
}

void LoggerCL::nextStep()
{
    LoggerCL::stepCnt++;
}

void LoggerCL::clear()
{
    LoggerCL::storedVars.clear();
}

} // END OF Namespace::Drops
