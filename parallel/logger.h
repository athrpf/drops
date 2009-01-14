/**
 * File: logger.h
 * Content: Definition of the logger-class. This class is used to log certain values like
 *          e.g. timespan measured during discretization etc.
 * Author:  Timo Henrich , SC RWTH Aachen
 * Version: 1.0
 *
 */

#ifndef DROPS_LOGGER_H
#define DROPS_LOGGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

#include "geom/multigrid.h"
#include "geom/topo.h"
#include "misc/utils.h"
#include "misc/problem.h"

// "Log" is used to toggle whether the logger is active or not.
#ifdef _LOG
	#define DROPS_LOGGER_SETVALUE(A,B)		DROPS::LoggerCL::setValue(A,B);
	#define DROPS_LOGGER_ADDVALUE(A,B)		DROPS::LoggerCL::addValue(A,B);
	#define DROPS_LOGGER_GETVALUE(A)		DROPS::LoggerCL::getValue(A);
	#define DROPS_LOGGER_WRITEOUT(A,B)		DROPS::LoggerCL::writeOut(A,B);
	#define DROPS_LOGGER_NEXTSTEP()			DROPS::LoggerCL::nextStep();
	#define DROPS_LOGGER_CLEAR()			DROPS::LoggerCL::clear();
#else
	#define DROPS_LOGGER_SETVALUE(A,B)	
	#define DROPS_LOGGER_ADDVALUE(A,B)	
	#define DROPS_LOGGER_GETVALUE(A)	
	#define DROPS_LOGGER_WRITEOUT(A,B)	
	#define DROPS_LOGGER_NEXTSTEP()		
	#define DROPS_LOGGER_CLEAR()		
#endif

namespace DROPS
{
    
class LoggerCL
{
    private:
    // Contains all values stored by setValue-Method
	static std::map<std::string,std::vector<double> > storedVars;
	// Current step in simulation
	static Uint stepCnt;
        
    public:

	/**
	* Stores a measured value for the current timestep under a
	* a specific name.
	* If a value is only set in once during a simulation it is treated like a
	* global variable. It appears in the column named "varName" in the log-file.
	* If a value is set multiple time in differen timesteps, a column is created for
	* each timestep name "varName[stepCnt]" (without bracket).
	*
	* \param varName A name for the value to be stored.
	* \param varValue The value to be store.
	*/
	static void setValue(const std::string& varName, const double& varValue);

	/**
	* Increments the current value of the logged variable "varName" by the passed
	* value.
	* If varName was not set before, it is initialized with zero.
	*
	* \param varName The name of the target log-entry.
	* \param varValue The increment value.
	*/
	static void addValue(const std::string& varName, const double& varValue);	

	/**
	* Returns the value stored under a specific varName.
	*
	* \return The value stored under the given varName. If no value was stored, zero is returned.
	*/
	static double getValue(const std::string& varName);	
	
	/**
	* Writes out the currently logged variables.
	* Note: The first column in output-file is always the number of used processors during simulation.
	*		The following columns are ordered in aplhabetic sequence according to their varNames.
	*
	* \param fileName A prefix-name for the logging-file. The logged-data is written into a file name
	*		 [PREFIX][NUM OF PROCS]
	* \param printHeader Toggles whether a comment with all varnames is written into the file or not.
	*/
	static void writeOut(const std::string& fileName,bool printHeader=false);
	
	/**
	* Notices the logger that a new simulation step has been reached.	
	*/
	static void nextStep();

	/*
	* Clears up the internal cache.
	*/
    static void clear();
            
}; // END OF CLASS

}

#endif
