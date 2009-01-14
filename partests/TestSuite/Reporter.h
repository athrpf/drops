//***************************************************************************
// TestSuite.cpp 															*
// Content: this classes are used by the testsuite to generate				* 	
//			reports															*
//			Timo Henrich , SC RWTH Aachen                               	*
// Version: 0.1                                                             *
// Date:                                                                    *
// Begin:   4. September 2007                                               *
//***************************************************************************
/// \author Timo Henrich 
/// \file Reporter.cpp


#include <fstream> 
#include <iostream>
#include <list>
#include <map>
#include <time.h>
#include <deque>
#include <queue>
using namespace std;

#ifndef __REPORTER_H__
#define __REPORTER_H__

#define REPORTER_RES_CODE_OK 1
#define REPORTER_RES_CODE_ERR_EXECUTION 2
#define REPORTER_RES_CODE_ERR_OUTPUT	3
#define REPORTER_RES_CODE_ERR_TESTS		4 


/**
 *  This class represents a single results of a test.
 *
 * 
 */ 
class ReporterTest{
	
 private: 
	map<string,string> settings; // list of all settings used for this test
	map<string,string> checkValues; // list of all checked values
	
	int usedProc;  // number of used processors
	int resultCode; // overall result-code
	
	vector<string> ResultMsg ;
	
 public: 
    /** class constructor
    *
    * \param number of uses processors
    * \param the result-code of this test.
    *
    */
 	ReporterTest(int usedProc,int resultCode) { 		
 		this->usedProc = usedProc; 
  		this->resultCode = resultCode ;
 	} // END OF FUNCTION	
	
	/**
	* Sets another number of used processors
	*/
	void setUsedProc(int usedProc);
	
	/**
	 * Sets another result-code.
	 */
	void setResultCode(int resultCode);

       /**
       * Get current result-code
       */
        int getResultCode();
	
	/**
	 * Inserts a set of new settings data
	 * 
	 * \param Name of this setting
	 * \param Value of this setting
	 */
	void addSettings(const string name,const string value);

	/**
	 * Inserts a set of new settings data
	 * 
	 * \param Name of this setting
	 * \param Value of this setting
	 */	
	void setSettings(map<string,string> settings);
	
	/**
	 * Inserts data about a check performed on a special dataset.
	 * 
	 * \param Name of the checked field
	 * \param Custom message explaining the testresult
	 */
	void addCheckValues(const string label,const string value);
	
	/**
	 * Inserts a general result-message explaining an error.
	 * 
	 * \param A const string containing a message
	 */	
	void addResultMsg(const string msg);
		
	/**
	 *  Dumps the whole data of this object to user-screen.
	 * 
	 */
	void dump ();
	
	
	/**
	 * Writes the whole data of this object to the file ,specified by the given 
	 * file-handle
	 * 
	 * \param A file-handle
	 * 
	 */
	void writeToHTMLFile(ofstream  &pTarget);
		
};  // END OF CLASS
 
/**
 *  This class represents all tests performed on a special program.
 * 
 */
class ReporterProgram  {

private :
  
        string label ;  // label of current program
        string location; // location of the uses executable
        string paramTpl ;  // used param-template
        string verifier;  // location of the program used to verify the tested programs output
        string verifierTpl; // location of the used param-template
  
  
	vector <ReporterTest*> testResults ; 	  // results of each test
	
public: 

	/**
	 * The class-constructor
	 * 
	 * \param A label describing the program
	 * \param The location of the programs-executable
	 * \param The used template for configuration-file generation
	 */
        ReporterProgram(const string label ,const string location ,const string paramTpl,const string  verifier = "" ,const string verifierTpl = "" ) {
                this->label = label; 
		this->location = location;		
		this->paramTpl = paramTpl;		
                this->verifier=     verifier;
                this->verifierTpl = verifierTpl ;                		
	} // END OF FUNCTION
	
	/** class destructor
	*	
	* 	Cleans the memory
	* 
	*/
	~ReporterProgram() {		
	} // END OF FUNCTION
	
	/**
	 * Adds a new set of test-data to this program-report.
	 */
	void addReporterTest(ReporterTest * pRepHandle);
	
	/**
	 * Dumps the whole data of this object to user-screen
	 * 
	 */
	void dump();
	
	/**
	 * Returns the label of the current object
	 */
	const string getLabel();


	/**
	 * Writes the whole data of this object to the file ,specified by the given 
	 * file-handle
	 * 
	 * \param A file-handle
	 * 
	 */
	void writeToHTMLFile(ofstream  &pTarget);
	
}; // END OF CLASS


typedef vector<ReporterProgram*> Reporter;

#endif