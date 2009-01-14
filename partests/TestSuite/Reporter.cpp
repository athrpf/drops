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
#include "./Reporter.h"
using namespace std;
        
	/**
	* Sets another number of used processors
	*/
	void ReporterTest::setUsedProc(int usedProc) {
		this->usedProc = usedProc; 		
	} // END OF FUNCTION
	
	/**
	 * Sets another result-code.
	 */
        void ReporterTest::setResultCode(int resultCode) {
		this->resultCode = resultCode ;							
	} // END OF FUNCTION

       /**
       * Get current result-code
       */
        int ReporterTest::getResultCode() {
          return this->resultCode;
        } // END OF FUNCTION
        
	/**
	 * Inserts a set of new settings data
	 * 
	 * \param Name of this setting
	 * \param Value of this setting
	 */
        void ReporterTest::addSettings(const string name,const string value) {
		this->settings.insert(map<string,string>::value_type(name,value));
	} // END OF FUNCTION

	/**
	 * Inserts a set of new settings data
	 * 
	 * \param Name of this setting
	 * \param Value of this setting
	 */	
        void ReporterTest::setSettings(map<string,string> settings) {
		this->settings = settings; 	
	} // END OF FUNCTION
	
	/**
	 * Inserts data about a check performed on a special dataset.
	 * 
	 * \param Name of the checked field
	 * \param Custom message explaining the testresult
	 */
        void ReporterTest::addCheckValues(const string label,const string value) {
		this->checkValues.insert(map<string,string>::value_type(label,value));
	} // END OF FUNCTION
	
	/**
	 * Inserts a general result-message explaining an error.
	 * 
	 * \param A const string containing a message
	 */	
        void ReporterTest::addResultMsg(const string msg) {
		//string tmp(msg);
                this->ResultMsg.push_back(msg);
	} // END OF FUNCTION	
		
	/**
	 *  Dumps the whole data of this object to user-screen.
	 * 
	 */
        void ReporterTest::dump () {
		cout << endl << "************************" <<endl;
		cout << " Result-Code: " << this->resultCode <<endl;
		cout << " Processors:	" << this->usedProc <<endl;
		cout << " Used Settings : "<< endl;

		   // Iterate over all settings
		   for( map<string,string>::iterator ii=this->settings.begin(); ii!=this->settings.end(); ++ii)
		   {
		       cout << " \t "<< (*ii).first << ": " << (*ii).second << endl;
		   } // END FOR		
	   cout << endl ;	   
	   cout << " Checked Values " << endl ;
	   
		   // Iterate over all checked fields
		   for( map<string,string>::iterator ii=this->checkValues.begin(); ii!=this->checkValues.end(); ++ii)
		   {
		       cout << " \t "<< (*ii).first << ": " << (*ii).second << endl;
		   } // END FOR				   
	   cout << endl ;
	   	   
		   // List all messages
	   cout << " Result-Messages " << endl ;
	       vector<string>::const_iterator cii;
		   for(cii=this->ResultMsg.begin(); cii!=this->ResultMsg.end(); cii++)
		   {
		      cout << "\t" << *cii << endl;
		   } // END FOR		   
		 		
	} // END OF FUNCTION	
	
	
	/**
	 * Writes the whole data of this object to the file ,specified by the given 
	 * file-handle
	 * 
	 * \param A file-handle
	 * 
	 */
        void ReporterTest::writeToHTMLFile(ofstream  &pTarget) {
		
 			pTarget << "<tr><td class=\"testContent\">" << endl ;
  		    pTarget << "<b>Used Processors:</b> " << this->usedProc  << " <br/>" << endl;
  		    pTarget << "<b>Result-Code:</b> "<< this->resultCode << " <br>" << endl ; 
			pTarget << "<b>Settings:</b>" <<endl;		
		
		   // Iterate over all settings
		   pTarget << "<table class=\"settingsTable\"> " << endl;
		   for( map<string,string>::iterator ii=this->settings.begin(); ii!=this->settings.end(); ++ii)
		   {
		       pTarget  << "<tr><td>  "<< (*ii).first << "</td><td>" << ": " << (*ii).second <<"</td></tr>" <<endl;
		   } // END FOR		
	   	   pTarget << "</table>"<< endl ;		
	   	   
	   	   if(!this->checkValues.empty()) {
			   	   pTarget << "<b>Checked Values:</b>" << endl ; 	   	   	   	   
				   // Iterate over all checked fields
				   pTarget << "<table cellspacing=\"2\" border=\"1\" cellpadding=\"2\" class=\"checkValuesTable\">" <<endl;		   
				   for( map<string,string>::iterator ii=this->checkValues.begin(); ii!=this->checkValues.end(); ++ii)
				   {
				       pTarget << "<tr><td valign=\"top\" width=\"200px\">"<< (*ii).first << "</td><td> " << (*ii).second <<"</td></tr>"<< endl;
				   } // END FOR		   	   
				   pTarget << "</table>" <<endl ;
	   	   } // END IF 
	   	   
		   // List all messages
	   	   pTarget << "<b>Result-Messages</b> " << endl ;	
	   	   
	   	   if(this->resultCode == REPORTER_RES_CODE_OK ) {       
	       	pTarget << "<table class=\"testPassed\">" << endl ;
	       	pTarget << "<tr><td>Test passed.</td></tr>" << endl;
	       	pTarget << "</table></td></tr>";
	   	   } else {
		   	   	//pTarget << "<table class=\"testFailed\">" << endl ;	
		   	   
		   	   pTarget << "<textarea readonly=\"readonly\" class=\"testFailed\">";
		   	   
		       vector<string>::const_iterator cii; 
			   for(cii=this->ResultMsg.begin(); cii!=this->ResultMsg.end(); cii++)
			   {
			     // pTarget << "<tr><td>" << *cii << "</td></td>"<< endl;
			     pTarget << *cii << "\n" << endl;
			   } // END FOR		
			   
			   	   	 
			   pTarget << "</textarea>";
			   
			   //pTarget << "</table></td></tr>";  	
	   	   }
		
		
	} // END OF FUNCTION
		         
	
	/**
	 * Adds a new set of test-data to this program-report.
	 */
        void ReporterProgram::addReporterTest(ReporterTest * pRepHandle) {
		this->testResults.push_back(pRepHandle);
	} // END OF FUNCTION
	
	/**
	 * Dumps the whole data of this object to user-screen
	 * 
	 */
        void ReporterProgram::dump() {
		cout << "#################################################################" <<endl;
		cout << "\t\t"<< this->label  <<endl;
		cout << "#################################################################" <<endl;
		cout << " Location: " << this->location << endl;
		cout << " Conf-Templ.:"<< this->paramTpl <<endl;
                if(this->verifier.length()>0)                
                  cout <<" Verifier: " << this->verifier << endl;
                if(this->verifierTpl.length()>0 )
                  cout <<" Conf-Templ: " << this->verifierTpl << endl;
                                
		cout << endl;
		cout << "Testresults:" <<endl;
		
	       vector<ReporterTest*>::const_iterator cii;
		   for(cii=this->testResults.begin(); cii!=this->testResults.end(); cii++)
		   {
		      (*cii)->dump();
		   } // END FOR		 		
		
			
	} // END OF FUNCTION
	
	/**
	 * Returns the label of the current object
	 */
        const string ReporterProgram::getLabel() {
			return this->label ;
	} // END OF FUNCTION


	/**
	 * Writes the whole data of this object to the file ,specified by the given 
	 * file-handle
	 * 
	 * \param A file-handle
	 * 
	 */
        void ReporterProgram::writeToHTMLFile(ofstream  &pTarget) {
	
	// Write common program information	
	pTarget << "<table width=\"95%\" align=\"right\" cellspacing=\"0\" cellpadding=\"3\" border=\"0\">"<<endl;
  	pTarget << "<tr><td class=\"programHeader\">" << this->label << "</td></tr>" << endl ;  			
	pTarget << "<tr><td class=\"programContent\">" << endl;
  	pTarget << "<!-- Common program information -->" << endl ;
    pTarget << "<b>Executable:</b>"<<  this->location << "<br/>" << endl;
    pTarget << "<b>Configuration Template:</b>" << this->paramTpl << "<br/>" << endl ;
    if(this->verifier.length()>0 )
      pTarget << "<b>Verifier:</b>"<<  this->verifier << "<br/>" << endl;
    if(this->verifierTpl.length()>0 )
      pTarget << "<b>Verifier Configuration:</b>"<<  this->verifierTpl << "<br/>" << endl;
      
	pTarget << "<!-- Test Toc -->" <<endl ; 
	pTarget << "<h4>Performed Tests</h4><ul>" << endl ;
		   // Create TOC
		   for(unsigned int i = 1 ; i <= this->testResults.size(); i++ )
		   {
		   		pTarget << "<li><a href=\"#"<< this->label << i << "\">Test " << i << "</a></li>" << endl;
		   } // END FOR		
						  			  			
  	pTarget << "</ul>" << endl ;
  	
  		  int testCnt = 1 ;
  		  
  		  // Iterate over all performed tests an write their results to output-file
	       vector<ReporterTest*>::const_iterator cii;
	       pTarget << "<table width=\"95%\" align=\"right\" cellspacing=\"0\" cellpadding=\"3\" border=\"0\"> " << endl;
		   for(cii=this->testResults.begin(); cii!=this->testResults.end(); cii++)
		   {
		   	  		                             
                          if((*cii)->getResultCode()==1)                            
		   	    pTarget << "<tr><td class=\"testHeader\">Test "<< testCnt ;
                          else
                            pTarget << "<tr><td class=\"testHeaderErr\">Test "<< testCnt ;
                          
		   	  pTarget <<  "<a id=\"" << this->label << testCnt <<"\" />" << endl ;
		   	  pTarget <<"</td></tr><tr>" << endl ;
   	 
		   	  
		      (*cii)->writeToHTMLFile(pTarget);
		      testCnt ++;
		      pTarget << "<!-- emptry row --><tr><td>&nbsp;</td></tr>" << endl ;
		   } // END FOR
		   
		   pTarget << " </table>" << endl ;		  		  
  				
  	pTarget << "<!-- emptry row --><tr><td>&nbsp;</td></tr>" << endl ; 		
	pTarget << "</td></tr>" << endl ;
	pTarget << "</table>"<<endl ;  		    		
	}