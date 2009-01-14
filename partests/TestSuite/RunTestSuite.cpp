//***************************************************************************
// RunTestSuite.cpp 			                                                  *
// Content: This application performs a test defined in the given config.xml*
// Version: 0.1                                                             *
// Date:                                                                    *
// Begin:   30. August 2007                                                 *
//***************************************************************************
/// \author Timo Henrich 
/// \file RunTestSuite.cpp 
/// \brief Testing the proper function of drops.

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <libgen.h>
#include <map>
#include <regex.h>
#include <math.h>
#include <stdlib.h>
#include "./TestSuite.h"
#include "./Reporter.h"
#include <string>
#include "./tinyxml.h"
#include <dirent.h>

using namespace std;



// Basedir of RunTestSuite executable.
string currentBaseDir;

/**
* Returns a string containing the given count of characters starting on the end of str.
*
*/
string getRSubStr(string str,int length)
{
  // ignore negative values
  if(length < 0 ) 
    length = 0 ;
    
   string ret("");
   
   for(int x = str.length()-length ; x >= 0 && x < str.length() ; x++)
   {      
     ret+= str[x];
     
   }
   
   return ret;
}


/**
* This function creates an output based on the generated test data. 
* Depending on the users settings, the output is written to an html-file or
* straight to the screen.
*
* \param A reference to an Reporter Objekt
* \oaram A string with the location of the output-file
*/
void PrepareTestOutput(Reporter* pResult,ofstream &target){
  // Print the test-results
  if(pResult != NULL ) {
  	
 	Reporter pTemp = *pResult;	
 	
  	
	
	// No output-file defined ? Then print the report to the user-screen
	if(!target.is_open()) {

	       Reporter::const_iterator cii;
		   for(cii=pTemp.begin(); cii!=pTemp.end(); cii++)
		   {
		       (*cii)->dump();
		   } // END FOR		 		
		  	
	} 
	// If output-file defined, generate HTML-report and store in this file.
	else {
	
		//target << " Writing results to '"<< outputFile << "' " <<endl;
		
			  // Open-Templates

                string topTplPath = currentBaseDir + "reportTop.tpl";
                string bottomTplPath = currentBaseDir + "reportBtn.tpl";
                
                ifstream tplTop(topTplPath.c_str(),ifstream::in);
                ifstream tplBtn(bottomTplPath.c_str(),ifstream::in);
			  			
			  // Check for errors
			  			  
			   
			  if (!tplTop || !tplBtn)  {                        
			    cerr << "Output generation failed! Templates not found.";
			    exit(-1) ; 			    
			  } // END IF 		

			  // Copy template (TOP) to output-file				  
			  char c ;
			  while((tplTop.get(c))) { target.put(c); } 
		 
		 	  // Datum auslesen  
			  struct tm *localT;
			  time_t currTime;
			
			  time(&currTime);
			  localT = localtime(&currTime);

     		  target << " Report generated on " << (localT->tm_mon+1) << "/" << (localT->tm_mday) << "/" << (localT->tm_year+1900) << " at " << localT->tm_hour << ":" << localT->tm_min << ":" << localT->tm_sec <<endl ;
		 
			  target << "<table width=\"80%\" cellpadding=\"3\" cellspacing=\"0\" border=\"0\" align=\"center\">" << endl;
			  target << "<!-- the document toc --><tr>" << endl ;
     		  target << "<td>" << endl;
     		  target << "<h2>Reports by tested programs:</h2>" <<  endl ;    				        		       	
       		  target << "<ul>" << endl ; 
		
			  // Create repots in HTML for all tested programs
				       Reporter::const_iterator cii;
					
					   // Iterate over all programs to get their names
					   for(cii=pTemp.begin(); cii!=pTemp.end(); cii++)
					   {				       
				     		target << "<li><a href=\"#"<<(*cii)->getLabel()<<"\">" <<(*cii)->getLabel() << "</a></li>" << endl;
					   } // END FOR
					   
			  target << "</ul></td></tr>" << endl << " <!-- every program gets on row --> " << endl;
					   
					   // Iterate over all programs once again to write
					   // their reports to the output-file
					   for(cii=pTemp.begin(); cii!=pTemp.end(); cii++)
					   {
					   		target << "<tr><td>";
					   		target <<  "<a id=\"" << (*cii)->getLabel() <<"\" />" << endl ;
					       (*cii)->writeToHTMLFile(target);
					       target << "</td></tr>";
					   } // END FOR					

			  // Copy template (BOTTOM) to output-file
			  
			  target << "</table>";			  
			  while((tplBtn.get(c))) { target.put(c); } 				
		
		
		
	}

                                                                                                          
  }  // END IF 
} // END OF FUNCTION

// Application entry-point//
int main(int argc, char * argv[]) {

  if(argc < 3 || strcmp(argv[1],"--help") == 0 )  { 
  	cout << "TestSuite for Drops"<<endl;
  	cout << "This program is used to check proper function of drops-applications." << endl <<endl;   	 
  	
  	cout << " Usage: " << argv[0] << "--dir dirname --file filename1 --out filename1 --help" << endl;
  	
  	cout << " Options:"<<endl;
  	cout << "\t--dir\t Performs tests as defined in the XML-files in the supplied directory "<<endl;
  	cout << "\t--file\t Performs tests as defined in the given XML-file "<<endl;
  	cout << "\t--out\t Writes all error and logging messages into the supplied file. "<<endl<<endl;
  	cout << " For more Informationen read the readme.txt file " << endl;  	
  	cout << endl;  	
  	return 1; 
  } // END IF 

  
  // detect current basedir
  currentBaseDir.append(dirname(argv[0]));
  currentBaseDir.append("/");
  
  // parameter indicates whether we test a single file or not.
  bool isLocationFile = true;
  // read command-line parameters
  string testLocation(argv[2]);
  
  if(strcmp(argv[1],"--dir")==0)
  {
    isLocationFile =false;            
    if(testLocation[(testLocation.length()-1)] == '/')
          testLocation =testLocation.substr(0,testLocation.length()-1);    
  }
  
    
    
// Open output-File    
  ofstream target;
  if(argc==5)
  {              
    target.open(argv[4]);
  } 
  
  cout << " Parameters: " ;
  if(isLocationFile)
    cout << " Location-Type: File " << endl;
   else
    cout << " Location-Type: Dir" << endl;
  
   cout <<" Location: "<< testLocation << endl;
   //cout <<" Output: " << argv[4] << endl<<endl;
	
// Perform test(s) for single file or for a whole directory	

  // Single File
  if(isLocationFile)
  {
         cout << " Perform Test-File " << testLocation <<endl;         

        // Create a TestSuite object
        TestSuite currentTest(testLocation);
        
        Reporter* pResult = currentTest.performTest();
        
            cout << "################################################# " << endl;
            cout << " Test finished !" <<endl;
            cout << "################################################# " << endl;
            
        PrepareTestOutput(pResult,target);		  
  }  
  // Directory  
  else
  {
      DIR *hdir;
      struct dirent *entry;

      hdir = opendir(testLocation.c_str());
            
      if(hdir == NULL)
      {
      
        cerr<<"Error reading directory: " << testLocation << endl ;
        exit(-1);
      }      
  
  
     for(struct dirent *entry=readdir(hdir);entry!=NULL ;entry = readdir(hdir))
     {           
          /// \todo Only read files with extension '.xml'
          if(strcmp(entry->d_name,".")==0 || strcmp(entry->d_name,"..")==0)
            continue;
          
         

        
        string tmp(testLocation);
        tmp.append("/");
        tmp.append(entry->d_name);        
        // Skip symbolic-links ('.','..') and all file without '.xml' extension.
                                    
        if(getRSubStr(tmp,1)=="." || getRSubStr(tmp,2)==".." || getRSubStr(tmp,4)!=".xml")        
          continue;
        
        // Create a TestSuite object
        TestSuite currentTest(tmp);                
        
        cout << " Perform Test-File " << testLocation <<endl;         
        
        Reporter* pResult = currentTest.performTest();
        
            cout << "################################################# " << endl;
            cout << " Test finished !" <<endl;
            cout << "################################################# " << endl;
            
            PrepareTestOutput(pResult,target);		                                                             
      };
      closedir(hdir);     
  }
    
  return 1; 	
} // END OF FUNCTION: main
