//***************************************************************************
// TestSuite.cpp                                                            *
// Content: this class perform some test as defined in the given config-fiel*
//                      and prints the results                              *
//                      Timo Henrich , SC RWTH Aachen                       *            
// Version: 0.1                                                             *
// Date:                                                                    *
// Begin:   30. August 2007                                                 *
// Main class of the application. Is used to perform tests and generate     *
// reports to tty or files                                                  *
//                                                                          *
//***************************************************************************

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <map>
#include <regex.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include "./tinyxml.h"
#include "./TestSuite.h"
#include "./Reporter.h"


using namespace std;


    // Helper functions. Used to perform a case insensitive string search
    bool ci_equal(char ch1, char ch2)
    {
    return toupper((unsigned char)ch1) == toupper((unsigned char)ch2);
    }

    size_t ci_find(const string& str1, const string& str2)
    {    
    string::const_iterator pos = search(str1. begin ( ), str1. end ( ), str2.
    begin ( ), str2. end ( ), ci_equal);
    if (pos == str1. end ( ))
    return string::npos;
    else
    return pos - str1. begin ( );
    }



     /**
     * Sets a new configuration-file
     *
     * \param The path to an configuration-file
     */
    void TestSuite::setConfigFile(const char * configFile) {
      if(configFile != NULL ) {
      
        this->configFile = configFile ; 
        /*
        delete(this->configFile);
        this->configFile = new char[strlen(configFile)];
        strcpy(this->configFile,configFile);
        */
      }  // END IF
    } // END OF FUNCTION

    /**
     * \short Checks a given string (the output of a executed program) for the occurence of certain phrases.
     * \param &programOutput A reference to the output string.
     * \param pattern An Array with patterns to be searched for.
     * \return A negative integer, if no pattern was found in the output-string.
     *        Otherwise a positive integer, indication the index of the found phrase in the passed vector.
     */
    int TestSuite::checkOutputForPhrase(const string &programOutput,const vector<string> &patterns)
    {      
      for(unsigned int i = 0 ; i < patterns.size() ; i++)
      {
      
        
        //if(programOutput.find(patterns[i])!=string::npos)               
        if(ci_find(programOutput,patterns[i])!=string::npos)
        {
          return i;
        }
      } // END FOR
      return -1;
    }    
        
                        
     /**
     *  Performs a system test like specified in the given config-file.
     *
     * If successfull,a pointer to a vector of reports is returned. In case of
     * error 'NULL' is returned.
     */
    Reporter* TestSuite::performTest() {
                
      this->lastReport = new Reporter();
                  
      cout << "################################################# " << endl;
      cout << " Perform Test " <<endl;
      cout << "################################################# " << endl;
                                
      // Try to read given configuration-file
      TiXmlDocument doc(this->configFile);
                                
      bool loadOkay = doc.LoadFile();
      // check whether file could be load or not.
      if ( !loadOkay )
      {
        cout << "Could not load the config-file '"<< configFile << "'. Error='"<<doc.ErrorDesc() <<"'"<<endl ;
        return NULL;
      } // END IF
                                
      // Compile the regex used for mode parsing
      regex_t modeCompiledExp ;
      regcomp (&modeCompiledExp, "([a-zA-Z0-9]*)(\\([.0-9]*\\)){0,1}", REG_EXTENDED);
                                
      // Get root element of configuration-file
      TiXmlElement  * docRoot = doc.RootElement();
      int testSetCnt = 1 ;
                                                                
      // iterate over all defined programs

      const char *  pMpiRunner = docRoot->Attribute("mprunner");
                                                                                              
      for( TiXmlElement* aElement = docRoot->FirstChildElement("program"); aElement; aElement = aElement->NextSiblingElement("program") ){
                                                                                
        // read common data

        const char * pProgLabel             =  aElement->Attribute("label");
        const char * pProgConfigTpl         =  aElement->Attribute("paramTpl");
        const char * pExecPath              =  aElement->Attribute("exec");
        char * pErrorKeywords         =(char*) aElement->Attribute("errorKeywords");
               
        string pVerifierPath("");
        string pVerifierConfigTpl("");
        
        if(aElement->Attribute("verifier")!=NULL)
          pVerifierPath.append(aElement->Attribute("verifier"));
        
        if(aElement->Attribute("verifierParamTpl")!=NULL)
          pVerifierConfigTpl.append(aElement->Attribute("verifierParamTpl"));
        
        // Define error-phrases
        // If any of the following phrases occures in the programs or the verifiers output,
        // TestSuite supposes there is something wrong.
        
          vector<string> errorPatterns;    
          
          // If no error-keywords are supplied, the standard-set is used
          if(pErrorKeywords==NULL)   
          {                     
            errorPatterns.push_back("fehler");                    
            errorPatterns.push_back("error");          
            errorPatterns.push_back("Errno");          
            errorPatterns.push_back("exit");
          }
          // otherwise use the keywords defined by user
          else
          {
          char * tmpKeyWord = strtok (pErrorKeywords," ;,");          
            while (tmpKeyWord != NULL)
            {                                    
               errorPatterns.push_back(tmpKeyWord);            
              // Parse string to get next keyword
              tmpKeyWord = strtok (NULL, " ;,");          
            }
          } // END IF                                            
        // Create a new section for the overall-report.
        ReporterProgram  * pCurrentProgReport = new ReporterProgram(pProgLabel,pExecPath,pProgConfigTpl,pVerifierPath,pVerifierConfigTpl);        
        this->lastReport->push_back(pCurrentProgReport);
                                                                                         
        if(pProgLabel == NULL ) {
          cerr << " Missing 'label' in 'program'-tag " << endl ;
          return NULL;
        } // END IF
                                        
        if(pExecPath == NULL ) {
          cerr << " Missing 'exec' in 'program'-tag " << endl;
          return NULL;
        }
                                                                                      
        int testCnt = 1 ;

        // Perform all defined tests for this program
        for( TiXmlElement* bElement = aElement->FirstChildElement("testset"); bElement; bElement = bElement->NextSiblingElement("testset") ){
                                                
          cout<< endl << "\tTestset " << testSetCnt << endl ;
          testSetCnt ++;
                                                   // read processor attribute 
           char *  pProcStr =(char*) bElement->Attribute("processors");
		  
		   
		  char  * pVeriProcStr=NULL;

		   if(bElement->Attribute("verifierProcessors")!=NULL)
			    pVeriProcStr= ((char*) bElement->Attribute("verifierProcessors"));
		   
                                                                                                                                                     
                                                  // create temp config-file
          if(pProgConfigTpl !=NULL )
            cout << "\t\tPreparing configuration-file" << endl ;          
          else      
            cout << "\t\tNo configuration-file defined." << endl ;
                                                                                                                                                                                                                                                                                                                                   
          // Generate config-files
          string pTempConfig, verififierConfig;
          TiXmlElement * config = bElement->FirstChildElement("config");
                             
          if(pProgConfigTpl != NULL )
            this->buildConfigFromTpl(bElement,pProgConfigTpl,pTempConfig );
                  
          if(!pVerifierConfigTpl.empty())
          {           
            this->buildConfigFromTpl(bElement,pVerifierConfigTpl,verififierConfig);            
          }                                                                                                              
          
          
          
          // Read all settings
          map<string,string> settings;
          for( TiXmlElement* sett = config->FirstChildElement("setting"); sett; sett = sett->NextSiblingElement("setting") ){
            string pSettLabel(sett->Attribute("label"));
            string pSettValue(sett->GetText());
            settings.insert(map<string,string>::value_type(pSettLabel,pSettValue));
          }
                                             
          //  Perform tests for all given processor numbers
          char * pch = strtok (pProcStr," ;,");
		  // If verifier is used, parse given processor-numbers
		  char * pvch = strtok(pVeriProcStr,";,");
          bool serFlag = true;
          while (pch != NULL || pProcStr== NULL&&serFlag  )
          {
            serFlag = false;
            int execResult=0;
            short int execLevel = 0 ;
            bool skipStep = false;
           
            // Prepare object that holds all logging-information
            int proc;
            if(pProcStr==NULL)
              proc =  1;
            else
              proc = atoi(pch);

            // Check if an mpi-runner is set.
            if(pProcStr!=NULL && pMpiRunner==NULL)
            {
              cerr << "Performing tests with mpi on different processors requires the 'mprunner' attribute to be set. "<<endl;
              return NULL;
            }            
            
            ReporterTest * pCurrentTest = new ReporterTest(proc,1);
            pCurrentTest->setSettings(settings);
            pCurrentProgReport->addReporterTest(pCurrentTest);

            // Keeps output of the last program executed.                     
            string programOutput("");
            
            //------------------------------------------------
            //  Test program
            //
            //------------------------------------------------
            
            // Prepare execution of program to test.
            if(pProcStr == NULL ) 
              cout << "\t\tPerform serial test (Test "<< testCnt <<"): ";
            else
              cout << "\t\tPerform test using " << pch <<" processor(s) (Test "<< testCnt <<"): ";

            cout<<flush;
            testCnt ++ ;                                                  
            string commandLine("");

                                   
            // If parallel tests is intendet, then set
            // mpi-prefix.            
            if(pProcStr != NULL ) {
              commandLine.append(pMpiRunner);
              commandLine.append(" -np ");
              commandLine.append(pch);
              commandLine.append(" ");
            } // END IF 
            
              commandLine.append( pExecPath);
              commandLine.append(" ");
              commandLine.append(pTempConfig);
              commandLine.append(" > /tmp/output.tmp 2>&1  \n");
            
            // Execute program            
            execResult = system(commandLine.c_str());
                                     
            this->readFileToBuffer("/tmp/output.tmp",programOutput);                        
                       
            unlink ("/tmp/output.tmp");;
            int foundPattern = checkOutputForPhrase(programOutput,errorPatterns);
            
            if(execResult == 0 && foundPattern < 0 )
              cout << " done. " <<endl;
            else
            {
              string msg = "Execution of program failed. ";
              msg.append("\nFound :");
              msg.append(errorPatterns[foundPattern]);
              msg.append("\nOutput follows:");
                     
              pCurrentTest->addResultMsg(msg);              
              pCurrentTest->setResultCode(REPORTER_RES_CODE_ERR_EXECUTION);              
              pCurrentTest->addResultMsg(programOutput);              
              cout << " failed " << endl;
              skipStep = true;
            }

              // Check if output-contains any error-messages
            cout << "\t\t\t Aufruf: " << commandLine << "\n";
            commandLine.erase();                                                        
                                                    

            //------------------------------------------------
            // Verify output (optional)
            //
            //------------------------------------------------            
             // Depending on the execution-result additional tests are
             // performed
            // if verifier program was given it will be executed now and
            // and the following test are proceded on its output.
            // if not, we check the normal program-output.

			// If no special proc-number is set for verfifier-run, run the the program
			// with the number of procs used before.

			/// \todo What should happen if less proc-numbers are defined than for the normal program-call?
			if(pvch==NULL)
				pvch = pch;
                        
            if(!pVerifierPath.empty() && !skipStep)
            {
                    // Increase execution level, so we can detetect which program
                    // causes the test to aboard.
                    execLevel++ ;
                    
                    cout << "\t\tRun verifier:" ;
                    if(pProcStr != NULL ) {
                      commandLine.append(pMpiRunner);
                      commandLine.append(" -np ");					 
                      commandLine.append(pvch);
                      commandLine.append(" ");
                    } 

                    commandLine.append( pVerifierPath);
                    commandLine.append(" ");
                    commandLine.append(verififierConfig);

                    char pTempConfig[] = " > /tmp/output.tmp 2>&1  ";
                    mkstemp(pTempConfig);                    
                    commandLine.append(pTempConfig);
                                  

                    execResult = system(commandLine.c_str());
                    
                                  
                    // Load output from last program-execution ( e.g the testet program or the verifier ) 
                    this->readFileToBuffer("/tmp/output.tmp" , programOutput );

                     unlink("/tmp/output.tmp");
                    // If any error occured during program-execution, aboard this test
                    // an log message.
                    // Check whether for error-messages in output

                    bool foundError = false;
                    // Check if output-contains any error-messages        
                    int foundPattern =  checkOutputForPhrase(programOutput,errorPatterns);
            
                    if(execResult != 0 || foundPattern >= 0 )
                    {
                      string msg = "Execution of verifier failed. ";
                      msg.append("\nFound :");
                      msg.append(errorPatterns[foundPattern]);
                      msg.append("\nOutput follows:");                    
                      pCurrentTest->addResultMsg(msg);
                      pCurrentTest->setResultCode(REPORTER_RES_CODE_ERR_EXECUTION);
                      pCurrentTest->addResultMsg(programOutput);
                      skipStep= true;
                      cout<< " failed. \n";
                    }
                    else
                    {
                      cout<< " done. \n";                      
                    } // END IF 
                    
                    cout << "\t\t\t Aufruf: " << commandLine << "\n";
                    commandLine.erase();
            } // END IF : Verifier

           
            //------------------------------------------------
            // Check output against defined rules
            //
            //------------------------------------------------                                               
            if(!skipStep)
            {
              cout << "\t\tCheck output." ;
                // Iterate over all values an check for
                // the right values in progam output
                TiXmlElement * config = bElement->FirstChildElement("values");
                for( TiXmlElement* sett = config->FirstChildElement("value"); sett; sett = sett->NextSiblingElement("value") )
                {

                  string result,pMode,pParams,modeStr;

                  // found value
                  if(this->getValueByRegex(sett->Attribute("regExp"), programOutput,result)) {

                    // Perform test on found value like defined in config-file
                    size_t matchCnt = 3 ;
                    regmatch_t matchPtr[3];
                    // Do the query
                    modeStr = string(sett->Attribute("mode"));
                    regexec (&modeCompiledExp,modeStr.c_str(), matchCnt, &matchPtr[0], 0);


                    pMode =  modeStr.substr(matchPtr[1].rm_so,matchPtr[1].rm_eo-matchPtr[1].rm_so) ;
                    // If found, parse mode-params
                    if(matchPtr[2].rm_so > 0 ) {
                      pParams = modeStr.substr(matchPtr[2].rm_so,matchPtr[2].rm_eo-matchPtr[2].rm_so) ;
                      pParams.erase(0,1);
                      pParams.erase(pParams.size()-1,1);
                    } // END IF


                    // Convert values
                    double a = atof(result.c_str());
                    double b = atof(sett->GetText());

                    // this string is appended to all messages..
                    string extender("( Output:");
                    extender.append(result);
                    extender.append(" | Expected:");
                    extender.append(sett->GetText());
                    extender.append(" )");

                    string msg("");


                    char convBuffer[255];
                    bool failedTest= false;
                    // Check what test to perfrom on the found value

                    // Value has to be (nearly) equal.
                    if(pMode.compare("equals") ==0) {

                      double max = b ;
                      if(fabs(a) > fabs(b) ) {
                        max= fabs(a) ;
                      } else {
                        max = fabs(b);
                      }// END IF

                      double tolerance = atof(pParams.c_str());
					  
					  double deviation=0;
					  if(max!=0)
							deviation = (fabs(a-b) / max);
							  
                      // check wether the value is in tolerance-range
                      if(deviation <= tolerance) {
                        msg.append("Test passed. Program-output equals expected output with tolerance.");
                        msg.append(extender);
                        
                        pCurrentTest->addCheckValues(sett->Attribute("label"),msg );
                      } else {
                        msg.append("Test failed: program-output ");
                        if(tolerance == 0 ) {
                          msg.append(" doesn't equal expected output ");
                        } else {
                          msg.append(" is not within the tolerance-range: ");
                          sprintf(convBuffer," (tolerance %f )",tolerance);
                          msg.append(convBuffer);

                        }

                        msg.append(extender);                        
                        pCurrentTest->addCheckValues(sett->Attribute("label"),msg);                    
                        failedTest = true;
                      } // END IF


                    } else {
                      // Value has to be below a certain value
                      if(pMode.compare("less") ==0) {

                        if ( a < b ) {
                          msg.append("Test passed (program-output is less than expected value)");
                          msg.append(extender);
                          
                          pCurrentTest->addCheckValues(sett->Attribute("label"),msg);                      
                        } else {

                          msg.append("Test failed: program-output is bigger or equal than expected output.");
                          msg.append(extender);
                          
                          pCurrentTest->addCheckValues(sett->Attribute("label"),msg);
                          failedTest = true;
                          
                        } // END IF


                      } else
                        // Value has to be below or equal to a certain value
                        if (pMode.compare("eless") ==0) {

                          if ( a <= b ) {
                            msg.append("Test passed (program-output is less or equal than expected value)");
                            msg.append(extender);
                            
                            pCurrentTest->addCheckValues(sett->Attribute("label"),msg);

                          } else {
                            msg.append("Test failed: program-output is bigger than expected output.");
                            msg.append(extender);
                            
                            pCurrentTest->addCheckValues(sett->Attribute("label"),msg);
                       
                            failedTest = true;
                            
                          } // END IF

                        } else
                          // Value has to be above a certain value
                          if (pMode.compare("more") ==0) {

                            if ( a > b ) {

                              msg.append("Test passed (program-output is bigger than expected value)");
                              msg.append(extender);

                              pCurrentTest->addCheckValues(sett->Attribute("label"),msg);

                            } else {
                              msg.append("Test failed: program-output is smaller or equal than expected output.");
                              msg.append(extender);
                              
                              pCurrentTest->addCheckValues(
                              sett->Attribute("label"),msg);
                              
                              failedTest = true;
                            } // END IF

                          } else
                            // Value has to be above or equal to a certain value
                            if (pMode.compare("emore") ==0) {

                              if ( a >= b ) {

                                msg.append("Test passed (program-output is bigger or equal than expected value)");
                                msg.append(extender);

                                pCurrentTest->addCheckValues(
                                sett->Attribute("label"),msg);

                              } else {

                                msg.append("Test failed: program-output is smaller than expected output.");
                                msg.append(extender);
                                
                                pCurrentTest->addCheckValues(
                                sett->Attribute("label"),msg);
                                failedTest = true;
                              } // END IF

                            }
                            // No valid Mode.
                            else {
                              pCurrentTest->setResultCode(REPORTER_RES_CODE_ERR_EXECUTION);
                              pCurrentTest->addResultMsg("Found invalid mode-code in configuration-file.");
                              pCurrentTest->addResultMsg(sett->Attribute("mode"));
                            } // END IF
                    } // END IF

                    if(failedTest) {                      
                      pCurrentTest->setResultCode(REPORTER_RES_CODE_ERR_TESTS);
                      pCurrentTest->addResultMsg("Value failed test:");
                      pCurrentTest->addResultMsg(sett->Attribute("label"));                      
                    } // END IF
                  }
                  
                  // value not found - generate message ...
                  else {
                    pCurrentTest->addCheckValues(
                        sett->Attribute("label"),"No value found in program-output"
                                                );

                    pCurrentTest->setResultCode(REPORTER_RES_CODE_ERR_TESTS);
                    pCurrentTest->addResultMsg("The following label doesn't fit the defined requirements ");
                    pCurrentTest->addResultMsg(sett->Attribute("label"));
                    
                  } // END IF
                } // END FOR : apply settings
                             
                cout<<" done \n";
            } // END IF : execRes != 0 

                                                                                                                                                  
            // Parse proc-string to get number of processors for next step.
            pch = strtok (NULL, " ;,");			
			pvch = strtok (NULL, " ;,");
			
            // Remove temporary output
            unlink ("./output.tmp");
            cout<<"\n";
          } // END WHILE
                                                                                        
          // remove the temporary conf-file for this testset
          pTempConfig.erase();
        } // END FOR : inner loop
                                        
      }       // END FOR      : outer loop
                
      return this->lastReport;
    } // END OF FUNCTION
                        
    /**
    *  
    *  This function searches for an value defined by the given
    * regular expression in the given text.
    * The regex is expected to contain at least one sub-expression
    *            
    * If no matching value is found, false is returned. Otherwise true.
    * If the regexp contains errors, an error-message is displayed.
    */
    bool TestSuite::getValueByRegex(const char * regExp,string &searchStr,string &pReturn) {
                                        

      pReturn.erase();      
      size_t matchCnt = 2 ;
      regmatch_t matchPtr[2];
                        
      regex_t compiledExp ;
      // Compile the regex
      regcomp (&compiledExp, regExp, REG_EXTENDED);
                                                                                        
      // Check for sub-expression
      if(compiledExp.re_nsub < 1 ) {
        cerr << "getValueByRegex: " << endl;
        cerr << "the regular expression has to contain at least one sub-expression: " << endl;
        cerr << regExp << endl;
        return false;
      } // END IF
                                                                                                                                                                

      // Do the query
       int resCode = regexec (&compiledExp,searchStr.c_str(), matchCnt, &matchPtr[0], 0);
                        
      // Handle errors
      if(resCode == REG_NOMATCH ) {
        return false;
      } else if(resCode != 0 ) {
        char buffer[1024];
        cerr << "getValueByRegex: " << endl;
        cerr << regExp << endl;
        regerror (resCode, &compiledExp, &buffer[0], sizeof(buffer));
        cerr << buffer << endl;
        return false ;
      }
                                                                              
      pReturn+=  searchStr.substr(matchPtr[1].rm_so,matchPtr[1].rm_eo-matchPtr[1].rm_so);
                                
      return true;
        
    } // END OF FUNCTION
 
        /**
     * \short This functions reads the whole given file into a string-object.
     *
     * \param pFileToOpen Location of the file that should be read.
     * \param target A pointer to an object of class string.
     *                The whole file-content is appended to this string.
     *
     * \return boolean. True if file could be readed sucessfull. False if not.
         */
    bool TestSuite::readFileToBuffer(const string pFileToOpen  ,string &rString ) {
          
      std::ifstream listfile (pFileToOpen.c_str(), std::ios_base::in);
      std::string line;
      rString.clear();
      rString.erase();             
      while(!listfile.eof() && !listfile.fail()) {
          std::getline(listfile, line);                    
          rString.append(line);
          rString.append("\n");           
      } // END OF WHILE      
      listfile.close();                 
      return true;
    } // END OF FUNCTION

    /**
     *
     *
     */
    void TestSuite::buildConfigFromTpl(TiXmlElement* bElement , const string pSourceFile,string &generatedFile)
    {

                                                  // Open template for configuration-file
      string pTplConfigStr("");

      if(!this->readFileToBuffer(pSourceFile,pTplConfigStr)) {
        cerr << " Error reading configuration-template '"<<pSourceFile<<"' "<<endl;
        return ;
      } // END IF
            
      FILE  *writer = NULL ;
      char pTempConfig[] = "/tmp/fileXXXXXX";
      mkstemp(pTempConfig);
      generatedFile = pTempConfig;
          
      writer = fopen(pTempConfig,"w");
      if(writer)
      {
                                                      // replace all knwon
                                                      // settings and write the new line into tmp-config
        string newString(pTplConfigStr.c_str());

            // Apply all settings to config-template..
        map<string,string> settings;
            
        TiXmlElement * config = bElement->FirstChildElement("config");
        for( TiXmlElement* sett = config->FirstChildElement("setting"); sett; sett = sett->NextSiblingElement("setting") ){
          int pos = -1;
          string pSettLabel(sett->Attribute("label"));
          string pSettValue(sett->GetText());
          settings.insert(map<string,string>::value_type(pSettLabel,pSettValue));
                                                                            // replace all occurences in current line
          while ((pos = newString.find(pSettLabel)) && pos !=-1) {
            newString.erase(pos, pSettLabel.length());
            newString.insert(pos, pSettValue);
          } // END FOR
        } // END FOR : apply settings

        fputs(newString.c_str(),writer);
        fclose(writer);

      } // END IF : Prepate configuration file                   
    } // END OF FUNCTION