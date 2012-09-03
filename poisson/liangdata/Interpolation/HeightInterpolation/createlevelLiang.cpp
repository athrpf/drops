#include "../PERIODICADDEDphysicaldata.hpp"
#include <fstream>
#include <sstream>
static PhysicalData pd;

 double* createlevelLiang() {

  double* levelLiang; //[4141]

  // allocate
  int NumCoordsLiang = (pd.NX-1)*pd.NY;
  levelLiang = new double[NumCoordsLiang];


  std::string level_filename = "LevelSetLiangData.txt";
  std::ifstream levelfile;
  levelfile.open(level_filename.c_str(), std::fstream::in);
  double curr_level;
  for (int k=0; k<NumCoordsLiang; k++) {
    levelfile >> curr_level;
    levelLiang[k] = curr_level;
  }
  

  

return levelLiang;


}




 
 
 
 
