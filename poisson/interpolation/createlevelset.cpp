#include "periodicdata.hpp"
#include <fstream>
#include <sstream>
#include "createlevelset.hpp"
extern PeriodicData pd;

 double* createlevelLiang() {

  double* levelLiang; //[4141]

  // allocate
  int NumCoordsLiang = (pd.NX-1)*pd.NY;
  levelLiang = new double[NumCoordsLiang];

  //Please mind: The following relative path to the source-data "LevelSetLiangData.txt" referrs to the
  //directory where fullDNS.cpp is located.
  std::string level_filename = "./liangdata/Interpolation/DataForPoissonCoeff/LevelSetLiangData.txt";
  std::ifstream levelfile;
  levelfile.open(level_filename.c_str(), std::fstream::in);
  double curr_level;
  for (int k=0; k<NumCoordsLiang; k++) {
    levelfile >> curr_level;
    levelLiang[k] = curr_level;
  }




return levelLiang;


}
