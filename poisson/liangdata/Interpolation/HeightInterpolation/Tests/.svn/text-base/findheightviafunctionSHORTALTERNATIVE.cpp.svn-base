#include "../../PERIODICADDEDphysicaldata.hpp"
#include "../createlevelLiang.hpp"
#include "../DiscreteLevelSetToDiscreteHeight.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
static PhysicalData pd;

int main() {

  //Create an array of Liang's level-set-values, which is converted to an array of discrete height-values by use of the 
  //function "DiscreteLevelSetToDiscreteHeight":
  createlevelLiang();
  DiscreteLevelSetToDiscreteHeight(createlevelLiang());
    
  for (int k=0; k<pd.NX; k++) {
    std::cout<< DiscreteLevelSetToDiscreteHeight(createlevelLiang())[k] << std::endl;
  }

return 0;

}




 
 
 
 
