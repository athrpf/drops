//#include <fstream>
//#include <stdlib.h>
//#include <string.h>
//#include <sstream>
//#include <math.h>
//#include<boost/lexical_cast.hpp>


  //void CreateRawData(int, double, int, double);
  
  
  void CreateRawData(int tsteps, double delta_t, int xsteps, double incx_spline){  
    
    //Raw-data-file-output:
    double x_temp; double t_temp;
    for(int i=0; i<tsteps; i++){
       std::string number = boost::lexical_cast<std::string>(i);
       std::string h_filename = "./hrawdata/h" + number + ".txt";
       std::string q_filename = "./qrawdata/q" + number + ".txt";
       std::ofstream hfile;
       std::ofstream qfile;
       hfile.open(h_filename.c_str());
       qfile.open(q_filename.c_str());
       t_temp = static_cast<double>(i)*delta_t;
       for(int k=0; k<xsteps; k++){
           x_temp = static_cast<double>(k)*incx_spline;
           hfile << sin(x_temp+t_temp) + 2. << std::endl;
           qfile << (-1.)*sin(x_temp+t_temp) + 2. << std::endl;
       }
       hfile.close();
       qfile.close();
   }
  }
  
  //int main(){
  
  //CreateRawData(10, 1.e-5, 10, 0.000209);
  
  //return 0;
  
  //}
  
