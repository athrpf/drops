#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include<boost/lexical_cast.hpp>


 
int main() {

    int tsteps = 10;
    int xsteps = 100;
    double delta_t = 0.1;
    double incx_spline= 0.5;
    double x_temp; 
    double t_temp;

    //Raw-data-file-output:  
    
    for(int i=0; i<tsteps; i++){
        std::string number = boost::lexical_cast<std::string>(i);
        std::string h_filename = "h" + number + ".txt";
        std::string q_filename = "q" + number + ".txt";
        std::ofstream hfile;
        std::ofstream qfile;
        hfile.open(h_filename.c_str());
        qfile.open(q_filename.c_str());
        t_temp = static_cast<double>(i)*delta_t;
        for(int k=0; k<xsteps; k++){
            x_temp = static_cast<double>(k)*incx_spline;
            hfile << k << " i is "<< number << std::endl;
            qfile << k << " i is "<< number << std::endl;
        }
        hfile.close();
        qfile.close();
    }
return 0;
}
