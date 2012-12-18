#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include<boost/lexical_cast.hpp>
#include "../liangdata/Interpolation/boxinterpLiangData.cpp"





int main() {


double x=0.1;
double t=3.e-4;
double delta_t=1.e-5;
    
    int tsteps=10; int xsteps=10; double incx_spline=0.209; 
    
    gsl_spline ** hspline, ** qspline; gsl_interp_accel * acc;
    
    hspline = new gsl_spline*[tsteps];
    qspline = new gsl_spline*[tsteps];
    acc = gsl_interp_accel_alloc();
    
    double * xd = new double[xsteps];
    for(int i=0; i<xsteps; i++){
            xd[i]=static_cast <double>(i)*incx_spline;
    }
    
    double * h_raw = new double[xsteps];
    double * q_raw = new double[xsteps];
     
    //Raw-data-file-output:
    double x_temp; double t_temp;
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
           hfile << sin(x_temp+t_temp) + 2. << std::endl;
           qfile << (-1.)*sin(x_temp+t_temp) + 2. << std::endl;
       }
       hfile.close();
       qfile.close();
   }
   
    //Raw-data-file-input, create splines:
    double x_temp2; double t_temp2;
    double curr_h; double curr_q;
    for(int i=0; i<tsteps; i++){
       std::string number = boost::lexical_cast<std::string>(i);
       std::string h_filename = "h" + number + ".txt";
       std::string q_filename = "q" + number + ".txt";
       std::ifstream hfile;
       std::ifstream qfile;
       hfile.open(h_filename.c_str(), std::fstream::in);
       qfile.open(q_filename.c_str(), std::fstream::in);
       t_temp2 = static_cast<double>(i)*delta_t;
       for(int k=0; k<xsteps; k++){
           
           hfile >> curr_h;
           qfile >> curr_q;
           
           h_raw[k] = curr_h;
           q_raw[k] = curr_q;
           
           x_temp2 = static_cast<double>(k)*incx_spline;
           std::cout <<"===============================" << std::endl;
           std::cout << curr_h - (sin(x_temp2+t_temp2) + 2.) << std::endl;
           std::cout << curr_q - ((-1.)*sin(x_temp2+t_temp2) + 2.) << std::endl;
           std::cout <<"===============================" << std::endl;
       }
       
       hfile.close();
       qfile.close();
       
       hspline[i] = gsl_spline_alloc(gsl_interp_cspline, xsteps);
       gsl_spline_init(hspline[i], xd, h_raw, xsteps);
        
       qspline[i] = gsl_spline_alloc(gsl_interp_cspline, xsteps);
       gsl_spline_init(qspline[i], xd, q_raw, xsteps);
       
   }
          
    delete[] h_raw;
    delete[] q_raw;
    
    
    
    gsl_spline ** hspline2, ** qspline2; gsl_interp_accel * acc2;
    
    hspline2 = new gsl_spline*[tsteps];
    qspline2 = new gsl_spline*[tsteps];
    acc2 = gsl_interp_accel_alloc();
    
    
    double * h_raw2 = new double[xsteps];
    double * q_raw2 = new double[xsteps];
    
    //double x_temp2; double t_temp2;
    for(int i=0; i<tsteps; i++){
        t_temp2 = static_cast<double>(i)*delta_t;
        for(int k=0; k<xsteps; k++){
            x_temp2 = static_cast<double>(k)*incx_spline;
            h_raw2[k] = sin(x_temp2+t_temp2) + 2.;
            q_raw2[k] = (-1.)*sin(x_temp2+t_temp2) + 2.;
        }
        hspline2[i] = gsl_spline_alloc(gsl_interp_cspline, xsteps);
        gsl_spline_init(hspline2[i], xd, h_raw2, xsteps);
        
        qspline2[i] = gsl_spline_alloc(gsl_interp_cspline, xsteps);
        gsl_spline_init(qspline2[i], xd, q_raw2, xsteps);
     }
      
     delete[] h_raw2;
     delete[] q_raw2;
   
  
    //int nt = boxnumber(t, delta_t, 0.);
    
   // double h = gsl_spline_eval(hspline[nt], x, acc);
    //double q = gsl_spline_eval(qspline[nt], x, acc);
    
   // double h2 = gsl_spline_eval(hspline2[nt], x, acc2);
    //double q2 = gsl_spline_eval(qspline2[nt], x, acc2);
    
   // std::cout << "h-h2 = " << h-h2 << std::endl;
   // std::cout << "q-q2 = " << q-q2 << std::endl;
    
    return 0;
}
