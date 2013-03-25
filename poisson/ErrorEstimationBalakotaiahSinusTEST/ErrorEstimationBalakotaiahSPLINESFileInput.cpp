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
#include <iomanip>


class HQUV{

private:

double Re; double We; double cotbeta; double x; double t; 

double inct_spline; double incx_spline; int xsteps;//time-increment splines, spacial 
        //increment of splines, #data-points in x-direction (depends on the size of the 
        //files in the raw-data-folders ./hrawdata and ./qrawdata)

double h; double h_x; double h_xx; double h_xxx; double h_xxxx;

double q; double q_x; double q_xx; double q_xxx; double q_xxxx;

public:

HQUV(double X, double T, int XSTEPS, double INCX, double INCT, double RE, double WE, double COTBETA); //Konstruktor
double get_Re(void) { return Re;}
double get_We(void) { return We;}
double get_cotbeta(void) { return cotbeta;}
double get_x(void) { return x;}
double get_t(void) { return t;} 


double get_inct_spline(void) { return inct_spline;}
double get_incx_spline(void) { return incx_spline;}
int get_xsteps(void) { return xsteps;}

double h0(void) {return h;}
double h1(void) {return h_x;}
double h2(void) {return h_xx;}
double h3(void) {return h_xxx;}
double h4(void) {return h_xxxx;}

double q0(void) {return q;}
double q1(void) {return q_x;}
double q2(void) {return q_xx;}
double q3(void) {return q_xxx;}
double q4(void) {return q_xxxx;}

double U1_0(void); double U1_1(void); double U1_2(void); double U1_3(void);
double U2_0(void); double U2_1(void); double U2_2(void); double U2_3(void);

double V2_0(void); double V2_1(void); double V2_2(void); double V2_3(void);
double V3_0(void); double V3_1(void); double V3_2(void); double V3_3(void);
};
///////////////////////Konstruktor definieren:
HQUV::HQUV(double X, double T, int XSTEPS, double INCX, double INCT, double RE, double WE, double COTBETA){
    
    x=X; t=T;  
    
    xsteps=XSTEPS; incx_spline=INCX; inct_spline=INCT; 
    
    Re=RE; We=WE; cotbeta=COTBETA;
    
    gsl_spline * hspline, * qspline; gsl_interp_accel * acc;
    
    hspline = new gsl_spline;
    qspline = new gsl_spline;
    
    acc = gsl_interp_accel_alloc();
    
    double * xd = new double[xsteps];
    for(int i=0; i<xsteps; i++){
            xd[i]=static_cast <double>(i)*incx_spline;
    }
    
    double * h_raw = new double[xsteps];
    double * q_raw = new double[xsteps];
  
       
    //Raw-data-file-input, create splines:
    
    int nt = boxnumber(t, inct_spline, 0.);
    double curr_h1; double curr_q1; double curr_h2; double curr_q2;
 
       std::string number1 = boost::lexical_cast<std::string>(nt);
       std::string number2 = boost::lexical_cast<std::string>(nt+1);
       std::string h1_filename = "./hrawdata/h" + number1 + ".txt";
       std::string h2_filename = "./hrawdata/h" + number2 + ".txt";
       std::string q1_filename = "./qrawdata/q" + number1 + ".txt";
       std::string q2_filename = "./qrawdata/q" + number2 + ".txt";
       std::ifstream h1file;
       std::ifstream h2file;
       std::ifstream q1file;
       std::ifstream q2file;
       h1file.open(h1_filename.c_str(), std::fstream::in);
       h2file.open(h2_filename.c_str(), std::fstream::in);
       q1file.open(q1_filename.c_str(), std::fstream::in);
       q2file.open(q2_filename.c_str(), std::fstream::in);
       
       for(int k=0; k<xsteps; k++){
           
           h1file >> curr_h1;
           h2file >> curr_h2;
           q1file >> curr_q1;
           q2file >> curr_q2;
           
           //linear interpolation of raw-data:
           h_raw[k] = ((curr_h2 - curr_h1)/inct_spline)*(t - nt*inct_spline) + curr_h1;
           q_raw[k] = ((curr_q2 - curr_q1)/inct_spline)*(t - nt*inct_spline) + curr_q1;
         
       }
       
       h1file.close();
       h2file.close();
       q1file.close();
       q2file.close();
       
       hspline = gsl_spline_alloc(gsl_interp_cspline, xsteps);
       gsl_spline_init(hspline, xd, h_raw, xsteps);
        
       qspline = gsl_spline_alloc(gsl_interp_cspline, xsteps);
       gsl_spline_init(qspline, xd, q_raw, xsteps);
          
    delete[] h_raw;
    delete[] q_raw;
    delete[] xd;
  
   
   
    double dx = incx_spline/(6.);
    
    h = gsl_spline_eval(hspline, x, acc);
    h_x = gsl_spline_eval_deriv(hspline, x, acc);
    h_xx = gsl_spline_eval_deriv2(hspline, x, acc);
    double h_xx_dx = gsl_spline_eval_deriv2(hspline, x + dx, acc);
    h_xxx = (h_xx_dx - h_xx)/dx; //gsl_spline_eval_deriv3(hspline, x, acc);
    h_xxxx = 0.;

    q = gsl_spline_eval(qspline, x, acc);
    q_x = gsl_spline_eval_deriv(qspline, x, acc);
    q_xx = gsl_spline_eval_deriv2(qspline, x, acc);
    double q_xx_dx = gsl_spline_eval_deriv2(qspline, x + dx, acc);
    q_xxx = (q_xx_dx - q_xx)/dx;//gsl_spline_eval_deriv3(qspline, x, acc);
    q_xxxx = 0.;    
    
}
/////////////////////////////////////////////
double HQUV::U1_0(void){
    return 3.*q/(h*h);
}
double HQUV::U2_0(void){
    return (-1.5)*q/(h*h*h);
}
double HQUV::U1_1(void){
    return 3.*q_x/(h*h)-6.*(q*h_x)/(h*h*h);
}
double HQUV::U2_1(void){
    return (-1.5)*q_x/(h*h*h) + 4.5*(q*h_x)/(h*h*h*h);
}
double HQUV::U1_2(void){
    double a = 3.*q_xx/(h*h) - 12.*(q_x*h_x)/(h*h*h);
    double b = 18.*(q*h_x*h_x)/(h*h*h*h) - 6.*(q*h_xx)/(h*h*h);
    return a + b; 
}
double HQUV::U2_2(void){
    double a = (-1.5)*q_xx/(h*h*h) + 9.*(q_x*h_x)/(h*h*h*h);
    double b = (-18.)*(q*h_x*h_x)/(h*h*h*h*h) + 4.5*(q*h_xx)/(h*h*h*h);
    return a + b; 
}
double HQUV::U1_3(void){
    double a = 3.*q_xxx/(h*h) - 18.*(q_xx*h_x)/(h*h*h) + 54.*(q_x*h_x*h_x)/(h*h*h*h);
    double b = (-18.)*(q_x*h_xx)/(h*h*h) - 72.*(q*h_x*h_x*h_x)/(h*h*h*h*h);
    double c = 54.*(q*h_x*h_xx)/(h*h*h*h) - 6.*(q*h_xxx)/(h*h*h); 
    return  a + b + c;
}
double HQUV::U2_3(void){
    double a = -1.5*q_xxx/(h*h*h) + 13.5*(q_xx*h_x)/(h*h*h*h) - 54.*(q_x*h_x*h_x)/(h*h*h*h*h);
    double b = 13.5*(q_x*h_xx)/(h*h*h*h) + 90.*(q*h_x*h_x*h_x)/(h*h*h*h*h*h);
    double c = 4.5*(q*h_xxx)/(h*h*h*h) - 54.*(q*h_x*h_xx)/(h*h*h*h*h);
    return a + b + c;
}
//////////////////////////////////////
double HQUV::V2_0(void){
    return 3.*(h_x*q)/(h*h*h) - 1.5*q_x/(h*h);
}
double HQUV::V3_0(void){
    return 0.5*q_x/(h*h*h) - 1.5*(h_x*q)/(h*h*h*h);
}
double HQUV::V2_1(void){
    return 6.*(q_x*h_x)/(h*h*h) - 9.*(q*h_x*h_x)/(h*h*h*h) + 3.*(q*h_xx)/(h*h*h) - 1.5*q_xx/(h*h);
}
double HQUV::V3_1(void){
    return 0.5*(q_xx)/(h*h*h) - 3.*(q_x*h_x)/(h*h*h*h) + 6.*(q*h_x*h_x)/(h*h*h*h*h) - 1.5*(q*h_xx)/(h*h*h*h);
}
double HQUV::V2_2(void){
    double a = 9.*(q_xx*h_x)/(h*h*h) - 27.*(q_x*h_x*h_x)/(h*h*h*h) + 9.*(q_x*h_xx)/(h*h*h);
    double b = 36.*(q*h_x*h_x*h_x)/(h*h*h*h*h) - 27.*(q*h_x*h_xx)/(h*h*h*h);
    double c = 3.*(q*h_xxx)/(h*h*h) - 1.5* q_xxx/(h*h);
    return a + b + c;
}
double HQUV::V3_2(void){
    double a = 0.5*(q_xxx)/(h*h*h) - 4.5*(q_xx*h_x)/(h*h*h*h) + 18.*(q_x*h_x*h_x)/(h*h*h*h*h);
    double b = -4.5*(q_x*h_xx)/(h*h*h*h) - 30.*(q*h_x*h_x*h_x)/(h*h*h*h*h*h);
    double c = 18.*(q*h_x*h_xx)/(h*h*h*h*h) -1.5*(q*h_xxx)/(h*h*h*h);
    return a + b + c;
}
double HQUV::V2_3(void){
    double a = 12.*(q_xxx*h_x)/(h*h*h) - 54.*(q_xx*h_x*h_x)/(h*h*h*h) + 18.*(q_xx*h_xx)/(h*h*h);
    double b = 144.*(q_x*h_x*h_x*h_x)/(h*h*h*h*h) - 108.*(q_x*h_x*h_xx)/(h*h*h*h);
    double c = 12.*(q_x*h_xxx)/(h*h*h) - 180.*(q*h_x*h_x*h_x*h_x)/(h*h*h*h*h*h);
    double d = 216.*(q*h_x*h_x*h_xx)/(h*h*h*h*h) - 27.*(q*h_xx*h_xx)/(h*h*h*h);
    double e = -36.*(q*h_x*h_xxx)/(h*h*h*h) + 3.*(q*h_xxxx)/(h*h*h) -1.5*(q_xxxx)/(h*h);
    return a + b + c + d + e;
}
double HQUV::V3_3(void){ 
    double a = 0.5*(q_xxxx)/(h*h*h) - 6.*(q_xxx*h_x)/(h*h*h*h) + 36.*(q_xx*h_x*h_x)/(h*h*h*h*h);
    double b = -9.*(q_xx*h_xx)/(h*h*h*h) - 120.*(q_x*h_x*h_x*h_x)/(h*h*h*h*h*h);
    double c = 72.*(q_x*h_x*h_xx)/(h*h*h*h*h) - 6.*(q_x*h_xxx)/(h*h*h*h);
    double d = 180.*(q*h_x*h_x*h_x*h_x)/(h*h*h*h*h*h*h) - 180.*(q*h_x*h_x*h_xx)/(h*h*h*h*h*h);
    double e = 18.*(q*h_xx*h_xx)/(h*h*h*h*h) + 24.*(q*h_x*h_xxx)/(h*h*h*h*h) - 1.5*(q*h_xxxx)/(h*h*h*h);
    return a + b + c + d + e;
}
/////////////////////////////////////////
double psurf_x(double h0, double h1, double h2, double h3, 
               double U1_0, double U1_1, double U1_2, 
               double U2_0, double U2_1, double U2_2, 
                            double V2_1, double V2_2,
                            double V3_1, double V3_2, double R, double W){
                       
    double H0_2 = h0*h0; double H0_3 = H0_2*h0; double H1_2 = h1*h1; 
    
    double a = (((-8.)*h1)/(1.+H1_2))*(U1_1 + 2.*U2_1*h0 + 2.*U2_0*h1 + V2_2*H0_2 + 2.*V2_1*h0*h1 + V3_2*H0_3 + 3.*V3_1*H0_2*h1);
    double b = (((-8.)*h2)/(1. + H1_2))*(U1_0 + 2.*U2_0*h0 + V2_1*H0_2 + V3_1*H0_3);
    double c = ((16.*H1_2*h2)/((1.+H1_2)*(1.+H1_2)))*(U1_0 + 2.*U2_0*h0 + V2_1*H0_2 + V3_1*H0_3);
    double d = ((-8.)*(1-H1_2)/(1.+H1_2))*(U1_2*h0 + U1_1*h1 + U2_2*H0_2 + 2.*U2_1*h0*h1);
    double e = (16.*(h1*h2)/(1.+H1_2))*(U1_1*h0 + U2_1*H0_2);
    double f = (16.*((1.-H1_2)*h1*h2)/((1.+H1_2)*(1.+H1_2)))*(U1_1*h0 + U2_1*H0_2);
    double g = R*W*((3.*h2*h2*h1)/((1.+H1_2)*(1.+H1_2)*sqrt(1.+H1_2))  -   h3/((1.+H1_2)*sqrt(1.+H1_2)));
    return a + b + c + d + e  + f + g;        
}
double diff_t(double A_t, double A_tplusdt, double dt){
    return (A_tplusdt - A_t)/dt;
}
double pbulk_x(double U1_0, double U1_1, 
               double U2_0, double U2_1, 
               double V2_0, double V2_1, double V2_2, double V2_3, double dtV2_0, double dtV2_1,
               double V3_0, double V3_1, double V3_2, double V3_3, double dtV3_0, double dtV3_1,
               double h0, double h1, double y, double R, double deltat, double cb){
               
    double H0_2 = h0*h0; double H0_3 = H0_2*h0; double H0_4 = H0_2*H0_2; double H0_5 = H0_4*h0; double H0_6 = H0_3*H0_3;         
               
    double a = y*8.*V2_1 + y*y*12.*V3_1 + y*y*y*(1./3.)*(4.*V2_3 - R*diff_t(V2_1,dtV2_1,deltat));
    double b = y*y*y*y*(V3_3 - 0.25*R*diff_t(V3_1,dtV3_1,deltat) - R*V2_0*V2_1 - 0.25*R*U1_0*V2_2 - 0.25*R*U1_1*V2_1);
    double c = (-0.2)*R*y*y*y*y*y*(5.*V2_0*V3_1 + 5.*V2_1*V3_0 + U1_1*V3_1 + U1_0*V3_2 + U2_1*V2_1 + U2_0*V2_2);
    double d = (-1./6.)*R*y*y*y*y*y*y*(U2_0*V3_2 + U2_1*V3_1 + 6.*V3_0*V3_1);
    double e = (1./3.)*R*diff_t(V2_1,dtV2_1,deltat)*H0_3 + (1./4.)*R*diff_t(V3_1,dtV3_1,deltat)*H0_4;              
    double f = 5.*R*V2_0*V3_0*H0_4*h1 - 4.*V2_2*H0_2*h1 - 4.*V3_2*H0_3*h1 + (1./6.)*R*U2_0*V3_2*H0_6;
    double g = 0.2*R*U1_0*V3_2*H0_5 + 0.2*R*U2_0*V2_2*H0_5;
    double h = 0.25*R*U1_0*V2_2*H0_4 - V3_3*H0_4 - (4./3.)*V2_3*H0_3;
    double i = R*U2_0*V3_1*H0_5*h1 + R*U1_0*V3_1*H0_4*h1 + R*U2_0*V2_1*H0_4*h1 + R*U1_0*V2_1*H0_3*h1;
    double j = 12.*cb*h1 - 12.*V3_1*H0_2 - 8.*V2_1*h0 - 8.*V2_0*h1 + 2.*R*V2_0*V2_0*H0_3*h1;
    double k = 3.*R*V3_0*V3_0*H0_5*h1 + 0.25*R*U1_1*V2_1*H0_4 + 0.2*R*U2_1*V2_1*H0_5;
    double l = 0.2*R*U1_1*V3_1*H0_5 + (1./6.)*R*U2_1*V3_1*H0_6 + R*V2_0*V2_1*H0_4;
    double m = R*diff_t(V3_0, dtV3_0, deltat)*H0_3*h1 + R*diff_t(V2_0, dtV2_0, deltat)*H0_2*h1;
    double n = R*V3_0*V3_1*H0_6 + R*V2_0*V3_1*H0_5 + R*V2_1*V3_0*H0_5 - 24.*V3_0*h0*h1;
    return a + b + c + d + e + f + g + h + i + j + k + l + m + n;          
}
double error(double U1_0, double U1_1, double U1_2, double U1_3, double dtU1_0,
              double U2_0, double U2_1, double U2_2, double U2_3, double dtU2_0,
              double V2_0, double V2_1, double V2_2, double V2_3, double dtV2_0, double dtV2_1,
              double V3_0, double V3_1, double V3_2, double V3_3, double dtV3_0, double dtV3_1,
              double h0, double h1, double h2, double h3, double y, double R, double W, double deltat, double cb){
              
     double a = R*(diff_t(U1_0, dtU1_0, deltat)*y + diff_t(U2_0, dtU2_0, deltat)*y*y);       
     double b = R*(U1_0*y + U2_0*y*y)*(U1_1*y + U2_1*y*y) + R*(V2_0*y*y + V3_0*y*y*y)*(U1_0 + 2*U2_0*y);        
     double c = -12. -4.*(U1_2*y + U2_2*y*y) -8.*U2_0;
     double d = pbulk_x(U1_0,U1_1,U2_0,U2_1,V2_0,V2_1,V2_2,V2_3,dtV2_0,dtV2_1,V3_0,V3_1,V3_2,V3_3,dtV3_0,dtV3_1,h0,h1,y,R,deltat, cb);      
     double e = psurf_x(h0, h1, h2, h3, U1_0, U1_1, U1_2, U2_0, U2_1, U2_2, V2_1, V2_2, V3_1,  V3_2, R, W);
     return a + b + c + d + e;        
}
                              

double ErrorBalakotaiahHQSplines(double X, double Y, double T, int XSTEPS, double INCX, double INCT, double RE, double WE, double COTBETA){

    double dt=1.e-5; //time-increment for temopral differentiation. Do not mix it up with the time-increment of the raw data.
    //HQUV::HQUV(double X, double T, int XSTEPS, double INCX, double INCT, double RE, double WE, double COTBETA)
    HQUV t(X, T, XSTEPS, INCX, INCT, RE, WE, COTBETA);
    HQUV t_dt(X, T + dt, XSTEPS, INCX, INCT, RE, WE, COTBETA);
    double retval; 
   
    retval = error(t.U1_0(), t.U1_1(), t.U1_2(), t.U1_3(), t_dt.U1_0(),
                  t.U2_0(), t.U2_1(), t.U2_2(), t.U2_3(), t_dt.U2_0(),
                  t.V2_0(), t.V2_1(), t.V2_2(), t.V2_3(), t_dt.V2_0(), t_dt.V2_1(),
                  t.V3_0(), t.V3_1(), t.V3_2(), t.V3_3(), t_dt.V3_0(), t_dt.V3_1(),
                  t.h0(), t.h1(), t.h2(), t.h3(), Y, RE, WE,dt, COTBETA);
                  
    return retval;

}


int main() {   

   
    double errorBala;
    std::string error_filename = "errorBalaspline.txt";
    std::ofstream errorfile;
    errorfile.precision(16);
    errorfile.open(error_filename.c_str());
    
    double newincx=0.01;
    double X;
    const double Y = 0.005;
    double T=4.e-3;
    int XSTEPS = 1000;
    double INCX=0.000209;
    double INCT=1.e-5;
    double RE=1.; 
    double WE=1.; 
    double COTBETA=0.;

    for(int k=0; k<20; k++){  

        X=newincx*static_cast<double>(k);
        
        errorBala=ErrorBalakotaiahHQSplines(X, Y, T, XSTEPS, INCX, INCT, RE, WE, COTBETA);

        std::cout << "For k= " << k << "the error in x-momentum-balance is: " << errorBala << std::endl; 
        errorfile << errorBala << std::endl;
    
    }        

errorfile.close();      

return 0;
}
