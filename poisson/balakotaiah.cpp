/****************************************
 *  Balakotaiah:                        \n*
 *   - instationary setup             \n*
 *   - constant diffusion             \n*
 *   - convection                     \n*
 *   - no reaction                    \n*
 ****************************************/
#include "poisson/poissonCoeff.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

extern DROPS::ParamCL P;

namespace Balakotaiah
{
  int getn(std::string dfile)
  {
    int n=0;
    double tmp;
    std::ifstream fd;
    fd.open(dfile.c_str());
    if(fd.fail())
      std::cout << dfile << " not available" << std::endl;

    while(!fd.eof()){
      if(fd.good()){
	fd >> tmp;
	n++;
      }
    }

    fd.close();
    return n;
  }

  int countTimesteps(const int starttime, const int stoptime, const int timeinc)
  {
    int tmp = 0;
    for(int i=starttime; i<=stoptime; i+=timeinc)
      tmp++;

    return tmp;
  }

  void getData(double *d, std::string dfile, std::string time)
  {
    std::ifstream fd;
    std::string file;
    file = dfile + time + ".dat";
    fd.open(file.c_str());
    int i=0;

    while(!fd.eof()){
      if(fd.good()){
	fd >> d[i];
	i++;
      }
    }

    fd.close();
  }

  static bool     init=false;

  static int      n;

  static gsl_spline ** hspline,
		   ** qspline;

  static gsl_interp_accel * acc;

  static std::string hfile,
		   qfile;

  static int  starttime,
		   stoptime,
		   timeinc;

  void initialize()
  {
    int id = omp_get_thread_num();
    //std::cout << "This is thread " << id << std::endl;
    if (id==0) { // only the master thread does this
      std::string file = P.get<std::string>("andreasOpt.DataDir");
      hfile = file + "/longh";
      qfile = file + "/longq";
      starttime = P.get<int>("andreasOpt.starttime");
      stoptime = P.get<int>("andreasOpt.stoptime");
      timeinc = P.get<int>("Time.StepSize");
      std::stringstream out;
      out << starttime;
      file=hfile+out.str()+".dat";
      n=getn(file);

      int tsteps = countTimesteps(starttime,stoptime,timeinc);
      acc = gsl_interp_accel_alloc();
      hspline = new gsl_spline*[tsteps];
      qspline = new gsl_spline*[tsteps];

      double * h,
	* q,
	* xa;

      h = new double[n];
      q = new double[n];
      xa = new double[n];

      for(int j=0; j<n; j++)
	xa[j] = j;

      std::cout << "Getting Data... ";
      for(int i=0; i<tsteps; i++){
	std::stringstream time;
	time << starttime + i*timeinc;
	getData(h, hfile, time.str());
	hspline[i] = gsl_spline_alloc(gsl_interp_cspline, n);
	gsl_spline_init(hspline[i], xa, h, n);
	getData(q, qfile, time.str());
	qspline[i] = gsl_spline_alloc(gsl_interp_cspline, n);
	gsl_spline_init(qspline[i], xa, q, n);
      }
      std::cout << " done" << std::endl;


      init=true;

      delete[] h;
      delete[] q;
      delete[] xa;
#pragma omp barrier
    }
  }

  double get_h(int t, double x)
  {
    if(!init)
      initialize();
    int nt = (int)((t+0.1)/timeinc);
    return gsl_spline_eval(hspline[nt], x, acc);
  }

  void uvh(double* res, int t, double x, double y)
  {
    double          hres[2],
      qres[2];
    if(!init)
      initialize();

    int nt = (int)((t+0.1)/timeinc);
    hres[0] = gsl_spline_eval(hspline[nt], x, acc);
    hres[1] = gsl_spline_eval_deriv(hspline[nt], x, acc);
    qres[0] = gsl_spline_eval(qspline[nt], x, acc);
    qres[1] = gsl_spline_eval_deriv(qspline[nt], x, acc);

    double h1 = hres[0];
    double h2 = h1*h1;
    double h3 = h2*h1;
    double h4 = h2*h2;
    double y2 = y*y;
    double y3 = y2*y;
    res[0]=3*qres[0]/h3*(hres[0]*y-y2/2);
    res[1]=-(1.5*qres[0]/h4*hres[1]-0.5*qres[1]/h3)*y3
      -(1.5*qres[1]/h2-3*qres[0]/h3*hres[1])*y2;
  }
  /*****************************************************************************************************************/

  double Interface( const DROPS::Point3DCL& p, double t)
  {
    //std::cout << "t: " << t  << "x: " << p[0] << "y: " << p[1] << std::endl;
    return get_h((int)(t+0.1), p[0]);
  }

  /// \brief Reaction: no reaction
  double Reaction(const DROPS::Point3DCL&, double){
    return 0.0;
  }

  /// \brief Diffusion
  double Diffusion(const DROPS::Point3DCL&, double){
    return 1.e-6;
  }

  /// \brief Initial value
  double InitialValue( const DROPS::Point3DCL& , double){
    return 1e-5;
  }

  DROPS::Point3DCL Flowfield(const DROPS::Point3DCL& p, double t){
    double res[2];
    //std::cout << "t: " << t << std::endl;
    uvh(res,(int)(t+0.1),p[0],p[1]);
    DROPS::Point3DCL v(0.);
    v[0] = res[0];
    v[1] = res[1];
    v[2] = 0.;
    return v;
  }

  /// \brief Solution
  double Solution( const DROPS::Point3DCL& p, double)
  {
    if(p[0] < 15)
      return 1e-4/15*p[0];
    else
      return 1e-4;
  }

  /// \brief Right-hand side
  double Source(const DROPS::Point3DCL&, double){
    return 0.;
  }
  static DROPS::RegisterScalarFunction regscaq("Balakotaiah_Reaction",     DROPS::instat_scalar_fun_ptr(Reaction)     );
  static DROPS::RegisterScalarFunction regscaf("Balakotaiah_Source",       DROPS::instat_scalar_fun_ptr(Source)       );
  static DROPS::RegisterScalarFunction regscaa("Balakotaiah_Diffusion",    DROPS::instat_scalar_fun_ptr(Diffusion)    );
  static DROPS::RegisterScalarFunction regscaint("Balakotaiah_Interface",  DROPS::instat_scalar_fun_ptr(Interface)    );
  static DROPS::RegisterScalarFunction regscas("Balakotaiah_Solution",     DROPS::instat_scalar_fun_ptr(Solution)     );
  static DROPS::RegisterVectorFunction regscav("Balakotaiah_Velocity",     DROPS::instat_vector_fun_ptr(Flowfield)    );
  static DROPS::RegisterScalarFunction regscai("Balakotaiah_InitialValue", DROPS::instat_scalar_fun_ptr(InitialValue) );
}//end of namespace
