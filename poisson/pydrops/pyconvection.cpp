#include <boost/python.hpp>
#include "misc/params.h"

#include "pyconnect.h"
#include "convection_diffusion.cpp"

bool check_dimension(const PdeFunction& f, int Nx, int Ny, int Nz, int Nt)
{
  return f.get_dimensions(Nx, Ny, Nz, Nt);
}

bool check_dimensions(int Nx, int Ny, int Nz, int Nt, const PdeFunction& C0, const PdeFunction& b_in, const PdeFunction& b_interface, const PdeFunction& source, const PdeFunction& Dw)
{
  if (check_dimension(C0, Nx, Ny, Nz, 1) &&
      check_dimension(b_in, 1, Ny, Nz, Nt) &&
      check_dimension(b_interface, Nx, 1, Nz, Nt) &&
      check_dimension(source, Nx, Ny, Nz, Nt) &&
      check_dimension(Dw, Nx, Ny, Nz, Nt)) {
    return true;
  }
  return false;
}

bool py_convection_diffusion(PdeFunction& C0, PdeFunction& b_in, PdeFunction& b_interface, PdeFunction& source, PdeFunction& Dw)
{
  try {
    using namespace DROPS;

    std::ifstream param;
    std::cout << "Using default parameter file: poissonex1.json\n";
    param.open( "poissonex1.json");
    /*    else
      param.open( argv[1]);
    if (!param){
      std::cerr << "error while opening parameter file\n";
      return 1;
      }*/
    param >> P;
    param.close();
    std::cout << P << std::endl;

    // set up data structure to represent a poisson problem
    // ---------------------------------------------------------------------
    std::cout << line << "Set up data structure to represent a Poisson problem ...\n";
    int Nx, Ny, Nz, Ns, Nt, N;

    Nx = P.get<int>("DomainCond.nx")+1;
    Ny = P.get<int>("DomainCond.ny")+1;
    Nz = P.get<int>("DomainCond.nz")+1;
    Ns = Nx*Ny*Nz;
    Nt = P.get<int>("Time.NumSteps")+1;
    N  = Ns*Nt;

    if (check_dimensions(Nx,Ny,Nz,Nt,C0,b_in,b_interface,source,Dw)) {
      convection_diffusion(P, C0.get_data(), b_in.get_data(), b_interface.get_data(), source.get_data(), Dw.get_data(), NULL);
    }
    else {
      return false;
    }
    return true;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


BOOST_PYTHON_MODULE(drops)
{
    using namespace boost::python;
    def("convection_diffusion", py_convection_diffusion);
    //boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    //def("setArray", &setArray);

    class_<PdeFunction>("PdeFunction", init<int,int,int,int>(args("Nx","Ny","Nz","Nt"), "__init__ docstring"))
      .def(init<int, int,int,int>())
      .def("set_value", &PdeFunction::set_value);
}
