#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/shared_ptr.hpp>

#include "misc/params.h"

#include "pyconnect.h"
#include "convection_diffusion.cpp"
#include "pypdefunction.h"

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


using namespace boost::python::numeric;

bool numpy_convection_diffusion(array& C0, array& b_in, array& source, array& Dw, array& b_interface, double uN, double Dmol,
				double lx, double ly, double lz,
				double dt, double theta, bool flag_pr, bool flag_bc, bool flag_supg) {
  using namespace boost::python;
  tuple t = extract<tuple>(source.attr("shape"));
  int nx = extract<int>(t[0]);
  int ny = extract<int>(t[1]);
  int nz = extract<int>(t[2]);
  int nt = extract<int>(t[3]);
  P.put<int>("DomainCond.nx", nx);
  P.put<int>("DomainCond.ny", ny);
  P.put<int>("DomainCond.nz", nz);
  P.put<int>("DomainCond.lx", lx);
  P.put<int>("DomainCond.ly", ly);
  P.put<int>("DomainCond.lz", lz);
  P.put<int>("PoissonCoeff.Reaction", 0.0);
  P.put<int>("Poisson.SolutionIsKnown",0);
  P.put<int>("DomainCond.RefineSteps", 0);
  //std::stringstream MeshFile;
  //MeshFile << lx << "x" << ly << "x" << lz << "@" << nx << "x" << ny << "x" << nz;
  //P.put<std::string>("DomainCond.MeshFile",MeshFile.str());
  P.put<std::string>("DomainCond.BoundaryType","0!2!2!0!2!2");
  P.put<int>("DomainCond.GeomType", 1);

  P.put<double>("PoissonCoeff.Dmol", Dmol);
  P.put<bool>("PoissonCoeff.Stabilization", flag_supg);

  P.put<int>("Time.NumSteps", nt);
  P.put<double>("Time.StepSize", dt);
  P.put<int>("Time.Scheme", 1);
  P.put<double>("Time.Theta", 1.0);
  P.put<bool>("Time.Convection", true);

  /* Convert numpy arrays to PyPdeFunctions */
  typedef const PdeFunction* PdeFunPtr;
  PdeFunPtr C0f(new PyPdeBoundaryFunction(C0,3));
  PdeFunPtr b_inf(new PyPdeBoundaryFunction(b_in,0));
  PdeFunPtr b_interfacef(new PyPdeBoundaryFunction(b_interface,1));
  PdeFunPtr sourcef(new PyPdeFunction(source));
  PdeFunPtr Dwf(new PyPdeFunction(Dw));
  convection_diffusion(P, C0f, b_inf, b_interfacef, sourcef, Dwf, NULL);
  return true;
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
      convection_diffusion(P, &C0, &b_in, &b_interface, &source, &Dw, NULL);
    }
    else {
      return false;
    }
    return true;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
  return false;
}


BOOST_PYTHON_MODULE(drops)
{
    using namespace boost::python;
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    def("convection_diffusion", numpy_convection_diffusion);
    //def("get_array", &get_array);
    //def("setArray", &setArray);

    class_<PyPdeFunction>("PyPdeFunction", init<numeric::array&>(args("x"), "__init__ docstring"))
      .def(init<numeric::array&>())
      .def("at",&PyPdeFunction::at);
}
