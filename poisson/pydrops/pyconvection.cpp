#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include "misc/params.h"
#include <boost/shared_ptr.hpp>

#include "pyconnect.h"
#include "convection_diffusion.cpp"

class PyPdeFunction : public PdeFunction {
private:
  int nx, ny, nz, nt;
  boost::python::numeric::array data;
public:

  PyPdeFunction(boost::python::numeric::array& data_) : data(data_)
  {
    using namespace boost::python;
    tuple t = extract<tuple>(data.attr("shape"));
    assert(len(t)==4);
    nx = boost::python::extract<int>(t[0]);
    ny = boost::python::extract<int>(t[1]);
    nz = boost::python::extract<int>(t[2]);
    nt = boost::python::extract<int>(t[3]);
  }

  virtual ~PyPdeFunction(){}

  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const {
    bool retval = false;
    if (Nx==nx && Ny==ny && Nz==nz && Nt==nt) {
      retval = true;
    }
    Nx = nx; Ny = ny; Nz = nz; Nt = nt;
    return retval;
  }

  virtual double operator()(int ix, int iy, int iz, int it) const {
    assert (ix>0 && ix<nx);
    assert (iy>0 && iy<ny);
    assert (iz>0 && iz<nz);
    assert (it>0 && it<nt);
    using namespace boost::python;
    tuple t = make_tuple<int,int,int,int>(ix,iy,iz,it);
    return extract<double>(data[t]);
  }

  double at(int ix, int iy, int iz, int it) const {
    assert (ix>0 && ix<nx);
    assert (iy>0 && iy<ny);
    assert (iz>0 && iz<nz);
    assert (it>0 && it<nt);
    using namespace boost::python;
    tuple t = make_tuple<int,int,int,int>(ix,iy,iz,it);
    return extract<double>(data[t]);
  }
};
/*
bool get_array(boost::python::numeric::array& a) {
  std::cout << "got an array\n" ;
  std::cout << "shape[0] = " << boost::python::extract<int>(data.getshape()[0]);
  std::cout << boost::python::extract<double>(a[boost::python::make_tuple(0,0,0,0)]) << std::endl;
  return true;
}
*/

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

bool numpy_convection_diffusion(array& C0, array& b_in, array& b_interface, array& source, array& Dw, double uN, double Dmol,
			   double lx, double ly, double lz,
			   double dt, double theta, bool flag_pr, bool flag_bc, bool flag_supg) {
  //P.put<int>("DomainCond.nx", nx);
//  P.put<int>("DomainCond.ny", ny);
//  P.put<int>("DomainCond.nz", nz);
//  P.put<int>("DomainCond.lx", lx);
//  P.put<int>("DomainCond.ly", ly);
//  P.put<int>("DomainCond.lz", lz);
  //  P.put<int>("Poisson.SolutionIsKnown",0);
  using namespace boost::python;
  tuple t = extract<tuple>(source.attr("shape"));
  int nx = extract<int>(t[0]);
  int ny = extract<int>(t[1]);
  int nz = extract<int>(t[2]);
  int nt = extract<int>(t[3]);
  P.put<int>("DomainCond.RefineSteps", 0);
  std::stringstream MeshFile;
  MeshFile << lx << "x" << ly << "x" << lz << "@" << nx << "x" << ny << "x" << nz;
  P.put<std::string>("DomainCond.MeshFile",MeshFile.str());
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
  PdeFunPtr C0f(new PyPdeFunction(C0));
  PdeFunPtr b_inf(new PyPdeFunction(b_in));
  PdeFunPtr b_interfacef(new PyPdeFunction(b_interface));
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
