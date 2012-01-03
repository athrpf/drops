#include "Python.h"
#include "numpy/arrayobject.h"	/* NumPy header */

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
DROPS::ParamCL P;

array numpy_convection_diffusion(array& C0, array& b_in, array& source, array& Dw, array& b_interface, string json_filename) {
  // 1. Read parameter file
  std::ifstream param;
  param.open(json_filename.c_str());
  param >> P;
  param.close();
  using namespace boost::python;
  tuple t = extract<tuple>(source.attr("shape"));
  int nx = extract<int>(t[0]);
  int ny = extract<int>(t[1]);
  int nz = extract<int>(t[2]);
  int nt = extract<int>(t[3]);
  P.put<std::string>("PoissonCoeff.Reaction", "Zero");
  P.put<int>("Poisson.SolutionIsKnown",0);
  P.put<int>("DomainCond.RefineSteps", 0);
  P.put<std::string>("DomainCond.BoundaryType","2!21!21!2!21!21");
  P.put<int>("DomainCond.GeomType", 1);

  // read from meshfile string
  std::string meshfile = P.get<std::string>("DomainCond.MeshFile"), delim("x@");
  size_t idx;
  while ((idx = meshfile.find_first_of(delim)) != std::string::npos )
    meshfile[idx] = ' ';
  std::istringstream brick_info(meshfile);
  int nx_mesh, ny_mesh, nz_mesh;
  double lx_mesh, ly_mesh, lz_mesh;
  brick_info >> lx_mesh >> ly_mesh >> lz_mesh >> nx_mesh >> ny_mesh >> nz_mesh;
  assert(nx_mesh==nx-1);
  assert(ny_mesh==ny-1);
  assert(nz_mesh==nz-1);

  int nt_test = P.get<int>("Time.NumSteps"); // again, number of intervals
  assert(nt_test==nt-1);

  std::string ofilename = P.get<std::string>("Err.Output", "outfile.out");
  std::ofstream outfile;  outfile.open(ofilename.c_str());

  /* Print out parameters */
  outfile << P << std::endl;

  /* Convert numpy arrays to PyPdeFunctions */
  PdeFunction::ConstPtr C0f(new PyPdeBoundaryFunction(&C0,3));
  PdeFunction::ConstPtr b_inf(new PyPdeBoundaryFunction(&b_in,0));
  PdeFunction::ConstPtr b_interfacef(new PyPdeBoundaryFunction(&b_interface,1));
  PdeFunction::ConstPtr sourcef(new PyPdeFunction(&source));
  PdeFunction::ConstPtr Dwf(new PyPdeFunction(&Dw));

  if (!check_dimensions(nx,ny,nz,nt,*C0f,*b_inf,*b_interfacef,*sourcef,*Dwf)) {
    std::cerr <<"Error in setting up DROPS: Wrong dimensions in inputs!\n";
    abort();
  }
  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = nx; c_sol_dim[1] = ny; c_sol_dim[2] = nz; c_sol_dim[3] = nt;
  PyArrayObject* newarray = (PyArrayObject*) PyArray_New(&PyArray_Type, 4, c_sol_dim, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
  delete[] c_sol_dim;
  boost::python::object obj(handle<>((PyObject*)newarray));
  double* solution_ptr = (double*)newarray->data;
  for (int k=0; k<nx*ny*nz*nt; ++k) { solution_ptr[k] = -1.2345;} // for testing purposes
  // if the problem is time-dependent, the first time step is set to the initial value
  bool adjoint = !(std::string("IA1Adjoint").compare(P.get<std::string>("PoissonCoeff.IAProb")));

  if (nt>1) {
    double* c0_ptr = NULL;
    if (!adjoint)
      c0_ptr = solution_ptr;
    else
      c0_ptr = solution_ptr + (nt-1)*nx*ny*nz;
    for (int ix=0; ix<nx; ++ix) for (int iy=0; iy<ny; ++iy) for (int iz=0; iz<nz; ++iz) solution_ptr[ix+nx*iy+nx*ny*iz] = C0f->operator()(ix, iy, iz,0);
  }
  array solution = extract<boost::python::numeric::array>(obj);
  double* solution_without_intitial_ptr = (adjoint | nt<2) ? solution_ptr : solution_ptr+nx*ny*nz;
  convection_diffusion(outfile, P, C0f, b_inf, b_interfacef, sourcef, Dwf, solution_without_intitial_ptr);
  outfile.close();
  return solution;
}

#include "prepy_product.cpp"

void drops_err_translator(const DROPS::DROPSErrCL& err)
{
  std::stringstream ss;
  ss << "DROPSErr: ";
  err.what(ss);
  PyErr_SetString(PyExc_UserWarning, ss.str().c_str());
}

void pydrops_err_translator(const PyDropsErr& err)
{
  std::string filename = err.write_err_to_file();
  std::string msg;
  msg += "DROPSErr: The following error occured during a call to DROPS:\n" + err.msg + "\nThe state of DROPS was written to " + filename + "\n";
  PyErr_SetString(PyExc_UserWarning, msg.c_str());
}

BOOST_PYTHON_MODULE(drops)
{
  import_array();		/* Initialize the Numarray module. */
  using namespace boost::python;
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("convection_diffusion", numpy_convection_diffusion);
  //register_exception_translator<DROPS::DROPSErrCL>(&drops_err_translator);
  register_exception_translator<PyDropsErr>(&pydrops_err_translator);
  def("setup_scalar_product_matrices", setup_sp_matrices);
  def("scalar_product", numpy_scalar_product);

  /*  class_<PyPdeFunction>("PyPdeFunction", init<numeric::array&>(args("x"), "__init__ docstring"))
    .def(init<numeric::array&>())
    .def("at",&PyPdeFunction::at);*/
}
