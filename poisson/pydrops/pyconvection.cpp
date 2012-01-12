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

enum ProblemType {
  IA1Direct=1,
  IA1Adjoint=2,
  IA2Direct=3,
  IA2Adjoint=4,
  IA2Sensitivity=5,
  IA2Gradient=6
};

static std::map<std::string, ProblemType> problem_type_map;


void initialize_problem_type_map() {
  problem_type_map["IA1Direct"] = IA1Direct;
  problem_type_map["IA1Adjoint"] = IA1Adjoint;
  problem_type_map["IA2Direct"] = IA2Direct;
  problem_type_map["IA2Adjoint"] = IA2Adjoint;
  problem_type_map["IA2Sensitivity"] = IA2Sensitivity;
  problem_type_map["IA2Gradient"] = IA2Gradient;
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

/**
 *  Solve a stationary convection diffusion problem from python.
 *
 *  \param b_in: boundary condition at inlet
 *  \param source: Depends on the problem to be solved:
 *          Direct: source term (rhs)
 *          Sensitivity: delta (perturbance)
 *          Adjoint: error in u-u_m (rhs)
 *          Gradient: psi (solution of adjoint problem
 *  \param presol: Depends on the problem to be solved:
 *          Direct: unused
 *          Sensitivity: solution to the corresponding direct problem
 *          Adjoint: unused
 *          Gradient: solution to the corresponding direct problem
 *  \param Dw: effective diffusion coefficient
 *  \param b_interface: boundary condition at surface
 *  \param json_filename: filename for json file to be used for parameters
 */
array numpy_stationary_convection_diffusion(array& b_in, array& source, boost::python::object& presol_obj, array& Dw, array& b_interface, string json_filename) {
  using namespace boost::python;
  initialize_problem_type_map();
  PdeFunction::ConstPtr presolf;
  if (presol_obj.ptr()==object().ptr()) {
    std::cout << "presol is None\n";
    // None
  } else {
    array presol = extract<boost::python::numeric::array>(presol_obj);
    presolf = PdeFunction::ConstPtr(new PyPdeFunction(&presol));
  }

  // read json file
  std::ifstream param;
  param.open(json_filename.c_str());
  param >> P;
  param.close();
  tuple t = extract<tuple>(source.attr("shape"));
  int nx = extract<int>(t[0]);
  int ny = extract<int>(t[1]);
  int nz = extract<int>(t[2]);
  //int nt = extract<int>(t[3]);
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

  /* Open output stream and print out parameters */
  std::string ofilename = P.get<std::string>("Err.Output", "outfile.out");
  std::ofstream outfile;  outfile.open(ofilename.c_str());
  outfile << P << std::endl;

  /* Convert numpy arrays to PyPdeFunctions */
  PdeFunction::ConstPtr b_inf(new PyPdeBoundaryFunction(&b_in,0));
  PdeFunction::ConstPtr b_interfacef(new PyPdeBoundaryFunction(&b_interface,1));
  PdeFunction::ConstPtr sourcef(new PyPdeFunction(&source));
  PdeFunction::ConstPtr Dwf(new PyPdeFunction(&Dw));

  //if (!check_dimensions(nx, ny, nz, 1, *C0f, *b_inf, *b_interfacef, *sourcef*, *Dwf)) {

  // set up solution array
  npy_intp* c_sol_dim = new npy_intp[4];
  c_sol_dim[0] = nx; c_sol_dim[1] = ny; c_sol_dim[2] = nz; c_sol_dim[3] = 1;
  PyArrayObject* newarray = (PyArrayObject*) PyArray_New(&PyArray_Type, 4, c_sol_dim, PyArray_DOUBLE, NULL, NULL, 0, NPY_F_CONTIGUOUS, NULL);
  delete[] c_sol_dim;
  boost::python::object obj(handle<>((PyObject*)newarray));
  double* solution_ptr = (double*)newarray->data;
  for (int k=0; k<nx*ny*nz; ++k) { solution_ptr[k] = -1.2345;} // for testing purposes
  std::string problem_type = P.get<std::string>("PoissonCoeff.IAProb");
  array solution = extract<boost::python::numeric::array>(obj);
  if (problem_type_map.find(problem_type)==problem_type_map.end()) {
    outfile.close();
    throw (PyDropsErr("The provided problem type (PoissonCoeff.IAProb) does not exist or is not an IA2 problem."));
  }
  switch (problem_type_map[std::string("IA2Direct")]) {
  case IA2Direct:
    std::cout << "solving the direct problem\n";
    CoefEstimation(outfile, P, b_inf, b_interfacef, sourcef, PdeFunction::ConstPtr(), PdeFunction::ConstPtr(), Dwf, solution_ptr);
    std::cout << "Done solving the direct problem\n";
    break;
  case IA2Adjoint:
    CoefEstimation(outfile, P, b_inf, b_interfacef, sourcef, PdeFunction::ConstPtr(), PdeFunction::ConstPtr(), Dwf, solution_ptr);
    break;
  case IA2Sensitivity:
    CoefEstimation(outfile, P, b_inf, b_interfacef, PdeFunction::ConstPtr(), presolf, sourcef, Dwf, solution_ptr);
    break;
  case IA2Gradient:
    CoefEstimation(outfile, P, b_inf, b_interfacef, PdeFunction::ConstPtr(), presolf, sourcef, Dwf, solution_ptr);
    break;
  default:
    outfile.close();
    throw (PyDropsErr("The provided problem type (PoissonCoeff.IAProb) does not exist or is not an IA2 problem."));
  }
  outfile.close();
  return solution;
}

/**
 *  Solve a non-stationary convection diffusion problem from python.
 *
 *
 */
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

  /* Open output stream and print out parameters */
  std::string ofilename = P.get<std::string>("Err.Output", "outfile.out");
  std::ofstream outfile;  outfile.open(ofilename.c_str());
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
    for (int ix=0; ix<nx; ++ix) for (int iy=0; iy<ny; ++iy) for (int iz=0; iz<nz; ++iz) c0_ptr[ix+nx*iy+nx*ny*iz] = C0f->operator()(ix, iy, iz,0);
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
  def("stationary_cd", numpy_stationary_convection_diffusion);
  //register_exception_translator<DROPS::DROPSErrCL>(&drops_err_translator);
  register_exception_translator<PyDropsErr>(&pydrops_err_translator);

  class_<PyScalarProductConnector>("ScalarProductConnector", "The ScalarProductConnector is an interface for computing the scalar product on a specific domain.", init<int, int, int, int, double, double, double, double, bool>(args("Nx","Ny","Nz","Nt","lx","ly","lz","tmax","h1"), "The constructor takes as arguments the number of grid points and the lengths of each side of the box."))
    .def("scalar_product", &PyScalarProductConnector::numpy_scalar_product);

  /*  class_<PyPdeFunction>("PyPdeFunction", init<numeric::array&>(args("x"), "__init__ docstring"))
    .def(init<numeric::array&>())
    .def("at",&PyPdeFunction::at);*/
}
