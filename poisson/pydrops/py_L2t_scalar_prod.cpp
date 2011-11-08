//**************************************************************************
// File:    pyL2t_scalar_product.cpp                          	           *
// Content: Python interface Calculation of the scalar product             *
//          in L2(\Omega X [t_0,t_f])                                      *
// Author:  Maka Karalashvili, Hans Pirnay				   *
//          RWTH Aachen                                                    *
// Version: 0.1                                                            *
// History: begin - Aug, 4 2011                                            *
//**************************************************************************

#include "Python.h"
#include "numpy/arrayobject.h"	/* NumPy header */

#include "poisson/poisson.h"
#include "poisson/integrTime.h"

#include <fstream>

#include "py_utils.hpp"
#include "drops_utils.hpp"
#include "py_param.hpp"

/**
   This code calculates the scalar product of two functions u and v over Tetrahedras
   The program can be called as a mex-function from matlab with the following syntax:

   [scalar_prod] = L2t_scalar_prod(u, v, xl, yl, zl, nx, ny, nz, dt, nt);

   scalar parameters:
   ------------------
   xl, yl, zl:     length of brick in x-/y-/z-direction [in mm]
   nx, ny, nz:     number of intervals in x-/y-/z-direction, i.e. Nxyz = (nx+1) x (ny+1) x (nz+1) grid points
   dt, nt:         length and number of time steps
   flag:           the flag of shoice of the quadrature for the time integration: flag=0.5 - the trapezium, flag = 1 - the quadrant

   matrix parameters:
   ------------------
   u:   Nxyz x nt        function in \Omega X [t_0,t_f]
   v:   Nxyz x nt        function in \Omega X [t_0,t_f]
**/

/*
class ParamCL
{
public:
  double lx, ly, lz;
  int    nx, ny, nz, nt; // n=number of intervalls
  double dt;

  ParamCL()
    : lx(1.0), ly(1.0), lz(1.0), nx(8), ny(2), nz(2), nt(50), // in mm
      dt(0.02)
  {}

  }*/

InstatParamCL C_L2SP;

class MatlabConnectCL
{ // holds the Matlab input matrices and Matlab output parameters

  typedef std::pair<float, float> f_pair;
  typedef std::pair<float, f_pair> cmp_key;
  //
  typedef std::map<cmp_key, DROPS::FaceCL*>  FACE_MAP;
  typedef std::map<cmp_key, DROPS::TetraCL*> TETRA_MAP;

private:
  int Nx, Ny, Nz, Nxy, Nxz, Nyz, Nxyz, Nt; // N=number of points
  double dx, dy, dz;

  const double *u, *v;   // input  matrices: u,v
  double scalar_prod;   // output scalar: the scalar product (1 x 1)
  //
  static DROPS::MultiGridCL* _MG;
  static TETRA_MAP tetra_map;
  typedef TETRA_MAP::const_iterator ti;

  // c-style array indexing for full-space vectors
  int GetNum( const DROPS::Point3DCL& p, double t) const
  {
    int retval = (rd(p[0]/dx)*Ny*Nz*Nt + rd(p[1]/dy)*Nz*Nt + rd(p[2]/dz)*Nt + rd(t/C_L2SP.dt_));
    if (retval <0 || retval >=Nx*Ny*Nz*Nt) {
      std::cout << "ERROR in GetNum: retval = " << retval << ", (" << p[0] << ","<< p[1] << ","<< p[2] << ","<< t << ")\nError Info: C_L2SP.dt_ = " << C_L2SP.dt_ << "\nNt = " << Nt << "\n";
      throw (retval);
    }
    return retval;
  }
  /*
    int GetNum( const DROPS::Point3DCL& p, double t=0.) const
    {
    return (rd(p[2]/dz)*Nxy + rd(p[1]/dy)*Nx + rd(p[0]/dx) + rd(t/C_L2SP.dt)*Nxyz);
    }
  */

public:
  MatlabConnectCL()
  {
    tetra_map.clear();
  }
  void SetMG( DROPS::MultiGridCL* MG) { _MG= MG; }
  //
  static void ClearMaps() {
    tetra_map.clear();
  }
  //
  static void setTetraMap()
  {
    //mexPrintf("L2: SETTING TETRA MAP\n");
    DROPS::Uint lvl= _MG->GetLastLevel();
    for(DROPS::MultiGridCL::TriangTetraIteratorCL
	  tit=_MG->GetTriangTetraBegin(lvl), tend=_MG->GetTriangTetraEnd(lvl);
        tit != tend; ++tit) {
      DROPS::TetraCL& tetra = *tit;

      DROPS::Point3DCL bc = DROPS::GetBaryCenter(tetra);
      f_pair pr= std::make_pair(rnd(bc[2]), rnd(bc[1]));
      cmp_key key= std::make_pair(rnd(bc[0]), pr);

      tetra_map[key]= &tetra;
    }
    //mexPrintf("L2: TETRA MAP SET\n");
  }
  static void DumpTetraMap()
  {
    FILE* f=fopen("tetra.map","w");
    //mexPrintf("DUMPING TETRA MAP\n");
    for (ti p= tetra_map.begin(); p!= tetra_map.end(); p++) {

      DROPS::Point3DCL bc = DROPS::Point3DCL(0.);
      bc[0]=p->first.first;
      bc[1]=p->first.second.second;
      bc[2]=p->first.second.first;

      f_pair pr= std::make_pair(bc[2], bc[1]);
      cmp_key key= std::make_pair(bc[0], pr);

      DROPS::TetraCL* tetra = tetra_map[key];
      fprintf(f,"key: %f %f %f \n",bc[0],bc[1],bc[2]);

      fprintf(f,"tetra:\n");
      DROPS::Point3DCL p = tetra->GetVertex(0)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = tetra->GetVertex(1)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = tetra->GetVertex(2)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = tetra->GetVertex(3)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
    }
    fclose(f);
    //mexPrintf("END DUMP TETRA MAP\n");
  }
  //functions
  double GetU( const DROPS::Point3DCL& p, double t=0.) const
  {
    double ret;

    f_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    DROPS::TetraCL* tetra= tetra_map[key];

    if (tetra == NULL) {//non-barycenter
      ret=u[GetNum(p,t)];
    }else {
      ret = 0.25*(u[GetNum(tetra->GetVertex(0)->GetCoord(),t)]+u[GetNum(tetra->GetVertex(1)->GetCoord(),t)]+
		  u[GetNum(tetra->GetVertex(2)->GetCoord(),t)]+u[GetNum(tetra->GetVertex(3)->GetCoord(),t)]) ;
    }

    return ret;
  };
  //
  double GetV( const DROPS::Point3DCL& p, double t=0.) const
  {
    double ret;

    f_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    DROPS::TetraCL* tetra= tetra_map[key];

    if (tetra == NULL) {//non-barycenter
      ret=v[GetNum(p,t)];
    }else {
      ret = 0.25*(v[GetNum(tetra->GetVertex(0)->GetCoord(),t)]+v[GetNum(tetra->GetVertex(1)->GetCoord(),t)]+
		  v[GetNum(tetra->GetVertex(2)->GetCoord(),t)]+v[GetNum(tetra->GetVertex(3)->GetCoord(),t)]) ;
    }

    return ret;
  };
  //
  //initialisation
  void Init( const InstatParamCL& P, const double* u_ptr, const double* v_ptr) {
    Nx= P.nx_+1;Ny= P.ny_+1; Nz=P.nz_+1;
    Nyz=Ny*Nz; Nxy=Nx*Ny; Nxz=Nx*Nz;
    Nxyz= Nxy*Nz;
    Nt = C_L2SP.nt_+1;
    dx= P.lx_/P.nx_; dy= P.ly_/P.ny_; dz= P.lz_/P.nz_;

    u =  u_ptr;
    v =  v_ptr;
  }
} MC_L2SP;

DROPS::MultiGridCL* MatlabConnectCL::_MG= NULL;
MatlabConnectCL::TETRA_MAP MatlabConnectCL::tetra_map;

//
namespace DROPS
{
  static double fu(const DROPS::Point3DCL& p, double t=0.)
  {
    return MC_L2SP.GetU(p,t);
  }
  static double fv(const DROPS::Point3DCL& p, double t=0.)
  {
    return MC_L2SP.GetV(p,t);
  }

  void L2tScalarProduct(MultiGridCL& MG, double* l2)
  {
    double u, v, t;
    Uint lvl=MG.GetLastLevel();
    //
    for (int step=0;step<=C_L2SP.nt_;step++) {
      t=step*C_L2SP.dt_;
      //
      double L2_t=0.0;
      for (MultiGridCL::TriangTetraIteratorCL sit(MG.GetTriangTetraBegin(lvl)),
	     send(MG.GetTriangTetraEnd(lvl)); sit!=send; ++sit){

	// why *6 ?
	double absdet= sit->GetVolume()*6.,
	  sum= 0;

	for(Uint i=0; i<4; ++i)
	  {
	    u= fu(sit->GetVertex(i)->GetCoord(),t);
	    v= fv(sit->GetVertex(i)->GetCoord(),t);
	    sum+= u*v;
	  }
	sum/= 120; // vertex-multiplier?
	u= fu(GetBaryCenter(*sit),t);
	v= fv(GetBaryCenter(*sit),t);
	// 2/15 is barycenter-multiplier?
	sum+= 2./15. * u*v;
	L2_t+= sum*absdet;
      }
      //
      l2[step] = L2_t;
    }
  }

} // end of namespace DROPS

static double L2t_scalar_prod(double* l2_scalar_prod)
{ // create Multigrid  and call ScalarProduct(...)
  try
    {
      DROPS::Point3DCL null(0.0);
      DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
      e1[0]= C_L2SP.lx_;
      e2[1]= C_L2SP.ly_;
      e3[2]= C_L2SP.lz_;

      DROPS::BrickBuilderCL brick(null, e1, e2, e3, C_L2SP.nx_, C_L2SP.ny_, C_L2SP.nz_);
      DROPS::MultiGridCL MG(brick);

      //prepare MC_L2SP
      MC_L2SP.SetMG(&MG);
      MatlabConnectCL::ClearMaps();
      //MatlabConnectCL::setFaceMap();
      MatlabConnectCL::setTetraMap();

      DROPS::L2tScalarProduct(MG, l2_scalar_prod);
    }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
  return -1.;
}

//the mexFunction
static PyObject* drops_L2t_scalar_prod(PyObject* self, PyObject* args)
{
  PyArrayObject* u=NULL;
  PyArrayObject* v=NULL;
  double lx, ly, lz, tmax;
  // get the input arguments
  int pyarg_retval = PyArg_ParseTuple(args, "O!O!dddd",
				      &PyArray_Type, &u,
				      &PyArray_Type, &v,
				      &lx, &ly, &lz, &tmax
				      );
  if (!pyarg_retval) {
    PyErr_SetString(PyExc_TypeError, "The input arguments are not correct.\nThe correct function call is sp = L2t_scalar_prod(u,v,lx,ly,lz,tmax)");
    return NULL;
  }

  C_L2SP.lx_ = lx;
  C_L2SP.ly_ = ly;
  C_L2SP.lz_ = lz;

  // dimensions
  int nx_u, nx_v;
  int ny_u, ny_v;
  int nz_u, nz_v;
  int nt_u, nt_v;
  nx_u = u->dimensions[0];
  nx_v = v->dimensions[0];
  ny_u = u->dimensions[1];
  ny_v = v->dimensions[1];
  nz_u = u->dimensions[2];
  nz_v = v->dimensions[2];
  nt_u = u->dimensions[3];
  nt_v = v->dimensions[3];
  if (nx_u!=nx_v || ny_u!=ny_v || nz_u!=nz_v || nt_u!=nt_v) {
    PyErr_SetString(PyExc_ValueError, "Vectors u and v do not have the same dimensions");
  }
  C_L2SP.nx_ = nx_u-1;
  C_L2SP.ny_ = ny_u-1;
  C_L2SP.nz_ = nz_u-1;
  C_L2SP.nt_ = nt_u-1;

  if (tmax==0.0) {
    C_L2SP.dt_ = 1.0;
  } else {
    C_L2SP.dt_ = tmax/C_L2SP.nt_;
  }

  //
  // Set the input matrices and output parameters.
  MC_L2SP.Init(C_L2SP, (const double*)u->data, (const double*)v->data);

  npy_intp l2_length = nt_u;
  PyArrayObject* l2 = (PyArrayObject*) PyArray_SimpleNew(1, &l2_length, PyArray_DOUBLE);
  double* l2_ptr    = (double*) l2->data;
  // Call the subroutine.
  L2t_scalar_prod(l2_ptr);

  PyObject* retval = Py_BuildValue("N",	PyArray_Return(l2));
  return retval;
}


static PyMethodDef DropsMethods[] = {
    {"l2_scalar",  drops_L2t_scalar_prod, METH_VARARGS,
     "Compute the scalar product."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initscalarprod(void)
{
    (void) Py_InitModule("scalarprod", DropsMethods);
    import_array();		/* Initialize the Numarray module. */
    /* A segfault will occur if I use numarray without this.. */
}
