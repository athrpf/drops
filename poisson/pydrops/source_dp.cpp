//********************************************************************************
// File:    source_dp.cpp                                                         *
// Content: poisson with convection and SUPG stabilization + interface for matlab *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers, Liang Zhang*
//          IGPM RWTH Aachen
//			Maka Karalashvili                                                     *
//          AVT.PT
// Version: 0.1                                                            *
// History: begin - Jul, 01 2008 modified 06.2011                                 *
//********************************************************************************

#include "mex.h"

#include "poisson/poisson.h"
#include "num/solver.h"
#include "integrTime.h"
#include "out/output.h"
#include <fstream>
#include <string>
#include <sstream>
//#include "boost/lexical_cast.hpp"
//using boost::lexical_cast;
using std::string;


const int SensFlagC= 64, AdjFlagC=32, DirBCFlagC= 16;
/**
   This code has the same functionality as poisson/ipfilm.cpp, but is extended by a matlab interface.
   It solves a convection-diffusion equation on a brick-shaped domain with certain bnd conditions
   modeling the (effective) heat/mass transfer in a flat film. Bnd conditions:
   - Dirichlet bc at x=0 (inflow)
   - Natural   bc at interface (Wall for heat and free surface for mass)
   - all other bnds: homogeneous natural bc

   The program can be called as a mex-function from matlab with the following syntax:

   [MaxIter, Csol] = source_dp( C0, B_in, F, alpha, B_interface, uN, amol, xl, yl, zl, nx, ny, nz, dt, nt, Theta, Tol, Iter,
   Flag_pr, Flag_bc, Flag_SUPG)

   scalar parameters:
   ------------------
   xl, yl, zl:     length of brick in x-/y-/z-direction [in mm]
   nx, ny, nz:     number of intervals in x-/y-/z-direction, i.e. Nxyz = (nx+1) x (ny+1) x (nz+1) grid points
   dt, nt:         length and number of time steps
   uN:             maximum velocity for Nusselt
   Theta:          parameter of one-step theta-scheme, controlling implicitness of time discretization
   Tol, Iter:      stopping criterion for the iterative solver (GMRES)
   a:          	  diffusion parameter a = lambda / (rho c) [in SI]
   //hflux:          bnd condition for heating: hflux = qh / (-lambda) [in SI]
   //                                       or  hflux = T_wall [K] if Dirichlet bc is used (cf. Flag)        //??
   Flag_pr:        to control if you solve the adjoint problem;
   Flag_bc:        to control the boundary type of the interface(in mass transport;
   Flag_SUPG:      to turn on stabilization, 1 with stabilization, o without.
   MaxIter (output):        number of iterations spent in the iterative solver (maximum over all time steps)

   matrix parameters:
   ------------------
   C0:     Nxyz x 1         initial temperature distribution T0(x,y,z)
   B_in:   Nyz  x Nt        temperature at inflow (x=0): T_in(y,z,t)
   F:      Nxyz x Nt        rhs term F(x,y,z,t)
   B_interface:             Interface boundary condition

   Csol (output):   Nxyz x nt        solution of the conv-diff problem,
   temperature distribution T(x,y,z,t) omitting initial time step

   Here Nt = nt+1,    Nyz = (ny+1) x (nz+1),    Nxyz = (nx+1) x Nyz.

   As one can see, the columns of the matrices correspond to the different time steps.
   Note, that all the values have to be passed in entities W, mm
**/

extern void _main();

inline int rd( double d) { return static_cast<int>( d+0.5); }                   // rounding

class ParamCL
{
public:
  double lx_, ly_, lz_;
  int    nx_, ny_, nz_, nt_;                                                      // n=number of intervals
  double dt_, hflux_, rho_, uN_, cp_, a_mol_, a_;
  std::string EnsDir_, EnsCase_;

  ParamCL()
    : lx_(180), ly_(0.8), lz_(0.3), nx_(8), ny_(2), nz_(2), nt_(50),                 // in mm
      dt_(0.02), hflux_(350), rho_(912), cp_(1540), a_mol_(0.118/1540/912), a_(1.0),// in SI
      EnsDir_("ensight"), EnsCase_("FilmTemp")
  {}

} C;

class MatlabConnectCL
{ // holds the Matlab input matrices and Matlab output parameters
  typedef std::pair<double, double> d_pair;
  typedef std::pair<double, d_pair> cmp_key;
  //
  typedef std::map<cmp_key, DROPS::FaceCL*>  FACE_MAP;
  typedef std::map<cmp_key, DROPS::TetraCL*> TETRA_MAP;

private:
  int Nx_, Ny_, Nz_, Nxy_, Nyz_, Nxz_, Nxyz_; // N=number of points
  double dx_, dy_, dz_;
  double a_mol_;

  const double *C0_, *B_in_, *B_Int_, *F_, // initial+boundary+rhs function,
    *a_;                                     // diffusion parameter as a function,
  double* C3D_,                                   // output matrices: temp solution (Nxyz x nt),
    *MaxIter_;                               // max. iterations of solver (1 x 1)
  //helper maps for barycenters
  static DROPS::MultiGridCL* MG_;

  static FACE_MAP face_map_;
  static TETRA_MAP tetra_map_;

  typedef FACE_MAP::const_iterator fi;
  typedef TETRA_MAP::const_iterator ti;


  static double rnd(double d) {// rounding four digits after comma
    int i_d = (int)(d*10000);
    double d_d = (double)i_d/10000;
    return d_d;
  }
  //
  int GetNum( const DROPS::Point3DCL& p, double t, int seg) const
  {
    if (seg == 0 || seg == 1){//yoz
      return (rd(p[2]/dz_)*Ny_ + rd(p[1]/dy_) + rd(t/C.dt_)*Ny_*Nz_);
    }
    if (seg == 2 || seg == 3){//xoz
      return (rd(p[2]/dz_)*Nx_ + rd(p[0]/dx_) + rd(t/C.dt_)*Nx_*Nz_);
    }
    if (seg == 4 || seg == 5){//xoy
      return (rd(p[1]/dy_)*Nx_ + rd(p[0]/dx_) + rd(t/C.dt_)*Nx_*Ny_);
    }
  }
  //
  int GetNum( const DROPS::Point3DCL& p, double t=0.) const
  {
    return (rd(p[2]/dz_)*Nxy_ + rd(p[1]/dy_)*Nx_ + rd(p[0]/dx_) + rd(t/C.dt_)*Nxyz_);
  }
public:
  MatlabConnectCL()
  {
    Nx_=Ny_=Nz_=Nxy_=Nyz_=Nxz_=Nxyz_=-1;
    dx_=dy_=dz_=0.0;

    face_map_.clear();
    tetra_map_.clear();
  }
  void SetMG( DROPS::MultiGridCL* MG) { MG_= MG; }
  //
  static void ClearMaps() {
    face_map_.clear();
    tetra_map_.clear();
  }
  //
  static void setFaceMap()
  {
    //mexPrintf("SETTING FACE MAP\n");
    DROPS::Uint lvl= MG_->GetLastLevel();
    for(DROPS::MultiGridCL::TriangFaceIteratorCL
	  fit=MG_->GetTriangFaceBegin(lvl), fend=MG_->GetTriangFaceEnd(lvl);
      	fit != fend;++fit) {
      DROPS::FaceCL& face = *fit;

      DROPS::Point3DCL bc = DROPS::GetBaryCenter(face);
      d_pair pr= std::make_pair(rnd(bc[2]), rnd(bc[1]));
      cmp_key key= std::make_pair(rnd(bc[0]), pr);

      face_map_[key]= &face;
    }
    //mexPrintf("FACE MAP SET\n");
  }
  static void DumpFaceMap()
  {
    FILE* f=fopen("face.map","w");
    //mexPrintf("DUMPING FACE MAP\n");
    for (fi p= face_map_.begin(); p!= face_map_.end(); p++) {

      DROPS::Point3DCL bc = DROPS::Point3DCL(0.);
      bc[0]=p->first.first;
      bc[1]=p->first.second.second;
      bc[2]=p->first.second.first;

      d_pair pr= std::make_pair(bc[2], bc[1]);
      cmp_key key= std::make_pair(bc[0], pr);

      DROPS::FaceCL* face = face_map_[key];

      fprintf(f,"bc: %f %f %f \n",bc[0],bc[1],bc[2]);

      fprintf(f,"face:\n");
      DROPS::Point3DCL p = face->GetVertex(0)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = face->GetVertex(1)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = face->GetVertex(2)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
    }
    fprintf(f," %d",(int)face_map_.size());
    fclose(f);
    mexPrintf("END DUMP FACE MAP\n");
  }
  //
  static void setTetraMap()
  {
    mexPrintf("SETTING TETRA MAP\n");
    DROPS::Uint lvl= MG_->GetLastLevel();
    for(DROPS::MultiGridCL::TriangTetraIteratorCL
	  tit=MG_->GetTriangTetraBegin(lvl), tend=MG_->GetTriangTetraEnd(lvl);
        tit != tend; ++tit) {
      DROPS::TetraCL& tetra = *tit;

      DROPS::Point3DCL bc = DROPS::GetBaryCenter(tetra);
      d_pair pr= std::make_pair(rnd(bc[2]), rnd(bc[1]));
      cmp_key key= std::make_pair(rnd(bc[0]), pr);

      tetra_map_[key]= &tetra;
    }
    mexPrintf("TETRA MAP SET\n");
  }
  static void DumpTetraMap()
  {
    FILE* f=fopen("tetra.map","w");
    mexPrintf("DUMPING TETRA MAP\n");
    for (ti p= tetra_map_.begin(); p!= tetra_map_.end(); p++) {

      DROPS::Point3DCL bc = DROPS::Point3DCL(0.);
      bc[0]=p->first.first;
      bc[1]=p->first.second.second;
      bc[2]=p->first.second.first;

      d_pair pr= std::make_pair(bc[2], bc[1]);
      cmp_key key= std::make_pair(bc[0], pr);

      DROPS::TetraCL* tetra = tetra_map_[key];
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
    mexPrintf("END DUMP TETRA MAP\n");
  }
  //
  double GetInitial( const DROPS::Point3DCL& p) const
  {
    return C0_[GetNum(p)];
  };
  //boundary functions
  //x=0;
  double GetInflow( const DROPS::Point3DCL& p, double t) const
  {
    return B_in_[GetNum(p,t,0)];
  };
  //y=0: if neumann condition is active
  double GetInterfaceFlux( const DROPS::Point3DCL& p, double t) const
  {
    double ret;

    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);

    DROPS::FaceCL* face= face_map_[key];

    if (face == NULL) {//non-barycenter
      ret= B_Int_[GetNum(p,t,3)];
    } else {
      ret= 1./3.*(B_Int_[GetNum(face->GetVertex(0)->GetCoord(),t,3)]+B_Int_[GetNum(face->GetVertex(1)->GetCoord(),t,3)]+B_Int_[GetNum(face->GetVertex(2)->GetCoord(),t,3)]);
    }
    return ret;
  };
  //y=0: if dirichlet condition is active
  double GetInterfaceValue( const DROPS::Point3DCL& p, double t) const
  {
    return B_Int_[GetNum(p,t,3)];
  };
  //rhs
  double GetRhs( const DROPS::Point3DCL& p, double t) const
  {
    double ret;

    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    DROPS::TetraCL* tetra= tetra_map_[key];

    if (tetra == NULL) {//non-barycenter
      ret=F_[GetNum(p,t)];
    }else {
      ret = 0.25*(F_[GetNum(tetra->GetVertex(0)->GetCoord(),t)]+F_[GetNum(tetra->GetVertex(1)->GetCoord(),t)]+
		  F_[GetNum(tetra->GetVertex(2)->GetCoord(),t)]+F_[GetNum(tetra->GetVertex(3)->GetCoord(),t)]) ;
    }

    return ret;
  };
  //coefficient functions
  double GetA( const DROPS::Point3DCL& p, double t) const
  {
    double ret;

    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    DROPS::TetraCL* tetra= tetra_map_[key];

    if (tetra == NULL) {//non-barycenter
      ret=a_[GetNum(p,t)]+a_mol_;
    }else {
      ret = 0.25*(a_[GetNum(tetra->GetVertex(0)->GetCoord(),t)]+a_[GetNum(tetra->GetVertex(1)->GetCoord(),t)]+
		  a_[GetNum(tetra->GetVertex(2)->GetCoord(),t)]+a_[GetNum(tetra->GetVertex(3)->GetCoord(),t)])  + a_mol_;
    }
    return ret;
  };

  template<class P1EvalT>
  void SetSol3D( const P1EvalT& sol, double t)
  {
    const int num= (rd(t/C.dt_)-1)*Nxyz_; // omit initial time step in output
    double *out= C3D_+num;

    DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
      {
	//mexPrintf("val: %f\n", sol.val(*sit));
	out[GetNum( sit->GetCoord())]= sol.val( *sit);
      }

  }

  //Check the input matrices
  void Init( const ParamCL& P, int Flag_bc, const double* C0, const double* B_in, const double* F, const double* a, const double* B_Int, double* maxiter, double* c_sol)
  {
    Nx_= P.nx_+1;Ny_= P.ny_+1; Nz_=P.nz_+1;
    Nyz_=Ny_*Nz_; Nxy_=Nx_*Ny_; Nxz_=Nx_*Nz_;
    Nxyz_= Nxy_*Nz_;
    dx_= P.lx_/P.nx_; dy_= P.ly_/P.ny_; dz_= P.lz_/P.nz_;

    // Save the matrix input arguments.
    C0_   = C0;
    B_in_ = B_in;
    F_    = F;
    a_    = a;
    a_mol_ = P.a_mol_;
    B_Int_= B_Int;
    if  (!(Flag_bc & DirBCFlagC)) {
      mexPrintf("solving with Neumann Cond.\n");
    } else {
      mexPrintf("solving with Dirichlet Cond.\n");
    }


    // Set the output pointer to the output arguments.
    MaxIter_ = maxiter;
    C3D_ = c_sol;
  }
} MC;

DROPS::MultiGridCL* MatlabConnectCL::MG_= NULL;
MatlabConnectCL::FACE_MAP MatlabConnectCL::face_map_;
MatlabConnectCL::TETRA_MAP MatlabConnectCL::tetra_map_;


// Solve the convection-diffusion equation
// du/dt - nu*laplace u + Vel grad u + q*u = f

class PoissonCoeffCL
{
public:
  static double alpha(const DROPS::Point3DCL& p, double t)
  {
    return MC.GetA(p,t);
  }
  static double f(const DROPS::Point3DCL& p, double t)
  {
    return MC.GetRhs(p,t);
  }
  //
  static double Initial(const DROPS::Point3DCL& p, double t)
  {
    return MC.GetInitial(p);
  }
  //This velocity profile is modified for the numerical testcase (as in SISC Paper)
  static DROPS::Point3DCL Vel(const DROPS::Point3DCL& p, double)
  {
    DROPS::Point3DCL ret;
    const double d= p[1]/C.ly_;
    ret[0]= C.uN_*(2-d)*d; // Nusselt
    return ret;
  }
  //Only used for Stabilization for flat film case
  static double h_Value()
  {//mesh size in flow direction
    double h=C.lx_/C.nx_;
    return h;
  }
  static double PecNum(const DROPS::Point3DCL& p, double t)
  {//Peclet Number
    double Pec=0.;
    Pec=fabs(Vel(p, t)[0])*h_Value()/(2.*alpha(p,t));
    return Pec;
  }
  static double Sta_Coeff(const DROPS::Point3DCL& p, double t)
  {//Stabilization coefficient
    if (PecNum(p,t)<=1)
      return 0.0;
    else
      return h_Value()/(2.*fabs(Vel(p, t)[0]))*(1.-1./PecNum(p, t));
  }
};

double Zero(const DROPS::Point3DCL&, double) { return 0.0; }
double Inflow(const DROPS::Point3DCL& p, double t) { return MC.GetInflow(p,t); }
double InterfaceFlux(const DROPS::Point3DCL& p, double t) { return MC.GetInterfaceFlux(p,t); }
double InterfaceValue(const DROPS::Point3DCL& p, double t) { return MC.GetInterfaceValue(p,t); }


namespace DROPS
{

  template<class Coeff>
  void Strategy(InstatPoissonP1CL<Coeff>& Poisson, double theta, double tol, int maxiter, int Flag_SUPG)
  {
    typedef InstatPoissonP1CL<Coeff> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();
    IdxDescCL& idx= Poisson.idx;
    VecDescCL& x= Poisson.x;
    VecDescCL& b= Poisson.b;
    MatDescCL& A= Poisson.A;
    MatDescCL& M= Poisson.M;
    MatDescCL& U= Poisson.U;

    //SUPG
    MatDescCL& M_SD=Poisson.M_SD;
    MatDescCL& U_SD=Poisson.U_SD;
    VecDescCL& vM_SD=Poisson.vM_SD;
    VecDescCL& vU_SD=Poisson.vU_SD;




    idx.Set(1);
    // erzeuge Nummerierung zu diesem Index
    Poisson.CreateNumbering(MG.GetLastLevel(), &idx);

    // Vektoren mit Index idx
    x.SetIdx(&idx);
    b.SetIdx(&idx);
    // Matrizen mit Index idx (Zeilen und Spalten)
    A.SetIdx(&idx, &idx);
    M.SetIdx(&idx, &idx);
    U.SetIdx(&idx, &idx);

    //SUPG
    vM_SD.SetIdx(&idx);
    vU_SD.SetIdx(&idx);
    M_SD.SetIdx(&idx, &idx);
    U_SD.SetIdx(&idx, &idx);

    mexPrintf("Number of Unknowns: %d\n", x.Data.size());
    mexPrintf("Theta: %g\n", theta);
    mexPrintf("Tolerance GMRES: %g\t", tol);
    mexPrintf("max. Num. GMRES-Iterations: %d\n", maxiter);
    /*  Point3DCL pIF; pIF[1]=C.ly; // Punkt auf Phasengrenze
	const double vel= norm( PoissonCoeffCL::Vel(pIF,0)),
	cfl= vel*C.dt/(C.lx_/C.nx),
	Re= (9.81*1e3)*std::pow(C.ly,3)/(3*C.nu*C.nu), //
	Pr= C.nu/C.a_mol;                              //??
	mexPrintf("Geschwindigkeit Phasengrenze: %g [mm/s]\tentspricht CFL = %g\n", vel, cfl);
	mexPrintf("Geschwindigkeit Phasengrenze: %g [mm/s]\tentspricht CFL = %g\nRe = %g, Pr = %g\n", vel, cfl, Re, Pr);
	mexPrintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");*/

    // stationaerer Anteil
    Poisson.SetupInstatSystem(A, M, Poisson.t);


    //parameters about solver?
    SSORPcCL pc(1.0);
    typedef GMResSolverCL<SSORPcCL> SolverT;
    SolverT solver(pc, 500, maxiter, tol);

    // Zeitdiskretisierung mit one-step-theta-scheme
    // theta=1 -> impl. Euler
    // theta=0.5 -> Crank-Nicholson
    InstatPoissonThetaSchemeCL<MyPoissonCL, SolverT>
      ThetaScheme(Poisson, solver, theta, true, Flag_SUPG==1);  //first bool(convection), second bool(stabilization)
    ThetaScheme.SetTimeStep(C.dt_);

    int MaxIter= 0;
    double MaxRes= 0., average= 0.;
    Poisson.Init(x, PoissonCoeffCL::Initial, 0);

    for (int step=1;step<=C.nt_;step++)
      {
	ThetaScheme.DoStep(x);
	average+= solver.GetIter();
	if (MaxIter<=solver.GetIter())
	  MaxIter= solver.GetIter();
	if (MaxRes<=solver.GetResid())
	  MaxRes= solver.GetResid();
	MC.SetSol3D(Poisson.GetSolution(), Poisson.t);
      }
    average/= C.nt_;
    mexPrintf("Num. Iterations in average: %g\n", average);
    mexPrintf("max. Num. Iterations: %d\n", MaxIter);
    mexPrintf("max. res. Norm: %g\n", MaxRes);
    mexPrintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    A.Reset(); M.Reset(); U.Reset();
    b.Reset(); U_SD.Reset(); M_SD.Reset();
    vU_SD.Reset(); vM_SD.Reset();
  }

} // end of namespace DROPS

// Create Geometry, boundary condition, then call strategy
static void source_dp( double theta, double tol, int iter, int Flag_pr, int Flag_bc, int Flag_SUPG)
{ // create Problem (Multigrid, BndData and Coeff) and call Strategy(...)
  try
    {
      DROPS::Point3DCL null(0.0);
      DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
      e1[0]= C.lx_;
      e2[1]= C.ly_;
      e3[2]= C.lz_;

      typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL>
	InstatPoissonOnBrickCL;
      typedef InstatPoissonOnBrickCL MyPoissonCL;

      DROPS::BrickBuilderCL brick(null, e1, e2, e3, C.nx_, C.ny_, C.nz_);

      mexPrintf("Delta t = %g\n", C.dt_);
      mexPrintf("Num. of timesteps = %d\n", C.nt_);

      const bool InterfaceBCType= !(Flag_bc & DirBCFlagC);
      const bool isneumann[6]=
	{ false, true,              // Gamma_in, Gamma_out
	  true,  InterfaceBCType,   // Gamma_h (wall), Gamma_r (surface)
	  true,  true };            // Gamma_r, Gamma_r
      //Interface function creation
      const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
	{ &Inflow, &Zero, &Zero, InterfaceBCType ? &InterfaceFlux : &InterfaceValue, &Zero, &Zero};

      DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
      MyPoissonCL prob(brick, PoissonCoeffCL(), bdata, Flag_pr & AdjFlagC);           //Adjoint problem
      DROPS::MultiGridCL& mg = prob.GetMG();

      //prepare MC
      MC.SetMG(&mg);
      MatlabConnectCL::ClearMaps();
      MatlabConnectCL::setFaceMap();
      MatlabConnectCL::setTetraMap();

      DROPS::Strategy(prob, theta, tol, iter, Flag_SUPG);

      return;
    }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}

//mexFunction gets parameters;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double Tol, theta;
  int Iter, Flag_pr, Flag_bc, Flag_SUPG;

  // Check for proper number of arguments.
  if(nrhs!=21)
    mexErrMsgTxt("(C0, b_in, F, alpha, b_interface, uN, amol, xl, yl, zl, nx, ny, nz, dt, nt, Theta, Tol, Iter, Flag_pr, Flag_bc, Flag_SUPG) as input required.");
  if(nlhs!=2)
    mexErrMsgTxt("Concentration field (3D instat.) and maximum number of GMRES-iterations as output required.");

  // Check to make sure the last input arguments are scalar.
  for (int index=5; index<nrhs; index++)
    if(!mxIsDouble(prhs[index]) ||
       mxGetN(prhs[index])*mxGetM(prhs[index])!=1)
      {
    	mexPrintf("Error in %dth input argument!\n", index+1);
        mexErrMsgTxt("Input must be a scalar.");
      }

  // Get the scalar input arguments.
  C.uN_ = mxGetScalar(prhs[5]);           //before was vicosity
  C.a_mol_ = mxGetScalar(prhs[6]);
  C.lx_ = mxGetScalar(prhs[7]);
  C.ly_ = mxGetScalar(prhs[8]);
  C.lz_ = mxGetScalar(prhs[9]);
  C.nx_ = rd( mxGetScalar(prhs[10]));
  C.ny_ = rd( mxGetScalar(prhs[11]));
  C.nz_ = rd( mxGetScalar(prhs[12]));
  C.dt_ = mxGetScalar(prhs[13]);
  C.nt_ = rd( mxGetScalar(prhs[14]));
  theta = mxGetScalar(prhs[15]);
  Tol =  mxGetScalar(prhs[16]);
  Iter = rd( mxGetScalar(prhs[17]));
  Flag_pr = rd( mxGetScalar(prhs[18]));
  Flag_bc = rd( mxGetScalar(prhs[19]));
  Flag_SUPG=rd( mxGetScalar(prhs[20]));
  if (Flag_pr) {
    C.uN_ = - C.uN_;
  }

  int Nx   = C.nx_+1;
  int Ny   = C.ny_+1;
  int Nz   = C.nz_+1;
  int Nxz  = Nx*Nz;
  int Nyz  = Ny*Nz;
  int Nxyz = Nx*Ny*Nz;

  // Check the dimensions of the input matrices.
  if (mxGetM(prhs[0]) != Nxyz || mxGetN(prhs[0])!= 1) {
    std::stringstream ss;
    ss << "Initial value matrix has wrong dimensions! Should be " << Nxyz << " x " << "1, but is " << mxGetM(prhs[0]) << " x " << mxGetN(prhs[0]);
    mexErrMsgTxt(ss.str().c_str());
  }
  if (mxGetM(prhs[1]) != Nyz  || mxGetN(prhs[1])!= C.nt_+1)
    mexErrMsgTxt("Input T_in has wrong dimensions.");
  if (mxGetM(prhs[2]) != Nxyz || mxGetN(prhs[2])!= C.nt_+1) {
    std::stringstream ss;
    ss << "F (input #" << 3 << ") has wrong dimensions! Should be " << Nxyz << " x " << C.nt_+1 << ", but is " << mxGetM(prhs[2]) << " x " << mxGetN(prhs[2]);
    mexErrMsgTxt(ss.str().c_str());
  }
  if (mxGetM(prhs[3]) != Nxyz || mxGetN(prhs[3])!= C.nt_+1){
    std::stringstream ss;
    ss << "a (input #" << 4 << ") has wrong dimensions! Should be " << Nxyz << " x " << C.nt_+1 << ", but is " << mxGetM(prhs[3]) << " x " << mxGetN(prhs[3]);
    mexErrMsgTxt(ss.str().c_str());
  }
  if (mxGetM(prhs[4]) != Nxz || mxGetN(prhs[4])!= C.nt_+1){
    std::stringstream ss;
    ss << "Interface-boundary (input #" << 5 << ") has wrong dimensions! Should be " << Nxz << " x " << C.nt_+1 << ", but is " << mxGetM(prhs[4]) << " x " << mxGetN(prhs[4]);
    mexErrMsgTxt(ss.str().c_str());
  }

  // Set the input matrices and output parameters.
  const double* C0    = mxGetPr(prhs[0]);
  const double* B_in  = mxGetPr(prhs[1]);
  const double* F     = mxGetPr(prhs[2]);
  const double* a     = mxGetPr(prhs[3]);
  const double* B_Int = mxGetPr(prhs[4]);

  // Allocate memory for output arguments.
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(Nxyz, C.nt_, mxREAL); // w/o initial time step
  double* maxiter = mxGetPr(plhs[0]);
  double* c_sol   = mxGetPr(plhs[1]);

  MC.Init(C, Flag_bc, C0, B_in, F, a, B_Int, maxiter, c_sol);

  // Call the subroutine.
  source_dp(theta, Tol, Iter, Flag_pr, Flag_bc, Flag_SUPG);

  return;
}
