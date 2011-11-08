//********************************************************************************
// File:    coeff_dp_stat.cpp                                                     *
// Content: poisson + interface for matlab *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers, Liang Zhang*
//          IGPM RWTH Aachen
//			Maka Karalashvili                                                     *
//          AVT.PT
// Version: 0.1                                                                   *
// History: begin - Jul, 01 2008 modified 07.2011                                 *
//********************************************************************************

#include "mex.h"

#include "poisson/poisson.h"
#include "num/solver.h"
#include "poisson/integrTime.h"
#include "out/output.h"
#include <fstream>

const int SensFlagC= 64, AdjFlagC=32, DirBCFlagC= 16;

/**
This code solves a steady-state equation(for fixed time moments) on a brick-shaped domain with certain bnd conditions
modeling the step 2 of (effective) mass transfer in a flat film , and is extended by a matlab interface. Bnd conditions:
- Dirichlet bc at x=0 (inflow)
- Dirichlet bc at y=y_top (Henry)
- all other bnds: homogeneous natural bc

The program can be called as a mex-function from matlab with the following syntax:

  [MaxIter, Csol] = coeff_dp_stat(B_in, a, F, B_interface, xl, yl, zl, nx, ny, nz, Tol, Iter, Flag_pr, Flag_bc, a_tilde, sol_dp)

scalar parameters:
------------------
  xl, yl, zl:     length of brick in x-/y-/z-direction [in mm]
  nx, ny, nz:     number of intervals in x-/y-/z-direction, i.e. Nxyz = (nx+1) x (ny+1) x (nz+1) grid points
  dt, nt:         length and number of time steps
  Tol, Iter:      stopping criterion for the iterative solver (GMRES)
  a:              the wavy transport coeff. (here dummy)
  B_interface:    bnd condition for Henry condition, which is zero in the second step of ia
  Flag_pr:        used for the sensitivity problem if Flag & SensFlagC, then the sensitivity operator is discretized(a_tilde, sol_dp parameters come in play).
  Flag_bc:		  used for the type of bc.

  MaxIter (output):        number of iterations spent in the iterative solver (maximum over all time steps)

matrix parameters:
------------------
  C_in:   Nyz  x Nt        temperature at inflow (x=0): T_in(y,z,t_i)
  F:      Nxyz x Nt        rhs term F(x,y,z,t_i)
  a_tilde:     Nxyz x Nt        the \tilde{a} in Sensitivity Operator
  sol_dp:      Nxyz x Nt        solution of the direct problem

  Csol (output):   Nxyz x nt        solution of the conv-diff problem,
                                    concentration distribution C(x,y,z,t_i) omitting initial time step

Nyz = (ny+1) x (nz+1),    Nxyz = (nx+1) x Nyz.

**/

extern void _main();

inline int rd( double d) { return static_cast<int>(d+0.5); } // rounding
//
class ParamCL
{
  public:
    double lx_, ly_, lz_;
    int    nx_, ny_, nz_; // n=number of intervalls
    double hflux_, rho_, nu_, cp_, a_mol_;
    std::string EnsDir_, EnsCase_;

    ParamCL()
      : lx_(180), ly_(0.8), lz_(0.3), nx_(8), ny_(2), nz_(2),  // in mm
        hflux_(350), rho_(912), nu_(4.7e-6), cp_(1540), a_mol_(0.118/1540/912),// in SI
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
    int Nx_, Ny_, Nz_, Nxy_, Nxz_, Nyz_, Nxyz_; // N=number of points
    double dx_, dy_, dz_;

    double *C_in_, *a_, *F_, *B_int_;   // input  matrices: inflow temp (Nyz x nt_i),
    							// a(Nxyz nt_i), rhs (Nxyz nt_i), qh
    double *a_tilde_, *sol_dp_; //for sensitivity operator
    double *C3D_, *MaxIter_;   // output matrices: temp solution (Nxyz x nt_i),
                             //     max. iterations of solver (1 x 1)
    //helper maps for barycenters
    static DROPS::MultiGridCL* MG_;

    static FACE_MAP face_map_;
    static TETRA_MAP tetra_map_;

    typedef FACE_MAP::const_iterator  fi;
    typedef TETRA_MAP::const_iterator ti;

    static double rnd(double d) {// rounding four digits after comma
        int i_d = (int)(d*10000);
        double d_d = (double)i_d/10000;
        return d_d;
     }
    int GetNum( const DROPS::Point3DCL& p,int seg) const
    {
    	if (seg == 0 || seg == 1){//yoz
		   	 return (rd(p[2]/dz_)*Ny_ + rd(p[1]/dy_));
    	}
    	if (seg == 2 || seg == 3){//xoz
 	   	     return (rd(p[2]/dz_)*Nx_ + rd(p[0]/dx_));
    	}
    	if (seg == 4 || seg == 5){//xoy
		   	 return (rd(p[1]/dy_)*Nx_ + rd(p[0]/dx_));
    	}
    }

    int GetNum( const DROPS::Point3DCL& p) const
    {
        return rd(p[2]/dz_)*Nxy_ + rd(p[1]/dy_)*Nx_ + rd(p[0]/dx_);
    }

  public:
	MatlabConnectCL()
	{
		Nx_=Ny_=Nz_=Nxy_=Nyz_=Nxz_=Nxyz_=-1;
		dx_=dy_=dz_=0.0;
        //
        face_map_.clear();
        tetra_map_.clear();

	}
    void SetMG( DROPS::MultiGridCL* MG) { MG_= MG; }
    //
    static void ClearMaps() {
        face_map_.clear();
        tetra_map_.clear();
    }
    static void setFaceMap()
    {
//      mexPrintf("SETTING FACE MAP\n");
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
//      mexPrintf("FACE MAP SET\n");
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
        fprintf(f," %d",face_map_.size());
      fclose(f);
      //mexPrintf("END DUMP FACE MAP\n");
    }
    //
    static void setTetraMap()
    {
      //mexPrintf("SETTING TETRA MAP\n");
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
      //mexPrintf("TETRA MAP SET\n");
    }
    static void DumpTetraMap()
    {
      FILE* f=fopen("tetra.map","w");
      //mexPrintf("DUMPING TETRA MAP\n");
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
      //mexPrintf("END DUMP TETRA MAP\n");
    }

	//boundary functions
	//x=0
    double GetInflow( const DROPS::Point3DCL& p) const
    {
      return C_in_[GetNum(p,0)];
    };
    //y=0: is neumann condition is active
    double GetInterfaceFlux( const DROPS::Point3DCL& p) const
    {
    	double ret;
    	int seg=3;

        d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
        cmp_key key= std::make_pair(rnd(p[0]), pr);

        DROPS::FaceCL* face= face_map_[key];

        if (face == NULL) {//non-barycenter
			ret= B_int_[GetNum(p,seg)];
	    } else {
	      	ret= 1./3.*(B_int_[GetNum(face->GetVertex(0)->GetCoord(),seg)]+B_int_[GetNum(face->GetVertex(1)->GetCoord(),seg)]+B_int_[GetNum(face->GetVertex(2)->GetCoord(),seg)]);
	    }

		return ret;
    };
	//y=0: if dirichlet condition is active
    double GetInterfaceValue( const DROPS::Point3DCL& p) const
    {
      return B_int_[GetNum(p,3)];
    };
    //coefficient functions
    double GetA( const DROPS::Point3DCL& p) const
    {
        double ret;

        d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
        cmp_key key= std::make_pair(rnd(p[0]), pr);
        DROPS::TetraCL* tetra= tetra_map_[key];

        if (tetra == NULL) {//non-barycenter
            ret=a_[GetNum(p)];
        }else {
            ret = 0.25*(a_[GetNum(tetra->GetVertex(0)->GetCoord())]+a_[GetNum(tetra->GetVertex(1)->GetCoord())]+
                        a_[GetNum(tetra->GetVertex(2)->GetCoord())]+a_[GetNum(tetra->GetVertex(3)->GetCoord())]) ;
        }

      return ret;
    };
    //rhs
    double GetRhs( const DROPS::Point3DCL& p) const
    {
        double ret;

        d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
        cmp_key key= std::make_pair(rnd(p[0]), pr);
        DROPS::TetraCL* tetra= tetra_map_[key];

        if (tetra == NULL) {//non-barycenter
            ret=F_[GetNum(p)];
        }else {
            ret = 0.25*(F_[GetNum(tetra->GetVertex(0)->GetCoord())]+F_[GetNum(tetra->GetVertex(1)->GetCoord())]+
                        F_[GetNum(tetra->GetVertex(2)->GetCoord())]+F_[GetNum(tetra->GetVertex(3)->GetCoord())]) ;
        }

        return ret;
    };
    //
    double GetATilde( const DROPS::Point3DCL& p) const
    {
      	if (a_tilde_ != NULL) {
            double ret;

            d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
            cmp_key key= std::make_pair(rnd(p[0]), pr);
            DROPS::TetraCL* tetra= tetra_map_[key];

            if (tetra == NULL) {//non-barycenter
                ret=a_tilde_[GetNum(p)];
            }else {
                ret = 0.25*(a_tilde_[GetNum(tetra->GetVertex(0)->GetCoord())]+a_tilde_[GetNum(tetra->GetVertex(1)->GetCoord())]+
                            a_tilde_[GetNum(tetra->GetVertex(2)->GetCoord())]+a_tilde_[GetNum(tetra->GetVertex(3)->GetCoord())]) ;
            }

            return ret;
      	}
    };
    double GetSolDP( const DROPS::Point3DCL& p) const
    {
      	if (sol_dp_!= NULL) {
            double ret;

            d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
            cmp_key key= std::make_pair(rnd(p[0]), pr);
            DROPS::TetraCL* tetra= tetra_map_[key];

            if (tetra == NULL) {//non-barycenter
                ret=sol_dp_[GetNum(p)];
            }else {
                ret = 0.25*(sol_dp_[GetNum(tetra->GetVertex(0)->GetCoord())]+sol_dp_[GetNum(tetra->GetVertex(1)->GetCoord())]+
                            sol_dp_[GetNum(tetra->GetVertex(2)->GetCoord())]+sol_dp_[GetNum(tetra->GetVertex(3)->GetCoord())]) ;
            }

            return ret;
      	}
    };
    //
  	template<class P1EvalT>
  	void SetSol3D( const P1EvalT& sol)
  	{
  		double *out= C3D_;

  		DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
  		{
  			out[GetNum( sit->GetCoord())]= sol.val( *sit);
  		}
  	}
    //initialisation
  	void Init( const ParamCL& P, int Flag_pr, int Flag_bc, mxArray *plhs[], const mxArray *prhs[])
  	{
    	Nx_= P.nx_+1;Ny_= P.ny_+1; Nz_=P.nz_+1;
    	Nyz_=Ny_*Nz_; Nxy_=Nx_*Ny_; Nxz_=Nx_*Nz_;
    	Nxyz_= Nxy_*Nz_;
    	dx_= P.lx_/P.nx_; dy_= P.ly_/P.ny_; dz_= P.lz_/P.nz_;

		// Check to make sure the first input arguments are double matrices.
	    if (mxGetPi(prhs[0])!=NULL)
            mexErrMsgTxt("Input C_in must be a double matrix.");
	    if (mxGetPi(prhs[1])!=NULL)
	  	    mexErrMsgTxt("Input   must be a double matrix.");
	    if (mxGetPi(prhs[2])!=NULL)
	  	    mexErrMsgTxt("Input F must be a double matrix.");
		if (mxGetPi(prhs[3])!=NULL)
	  	    mexErrMsgTxt("Input interface_bc  must be a double matrix.");

	  	if (SensFlagC & Flag_pr)
	  	{
		    if(mxGetPi(prhs[14])!=NULL)
		  	    mexErrMsgTxt("Input _tilde must be a double matrix.");
		    if(mxGetPi(prhs[15])!=NULL)
		  	    mexErrMsgTxt("Input sol_dp must be a double matrix.");
	  	}

	 	// Check the dimensions of the input matrices.
	 	if (mxGetM(prhs[0]) != Nyz_  || mxGetN(prhs[0])!= 1)
	 	    mexErrMsgTxt("Input C_in has wrong dimensions.");
	 	if (mxGetM(prhs[1]) != Nxyz_ || mxGetN(prhs[1])!= 1)
	 	    mexErrMsgTxt("Input  has wrong dimensions.");
	 	if (mxGetM(prhs[2]) != Nxyz_ || mxGetN(prhs[2])!= 1)
	 	    mexErrMsgTxt("Input F has wrong dimensions.");
	 	if (mxGetM(prhs[3]) != Nxz_ || mxGetN(prhs[3])!= 1)
	 	    mexErrMsgTxt("Input interface_bc has wrong dimensions.");

	  	if (SensFlagC & Flag_pr)
	  	{
		 	if (mxGetM(prhs[14]) != Nxyz_ || mxGetN(prhs[14])!= 1)
		 	    mexErrMsgTxt("Input a_tilde has wrong dimensions.");
		 	if (mxGetM(prhs[15]) != Nxyz_ || mxGetN(prhs[15])!= 1)
		 	    mexErrMsgTxt("Input sol_dp has wrong dimensions.");
	  	}

	 	// Get the matrix input arguments.
	 	C_in_=  mxGetPr(prhs[0]);
	 	a_= mxGetPr(prhs[1]);
	 	F_=     mxGetPr(prhs[2]);
	 	if  (!(Flag_bc & DirHeatBCFlagC)) {
			//mexPrintf("solving with Neumann Cond.\n");
	 		hflux= mxGetPr(prhs[3]);
	 		Twall= NULL;
	 	} else {
			//mexPrintf("solving with Dirichlet Cond.\n");
	 		hflux= NULL;
	 		Twall= mxGetPr(prhs[3]);
	 	}
	 	//
	  	if (SensFlagC & Flag_pr)
	  	{
			//mexPrintf("solving SP\n");
		 	a_tilde_= mxGetPr(prhs[14]);
		 	sol_dp_= mxGetPr(prhs[15]);
	  	} else
	  	{
	  		a_tilde_ = NULL;
	  		sol_dp_ = NULL;
	  	}
	 	// Allocate memory for output arguments.
	 	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	 	plhs[1] = mxCreateDoubleMatrix(Nxyz_, 1, mxREAL); // w/o initial time step

	 	// Set the output pointer to the output arguments.
	 	MaxIter_ = mxGetPr(plhs[0]);
	 	C3D_     = mxGetPr(plhs[1]);
  	}
} MC;

DROPS::MultiGridCL* MatlabConnectCL::MG_= NULL;
MatlabConnectCL::FACE_MAP MatlabConnectCL::face_map_;
MatlabConnectCL::TETRA_MAP MatlabConnectCL::tetra_map_;


// Solve the steady-state poiss. equation
// -div(a*grad u )  = f

class PoissonCoeffCL
{
  public:
    static double alpha(const DROPS::Point3DCL& p, double)
    {
      	return MC.GetA(p);
    }
    static double f(const DROPS::Point3DCL& p, double)
    {
    	return MC.GetRhs(p);
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL&, double)
    {
        return DROPS::Point3DCL(0.);
    }
   //Only used for Stabilization for flat film case
    static double h_Value()
    {//mesh size in flow direction
        double h=C.lx_/C.nx_;
        return h;
    }
    static double PecNum(const DROPS::Point3DCL& p, double t)
    { //Peclet Number
        double Pec=0.;
        Pec=Vel(p, t)[0]*h_Value()/(2.*alpha(p,t));
        return Pec;
    }
    static double Sta_Coeff(const DROPS::Point3DCL& p, double t)
    {//Stabilization coefficient
        if (PecNum(p,t)<=1)
          return 0.0;
        else
          return h_Value()/(2.*Vel(p, t)[0])*(1.-1./PecNum(p, t));
    }
};
class GradSrcCL //for sensitivity operator
{
  public:
	static double a_tilde(const DROPS::Point3DCL& p, double)
	{
		return MC.GetATilde(p);
	}
	static double sol_dp(const DROPS::Point3DCL& p, double)
	{
		return MC.GetSolDP(p);
	}
};
//
double Zero(const DROPS::Point3DCL&, double) { return 0.0; }
double Inflow(const DROPS::Point3DCL& p, double) { return MC.GetInflow(p); }
double InterfaceFlux(const DROPS::Point3DCL& p, double) { return MC.GetInterfaceFlux(p); }
double InterfaceValue(const DROPS::Point3DCL& p, double) { return MC.GetInterfaceValue(p); }

namespace DROPS
{

template<class Coeff>
void Strategy(InstatPoissonP1CL<Coeff>& Poisson, int Flag_pr, int Flag_bc,
  			  double tol, int maxiter)
{
  typedef InstatPoissonP1CL<Coeff> MyPoissonCL;

  MultiGridCL& MG= Poisson.GetMG();
  IdxDescCL& idx= Poisson.idx;
  VecDescCL& x= Poisson.x;
  VecDescCL& b= Poisson.b;
  MatDescCL& A= Poisson.A;
  MatDescCL& M= Poisson.M;

  VecDescCL cplA, dummy;

  idx.Set(1, 0, 0, 0);
  // vreate Numbering for this Index
  Poisson.CreateNumbering(MG.GetLastLevel(), &idx);
  // Apply numbering to the vectors
  x.SetIdx( &idx);
  b.SetIdx( &idx);
  cplA.SetIdx( &idx);
  dummy.SetIdx( &idx);
  // Apply numbering to the matrices - rows&columns
  A.SetIdx(&idx, &idx);
  M.SetIdx(&idx, &idx);
  //
  mexPrintf("Number of Unknowns: %d\n", x.Data.size());
  mexPrintf("Tolerance: %g\t", tol);
  mexPrintf("max. Num. of Iterations: %d\n", maxiter);
  //
  Poisson.SetupInstatSystem(A, M, Poisson.t);
  if (Flag_pr & SensFlagC)
  {
	  mexPrintf("setting GradSrc for SP\n");
	  Poisson.SetupInstatRhs( cplA, dummy, 0, dummy, 0);
	  //NOTE: this is obsolette here as hflux=0 in the second step of
	  //the incremental identification
	  if (Flag_bc & DirBCFlagC) {
		  Poisson.SetupGradSrc(b,&GradSrcCL::sol_dp,&GradSrcCL::a_tilde, NULL);
	  }
	  else {
		  Poisson.SetupGradSrc(b,&GradSrcCL::sol_dp,&GradSrcCL::a_tilde, &InterfaceFlux);
	  }
  } else
  {
  	Poisson.SetupInstatRhs( cplA, dummy, 0, b, 0);
  }
  b.Data+= cplA.Data;

  //initialise solver
  SSORPcCL pc(1.0);
  PCG_SsorCL solver(pc, maxiter, tol);
  //solve
  solver.Solve( A.Data, x.Data, b.Data);
  MC.SetSol3D(Poisson.GetSolution());

  mexPrintf("Num. of Iterations max.: %d\n", solver.GetIter());
  mexPrintf("max. Residual Norm: %g\n", solver.GetResid());

  A.Reset(); M.Reset();
  b.Reset();
}

} // end of namespace DROPS


static void coeff_ip_stat( double tol, int iter, int Flag_bc, int Flag_pr)
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
    //setup boundary
    const bool InterfaceBCType= !(Flag_bc & DirBCFlagC);
    if (!InterfaceBCType) {
    	mexPrintf("Solving with Dirichlet boundary");
    else
    	mexPrintf("Solving with Neumann boundary");
    }
    //
    const bool isneumann[6]=
      { false, true,   // Gamma_in, Gamma_out
        true, InterfaceBCType,     // Gamma_h (wall), Gamma_r (surface)
        true, true };  // Gamma_r, Gamma_r
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
	{ &Inflow, &Zero, &Zero, InterfaceBCType ? &InterfaceFlux : &InterfaceValue, &Zero, &Zero};

    DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();

    //prepare MC
    MC.SetMG(&mg);
    MatlabConnectCL::ClearMaps();
    MatlabConnectCL::setFaceMap();
    MatlabConnectCL::setTetraMap();

    DROPS::Strategy(prob, Flag_pr, Flag_bc, tol, iter);

    return;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}
//the mexFunction
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double Tol;
  int Iter, Flag_pr, Flag_bc;

  // Check for proper number of arguments.
  if(nrhs>16)
    mexErrMsgTxt("wrong number of inputs");
  if(nlhs!=2)
    mexErrMsgTxt("Concentration field (3D instat.) and maximum number of iterations as output required.");

  int chkIdx = nrhs;
  if (nrhs == 16) chkIdx-=2;
	// Check to make sure the last input arguments are scalar.
  for (int index=4; index<chkIdx; index++)
    if(!mxIsDouble(prhs[index]) ||
      mxGetN(prhs[index])*mxGetM(prhs[index])!=1)
    {
    	mexPrintf("Error in %dth input argument!\n", index+1);
        mexErrMsgTxt("Input must be a scalar.");
    }


  // Get the scalar input arguments.
  C.lx_ = mxGetScalar(prhs[4]);
  C.ly_ = mxGetScalar(prhs[5]);
  C.lz_ = mxGetScalar(prhs[6]);
  C.nx_ = rd( mxGetScalar(prhs[7]));
  C.ny_ = rd( mxGetScalar(prhs[8]));
  C.nz_ = rd( mxGetScalar(prhs[9]));
  Tol =  mxGetScalar(prhs[10]);
  Iter = rd( mxGetScalar(prhs[11]));
  Flag_pr = rd( mxGetScalar(prhs[12]));
  Flag_bc = rd( mxGetScalar(prhs[13]));

  // Set the input matrices and output parameters.
  MC.Init(C, Flag_pr, Flag_bc, plhs, prhs);

  // Call the subroutine.
  coeff_ip_stat(Tol, Iter, Flag_bc, Flag_pr);

  return;
}
