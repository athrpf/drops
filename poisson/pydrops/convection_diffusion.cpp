#include "poisson/poisson.h"
#include "num/solver.h"
#include "integrTime.h"
#include "out/output.h"
#include <fstream>
#include <string>
#include <sstream>

// Solve the convection-diffusion equation
// du/dt - nu*laplace u + Vel grad u + q*u = f
double Zero(const DROPS::Point3DCL&, double) { return 0.0; }

class PythonConnectCL
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
  static double GetInitial( const DROPS::Point3DCL& p) const
  {
    return C0_[GetNum(p)];
  };
  //boundary functions
  //x=0;
  static double GetInflow( const DROPS::Point3DCL& p, double t) const
  {
    return B_in_[GetNum(p,t,0)];
  };
  //y=0: if neumann condition is active
  static double GetInterfaceFlux( const DROPS::Point3DCL& p, double t) const
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
  static double GetInterfaceValue( const DROPS::Point3DCL& p, double t) const
  {
    return B_Int_[GetNum(p,t,3)];
  };
  //rhs
  static double GetRhs( const DROPS::Point3DCL& p, double t) const
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
  static double GetA( const DROPS::Point3DCL& p, double t) const
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

//the main function
void convection_diffusion(DROPS::ParamCL& param, double* C0, double* b_in, double* b_interface, double* source, double* Diff);
{
    try
    {
        DROPS::Point3DCL null(0.0);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= param.get<double>("Geom.lx");
        e2[1]= param.get<double>("Geom.ly");
        e3[2]= param.get<double>("Geom.lz");
        
        //create geometry
        DROPS::PoissonBndDataCL* bdata = 0;
        //In poissonP1.cpp we use builddomain function
        DROPS::BrickBuilderCL brick(null, e1, e2, e3, param.get<int>("Geom.nx"), param.get<int>("Geom.nx"),param.get<int>("Geom.nx"),);

        const bool isneumann[6]=
        { false, true,              // inlet, outlet
          true,  false,             // wall, interface
          true,  true };            //  
          
        DROPS::PoisssonCoeffCL<DROPS::ParamCL> PoissonCoef(param);
        const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
        { param.c_in_, &Zero, &Zero, param.c_surface_, &Zero, &Zero};

        DROPS::PoissonBndDataCL bdata(6, isneumann, bnd_fun);
        // Setup the problem
        DROPS::PoissonP1CL<DROPS::PoissonCoeffCL<DROPS::ParamCL> > prob( *mg, PoissonCoef, *bdata);
        DROPS::MultiGridCL& mg = prob.GetMG();

        // Refine the grid
        // Create new tetrahedra
        for ( int ref=1; ref <= param.get<int>("Geom.Refinement"); ++ref){
            std::cout << " refine (" << ref << ")\n";
            DROPS::MarkAll( *mg);
            mg->Refine();
        }
        mg->SizeInfo(cout);

        // Solve the problem
        DROPS::Strategy( prob);
        std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
        delete mg;
        delete bdata;
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }

  return;
}
