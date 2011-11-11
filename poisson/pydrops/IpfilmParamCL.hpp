#include "../poisson.h"

 class IpfilmParamCL
{
public:
  //Geometry
  double lx_, ly_, lz_;
  int    nx_, ny_, nz_;
  //Refinement steps
  int    refinement_;
  //Instationary
  int    nt_;
  double dt_;
  //Solver
  int max_iter_;
  double tol_;
  //Convection, SUPG, solution is known
  bool convection_, SUPG_, SolutionIsKnown_, ALE_;
  //Initial function
  instat_scalar_fun_ptr init_;
  //Analytical solution
  instat_scalar_fun_ptr solution_;
  //film surface function
  instat_scalar_fun_ptr film_height_;
  //inlet boundary condition
  instat_scalar_fun_ptr c_in_;
  //interface boundary condition
  instat_scalar_fun_ptr c_surface_;
  //diffusion parameter
  instat_scalar_fun_ptr diffusion_parameter_;
  //velocity field
  instat_vector_fun_ptr velocity_;
  
  IpfilmParamCL()
    : lx_(1.), ly_(1.), lz_(1.), nx_(4), ny_(4), nz_(4), refinement_(1),
      nt_(10), dt_{0.01}, 
      max_iter_(1000), tol_(1.0e-8),
      convection_(false), SUPG_(false), SolutionIsKnown_(false), ALE_(false);
      init_= NULL, solution_ = NULL, film_height_ = NULL, 
      c_in_ = NULL, c_surface_ = NULL, diffusion_parameter_ =NULL,
      velocity_ = NULL 
  {}
  
 void SetupCoeffGeom(double lx, double ly, double lz, int nx, int ny, int nz)
 { lx_ = lx; ly_ = ly; lz_ = lz; nx = nx_; ny = ny_; nz = nz_;}
 
 void SetupRefinement(int refinement) {refinement_ = refinement;}
 
 void SetupCoeffInstat(int nt, double dt) { nt_= nt; dt_=dt;}
 
 void SetupCoeffSolver(int max_iter, double tol){ max_iter_ = max_iter; tol_= tol;}
 
 void Ifconvection(bool convection){ convection_ = convection;}
 
 void IfSUPG(bool SUPG) { SUPG_ = SUPG; }
 
 void SolutionIsKnown(bool SolutionIsKnown = false){SolutionIsKnown_ = SolutionIsKnown;}
 
 void IfALE(bool ALE = false) {ALE_ = ALE;}
 
 void SetupInit(instat_scalar_fun_ptr init){init_ = init;}
 
 void SetupSolution(instat_scalar_fun_ptr solution){solution_ = solution;}
 
 void SetupBoundaryConditions(instat_scalar_fun_ptr c_in, instat_scalar_fun_ptr c_surface) {c_in_ = c_in; c_surface_ = c_surface;}
 
 void SetupPoissonCoeff(instat_scalar_fun_ptr diffusion_parameter, instat_vector_fun_ptr velocity) { diffusion_parameter = diffusion_parameter; velocity_ = velocity;)
} 