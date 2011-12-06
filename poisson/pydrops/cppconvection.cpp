#include "misc/params.h"

#include "pyconnect.h"
#include "convection_diffusion.cpp"


int main(int argc, char** argv)
{
  try {
    using namespace DROPS;

    std::ifstream param;
    if (argc!=2){
      std::cout << "Using default parameter file: poissonex1.json\n";
      param.open( "poissonex1.json");
    }
    else
      param.open( argv[1]);
    if (!param){
      std::cerr << "error while opening parameter file\n";
      return 1;
    }
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
    PdeFunction C0(Nx,Ny,Nz,1);
    C0.set_value(1.0);
    PdeFunction b_in(1,Ny,Nz,Nt);
    b_in.set_value(1.0);
    PdeFunction b_interface(Nx,1,Nz,Nt);
    b_interface.set_value(2.0);
    PdeFunction source(Nx,Ny,Nz,Nt);
    source.set_value(0.5);
    PdeFunction Dw(Nx,Ny,Nz,Nt);
    Dw.set_value(0.001);
    convection_diffusion(P, C0.get_data(), b_in.get_data(), b_interface.get_data(), source.get_data(), Dw.get_data(), NULL);

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
