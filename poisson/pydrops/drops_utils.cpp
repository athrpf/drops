
#include "drops_utils.hpp"

#include "functors.hpp"

MassTransferBrick::MassTransferBrick(int nx, int ny, int nz,
				     double lx, double ly, double lz,
				     DROPS::PoissonBndDataCL::bnd_val_fun inlet_concentration,
				     DROPS::PoissonBndDataCL::bnd_val_fun interface_concentration)
  :
  nx_(nx), ny_(ny), nz_(nz),
  lx_(lx), ly_(ly), lz_(lz),
  brick_(NULL), bdata_(NULL)
{
  using namespace DROPS;
  Point3DCL null(0.0);
  Point3DCL e1(0.0), e2(0.0), e3(0.0);
  e1[0]= lx;
  e2[1]= ly;
  e3[2]= lz;
  brick_ = new BrickBuilderCL(null, e1, e2, e3, nx, ny, nz);

  const bool isneumann[6]=
    { false, true,     // Gamma_in, Gamma_out
      true, false,     // Gamma_h (wall), Gamma_r (surface)
      true, true };    // Gamma_r, Gamma_r

  Zero zero;
  const DROPS::PoissonBndDataCL::bnd_val_fun bnd_fun[6]=
	{inlet_concentration, zero, zero, interface_concentration, zero, zero};

  bdata_ = new DROPS::PoissonBndDataCL(6, isneumann, bnd_fun);
}

MassTransferBrick::~MassTransferBrick()
{
  delete brick_;
  delete bdata_;
}
