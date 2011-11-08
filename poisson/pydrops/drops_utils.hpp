
#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "geom/builder.h"
#include "poisson/poisson.h"
#include "boost/function.hpp"

const int AdjFlagC=32, DirBCFlagC= 16;

/** Rounding
 *  Used for rounding to full integers, used for finding the index corresponding
 *  to a double in an array
 */
inline int rd( double d) { return static_cast<int>( d+0.5); }


//Used to setup computational domain and setup boundary conditions.
class MassTransferBrick
{
private:
  int nx_, ny_, nz_;
  double lx_, ly_, lz_;
  DROPS::BrickBuilderCL* brick_;
  DROPS::PoissonBndDataCL* bdata_;

public:
  MassTransferBrick(int nx, int ny, int nz,
		    double lx, double ly, double lz,
		    DROPS::PoissonBndDataCL::bnd_val_fun inlet_concentration,
		    DROPS::PoissonBndDataCL::bnd_val_fun interface_concentration);
  ~MassTransferBrick();

  DROPS::BrickBuilderCL& get_brick() {return *brick_;}
  DROPS::PoissonBndDataCL& get_bdata() {return *bdata_;}
};

#endif
