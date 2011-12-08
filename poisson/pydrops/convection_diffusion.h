#ifndef __CONVECTION_DIFFUSION_H__
#define __CONVECTION_DIFFUSION_H__

#include "misc/params.h"

void convection_diffusion(DROPS::ParamCL& P, const double* C0, const double* b_in, const double* b_interface, const double* source, const double* Dw, double* C_sol);

#endif
