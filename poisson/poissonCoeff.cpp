/// \file poissonCoeff.cpp
/// \brief boundary and source functions for the poisson-type problems
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "misc/bndmap.h"
#include "poisson/params.h"

//========================================================================
//                       Functions for poissonP1.cpp
//========================================================================

double Heat(const DROPS::Point3DCL&, double)
{
	extern DROPS::ParamPoissonProblemCL C;
	return C.exp_Heat/C.exp_Lambda*1e-3;
}

static DROPS::RegisterScalarFunction regscaheat("Heat", Heat);

