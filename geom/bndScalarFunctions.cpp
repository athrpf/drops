/// \file bndScalarFunctions.cpp
/// \brief collections of scalar boundary functions
/// \author LNM RWTH Aachen: Martin Horsky; SC RWTH Aachen:
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
#include "misc/container.h"
#include "poisson/params.h"
#include "misc/bndmap.h"

#ifndef BNDSCALARFUNCTIONS_H_
#define BNDSCALARFUNCTIONS_H_

namespace DROPS
{
//========================================================================
//                         General Functions
//========================================================================

/// returning zero
double Zero( const Point3DCL&, double) { return 0.; }

/// boundary description of a neumann problem
// uses constant function f = (-1)^seg *4.0
template<int sel>
double NeuConst( const Point3DCL& , double ){return std::pow(-1,sel)*4.0; }

/// boundary description of a neumann problem
// uses exp-function f =  e^(t)* e^(px + py + pz)
template<int sel>
double NeuExp( const Point3DCL& p, double t){return std::pow(-1,sel)*std::exp(t)*std::exp(p[0]+p[1]+p[2]); }

/// boundary description of a neumann problem
// uses polynomial function
double NeuPoly( const Point3DCL& p, double ){return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]);}

//========================================================================
//                       Functions for poissonP1.cpp
//========================================================================

double Heat(const Point3DCL&, double)
{
	extern ParamPoissonProblemCL C;
	return C.exp_Heat/C.exp_Lambda*1e-3;
}

//========================================================================
//                   Registrierung der Funktionen
//========================================================================
static RegisterScalarFunction regscazero("Zero", Zero);
static RegisterScalarFunction regscaconstpos("NeuConstPos", NeuConst<0>);
static RegisterScalarFunction regscaconstneg("NeuConstNeg", NeuConst<1>);
static RegisterScalarFunction regscaexppos("NeuExpPos", DROPS::NeuExp<0>);
static RegisterScalarFunction regscaexpneg("NeuExpNeg", DROPS::NeuExp<1>);
static RegisterScalarFunction regscapoly("NeuPoly", NeuPoly);
static RegisterScalarFunction regscaheat("Heat", Heat);


}//end namespace DROPS
#endif /* BNDSCALARFUNCTIONS_H_ */
