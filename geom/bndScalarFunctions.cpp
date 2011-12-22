/// \file bndScalarFunctions.cpp
/// \brief collections of general scalar functions (like zero, one, etc..). No problem-specific functions!
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
/// time independent function return zero
double tid_Zero( const Point3DCL&) { return 0.; }
/// returning one
double One( const Point3DCL&, double) { return 1.; }
/// time independent function return one
double tid_One( const Point3DCL&) { return 1.; }


//========================================================================
//                   Registrierung der Funktionen
//========================================================================
/*  static RegisterScalarFunction regscazero("Zero", Zero);
  static RegisterScalarFunction regtidscazero("Zero", tid_Zero);
  static RegisterScalarFunction regscaone("One", One);
  static RegisterScalarFunction regtidscaone("One", tid_One);*/
//boost function
  static RegisterScalarFunction regscazero("Zero", instat_scalar_fun_ptr(Zero));
  static RegisterScalarFunction regtidscazero("Zero", scalar_fun_ptr(tid_Zero));
  static RegisterScalarFunction regscaone("One", instat_scalar_fun_ptr(One));
  static RegisterScalarFunction regtidscaone("One", scalar_fun_ptr(tid_One));
}//end namespace DROPS
#endif /* BNDSCALARFUNCTIONS_H_ */
