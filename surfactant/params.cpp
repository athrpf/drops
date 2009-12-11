/// \file
/// \brief parameters for  surfactant PDE on the interface.
/// \author LNM RWTH Aachen: Joerg Grande

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

#include "surfactant/params.h"

namespace DROPS
{

void ParamSurfactantCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt( tm_NumSteps,    "NumSteps");
    rp_.RegDouble( tm_StepSize, "StepSize");
    rp_.EndGroup();

    rp_.BeginGroup("Levelset");
    rp_.RegInt( lvs_Iter,    "Iter");
    rp_.RegDouble( lvs_Tol,  "Tol");
    rp_.RegDouble( lvs_Theta,"Theta");
    rp_.RegDouble( lvs_SD,   "SD");
    rp_.RegDouble( lvs_CurvDiff,  "CurvDiff");
    rp_.RegInt( lvs_VolCorr,      "VolCorrection");
    rp_.EndGroup();

    rp_.BeginGroup("Reparam");
    rp_.RegInt( rpm_Freq,      "Freq");
    rp_.RegInt( rpm_Method,    "Method");
    rp_.EndGroup();

    rp_.BeginGroup("AdaptRef");
    rp_.RegInt( ref_Freq,       "Freq");
    rp_.RegInt( ref_FinestLevel,"FinestLevel");
    rp_.RegDouble( ref_Width,   "Width");
    rp_.EndGroup();

    rp_.BeginGroup("Exp");
    rp_.RegCoord( exp_Radius,     "RadDrop");
    rp_.RegCoord( exp_PosDrop,      "PosDrop");
    rp_.RegCoord( exp_Velocity,   "Velocity");
    rp_.EndGroup();

    rp_.BeginGroup("SurfTransp");
    rp_.RegDouble ( surf_Theta,    "Theta");
    rp_.RegInt    ( surf_Iter,     "Iter");
    rp_.RegDouble ( surf_Tol,      "Tol");
    rp_.RegDouble ( surf_OmitBound,"OmitBound");
    rp_.RegDouble ( surf_Visc,     "Visc");
    rp_.EndGroup();

    rp_.RegInt( TestCase,     "TestCase");
    rp_.RegInt( cdiv,         "InitialDivisions");
    rp_.RegString( EnsCase,   "EnsightCase");
    rp_.RegString( EnsDir,    "EnsightDir");
}

} // end of namespace DROPS



