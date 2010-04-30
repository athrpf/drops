/// \file stokesparams.cpp
/// \brief parameters for two-phase flow problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande,Eva Loch, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#include "stokes/params.h"

namespace DROPS
{
void ParamStokesMaterialDataCL::RegisterParams()
{
    rp_.BeginGroup("Mat");
    rp_.RegDouble( mat_Dens,    "Dens");
    rp_.RegDouble( mat_Visc,    "Visc");
    rp_.EndGroup();
}

void ParamStokesExperimentalDataCL::RegisterParams()
{
    rp_.BeginGroup("Exp");
    rp_.RegCoord(  exp_Gravity,    "Gravity");
    rp_.RegDouble( exp_InflowVel,  "InflowVel");
    rp_.RegDouble( exp_RadInlet,   "RadInlet");
    rp_.RegInt(    exp_FlowDir,    "FlowDir");
    rp_.RegDouble( exp_InflowFreq, "InflowFreq");
    rp_.RegDouble( exp_InflowAmpl, "InflowAmpl");
    rp_.EndGroup();
}

#ifndef DROPS_WIN
void ParamErrCL::RegisterParams()
{
	rp_.BeginGroup("Err");
	rp_.RegInt( err_DoErrorEstimate,    "DoErrorEstimate");
	rp_.RegDouble( err_RelReduction,    "RelReduction");
	rp_.RegDouble( err_MinRatio,        "MinRatio");
	rp_.RegDouble( err_Threshold,       "Threshold");
	rp_.RegDouble( err_Meas,            "Meas");
    rp_.RegInt( err_DoMark,             "DoMark");
    rp_.RegInt( err_NumRef,             "NumRef");
	rp_.EndGroup();
}
#endif

void ParamMiscCL::RegisterParams()
{
	rp_.BeginGroup("Misc");
	rp_.RegDouble( misc_Omega,    "Omega");
	rp_.RegDouble( misc_Tau,      "Tau");
	rp_.RegInt( misc_ModifyGrid,  "ModifyGrid");
	rp_.EndGroup();
}


void ParamStokesProblemCL::RegisterParams()
{

}

} // end of namespace DROPS
