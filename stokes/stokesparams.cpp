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
    rp_.RegDouble( exp_InflowVel,  "InflowVel",  0.0);
    rp_.RegDouble( exp_RadInlet,   "RadInlet",   1.0);
    rp_.RegInt(    exp_FlowDir,    "FlowDir",    1);
    rp_.RegDouble( exp_InflowFreq, "InflowFreq", 0.0);
    rp_.RegDouble( exp_InflowAmpl, "InflowAmpl", 0.0);
    rp_.EndGroup();
}

void ParamMiscCL::RegisterParams()
{
	rp_.BeginGroup("Misc");
	rp_.RegDouble( misc_Omega,    "Omega");
	rp_.RegDouble( misc_Tau,      "Tau");
	rp_.RegInt( misc_ModifyGrid,  "ModifyGrid", 0);
	rp_.RegDouble( misc_MarkLower,"MarkLower",  0.0);
	rp_.EndGroup();
}

void ParamStokesCoeffCL::RegisterParams()
{
	rp_.BeginGroup("StokesCoeff");
	rp_.RegString( stc_Reaction,	"Reaction");
	rp_.RegString( stc_Source,	    "Source");
	rp_.RegString( stc_Solution_Vel,    "Solution_Vel");
	rp_.RegString( stc_Solution_DVel,    "Solution_DVel");
	rp_.RegString( stc_Solution_Pr, "Solution_Pr");
	rp_.EndGroup();
}


void ParamStokesProblemCL::RegisterParams()
{

}

} // end of namespace DROPS
