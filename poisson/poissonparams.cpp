/// \file lsetparams.cpp
/// \brief parameters for two-phase flow problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande; SC RWTH Aachen: Oliver Fortmeier

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

#include "poisson/params.h"

namespace DROPS
{

void ParamPoissonCL::RegisterParams()
{
    rp_.BeginGroup( "Poisson");
	rp_.RegInt(    pos_PcIter,			"PcIter");
	rp_.RegDouble( pos_PcTol,			"PcTol");
	rp_.RegInt(    pos_Iter,			"Iter");
	rp_.RegInt(    pos_Restart,			"Restart");
	rp_.RegInt(    pos_RelativeErr,		"RelativeErr");
	rp_.RegDouble( pos_Tol,				"Tol");
	rp_.RegInt(    pos_Method,			"Method");
	rp_.RegInt(    pos_SmoothingSteps,	"SmoothingSteps");
	rp_.RegInt(    pos_NumLvl,			"NumLvl");
	rp_.RegInt(    pos_SolutionIsKnown, "SolutionIsKnown");
	rp_.RegDouble( pos_Relax,			"Relax");
    rp_.EndGroup();
}

void ParamTimeDiscPoissonCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt(    tm_NumSteps,  "NumSteps");
    rp_.RegDouble( tm_StepSize,  "StepSize");
    rp_.RegInt(    tm_Scheme,    "Scheme");
    rp_.RegDouble( tm_Theta,	 "Theta");
    rp_.RegDouble( tm_Nu,		 "Nu");
    rp_.RegInt( tm_Convection,   "Convection");
    rp_.EndGroup();
}

void ParamExperimentalDataPoissonCL::RegisterParams()
{
    rp_.BeginGroup("Exp");
    rp_.RegDouble( exp_Heat,	   "Heat");
    rp_.RegDouble( exp_Rho,		   "Rho");
    rp_.RegDouble( exp_Mu,		   "Mu");
    rp_.RegDouble( exp_Cp,		   "Cp");
    rp_.RegDouble( exp_Lambda,	   "Lambda");
    rp_.EndGroup();
}

void ParamErrCL::RegisterParams()
{
	rp_.BeginGroup("Err");
	rp_.RegInt( err_DoErrorEstimate, 	"DoErrorEstimate");
	rp_.RegDouble( err_RelReduction,	"RelReduction");
	rp_.RegDouble( err_MinRatio,		"MinRatio");
	rp_.RegDouble( err_Threshold,		"Threshold");
	rp_.RegDouble( err_Meas,			"Meas");
    rp_.RegInt( err_DoMark,				"DoMark");
    rp_.RegInt( err_NumRef,				"NumRef");
	rp_.EndGroup();
}

void ParamPoissonProblemCL::RegisterParams()
{

}

} // end of namespace DROPS
