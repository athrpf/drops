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

#include "levelset/params.h"

namespace DROPS
{

void ParamEnsightCL::RegisterParams()
{
    rp_.BeginGroup( "Ensight" );
    rp_.RegInt(     ens_EnsightOut, "EnsightOut" );
    rp_.RegString(  ens_EnsCase,    "EnsCase" );
    rp_.RegString(  ens_EnsDir,     "EnsDir" );
    rp_.RegString(  ens_GeomName,   "GeomName");
    rp_.RegInt(     ens_MasterOut,  "MasterOut");
    rp_.RegInt(     ens_Binary,     "Binary");
    rp_.EndGroup();
}

void ParamVTKCL::RegisterParams()
{
    rp_.BeginGroup( "VTK" );
    rp_.RegInt(     vtk_VTKOut,    "VTKOut" );
    rp_.RegString(  vtk_VTKDir,    "VTKDir" );
    rp_.RegString(  vtk_VTKName,   "VTKName");
    rp_.RegInt(     vtk_Binary,    "Binary");
    rp_.EndGroup();
}

void ParamInfoCL::RegisterParams()
{
    rp_.BeginGroup( "Info" );
    rp_.RegInt(     inf_PrintSize,   "PrintSize");
    rp_.RegInt(     inf_PrintNumUnk, "PrintNumUnk");
    rp_.RegInt(     inf_CheckMG,     "CheckMG");
    rp_.EndGroup();
}

void ParamQuadCL::RegisterParams()
{
    rp_.BeginGroup("QuadrilateralGrid");
    rp_.RegInt(   qlg_Quad,         "Quad");
    rp_.RegInt(   qlg_GridX,        "GridX");
    rp_.RegInt(   qlg_GridY,        "GridY");
    rp_.RegInt(   qlg_GridZ,        "GridZ");
    rp_.RegCoord( qlg_Stepsize,     "Stepsize");
    rp_.RegCoord( qlg_Barycenter,   "Barycenter");
    rp_.RegCoord( qlg_Rotation,     "Rotation");
    rp_.RegString(qlg_FileName, "FileName");
    rp_.EndGroup();
}

void ParamMGSerCL::RegisterParams()
{
    rp_.BeginGroup( "Restart");
    rp_.RegString( rst_Outputfile,    "Outputfile");
    rp_.RegString( rst_Inputfile,     "Inputfile");
    rp_.RegInt(    rst_Serialization, "Serialization");
    rp_.RegInt(    rst_Overwrite,     "Overwrite");
    rp_.RegInt(    rst_Binary,        "Binary");
    rp_.EndGroup();
}

void ParamBrickCL::RegisterParams()
{
    rp_.BeginGroup( "Brick");
    rp_.RegInt(   brk_BasicRefX,  "BasicRefX");
    rp_.RegInt(   brk_BasicRefY,  "BasicRefY");
    rp_.RegInt(   brk_BasicRefZ,  "BasicRefZ");
    rp_.RegCoord( brk_orig,       "orig");
    rp_.RegCoord( brk_dim,        "dim");
    rp_.EndGroup();
}

void ParamReparamCL::RegisterParams()
{
    rp_.BeginGroup( "Reparam");
    rp_.RegInt(    rpm_Freq,       "Freq");
    rp_.RegInt(    rpm_Method,     "Method");
    rp_.RegDouble( rpm_NarrowBand, "NarrowBand");
    rp_.RegDouble( rpm_MinGrad,    "MinGrad");
    rp_.RegDouble( rpm_MaxGrad,    "MaxGrad");
    rp_.EndGroup();
}

void ParamAdaptRefCL::RegisterParams()
{
    rp_.BeginGroup("AdaptRef");
    rp_.RegInt(    ref_Freq,           "Freq");
    rp_.RegInt(    ref_FinestLevel,    "FinestLevel");
    rp_.RegInt(    ref_CoarsestLevel,  "CoarsestLevel");
    rp_.RegDouble( ref_Width,          "Width");
    rp_.RegInt(    ref_LoadBalStrategy,"LoadBalStrategy" );
    rp_.RegInt(    ref_Partitioner,    "Partitioner");
    rp_.EndGroup();
}

void ParamStokesCL::RegisterParams()
{
    rp_.BeginGroup("Stokes");
    rp_.RegInt(    stk_StokesMethod, "StokesMethod");
    rp_.RegDouble( stk_InnerTol,     "InnerTol");
    rp_.RegDouble( stk_OuterTol,     "OuterTol");
    rp_.RegInt(    stk_InnerIter,    "InnerIter");
    rp_.RegInt(    stk_OuterIter,    "OuterIter");
    rp_.RegInt(    stk_PcAIter,      "PcAIter");
    rp_.RegDouble( stk_PcATol,       "PcATol");
    rp_.RegDouble( stk_PcSTol,       "PcSTol");
    rp_.RegDouble( stk_XFEMStab,     "XFEMStab");
    rp_.RegDouble( stk_Theta,        "Theta");
    rp_.EndGroup();
}

void ParamNavStokesCL::RegisterParams()
{
    rp_.BeginGroup( "NavStokes");
    rp_.RegDouble( ns_Nonlinear,  "Nonlinear");
    rp_.RegDouble( ns_Tol,        "Tol");
    rp_.RegDouble( ns_Reduction,  "Reduction");
    rp_.RegInt(    ns_Iter,       "Iter");
    rp_.EndGroup();
}

void ParamTimeDiscCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt(    tm_NumSteps,  "NumSteps");
    rp_.RegDouble( tm_StepSize,  "StepSize");
    rp_.RegInt(    tm_Scheme,    "Scheme");
    rp_.EndGroup();
}

void ParamLevelSetCL::RegisterParams()
{
    rp_.BeginGroup("Levelset");
    rp_.RegInt(    lvs_Iter,           "Iter");
    rp_.RegDouble( lvs_Tol,            "Tol");
    rp_.RegDouble( lvs_SD,             "SD");
    rp_.RegDouble( lvs_CurvDiff,       "CurvDiff");
    rp_.RegInt(    lvs_VolCorrection,  "VolCorrection");
    rp_.RegDouble( lvs_Theta,          "Theta");
    rp_.EndGroup();
}

void ParamCouplingCL::RegisterParams()
{
    rp_.BeginGroup("Coupling");
    rp_.RegInt(    cpl_Iter,        "Iter");
    rp_.RegDouble( cpl_Tol,         "Tol");
    rp_.RegDouble( cpl_Stab,        "Stab");
    rp_.RegDouble( cpl_Projection,  "Projection");
    rp_.EndGroup();
}

void ParamMaterialDataCL::RegisterParams()
{
    rp_.BeginGroup("Mat");
    rp_.RegDouble( mat_DensDrop,    "DensDrop");
    rp_.RegDouble( mat_ViscDrop,    "ViscDrop");
    rp_.RegDouble( mat_DensFluid,   "DensFluid");
    rp_.RegDouble( mat_ViscFluid,   "ViscFluid");
    rp_.RegDouble( mat_SmoothZone,  "SmoothZone");
    rp_.EndGroup();
}

void ParamExperimentalDataCL::RegisterParams()
{
    rp_.BeginGroup("Exp");
    rp_.RegCoord(  exp_RadDrop,    "RadDrop");
    rp_.RegCoord(  exp_PosDrop,    "PosDrop");
    rp_.RegCoord(  exp_Gravity,    "Gravity");
    rp_.RegDouble( exp_InflowVel,  "InflowVel");
    rp_.RegDouble( exp_RadInlet,   "RadInlet");
    rp_.RegInt(    exp_FlowDir,    "FlowDir");
    rp_.RegDouble( exp_InflowFreq, "InflowFreq");
    rp_.RegDouble( exp_InflowAmpl, "InflowAmpl");
    rp_.EndGroup();
}

void ParamSurfaceTensionCL::RegisterParams()
{
    rp_.BeginGroup("SurfTens");
    rp_.RegInt(    sft_VarTension,  "VarTension");
    rp_.RegDouble( sft_SurfTension, "SurfTension");
    rp_.RegDouble( sft_JumpWidth,   "JumpWidth");
    rp_.RegDouble( sft_RelPos,      "RelPos");
    rp_.RegDouble( sft_DirtFactor,  "DirtFactor");
    rp_.EndGroup();
}

void ParamTransportCL::RegisterParams()
{
    rp_.BeginGroup("Transp");
    rp_.RegInt    ( trp_DoTransp,        "DoTransp");
    rp_.RegDouble ( trp_Theta,           "Theta");
    rp_.RegInt    ( trp_Iter,            "Iter");
    rp_.RegDouble ( trp_Tol,             "Tol");
    rp_.RegDouble ( trp_DiffPos,         "DiffPos");
    rp_.RegDouble ( trp_DiffNeg,         "DiffNeg");
    rp_.RegDouble ( trp_HPos,            "HPos");
    rp_.RegDouble ( trp_HNeg,            "HNeg");
    rp_.RegDouble ( trp_IniCPos,         "IniCPos");
    rp_.RegDouble ( trp_IniCNeg,         "IniCNeg");
    rp_.RegDouble ( trp_NitschePenalty,  "NitschePenalty");
    rp_.RegDouble ( trp_NitscheXFEMStab, "NitscheXFEMStab");
    rp_.EndGroup();
}

void ParamSurfactantTransportCL::RegisterParams()
{
    rp_.BeginGroup("SurfTransp");
    rp_.RegInt    ( surf_DoTransp, "DoTransp");
    rp_.RegDouble ( surf_Theta,    "Theta");
    rp_.RegInt    ( surf_Iter,     "Iter");
    rp_.RegDouble ( surf_Tol,      "Tol");
    rp_.RegDouble ( surf_OmitBound,"OmitBound");
    rp_.RegDouble ( surf_Visc,     "Visc");
    rp_.EndGroup();
}

void ParamDomainCondCL::RegisterParams()
{
    rp_.BeginGroup("DomainCond");
    rp_.RegInt(    dmc_InitialCond,   "InitialCond");
    rp_.RegString( dmc_InitialFile,   "InitialFile");
    rp_.RegString( dmc_MeshFile,      "MeshFile");
    rp_.RegInt(    dmc_GeomType,      "GeomType");
    rp_.RegInt(    dmc_BoundaryType,  "BoundaryType");
    rp_.EndGroup();
}

void ParamMesszelleNsCL::RegisterParams()
{

}

void ParamFilmCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt(    tm_NumSteps,       "NumSteps");
    rp_.RegDouble( tm_StepSize,       "StepSize");
    rp_.EndGroup();

    rp_.BeginGroup("Stokes");
    rp_.RegInt(    stk_InnerIter,    "InnerIter");
    rp_.RegInt(    stk_OuterIter,    "OuterIter");
    rp_.RegDouble( stk_InnerTol,     "InnerTol");
    rp_.RegDouble( stk_OuterTol,     "OuterTol");
    rp_.RegInt(    stk_StokesMethod, "StokesMethod");
    rp_.RegInt(    stk_PcAIter,      "PcAIter");
    rp_.RegDouble( stk_PcATol,       "PcATol");
    rp_.RegDouble( stk_PcSTol,       "PcSTol");
    rp_.RegDouble( stk_XFEMStab,     "XFEMStab");
    rp_.RegDouble( stk_Theta,        "Theta");
    rp_.EndGroup();

    rp_.BeginGroup( "NavStokes");
    rp_.RegDouble( ns_Nonlinear,  "Nonlinear");
    rp_.RegDouble( ns_Tol,        "Tol");
    rp_.RegDouble( ns_Reduction,  "Reduction");
    rp_.RegInt(    ns_Iter,       "Iter");
    rp_.EndGroup();

    rp_.BeginGroup("Levelset");
    rp_.RegInt(    lvs_Iter,           "Iter");
    rp_.RegDouble( lvs_Tol,            "Tol");
    rp_.RegDouble( lvs_SD,             "SD");
    rp_.RegDouble( lvs_CurvDiff,       "CurvDiff");
    rp_.RegInt(    lvs_VolCorrection,  "VolCorrection");
    rp_.RegDouble( lvs_Theta,          "Theta");
    rp_.EndGroup();

    rp_.BeginGroup("Coupling");
    rp_.RegInt(    cpl_Iter,  "Iter");
    rp_.RegDouble( cpl_Tol,   "Tol");
    rp_.EndGroup();

    rp_.BeginGroup("Reparam");
    rp_.RegInt(    rpm_Freq,      "Freq");
    rp_.RegInt(    rpm_Method,    "Method");
    rp_.EndGroup();

    rp_.BeginGroup("AdaptRef");
    rp_.RegInt(    ref_Freq,         "Freq");
    rp_.RegInt(    ref_FinestLevel,  "FinestLevel");
    rp_.RegInt(    ref_CoarsestLevel,"CoarsestLevel");
    rp_.RegDouble( ref_Width,        "Width");
    rp_.RegInt(    ref_LoadBalStrategy,"LoadBalStrategy" );
    rp_.RegInt(    ref_Partitioner,    "Partitioner");
    rp_.EndGroup();

    rp_.BeginGroup("Mat");
    rp_.RegDouble( mat_DensFluid,    "DensFluid");
    rp_.RegDouble( mat_ViscFluid,    "ViscFluid");
    rp_.RegDouble( mat_DensGas,      "DensGas");
    rp_.RegDouble( mat_ViscGas,      "ViscGas");
    rp_.RegDouble( mat_SmoothZone,   "SmoothZone");
    rp_.RegDouble( mat_SurfTension,  "SurfTension");
    rp_.EndGroup();

    rp_.BeginGroup("Exp");
    rp_.RegDouble( exp_Thickness,  "Thickness");
    rp_.RegCoord(  exp_Gravity,    "Gravity");
    rp_.RegDouble( exp_PumpFreq,   "PumpFreq");
    rp_.RegDouble( exp_PumpAmpl,   "PumpAmpl");
    rp_.RegDouble( exp_Ampl_zDir,  "Ampl_zDir");
    rp_.EndGroup();

    // miscellaneous

    rp_.RegInt(    mcl_InitialCond,          "InitialCond");
    rp_.RegString( mcl_EnsightCase,          "EnsightCase");
    rp_.RegString( mcl_EnsightDir,           "EnsightDir");
    rp_.RegString( mcl_InitialFile,          "InitialFile");
    rp_.RegString( mcl_BndCond,              "BndCond");
    rp_.RegCoord(  mcl_MeshSize,             "MeshSize");
    rp_.RegCoord(  mcl_MeshResolution,       "MeshResolution");
    rp_.RegString( mcl_SerializationFile,    "SerializationFile");
    rp_.RegString( mcl_DeserializationFile,  "DeserializationFile");
}

} // end of namespace DROPS
