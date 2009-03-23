/// \file
/// \brief parameters for two-phase flow problems.

#include "levelset/params.h"

namespace DROPS
{

void ParamEnsightCL::RegisterParams()
{
    rp_.BeginGroup( "Ensight" );
    rp_.RegInt(     ensight,   "EnsightOut" );
    rp_.RegString(  EnsCase,   "EnsCase" );
    rp_.RegString(  EnsDir,    "EnsDir" );
    rp_.RegString(  geomName,  "GeomName");
    rp_.RegInt(     masterOut, "MasterOut");
    rp_.RegInt(     binary,    "Binary");
    rp_.EndGroup();
}

void ParamVTKCL::RegisterParams()
{
    rp_.BeginGroup( "VTK" );
    rp_.RegInt(     vtk,       "VTKOut" );
    rp_.RegString(  vtkDir,    "VTKDir" );
    rp_.RegString(  vtkName,   "VTKName");
    rp_.RegInt(     vtkBinary, "Binary");
    rp_.EndGroup();
}

void ParamInfoCL::RegisterParams()
{
    rp_.BeginGroup( "Info" );
    rp_.RegInt(     printSize,   "PrintSize");
    rp_.RegInt(     printNumUnk, "PrintNumUnk");
    rp_.RegInt(     checkMG,     "CheckMG");
    rp_.EndGroup();
}

void ParamQuadCL::RegisterParams()
{
    rp_.BeginGroup("QuadrilateralGrid");
    rp_.RegInt(   quad,         "Quad");
    rp_.RegInt(   gridX,        "GridX");
    rp_.RegInt(   gridY,        "GridY");
    rp_.RegInt(   gridZ,        "GridZ");
    rp_.RegCoord( stepsize,     "Stepsize");
    rp_.RegCoord( barycenter,   "Barycenter");
    rp_.RegCoord( rotation,     "Rotation");
    rp_.RegString(quadFileName, "FileName");
    rp_.EndGroup();
}

void ParamMGSerCL::RegisterParams()
{
    rp_.BeginGroup( "Restart");
    rp_.RegString( ser_dir,               "Outputfile");
    rp_.RegString( deserialization_file,  "Inputfile");
    rp_.RegInt(serialization,             "Serialization");
    rp_.RegInt(overwrite,                 "Overwrite");
    rp_.EndGroup();
}

void ParamBrickCL::RegisterParams()
{
    rp_.BeginGroup( "Brick");
    rp_.RegInt(   basicref_x, "BasicRefX");
    rp_.RegInt(   basicref_y, "BasicRefY");
    rp_.RegInt(   basicref_z, "BasicRefZ");
    rp_.RegCoord( orig,       "orig");
    rp_.RegCoord( dim,        "dim");
    rp_.EndGroup();
}

void ParamReparamCL::RegisterParams()
{
    rp_.BeginGroup( "Reparam");
    rp_.RegInt(    RepFreq,   "Freq");
    rp_.RegInt(    RepMethod, "Method");
    rp_.RegDouble( NarrowBand,"NarrowBand");
    rp_.RegDouble( MinGrad,   "MinGrad");
    rp_.RegDouble( MaxGrad,   "MaxGrad");
    rp_.EndGroup();
}

void ParamAdaptRefCL::RegisterParams()
{
    rp_.BeginGroup("AdaptRef");
    rp_.RegInt(    ref_freq,       "Freq");
    rp_.RegInt(    ref_flevel,     "FinestLevel");
    rp_.RegDouble( ref_width,      "Width");
    rp_.RegInt(    refineStrategy, "RefineStrategy" );
    rp_.EndGroup();
}

void ParamStokesCL::RegisterParams()
{
    rp_.BeginGroup("Stokes");
    rp_.RegInt( StokesMethod, "StokesMethod");
    rp_.RegDouble( inner_tol, "InnerTol");
    rp_.RegDouble( outer_tol, "OuterTol");
    rp_.RegInt( inner_iter,   "InnerIter");
    rp_.RegInt( outer_iter,   "OuterIter");
    rp_.RegInt( pcA_iter,     "PcAIter");
    rp_.RegDouble( pcA_tol,   "PcATol");
    rp_.RegDouble( pcS_tol,   "PcSTol");
    rp_.RegDouble(XFEMStab,   "XFEMStab");
    rp_.EndGroup();
}

void ParamNavStokesCL::RegisterParams()
{
    rp_.BeginGroup( "NavStokes");
    rp_.RegDouble( nonlinear,  "Nonlinear");
    rp_.RegDouble( ns_tol,     "Tol");
    rp_.RegDouble( ns_red,     "Reduction");
    rp_.RegInt( ns_iter,       "Iter");
    rp_.EndGroup();
}

void ParamTimeDiscCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt( num_steps,    "NumSteps");
    rp_.RegDouble( dt,        "StepSize");
    rp_.RegDouble( theta,     "Theta");
    rp_.RegInt(    scheme,    "Scheme");
    rp_.EndGroup();
}

void ParamLevelSetCL::RegisterParams()
{
    rp_.BeginGroup("Levelset");
    rp_.RegInt( lset_iter,    "Iter");
    rp_.RegDouble( lset_tol,  "Tol");
    rp_.RegDouble( lset_SD,   "SD");
    rp_.RegDouble( CurvDiff,  "CurvDiff");
    rp_.RegInt( VolCorr,      "VolCorrection");
    rp_.EndGroup();
}

void ParamCouplingCL::RegisterParams()
{
    rp_.BeginGroup("Coupling");
    rp_.RegInt(    cpl_iter,  "Iter");
    rp_.RegDouble( cpl_tol,   "Tol");
    rp_.RegDouble( cpl_stab,  "Stab");
    rp_.RegDouble( cpl_proj,  "Projection");
    rp_.EndGroup();
}

void ParamMaterialDataCL::RegisterParams()
{
    rp_.BeginGroup("Mat");
    rp_.RegDouble( rhoD,      "DensDrop");
    rp_.RegDouble( muD,       "ViscDrop");
    rp_.RegDouble( rhoF,      "DensFluid");
    rp_.RegDouble( muF,       "ViscFluid");
    rp_.RegDouble( sm_eps,    "SmoothZone");
    rp_.EndGroup();
}

void ParamExperimentalDataCL::RegisterParams()
{
    rp_.BeginGroup("Exp");
    rp_.RegCoord( Radius,     "RadDrop");
    rp_.RegCoord( Mitte,      "PosDrop");
    rp_.RegCoord( g,          "Gravity");
    rp_.RegDouble( Anstroem,  "InflowVel");
    rp_.RegDouble( r_inlet,   "RadInlet");
    rp_.RegInt( flow_dir,     "FlowDir");
    rp_.RegDouble( inflow_freq, "InflowFreq");
    rp_.RegDouble( inflow_ampl, "InflowAmpl");
    rp_.EndGroup();
}

void ParamSurfaceTensionCL::RegisterParams()
{
    rp_.BeginGroup("SurfTens");
    rp_.RegInt( st_var,       "VarTension");
    rp_.RegDouble( sigma,     "SurfTension");
    rp_.RegDouble( st_jumpWidth, "JumpWidth");
    rp_.RegDouble( st_relPos, "RelPos");
    rp_.RegDouble( st_red,    "DirtFactor");
    rp_.EndGroup();
}

void ParamTransportCL::RegisterParams()
{
    rp_.BeginGroup("Transp");
    rp_.RegInt    ( transp_do,      "DoTransp");
    rp_.RegDouble ( transp_theta,   "Theta");
    rp_.RegInt    ( transp_iter,    "Iter");
    rp_.RegDouble ( transp_tol,     "Tol");
    rp_.RegDouble ( transp_diffPos, "DiffPos");
    rp_.RegDouble ( transp_diffNeg, "DiffNeg");
    rp_.RegDouble ( transp_H,       "H");
    rp_.RegDouble ( transp_cPos,    "IniCPos");
    rp_.RegDouble ( transp_cNeg,    "IniCNeg");
    rp_.EndGroup();
}

void ParamDomainCondCL::RegisterParams()
{
    rp_.BeginGroup("DomainCond");
    rp_.RegInt( IniCond,      "InitialCond");
    rp_.RegString( IniData,   "InitialFile");
    rp_.RegString( meshfile,  "MeshFile");
    rp_.RegInt( GeomType,     "GeomType");
    rp_.RegInt( bnd_type,     "BoundaryType");
    rp_.EndGroup();
}

void ParamMesszelleNsCL::RegisterParams()
{

}

void ParamFilmCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt( num_steps,    "NumSteps");
    rp_.RegDouble( dt,        "StepSize");
    rp_.RegDouble( theta,     "ThetaStokes");
    rp_.RegDouble( lset_theta,"ThetaLevelset");
    rp_.EndGroup();

    rp_.BeginGroup("Stokes");
    rp_.RegInt( inner_iter,   "InnerIter");
    rp_.RegInt( outer_iter,   "OuterIter");
    rp_.RegDouble( inner_tol, "InnerTol");
    rp_.RegDouble( outer_tol, "OuterTol");
    rp_.RegInt( StokesMethod, "StokesMethod");
    rp_.RegInt( pcA_iter,     "PcAIter");
    rp_.RegDouble( pcA_tol,   "PcATol");
    rp_.RegDouble( pcS_tol,   "PcSTol");
    rp_.RegDouble(XFEMStab,   "XFEMStab");
    rp_.EndGroup();

    rp_.BeginGroup( "NavStokes");
    rp_.RegDouble( nonlinear,  "Nonlinear");
    rp_.RegDouble( ns_tol,     "Tol");
    rp_.RegDouble( ns_red,     "Reduction");
    rp_.RegInt( ns_iter,       "Iter");
    rp_.EndGroup();

    rp_.BeginGroup("Levelset");
    rp_.RegInt( lset_iter,    "Iter");
    rp_.RegDouble( lset_tol,  "Tol");
    rp_.RegDouble( lset_SD,   "SD");
    rp_.RegDouble( CurvDiff,  "CurvDiff");
    rp_.RegInt( VolCorr,      "VolCorrection");
    rp_.EndGroup();

    rp_.BeginGroup("Coupling");
    rp_.RegInt( cpl_iter,     "Iter");
    rp_.RegDouble( cpl_tol,   "Tol");
    rp_.EndGroup();

    rp_.BeginGroup("Reparam");
    rp_.RegInt( RepFreq,      "Freq");
    rp_.RegInt( RepMethod,    "Method");
    rp_.EndGroup();

    rp_.BeginGroup("AdaptRef");
    rp_.RegInt( ref_freq,     "Freq");
    rp_.RegInt( ref_flevel,   "FinestLevel");
    rp_.RegDouble( ref_width, "Width");
    rp_.EndGroup();

    rp_.BeginGroup("Mat");
    rp_.RegDouble( rhoF,      "DensFluid");
    rp_.RegDouble( muF,       "ViscFluid");
    rp_.RegDouble( rhoG,      "DensGas");
    rp_.RegDouble( muG,       "ViscGas");
    rp_.RegDouble( sm_eps,    "SmoothZone");
    rp_.RegDouble( sigma,     "SurfTension");
    rp_.EndGroup();

    rp_.BeginGroup("Exp");
    rp_.RegDouble( Filmdicke, "Thickness");
    rp_.RegCoord( g,          "Gravity");
    rp_.RegDouble( PumpFreq,  "PumpFreq");
    rp_.RegDouble( PumpAmpl,  "PumpAmpl");
    rp_.RegDouble( AmplZ,     "Ampl_zDir");
    rp_.EndGroup();

    rp_.RegInt( IniCond,      "InitialCond");
    rp_.RegString( EnsCase,   "EnsightCase");
    rp_.RegString( EnsDir,    "EnsightDir");
    rp_.RegString( IniData,   "InitialFile");
    rp_.RegString( BndCond,   "BndCond");
    rp_.RegCoord( mesh_size,  "MeshSize");
    rp_.RegCoord( mesh_res,   "MeshResolution");
    rp_.RegString( serialization_file,    "SerializationFile");
    rp_.RegString( deserialization_file,  "DeserializationFile");
}

} // end of namespace DROPS



