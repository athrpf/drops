/***************************************************************************
*  File:    params.cpp                                                     *
*  Content: Reading of Parameterfiles for parallel programms               *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           10.11.2005                                             *
*  last modified:   10.11.2005                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file partests/params.cpp

#include "partests/params.h"

namespace DROPS
{

void ParamLoadBalCL::RegisterParams()
{
    rp_.BeginGroup( "LoadBalancing" );
    rp_.RegInt(     refineStrategy, "RefineStrategy" );
    rp_.RegDouble(  quality       , "Qualtity" );
    rp_.EndGroup();
}

void ParamEnsightCL::RegisterParams()
{
    rp_.BeginGroup( "Ensight" );
    rp_.RegInt(     ensight,   "EnsightOut" );
    rp_.RegString(  ensCase,   "EnsCase" );
    rp_.RegString(  ensDir,    "EnsDir" );
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
    rp_.RegInt(     vtkBinary,    "Binary");
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
    rp_.BeginGroup( "MGSerialization");
    rp_.RegInt(serialization, "Serialization");
    rp_.RegInt(overwrite, "Overwrite");
    rp_.RegString(ser_dir, "SerDir");
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
    rp_.RegInt(   RepFreq,   "Freq");
    rp_.RegInt(   RepMethod, "Method");
    rp_.RegDouble(NarrowBand,"NarrowBand");
    rp_.EndGroup();
}

void ParamApdaptRefCL::RegisterParams()
{
    rp_.BeginGroup("AdaptRef");
    rp_.RegInt( ref_freq,     "Freq");
    rp_.RegInt( ref_flevel,   "FinestLevel");
    rp_.RegDouble( ref_width, "Width");
    rp_.EndGroup();
}

void ParamParRefCL::RegisterParams()
{
    rp_.BeginGroup( "Refining");
    rp_.RegInt(  init_cond , "InitCond");
    rp_.RegInt(  refined,    "Refined");
    rp_.RegInt(  markall, "MarkAll");
    rp_.RegInt(  markdrop, "MarkDrop");
    rp_.RegInt(  markcorner, "MarkCorner");
    rp_.RegInt(  markingproc, "MarkingProc");
    rp_.RegInt(  Strategy, "Strategy");
    rp_.EndGroup();

    rp_.BeginGroup( "Coarsening");
    rp_.RegInt( coarsedrop, "CoarseDrop");
    rp_.RegInt( coarseall, "CoarseAll");
    rp_.RegInt( unmarkingproc, "UnMarkingProc");
    rp_.EndGroup();

    rp_.BeginGroup( "LoadBalancing");
    rp_.RegInt( refineStrategy, "RefineStrategy");
    rp_.RegInt( coarseStrategy, "CoarseStrategy");
    rp_.RegInt( middleMig, "MiddleMig");
    rp_.EndGroup();

    rp_.BeginGroup( "Misc");
    rp_.RegInt( printSize, "PrintSize");
    rp_.RegInt( printPMG, "PrintPMG");
    rp_.RegInt( printGEO, "PrintGEO");
    rp_.RegInt( printTime, "PrintTime");
    rp_.RegInt( checkRef, "CheckAfterRef");
    rp_.RegInt( checkMig, "CheckAfterMig");
    rp_.RegInt( checkDDD, "CheckDDD");
    // rp_.RegInt( exitOnMGErr, "ExitOnMGErr");
    rp_.RegString(init_pre, "InitPrefix");
    rp_.EndGroup();
}


void ParamParStokesCL::RegisterParams()
{
    rp_.BeginGroup( "Stokes");
    rp_.RegDouble( nu, "nu");
    rp_.EndGroup();

    rp_.BeginGroup( "Refining");
    rp_.RegInt( basicref_x, "BasicRefX");
    rp_.RegInt( basicref_y, "BasicRefY");
    rp_.RegInt( basicref_z, "BasicRefZ");
    rp_.RegDouble( dx, "dx");
    rp_.RegDouble( dy, "dy");
    rp_.RegDouble( dz, "dz");
    rp_.RegInt( refall, "RefAll");
    rp_.EndGroup();

    rp_.BeginGroup( "LoadBalancing");
    rp_.RegInt( refineStrategy, "RefineStrategy");
    rp_.EndGroup();

    rp_.BeginGroup( "Solver");
    rp_.RegDouble( relax, "Relax");
    rp_.RegInt(    pc_iter,   "PCIter");
    rp_.RegDouble( pc_rel_tol,"PCRelTol");
    rp_.RegInt( inner_iter, "InnerIter");
    rp_.RegInt( outer_iter, "OuterIter");
    rp_.RegDouble( inner_tol, "InnerTol");
    rp_.RegDouble( outer_tol, "OuterTol");
    rp_.RegInt( restart, "Restart");
    rp_.RegInt( relative, "Relative");
    rp_.RegDouble( reduction, "Reduction");
    rp_.RegInt( accur, "Accur");
    rp_.EndGroup();

    rp_.BeginGroup( "Misc");
    rp_.RegInt( printInfo, "PrintInfo");
    rp_.EndGroup();
}

void ParamParInstatStokesCL::RegisterParams()
{
    rp_.BeginGroup( "ExpData");
    rp_.RegDouble(inflowVel, "InflowVel");
    rp_.RegCoord( g,         "g");
    rp_.RegDouble(frequence, "Frequence");
    rp_.RegDouble(ampl,      "Amplitude");
    rp_.EndGroup();

    rp_.BeginGroup( "Time");
    rp_.RegInt(   timesteps, "TimeSteps");
    rp_.RegDouble(theta,     "Theta");
    rp_.RegDouble(stepsize,  "StepSize");
    rp_.EndGroup();
}

void ParamParPoissonCL::RegisterParams()
{
    rp_.BeginGroup( "Poisson");
    rp_.RegDouble( nu, "nu");
    rp_.EndGroup();

    rp_.BeginGroup( "Refining");
    rp_.RegInt( refall, "RefAll");
    rp_.RegInt( markdrop, "MarkDrop");
    rp_.RegInt( markcorner, "MarkCorner");
    rp_.RegInt( adaptiv,  "Adpativ");
    rp_.EndGroup();

    rp_.BeginGroup( "LoadBalancing");
    rp_.RegInt( refineStrategy, "RefineStrategy");
    rp_.RegInt( transferUnks,   "TransferUnknowns");
    rp_.EndGroup();

    rp_.BeginGroup( "Solver");
    rp_.RegInt(    solver,      "Solver");
    rp_.RegInt(    precond,     "Precond");
    rp_.RegDouble( relax,       "Relax");
    rp_.RegInt(    pciter,      "PCIter");
    rp_.RegDouble( pctol,       "PCTol");
    rp_.RegInt(    iter,        "Iteration");
    rp_.RegDouble( tol,         "Tol");
    rp_.RegInt(    restart,     "Restart");
    rp_.RegInt(    useMGS,      "UseMGS");
    rp_.RegInt(    relative,    "Relative");
    rp_.RegInt(    accur,       "Accur");
    rp_.RegInt(    modified,    "Modified");
    rp_.RegInt(    preCondMeth, "PreCondMeth");

    rp_.EndGroup();

    rp_.BeginGroup( "Misc");
    rp_.RegInt( printGEO,      "PrintGEO");
    rp_.RegInt( printMG,       "PrintMG");
    rp_.RegInt( printSize,     "PrintSize");
    rp_.RegInt( check,         "CheckMG");
    rp_.RegInt( printUnknowns, "PrintUnknowns");
    rp_.EndGroup();

    rp_.BeginGroup( "Ensight");
    rp_.RegInt(   ensight, "ensight");
    rp_.RegString(EnsCase, "EnsCase");
    rp_.RegString(EnsDir, "EnsDir");
    rp_.RegString(geomName, "GeomName");
    rp_.RegString(varName, "VarName");
    rp_.EndGroup();
}

void ParamParExchangeCL::RegisterParams()
{
    rp_.BeginGroup( "Unknowns");
    rp_.RegInt( numsV1, "NumsV1");
    rp_.RegInt( numsE1, "NumsE1");
    rp_.RegInt( numsT1, "NumsT1");
    rp_.RegInt( numsV2, "NumsV2");
    rp_.RegInt( numsE2, "NumsE2");
    rp_.RegInt( numsT2, "NumsT2");
    rp_.EndGroup();

    rp_.BeginGroup( "Refining");
    rp_.RegInt( basicref_x, "BasicRefX");
    rp_.RegInt( basicref_y, "BasicRefY");
    rp_.RegInt( basicref_z, "BasicRefZ");
    rp_.RegDouble( dx, "dx");
    rp_.RegDouble( dy, "dy");
    rp_.RegDouble( dz, "dz");
    rp_.RegInt( refall, "RefAll");
    rp_.RegInt( refineStrategy,  "RefineStrategy");
    rp_.RegInt( migEveryTime, "MigEveryTime");
    rp_.EndGroup();

    rp_.BeginGroup( "Misc");
    rp_.RegInt( printMG,  "PrintMG");
    rp_.RegInt( printEx,  "PrintEx");
    rp_.RegInt( checkMG,  "CheckMG");
    rp_.RegInt( timeMeas, "TimeMeas");
    rp_.RegInt( tests, "Tests");
    rp_.RegInt( multruns, "MultRuns");
    rp_.RegInt( printMsgSize, "PrintMsgSize");
    rp_.EndGroup();
}

void ParamParNavStokesCL::RegisterParams()
{
    rp_.BeginGroup("NavierStokes");
    rp_.RegDouble( reduction, "Reduction");
    rp_.RegInt(    nav_iter,  "Iter");
    rp_.RegDouble( nav_tol,   "Tol");
    rp_.RegInt(    markTop,   "RefTop");
    rp_.EndGroup();
}

void ParamParBrickFlowCL::RegisterParams()
{
}

void ParParamMesszelleNsCL::RegisterParams()
{
}

void ParamParFilmCL::RegisterParams()
{
    rp_.BeginGroup("Geometry");
    rp_.RegCoord(  Dim,    "Dim");
    rp_.RegInt(    nx,     "nx");
    rp_.RegInt(    ny,     "ny");
    rp_.RegInt(    nz,     "nz");
    rp_.RegInt(    Refine, "Refine");
    rp_.EndGroup();

    rp_.BeginGroup("Material");
    rp_.RegDouble( KineticViscosity,  "KineticViscosity");
    rp_.RegDouble( HeatConductivity,  "HeatConductivity");
    rp_.RegDouble( InflowTemperature, "InflowTemperature");
    rp_.RegDouble( WallTemperature,   "WallTemperature");
    rp_.RegCoord ( g,                 "g");
    rp_.EndGroup();

    rp_.BeginGroup("Solver");
    rp_.RegDouble( Relax,   "Relax");
    rp_.RegInt(    Restart, "Restart");
    rp_.RegInt(    Iter,    "Iter");
    rp_.RegDouble( Tol,     "Tol");
    rp_.RegInt(    UseMGS,  "UseMGS");
    rp_.EndGroup();

    rp_.BeginGroup( "Ensight");
    rp_.RegInt(   Ensight,  "Ensight");
    rp_.RegString(EnsCase,  "EnsCase");
    rp_.RegString(EnsDir,   "EnsDir");
    rp_.RegString(geomName, "GeomName");
    rp_.RegString(varName,  "VarName");
    rp_.EndGroup();
}

void ParamParSerCL::RegisterParams()
{
    rp_.BeginGroup( "Refining");
    rp_.RegInt( markall,     "MarkAll");
    rp_.RegInt( markdrop,    "MarkDrop");
    rp_.RegInt( markcorner,  "MarkCorner");
    rp_.RegInt( markingproc, "MarkingProc");
    rp_.EndGroup();

    rp_.BeginGroup( "Misc");
    rp_.RegInt( mode,     "Mode");
    rp_.RegInt( unknowns, "Unknowns");
    rp_.EndGroup();
}

void ParamParFilmStokesCL::RegisterParams()
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
    rp_.EndGroup();

    rp_.RegInt( IniCond,      "InitialCond");
    rp_.RegString( IniData,   "InitialFile");
}

} // end of namespace DROPS
