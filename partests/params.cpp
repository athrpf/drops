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
    rp_.RegDouble(  quality       , "Qualtity" );
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
    rp_.RegInt( poisson_printTime, "PrintTime");
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

} // end of namespace DROPS
