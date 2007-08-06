/// \file
/// \brief parameters for two-phase flow problems.

#include "levelset/params.h"

namespace DROPS
{

void ParamMesszelleCL::RegisterParams()
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
    rp_.RegDouble( cpl_stab,  "Stab");
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
    rp_.RegDouble( rhoD,      "DensDrop");
    rp_.RegDouble( muD,       "ViscDrop");
    rp_.RegDouble( rhoF,      "DensFluid");
    rp_.RegDouble( muF,       "ViscFluid");
    rp_.RegDouble( sm_eps,    "SmoothZone");
    rp_.RegDouble( sigma,     "SurfTension");
    rp_.EndGroup();

    rp_.BeginGroup("Exp");
    rp_.RegCoord( Radius,     "RadDrop");
    rp_.RegCoord( Mitte,      "PosDrop");
    rp_.RegCoord( g,          "Gravity");
    rp_.RegDouble( Anstroem,  "InflowVel");
    rp_.RegDouble( r_inlet,   "RadInlet");
    rp_.RegInt( flow_dir,     "FlowDir");
    rp_.EndGroup();

    rp_.RegDouble(XFEMStab,   "XFEMStab");
    rp_.RegInt( IniCond,      "InitialCond");
    rp_.RegInt( num_ref,      "NumRef");
    rp_.RegString( EnsCase,   "EnsightCase");
    rp_.RegString( EnsDir,    "EnsightDir");
    rp_.RegString( IniData,   "InitialFile");
    rp_.RegString( meshfile,  "MeshFile");
}

void ParamMesszelleNsCL::RegisterParams()
{
    rp_.BeginGroup( "NavStokes");
    rp_.RegInt( scheme,        "Scheme");
    rp_.RegDouble( nonlinear,  "Nonlinear");
    rp_.RegDouble( stat_nonlinear, "NonlinearStat");
    rp_.RegDouble( stat_theta, "ThetaStat");
    rp_.RegDouble( ns_tol,     "Tol");
    rp_.RegDouble( ns_red,     "Reduction");
    rp_.RegInt( ns_iter,       "Iter");
    rp_.EndGroup();
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
    rp_.EndGroup();

    rp_.RegInt( IniCond,      "InitialCond");
    rp_.RegString( EnsCase,   "EnsightCase");
    rp_.RegString( EnsDir,    "EnsightDir");
    rp_.RegString( IniData,   "InitialFile");
    rp_.RegString( BndCond,   "BndCond");
    rp_.RegCoord( mesh_size,  "MeshSize");
    rp_.RegCoord( mesh_res,   "MeshResolution");
}

} // end of namespace DROPS



