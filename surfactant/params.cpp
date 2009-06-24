/// \file
/// \brief parameters for  surfactant PDE on the interface.

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



