/// \file
/// \brief parameters for  surfactant PDE on the interface.

#include "surfactant/params.h"

namespace DROPS
{

void ParamSurfactantCL::RegisterParams()
{
    rp_.BeginGroup("Time");
    rp_.RegInt( num_steps,    "NumSteps");
    rp_.RegDouble( dt,        "StepSize");
    rp_.RegDouble( theta_surf,"ThetaSurf");
    rp_.RegDouble( lset_theta,"ThetaLevelset");
    rp_.EndGroup();

    rp_.BeginGroup("Levelset");
    rp_.RegInt( lset_iter,    "Iter");
    rp_.RegDouble( lset_tol,  "Tol");
    rp_.RegDouble( lset_SD,   "SD");
    rp_.RegDouble( CurvDiff,  "CurvDiff");
    rp_.RegInt( VolCorr,      "VolCorrection");
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
    rp_.RegDouble( muI,       "ViscIface");
    rp_.EndGroup();

    rp_.BeginGroup("Exp");
    rp_.RegCoord( Radius,     "RadDrop");
    rp_.RegCoord( Mitte,      "PosDrop");
    rp_.RegCoord( Velocity,   "Velocity");
    rp_.EndGroup();

    rp_.RegInt( TestCase,     "TestCase");
    rp_.RegInt( cdiv,         "InitialDivisions");
    rp_.RegString( EnsCase,   "EnsightCase");
    rp_.RegString( EnsDir,    "EnsightDir");
}

} // end of namespace DROPS



