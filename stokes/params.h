/// \file params.h
/// \brief parameters for two-phase flow and other problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Eva Loch, Volker Reichelt, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_STOKES_PARAMS_H
#define DROPS_STOKES_PARAMS_H

#include "misc/params.h"
#include "levelset/params.h"
#include "poisson/params.h"

namespace DROPS
{

/// \brief Parameter class for the material data
class ParamStokesMaterialDataCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name Material data
  ///

  //@{
    double mat_Dens,                        ///< density
           mat_Visc;                        ///< dynamic viscosity

  //@}
  public:
    ParamStokesMaterialDataCL()                             { RegisterParams(); }
    ParamStokesMaterialDataCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the experimental data
class ParamStokesExperimentalDataCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  ///\name Experimental setup
  //@{
    double    exp_InflowVel,                       ///< max. inflow velocity (parabolic profile)
              exp_RadInlet;                        ///< radius at inlet of measuring device
    int       exp_FlowDir;                         ///< flow direction (x/y/z = 0/1/2)
    Point3DCL exp_Gravity;                         ///< gravity
    double    exp_InflowFreq,                      ///< inflow frequence
              exp_InflowAmpl;                      ///< inflow amplitude
  //@}
  public:
    ParamStokesExperimentalDataCL()                             { RegisterParams(); }
    ParamStokesExperimentalDataCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for miscellaneous
class ParamMiscCL : public virtual ParamBaseCL
{
  protected:
	void RegisterParams();

  public:
  /// \name miscellaneous
  //@{
    double      misc_Omega;                            ///< damping factor for SSOR
    double      misc_Tau;                              ///< factor for mass matrix in Uzawa method
    int         misc_ModifyGrid;
    double      misc_MarkLower;						   ///< range from 0 to the maximum value of Y coordinate, mark all tetrahedron (y<MarkLower) for further refinement
  //@}

  public:
      ParamMiscCL()                             { RegisterParams(); }
      ParamMiscCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }

};

/// \brief Parameter class for Stokes coefficents
class ParamStokesCoeffCL : public virtual ParamBaseCL
{
  protected:
	void RegisterParams();

  public:
  /// \name StokesCoeff
  //@{
    std::string    		stc_Reaction;
    std::string    		stc_Source;
    std::string    		stc_Solution_Vel;
    std::string    		stc_Solution_DVel;
    std::string    		stc_Solution_Pr;
    double              ns_Nonlinear;

  //@}

  public:
      ParamStokesCoeffCL()                             { RegisterParams(); }
      ParamStokesCoeffCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }

};

/// \brief Parameter class for the poisson case
class ParamStokesProblemCL:
	public ParamStokesCL,
	public ParamEnsightCL,
	public ParamTimeDiscCL,
	public ParamDomainCondCL,
	public ParamAdaptRefCL,
	public ParamMiscCL,
	public ParamStokesCoeffCL,
	public ParamStokesExperimentalDataCL,
	public ParamStokesMaterialDataCL,
	public ParamErrCL
{
  private:
    void RegisterParams();
  public:
    ParamStokesProblemCL() { RegisterParams(); }
    ParamStokesProblemCL( const std::string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamStokesCL::rp_.ReadParams( file);
        ParamEnsightCL::rp_.ReadParams( file);
        ParamTimeDiscCL::rp_.ReadParams( file);
        ParamDomainCondCL::rp_.ReadParams( file);
        ParamErrCL::rp_.ReadParams( file);
    }
};


} // end of namespace DROPS

#endif
