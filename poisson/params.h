/// \file params.h
/// \brief parameters for two-phase flow and other problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_POISSON_PARAMS_H
#define DROPS_POISSON_PARAMS_H

#include "misc/params.h"
#include "levelset/params.h"

namespace DROPS
{

///\brief Parameter class for Poisson equation
class ParamPoissonCL : public virtual ParamBaseCL
{
  protected:
	void RegisterParams();

  public:
  /// \name parameter for the Poisson equation
  //@{
	int    pos_PcIter;
	double pos_PcTol;
	int    pos_Iter;
	int    pos_Restart;
	int    pos_RelativeErr;
	double pos_Tol;
	int    pos_Method;
	int    pos_SmoothingSteps;
	int    pos_NumLvl;
	int    pos_SolutionIsKnown;
	double pos_Relax;

  //@}
  public:
    ParamPoissonCL()                             { RegisterParams(); }
    ParamPoissonCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }

};

/// \brief Parameter class for time discretization
class ParamTimeDiscPoissonCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  /// \name parameter for the time discretization
  //@{
    double tm_StepSize;                            ///< time step size
    int    tm_NumSteps;                            ///< number of timesteps
    int    tm_Scheme;                              ///< not used a.t.m.
    double tm_Theta;
    double tm_Nu;
    int    tm_Convection;

  //@}
  public:
    ParamTimeDiscPoissonCL()                             { RegisterParams(); }
    ParamTimeDiscPoissonCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the experimental data
class ParamExperimentalDataPoissonCL : public virtual ParamBaseCL
{
  protected:
    void RegisterParams();

  public:
  ///\name Experimental setup
  //@{
    double exp_Heat,
           exp_Rho,
           exp_Mu,
           exp_Cp,
           exp_Lambda;
  //@}
  public:
    ParamExperimentalDataPoissonCL()                             { RegisterParams(); }
    ParamExperimentalDataPoissonCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for error estimator
class ParamErrCL : public virtual ParamBaseCL
{
  protected:
	void RegisterParams();

  public:
  /// \name error estimator
  //@{
    int    		err_DoErrorEstimate;
    double 		err_RelReduction;
    double 		err_MinRatio;
    double 		err_Threshold;
    double 		err_Meas;
    int    		err_DoMark;
    int    		err_NumRef;
  //@}

  public:
      ParamErrCL()                             { RegisterParams(); }
      ParamErrCL( const std::string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

/// \brief Parameter class for the poisson case
class ParamPoissonProblemCL:
	    public ParamPoissonCL,
        public ParamEnsightCL,
        public ParamTimeDiscPoissonCL,
        public ParamExperimentalDataPoissonCL,
        public ParamDomainCondCL,
        public ParamErrCL,
        public ParamVTKCL

{
  private:
    void RegisterParams();
  public:
    ParamPoissonProblemCL() { RegisterParams(); }
    ParamPoissonProblemCL( const std::string& filename) {
        RegisterParams();
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamPoissonCL::rp_.ReadParams( file);
        ParamTimeDiscPoissonCL::rp_.ReadParams( file);
        ParamEnsightCL::rp_.ReadParams( file);
        ParamExperimentalDataPoissonCL::rp_.ReadParams( file);
        ParamDomainCondCL::rp_.ReadParams( file);
        ParamErrCL::rp_.ReadParams( file);
        ParamVTKCL::rp_.ReadParams( file);
    }
};
} // end of namespace DROPS

#endif





