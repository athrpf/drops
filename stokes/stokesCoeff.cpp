/// \file stokesCoeff.cpp
/// \brief  solution and rhs functions for stokes-problems
/// \author LNM RWTH Aachen: Eva Loch, Yuanjun Zhang

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

#include "misc/bndmap.h"
#include "stokes/params.h"

using namespace std;

extern DROPS::ParamStokesProblemCL C_Stokes;

double ScaZero(const DROPS::Point3DCL&, double =0)
{
  return 0.0;
}

DROPS::SMatrixCL<3, 3> MatZero(const DROPS::Point3DCL& , double =0)
{
  return DROPS::SMatrixCL<3, 3> (0.0);
}

DROPS::SVectorCL<3> VecZero(const DROPS::Point3DCL& , double =0)
{
  return DROPS::SVectorCL<3> (0.0);
}


static DROPS::RegisterScalarFunction regscazero("ScaZero",        ScaZero);
static DROPS::RegisterVectorFunction regveczero("VecZero",        VecZero);
static DROPS::RegisterMatrixFunction regmatdvelsolution("MatZero",MatZero );

// Default parameters
namespace stokes {

DROPS::Point3DCL Source( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.0);
    ret[2]= 3.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
    return ret;
}

DROPS::SVectorCL<3> VelSolution( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    return ret;
}

DROPS::SMatrixCL<3, 3> DVelSolution(const DROPS::Point3DCL& p, double)
{
    DROPS::SMatrixCL<3, 3> ret;
    ret(0,0)= std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret(0,1)= std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret(0,2)= std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;

    ret(1,0)=   std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret(1,1)=   std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret(1,2)= - std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;

    ret(2,0)= -2.*std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    ret(2,1)=  2.*std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;
    ret(2,2)= -2.*std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    return ret;
}

double PrSolution( const DROPS::Point3DCL& p, double)
{
    return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]) - 0.125208551608365;
}

static DROPS::RegisterVectorFunction regvecsource("Source_stokes",     Source );
static DROPS::RegisterVectorFunction regvecvelsolution("VelSol_stokes",   VelSolution );
static DROPS::RegisterMatrixFunction regmatdvelsolution("DVelSol_stokes", DVelSolution );
static DROPS::RegisterScalarFunction regscaprsolution("PrSol_stokes",      PrSolution);

}//end of namespace



// File : FSdrops
namespace FSdrops{

DROPS::Point3DCL Source_stat( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.0);
    ret[2]= 3.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
    return ret;
}

DROPS::SVectorCL<3> VelSolution_stat( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    return ret;
}

double PrSolution_stat( const DROPS::Point3DCL& p, double)
{
	return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]);
}

static DROPS::RegisterVectorFunction regvecsource("Source_FSdrops_stat",         Source_stat );
static DROPS::RegisterVectorFunction regvecvelsolution("VelSol_FSdrops_stat",    VelSolution_stat );
//static DROPS::RegisterMatrixFunction regmatdvelsolution("DVelSol_FSdrops_stat",  MatZero );
static DROPS::RegisterScalarFunction regscaprsolution("PrSol_FSdrops_stat",      PrSolution_stat);


DROPS::SVectorCL<3> helper( const DROPS::Point3DCL& p)
{
	DROPS::SVectorCL<3> ret;
    ret[0]=     std::cos(p[0])*std::sin(p[1])*std::sin(p[2]);
    ret[1]=     std::sin(p[0])*std::cos(p[1])*std::sin(p[2]);
    ret[2]= -2.*std::sin(p[0])*std::sin(p[1])*std::cos(p[2]);
    return ret;
}

DROPS::Point3DCL Source_instat( const DROPS::Point3DCL& p, double t)
{
	return (1 + 3*t)*helper(p);
}

DROPS::SVectorCL<3> VelSolution_instat( const DROPS::Point3DCL& p, double t)
{
	return t*helper(p);
}

static DROPS::RegisterVectorFunction reginvecsource("Source_FSdrops_instat",     Source_instat );
static DROPS::RegisterVectorFunction reginvecvelsolution("VelSol_FSdrops_instat",VelSolution_instat );

DROPS::Point3DCL Source_instat2( const DROPS::Point3DCL& p, double t)
{
	DROPS::SVectorCL<3> ret;
	ret[0]=       2/3*t*std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
	ret[1]=      -2/3*t*std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
	ret[2]= (4/3+3*t)*t*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
	return ret;
}

DROPS::SVectorCL<3> VelSolution_instat2( const DROPS::Point3DCL& p, double t)
{
	DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])*t*t;
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])*t*t;
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*t;
    return ret/3.;
}

double PrSolution_instat2( const DROPS::Point3DCL& p, double t)
{
	return std::cos(p[0])*std::sin(p[1])*std::sin(p[2])*t*t;
}

static DROPS::RegisterVectorFunction regin2vecsource("Source_FSdrops_instat2",     Source_instat2 );
static DROPS::RegisterVectorFunction regin2vecvelsolution("VelSol_FSdrops_instat2",VelSolution_instat2 );
static DROPS::RegisterScalarFunction regin2scaprsolution("PrSol_FSdrops_instat2",  PrSolution_instat2);

DROPS::Point3DCL Source_instat3( const DROPS::Point3DCL& , double )
{
	return DROPS::SVectorCL<3>(-2.);
}

DROPS::SVectorCL<3> VelSolution_instat3( const DROPS::Point3DCL& , double t)
{
	return DROPS::SVectorCL<3>(-t);
}

double PrSolution_instat3( const DROPS::Point3DCL& p, double)
{
	return -(p[0]+p[1]+p[2]);
}

static DROPS::RegisterVectorFunction regin3vecsource("Source_FSdrops_instat3",     Source_instat3 );
static DROPS::RegisterVectorFunction regin3vecvelsolution("VelSol_FSdrops_instat3",VelSolution_instat3 );
static DROPS::RegisterScalarFunction regin3scaprsolution("PrSol_FSdrops_instat3",  PrSolution_instat3);



}//end of namespace


// File : drivcav.cpp
namespace drivcav{

DROPS::SVectorCL<3> Inflow1( const DROPS::Point3DCL& , double)
{
	DROPS::SVectorCL<3> ret(0.);
	ret[0]= 1.;
	return ret;
}

static DROPS::RegisterVectorFunction reginflow1("Inflow_drivcav_v1",Inflow1 );

DROPS::SVectorCL<3> Inflow2( const DROPS::Point3DCL& p, double)
{
    const double st = 0.1;
    const DROPS::SVectorCL<3> ret= DROPS::std_basis<3>(1);
    const double d0= std::fabs(p[0]-.5);
    const double d1= std::fabs(p[1]-.5);
    const double m= std::max(d0, d1);
    return (.5-st<m) ? ((.5-m)/st)*ret : ret;
}

static DROPS::RegisterVectorFunction reginflow2("Inflow_drivcav_v2",Inflow2 );

}//end of namespace

// File : isdrops.cpp
namespace isdrops{

DROPS::Point3DCL Source( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret;
    ret[0]= 2./3.*t*std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
    ret[1]= -2./3.*t*std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
    ret[2]= std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*(4./3. + 3.*t);
    return ret;
}

DROPS::SVectorCL<3> VelSolution( const DROPS::Point3DCL& p, double t)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])*t*t;
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])*t*t;
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*t*t;
    return ret/3.;
}

double PrSolution( const DROPS::Point3DCL& p, double t)
{
	 return std::cos(p[0])*std::sin(p[1])*std::sin(p[2])*t*t
	        -(std::sin( 1.) -2.*std::sin( 1.)*std::cos( 1.)
	        + std::sin( 1.)*std::pow( std::cos( 1.), 2))*t*t; // (...)==0.1778213062
}

static DROPS::RegisterVectorFunction regvecsource("Source_isdrops",     Source );
static DROPS::RegisterVectorFunction regvecvelsolution("VelSol_isdrops",   VelSolution );
static DROPS::RegisterScalarFunction regscaprsolution("PrSol_isdrops",      PrSolution);

}// end of namespace

// File : MGsdropsP2.cpp
namespace MGsdropsP2{
DROPS::Point3DCL Source( const DROPS::Point3DCL& p, double)
{
	 DROPS::SVectorCL<3> ret(0.0);
	 ret[2]= 3.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
	 return ret;
}

DROPS::SVectorCL<3> VelSolution( const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    return ret;
}

DROPS::SMatrixCL<3, 3> DVelSolution(const DROPS::Point3DCL& p, double)
{
    DROPS::SMatrixCL<3, 3> ret;
    ret(0,0)= std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret(0,1)= std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret(0,2)= std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;

    ret(1,0)=   std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret(1,1)=   std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret(1,2)= - std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;

    ret(2,0)= -2.*std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    ret(2,1)=  2.*std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;
    ret(2,2)= -2.*std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    return ret;
}

double PrSolution( const DROPS::Point3DCL& p, double)
{
    return std::cos(p[0])*std::sin(p[1])*std::sin(p[2]) - 0.125208551608365;
}

static DROPS::RegisterVectorFunction regvecsource("Source_MGsdropsP2",     Source );
static DROPS::RegisterVectorFunction regvecvelsolution("VelSol_MGsdropsP2",   VelSolution );
static DROPS::RegisterScalarFunction regscaprsolution("PrSol_MGsdropsP2",      PrSolution);
static DROPS::RegisterMatrixFunction regmatdvelsolution("DVelSol_MGsdropsP2", DVelSolution );
}// end of namespace

namespace stsdrops {

DROPS::SVectorCL<3> VelSolution( const DROPS::Point3DCL& p, double)
{
	DROPS::SVectorCL<3> ret;
	ret[0]=    std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
	ret[1]=  - std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
	ret[2]= 2.*std::cos(p[0])*std::sin(p[1])*std::cos(p[2]);
	return ret/3.;
}

DROPS::SMatrixCL<3, 3> DVelSolution(const DROPS::Point3DCL& p, double)
{
    DROPS::SMatrixCL<3, 3> ret;
    ret(0,0)= std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret(0,1)= std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret(0,2)= std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;

    ret(1,0)=   std::sin(p[0])*std::cos(p[1])*std::sin(p[2])/3.;
    ret(1,1)=   std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    ret(1,2)= - std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;

    ret(2,0)= -2.*std::sin(p[0])*std::sin(p[1])*std::cos(p[2])/3.;
    ret(2,1)=  2.*std::cos(p[0])*std::cos(p[1])*std::cos(p[2])/3.;
    ret(2,2)= -2.*std::cos(p[0])*std::sin(p[1])*std::sin(p[2])/3.;
    return ret;
}

double PrSolution( const DROPS::Point3DCL& p, double)
{
	return std::cos(p[0])*std::sin(p[1])*std::sin(p[2])
	               -(std::sin( 1.) -2.*std::sin( 1.)*std::cos( 1.) + std::sin( 1.)*std::pow( std::cos( 1.), 2)); // (...)==0.1778213062
}

double ConstQ(const DROPS::Point3DCL&, double =0)
{
    return C_Stokes.mat_Dens;
}

DROPS::Point3DCL Source( const DROPS::Point3DCL& p, double)
{
    const double g= C_Stokes.mat_Dens;
    DROPS::SVectorCL<3> ret;
    ret[0]= g/3.*std::sin(p[0])*std::sin(p[1])*std::sin(p[2]);
    ret[1]= -g/3.*std::cos(p[0])*std::cos(p[1])*std::sin(p[2]);
    ret[2]= std::cos(p[0])*std::sin(p[1])*std::cos(p[2])*(2./3.*g + 3.);
    return ret;
}

static DROPS::RegisterVectorFunction regvecsource("Source_stsdrops",     Source );
static DROPS::RegisterVectorFunction regvecvelsolution("VelSol_stsdrops",   VelSolution );
static DROPS::RegisterScalarFunction regscaprsolution("PrSol_stsdrops",      PrSolution);
static DROPS::RegisterMatrixFunction regmatdvelsolution("DVelSol_stsdrops", DVelSolution );
static DROPS::RegisterScalarFunction regscaconstq("ConstQ",      ConstQ);


}// end of namespace
