/// \file transportCoeff.cpp
/// \brief boundary and source functions for the transport Problem
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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
#include "levelset/params.h"

double Initialcneg (const DROPS::Point3DCL& , double )
{
    extern DROPS::ParamMesszelleNsCL C;
    return C.trp_IniCNeg;
}

double Initialcpos (const DROPS::Point3DCL& , double )
{
    extern DROPS::ParamMesszelleNsCL C;
    return C.trp_IniCPos;
}

double Reaction (const DROPS::Point3DCL& , double )
{
  return 0.0;
}

double ZeroFct (const DROPS::Point3DCL& ,  double )
{  
    return 0.;
}

double tid_ZeroFct (const DROPS::Point3DCL&)
{  
    return 0.;
}

double Rhs (const DROPS::Point3DCL& ,  double )
{  
    return 0.;
}

double Dirichlet (const DROPS::Point3DCL& p, double )
{
/*  static double x0 = C.exp_PosDrop[0];
  static double y0 = C.exp_PosDrop[1];
  static double R = C.exp_RadDrop[0];  */
//  double x = p[0];
  double y = p[1];  
  
  int c= (int) (y * 2.0);
  switch(c)
  {
    case 0: return 0; break;
    case 1: return (y-0.5)*(1-y); break;
    case 2: return -(y-1)*(1.5-y); break;
/*    case 1:
    case 2:
            return (y-0.5)*(1.5-y);
*/    case 3: 
    default:
        return 0;
        break; 
  }
}

double DirichletConst (const DROPS::Point3DCL& p, double )
{
  extern DROPS::ParamMesszelleNsCL C;
  double y=p[1]-C.exp_PosDrop[1];
  if (y>0)
    return C.trp_IniCPos;    
  else
    return C.trp_IniCNeg;    
}

double DirichletPos (const DROPS::Point3DCL&, double )
{
  extern DROPS::ParamMesszelleNsCL C;
  return C.trp_IniCPos;  
}

double DirichletConstt (const DROPS::Point3DCL& p, double )
{
  extern DROPS::ParamMesszelleNsCL C;
  double y=p[1]-C.exp_PosDrop[1];
  if (y>0)
    return C.trp_IniCPos;    
  else
    return C.trp_IniCNeg * C.trp_HNeg / C.trp_HPos;    
}

DROPS::SVectorCL<3> PotentialFlowfield (const DROPS::Point3DCL& p, double )
{  
    extern DROPS::ParamMesszelleNsCL C;  
    DROPS::SVectorCL<3> ret(0.);
    double x=p[0]-C.exp_PosDrop[0]; double y=p[1]-C.exp_PosDrop[1];
    double r2 = x*x+y*y;
    double r4 = r2*r2;
    double R=C.exp_RadDrop[0];
    static bool test=true;
    if(test) { std::cout << "R = " << R << std::endl; test=false;}
    double R2 = R*R;
    if (R2>=r2){
      ;//nothing, v=0
    }else{
      ret[0] = 1 + R2*(y*y-x*x)/r4;
      ret[1] = -2*R2*(y*x)/r4;
    }
//    std::cout << ret[0] << " " << ret[1] << std::endl;
    return ret;
}

template<int i>
DROPS::SVectorCL<3> StraightFlowfield (const DROPS::Point3DCL&, double )
{  
    DROPS::SVectorCL<3> ret(0.);
    ret[i] = 1.0;
    return ret;
}


double cylinderdistance( const DROPS::Point3DCL& p)
{
    extern DROPS::ParamMesszelleNsCL C;
    double x=p[0]-C.exp_PosDrop[0]; double y=p[1]-C.exp_PosDrop[1];
    double R=C.exp_RadDrop[0];
    static bool test=true;
    if(test) { std::cout << "R = " << R << std::endl; test=false;}
    return std::sqrt(std::abs(x*x+y*y)) - R;
}

template<int i>
double planedistance( const DROPS::Point3DCL& p)
{
    extern DROPS::ParamMesszelleNsCL C;
    double x=p[i]-C.exp_PosDrop[i]; 
    return x;
}


static DROPS::RegisterScalarFunction regsca_inicneg("IniCnegFct", Initialcneg);
static DROPS::RegisterScalarFunction regsca_inicpos("IniCposFct", Initialcpos);
static DROPS::RegisterScalarFunction regsca_reaction("ReactionFct", Reaction);
static DROPS::RegisterScalarFunction regsca_zero("ZeroFct", ZeroFct);
static DROPS::RegisterScalarFunction regsca_tidzero("ZeroFct", tid_ZeroFct);
static DROPS::RegisterScalarFunction regsca_Rhs("Rhs", Rhs);
static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletPos);
static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletPos);
//static DROPS::RegisterScalarFunction regsca_dir("Dirichlet", DirichletConst);
//static DROPS::RegisterScalarFunction regsca_dirt("Dirichlett", DirichletConstt);
static DROPS::RegisterScalarFunction regsca_cdist("cylinderdistance", cylinderdistance);
static DROPS::RegisterScalarFunction regsca_pdistx("planedistancex", planedistance<0>);
static DROPS::RegisterScalarFunction regsca_pdisty("planedistancey", planedistance<1>);
static DROPS::RegisterScalarFunction regsca_pdistz("planedistancez", planedistance<2>);
static DROPS::RegisterVectorFunction regvec_potflow("potflow", PotentialFlowfield);
static DROPS::RegisterVectorFunction regvec_strflowx("straightflowx", StraightFlowfield<0>);
static DROPS::RegisterVectorFunction regvec_strflowy("straightflowy", StraightFlowfield<1>);
static DROPS::RegisterVectorFunction regvec_strflowz("straightflowz", StraightFlowfield<2>);


