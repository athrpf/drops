//**************************************************************************
// File:    mzelle_hdr.h                                                   *
// Content: header for all main programs simulating the measuring cell     *
// Author:  Sven Gross, Joerg Grande, IGPM RWTH Aachen                     *
//**************************************************************************

#ifndef DROPS_MZELLE_HDR_H
#define DROPS_MZELLE_HDR_H

#include "geom/multigrid.h"
#include "levelset/params.h"

// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0


class ZeroFlowCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    ZeroFlowCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( C.rhoD, C.rhoF ), DROPS::H_sm, C.sm_eps),
        mu(  DROPS::JumpCL( C.muD,  C.muF),   DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma), g( C.g)    {}
};

class DimLessCoeffCL
{
// \Omega_1 = Tropfen,    \Omega_2 = umgebendes Fluid
  public:
    static DROPS::Point3DCL f(const DROPS::Point3DCL&, double)
        { DROPS::Point3DCL ret(0.0); return ret; }
    const DROPS::SmoothedJumpCL rho, mu;
    const double SurfTens;
    const DROPS::Point3DCL g;

    DimLessCoeffCL( const DROPS::ParamMesszelleCL& C)
      : rho( DROPS::JumpCL( 1., C.rhoF/C.rhoD ), DROPS::H_sm, C.sm_eps),
        mu ( DROPS::JumpCL( 1., C.muF/C.muD),    DROPS::H_sm, C.sm_eps),
        SurfTens( C.sigma/C.rhoD), g( C.g)    {}
};

class EllipsoidCL
{
  private:
    static DROPS::Point3DCL Mitte_;
    static DROPS::Point3DCL Radius_;
    
  public:
    EllipsoidCL( const DROPS::Point3DCL& Mitte, const DROPS::Point3DCL& Radius)
    { Init( Mitte, Radius); }
    static void Init( const DROPS::Point3DCL& Mitte, const DROPS::Point3DCL& Radius)
    { Mitte_= Mitte;    Radius_= Radius; }
    static double DistanceFct( const DROPS::Point3DCL& p)
    {
        DROPS::Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static double GetVolume() { return 4./3.*M_PI*Radius_[0]*Radius_[1]*Radius_[2]; }
};

DROPS::Point3DCL EllipsoidCL::Mitte_;
DROPS::Point3DCL EllipsoidCL::Radius_;


double sigma;
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } 

double One( const DROPS::Point3DCL&) { return 1.; }

#endif
