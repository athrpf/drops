//**************************************************************************
// File:    mzelle_hdr.h                                                   *
// Content: header for all main programs simulating the measuring cell     *
// Author:  Sven Gross, Joerg Grande, IGPM RWTH Aachen                     *
//**************************************************************************

#ifndef DROPS_MZELLE_HDR_H
#define DROPS_MZELLE_HDR_H

#include "geom/multigrid.h"
#include "levelset/params.h"
#include "levelset/levelset.h"
#include "num/discretize.h"

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


class InterfaceInfoCL
{
  public:
    DROPS::Point3DCL bary, min, max;
    double maxGrad, Vol, h_min, h_max;

    void Update (const DROPS::LevelsetP2CL& ls) {
        ls.GetInfo( maxGrad, Vol, bary, min, max);
        std::pair<double, double> h= h_interface( ls.GetMG().GetTriangEdgeBegin( ls.Phi.RowIdx->TriangLevel), ls.GetMG().GetTriangEdgeEnd( ls.Phi.RowIdx->TriangLevel), ls.Phi);
        h_min= h.first; h_max= h.second;
    }
} IFInfo; 

double eps=5e-4, // Sprungbreite
    lambda=1.5, // Position des Sprungs zwischen Oberkante (lambda=0) und Schwerpunkt (lambda=1)
    sigma_dirt_fac= 0.8; // gesenkte OFspannung durch Verunreinigungen im unteren Teil des Tropfens 
double sigma;

double sm_step(const DROPS::Point3DCL& p) 
{
    double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
        y= p[1] - y_mid;
    if (y > eps) return sigma;
    if (y < -eps) return sigma_dirt_fac*sigma;
    const double z=y/eps*M_PI/2.; 
    return sigma_dirt_fac*sigma + (sigma - sigma_dirt_fac*sigma) * (std::sin(z)+1)/2;
} 

DROPS::Point3DCL grad_sm_step (const DROPS::Point3DCL& p) 
{
    double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
        y= p[1] - y_mid;
    DROPS::Point3DCL ret;
    if (y > eps) return ret;
    if (y < -eps) return ret;
    const double z=y/eps*M_PI/2.; 
    ret[1]= (sigma - sigma_dirt_fac*sigma) * std::cos(z)/eps*M_PI/4;
    return ret;
} 

double lin(const DROPS::Point3DCL& p)
{
    const double y_top= IFInfo.max[1], 
                 y_bot= IFInfo.bary[1],
                 y_slope= sigma*(1 - sigma_dirt_fac)/(y_top - y_bot);
    return sigma + (p[1] - y_top)*y_slope; 
}

DROPS::Point3DCL grad_lin (const DROPS::Point3DCL&) 
{ 
    const double y_top= IFInfo.max[1], 
                 y_bot= IFInfo.bary[1],
                 y_slope= sigma*(1 - sigma_dirt_fac)/(y_top - y_bot);
    DROPS::Point3DCL ret;
    ret[1]= y_slope;
    return ret; 
}

double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } 
DROPS::Point3DCL gsigma (const DROPS::Point3DCL&, double) { return DROPS::Point3DCL(); } 

double sigma_step(const DROPS::Point3DCL& p, double) { return sm_step( p); } 
DROPS::Point3DCL gsigma_step (const DROPS::Point3DCL& p, double) { return grad_sm_step( p); } 

double One( const DROPS::Point3DCL&) { return 1.; }

#endif
