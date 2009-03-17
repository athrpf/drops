//**************************************************************************
// File:    two_phase_hdr.h                                                *
// Content: header for all parallel main programs simulating a             *
//          two-phase flow                                                 *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   03. September 2008                                             *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file two_phase_hdr.h
/// \brief Header for all parallel main programs simulating a two-phase flow

#ifndef DROPS_TWO_PHASE_HDR_H
#define DROPS_TWO_PHASE_HDR_H

#include "geom/multigrid.h"
#include "misc/utils.h"
#include "num/discretize.h"
#include "stokes/stokes.h"
#include "partests/params.h"
#include "levelset/levelset.h"
#include "levelset/mzelle_hdr.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

namespace DROPS{
/// \brief Class for describing surface tension forces
/** This class needs access to an InterfaceInfoCL named as IFInfo.*/
class SurfaceTensionCL
{
  public:
    static double eps;              ///< depth of the jump
    static double lambda;           ///< position of the jump (top: lambda=0, barycenter lambda=1)
    static double sigma;            ///< surface tension at the top of the drop
    static double sigma_dirt_fac;   ///< factor of surface tension at the bottom of the drop

    /// \brief Constant surface tension
    static double sigmaf (const Point3DCL&, double) {
        return sigma;
    }

    /// \brief Gradient of constant surface tension
    static Point3DCL gsigma (const Point3DCL&, double) {
        return Point3DCL();
    }

    /// \brief Surface tension modeled by a step
    static double sigma_step(const Point3DCL& p, double)
    {
        double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
               y    = p[1] - y_mid;
        if (y > eps) return sigma;
        if (y < -eps) return sigma_dirt_fac*sigma;
        const double z=y/eps*M_PI/2.;
        return sigma_dirt_fac*sigma + (sigma - sigma_dirt_fac*sigma) * (std::sin(z)+1)/2;
    }

    /// \brief Gradient of surface tension modeled by a step
    static Point3DCL gsigma_step (const Point3DCL& p, double)
    {
        double y_mid= lambda*IFInfo.bary[1] + (1-lambda)*IFInfo.max[1], // zwischen Tropfenschwerpunkt und Oberkante
               y    = p[1] - y_mid;
        Point3DCL ret;
        if (y > eps) return ret;
        if (y < -eps) return ret;
        const double z=y/eps*M_PI/2.;
        ret[1]= (sigma - sigma_dirt_fac*sigma) * std::cos(z)/eps*M_PI/4;
        return ret;
    }
};

/// \brief Create all numberings for stokes DOF and level-set DOF and assign them
template<typename StokesT, typename LevelsetT>
void CreateIdxAndAssignIdx(StokesT& Stokes, LevelsetT& lset, const MultiGridCL& mg)
{
    // Create numbering
    Stokes.CreateNumberingVel( mg.GetLastLevel(), &Stokes.vel_idx);
    Stokes.CreateNumberingPr ( mg.GetLastLevel(), &Stokes.pr_idx);
    lset.CreateNumbering     ( mg.GetLastLevel(), &lset.idx);

    // Tell matrices and vectors about the numbering
    Stokes.v.SetIdx(   &Stokes.vel_idx);
    Stokes.p.SetIdx(   &Stokes.pr_idx);
    Stokes.b.SetIdx(   &Stokes.vel_idx);
    Stokes.c.SetIdx(   &Stokes.pr_idx);
    Stokes.A.SetIdx(   &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.B.SetIdx(   &Stokes.pr_idx,  &Stokes.vel_idx);
    Stokes.M.SetIdx(   &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.N.SetIdx(   &Stokes.vel_idx, &Stokes.vel_idx);
    Stokes.prM.SetIdx( &Stokes.pr_idx,  &Stokes.pr_idx);
    Stokes.prA.SetIdx( &Stokes.pr_idx,  &Stokes.pr_idx);
    lset.Phi.SetIdx(   &lset.idx);
}



/*******************************************************************************
*    G L O B A L   A N D   S T A T I C   V A R I A B L E S                     *
*******************************************************************************/


const char line[] ="------------------------------------------------------------";

// Init of static members
double SurfaceTensionCL::eps           = 5e-4;
double SurfaceTensionCL::lambda        = 1.5;
double SurfaceTensionCL::sigma         = 0.0;
double SurfaceTensionCL::sigma_dirt_fac= 0.8;

}       // end of namespace

#endif
