/// \file twophasedrops.cpp
/// \brief flow in measurement cell or brick, 
///        basically a copy of twophasedrops.cpp with extensions for mass transport with NitscheXFEM
/// \author LNM RWTH Aachen: Hieu Nguyen, Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Christoph Lehrenfeld; SC RWTH Aachen: Oliver Fortmeier

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

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "navstokes/instatnavstokes2phase.h"
#include "stokes/integrTime.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "levelset/coupling.h"
#include "misc/params.h"
#include "levelset/adaptriang.h"
#include "levelset/mzelle_hdr.h"
#include "levelset/surfacetension.h"
#include "transport/transportNitsche.h"
#include "surfactant/ifacetransp.h"
#include "levelset/twophaseutils.h"
#include "misc/bndmap.h"

#include "num/stokessolverfactory.h"
#ifndef _PAR
#include "num/stokessolver.h"
#else
#include "num/parstokessolver.h"
#include "parallel/loadbal.h"
#include "parallel/parmultigrid.h"
#endif
#include <fstream>
#include <sstream>

/// \todo serialization of mass transp. data
/// \todo boundary values for concentration should get more elegant (like poisson perhaps?)
/// \todo solutiononpart-output for ensight AND vtk
/// \todo surfacetension, varSurfaceTension ... flags, output and cases!
DROPS::ParamCL P;
// rho*du/dt - mu*laplace u + Dp = f + rho*g - okn
//                        -div u = 0
//                             u = u0, t=t0


class InflowLset{
  private:
    static DROPS::scalar_fun_ptr fun;
  public:
    InflowLset(DROPS::scalar_fun_ptr ffun){Init(ffun);}
    static void Init(DROPS::scalar_fun_ptr ffun) {fun = ffun;}
    static double Fct( const DROPS::Point3DCL& p, double){
      return fun(p);
    }
};
DROPS::scalar_fun_ptr InflowLset::fun;



typedef DROPS::BndDataCL<> cBndDataCL;
/// \name Initial data for transport equation
//@{
typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;

/// \name Initial data and rhs for surfactant transport
//@{
const double a( -13./8.*std::sqrt( 35./M_PI));
double surf_rhs (const DROPS::Point3DCL& p, double)
{
    return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
double surf_sol (const DROPS::Point3DCL& p, double)
{
    return 1. + std::sin( atan2( p[0] -P.get<DROPS::Point3DCL>("Exp.PosDrop")[0], p[2] -P.get<DROPS::Point3DCL>("Exp.PosDrop")[2]));
}
//@}

namespace DROPS // for Strategy
{
 
void InitVel (const MultiGridCL& MG, VecDescCL& v, instat_vector_fun_ptr vf)
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
 
    DROPS_FOR_TRIANG_CONST_VERTEX( MG, lvl, it)
        if (it->Unknowns.Exist( idx))
            DoFHelperCL<Point3DCL, VectorCL>::set( v.Data, it->Unknowns( idx),
                vf( it->GetCoord(), 0.));
 
    DROPS_FOR_TRIANG_CONST_EDGE( MG, lvl, it)
        if (it->Unknowns.Exist( idx))
            DoFHelperCL<Point3DCL, VectorCL>::set( v.Data, it->Unknowns( idx),
                vf( BaryCenter( it->GetVertex( 0)->GetCoord(), it->GetVertex( 1)->GetCoord()), 0.));
}


void  OnlyTransportStrategy( MultiGridCL& MG, LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)    // do just the transport stuff
{
    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & tdscalarmap = InScaMap::getInstance();
    ScaMap & scalarmap = ScaMap::getInstance();
    instat_vector_fun_ptr Flowfield = tdvectormap[P.get("Transp.Flow", std::string("ZeroVel"))];
    instat_scalar_fun_ptr Reaction = tdscalarmap["ReactionFct"];
    instat_scalar_fun_ptr Rhs = tdscalarmap["Rhs"];
    instat_scalar_fun_ptr Initialcneg = tdscalarmap["IniCnegFct"];
    instat_scalar_fun_ptr Initialcpos = tdscalarmap["IniCposFct"];
    scalar_fun_ptr distance = scalarmap[P.get("Transp.Levelset", std::string("Ellipsoid"))];

    cBndDataCL *pBnd_c, *pBnd_ct;
    DROPS::BuildBoundaryData( &MG, pBnd_c,  P.get<std::string>("Transp.BoundaryType","21!2!21!21!21!21"), P.get<std::string>("Transp.BoundaryFncs","Zero!Dirichlet!Zero!Zero!Zero!Zero"));
    DROPS::BuildBoundaryData( &MG, pBnd_ct, P.get<std::string>("Transp.BoundaryType","21!2!21!21!21!21"), P.get<std::string>("Transp.BoundaryFncs_t","Zero!Dirichlett!Zero!Zero!Zero!Zero"));
    cBndDataCL & Bnd_c(*pBnd_c);
    cBndDataCL & Bnd_ct(*pBnd_ct); 
   
    DROPS::instat_scalar_fun_ptr sigmap = 0;
    SurfaceTensionCL sf( sigmap, Bnd_c);    
    
    LevelsetP2CL lset( MG, lsetbnddata, sf, P.get<double>("Levelset.SD"), P.get<double>("Levelset.CurvDiff"));
    // levelset wrt the previous time step:
    LevelsetP2CL oldlset( MG, lsetbnddata, sf, P.get<double>("Levelset.SD"), P.get<double>("Levelset.CurvDiff"));
    //Prolongate and Restrict solution vector levelset from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    LevelsetRepairCL oldlsetrepair( oldlset);
    adap.push_back( &oldlsetrepair);
    IdxDescCL* lidx= &lset.idx;
    // index wrt the interface at previous time step
    IdxDescCL* oldlidx= &oldlset.idx; 
    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    oldlset.CreateNumbering( MG.GetLastLevel(), oldlidx);
    oldlset.Phi.SetIdx( oldlidx);
    SetInitialLevelsetConditions( lset, MG, P);
    SetInitialLevelsetConditions( oldlset, MG, P);
    lset.Init( distance );
    oldlset.Init( distance);
    DisplayDetailedGeom( MG);
    const double Vol= lset.GetVolume(); //0.5 * 0.125 * M_PI; //EllipsoidCL::GetVolume();
    std::cout << "initial volume(abs value): " << lset.GetVolume() << std::endl;
    
    //VelocityContainer vel(Stokes.v,Stokes.GetBndData().Vel,MG);
    VelocityContainer vel(Flowfield);
    
    TransportP1XCL massTransp( MG, Bnd_c, Bnd_ct, vel, lsetbnddata, lset.Phi, oldlset.Phi,P,0,Reaction,Rhs);
    TransportXRepairCL transprepair(massTransp, MG.GetLastLevel());
    
    // index of the concentration wrt the interface at actual time step:
    MLIdxDescCL* cidx= &massTransp.idx;
    
    // index of the concentration wrt the interface at previous time step:
    MLIdxDescCL* cidx_old= &massTransp.oldidx; 
    
    //This following Vector c_out is responsable for the communication from concentration to surface tension. 
    //Before a new surface tension is computed c_out should be updated (via GetSolutionOnPart)
    //Important: This vector has to be kept in memory as long as the surface tension is computed!
    VecDescCL c_out;
    //c_in: s. c_out but c_in is only used for visualization until now
    VecDescCL c_in;
    IdxDescCL p1idx(P1_FE,Bnd_c,0);
    /// \todo for periodic stuff: matching function here
    p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, MG.GetBnd().GetMatchFun(), &lset.Phi,&lsetbnddata);   
    c_in.SetIdx( &p1idx);
    c_out.SetIdx( &p1idx);
    
    {
        adap.push_back(&transprepair);
        massTransp.CreateNumbering( MG.GetLastLevel(), cidx, cidx_old, lset.Phi, oldlset.Phi);
        massTransp.ct.SetIdx( cidx);
        massTransp.c.SetIdx( cidx);
        massTransp.oldct.SetIdx( cidx_old);
        std::cout << massTransp.ct.Data.size() << " concentration unknowns,\n";
 
        if (P.get<int>("DomainCond.InitialCond") != -1)
          massTransp.Init( Initialcneg, Initialcpos, 0);
        else
          ReadFEFromFile( massTransp.ct, MG, P.get<std::string>("DomainCond.InitialFile")+"concentrationTransf");
          
        massTransp.GetSolutionOnPart( c_out, true , false);
        massTransp.GetSolutionOnPart( c_in, false , false);
        double c_mean = massTransp.MeanDropConcentration();
        std::cout << "START:: Mean concentration in drop: " << std::setprecision(12) << c_mean <<"\n";        
    }
    
    LevelsetModifyCL lsetmod( P.get<int>("Reparam.Freq"), P.get<int>("Reparam.Method"), P.get<double>("Reparam.MaxGrad"), P.get<double>("Reparam.MinGrad"), P.get<int>("Levelset.VolCorrection"), Vol);

    // for serialization of geometry and numerical data
    if (P.get("Transp.DoTransp", 0))
      std::cout << "WARNING: mass transport data is not serialized, yet!" << std::endl;

    // Initialize Ensight6 output
    //Update c from ct
    //massTransp.TransformWithScaling(massTransp.ct, massTransp.c, 1.0/massTransp.GetHenry(true), 1.0/massTransp.GetHenry(false));
    std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));

    Ensight6OutCL ensight( P.get<std::string>("Ensight.EnsCase") + ".case", (P.get("Ensight.EnsightOut", 0) ? P.get<int>("Time.NumSteps")/P.get("Ensight.EnsightOut", 0)+1 : 0), P.get<int>("Ensight.Binary"), P.get<int>("Ensight.MasterOut"));
    ensight.Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(),   P.get<std::string>("Ensight.GeomName"),     ensf + ".geo", true));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
    if (P.get("Transp.DoTransp", 0)) {
        ensight.Register( make_Ensight6Scalar( massTransp.GetSolution( massTransp.ct,true),
                                                                        "TransConc",     ensf + ".ct",  true));
        ensight.Register( make_Ensight6P1XScalar( MG, lset.Phi, massTransp.ct, "XTransConcentration",   ensf + ".xconc", true));
    }
    
    // writer for vtk-format
    VTKOutCL vtkwriter(adap.GetMG(), "DROPS data", (P.get("VTK.VTKOut", 0) ? P.get<int>("Time.NumSteps")/P.get("VTK.VTKOut", 0)+1 : 0),
                P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"), P.get<int>("VTK.Binary"));

    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );

    if (P.get("Transp.DoTransp", 0)) {
        vtkwriter.Register( make_VTKScalar( massTransp.GetSolution( massTransp.ct,false), "TransConcentration") );
        vtkwriter.Register( make_VTKScalar( massTransp.GetSolution( c_out,false), "XConcentrationPos") );
        vtkwriter.Register( make_VTKScalar( massTransp.GetSolution( c_in,false), "XConcentrationNeg") );
    }

    if (P.get("Ensight.EnsightOut", 0))
        ensight.Write(0);
    if (P.get("VTK.VTKOut", 0))
        vtkwriter.Write(0);
        
//		massTransp.CheckSolution(Solutioncneg,Solutioncpos,0);
    double cmean_old = massTransp.MeanDropConcentration();

    for (int step= 1; step<=P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "============================================================ step " << std::setw(8) << step << "  /  " << std::setw(8) << P.get<int>("Time.NumSteps") << std::endl;
        double c_mean = massTransp.MeanDropConcentration();
        double t= P.get<double>("Time.StepSize") * step;
        std::cout << "Mean concentration in drop: " << c_mean <<"\n";
        if(step > 5 && std::abs(cmean_old - c_mean)/(std::abs(cmean_old) * P.get<double>("Time.StepSize")) < 1e-5 ){
          std::cout << "I think I found a stationary solution! " << std::endl;
          break;
        }
        cmean_old = c_mean;

        // grid modification
        bool doGridMod= P.get<int>("AdaptRef.Freq") && step%P.get<int>("AdaptRef.Freq") == 0;
        if (doGridMod) {
            adap.UpdateTriang( lset);
        }
      
        massTransp.DoStep( t);
        
        massTransp.GetSolutionOnPart( c_out, true , false);
        massTransp.GetSolutionOnPart( c_in, false , false);

        bool ensightoutnow = P.get("Ensight.EnsightOut", 0) && step%P.get("Ensight.EnsightOut", 0)==0;
        bool vtkoutnow = P.get("VTK.VTKOut", 0) && (step%P.get("VTK.VTKOut", 0)==0 || step < 20);
        if (ensightoutnow)
            ensight.Write(t);
        if (vtkoutnow)
            vtkwriter.Write(t);
    }
    std::cout << std::endl;

    delete pBnd_c;
    delete pBnd_ct;
}

void Strategy( InstatNavierStokes2PhaseP2P1CL& Stokes,  LsetBndDataCL& lsetbnddata, AdapTriangCL& adap)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();

    InVecMap & tdvectormap = InVecMap::getInstance();
    InScaMap & tdscalarmap = InScaMap::getInstance();
    ScaMap & scalarmap = ScaMap::getInstance();

    instat_vector_fun_ptr Flowfield = tdvectormap["ZeroVel"];
    instat_scalar_fun_ptr Reaction = tdscalarmap["ReactionFct"];
    instat_scalar_fun_ptr Rhs = tdscalarmap["Rhs"]; 
    instat_scalar_fun_ptr Initialcneg = tdscalarmap["IniCnegFct"];
    instat_scalar_fun_ptr Initialcpos = tdscalarmap["IniCposFct"];

//    scalar_fun_ptr distance = scalarmap["distance"];
    scalar_fun_ptr distance = scalarmap[P.get("Transp.Levelset", std::string("Ellipsoid"))];
    InflowLset::Init(distance);

    cBndDataCL *pBnd_c, *pBnd_ct;
    DROPS::BuildBoundaryData( &MG, pBnd_c,  P.get<std::string>("Transp.BoundaryType","2!2!2!2!2!2"), P.get<std::string>("Transp.BoundaryFncs","Dirichlet!Dirichlet!Dirichlet!Dirichlet!Dirichlet!Dirichlet"));
    DROPS::BuildBoundaryData( &MG, pBnd_ct, P.get<std::string>("Transp.BoundaryType","2!2!2!2!2!2"), P.get<std::string>("Transp.BoundaryFncs","Dirichlett!Dirichlett!Dirichlett!Dirichlett!Dirichlett!Dirichlett"));
    cBndDataCL & Bnd_c(*pBnd_c);
    cBndDataCL & Bnd_ct(*pBnd_ct); 
    
    typedef InstatNavierStokes2PhaseP2P1CL StokesProblemT;


    // initialization of surface tension
    sigma= Stokes.GetCoeff().SurfTens;
    //todo: weg oder fallunterscheidung einfuehren
    eps= P.get<double>("SurfTens.JumpWidth");    lambda= P.get<double>("SurfTens.RelPos");    sigma_dirt_fac= P.get<double>("SurfTens.DirtFactor");
    instat_scalar_fun_ptr sigmap  = 0;
    if (P.get<int>("SurfTens.VarTension"))
    {
        sigmap  = &sigma_step;
    }
    else
    {
        sigmap  = &sigmaf;
    }
    double cp=0., coeffC[5];
    //coefficients for ansatz of var. surface tension
//    coeffC[0]= 1.625; coeffC[1]= 0.0; coeffC[2]= 0.0; coeffC[3]= coeffC[4]= 0.;
    coeffC[0]= 1.625; coeffC[1]= -28.07768; coeffC[2]= 222.7858; coeffC[3]= coeffC[4]= 0.;
    
    SurfaceTensionCL sf( sigmap, Bnd_c);
    sf.SetCoeff(coeffC, cp);
    LevelsetP2CL lset( MG, lsetbnddata, sf, P.get<double>("Levelset.SD"), P.get<double>("Levelset.CurvDiff"));
    // levelset wrt the previous time step:
    LevelsetP2CL oldlset( MG, lsetbnddata, sf, P.get<double>("Levelset.SD"), P.get<double>("Levelset.CurvDiff"));
    //Prolongate and Restrict solution vector levelset from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    LevelsetRepairCL oldlsetrepair( oldlset);
    adap.push_back( &oldlsetrepair);
    //Prolongate and Restrict solution vector for the velocity from old mesh to new mesh after mesh adaptation:
    //always act on the same grid with possibly different interface position
    VelocityRepairCL velrepair( Stokes);
    adap.push_back( &velrepair);
    PressureRepairCL prrepair( Stokes, lset);
    adap.push_back( &prrepair);
    IdxDescCL* lidx= &lset.idx;
    // index wrt the interface at previous time step
    IdxDescCL* oldlidx= &oldlset.idx; 
    MLIdxDescCL* vidx= &Stokes.vel_idx;
    IdxDescCL old_vidx(vecP2_FE);
    MLIdxDescCL* pidx= &Stokes.pr_idx;

    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    VecDescCL old_v(&old_vidx);
    
    //Prolongate and Restrict solution vector old_v from old mesh to new mesh after mesh adaptation:
    DROPS::VelTranspRepairCL old_vrepair(old_v, MG, Stokes.GetBndData().Vel, old_vidx, lset.Phi, lsetbnddata, 0.);
    adap.push_back( &old_vrepair);
    oldlset.CreateNumbering( MG.GetLastLevel(), oldlidx);
    oldlset.Phi.SetIdx( oldlidx);
    
    if (P.get<int>("SurfTens.VarTension"))
        lset.SetSurfaceForce( SF_ImprovedLBVar);
    else
        lset.SetSurfaceForce( SF_ImprovedLB);

    if ( StokesSolverFactoryHelperCL().VelMGUsed(P))
        Stokes.SetNumVelLvl ( Stokes.GetMG().GetNumLevel());
    if ( StokesSolverFactoryHelperCL().PrMGUsed(P))
        Stokes.SetNumPrLvl  ( Stokes.GetMG().GetNumLevel());

    SetInitialLevelsetConditions( lset, MG, P);
    SetInitialLevelsetConditions( oldlset, MG, P);
    Stokes.CreateNumberingVel( MG.GetLastLevel(), vidx);
    Stokes.CreateNumberingPr(  MG.GetLastLevel(), pidx, 0, &lset);
    old_vidx.CreateNumbering( MG.GetLastLevel(), MG, Stokes.GetBndData().Vel);
    old_v.SetIdx  ( &old_vidx);

    Stokes.SetIdx();
    Stokes.v.SetIdx  ( vidx);
    Stokes.p.SetIdx  ( pidx);
    Stokes.InitVel( &Stokes.v, Flowfield);
    SetInitialConditions( Stokes, lset, MG, P);
    InitVel(MG, old_v, Flowfield);
   
    lset.Init( distance);
    oldlset.Init( distance);

    DisplayDetailedGeom( MG);
    DisplayUnks(Stokes, lset, MG);
   
    const double Vol= EllipsoidCL::GetVolume();
    std::cout << "initial volume: " << lset.GetVolume()/Vol << std::endl;
    double dphi= lset.AdjustVolume( Vol, 1e-9);
    std::cout << "initial volume correction is " << dphi << std::endl;
    lset.Phi.Data+= dphi;
    oldlset.Phi.Data+= dphi;
    std::cout << "new initial volume: " << lset.GetVolume()/Vol << std::endl;
   
    VelocityContainer vel(Stokes.v,Stokes.GetBndData().Vel,MG);
    //VelocityContainer vel(Flowfield);   
    TransportP1XCL massTransp( MG, Bnd_c, Bnd_ct, vel, lsetbnddata, lset.Phi, oldlset.Phi,P,0,Reaction,Rhs);
    TransportXRepairCL transprepair(massTransp, MG.GetLastLevel());
    
    // index of the concentration wrt the interface at actual time step:
    MLIdxDescCL* cidx= &massTransp.idx;
    
    // index of the concentration wrt the interface at previous time step:
    MLIdxDescCL* cidx_old= &massTransp.oldidx; 
    
    //This following Vector c_out is responsable for the communication from concentration to surface tension. 
    //Before a new surface tension is computed c_out should be updated (via GetSolutionOnPart)
    //Important: This vector has to be kept in memory as long as the surface tension is computed!
    VecDescCL c_out;
    //c_in: s. c_out but c_in is only used for visualization until now
    VecDescCL c_in;
    IdxDescCL p1idx(P1_FE,Bnd_c,0);
    /// \todo for periodic stuff: matching function here
    p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, 0, &lset.Phi,&lsetbnddata);   
    c_in.SetIdx( &p1idx);
    c_out.SetIdx( &p1idx);
    
    if (P.get("Transp.DoTransp", 0)) {
        adap.push_back(&transprepair);
        massTransp.CreateNumbering( MG.GetLastLevel(), cidx, cidx_old, lset.Phi, oldlset.Phi);
        massTransp.ct.SetIdx( cidx);
        massTransp.c.SetIdx( cidx);
        massTransp.oldct.SetIdx( cidx_old);
        if (P.get<int>("DomainCond.InitialCond") != -1)
          massTransp.Init( Initialcpos, Initialcpos, true);
        else
          ReadFEFromFile( massTransp.ct, MG, P.get<std::string>("DomainCond.InitialFile")+"concentrationTransf");
        
        std::cout << massTransp.ct.Data.size() << " concentration unknowns,\n";
 
        if (P.get<int>("SurfTens.VarTension")){
            massTransp.GetSolutionOnPart( c_out, true , false);
            massTransp.GetSolutionOnPart( c_in, false , false);
//            P1XtoP1 (*massTransp.c.RowIdx, massTransp.c.Data, p1idx, c_out.Data, c_in.Data, lset.Phi, MG);         
            sf.SetConcentration(&c_out);
            sf.SetInputMethod(Sigma_C);
            sf.SetTime(0.);
        }

        if (P.get<int>("DomainCond.InitialCond") != -1)
          massTransp.Init( Initialcneg, Initialcpos, true);
        else
          ReadFEFromFile( massTransp.ct, MG, P.get<std::string>("DomainCond.InitialFile")+"concentrationTransf");
        massTransp.GetSolutionOnPart( c_out, true , false);
        massTransp.GetSolutionOnPart( c_in, false , false);
//        P1XtoP1 (*massTransp.c.RowIdx, massTransp.c.Data, p1idx, c_out.Data, c_in.Data, lset.Phi, MG);         
          
        double c_mean = massTransp.MeanDropConcentration();
        std::cout << "START:: Mean concentration in drop: " << std::setprecision(12) << c_mean <<"\n";        
    }
    /// \todo rhs beruecksichtigen
    SurfactantcGP1CL surfTransp( MG, Stokes.GetBndData().Vel, P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"), &Stokes.v, lset.Phi, lset.GetBndData(), P.get<double>("Time.StepSize"), P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"), P.get<double>("SurfTransp.OmitBound"));
    InterfaceP1RepairCL surf_repair( MG, lset.Phi, lset.GetBndData(), surfTransp.ic);
    if (P.get("SurfTransp.DoTransp", 0))
    {
        adap.push_back( &surf_repair);
        surfTransp.idx.CreateNumbering( MG.GetLastLevel(), MG, &lset.Phi, &lset.GetBndData());
        std::cout << "Surfactant transport: NumUnknowns: " << surfTransp.idx.NumUnknowns() << std::endl;
        surfTransp.ic.SetIdx( &surfTransp.idx);
        surfTransp.Init( &surf_sol);
    }

    // Stokes-Solver
    StokesSolverFactoryCL<InstatNavierStokes2PhaseP2P1CL> stokessolverfactory(Stokes, P);
    StokesSolverBaseCL* stokessolver = stokessolverfactory.CreateStokesSolver();
//     StokesSolverAsPreCL pc (*stokessolver1, 1);
//     GCRSolverCL<StokesSolverAsPreCL> gcr(pc, C.stk_OuterIter, C.stk_OuterIter, C.stk_OuterTol, /*rel*/ false);
//     BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> >* stokessolver =
//             new BlockMatrixSolverCL<GCRSolverCL<StokesSolverAsPreCL> > (gcr);

    // Navier-Stokes-Solver
    NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>* navstokessolver = 0;
    if (P.get("NavStokes.Nonlinear", 0.0)==0.0)
        navstokessolver = new NSSolverBaseCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver);
    else
        navstokessolver = new AdaptFixedPtDefectCorrCL<InstatNavierStokes2PhaseP2P1CL>(Stokes, *stokessolver, P.get<int>("NavStokes.Iter"), P.get<double>("NavStokes.Tol"), P.get<double>("NavStokes.Reduction"));
    // Level-Set-Solver
#ifndef _PAR
    SSORPcCL ssorpc;
    GMResSolverCL<SSORPcCL>* gm = new GMResSolverCL<SSORPcCL>( ssorpc, 100, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"));
#else
    ParJac0CL jacparpc( *lidx);
    ParPreGMResSolverCL<ParJac0CL>* gm = new ParPreGMResSolverCL<ParJac0CL>
           (/*restart*/100, P.get<int>("Levelset.Iter"), P.get<double>("Levelset.Tol"), *lidx, jacparpc,/*rel*/true, /*acc*/ true, /*modGS*/false, LeftPreconditioning, /*parmod*/true);
#endif

    LevelsetModifyCL lsetmod( P.get<int>("Reparam.Freq"), P.get<int>("Reparam.Method"), P.get<double>("Reparam.MaxGrad"), P.get<double>("Reparam.MinGrad"), P.get<int>("Levelset.VolCorrection"), Vol);

    // Time discretisation + coupling
    TimeDisc2PhaseCL* timedisc= CreateTimeDisc(Stokes, lset, navstokessolver, gm, P, lsetmod);

    if (P.get<int>("Time.NumSteps") != 0){
        timedisc->SetTimeStep( P.get<double>("Time.StepSize"));
        timedisc->SetSchurPrePtr( stokessolverfactory.GetSchurPrePtr() );
    }

    if (P.get("NavStokes.Nonlinear", 0.0)!=0.0 || P.get<int>("Time.NumSteps") == 0) {
        stokessolverfactory.SetMatrixA( &navstokessolver->GetAN()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( navstokessolver->GetAN(), &Stokes.B.Data,
                                         &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
    }
    else {
        stokessolverfactory.SetMatrixA( &timedisc->GetUpperLeftBlock()->GetFinest());
            //for Stokes-MGM
        stokessolverfactory.SetMatrices( timedisc->GetUpperLeftBlock(), &Stokes.B.Data,
                                         &Stokes.M.Data, &Stokes.prM.Data, &Stokes.pr_idx);
    }

    UpdateProlongationCL PVel( Stokes.GetMG(), stokessolverfactory.GetPVel(), &Stokes.vel_idx, &Stokes.vel_idx);
    adap.push_back( &PVel);
    UpdateProlongationCL PPr ( Stokes.GetMG(), stokessolverfactory.GetPPr(), &Stokes.pr_idx, &Stokes.pr_idx);
    adap.push_back( &PPr);
    // For a two-level MG-solver: P2P1 -- P2P1X;
//     MakeP1P1XProlongation ( Stokes.vel_idx.NumUnknowns(), Stokes.pr_idx.NumUnknowns(),
//         Stokes.pr_idx.GetFinest().GetXidx().GetNumUnknownsStdFE(),
//         stokessolverfactory.GetPVel()->GetFinest(), stokessolverfactory.GetPPr()->GetFinest());

    std::ofstream* infofile = 0;
    IF_MASTER {
        infofile = new std::ofstream ((P.get<std::string>("VTK.VTKName","ns_transp")+".info").c_str());
    }
    IFInfo.Init(infofile);
    IFInfo.WriteHeader();

    if (P.get<int>("Time.NumSteps") == 0)
        SolveStatProblem( Stokes, lset, *navstokessolver);

    // for serialization of geometry and numerical data
    if (P.get("Transp.DoTransp", 0))
      std::cout << "WARNING: mass transport data is not serialized, yet!" << std::endl;
    TwoPhaseStoreCL<InstatNavierStokes2PhaseP2P1CL> ser(MG, Stokes, lset, 0, P.get<std::string>("Restart.Outputfile"), P.get<int>("Restart.Overwrite"), P.get<int>("Restart.Binary"));




    // Output-Registrations:
    Ensight6OutCL* ensight = NULL;
    if (P.get<int>("Ensight.EnsightOut",0)){
        // Initialize Ensight6 output
        std::string ensf( P.get<std::string>("Ensight.EnsDir") + "/" + P.get<std::string>("Ensight.EnsCase"));
        ensight = new Ensight6OutCL( P.get<std::string>("Ensight.EnsCase") + ".case", 
                                     P.get<int>("Time.NumSteps")/P.get("Ensight.EnsightOut", 0)+1,
                                     P.get<int>("Ensight.Binary"), P.get<int>("Ensight.MasterOut"));
        ensight->Register( make_Ensight6Geom      ( MG, MG.GetLastLevel(), P.get<std::string>("Ensight.GeomName"),
                                                    ensf + ".geo", true));
        ensight->Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));
        ensight->Register( make_Ensight6Scalar    ( Stokes.GetPrSolution(),  "Pressure",      ensf + ".pr",  true));
        ensight->Register( make_Ensight6Vector    ( Stokes.GetVelSolution(), "Velocity",      ensf + ".vel", true));
        ensight->Register( make_Ensight6Scalar    ( ScalarFunAsP2EvalCL( sigmap, 0., &MG, MG.GetLastLevel()),
                                                    "Surfaceforce",  ensf + ".sf",  true));

        if (P.get("Transp.DoTransp", 0)) {
            ensight->Register( make_Ensight6Scalar( massTransp.GetSolution(),"Concentration", ensf + ".c",   true));
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, massTransp.ct, 
                                                      "XTransConcentration",   ensf + ".xconc", true));

        }
        if (P.get("SurfTransp.DoTransp", 0)) {
            ensight->Register( make_Ensight6IfaceScalar( MG, surfTransp.ic,  "InterfaceSol",  ensf + ".sur", true));
        }

#ifndef _PAR
        if (Stokes.UsesXFEM())
            ensight->Register( make_Ensight6P1XScalar( MG, lset.Phi, Stokes.p, "XPressure",   ensf + ".pr", true));
#endif
        ensight->Write( Stokes.v.t);
    }

    // writer for vtk-format
    VTKOutCL * vtkwriter = NULL;
    if (P.get<int>("VTK.VTKOut",0)){
        vtkwriter = new VTKOutCL(adap.GetMG(), "DROPS data", 
                                 P.get<int>("Time.NumSteps")/P.get("VTK.VTKOut", 0)+1,
                                 P.get<std::string>("VTK.VTKDir"), P.get<std::string>("VTK.VTKName"), 
                                 P.get<int>("VTK.Binary"));
        vtkwriter->Register( make_VTKVector( Stokes.GetVelSolution(), "velocity") );
        vtkwriter->Register( make_VTKScalar( Stokes.GetPrSolution(), "pressure") );
        vtkwriter->Register( make_VTKScalar( lset.GetSolution(), "level-set") );

        if (P.get("Transp.DoTransp", 0)) {
            vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( massTransp.ct,false), "TransConcentration") );
            vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_out,false), "XConcentrationPos") );
            vtkwriter->Register( make_VTKScalar( massTransp.GetSolution( c_in,false), "XConcentrationNeg") );

        }

        if (P.get("SurfTransp.DoTransp", 0)) {
            vtkwriter->Register( make_VTKIfaceScalar( MG, surfTransp.ic,  "InterfaceSol"));
        }
        vtkwriter->Write(Stokes.v.t);
    }

    for (int step= 1; step<=P.get<int>("Time.NumSteps"); ++step)
    {
        std::cout << "============================================================ step " << step << std::endl;
        double c_mean = massTransp.MeanDropConcentration();
        double t= Stokes.v.t;
        std::cout << "Mean concentration in drop: " << c_mean <<"\n";
        
        IFInfo.Update( lset, Stokes.GetVelSolution());
//         IFInfo.Write(Stokes.t, c_mean);
        IFInfo.Write(Stokes.v.t);
        
        if (P.get("SurfTransp.DoTransp", 0)) surfTransp.InitOld();
        timedisc->DoStep( P.get<int>("Coupling.Iter"));
        
//         if (P.get("Transp.DoTransp", 0)) massTransp.DoStep( step*P.get<double>("Time.StepSize"));
        if (P.get("SurfTransp.DoTransp", 0)) {
            surfTransp.DoStep( step*P.get<double>("Time.StepSize"));
            BndDataCL<> ifbnd( 0);
            std::cout << "surfactant on \\Gamma: " << Integral_Gamma( MG, lset.Phi, lset.GetBndData(), make_P1Eval(  MG, ifbnd, surfTransp.ic)) << '\n';
        }
 
        // WriteMatrices( Stokes, step);

        // grid modification
        bool doGridMod= P.get<int>("AdaptRef.Freq") && step%P.get<int>("AdaptRef.Freq") == 0;
        if (doGridMod) {
            adap.UpdateTriang( lset);
        }
      
        if (P.get("Transp.DoTransp", 0)) {
            massTransp.DoStep( t);
            old_vrepair.SetTime(t);\
            
        //}

            if (P.get<int>("SurfTens.VarTension")){
                sf.SetConcentration(&c_out);
                sf.SetInputMethod(Sigma_C);
                sf.SetTime(t);
            }
        }
        
        timedisc->Update(); 
        
        //Update c from ct
//        massTransp.TransformWithScaling(massTransp.ct, massTransp.c, 1.0/massTransp.GetHenry(true), 1.0/massTransp.GetHenry(false));
		
		massTransp.GetSolutionOnPart( c_out, true , false);
		massTransp.GetSolutionOnPart( c_in, false , false);

/// \todo for periodic stuff: matching function here
//        p1idx.CreateNumbering( MG.GetLastLevel(), MG, Bnd_c, 0, &lset.Phi,&lsetbnddata);   
//        P1XtoP1 (*massTransp.c.RowIdx, massTransp.c.Data, p1idx, c_out.Data, c_in.Data, lset.Phi, MG);         
        
        if (ensight && step%P.get("Ensight.EnsightOut", 0)==0)
            ensight->Write( Stokes.v.t);
        if (vtkwriter && step%P.get("VTK.VTKOut", 0)==0)
            vtkwriter->Write(Stokes.v.t);
//        if (C.rst_Serialization && step%C.rst_Serialization==0)
//            ser.Write();
    }
    IFInfo.Update( lset, Stokes.GetVelSolution());
    IFInfo.Write(Stokes.v.t);
    std::cout << std::endl;
    delete timedisc;
    delete navstokessolver;
    delete stokessolver;
    delete gm;
    if (infofile) delete infofile;
//     delete stokessolver1;
    if (vtkwriter) delete vtkwriter;
    if (ensight) delete ensight;
    delete pBnd_c;
    delete pBnd_ct;

}

} // end of namespace DROPS

int main (int argc, char** argv)
{
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
#endif
  try
  {
#ifdef _PAR
    DROPS::ParMultiGridInitCL pmginit;
#endif
    std::ifstream param;
    if (argc!=2)
    {
        std::cout << "Using default parameter file: risingbutanoldroplet.json\n";
        param.open( "risingbutanoldroplet.json");
    }
    else
        param.open( argv[1]);
    if (!param)
    {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> P;
    param.close();
    std::cout << P << std::endl;

    DROPS::MultiGridCL* mg= 0;
    DROPS::StokesBndDataCL* bnddata= 0;
    DROPS::LsetBndDataCL* lsetbnddata= 0;

    double RadInlet = P.get<double>("Exp.RadInlet");
    DROPS::BuildDomain( mg, P.get<std::string>("DomainCond.MeshFile"), P.get<int>("DomainCond.GeomType"), P.get<std::string>("Restart.Inputfile"), RadInlet);
    P.put("Exp.RadInlet", RadInlet);
    DROPS::BuildBoundaryData( mg, bnddata, P.get<std::string>("DomainCond.BoundaryType"), P.get<std::string>("DomainCond.BoundaryFncs"));

    // todo: reasonable implementation needed
    std::string lsetbndtype = "98" /*NoBC*/, lsetbndfun = "Zero";
    for( size_t i= 1; i<mg->GetBnd().GetNumBndSeg(); ++i) {
        lsetbndtype += "!98";
        lsetbndfun  += "!Zero";
    }

    DROPS::BuildBoundaryData( mg, lsetbnddata, lsetbndtype, lsetbndfun);

    std::cout << "Generated MG of " << mg->GetLastLevel() << " levels." << std::endl;

    DROPS::EllipsoidCL::Init(P.get<DROPS::Point3DCL>("Exp.PosDrop"), P.get<DROPS::Point3DCL>("Exp.RadDrop"));
    DROPS::AdapTriangCL adap( *mg, P.get<double>("AdaptRef.Width"), P.get<int>("AdaptRef.CoarsestLevel"), P.get<int>("AdaptRef.FinestLevel"), ((P.get<std::string>("Restart.Inputfile") == "none") ? P.get<int>("AdaptRef.LoadBalStrategy") : -P.get<int>("AdaptRef.LoadBalStrategy")), P.get<int>("AdaptRef.Partitioner"));
    // If we read the Multigrid, it shouldn't be modified;
    // otherwise the pde-solutions from the ensight files might not fit.
    if (!P.get<int>("Transp.UseNSSol",1) || (P.get<std::string>("Restart.Inputfile") == "none")){
        DROPS::ScaMap & scalarmap = DROPS::ScaMap::getInstance();
        DROPS::scalar_fun_ptr distance = scalarmap[P.get("Transp.Levelset", std::string("Ellipsoid"))];
        adap.MakeInitialTriang( distance);
    }
    
    std::cout << DROPS::SanityMGOutCL(*mg) << std::endl;
#ifdef _PAR
    adap.GetLb().GetLB().SetWeightFnct(3);
    if (DROPS::ProcCL::Check( CheckParMultiGrid( adap.GetPMG())))
        std::cout << "As far as I can tell the ParMultigridCl is sane\n";
#endif

    if (P.get<int>("Transp.UseNSSol",1)){
      DROPS::InstatNavierStokes2PhaseP2P1CL prob( *mg, DROPS::TwoPhaseFlowCoeffCL(P), *bnddata, P.get<double>("Stokes.XFEMStab")<0 ? DROPS::P1_FE : DROPS::P1X_FE, P.get<double>("Stokes.XFEMStab"));
      Strategy( prob, *lsetbnddata, adap);    // do all the stuff
    }
    else
    {
      OnlyTransportStrategy( *mg, *lsetbnddata, adap);    // do just the transport stuff
    }

    delete mg;
    delete lsetbnddata;
    delete bnddata;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
