/// \file
/// \brief Solve a non-stationary convection-diffusion-equation on a moving interface
/// \author LNM RWTH Aachen: Joerg Grande

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

#include "surfactant/ifacetransp.h"
#include "misc/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/surfacetension.h"
#include "out/ensightOut.h"

#include <fstream>

using namespace DROPS;

DROPS::ParamCL P;
std::string ensf; // basename of the ensight files

DROPS::Point3DCL u_func (const DROPS::Point3DCL&, double)
{
    return P.get<DROPS::Point3DCL>("Exp.Velocity");
}

typedef DROPS::Point3DCL (*bnd_val_fun) (const DROPS::Point3DCL&, double);

DROPS::BndCondT bc[6]= {
    DROPS::DirBC, DROPS::DirBC,
    DROPS::DirBC, DROPS::DirBC,
    DROPS::DirBC, DROPS::DirBC
};

bnd_val_fun bf[6]= {
    &u_func, &u_func, &u_func, &u_func, &u_func, &u_func
};

double sphere_2 (const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL x( p - P.get<DROPS::Point3DCL>("Exp.PosDrop"));

    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (P.get<DROPS::Point3DCL>("Exp.PosDrop") + t*u_func(p, t)));

    return x.norm() - P.get<DROPS::Point3DCL>("Exp.RadDrop")[0];
}


// TestCase == 0: Sphere around 0, RadDrop 1, v == 0
// A right hand side from C.J. Heine...
const double a( -13./8.*std::sqrt( 35./M_PI));
double rhs0 (const DROPS::Point3DCL& p, double)
{
    return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
// ...and the corresponding solution (extended)
double sol0 (const DROPS::Point3DCL& p, double)
{
//    return p.norm_sq()/(12. + p.norm_sq())*rhs0( p, 0.);
    return 1. + p.norm_sq()/(12. + p.norm_sq())*rhs0( p, 0.);
}

double sol0t (const DROPS::Point3DCL& p, double t)
{
    const Point3DCL q( p - (P.get<DROPS::Point3DCL>("Exp.PosDrop") + t*u_func(p, t)));
    const double val( a*(3.*q[0]*q[0]*q[1] - q[1]*q[1]*q[1]));

//    return q.norm_sq()/(12. + q.norm_sq())*val;
    return 1. + q.norm_sq()/(12. + q.norm_sq())*val;
}

template<class DiscP1FunType>
double L2_error (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    const DiscP1FunType& discsol, DROPS::instat_scalar_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qsol, qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
            	triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    qdiscsol.assign(  *it, &triangle.GetBary( tri), discsol);
                    d+= DROPS::Quad5_2DCL<>( std::pow( qdiscsol - qsol, 2)).quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls, const BndDataCL<>& lsbnd,
    DROPS::instat_scalar_fun_ptr extsol, double t= 0.)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfaceTriangleCL triangle;
    DROPS::Quad5_2DCL<> qsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
    	triangle.Init( *it, ls, lsbnd);
        if (triangle.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
            	triangle.ComputeForChild( ch);
                for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &triangle.GetBary( tri), extsol, t);
                    d+= DROPS::Quad5_2DCL<>( qsol*qsol).quad( triangle.GetAbsDet( tri));
                }
            }
        }
    }
    return std::sqrt( d);
}

void LinearLSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, DROPS::scalar_fun_ptr d)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord());

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= ls.Data[it->Unknowns( idx)]=
            0.5*(ls.Data[it->GetVertex( 0)->Unknowns( idx)] + ls.Data[it->GetVertex( 1)->Unknowns( idx)]);
}

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
}

void InitVel ( const MultiGridCL& mg, VecDescCL* vec, BndDataCL<Point3DCL>& Bnd, instat_vector_fun_ptr LsgVel, double t)
{
    VectorCL& lsgvel= vec->Data;
    const Uint lvl  = vec->GetLevel(),
               vidx = vec->RowIdx->GetIdx();

   DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel(sit->GetCoord(), t));
    }
    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, sit) {
        if (!Bnd.IsOnDirBnd( *sit))
            DoFHelperCL<Point3DCL, VectorCL>::set( lsgvel, sit->Unknowns( vidx),
                LsgVel( (sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2., t));
    }
    vec->t = t;
}

void P1Init (instat_scalar_fun_ptr icf, VecDescCL& ic, MultiGridCL& mg, double t)
{
    const Uint lvl= ic.GetLevel(),
               idx= ic.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {
        if (it->Unknowns.Exist( idx))
            ic.Data[it->Unknowns( idx)]= icf( it->GetCoord(), t);
    }
}


void Strategy (DROPS::MultiGridCL& mg, DROPS::AdapTriangCL& adap, DROPS::LevelsetP2CL& lset,
    DROPS::Ensight6OutCL& ensight)
{
    using namespace DROPS;

    LSInit( mg, lset.Phi, &sphere_2move, 0.);

    BndDataCL<> lsetbnd2( 6);
    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    DROPS::LevelsetP2CL lset2( mg, lsetbnd2, sf, P.get<double>("Levelset.Theta"), P.get<double>("Levelset.SD")); // Only for output
    lset2.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset2.Phi.SetIdx( &lset2.idx);
    LSInit( mg, lset2.Phi, &sphere_2move, 0.);

    const double Vol= lset.GetVolume();
    std::cout << "droplet volume: " << Vol << std::endl;

    BndDataCL<Point3DCL> Bnd_v( 6, bc, bf);
    IdxDescCL vidx( vecP2_FE);
    vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( mg, &v, Bnd_v, u_func, 0.);

    lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.StepSize"));

    SurfactantcGP1CL timedisc( mg, Bnd_v, P.get<double>("SurfTransp.Theta"), P.get<double>("SurfTransp.Visc"), &v, lset.Phi, lset.GetBndData(), P.get<double>("Time.StepSize"), P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"), P.get<double>("SurfTransp.OmitBound"));

    LevelsetRepairCL lsetrepair( lset);
    adap.push_back( &lsetrepair);
    InterfaceP1RepairCL ic_repair( mg, lset.Phi, lset.GetBndData(), timedisc.ic);
    adap.push_back( &ic_repair);
    LevelsetRepairCL lset2repair( lset2);
    adap.push_back( &lset2repair);

    // Init Interface-Sol
    timedisc.idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << timedisc.idx.NumUnknowns() << std::endl;
    timedisc.ic.SetIdx( &timedisc.idx);
    timedisc.Init( &sol0);
    // timedisc.Update();

    // Additional Ensight6-variables
    ensight.Register( make_Ensight6IfaceScalar( mg, timedisc.ic,                 "InterfaceSol", ensf + ".sur", true));
    ensight.Register( make_Ensight6Vector(      make_P2Eval( mg, Bnd_v, v),      "Velocity",     ensf + ".vel",  true));
    ensight.Register( make_Ensight6Scalar(      lset2.GetSolution(),             "Levelset2",    ensf + ".scl2", true));
    ensight.Register( make_Ensight6Scalar( ScalarFunAsP2EvalCL( sol0t, 0., &mg), "TrueSol",      ensf + ".sol", true));
    ensight.Write( 0.);

    // timedisc.SetTimeStep( P.get<double>("Time.StepSize"), P.get<double>("SurfTransp.Theta"));
//    std::cout << "L_2-error: " << L2_error( mg, lset.Phi, timedisc.GetSolution(), &sol0t, 0.)
//              << " norm of true solution: " << L2_norm( mg, lset.Phi, &sol0t, 0.)
//              << std::endl;
    BndDataCL<> ifbnd( 0);
    std::cerr << "initial surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';

    for (int step= 1; step <= P.get<int>("Time.NumSteps"); ++step) {
        std::cout << "======================================================== step " << step << ":\n";

        timedisc.InitOld();
        LSInit( mg, lset.Phi, &sphere_2move, step*P.get<double>("Time.StepSize"));
        timedisc.DoStep( step*P.get<double>("Time.StepSize"));
        std::cerr << "surfactant on \\Gamma: " << Integral_Gamma( mg, lset.Phi, lset.GetBndData(), make_P1Eval(  mg, ifbnd, timedisc.ic)) << '\n';

        //lset2.DoStep();
//        VectorCL rhs( lset2.Phi.Data.size());
//        lset2.ComputeRhs( rhs);
//        lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v, step*P.get<double>("Time.StepSize")));
//        lset2.SetTimeStep( P.get<double>("Time.StepSize"));
//        lset2.DoStep( rhs);

        std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        if (P.get("Levelset.VolCorr", 0)) {
            double dphi= lset.AdjustVolume( Vol, 1e-9);
            std::cout << "volume correction is " << dphi << std::endl;
            lset.Phi.Data+= dphi;
            std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
        }
        //if (C.rpm_Freq && step%C.rpm_Freq==0) { // reparam levelset function
            // lset.ReparamFastMarching( C.rpm_Method);
        if (P.get("AdaptRef.Freq", 0) != 0 && step%P.get("AdaptRef.Freq", 0) == 0) {
            if (P.get("AdaptRef.Freq", 0) != 0) {
                adap.UpdateTriang( lset);
                vidx.DeleteNumbering( mg);
                vidx.CreateNumbering( mg.GetLastLevel(), mg, Bnd_v);
                v.SetIdx( &vidx);
                InitVel( mg, &v, Bnd_v, u_func, step*P.get<double>("Time.StepSize"));
                LSInit( mg, lset.Phi, &sphere_2move, step*P.get<double>("Time.StepSize"));
                timedisc.Update();

                lset2.SetupSystem( make_P2Eval( mg, Bnd_v, v), P.get<double>("Time.StepSize"));
            }
            std::cout << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            if (P.get("Levelset.VolCorr", 0)) {
                double dphi= lset.AdjustVolume( Vol, 1e-9);
                std::cout << "volume correction is " << dphi << std::endl;
                lset.Phi.Data+= dphi;
                std::cout << "new rel. Volume: " << lset.GetVolume()/Vol << std::endl;
            }
        }
        ensight.Write( step*P.get<double>("Time.StepSize"));
//        std::cout << "L_2-error: " << L2_error( mg, lset.Phi, timedisc.GetSolution(), &sol0t, step*P.get<double>("Time.StepSize"))
//                  << " norm of true solution: " << L2_norm( mg, lset.Phi, &sol0t, step*P.get<double>("Time.StepSize"))
//                  << std::endl;
    }
    std::cout << std::endl;
}


int main (int argc, char* argv[])
{
  try {
    std::ifstream param;
    if (argc!=2) {
        std::cout << "Using default parameter file: surfactant.json\n";
        param.open( "surfactant.json");
    }
    else
        param.open( argv[1]);
    if (!param) {
        std::cout << "error while opening parameter file\n";
        return 1;
    }
    param >> P;
    param.close();
    std::cout << P << std::endl;
    ensf= P.get<std::string>("EnsightDir") + "/" + P.get<std::string>("EnsightCase"); // basename of the ensight files

    std::cout << "Setting up interface-PDE:\n";
    DROPS::BrickBuilderCL brick( DROPS::MakePoint3D( -2., -2., -2.),
                                 4.*DROPS::std_basis<3>( 1),
                                 4.*DROPS::std_basis<3>( 2),
                                 4.*DROPS::std_basis<3>( 3),
                                 P.get<int>("InitialDivisions"), P.get<int>("InitialDivisions"), P.get<int>("InitialDivisions"));
    DROPS::MultiGridCL mg( brick);
    DROPS::AdapTriangCL adap( mg, P.get<double>("AdaptRef.Width"), 0, P.get<int>("AdaptRef.FinestLevel"));
    adap.MakeInitialTriang( sphere_2);

    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

    DROPS::LevelsetP2CL lset( mg, lsbnd, sf);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    LinearLSInit( mg, lset.Phi, &sphere_2);
//    lset.Init( &sphere_2);

    // Initialize Ensight6 output
    std::string ensf( P.get<std::string>("EnsightDir") + "/" + P.get<std::string>("EnsightCase"));
    Ensight6OutCL ensight( P.get<std::string>("EnsightCase") + ".case", P.get<int>("Time.NumSteps") + 1);
    // ensight.Register( make_Ensight6Geom      ( mg, mg.GetLastLevel(),   "Geometrie",     ensf + ".geo", true));
    ensight.Register( make_Ensight6Geom      ( mg, mg.GetLastLevel(),   "Geometrie",     ensf + ".geo", false));
    ensight.Register( make_Ensight6Scalar    ( lset.GetSolution(),      "Levelset",      ensf + ".scl", true));

    if (P.get("TestCase", 0) > 0) { // Time dependent tests in Strategy
        Strategy( mg, adap, lset, ensight);
        return 0;
    }

    DROPS::IdxDescCL ifaceidx( P1IF_FE);
    ifaceidx.GetXidx().SetBound( P.get<double>("SurfTransp.OmitBound"));
    ifaceidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi, &lset.GetBndData());
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi, lset.GetBndData());
    std::cout << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi, lset.GetBndData(), P.get<double>("SurfTransp.Visc"));
    DROPS::MatrixCL L;
    L.LinComb( 1.0, A.Data, 1.0, M.Data);
    DROPS::VecDescCL b( &ifaceidx);
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, lset.GetBndData(), rhs0);

    DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    DROPS::WriteToFile( b.Data, "rhs_iface.txt", "rhs");

    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, P.get<int>("SurfTransp.Iter"), P.get<double>("SurfTransp.Tol"), true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data);
    std::cout << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    DROPS::WriteToFile( x.Data, "x_iface.txt", "solution");

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE);
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);
    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);
    DROPS::NoBndDataCL<> nobnd;
    ensight.Register( make_Ensight6Scalar( make_P1Eval( mg, nobnd, xext), "InterfaceSol", ensf + ".sur"));
    ensight.Write();

    double L2_err( L2_error( mg, lset.Phi, lset.GetBndData(), make_P1Eval( mg, nobnd, xext), &sol0));
    std::cout << "L_2-error: " << L2_err << std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
