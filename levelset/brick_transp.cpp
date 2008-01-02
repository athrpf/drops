//****************************************************************************
// File:    brick_trans.cpp                                                  *
// Content: test case: two-phase flow in square pipe with transport          *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen      *
//****************************************************************************

#include "geom/multigrid.h"
#include "geom/builder.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "levelset/params.h"
#include "levelset/mzelle_hdr.h"
#include "num/bndData.h"
#include "poisson/transport2phase.h"
#include <fstream>

DROPS::ParamMesszelleCL C;

DROPS::SVectorCL<3> Inflow (const DROPS::Point3DCL& p, double)
{
    DROPS::SVectorCL<3> ret(0.);
    const double x = p[0]*(2*C.r_inlet-p[0]) / (C.r_inlet*C.r_inlet),
                 z = p[2]*(2*C.r_inlet-p[2]) / (C.r_inlet*C.r_inlet);
    ret[1]= x * z * C.Anstroem;
    return ret;
}

double Initialcneg (const DROPS::Point3DCL&, double)
{
    return 1.0;
}

double Initialcpos (const DROPS::Point3DCL&, double)
{
    return 0.3; // 0.5;
}

typedef DROPS::BndDataCL<DROPS::Point3DCL> VelBndDataCL;
typedef VelBndDataCL::bnd_val_fun  vel_bnd_val_fun;
typedef DROPS::BndDataCL<> cBndDataCL;
typedef cBndDataCL::bnd_val_fun  c_bnd_val_fun;

const DROPS::BndCondT v_bc[6]= { DROPS::Dir0BC, DROPS::Dir0BC, DROPS::DirBC, DROPS::DirBC, DROPS::Dir0BC, DROPS::Dir0BC };
const vel_bnd_val_fun v_bfun[6]= { 0, 0, &Inflow, &Inflow, 0, 0 };

const DROPS::BndCondT c_bc[6]= { DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC, DROPS::DirBC };
const c_bnd_val_fun c_bfun[6]= { & Initialcpos,  & Initialcpos, & Initialcpos,& Initialcpos, & Initialcpos, & Initialcpos };

double D[2]= { /*pos. part*/ 10e-3, /*neg. part*/ 5e-3 };
double H= 0.5; // in the neg. part

namespace DROPS
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

typedef P2EvalCL<SVectorCL<3>, const VelBndDataCL, const VecDescCL> const_DiscVelSolCL;

void Strategy (MultiGridCL& MG)
{
    LevelsetP2CL lset( MG, &sigmaf, /*grad sigma*/ 0, C.lset_theta, C.lset_SD, -1, C.lset_iter, C.lset_tol, C.CurvDiff);
    IdxDescCL* lidx= &lset.idx;
    lset.CreateNumbering( MG.GetLastLevel(), lidx);
    lset.Phi.SetIdx( lidx);
    lset.Init( EllipsoidCL::DistanceFct);

    VelBndDataCL Bnd_v( 6, v_bc, v_bfun);
    IdxDescCL  vidx( vecP2_FE);
    CreateNumb( MG.GetLastLevel(), vidx, MG, Bnd_v);
    VecDescCL v( &vidx);
    InitVel( MG, v, &Inflow);

    cBndDataCL Bnd_c( 6, c_bc, c_bfun);
    TransportP1CL c( MG, Bnd_c, Bnd_v, /*theta*/ 0.5, D, H, &v, lset,
        /*t*/ 0., C.dt, C.outer_iter, C.outer_tol);
    IdxDescCL* cidx= &c.idx;
    c.CreateNumbering( MG.GetLastLevel(), cidx);
    c.ct.SetIdx( cidx);
    c.Init( &Initialcneg, &Initialcpos);
    c.Update();

    EnsightP2SolOutCL ensight( MG, lidx);
    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 datc = filename+".c" ,
                 datct = filename+".ct" ,
                 datvec= filename+".vel",
                 datscl= filename+".scl";
    ensight.CaseBegin( string(C.EnsCase+".case").c_str(), C.num_steps+1);
    ensight.DescribeGeom( "Messzelle", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "Concentration", datc, true);
    ensight.DescribeScalar( "TransConc", datct, true);
    ensight.DescribeVector( "Velocity", datvec, true);
    ensight.putGeom( datgeo);

    MG.SizeInfo( std::cerr);
    std::cerr << c.c.Data.size() << " concentration unknowns,\n";
    std::cerr << v.Data.size() << " velocity unknowns,\n";
    std::cerr << lset.Phi.Data.size() << " levelset unknowns.\n";

    const double Vol= EllipsoidCL::GetVolume();
    std::cerr << "rel. Volume: " << lset.GetVolume()/Vol << std::endl;

    ensight.putVector( datvec, const_DiscVelSolCL( &v, &Bnd_v, &MG, 0.), 0);
    ensight.putScalar( datc,  c.GetSolution(), 0);
    ensight.putScalar( datct,  c.GetSolution( c.ct), 0);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.Commit();

    c.SetTimeStep( C.dt);
    for (int step= 1; step <= C.num_steps; ++step) {
        std::cerr << "======================================================== Schritt " << step << ":\n";
        c.DoStep( step*C.dt);
        ensight.putScalar( datc, c.GetSolution(), step*C.dt);
        ensight.putScalar( datct, c.GetSolution( c.ct), step*C.dt);
        ensight.putVector( datvec, const_DiscVelSolCL( &v, &Bnd_v, &MG, step*C.dt), step*C.dt);
        ensight.putScalar( datscl, lset.GetSolution(), step*C.dt);
        ensight.Commit();
    }
    ensight.CaseEnd();
    std::cerr << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try {
    if (argc != 2) {
        std::cerr << "You have to specify one parameter:\n\t"
                  << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param) {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    DROPS::Point3DCL e1(0.), e2(0.), e3(0.), orig;
    e1[0]=2*C.r_inlet; e2[1]=1.0; e3[2]= 2*C.r_inlet;
    DROPS::BrickBuilderCL builder( orig, e1, e2, e3, 20, 20, 20);
    DROPS::MultiGridCL mg( builder);
    std::cerr << DROPS::SanityMGOutCL( mg) << std::endl;
    EllipsoidCL::Init( C.Mitte, C.Radius);

    Strategy( mg);    // do all the stuff

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
