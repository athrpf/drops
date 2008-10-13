#include "surfactant/ifacetransp.h"
#include "surfactant/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "out/ensightOut.h"

#include <fstream>

DROPS::ParamSurfactantCL C;

double sphere_2 (const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL x( p - C.Mitte);

    return x.norm() - C.Radius[0];
}

double sphere_2ls (const DROPS::Point3DCL& p, double)
{
    return sphere_2( p);
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

// TestCase == 0: Sphere around 0, radius 1, v == 0
// A right hand side from C.J. Heine...
const double a( -13./8.*std::sqrt( 35./M_PI));
double rhs0 (const DROPS::Point3DCL& p, double)
{
    return a*(3.*p[0]*p[0]*p[1] - p[1]*p[1]*p[1]);
}
// ...and the corresponding solution (extended)
double sol0 (const DROPS::Point3DCL& p, double)
{
    return p.norm_sq()/(12. + p.norm_sq())*rhs0( p, 0.);
}

typedef DROPS::P1EvalCL<double, const DROPS::NoBndDataCL<>, DROPS::VecDescCL> DiscP1FunT;

double L2_error (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls,
    DiscP1FunT& discsol, DROPS::instat_scalar_fun_ptr extsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfacePatchCL patch;
    DROPS::Quad5_2DCL<> qsol, qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri) {
                    qsol.assign( *it, &patch.GetBary( tri), extsol);
                    qdiscsol.assign(  *it, &patch.GetBary( tri), discsol);
                    d+= DROPS::Quad5_2DCL<>( std::pow( qdiscsol - qsol, 2)).quad( patch.GetFuncDet( tri));
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

int main (int argc, char* argv[])
{
  try {
    std::ifstream param;
    if (argc!=2) {
        std::cerr << "Using default parameter file: surfactant.param\n";
        param.open( "surfactant.param");
    }
    else
        param.open( argv[1]);
    if (!param) {
        std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    std::cerr << C << std::endl;

    std::cerr << "Setting up interface-PDE:\n";
    DROPS::BrickBuilderCL brick( DROPS::MakePoint3D( -2., -2., -2.),
                                 4.*DROPS::std_basis<3>( 1),
                                 4.*DROPS::std_basis<3>( 2),
                                 4.*DROPS::std_basis<3>( 3),
                                 C.cdiv, C.cdiv, C.cdiv);
    DROPS::MultiGridCL mg( brick);
    DROPS::AdapTriangCL adap( mg, C.ref_width, 0, C.ref_flevel);
    adap.MakeInitialTriang( sphere_2);

    DROPS::LevelsetP2CL lset( mg);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);
    LinearLSInit( mg, lset.Phi, &sphere_2);

    DROPS::IdxDescCL ifaceidx( DROPS::P1_FE);
    DROPS::CreateNumbOnInterface( mg.GetLastLevel(), ifaceidx, mg, lset.Phi);
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns << std::endl;

    DROPS::IdxDescCL ifacefullidx( DROPS::P1_FE, DROPS::NoBndDataCL<>());
    ifacefullidx.CreateNumbering( mg.GetLastLevel(), mg);

    DROPS::EnsightP2SolOutCL ensight( mg, &lset.idx);
    const std::string filename( C.EnsDir + "/" + C.EnsCase);
    const std::string datgeo=  filename+".geo",
                      datscl=  filename+".scl",
                     datsurf = filename+".sur";
    ensight.CaseBegin( std::string( C.EnsCase+".case").c_str(), 1);
    ensight.DescribeGeom( "Geometrie", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true);
    ensight.DescribeScalar( "InterfaceSol", datsurf,  true);
    ensight.putGeom( datgeo);

    DROPS::MatDescCL M( &ifaceidx, &ifaceidx);
    DROPS::SetupInterfaceMassP1( mg, &M, lset.Phi);
    std::cerr << "M is set up.\n";
    DROPS::MatDescCL A( &ifaceidx, &ifaceidx);
    DROPS::SetupLBP1( mg, &A, lset.Phi);
    DROPS::MatrixCL L;
    L.LinComb( 1.0, A.Data, 1.0, M.Data);
    DROPS::VecDescCL b( &ifaceidx);
    DROPS::SetupInterfaceRhsP1( mg, &b, lset.Phi, rhs0);

    DROPS::WriteToFile( M.Data, "m_iface.txt", "M");
    DROPS::WriteToFile( A.Data, "a_iface.txt", "A");
    DROPS::WriteToFile( b.Data, "rhs_iface.txt", "rhs");

    typedef DROPS::SSORPcCL SurfPcT;
    SurfPcT surfpc;
    typedef DROPS::PCGSolverCL<SurfPcT> SurfSolverT;
    SurfSolverT surfsolver( surfpc, 1000, 1e-6, true);

    DROPS::VecDescCL x( &ifaceidx);
    surfsolver.Solve( L, x.Data, b.Data);
    std::cerr << "Iter: " << surfsolver.GetIter() << "\tres: " << surfsolver.GetResid() << '\n';

    DROPS::WriteToFile( x.Data, "x_iface.txt", "solution");

    DROPS::VecDescCL xext( &ifacefullidx);
    DROPS::Extend( mg, x, xext);
    DROPS::NoBndDataCL<> nobnd;
    DiscP1FunT surfsol( &xext, &nobnd, &mg, 0.);
    ensight.putScalar( datscl,  lset.GetSolution(), 0);
    ensight.putScalar( datsurf, surfsol, 0);
    ensight.Commit();

    double L2_err( L2_error( mg, lset.Phi, surfsol, &sol0));
    std::cerr << "L_2-error: " << L2_err << std::endl;

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
