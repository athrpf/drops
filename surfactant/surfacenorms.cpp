//****************************************************************************
// File:    brickflow.cpp                                                    *
// Content: test case: one-phase flow in square pipe with oscillating inflow *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen      *
//****************************************************************************
#include "surfactant/ifacetransp.h"
#include "surfactant/params.h"
#include "geom/builder.h"
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "out/ensightOut.h"

#include <fstream>

DROPS::ParamSurfactantCL C;
std::string filename,
            datgeo,
            datscl,
            datvel,
            datsurf,
            datsol,
            datlset;

DROPS::Point3DCL u_func (const DROPS::Point3DCL&, double)
{
    return C.Velocity;
}

double sphere_2 (const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL x( p - C.Mitte);

    return x.norm() - C.Radius[0];
}

double sphere_2move (const DROPS::Point3DCL& p, double t)
{
    DROPS::Point3DCL x( p - (C.Mitte + t*u_func(p, t)));

    return x.norm() - C.Radius[0];
}

typedef double (*dist_funT) (const DROPS::Point3DCL&, double);

void LSInit (const DROPS::MultiGridCL& mg, DROPS::VecDescCL& ls, dist_funT d, double t)
{
    const DROPS::Uint lvl= ls.GetLevel(),
                      idx= ls.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( it->GetCoord(), t);

    DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it)
        ls.Data[it->Unknowns( idx)]= d( 0.5*(it->GetVertex( 0)->GetCoord() + it->GetVertex( 1)->GetCoord()), t);
}

#define DROPS_FOR_TETRA_INTERFACE_BEGIN( t, ls, p, n) \
    (p).Init( (t), (ls)); \
    if ((p).Intersects()) { /*We are at the phase boundary.*/ \
        for (int ch__= 0; ch__ < 8; ++ch__) { \
            (p).ComputeForChild( ch__); \
            for (int n= 0; n < (p).GetNumTriangles(); ++n) \

#define DROPS_FOR_TETRA_INTERFACE_END }}
 

template <typename DiscP1FunT>
double L2_norm (const DROPS::MultiGridCL& mg, const DROPS::VecDescCL& ls,
    const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = ls.GetLevel();
    DROPS::InterfacePatchCL patch;
    DROPS::Quad5_2DCL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        DROPS_FOR_TETRA_INTERFACE_BEGIN( *it, ls, patch, tri) {
            qdiscsol.assign(  *it, &patch.GetBary( tri), discsol);
            d+= DROPS::Quad5_2DCL<>( qdiscsol*qdiscsol).quad( patch.GetFuncDet( tri));
        }
        DROPS_FOR_TETRA_INTERFACE_END
    }
    return std::sqrt( d);
}

template <typename DiscP1FunT>
double L2_norm_Omega (const DROPS::MultiGridCL& mg, const DiscP1FunT& discsol)
{
    double d( 0.);
    const DROPS::Uint lvl = discsol.GetLevel();
    DROPS::Quad5CL<> qdiscsol;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        qdiscsol.assign(  *it, discsol);
        d+= DROPS::Quad5CL<>( qdiscsol*qdiscsol).quad( it->GetVolume()*6.);
    }
    return std::sqrt( d);
}

namespace DROPS // for Strategy
{

const int N= 11;

void Strategy (DROPS::MultiGridCL& mg, DROPS::LevelsetP2CL& lset)
{
    using namespace DROPS;

    LSInit( mg, lset.Phi, &sphere_2move, 1.);
    IdxDescCL fullidx( P1_FE);
    //DROPS::CreateNumbOnInterface( mg.GetLastLevel(), idx, mg, lset.Phi, C.surf_omit_bound);
    fullidx.CreateNumbering( mg.GetLastLevel(), mg);
    NoBndDataCL<> bnd;

    VecDescCL DV[N];
    for (int i= 0; i < N; ++i)
        DV[i].SetIdx( &fullidx);

    ReadEnsightP2SolCL reader( mg);
    reader.ReadScalar( "ensight/surfactant.sur1",   DV[0], bnd);
    reader.ReadScalar( "ensight/surfactant.sur2",   DV[1], bnd);
    reader.ReadScalar( "ensight/surfactant.sur4",   DV[2], bnd);
    reader.ReadScalar( "ensight/surfactant.sur8",  DV[3], bnd);
    reader.ReadScalar( "ensight/surfactant.sur16",  DV[4], bnd);
    reader.ReadScalar( "ensight/surfactant.sur32",  DV[5], bnd);
    reader.ReadScalar( "ensight/surfactant.sur64", DV[6], bnd);
    reader.ReadScalar( "ensight/surfactant.sur128", DV[7], bnd);
    reader.ReadScalar( "ensight/surfactant.sur256",  DV[8], bnd);
    reader.ReadScalar( "ensight/surfactant.sur512", DV[9], bnd);
    reader.ReadScalar( "ensight/surfactant.sur1024", DV[10], bnd);

    for (int i= 0; i < N - 1; ++i) {
        DV[i].Data-=DV[N - 1].Data;
//        std::cerr << norm( DV[i].Data) << ", ";
        std::cerr << L2_norm( mg, lset.Phi, make_P1Eval( mg, bnd, DV[i], 1.)) << ", ";
    }
    std::cerr << std::endl;
    for (int i= 0; i < N - 1; ++i) {
        std::cerr << L2_norm_Omega( mg, make_P1Eval( mg, bnd, DV[i], 1.)) << ", ";
    }
    std::cerr << std::endl;
}

} // end of namespace DROPS


int main (int argc, char** argv)
{
  try
  {
    if (argc!=2)
    {
        std::cerr << "You have to specify one parameter:\n\t"
                  << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
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
    lset.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset.Phi.SetIdx( &lset.idx);
    Strategy( mg, lset);
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
