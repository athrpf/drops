#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"

using namespace DROPS;

typedef double (*fun_ptr)(const SVectorCL<3>&);

enum  OutputModeT { SILENT, NOISY };


double f(const SVectorCL<3>& p)
{ return p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2] +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]; }

double g(const SVectorCL<3>& p)
{  return p[0] +10.*p[1] +100.*p[2]+1000.; }

double h(const SVectorCL<3>& p)
{  return std::sin(M_PI*p[0])*std::sin(M_PI*p[1])*std::sin(M_PI*p[2]); }

double g2(const DROPS::Point3DCL& p)
{
    return (-1.)*p[0];
}

void MarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte( 0.5);

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It( mg.GetTriangTetraBegin(maxLevel)),
             ItEnd( mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It) {
        if ( (GetBaryCenter( *It)-Mitte).norm() <= std::max( 0.2, 1.5*std::pow( It->GetVolume(), 1.0/3.0)) )
            It->SetRemoveMark();
    }
}


typedef NoBndDataCL<double> BndCL;

BndCL Bnd;

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns);
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &Bnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord()));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5));
    }
}


int CheckResult( DROPS::TetraCL& t, fun_ptr f, double* dof, OutputModeT om)
{
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    if (om!=SILENT) std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0.;
    for (int i=0; i<4; ++i) {
        diff= dof[i] - f( t.GetVertex( i)->GetCoord());
        if ( std::abs( diff) > vmaxdiff) { ++ret; vmaxdiff= std::abs( diff); v= t.GetVertex( i); }
    }
    if (om!=SILENT) std::cout << "\n\nEdges:" << std::endl;
    for (int i= 0; i<6; ++i) {
        const DROPS::EdgeCL* sit= t.GetEdge( i);
        diff = dof[i+4] - f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5);
        if ( std::abs(diff) > emaxdiff) { ++ret; emaxdiff= std::abs(diff); e= sit; }
    }
    if (om!=SILENT) {
        std::cout << "maximale Differenz Vertices: " << vmaxdiff << " auf\n";
        if (v) v->DebugInfo( std::cout);
        std::cout << "maximale Differenz Edges: " << emaxdiff << " auf\n";
        if (e) e->DebugInfo( std::cout);
        std::cout << std::endl;
    }
    return ret;
}


DROPS::Uint Rule(DROPS::Uint r)
{
    return r < 64 ? r : 127;
}

// Checks every possible tetra-modification.
int Test()
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting restriction on single tetra-combinations:\n";
    int ttt, ret= 0;
    for (DROPS::Uint i= 0; i<=64; ++i) {
        for (DROPS::Uint j= 0; j<=64; ++j) {
            DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
//            std::cout << Rule( i) << "\t-->\t" << Rule( j) << " ";
            DROPS::TetraBuilderCL tet( Rule( i), DROPS::std_basis<3>( 1),
                                                 DROPS::std_basis<3>( 2),
                                                 DROPS::std_basis<3>( 3),
                                                 DROPS::Point3DCL( 1.0));
            DROPS::MultiGridCL mg( tet);
            DROPS::TetraCL& t0= *mg.GetTriangTetraBegin( 0);
            
            DROPS::IdxDescCL i0, i1;
            i0.Set( 1,1,0,0); i0.TriangLevel= mg.GetLastLevel(); i0.NumUnknowns= 0;
            DROPS::CreateNumbOnVertex( i0.GetIdx(), i0.NumUnknowns, 1,
                                       mg.GetTriangVertexBegin( i0.TriangLevel),
                                       mg.GetTriangVertexEnd( i0.TriangLevel),
                                       Bnd);
            DROPS::CreateNumbOnEdge( i0.GetIdx(), i0.NumUnknowns, 1,
                                     mg.GetTriangEdgeBegin( i0.TriangLevel),
                                     mg.GetTriangEdgeEnd( i0.TriangLevel),
                                     Bnd);
            DROPS::VecDescCL v0, v1;
            v0.SetIdx( &i0);
            SetFun( v0, mg, f);
//            SetFun( v0, mg, g2);

            i1.Set( 1,1,0,0);
            i1.TriangLevel= 0; i1.NumUnknowns= 0;
            DROPS::CreateNumbOnVertex( i1.GetIdx(), i1.NumUnknowns, 1,
                                       mg.GetTriangVertexBegin( i1.TriangLevel),
                                       mg.GetTriangVertexEnd( i1.TriangLevel),
                                       Bnd);
            DROPS::CreateNumbOnEdge( i1.GetIdx(), i1.NumUnknowns, 1,
                                     mg.GetTriangEdgeBegin( i1.TriangLevel),
                                     mg.GetTriangEdgeEnd( i1.TriangLevel),
                                     Bnd);
            v1.SetIdx( &i1);
            DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun0( &v0, &Bnd, &mg);
            double dof[10];
            RestrictP2( t0, v0, Bnd, dof);
            ttt= CheckResult( t0, f, dof, NOISY);
            RestrictP2( t0, fun0, dof);
            ttt= CheckResult( t0, f, dof, NOISY);
//            ttt= CheckResult( t0, g2, dof, NOISY);
            ret+= ttt;
        }
    }
    return ret;
}




int main ()
{
  try {
    return  Test();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
