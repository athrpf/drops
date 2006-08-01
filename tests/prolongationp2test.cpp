#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"

using namespace DROPS;

namespace DROPS {

class BndCL
{
  public:
    typedef double bnd_type;

    inline bool IsOnDirBnd (const VertexCL&) const { return false; }
    inline bool IsOnNeuBnd (const VertexCL&) const { return false; }
    inline bool IsOnDirBnd (const EdgeCL&) const { return false; }
    inline bool IsOnNeuBnd (const EdgeCL&) const { return false; }

    static inline bnd_type GetDirBndValue (const VertexCL&)
        { throw DROPSErrCL("BndCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
    static inline bnd_type GetDirBndValue (const EdgeCL&)
        { throw DROPSErrCL("BndCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on edge."); }
} Bnd;

} // end of namespace DROPS


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


int TestProlongation()
{
    int ret= 0;
    DROPS::IdxDescCL i0, i1;
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting prolongation for P2-elements:\n";
    for (DROPS::Uint i=0; i<63; ++i) {
        DROPS::TetraBuilderCL brick( i);
        DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
        DROPS::MultiGridCL mg( brick);
        i0.Set( 1,1,0,0); i0.TriangLevel= 0; i0.NumUnknowns= 0;
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
        i1.Set( 1,1,0,0); i1.TriangLevel= mg.GetLastLevel() ; i1.NumUnknowns= 0;
        DROPS::CreateNumbOnVertex( i1.GetIdx(), i1.NumUnknowns, 1,
                                   mg.GetTriangVertexBegin( i1.TriangLevel),
                                   mg.GetTriangVertexEnd( i1.TriangLevel),
                                   Bnd);
        DROPS::CreateNumbOnEdge( i1.GetIdx(), i1.NumUnknowns, 1,
                                 mg.GetTriangEdgeBegin( i1.TriangLevel),
                                 mg.GetTriangEdgeEnd( i1.TriangLevel),
                                 Bnd);
        v1.SetIdx( &i1);
        DROPS::MatDescCL P;
        SetupP2ProlongationMatrix( mg, P, &i0, &i1);
//        TetraCL* t= &*mg.GetTetrasBegin( 0);
        for (int k= 0; k<10; ++k) {
            VectorCL v( 10u);
            v[k]= 1.0;
            VectorCL w= P.Data*v;

//            for (TetraCL::ChildPIterator cp= t->GetChildBegin(); cp != t->GetChildEnd(); ++cp) {
//                (*cp)->DebugInfo( std::cout);
//                std::cout << std::endl;
//            }
            for (MultiGridCL::TriangVertexIteratorCL vi= mg.GetTriangVertexBegin(); vi != mg.GetTriangVertexEnd(); ++vi) {
//                vi->DebugInfo( std::cout);
//                std::cout << std::endl;
                SVectorCL<3> c= vi->GetCoord();
                if ( std::abs( w[vi->Unknowns( i1.GetIdx())] - FE_P2CL::H(k, c[0], c[1], c[2])) > 1e-15) {
                    std::cout << "Inconsistency: rule= " << i << ", k= " << k << std::endl;
                    ++ret;
                }
//                std::cout << "coarse: " << vi->Unknowns( i0.GetIdx()) << "\tfine: " << vi->Unknowns( i1.GetIdx())
//                          << std::endl << std::endl;
            }
            for (MultiGridCL::TriangEdgeIteratorCL ei= mg.GetTriangEdgeBegin(); ei != mg.GetTriangEdgeEnd(); ++ei) {
//                ei->DebugInfo( std::cout);
//                std::cout << std::endl;
//                std::cout << "coarse: " << ei->Unknowns( i0.GetIdx()) << "\tfine: " << ei->Unknowns( i1.GetIdx())
//                          << std::endl << std::endl;
                SVectorCL<3> c= ConvexComb( 0.5, ei->GetVertex( 0)->GetCoord(), ei->GetVertex( 1)->GetCoord());
                if ( std::abs( w[ei->Unknowns( i1.GetIdx())] - FE_P2CL::H(k, c[0], c[1], c[2])) > 1e-15) {
                    std::cout << "Inconsistency: rule= " << i << ", k= " << k << std::endl;
                    ++ret;
                }
            }
        }
//        std::cout << "Rule: " << i << std::endl;
//        std::cout << P.Data << std::endl;
        DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllVertexBegin( i0.TriangLevel),
                                    mg.GetAllVertexEnd( i0.TriangLevel));
        DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllEdgeBegin( i0.TriangLevel),
                                    mg.GetAllEdgeEnd( i0.TriangLevel));
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllVertexBegin( i1.TriangLevel),
                                    mg.GetAllVertexEnd( i1.TriangLevel));
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllEdgeBegin( i1.TriangLevel),
                                    mg.GetAllEdgeEnd( i1.TriangLevel));
    }

    std::cout << "\n-----------------------------------------------------------------"
              << std::endl;
    return ret;
}


// Returns 0, iff everything seems ok.
int main ()
{
  try {
    int ret= TestProlongation();
    return ret;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
