#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"

using namespace DROPS;
double f(const SVectorCL<3>& p)
//{ return 42; }
{ return p[0] +10.*p[1] +100.*p[2]+1; }

void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

class BndCL
{
  public:
    typedef BndCL self;
    typedef double bnd_type;

    inline bool IsOnDirBnd (const VertexCL&) const { return false; }
    inline bool IsOnNeuBnd (const VertexCL&) const { return false; }
    inline bool IsOnDirBnd (const EdgeCL&) const { return false; }
    inline bool IsOnNeuBnd (const EdgeCL&) const { return false; }
    
    static inline bnd_type GetDirBndValue (const VertexCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
    static inline bnd_type GetDirBndValue (const EdgeCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
} Bnd;

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg)
{
    vd.Data.resize(vd.RowIdx->NumUnknowns);
    P1EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun(&vd, &Bnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl), theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit)
    {
        fun.SetDoF(*sit, f(sit->GetCoord()));
    }
}



int main (int argc, char** argv)
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;
//    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);
    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 3, 3, 3);
    DROPS::MultiGridCL mg(brick);
    MarkDrop(mg, 0);
//    MarkAll( mg);
    mg.Refine();

    IdxDescCL i0, i1;
    i0.Set(0, 1,0,0,0); i0.TriangLevel= 0; i0.NumUnknowns= 0;
    i1.Set(1, 1,0,0,0); i1.TriangLevel= 1; i1.NumUnknowns= 0;
    CreateNumbOnVertex(i0.Idx, i0.NumUnknowns, 1, mg.GetTriangVertexBegin(i0.TriangLevel), mg.GetTriangVertexEnd(i0.TriangLevel), Bnd );
    CreateNumbOnVertex(i1.Idx, i1.NumUnknowns, 1, mg.GetTriangVertexBegin(i1.TriangLevel), mg.GetTriangVertexEnd(i1.TriangLevel), Bnd );

    VecDescBaseCL<VectorCL> v0, v1;
    v0.SetIdx(&i0);
    v1.SetIdx(&i1);
    SetFun(v0, mg);
    
    P1EvalCL<double, BndCL, const VecDescBaseCL<VectorCL> > fun0(&v0, &Bnd, &mg);
    P1EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun1(&v1, &Bnd, &mg);
    Interpolate(fun1, fun0);
    std::cout << "Verts:" << std::endl;
    double diff;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(1), theend= mg.GetTriangVertexEnd(1); sit!=theend; ++sit)
    {
        diff= fun1.val(*sit) - f(sit->GetCoord()); 
        std::cout << diff << "\t";
        if (diff!=0.) return 1;
    }

    std::cerr << std::endl << DROPS::SanityMGOutCL(mg) << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
