#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"

using namespace DROPS;
double f(const SVectorCL<3>& p)
{ return p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2] +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]; }

void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

class BndCL
{
  public:
    typedef BndCL self;
    typedef double bnd_type;

    static inline bool IsOnDirBnd (const VertexCL&) { return false; }
    static inline bool IsOnNeuBnd (const VertexCL&) { return false; }
    static inline bool IsOnDirBnd (const EdgeCL&) { return false; }
    static inline bool IsOnNeuBnd (const EdgeCL&) { return false; }
    
    static inline bnd_type GetDirBndValue (const VertexCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
    static inline bnd_type GetDirBndValue (const EdgeCL&)
        { throw DROPSErrCL("StokesBndDataPrCL::GetDirBndValue: Attempt to use Dirichlet-boundary-conditions on vertex."); }
};

template<bool _OnBnd>
void CreateNumbOnVertex(const Uint idx, IdxT& counter, Uint NumUnknown,
                        const MultiGridCL::TriangVertexIteratorCL& begin,
                        const MultiGridCL::TriangVertexIteratorCL& end)
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence. Simplices on Dirichlet boundaries are skipped.
{
    if (NumUnknown == 0) return;
    for ( MultiGridCL::TriangVertexIteratorCL it=begin; it!=end; ++it)
    {
      if ( _OnBnd || !BndCL::IsOnDirBnd(*it) )
      {        
// TODO: Fix this!    immer Alloc?  oder je nachdem Resize?
          if ( !it->Unknowns.Exist() ) it->Unknowns.Init(2);
          if ( !it->Unknowns.Exist(idx) )
              it->Unknowns.Get()->Alloc(idx, NumUnknown);
          for (Uint i=0; i<NumUnknown; ++i) it->Unknowns(idx)[i] = counter++;
      }
    }
}

template<bool _OnBnd>
void CreateNumbOnEdge(const Uint idx, IdxT& counter, Uint NumUnknown,
                      const MultiGridCL::TriangEdgeIteratorCL& begin,
                      const MultiGridCL::TriangEdgeIteratorCL& end)
// allocates memory for the Unknown-indices on all simplices between begin and end 
// and numbers them in order of appearence. Simplices on Dirichlet boundaries are skipped.
{
    if (NumUnknown == 0) return;
    for ( MultiGridCL::TriangEdgeIteratorCL it=begin; it!=end; ++it)
    {
      if ( _OnBnd || !(BndCL::IsOnDirBnd(*it)) )
      {        
// TODO: Fix this!    immer Alloc?  oder je nachdem Resize?
          if ( !it->Unknowns.Exist() ) it->Unknowns.Init(2);
          if ( !it->Unknowns(idx) )
              it->Unknowns.Get()->Alloc(idx, NumUnknown);
          for (Uint i=0; i<NumUnknown; ++i) it->Unknowns(idx)[i] = counter++;
      }
    }
}

void SetFun(VecDescBaseCL<VectorBaseCL<double> >& vd, MultiGridCL& mg)
{
    vd.Data.resize(vd.RowIdx->NumUnknowns);
    BndCL bnd;
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorBaseCL<double> > > fun(&vd, &bnd, &mg);
    const Uint lvl= vd.RowIdx->TriangLevel;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl), theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit)
    {
        fun.SetDoF(*sit, f(sit->GetCoord()));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl), theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit)
    {
        fun.SetDoF(*sit, f((sit->GetVertex(0)->GetCoord() + sit->GetVertex(1)->GetCoord())/2.));
    }
}

int main (int argc, char** argv)
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);
    DROPS::MultiGridCL mg(brick);
    MarkAll(mg);
    mg.Refine();

    IdxDescCL i0, i1;
    i0.Set(0, 1,0,0,0); i0.TriangLevel= 0; i0.NumUnknowns= 0;
    i1.Set(1, 1,0,0,0); i1.TriangLevel= 1; i1.NumUnknowns= 0;
    CreateNumbOnVertex<true>(i0.Idx, i0.NumUnknowns, 1, mg.GetTriangVertexBegin(i0.TriangLevel), mg.GetTriangVertexEnd(i0.TriangLevel) );
    CreateNumbOnVertex<true>(i1.Idx, i1.NumUnknowns, 1, mg.GetTriangVertexBegin(i1.TriangLevel), mg.GetTriangVertexEnd(i1.TriangLevel) );
    VecDescBaseCL<VectorBaseCL<double> > v0, v1;
    v0.SetIdx(&i0);
    v1.SetIdx(&i1);
    SetFun(v0, mg);
    v1.Data.resize(v1.RowIdx->NumUnknowns);
    BndCL bnd;
    P1EvalCL<double, BndCL, const VecDescBaseCL<VectorBaseCL<double> > > fun0(&v0, &bnd, &mg);
    P1EvalCL<double, BndCL,VecDescBaseCL<VectorBaseCL<double> > > fun1(&v1, &bnd, &mg);
    Interpolate(fun1, fun0);
    std::cout << "Verts:" << std::endl;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(1), theend= mg.GetTriangVertexEnd(1); sit!=theend; ++sit)
    {
        std::cout << fun1.val(*sit) - f(sit->GetCoord()) << std::endl;
    }
    std::cout << "Edges:" << std::endl;
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(1), theend= mg.GetTriangEdgeEnd(1); sit!=theend; ++sit)
    {
        std::cout << fun1.val(*sit, .5) - f((sit->GetVertex(0)->GetCoord()+sit->GetVertex(1)->GetCoord())/2.) << std::endl;
    }

    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
