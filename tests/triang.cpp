#include "misc/utils.h"
#include "geom/multigrid.h"
#include "geom/builder.h"

using namespace DROPS;

enum  OutputModeT { SILENT, NOISY };

template <class SimplexT>
struct TriangFillCL;

template <class SimplexT>
class TriangCL
{
  public:
    typedef std::vector<SimplexT*> LevelCont;
    typedef SimplexT**             iterator;
    typedef const SimplexT**       const_iterator;

  private:
    mutable std::vector<LevelCont> triang_;
    MultiGridCL&                   mg_;

    int  StdIndex    (int lvl) const { return lvl >= 0 ? lvl : lvl + triang_.size(); }
    inline void MaybeCreate (int lvl) const;

  public:
    TriangCL (MultiGridCL& mg) : triang_( mg.GetNumLevel()), mg_( mg) {}

    void clear () { triang_.clear(); }
    Uint size  (int lvl= -1) const
        { MaybeCreate( StdIndex( lvl)); return triang_[StdIndex( lvl)].size(); }

    iterator begin (int lvl= -1)
        { MaybeCreate( StdIndex( lvl)); return &*triang_[StdIndex( lvl)].begin(); }
    iterator end   (int lvl= -1)
        { MaybeCreate( StdIndex( lvl)); return &*triang_[StdIndex( lvl)].end(); }
    const_iterator begin (int lvl= -1) const
        { MaybeCreate( StdIndex( lvl)); return const_cast<const_iterator>( &*triang_[StdIndex( lvl)].begin()); }
    const_iterator end   (int lvl= -1) const
        { MaybeCreate( StdIndex( lvl)); return const_cast<const_iterator>( &*triang_[StdIndex( lvl)].end()); }
};

template <class SimplexT>
  inline void
  TriangCL<SimplexT>::MaybeCreate(int lvl) const
{
    Assert ( StdIndex( lvl) >= 0 && StdIndex( lvl) < mg_.GetNumLevel(),
        DROPSErrCL( "TriangCL::MaybeCreate: Wrong level."), DebugContainerC);
    if (triang_.size() != mg_.GetNumLevel()) {
        triang_.clear();
        triang_.resize( mg_.GetNumLevel());
    }
    const int level= StdIndex( lvl);
    if (triang_[level].empty())
        TriangFillCL<SimplexT>::fill( mg_, triang_[level], level);
}

template <class SimplexT>
struct TriangFillCL
{
  static void
  fill (MultiGridCL& mg, typename TriangCL<SimplexT>::LevelCont& c, int lvl);
};

template <>
struct TriangFillCL<VertexCL>
{
    static void
    fill (MultiGridCL& mg, TriangCL<VertexCL>::LevelCont& c, int lvl) {
        for (MultiGridCL::VertexIterator it= mg.GetAllVertexBegin( lvl),
             theend= mg.GetAllVertexEnd( lvl); it != theend; ++it)
            if (it->IsInTriang( lvl)) c.push_back( &*it);
        TriangCL<VertexCL>::LevelCont tmp= c;
        c.swap( tmp);
    }
};
template <>
struct TriangFillCL<EdgeCL>
{
    static void
    fill (MultiGridCL& mg, TriangCL<EdgeCL>::LevelCont& c, int lvl) {
        for (MultiGridCL::EdgeIterator it= mg.GetAllEdgeBegin( lvl),
             theend= mg.GetAllEdgeEnd( lvl); it != theend; ++it)
            if (it->IsInTriang( lvl)) c.push_back( &*it);
        TriangCL<EdgeCL>::LevelCont tmp= c;
        c.swap( tmp);
    }
};
template <>
struct TriangFillCL<FaceCL>
{
    static void
    fill (MultiGridCL& mg, TriangCL<FaceCL>::LevelCont& c, int lvl) {
        for (MultiGridCL::FaceIterator it= mg.GetAllFaceBegin( lvl),
             theend= mg.GetAllFaceEnd( lvl); it != theend; ++it)
            if (it->IsInTriang( lvl)) c.push_back( &*it);
        TriangCL<FaceCL>::LevelCont tmp= c;
        c.swap( tmp);
    }
};
template <>
struct TriangFillCL<TetraCL>
{
    static void
    fill (MultiGridCL& mg, TriangCL<TetraCL>::LevelCont& c, int lvl) {
        for (MultiGridCL::TetraIterator it= mg.GetAllTetraBegin( lvl),
             theend= mg.GetAllTetraEnd( lvl); it != theend; ++it)
            if (it->IsInTriang( lvl)) c.push_back( &*it);
        TriangCL<TetraCL>::LevelCont tmp= c;
        c.swap( tmp);
    }
};


typedef  TriangCL<VertexCL> VertexTriangCL;
typedef  TriangCL<EdgeCL>   EdgeTriangCL;
typedef  TriangCL<FaceCL>   FaceTriangCL;
typedef  TriangCL<TetraCL>  TetraTriangCL;


void MarkDrop(DROPS::MultiGridCL& mg, int maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.2,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
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


Uint
Old( DROPS::MultiGridCL& mg)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nOld Iterators:\n";

    Uint tmp= 0;
    for (MultiGridCL::TriangVertexIteratorCL sit= mg.GetTriangVertexBegin(),
         end=mg.GetTriangVertexEnd(); sit != end; ++sit) {
        tmp+= sit->GetId().GetIdent();
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit= mg.GetTriangEdgeBegin(),
         end=mg.GetTriangEdgeEnd(); sit != end; ++sit) {
        ++tmp;
    }
    for (MultiGridCL::TriangTetraIteratorCL sit= mg.GetTriangTetraBegin(),
         end=mg.GetTriangTetraEnd(); sit != end; ++sit) {
        tmp+= sit->GetId().GetIdent();
    }
    return tmp;
}

Uint
New( DROPS::MultiGridCL&, const VertexTriangCL& vt, const EdgeTriangCL& et,
    const TetraTriangCL& tt)
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nNew Iterators:\n";

    Uint tmp= 0;
    for (VertexTriangCL::const_iterator sit= vt.begin(),
         end=vt.end(); sit != end; ++sit) {
        tmp+= (*sit)->GetId().GetIdent();
    }
    for (EdgeTriangCL::const_iterator sit= et.begin(),
         end=et.end(); sit != end; ++sit) {
        ++tmp;
    }
    for (TetraTriangCL::const_iterator sit= tt.begin(),
         end=tt.end(); sit != end; ++sit) {
        tmp+= (*sit)->GetId().GetIdent();
    }
    return tmp;
}

int main ()
{
  try {
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
				DROPS::std_basis<3>(2),
				DROPS::std_basis<3>(3),
				30, 30, 30);
    DROPS::MultiGridCL mg( brick);
    mg.SizeInfo( std::cout);
    MarkDrop( mg, -1);
    mg.Refine();
    MarkDrop( mg, -1);
    mg.Refine();
    MarkDrop( mg, -1);
    mg.Refine();
    mg.SizeInfo( std::cout);
   
    VertexTriangCL vt( mg);
    EdgeTriangCL et( mg);
    TetraTriangCL tt( mg);

    TimerCL time;
    time.Start();
    Uint q0= Old( mg);
    time.Stop();
    std::cout << "value: " << q0 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    time.Start();
    Uint q1= New( mg, vt, et, tt);
    time.Stop();
    std::cout << "value: " << q1 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    time.Start();
    Uint q2= New( mg, vt, et, tt);
    time.Stop();
    std::cout << "value: " << q2 << "\ttime: " << time.GetTime() << " seconds"
        << std::endl;
    time.Reset();
    std::cout << "tt.size: " << tt.size() << std::endl;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
