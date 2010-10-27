/// \file fastmarch.cpp
/// \brief reparametrization methods
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

//#define MINPACK

#include "levelset/fastmarch.h"
#include "num/solver.h"
#include <fstream>
#include <cstring>
#ifdef _OPENMP
#  include <omp.h>
#endif
#if __GNUC__ >= 4 && !defined(__INTEL_COMPILER)
#    include <tr1/unordered_map>
#else
#    include <map>
#endif
#include <set>
#ifdef MINPACK
#  include <cminpack.h>
#endif

namespace DROPS
{
#define COUNTMEM

// R E P A R A M  D A T A  C L
//----------------------------

ReparamDataCL::~ReparamDataCL()
{
    for (size_t i=0; i<perpFoot.size(); ++i)
        if ( perpFoot[i]!=0)
            delete perpFoot[i];
}

/** Initializes the data structures augmIdx_ and map_.*/
void ReparamDataCL::InitPerMap()
{
    if (!per)
        return;

    const int    lvl     = phi.GetLevel();
    augmIdx= new IdxDescCL( P2_FE, *bnd);
    augmIdx->CreateNumbering( lvl, mg);
    const Uint   augm_idx= augmIdx->GetIdx();
    const Uint   idx     = phi.RowIdx->GetIdx();
    const size_t size    = phi.Data.size();      // number of DOF for lset
    const size_t sizeAugm= augmIdx->NumUnknowns();

    std::vector<bool> ini( size, false);
    IdxT k= 0;
    map.resize( sizeAugm-size);

    DROPS_FOR_TRIANG_VERTEX( mg, lvl, it){
        const IdxT Nr= it->Unknowns( idx);
        if (ini[Nr]){
            it->Unknowns( augm_idx)= size + k;
            map[k++]= Nr;
        }
        else{ // touched for the first time
            ini[Nr]= true;
            it->Unknowns( augm_idx)= Nr;
        }
    }
    DROPS_FOR_TRIANG_EDGE( mg, lvl, it){
        const IdxT Nr= it->Unknowns( idx);
        if (ini[Nr]){
            it->Unknowns( augm_idx)= size + k;
            map[k++]= Nr;
        }
        else{ // touched for the first time
            ini[Nr]= true;
            it->Unknowns( augm_idx)= Nr;
        }
    }
}


/** Iterate over all vertices and edges of the finest triangulation and put their
    coordinates (or barycenter) in coord_.
    \pre for periodic boundaries InitPerMap has to be called first
*/
void ReparamDataCL::InitCoord()
{
    const size_t coord_size= per ? augmIdx->NumUnknowns() : phi.Data.size();
    const Uint idx         = per ? augmIdx->GetIdx()      : phi.RowIdx->GetIdx();
    const int lvl= phi.GetLevel();
    coord.resize( coord_size);


#pragma omp parallel for
    for ( int i=0; i<std::distance(mg.GetTriangVertexBegin(lvl), mg.GetTriangVertexEnd(lvl)); ++i ){
        MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl)+i;
        coord[ it->Unknowns( idx)]= it->GetCoord();
    }
#pragma omp parallel for
    for ( int i=0; i<std::distance(mg.GetTriangEdgeBegin(lvl), mg.GetTriangEdgeEnd(lvl)); ++i ){
        MultiGridCL::TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(lvl)+i;
        coord[ it->Unknowns( idx)]= GetBaryCenter( *it);
    }
}

// I N I T  Z E R O  C L
//----------------------

/** Free all memory, especially free memory of created augmented index*/
InitZeroCL::~InitZeroCL()
{
    if ( data_.augmIdx){
        data_.augmIdx->DeleteNumbering( data_.mg);
        delete data_.augmIdx;
        data_.augmIdx= 0;
    }
}

// I N I T  Z E R O  N O  M O D  C L
//----------------------------------

/** Iterate over all tetrahedra and check, if a child is intersected by the zero level
    of the level set function
*/
void InitZeroNoModCL::Perform()
{
    InterfaceTriangleCL patch;                              // check for intersection
    LocalNumbP2CL       n;                                  // local numbering of dof
    const RefRuleCL    RegRef= GetRefRule( RegRefRuleC);    // determine regular children
    const int lvl= base::data_.phi.GetLevel();

#pragma omp parallel for
    for ( int i=0; i<std::distance(data_.mg.GetTriangTetraBegin(lvl), data_.mg.GetTriangTetraEnd(lvl)); ++i ){
        MultiGridCL::TriangTetraIteratorCL it= data_.mg.GetTriangTetraBegin(lvl)+i;
        patch.Init( *it, base::data_.phi, *base::data_.bnd);
        if ( patch.Intersects()){                                       // tetra (*it) is intersected
            n.assign( *it, *base::data_.phi.RowIdx, BndDataCL<>(0));    // create local numbering
            for (Uint ch= 0; ch<MaxChildrenC; ++ch) {
                if (patch.ComputeForChild( ch)){                        // Child ch has an intersection
                    const ChildDataCL data= GetChildData( RegRef.Children[ch]);
                    // mark all vertices as finished
                    for ( int vert=0; vert<4; ++vert){
                        IdxT dof= data_.Map( n.num[ data.Vertices[vert]]);
#pragma omp critical
{
                        base::data_.typ[ dof]= data_.Finished;
                        base::data_.phi.Data[ dof]= std::abs(base::data_.phi.Data[ dof]);
}
                    }
                }
            }
        }
    }
}

// I N I T  Z E R O  E X A C T  C L
//---------------------------------

/// \brief The class provides, for each tetra, that is cut by the interface, a set of all levelset-unknowns, that are at most two tetras away. Only levelset-unknowns are considered, that belong to a tetra that is cut by the interface.
///
/// Consider the graph with level-set-unknowns (and their respective periodic unknowns) at the interface as vertices. The edges of the graph are tetras.
/// The class computes, for each vertex, the set of vertices that are at most two edges away.
/// The algorithm is breadth-first search. This mapping is inverted to provide the desired map.
/// The work is done in the constructor. The resulting map can be accessed via operator().
///
/// The parallelization of this class probably involves a lot of work.
class TetraNeighborCL
{
  public:
    typedef std::set<IdxT> IdxSetT;
#if __GNUC__ >= 4 && !defined(__INTEL_COMPILER)
    typedef std::tr1::unordered_map<const TetraCL*, IdxSetT> TetraToIdxMapT;
#else
    typedef std::map<const TetraCL*, IdxSetT> TetraToIdxMapT;
#endif
    typedef std::vector<std::pair<const TetraCL*, IdxSetT> > TetraToIdxVecT;

  private:
    typedef std::set<const TetraCL*> TetraSetT;
    typedef std::vector<TetraSetT>   IdxToTetraMapT;

    IdxToTetraMapT idx_to_tetra_; ///< Used to accumulate the reverse of the desired mapping; cleared at the end of the constructor.
    TetraToIdxMapT tetra_to_idx_; ///< A set of levelset-unknowns for each tetra, that is cut by the interface.

    ///
    /// \brief Updates the neighborhood information: traverse all outgoing edges (=tetras) of idx and add them to their end-vertices.
    void insert_neighbor_tetras (size_t idx, const TetraCL& t, Uint ls_idx, const IdxToTetraMapT& idx_to_tetra);
    /// \brief Traverse the edges of the graph and update the vertex-neighborhoods.
    void ComputeNeighborTetras (const MultiGridCL& mg, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const ReparamDataCL::perMapVecT& mapPerDoF);
    /// \brief Reverse the map produced by ComputeNeighborTetras().
    void ComputeIdxToTetraMap ();

  public:
    TetraNeighborCL (const MultiGridCL& mg, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const ReparamDataCL::perMapVecT& mapPerDoF)
    {
        ComputeNeighborTetras ( mg, ls, lsetbnd, mapPerDoF);
        ComputeIdxToTetraMap ();
        idx_to_tetra_.clear();
    }

    /// \brief Accessor to the map from levelset-unknowns to tetras.
    TetraToIdxMapT& operator() () { return tetra_to_idx_; }

    /// \brief Accessor to the map tetra_to_idx_ as a vector (used for OpenMP purpose)
    TetraToIdxVecT GetTetraToIdxAsVector() { return Map2Vec(tetra_to_idx_); }
};

void TetraNeighborCL::insert_neighbor_tetras (size_t idx, const TetraCL& t, Uint ls_idx,
    const IdxToTetraMapT& idx_to_tetra)
{
    IdxT idx_j;
    for (Uint j= 0; j < 4; ++j) {
        idx_j= t.GetVertex( j)->Unknowns( ls_idx);
        idx_to_tetra_[idx].insert( idx_to_tetra[idx_j].begin(), idx_to_tetra[idx_j].end());
    }
    for (Uint j= 0; j < 6; ++j) {
        idx_j= t.GetEdge( j)->Unknowns( ls_idx);
        idx_to_tetra_[idx].insert( idx_to_tetra[idx_j].begin(), idx_to_tetra[idx_j].end());
    }
}

void TetraNeighborCL::ComputeNeighborTetras (const MultiGridCL& mg, const VecDescCL& ls, const BndDataCL<>& lsetbnd, const ReparamDataCL::perMapVecT& mapPerDoF)
{
    const DROPS::Uint lvl= ls.GetLevel();
    const Uint ls_idx= ls.RowIdx->GetIdx();
    const IdxT ls_size= ls.Data.size();

    DROPS::InterfacePatchCL patch;
    IdxToTetraMapT idx_to_tetra1( ls_size);
    idx_to_tetra_.resize( ls_size + mapPerDoF.size());

    // Record for each index, which tetras are immediate neighbors. Do this only around the interface.
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls, lsetbnd);
        if (!patch.Intersects()) continue;

        for (Uint i= 0; i < 4; ++i)
            idx_to_tetra1[it->GetVertex( i)->Unknowns( ls_idx)].insert( &*it);
        for (Uint i= 0; i < 6; ++i)
            idx_to_tetra1[it->GetEdge( i)->Unknowns( ls_idx)].insert( &*it);
    }

    // Record now also the neighbors of neighbors. Do this only, if the idx belongs to a tetra that is cut by the interface.
    IdxT idx;
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        for (Uint i= 0; i < 4; ++i) {
            idx= it->GetVertex( i)->Unknowns( ls_idx);
            if (idx_to_tetra1[idx].empty()) continue; // The vertex has no tetra on an interface.
            insert_neighbor_tetras( idx, *it, ls_idx, idx_to_tetra1);
        }
        for (Uint i= 0; i < 6; ++i) {
            idx= it->GetEdge( i)->Unknowns( ls_idx);
            if (idx_to_tetra1[idx].empty()) continue; // The edge has no tetra on an interface.
            insert_neighbor_tetras( idx, *it, ls_idx, idx_to_tetra1);
        }
    }

    // for dof's on periodic boundaries the respective entries of idx_to_tetra_ are copied
    for (size_t i=0, s= mapPerDoF.size(); i<s; ++i) {
        const IdxT augm_dof= ls_size + i,
                        dof= mapPerDoF[i];
        idx_to_tetra_[augm_dof]= idx_to_tetra_[dof];
    }
}

void TetraNeighborCL::ComputeIdxToTetraMap ()
{
    for (size_t idx= 0; idx < idx_to_tetra_.size(); ++idx)
        for (TetraSetT::iterator sit= idx_to_tetra_[idx].begin(); sit != idx_to_tetra_[idx].end(); ++sit)
            tetra_to_idx_[*sit].insert( idx);
}

/// \brief Computes the distance d from a point p to the triangle tri_.
///
/// The latter is given via the world coordinates of its vertices in the constructor.
///
/// The base point l is the point in tri with distance d from p. Although currently not used,
/// the code to compute l is tested and functional. l is given in barycentric coordinates with respect to tri.
class ExactDistanceInitCL
{
  private:
    /// \brief Computes the distance d of a point p to the line through v0, v1.  If the base l is not in the segment [v0, v1], returns false.
    static bool dist (const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& p, double& d /*, Point2DCL& l*/);

    const Point3DCL* tri_; ///< The current triangle to which distances are computed; *no* copies of the points.
    QRDecompCL<3,2> qr_;

  public:
    ExactDistanceInitCL (const Point3DCL tri[3]);

    /// \brief Computes the distance d from a point p to the triangle tri_. The latter is given via the world coordinates of its vertices.
    void dist (/*const Point3DCL tri[3],*/ const Point3DCL& p, double& d /*, SVectorCL<3>& l*/);

    /// \brief Simple tests to check the correctness of local distance computations manually.
    static void Test ();
};

ExactDistanceInitCL::ExactDistanceInitCL (const Point3DCL tri[3])
    : tri_( tri)
{
    SMatrixCL<3,2> M( Uninitialized);
    M.col( 0, tri[1] - tri[0]);
    M.col( 1, tri[2] - tri[0]);
    qr_.GetMatrix()= M;
    qr_.prepare_solve();
}

void ExactDistanceInitCL::Test ()
{
    Point3DCL tri[3];
    tri[0]= std_basis<3>( 0);
    tri[1]= std_basis<3>( 1);
//    tri[2]= std_basis<3>( 2);
    tri[2]= tri[1]+1e-5*std_basis<3>( 2);
    ExactDistanceInitCL dist_to_tri( tri);

    Point3DCL p= std_basis<3>( 3);
    p[0]= .5;
    p[1]= -1e5;
//    p= -p;
    double d;
    Point3DCL l;
    dist_to_tri.dist( /*tri,*/ p, d /*, l*/);
    std::cerr.precision( 20);
    std::cerr << "p: " << p << " d: " << d << std::endl; // << " l: " << l << std::endl;
}

bool ExactDistanceInitCL::dist (const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& p, double& d /*, Point2DCL& l*/)
{
    const double ll= inner_prod( v1 - v0, p - v0)/(v1 - v0).norm_sq();
    if (ll >= 0. && 1. - ll >= 0.) {
        d= ((1. - ll)*v0 + ll*v1 - p).norm(); // l= MakePoint2D( 1. - ll, ll);
        return true;
    }
    return false;
}

void ExactDistanceInitCL::dist (/*const Point3DCL tri[3],*/ const Point3DCL& p, double& d /*, SVectorCL<3>& l*/)
{
//     SMatrixCL<3,2> M;
//     M.col( 0, tri[1] - tri[0]);
//     M.col( 1, tri[2] - tri[0]);
//     QRDecompCL<3, 2> qr( M);*/
    Point3DCL ll( p - tri_[0]);
    qr_.Solve( ll); // The first two components are the least-squares-solution, the last is the residual.

    if (ll[0] >= 0. && ll[1] >= 0. && 1. - ll[0] - ll[1] >= 0.) {
        d= std::fabs( ll[2]); // l= MakePoint3D( 1. - ll[0] - ll[1], ll[0], ll[1]);
        return;
    }

    // Compute the distance to the three boundary segments.
    double tmpd;
    // Point2DCL ll2;
    if (ll[1] < 0. && ExactDistanceInitCL::dist( tri_[0], tri_[1], p, tmpd /*, ll2*/)) {
        d= tmpd; // l= MakePoint3D( ll2[0], ll2[1], 0.);
        return;
    }
    if (ll[0] < 0. && ExactDistanceInitCL::dist( tri_[0], tri_[2], p, tmpd /*, ll2*/)) {
        d= tmpd; // l= MakePoint3D( ll2[0], 0., ll2[1]);
        return;
    }
    if (1. - ll[0] - ll[1] < 0. && ExactDistanceInitCL::dist( tri_[1], tri_[2], p, tmpd /*, ll2*/)) {
        d= tmpd; // l= MakePoint3D( 0., ll2[0], ll2[1]);
        return;
    }

    // Compute the distance to the vertices.
    d= std::min( (p - tri_[0]).norm(), std::min( (p - tri_[1]).norm(), (p - tri_[2]).norm()));
}

// I N I T  Z E R O  E X A C T  C L
//---------------------------------

/** Copy corners of the triangle and compute QR decomposition*/
InitZeroExactCL::DistanceTriangCL::DistanceTriangCL(const Point3DCL tri[3])
{
    std::copy( tri, tri+3, tri_);
    SMatrixCL<3,2> M( Uninitialized);
    M.col( 0, tri[1] - tri[0]);
    M.col( 1, tri[2] - tri[0]);
    qr_.GetMatrix()= M;
    qr_.prepare_solve();
}

/** Take data from the buffer to initialize this class
    The buffer must contain 9 doubles for tri_ and 10
    doubles for qr_
*/
InitZeroExactCL::DistanceTriangCL::DistanceTriangCL(const double* buff)
{
    for ( int t=0; t<3; ++t) for ( int d=0; d<3; ++d)
        tri_[t][d]= buff[t*3+d];
    qr_.Deserialize( buff+9);
}


bool InitZeroExactCL::DistanceTriangCL::dist (const Point3DCL& v0, const Point3DCL& v1, const Point3DCL& p, double& d, Point2DCL*)
{
    const double ll= inner_prod( v1 - v0, p - v0)/(v1 - v0).norm_sq();
    if (ll >= 0. && 1. - ll >= 0.) {
        d= ((1. - ll)*v0 + ll*v1 - p).norm(); // l= MakePoint2D( 1. - ll, ll);
        return true;
    }
    return false;
}


/** Compute the distance of point p to the triangle. If a pointer l is given,
    the corresponding perpednicular foot is given here on return.
    \param p       point, to which the distance is computed
    \param locPerp location of the perpendicular foot. 2 : inside triangle,
                   1 : on boundary segment, 0 : corner of triangle
    \param l       perpendicular foot
*/
double InitZeroExactCL::DistanceTriangCL::dist( const Point3DCL& p, byte* locPerp, Point3DCL* l)
{
    Point3DCL ll( p - tri_[0]);
    qr_.Solve( ll); // The first two components are the least-squares-solution, the last is the residual.

    // ll is located in triangle
    if (ll[0] >= 0. && ll[1] >= 0. && 1. - ll[0] - ll[1] >= 0.) {
        if ( l!=0){
            *l= MakePoint3D( 1. - ll[0] - ll[1], ll[0], ll[1]);
            *locPerp= 2;
        }
        return std::fabs( ll[2]);
    }

    // Compute the distance to the three boundary segments.
    double tmpd;
    Point2DCL* ll2= (l==0) ? 0 : new Point2DCL;
    // Point2DCL ll2;
    if (ll[1] < 0. && dist( tri_[0], tri_[1], p, tmpd, ll2)) {
        if ( l!=0){
            *l= MakePoint3D( (*ll2)[0], (*ll2)[1], 0.);
            *locPerp= 1;
        }
        return tmpd;
    }
    if (ll[0] < 0. && dist( tri_[0], tri_[2], p, tmpd, ll2)) {
        if ( l!=0){
            *l= MakePoint3D( (*ll2)[0], 0., (*ll2)[1]);
            *locPerp= 1;
        }
        return  tmpd;
    }
    if (1. - ll[0] - ll[1] < 0. && dist( tri_[1], tri_[2], p, tmpd, ll2)) {
        if ( l!=0){
            *l= MakePoint3D( 0., (*ll2)[0], (*ll2)[1]);
            *locPerp= 1;
        }
        return tmpd;
    }

    // Compute the distance to the vertices.
    if ( l!=0) *locPerp= 0;
    return std::min( (p - tri_[0]).norm(), std::min( (p - tri_[1]).norm(), (p - tri_[2]).norm()));
}

/** Collect all neighbor and neighbor-neighbor tetras of a level set dof*/
void InitZeroExactCL::InitDofToTetra()
{
    const Uint lvl= data_.phi.GetLevel();
    const Uint idx= data_.phi.RowIdx->GetIdx();
    const IdxT size= data_.phi.Data.size();
    const MultiGridCL& mg= data_.mg;

    DROPS::InterfacePatchCL patch;
    DofToTetraMapT dof_to_tetra1( size);
    dofToTetra_.resize( size + data_.map.size());

    // Record for each index, which tetras are immediate neighbors. Do this only around the interface.
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, data_.phi, *data_.bnd);
        if (!patch.Intersects()) continue;

        for (Uint i= 0; i < 4; ++i)
            dof_to_tetra1[it->GetVertex( i)->Unknowns( idx)].insert( &*it);
        for (Uint i= 0; i < 6; ++i)
            dof_to_tetra1[it->GetEdge( i)->Unknowns( idx)].insert( &*it);
    }

    // Record now also the neighbors of neighbors. Do this only, if the idx belongs to a tetra that is cut by the interface.
    IdxT dof;
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        for (Uint i= 0; i<10; ++i) {
            dof= (i<4) ? it->GetVertex( i)->Unknowns( idx) : it->GetEdge( i-4)->Unknowns( idx);
            if (dof_to_tetra1[dof].empty()) continue; // The vertex/edge has no tetra on an interface.
            insert_neighbor_tetras( dof, *it, idx, dof_to_tetra1);
        }
    }

    // for dof's on periodic boundaries the respective entries of idx_to_tetra_ are copied
    for (size_t i=0, s= data_.map.size(); i<s; ++i) {
        const IdxT augm_dof= size + i,
                        dof= data_.map[i];
        dofToTetra_[augm_dof]= dofToTetra_[dof];
    }
}

/** Put all neighbors of dof into the intersected tetra set, too */
void InitZeroExactCL::insert_neighbor_tetras (IdxT dof, const TetraCL& t, Uint idx,
    const DofToTetraMapT& dof_to_tetra)
{
    for (Uint j=0; j<10; ++j) {
        const IdxT dofJ=(j<4) ? t.GetVertex( j)->Unknowns( idx) : t.GetEdge( j-4)->Unknowns( idx);
        dofToTetra_[dof].insert( dof_to_tetra[dofJ].begin(), dof_to_tetra[dofJ].end());
    }
}

inline bool operator < (const InitZeroExactCL::DistTriangIndexHelperCL& a, const InitZeroExactCL::DistTriangIndexHelperCL& b)
{
    if ( a.tetra<b.tetra)
        return true;
    else if ( a.tetra==b.tetra){
        if ( a.childNum<b.childNum)
            return true;
        else if ( a.childNum==b.childNum){
            if ( a.triangNum<b.triangNum)
                return true;
        }
    }
    return false;
}

inline bool operator == (const InitZeroExactCL::DistTriangIndexHelperCL& a, const InitZeroExactCL::DistTriangIndexHelperCL& b)
{
    return (a.tetra==b.tetra && a.childNum==b.childNum && a.triangNum==b.triangNum);
}

/** Determine all triangles representing the interface */
void InitZeroExactCL::BuildDistTriang()
{
    const MultiGridCL&  mg= data_.mg;
    const Uint          lvl= data_.phi.GetLevel();
    InterfaceTriangleCL patch;
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it){
        patch.Init( *it, data_.phi, *data_.bnd);
        if ( !patch.Intersects()) continue;
        for ( int ch=0; ch<8; ++ch){
            if (!patch.ComputeForChild( ch)) continue; // Child ch has no intersection
            for (int tri= 0; tri < patch.GetNumTriangles(); ++tri) {
                distTriangPos_[ DistTriangIndexHelperCL(&(*it), ch, tri)]= distTriang_.size();
                distTriang_.push_back( DistanceTriangCL(&patch.GetPoint( tri)));
            }
        }
    }
}

/** Assign each dof located in the vicinity of the interface the position
    of close triangles in distTriang_. */
void InitZeroExactCL::AssociateTriangles()
{
    const RefRuleCL     RegRef= GetRefRule( RegRefRuleC);  // determine regular children
    InterfaceTriangleCL patch;
    dofToDistTriang_.resize( data_.phi.Data.size());

    for ( size_t dof=0; dof<dofToTetra_.size(); ++dof){
        for ( TetraSetT::const_iterator it= dofToTetra_[dof].begin(); it!=dofToTetra_[dof].end(); ++it){
            patch.Init( **it, data_.phi, *data_.bnd);
            Assert( patch.Intersects(), DROPSErrCL("InitZeroExactCL::AssociateTriangles: Tetra is not intersected"), DebugNumericC);
            for ( int ch=0; ch<8; ++ch){
                if (!patch.ComputeForChild( ch)) continue; // Child ch has no intersection
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri) {
                    const size_t posOfDistTriang= distTriangPos_[ DistTriangIndexHelperCL(&(**it), ch, tri)];
                    dofToDistTriang_[dof].push_back( posOfDistTriang);
                }
            }
        }
    }
}

void InitZeroExactCL::DisplayMem() const
{
    const size_t memPerDistTriang= 9+10;
    const size_t memPerDistHelper= 8+4+4;
    const size_t memPerTriangPos= 3*8+memPerDistHelper+8;
    const size_t memPerTetraInSet= 4*8;

    size_t memDofToTetra=0, memDofToDist=0;

    for ( size_t i=0; i<dofToTetra_.size(); ++i)
        memDofToTetra += memPerTetraInSet*dofToTetra_[i].size();
    for ( size_t i=0; i<dofToDistTriang_.size(); ++i)
        memDofToDist+= 8*dofToDistTriang_[i].size();

    const size_t memByte=   distTriang_.size()*memPerDistTriang     // all distance Triangles
                          + distTriangPos_.size()*memPerTriangPos   // find position of triangle
                          + memDofToTetra                           // dof -> {intersected tera}
                          + memDofToDist;                           // dof -> {distance Triangles}
    double mem=(double)memByte/1024/1024;
#ifdef _PAR
    mem= ProcCL::GlobalMax(mem);
#endif
    std::cout << " * Add. Memory for ExactDistance: " <<mem<< " MB\n";
}

void InitZeroExactCL::DetermineDistances()
{
#ifdef COUNTMEM
    DisplayMem();
#endif
#pragma omp parallel
{
    Point3DCL *perp   = data_.UsePerp() ? new Point3DCL() : 0;
    byte      *locPerp= data_.UsePerp() ? new byte() : 0;
    double distance, newDistance;
    const int augm_size= data_.per ? data_.augmIdx->NumUnknowns() : data_.phi.Data.size();
#pragma omp for
    for (int augm_dof=0; augm_dof<augm_size; ++augm_dof){
        const int dof= data_.Map(augm_dof);
        if ( !dofToDistTriang_[dof].empty()){
            distance= std::numeric_limits<double>::max();
            for ( size_t posTri=0; posTri<dofToDistTriang_[dof].size(); ++posTri){
                newDistance= distTriang_[ dofToDistTriang_[dof][posTri]].dist( data_.coord[augm_dof], locPerp, perp);
                if ( newDistance<distance){
                    if ( data_.UsePerp() && *locPerp==2)
                        data_.UpdatePerp( dof, distance, *perp);
                    distance= newDistance;
                }
            }
            if ( data_.typ[dof] != ReparamDataCL::Finished) {
                data_.phi.Data[dof]= distance;
                data_.typ[dof]= ReparamDataCL::Finished;
            }
            else{
                data_.phi.Data[dof]= std::min( data_.phi.Data[dof], distance);
            }
        }
    }
    if (perp)    delete perp;    perp=0;
    if (locPerp) delete locPerp; locPerp=0;
}
}

void InitZeroExactCL::Clean()
{
    distTriang_.clear();
    distTriangPos_.clear();
    dofToTetra_.clear();
    dofToDistTriang_.clear();
}

/** Note, that the correctness of this method depends on the geometry of the tetras. If
    the tetra-neighborhood of a vertex contains tetras of very different size, the
    distance can be of by O(h) without noticing.
*/
void InitZeroExactCL::Perform()
{
#ifdef _PAR
    InitDofToTetra();
    BuildDistTriang();
    AssociateTriangles();
    DetermineDistances();
#else
    // This is the implementation of J. Grande, which is much faster for the sequential version
    VecDescCL oldv( data_.phi);
    InterfaceTriangleCL patch;
    double dd;
    TetraNeighborCL tetra_to_idx( data_.mg, oldv, *data_.bnd, data_.map);
    for (TetraNeighborCL::TetraToIdxMapT::iterator it= tetra_to_idx().begin(); it!=tetra_to_idx().end(); ++it) {
        const TetraCL* t( it->first);
        patch.Init( *t, oldv, *data_.bnd);
        for (int ch= 0; ch < 8; ++ch) {
            if (!patch.ComputeForChild( ch)) continue; // Child ch has no intersection

            TetraNeighborCL::IdxSetT& idxset= it->second;
            for (int tri= 0; tri < patch.GetNumTriangles(); ++tri) {
                ExactDistanceInitCL dist_to_tri( &patch.GetPoint( tri));
                for (TetraNeighborCL::IdxSetT::iterator sit= idxset.begin(); sit!=idxset.end(); ++sit) {
                    dist_to_tri.dist( data_.coord[*sit], dd);
                    const int dof= data_.Map(*sit);
                    if ( data_.typ[dof] != ReparamDataCL::Finished) {
                        data_.phi.Data[dof]= dd;
                        data_.typ[dof]= ReparamDataCL::Finished;
                    }
                    else
                        data_.phi.Data[dof]= std::min( data_.phi.Data[dof], dd);
                }
            }
        }
    }
#endif
    Clean();
}

// P A R  I N I T  Z E R O  E X A C T  C L
//----------------------------------------

#ifdef _PAR

ReparamDataCL*                            ParInitZeroExactCL::actualData_     = 0;
ParInitZeroExactCL::ToSendDistTriangMapT* ParInitZeroExactCL::toSend_         = 0;
InitZeroExactCL::DofToDistTriangT*        ParInitZeroExactCL::distTriangs_    = 0;
std::map<int, size_t>*                    ParInitZeroExactCL::actualOffset_   = 0;
size_t                                    ParInitZeroExactCL::maxTriangsPerDOF_=0;

extern "C" int ExecGatherDistTriangVertexC(OBJT objp){
    return ParInitZeroExactCL::ExecGatherDistTriang<VertexCL>(objp);
}
extern "C" int ExecGatherDistTriangEdgeC(OBJT objp){
    return ParInitZeroExactCL::ExecGatherDistTriang<EdgeCL>(objp);
}
extern "C" int HandlerDistTriangGatherPosVertexC(OBJT objp, void* buf){
    return ParInitZeroExactCL::HandlerDistTriangGatherPos<VertexCL>(objp, buf);
}
extern "C" int HandlerDistTriangGatherPosEdgeC(OBJT objp, void* buf){
    return ParInitZeroExactCL::HandlerDistTriangGatherPos<EdgeCL>(objp, buf);
}
extern "C" int HandlerDistTriangScatterPosVertexC(OBJT objp, void* buf){
    return ParInitZeroExactCL::HandlerDistTriangScatterPos<VertexCL>(objp, buf);
}
extern "C" int HandlerDistTriangScatterPosEdgeC(OBJT objp, void* buf){
    return ParInitZeroExactCL::HandlerDistTriangScatterPos<EdgeCL>(objp, buf);
}

/** Use the DDD Handler  ExecGatherDistTriang, but we do not distinguish between
    neighbor processes to to lag of DDD functions to whom the data are communicated
 */
void ParInitZeroExactCL::GatherDistTriang()
{
    actualData_= &data_; toSend_= &toSendDistTriang_; distTriangs_= &dofToDistTriang_;
    maxTriangsPerDOF_=0;
    DynamicDataInterfaceCL::IFExecLocal( InterfaceCL<VertexCL>::GetIF(), &ExecGatherDistTriangVertexC);
    DynamicDataInterfaceCL::IFExecLocal( InterfaceCL<EdgeCL>::GetIF(),   &ExecGatherDistTriangEdgeC);
    actualData_= 0; toSend_= 0; distTriangs_= 0;
}

void ParInitZeroExactCL::CommunicateDistTriang()
{
    // fill send buffer
    const size_t dPerDistTriang= 3*3+10;
    VectorCL sendBuf( toSendDistTriang_.size()*dPerDistTriang);
    size_t i=0;
    for ( ToSendDistTriangMapT::const_iterator it=toSendDistTriang_.begin(); it!=toSendDistTriang_.end(); ++it, ++i){
        // Put tri_ into the buffer
        for ( int t=0; t<3; ++t) for ( int d=0; d<3; ++d)
                sendBuf[ i*dPerDistTriang+t*3+d]= distTriang_[*it].GetTri(t)[d];
        // Put QR decomposition into the buffer
        distTriang_[*it].GetQR().Serialize( Addr(sendBuf)+i*dPerDistTriang+9);
    }
    // send data to neighbors
    const ExchangeCL::ProcNumCT neighs= data_.phi.RowIdx->GetEx().GetNeighbors();
    std::vector<ProcCL::RequestT> req( data_.phi.RowIdx->GetEx().GetNumNeighs());
    size_t pos=0;
    for ( ExchangeCL::ProcNumCT::const_iterator it= neighs.begin(); it!=neighs.end(); ++it, ++pos){
        req[pos]= ProcCL::Isend( sendBuf, *it, 1501);
    }
    // receive data from neighbors
    ProcCL::StatusT stat;
    pos=0;
    VectorCL recvBuf;
    size_t receivedElements=0;
    for ( ExchangeCL::ProcNumCT::const_iterator it= neighs.begin(); it!=neighs.end(); ++it, ++pos){
        ProcCL::Probe(*it, 1501, stat);
        const int numData= ProcCL::GetCount<GIDT>(stat);
        recvBuf.resize( numData);
        ProcCL::Recv( recvBuf, *it, 1501);
        offset_[*it]= receivedElements;
        receivedElements+= numData/dPerDistTriang;
        // Handle Data
        for ( size_t i=0; i<numData/dPerDistTriang; ++i){
            distTriang_.push_back( InitZeroExactCL::DistanceTriangCL(Addr(recvBuf)+i*dPerDistTriang));
        }
    }
    ProcCL::WaitAll( req);
}

void ParInitZeroExactCL::AssociateTrianglesOnProcBnd()
{
    actualData_= &data_; toSend_=&toSendDistTriang_;
    distTriangs_=&dofToDistTriang_; actualOffset_= &offset_;
    maxTriangsPerDOF_= ProcCL::GlobalMax( maxTriangsPerDOF_);
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<VertexCL>::GetIF(), (maxTriangsPerDOF_+2)*sizeof(IdxT),
                   HandlerDistTriangGatherPosVertexC, HandlerDistTriangScatterPosVertexC);
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<EdgeCL>::GetIF(), (maxTriangsPerDOF_+2)*sizeof(IdxT),
                   HandlerDistTriangGatherPosEdgeC, HandlerDistTriangScatterPosEdgeC);
    actualData_=0; toSend_=0;
    distTriangs_=0; actualOffset_=0;
}

/** Call of base::Clean as well*/
void ParInitZeroExactCL::Clean()
{
    base::Clean();
    toSendDistTriang_.clear();
    offset_.clear();
}

void ParInitZeroExactCL::DisplayMem() const
{
    const size_t memPerSendDist= 3*8+8, memPerOffset=3*8+4+8;
    const size_t memByte=toSendDistTriang_.size()*memPerSendDist + offset_.size()*memPerOffset;
    const double mem= ProcCL::GlobalMax((double)memByte)/1024/1024;
    std::cout << " * Add. memory for parallel implementation of ExactDistance: " << mem << " MB" << std::endl;
}

void ParInitZeroExactCL::Perform()
{
    // same as sequential
    base::InitDofToTetra();
    base::BuildDistTriang();
    base::AssociateTriangles();

    // communication:
    GatherDistTriang();
    CommunicateDistTriang();
    AssociateTrianglesOnProcBnd();

    // determine distances
    DetermineDistances();

#ifdef COUNTMEM
    DisplayMem();
#endif

    Clean();
}

#endif

// I N I T  Z E R O  P 2  C L
//---------------------------

void InitZeroP2CL::BuildRepTetra()
{
    const MultiGridCL& mg= data_.mg;
    const Uint         lvl= data_.phi.GetLevel();
    LocalP1CL<Point3DCL> Gref[10];
    InterfacePatchCL patch;
    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it){
        patch.Init( *it, data_.phi, *data_.bnd);
        if ( patch.Intersects()){
            storePosOfTetra_[&*it]= tetras_.size();
            tetras_.push_back( RepTetra(*it, Gref, data_));
        }
    }
}

/** Assign each dof located in the vicinity of the interface the position
    of close tetras in tetras_. */
void InitZeroP2CL::AssociateTetras()
{
    dofToRepTetra_.resize( data_.phi.Data.size());
    for ( size_t dof=0; dof<base::dofToTetra_.size(); ++dof){
        for ( base::TetraSetT::const_iterator it= base::dofToTetra_[dof].begin(); it!=base::dofToTetra_[dof].end(); ++it){
            dofToRepTetra_[dof].push_back( storePosOfTetra_[*it]);
        }
    }
}

void InitZeroP2CL::Clean() {}
void InitZeroP2CL::DisplayMem() const {}

void InitZeroP2CL::DetermineDistances()
{
#ifdef MINPACK
    const int n=4, lwa=50;
    int info=0;
    double x[4], f[4], tol=1e-8, wa[50];
    for (size_t dof=0; dof<dofToRepTetra_.size(); ++dof){
        for ( size_t tetra=0; tetra<dofToRepTetra_[dof].size(); ++tetra){
            actualTetra_= &( tetras_[dofToRepTetra_[dof][tetra]]);
            std::copy( actualTetra_->baryCenter.begin(), actualTetra_->baryCenter.end(), x);
            x[3]= ( actualTetra_->baryCenter-data_.coord[dof]).norm();
            actualDOF_= dof;
            info= hybrd1 ( self::f_P2, this, 4, x, f, tol, wa, lwa);

            const Point3DCL p( x, x+3);
            if ( info!=1){
                std::cerr << "Warning InitZeroP2CL::DetermineDistances: MINPACK info " << info
                          << ", tol " << tol << ", p " << p << ", f(p,lambda)= "
                          << f[0] << ' ' << f[1] << ' ' << f[2] << ' ' << f[3] << std::endl;
            }
            BaryCoordCL pbary( actualTetra_->w2b(p)); // barycentric coordinates of p
            const bool perpInsideTetra= ( pbary[0]>=-DoubleEpsC && pbary[1]>=-DoubleEpsC && pbary[2]>=-DoubleEpsC && pbary[3]>=-DoubleEpsC);
            double newDist = perpInsideTetra ? (p-data_.coord[ dof]).norm() : data_.phi.Data[ dof];

            if (!perpInsideTetra){

            }

        }
    }
    actualTetra_= 0;
#else
    throw DROPSErrCL("InitZeroP2CL::DetermineDistances: Cannot use MINPACK");
#endif
}


int InitZeroP2CL::f_P2(void *this_class, int, const double *x, double *fvec, int)
{
    InitZeroP2CL* actual= static_cast<InitZeroP2CL*>(this_class);
    const Point3DCL p( x, x+3);                         // actual point p (given by optimizer)
    Point3DCL Gphi;                                     // gradient of phi at point p
    const double lambda= x[3];                          // value of lambda (given by optimizer)
    BaryCoordCL pbary= actual->actualTetra_->w2b(p);    // barycentric coordinates of p

    // Compute gradient of phi at point p
    for ( Uint i=0; i<10; ++i)
        Gphi += actual->actualTetra_->valPhi[i]*actual->actualTetra_->G[i](pbary);

    // Compute f
    fvec[0]= lambda*Gphi[0] + p[0] - actual->data_.coord[actual->actualDOF_][0];
    fvec[1]= lambda*Gphi[1] + p[1] - actual->data_.coord[actual->actualDOF_][1];
    fvec[2]= lambda*Gphi[2] + p[2] - actual->data_.coord[actual->actualDOF_][2];
    fvec[3]= P2EvalT::val( actual->actualTetra_->valPhi, pbary);
    return 1;
}

void InitZeroP2CL::Perform()
{
    base::InitDofToTetra();

}

// F A S T M A R C H I N G  C L
//-----------------------------

/** Iterate over all tetrahedra and store a child tetrahedra representation of each vertex*/
void FastmarchingCL::InitNeigh()
{
    neigh_.resize( data_.phi.Data.size());
    const Uint lvl= data_.phi.GetLevel();
    const Uint idx= data_.per ? data_.augmIdx->GetIdx() : data_.phi.RowIdx->GetIdx();
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);
    IdxT Numb[10];

    DROPS_FOR_TRIANG_TETRA( base::data_.mg, lvl, it){
        for ( int v=0; v<10; ++v){ // collect data on all DoF
            if (v<4)
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            else
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);
        }
        for ( Uint ch=0; ch<MaxChildrenC; ++ch){
            const ChildDataCL chdata= GetChildData( RegRef.Children[ch]);
            // build tetra
            ReprTetraT t;
            for ( Uint vert=0; vert<NumVertsC; ++vert)
                t[vert]= Numb[ chdata.Vertices[vert]];
            // store tetra
            for ( Uint vert=0; vert<NumVertsC; ++vert)
                neigh_[ base::data_.Map( t[vert])].push_back( t);
        }
    }
}

/** Iterate over all vertices and edges and check if a dof is marked as finished.
    If so, put all neighbor dof in the close set.
    \pre InitNeigh has to be called
*/
void FastmarchingCL::InitClose()
{
    std::set<IdxT> closeVerts;
    for (size_t dof=0; dof<data_.typ.size(); ++dof){
        if ( data_.typ[dof] == data_.Finished)
            for (Uint n = 0; n < neigh_[dof].size(); ++n)
                for (int j = 0; j < 4; ++j)
                    closeVerts.insert( neigh_[dof][n][j]);
    }
    std::set<IdxT>::const_iterator closeit;
    for ( closeit=closeVerts.begin(); closeit!=closeVerts.end(); ++closeit){
        Update( *closeit);
    }
#if DROPSDebugC&DebugNumC
    // Check if all values are >0
    for (size_t i=0; i<data_.phi.Data.size(); ++i){
        Assert(data_.phi.Data[i]>=0, DropsErrCL("FastmarchingCL::InitClose: Values must be greater than zero"), ~0);
    }
#endif
}

/** Iterate over all neighbors of NrI and check if distance has changed. If a neighbor
    is not in the close set, put this dof in this set. If NrI is already marked as finished
    then do nothing
    \param NrI dof to be updated
*/
void FastmarchingCL::Update( const IdxT NrI)
{
    const IdxT MapNrI = data_.Map(NrI);

    // Update all vertices that are not Finished
    if ( data_.typ[MapNrI] == data_.Finished)
        return;

    IdxT upd[3];
    double minval = ( data_.typ[MapNrI] == data_.Close) ? data_.phi.Data[MapNrI] : 1e99;

    // Update all neighbor vertices
    for ( Uint n = 0; n < neigh_[MapNrI].size(); ++n) {
        int num = 0;
        for ( int j = 0; j < 4; ++j) {
            const IdxT NrJ   = neigh_[MapNrI][n][j];
            const IdxT MapNrJ= data_.Map( NrJ);
            if ( data_.typ[MapNrJ] == data_.Finished) {
                upd[num++] = NrJ;
                minval = std::min(minval,
                        data_.phi.Data[ MapNrJ]+(data_.coord[NrJ]-data_.coord[NrI]).norm());
            }
        }
        minval = std::min(minval, CompValueProj(NrI, num, upd));
    }

    data_.phi.Data[MapNrI] = minval;
    if (data_.typ[MapNrI] != data_.Close) {
        close_.insert( DistIdxT( minval, MapNrI));
        data_.typ[MapNrI] = data_.Close;
    }
}

/** Compute the projection to an edge or an face
    \param Nr the dof to be computed
    \param num number of intersection points (2: edge, 3: face)
    \param upd corners of edge/face
    \return distance of Nr to edge/face
*/
double FastmarchingCL::CompValueProj( IdxT Nr, int num, const IdxT upd[3]) const
{
#ifdef _PAR
    if (data_.per)
        throw DROPSErrCL("FastmarchingCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif
    double val= 1e99;
    VectorBaseCL<Point3DCL>& coord_=data_.coord;

    switch (num){
        case 2:{ // projection on edge
            const Point3DCL a= coord_[upd[1]] - coord_[upd[0]];
            const Point3DCL b= coord_[  Nr  ] - coord_[upd[0]];
            double bary= inner_prod(a,b)/a.norm_sq();

            data_.Normalize(bary);
            const Point3DCL lotfuss= (1-bary)*coord_[upd[0]] + bary*coord_[upd[1]];
            const double y= (1-bary)*data_.phi.Data[ data_.Map(upd[0])] + bary*data_.phi.Data[ data_.Map(upd[1])];
            val= y + (lotfuss - coord_[Nr]).norm();
        }
        break;

        case 3:{ // projection of face
            const Point3DCL a= coord_[upd[1]] - coord_[upd[0]];
            const Point3DCL b= coord_[upd[2]] - coord_[upd[0]];
            const Point3DCL c= coord_[  Nr  ] - coord_[upd[0]];
            double bary1= inner_prod(a,c)/a.norm_sq(),
                   bary2= inner_prod(b,c)/b.norm_sq();
            data_.Normalize(bary1, bary2);
            const Point3DCL lotfuss= (1-bary1-bary2)*coord_[upd[0]] + bary1*coord_[upd[1]] + bary2*coord_[upd[2]];
            const double y= (1-bary1-bary2)*data_.phi.Data[data_.Map(upd[0])]
                   + bary1*data_.phi.Data[data_.Map(upd[1])]+bary2*data_.phi.Data[data_.Map(upd[2])];
            val= y + (lotfuss - coord_[Nr]).norm();
        }
    }
    return val;
}

/** While some vertices are still marked as Close, determine the distance of these
    vertices by the FMM
*/
void FastmarchingCL::DetermineDistances()
{
#ifdef COUNTMEM
    size_t elemClose=0, memPerClose=28;
    size_t elemRepTetra=0, memPerTetra=4*8;
    for ( size_t i=0; i<neigh_.size(); ++i)
        elemRepTetra+= neigh_[i].size();
#endif
    IdxT next;

    while ( !close_.empty()) {
#ifdef COUNTMEM
        elemClose= std::max( elemClose, close_.size());
#endif
        // remark: next < size_   =>   Map not needed for next
        next= close_.GetNearest().second;
        data_.typ[next] = data_.Finished;

        std::set<IdxT> neighVerts;
        for ( Uint n = 0; n < neigh_[next].size(); ++n) { // collect all neighboring verts in neighVerts
            for ( Uint i = 0; i < 4; ++i)
                neighVerts.insert(neigh_[next][n][i]);
        }
        for (std::set<IdxT>::const_iterator it = neighVerts.begin(), end = neighVerts.end(); it != end; ++it) { // update all neighboring verts, mark as Close
            Update( *it);
        }
        neigh_[next].clear(); // will not be needed anymore
    }

#ifdef COUNTMEM
    usedMem_=  elemClose*memPerClose                // elements in close
             + elemRepTetra*memPerTetra             // tetrahedra to represent neighbors
             + (8+3*8+1)* data_.phi.Data.size();   // value of phi + coord + type
    const double mem= (double)usedMem_/1024/1024;
    std::cout << " * Used memory for Fastmarching: " << mem << " MB" << std::endl;
#endif
}

/** Apply the FMM to a level set function*/
void FastmarchingCL::Perform()
{
#ifdef _PAR
    if (data_.per)
        throw DROPSErrCL("FastmarchingCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif
    InitNeigh();
    InitClose();
    DetermineDistances();
}

#ifdef _PAR

// F A S T M A R C H I N G  O N  M A S T E R  C L
//-----------------------------------------------
ReparamDataCL*    FastmarchingOnMasterCL::actualData_= 0;
std::vector<IdxT> FastmarchingOnMasterCL::globNumb_  = std::vector<IdxT>();
std::vector<IdxT> FastmarchingOnMasterCL::locNumb_   = std::vector<IdxT>();


/** Create the mapping local -> global and global -> local dof*/
void FastmarchingOnMasterCL::CreateGlobNumb()
{
    Comment("Create global numbering\n", DebugParallelNumC);
    const IdxT numDoF   = data_.phi.RowIdx->NumUnknowns();
    const IdxT idx      = data_.phi.RowIdx->GetIdx();
    const Uint lvl      = data_.phi.GetLevel();
    IdxT numExclusiveDoF= 0;
    globNumb_.resize(0); globNumb_.resize(numDoF, NoIdx);
    size_= numDoF;

    // If owning an exclusive DoF, put them into list, otherwise this dof is NoIdx
    DROPS_FOR_TRIANG_VERTEX( data_.mg, lvl, it){
        if ( it->IsExclusive(PrioHasUnk)){
            globNumb_[ it->Unknowns(idx)] = numExclusiveDoF++;
        }
    }
    DROPS_FOR_TRIANG_EDGE( data_.mg, lvl, it){
        if ( it->IsExclusive(PrioHasUnk)){
            globNumb_[ it->Unknowns(idx)] = numExclusiveDoF++;
        }
    }
    locNumb_.resize( numExclusiveDoF, NoIdx);

    // Get offset for each process
    exclusive_.resize( ProcCL::Size());
    ProcCL::Gather( (int)numExclusiveDoF, Addr( exclusive_), -1);
    offset_.resize(ProcCL::Size()+1);
    for (int i=0; i<ProcCL::Size(); ++i)
        offset_[i+1]  = offset_[i] + exclusive_[i];

    // Append offset to each global number (so it is unique)
    for (IdxT i=0; i<numDoF; ++i){
        if (globNumb_[i]!=NoIdx)
            globNumb_[i] += offset_[ ProcCL::MyRank()];
    }

    // Now, get global numbers of distributed DOF not stored on this process
    actualData_= &data_;        // make data static accessible, so DDD can access it
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<VertexCL>::GetIF(), sizeof(IdxT),
                   HandlerGlobDOFGatherVertexC, HandlerGlobDOFScatterVertexC );
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<EdgeCL>::GetIF(), sizeof(IdxT),
                   HandlerGlobDOFGatherEdgeC, HandlerGlobDOFScatterEdgeC );
    actualData_= 0;

    // in debug mode, check if everything is right
#if DROPSDebugC&DebugParallelNumC
    DROPS_FOR_TRIANG_VERTEX( data_.mg, lvl, it)
        if (globNumb_[ it->Unknowns(idx) ]==NoIdx)
            it->DebugInfo(std::cout);
    DROPS_FOR_TRIANG_EDGE( data_.mg, lvl, it)
        if (globNumb_[ it->Unknowns(idx) ]==NoIdx)
            it->DebugInfo(std::cout);

    for (IdxT i=0; i<numDoF; ++i)
        if (globNumb_[i]==NoIdx)
            throw DROPSErrCL("FastMarchCL::CreateGlobNumb: Not all local DoF numbers mapped on global DoF number");

    IdxT maxEntry=*std::max_element(globNumb_.begin(), globNumb_.end());
    if (maxEntry>(IdxT)offset_[ProcCL::Size()]){
        std::cerr << "["<<ProcCL::MyRank()<<"] max entry is "<<maxEntry<<" but there are only "
                  << offset_[ProcCL::Size()] << " global dof!" << std::endl;
        throw DROPSErrCL("FastMarchCL::CreateGlobNumb: max entry in globNumb is bigger than global size");
    }
#endif
}

extern "C" int HandlerGlobDOFGatherVertexC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerGlobDOFGather<VertexCL>(objp, buf);
}
extern "C" int HandlerGlobDOFGatherEdgeC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerGlobDOFGather<EdgeCL>(objp, buf);
}
extern "C" int HandlerGlobDOFScatterVertexC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerGlobDOFScatter<VertexCL>(objp,buf);
}
extern "C" int HandlerGlobDOFScatterEdgeC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerGlobDOFScatter<EdgeCL>(objp,buf);
}

/** this function uses the DDD-Interfaces of master vertices and edges */
void FastmarchingOnMasterCL::DistributeFinished()
{
    Comment("Distribute Finished\n", DebugParallelNumC);

    actualData_= &data_;
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<VertexCL>::GetIF(),  sizeof(CoupMarkValST),
                   HandlerFinishedGatherVertexC,   HandlerFinishedScatterVertexC );
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<EdgeCL>::GetIF(), sizeof(CoupMarkValST),
                   HandlerFinishedGatherEdgeC,   HandlerFinishedScatterEdgeC );
    actualData_=0 ;
}

extern "C" int HandlerFinishedGatherVertexC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerFinishedGather<VertexCL>(objp,buf);
}
extern "C" int HandlerFinishedGatherEdgeC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerFinishedGather<EdgeCL>(objp,buf);
}
extern "C" int HandlerFinishedScatterVertexC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerFinishedScatter<VertexCL>(objp,buf);
}
extern "C" int HandlerFinishedScatterEdgeC(OBJT objp, void* buf){
    return FastmarchingOnMasterCL::HandlerFinishedScatter<EdgeCL>(objp,buf);
}

/** Get all local data of types, coordinates, values and tetrahedra
    Here, as well, the mapping "send/receive element -> local number"
    is created
    \pre Global numbering must exists, i.e., CreateGlobalNumb() has
         to be called first
*/
void FastmarchingOnMasterCL::CollectLocalData(std::vector<byte>& typ, std::vector<double>& coord, std::vector<double>& values, std::vector<IdxT>& tetraList) const
{
    const Uint idx= data_.phi.RowIdx->GetIdx();
    const Uint lvl= data_.phi.GetLevel();

    // collect all types, coordinates and values on exclusive vertices
    // and edges, and as well, create mapping the mapping localNumb_ for
    // later receiving
    typ.clear();
    coord.clear();
    values.clear();
    locNumb_.clear();
    DROPS_FOR_TRIANG_VERTEX( data_.mg, lvl, it){
        if ( it->IsExclusive(PrioHasUnk)){
            const IdxT dof= it->Unknowns(idx);
            typ.push_back( data_.typ[ dof]);
            for ( int i=0; i<3; ++i)
                coord.push_back( data_.coord[ dof][i]);
            values.push_back( data_.phi.Data[ dof]);
            locNumb_.push_back( dof);
        }
    }
    DROPS_FOR_TRIANG_EDGE( data_.mg, lvl, it){
        if ( it->IsExclusive(PrioHasUnk)){
            const IdxT dof= it->Unknowns(idx);
            typ.push_back( data_.typ[ dof]);
            for ( int i=0; i<3; ++i)
                coord.push_back( data_.coord[ dof][i]);
            values.push_back( data_.phi.Data[ dof]);
            locNumb_.push_back( dof);
        }
    }

    // collect all tetrahedra
    tetraList.clear();
    IdxT Numb[10];
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);
    DROPS_FOR_TRIANG_TETRA( data_.mg, lvl, it){
        for ( int v=0; v<10; ++v){
            if (v<4)
                Numb[v]= globNumb_[ it->GetVertex(v)->Unknowns(idx)];
            else
                Numb[v]= globNumb_[ it->GetEdge(v-4)->Unknowns(idx)];
        }
        for ( Uint ch=0; ch<MaxChildrenC; ++ch){
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            for ( Uint vert= 0; vert<NumVertsC; ++vert)
                tetraList.push_back( Numb[ data.Vertices[vert]]);
        }
    }
}

/** At first, collect all local data such as values, types and tetrahedra.
    Afterwards, send these data to the master process, who is responsible
    for performing the FMM.
    The master process puts the received data into the member data_. This
    function also initializes the neighbor sets.
*/
void FastmarchingOnMasterCL::Collect()
{
    // Get local data
    std::vector<byte>   localTyp;
    std::vector<double> localCoord;
    std::vector<double> localValues;
    std::vector<IdxT>   localTetraList;
    CollectLocalData(localTyp, localCoord, localValues, localTetraList);

    // Gather all data at master process
    std::valarray<byte>   types    = ProcCL::Gatherv( localTyp,       master_);
    std::valarray<double> coord    = ProcCL::Gatherv( localCoord,     master_);
    std::valarray<double> values   = ProcCL::Gatherv( localValues,    master_);
    std::valarray<IdxT>   tetraList= ProcCL::Gatherv( localTetraList, master_);

    Assert( types.size()    ==  values.size(), DROPSErrCL("FastmarchingOnMasterCL::Collect: types size does not match"), DebugParallelNumC);
    Assert( coord.size()    ==3*values.size(), DROPSErrCL("FastmarchingOnMasterCL::Collect: coord size does not match"), DebugParallelNumC);

    if ( ProcCL::MyRank()==master_){
        // Transform data on master, i.e., put collected data into the member data_
        data_.phi.Data.resize( values.size());
        std::copy( Addr(values), Addr(values)+values.size(), Addr(data_.phi.Data));

        data_.typ.resize( types.size());
        std::copy( Addr(types), Addr(types)+types.size(), Addr(data_.typ));

        data_.coord.resize( values.size());
        for ( size_t i=0; i<data_.coord.size(); ++i) for ( int j=0; j<3; ++j)
            data_.coord[i][j]= coord[i*3+j];

        // Initialize list of neighbors
        base::neigh_.resize(0); base::neigh_.resize( values.size());
        for ( size_t t=0; t<tetraList.size(); t+=4){
            ReprTetraT tetra;
            for (int i=0; i<4; ++i)
                tetra[i]= tetraList[t+i];
            for (int i=0; i<4; ++i)
                base::neigh_[ tetra[i]].push_back(tetra);
        }
    }
}

/** After the master has performed the FMM, the values of the level set
    has to be send back to all other processes.
*/
void FastmarchingOnMasterCL::Distribute()
{
    // Send values back
    VectorCL recvBuf( exclusive_[ProcCL::MyRank()]);
    ProcCL::Scatterv( Addr(data_.phi.Data), Addr(exclusive_), Addr(offset_), Addr(recvBuf), size_, master_);

    // Put data back into data_.phi.Data
    data_.phi.Data.resize( size_);
    data_.phi.Data= 0.;
    for ( size_t i=0; i<recvBuf.size(); ++i)
        data_.phi.Data[ locNumb_[i]]=recvBuf[i];

    // Data on distributed DOF has just been send to one process, so accumulate the vector
    data_.phi.RowIdx->GetEx().Accumulate( data_.phi.Data);
}

/** Free memory of static members*/
void FastmarchingOnMasterCL::CleanUp()
{
    actualData_= 0;
    globNumb_.clear();
    locNumb_.clear();
}

void FastmarchingOnMasterCL::Perform()
{
    DistributeFinished();
    CreateGlobNumb();
    Collect();
    if (ProcCL::MyRank()==master_){
        base::InitClose();
        base::DetermineDistances();
    }
    Distribute();
    CleanUp();
}
#endif      // of _PAR

// D I R E C T  D I S T A N C E  C L
//----------------------------------

void DirectDistanceCL::InitFrontVector()
/** Count number of frontier vertices and perpendiculart feet, put
    the coordinates into the the vector front_ and the value on
    these points into the vector vals_
*/
{
    // Count number of frontier vertices and perpendiccular feet
    Uint numFront= 0;
    for ( size_t i=0; i<data_.typ.size(); ++i)
        if ( data_.typ[i]==ReparamDataCL::Finished)
            ++numFront;
    for ( size_t i=0; i<data_.perpFoot.size() && data_.UsePerp(); ++i)
        if ( data_.perpFoot[i]!=0)
            ++numFront;

    // Put all frontier vertices and perpendicular feet into a single array
    // and store distance of each frontier vertex
    front_.resize(0); front_.resize( 3*numFront);
    vals_.resize(0);  vals_.resize( numFront);
    size_t pos=0;
    size_t posDist=0;
    for ( size_t i=0; i<data_.typ.size(); ++i){
        if ( data_.typ[i]==ReparamDataCL::Finished){
            for ( int j=0; j<3; ++j){
                front_[ pos++]= data_.coord[i][j];
            }
            vals_[ posDist++]= data_.phi.Data[i];
        }
    }
    for ( size_t i=0; i<data_.perpFoot.size() && data_.UsePerp(); ++i){
        if ( data_.perpFoot[i]!=0){
            for ( int j=0; j<3; ++j){
                front_[ pos++]= (*data_.perpFoot[i])[j];
            }
            // we do not need to set distance to 0, because std::valarray is initialized by 0
        }
    }
#ifndef _PAR
    Uint numLsetUnk=data_.phi.Data.size();
    std::cout << " * Lset unk " << numLsetUnk
              << ", frontier " << posDist
              << ", perpFoot " << (numFront-posDist) << std::endl;
#endif
}

void DirectDistanceCL::BuildKDTree()
/** Take the elements out of front_ and build a kd-tree representing this set
    \pre InitFrontVector has to be called
 */
{
    kdTree_= new KDTreeCL<double>( Addr(front_), front_.size()/3, 3, true);
}

void DirectDistanceCL::DetermineDistances()
/** Iterate over all off-site vertices and assign shortest distance to a
    frontier vertex or perpendicular foot to phi
*/
{
#pragma omp parallel for
    for ( int dof=0; dof<(int)data_.phi.Data.size(); ++dof) {
        if ( data_.typ[dof]!=ReparamDataCL::Finished && data_.typ[dof]!=ReparamDataCL::Handled) {
            double newPhi= std::numeric_limits<double>::max();
            KDTreeResultVectorCL<double> e;
            std::vector<double> qv( data_.coord[dof].begin(), data_.coord[dof].end());
            kdTree_->GetNNearest( qv, numNeigh_, e);
            // check for smallest distance in neighbor vertices
            for ( int n=0; n<(int)e.size(); ++n) {
                newPhi= std::min( newPhi, std::sqrt( e[n].Distance)+vals_[ e[n].Index]);
            }
#pragma omp critical
            data_.phi.Data[dof]= newPhi;
        }
    }
}

void DirectDistanceCL::DisplayMem() const
{
    const size_t memFront= front_.size()*8, memVals= vals_.size()*8;
    const size_t memPerKDInterval=2*8;
    const size_t memPerKDNode= 4+3*8+2*4+2*memPerKDInterval+2*8;
    const size_t memKDTree= front_.size()*(8+memPerKDNode);
    const size_t memData= (8+3*8+1)* data_.phi.Data.size();
    const size_t memByte= memData+memFront+memVals+memKDTree;
    const double mem=double(memByte)/1024/1024;
    std::cout << " * Memory by Direct Distance " << mem << " MB.\n";
}

void DirectDistanceCL::Perform()
/** Perform all steps that are necessary to determine the value of the level set
    function on off-site vertices.
*/
{
    InitFrontVector();
    BuildKDTree();
    DetermineDistances();
    DisplayMem();
    delete kdTree_; kdTree_=0;
}

// P A R  D I R E C T  D I S T A N C E  C L
//-----------------------------------------

#ifdef _PAR
ParDirectDistanceCL::MultiFrontT ParDirectDistanceCL::onProc_=
    std::map<IdxT, std::list<ParDirectDistanceCL::TransferST> >();
ReparamDataCL* ParDirectDistanceCL::actualData_=0;

extern "C" int HandlerFrontierGatherVertexC(OBJT objp, void* buf){
    return ParDirectDistanceCL::HandlerFrontierGather<VertexCL>(objp,buf);
}
extern "C" int HandlerFrontierGatherEdgeC(OBJT objp, void* buf){
    return ParDirectDistanceCL::HandlerFrontierGather<EdgeCL>(objp,buf);
}
extern "C" int HandlerFrontierScatterVertexC(OBJT objp, void* buf){
    return ParDirectDistanceCL::HandlerFrontierScatter<VertexCL>(objp,buf);
}
extern "C" int HandlerFrontierScatterEdgeC(OBJT objp, void* buf){
    return ParDirectDistanceCL::HandlerFrontierScatter<EdgeCL>(objp,buf);
}

void ParDirectDistanceCL::CommunicateFrontierSetOnProcBnd()
{
    actualData_=&data_;
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<VertexCL>::GetIF(), sizeof(TransferST),
            HandlerFrontierGatherVertexC, HandlerFrontierScatterVertexC );
    DynamicDataInterfaceCL::IFExchange(InterfaceCL<EdgeCL>::GetIF(), sizeof(TransferST),
            HandlerFrontierGatherEdgeC, HandlerFrontierScatterEdgeC );
    actualData_=0;

    for ( MultiFrontT::iterator it=onProc_.begin(); it!=onProc_.end(); ++it){
        const IdxT dof= it->first;

        // Put own value into map
        TransferST tmp;
        tmp.value= data_.typ[dof]==ReparamDataCL::Finished ? data_.phi.Data[dof] : std::numeric_limits<double>::max();;
        tmp.procID= ProcCL::MyRank();
        tmp.perp= data_.perpFoot[ dof] ? *data_.perpFoot[dof] : Point3DCL(std::numeric_limits<double>::max());
        it->second.push_front( tmp);

        // search for process determined the smallest value of phi
        double minval= std::numeric_limits<double>::max();
        MultiFrontT::mapped_type::const_iterator minProc;
        for ( MultiFrontT::mapped_type::const_iterator pit= it->second.begin(); pit!=it->second.end(); ++pit){
            if ( pit->value<=minval){
                minval= pit->value;
                minProc= pit;
            }
        }

        data_.phi.Data[ dof]= minProc->value;
        if ( minProc->perp[0]<std::numeric_limits<double>::max()){
            if (data_.perpFoot[dof])
                *data_.perpFoot[dof]= minProc->perp;
            else
                data_.perpFoot[dof]= new Point3DCL( minProc->perp);
        }
        else{
            if (data_.perpFoot[dof]){
                delete data_.perpFoot[dof]; data_.perpFoot[dof]=0;
            }
        }
        data_.typ[dof]= minProc->procID==ProcCL::MyRank() ? ReparamDataCL::Finished : ReparamDataCL::Handled;
    }
}

void ParDirectDistanceCL::GatherFrontier()
/** Call of MPI function Allgatherv
    \pre To gather the frontier and vertices and perpendicular feet, this information must
    be collected by calling base::InitFrontVector
*/
{
    //call MPI_Allgatherv to gather values and frontier vertices
    VectorCL allVals ( ProcCL::Gatherv( base::vals_,  -1)),
             allFront( ProcCL::Gatherv( base::front_, -1));

    // local sets are not need any more ...
    base::vals_.resize( allVals.size()); vals_=allVals;
    base::front_.resize( allFront.size()); front_=allFront;

    Uint numLsetUnk= ProcCL::GlobalSum(data_.phi.Data.size());
    std::cout << " * Lset unk " << numLsetUnk
              << ", frontier and perpendicular feet " << vals_.size() << std::endl;
}

void ParDirectDistanceCL::CleanUp()
{
    onProc_.clear();
    actualData_=0;
}

void ParDirectDistanceCL::Perform()
{
    base::InitFrontVector();
    CommunicateFrontierSetOnProcBnd();
    GatherFrontier();
    base::BuildKDTree();
    base::DetermineDistances();
    delete kdTree_; kdTree_=0;
}
#endif

// R E P A R A M  C L
//-------------------

/** Pass all necessary information to reparametrization class
    \param mg       multi grid
    \param phi      values of level set function
    \param periodic periodic boundaries are used
    \param bnd      boundary conditions for periodic boundaries
*/
ReparamCL::ReparamCL( MultiGridCL& mg, VecDescCL& phi, bool gatherPerp, bool periodic, const BndDataCL<>* bnd)
  : data_( mg, phi, gatherPerp, periodic, bnd)
{ }

/** Clean everything up*/
ReparamCL::~ReparamCL()
{
    if (initZero_)
        delete initZero_;
    initZero_=0;
    if (propagate_)
        delete propagate_;
    propagate_=0;
}

/** Assign each dof the sign stored in old*/
void ReparamCL::RestoreSigns()
{
#pragma omp parallel for schedule(static)
    for ( int i=0; i<(int)data_.old.size(); ++i){
        if ( data_.old[i]<0){
            data_.phi.Data[i]*= -1.;
        }
    }
}

void ReparamCL::Perform()
{
    std::cout << "Reparametrize level set function\n";
#pragma omp parallel
{
    int numThreads=1;
#ifdef _OPENMP
    numThreads= omp_get_num_threads();
#endif
#pragma omp master
    std::cout << " * Using "<<numThreads<<" thread(s)" << std::endl;
}
#ifdef _PAR
    ParTimerCL timer, alltimer;
#else
    TimerCL timer, alltimer;
#endif

    timer.Reset();
    initZero_->Perform();
    timer.Stop();
    std::cout << " * Init frontier set by " << initZero_->GetName() << " took " << timer.GetTime() << " sec." << std::endl;

    timer.Reset();
    propagate_->Perform();
    timer.Stop();
    std::cout << " * Propagation by " << propagate_->GetName() << " took " << timer.GetTime() << " sec." << std::endl;
    RestoreSigns();

    alltimer.Stop();
    std::cout << " * Re-parametrization took " << alltimer.GetTime() << " sec." << std::endl;
}

// R E P A R A M  F A C T O R Y  C L
//----------------------------------


/** Construct the reparametrization class. For method the following arguments are valid
    \param mg       multi grid
    \param phi      values of level set function
    \param method   used methods to reparametrize a level set function
    \param periodic periodic boundaries are used
    \param bnd      boundary conditions for periodic boundaries
    \return pointer to a reparametrization class
*/
std::auto_ptr<ReparamCL> ReparamFactoryCL::GetReparam( MultiGridCL& mg,
        VecDescCL& phi, int method, bool periodic, const BndDataCL<>* bnd)
{
    int initMethod= method%10;
    int propMethod= method/10;
    std::auto_ptr<ReparamCL> reparam(new ReparamCL(mg, phi, propMethod==1, periodic, bnd));
    switch (initMethod) {
        case 0: {
            reparam->initZero_ = new InitZeroNoModCL( reparam->data_);
            break;
        }
        case 1: {
            reparam->initZero_ = new InitZeroP1CL<1>( reparam->data_);
            break;
        }
        case 2: {
            reparam->initZero_ = new InitZeroP1CL<0>( reparam->data_);
            break;
        }
        case 3: {
#ifdef _PAR
            reparam->initZero_ = new ParInitZeroExactCL( reparam->data_);
#else
            reparam->initZero_ = new InitZeroExactCL( reparam->data_);
#endif
            break;
        }
        default: {
            throw DROPSErrCL("ReparamFactoryCL::GetReparam: Unknown method for InitZero");
        }
    }
    switch (propMethod) {
        case 0: {
#ifdef _PAR
            reparam->propagate_ = new FastmarchingOnMasterCL( reparam->data_);
#else
            reparam->propagate_ = new FastmarchingCL( reparam->data_);
#endif
            break;
        }
        case 1: {
#ifdef _PAR
            reparam->propagate_ = new ParDirectDistanceCL( reparam->data_);
#else
            reparam->propagate_ = new DirectDistanceCL( reparam->data_);
#endif
            break;
        }
        default: {
            throw DROPSErrCL("ReparamFactoryCL::GetReparam: Unknown method for Propagate");
        }
    }
    return reparam;
}
} // end of namespace DROPS
