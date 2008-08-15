//**************************************************************************
// File:    multigrid.cpp                                                  *
// Content: classes that constitute the multigrid                          *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - August, 3 2000                                         *
//                                                                         *
// Remarks: We should use the const-qualifier to make it difficult to      *
//          accidentally change the multigrid structure from anywhere      *
//          outside of the multigrid algorithms.                           *
//          Thus the pointer to user data structures should probably be    *
//          a pointer to mutable.                                          *
//**************************************************************************

#include "geom/multigrid.h"
#include "num/solver.h"

namespace DROPS
{

//
// static members of TetraCL
//
SArrayCL<EdgeCL*, NumAllEdgesC> TetraCL::_ePtrs(static_cast<EdgeCL*>(0));
SArrayCL<FaceCL*, NumAllFacesC> TetraCL::_fPtrs(static_cast<FaceCL*>(0));

// E d g e C L


// Assumes that the next Level with respect to ep exists.
void EdgeCL::BuildMidVertex(VertContT& container, const BoundaryCL& Bnd)
// TODO: For nonlinear boundaries, we must treat the case, in which an edge lies in two
//       boundary-segments in a different manner: Project to the common "edge" of the
//       boundary-segments!
// TODO: Due to MeshReader-Boundaries, which have no 2D-parametrization, we calculate
//       the barycenter of the edge directly. This, of course, breaks nonlinear boundary
//       segments.
{
    const VertexCL* const vp0 ( GetVertex( 0));
    const VertexCL* const vp1 ( GetVertex( 1));

    if (IsOnBoundary()) {
        const BndIdxT bndidx= *GetBndIdxBegin();
        const BndPointCL& bndvert0= *std::find_if( vp0->GetBndVertBegin(), vp0->GetBndVertEnd(),
                                                   BndPointSegEqCL( bndidx));
        const BndPointCL& bndvert1= *std::find_if( vp1->GetBndVertBegin(), vp1->GetBndVertEnd(),
                                                   BndPointSegEqCL( bndidx));
        BndPairCL bndpair= Bnd.GetBndSeg( bndidx)->MidProject( bndvert0, bndvert1);
        // XXX: Revise this for nonlinear boundary-segments.
        bndpair.second= GetBaryCenter( *this);
        container.push_back( VertexCL( bndpair.second, GetLevel() + 1));
        SetMidVertex( &container.back());
        container.back().AddBnd( BndPointCL( bndidx, bndpair.first));
        if ( std::distance( GetBndIdxBegin(), GetBndIdxEnd()) == 2 ) {
            const BndIdxT bndidx= *(GetBndIdxBegin() + 1);
            const BndPointCL& bndvert0= *std::find_if( vp0->GetBndVertBegin(), vp0->GetBndVertEnd(),
                                                       BndPointSegEqCL( bndidx));
            const BndPointCL& bndvert1= *std::find_if( vp1->GetBndVertBegin(), vp1->GetBndVertEnd(),
                                                       BndPointSegEqCL( bndidx));
            BndPairCL bndpair1= Bnd.GetBndSeg( bndidx)->MidProject( bndvert0, bndvert1);
            container.back().AddBnd( BndPointCL( bndidx, bndpair1.first));
            container.back().BndSort();
//            Assert( bndpair.second == bndpair1.second, DROPSErrCL("BuildMidVertex: Projection leads to different 3D-coords."), ~0 );
        }
    }
    else {
        container.push_back( VertexCL( BaryCenter( vp0->GetCoord(), vp1->GetCoord()), GetLevel() + 1));
        SetMidVertex( &container.back());
    }
}

Point3DCL GetBaryCenter(const EdgeCL& e)
{
    return (e.GetVertex(0)->GetCoord() + e.GetVertex(1)->GetCoord() )*0.5;
}

// F a c e C L

// Assumptions:
// - a green child and its parent with a common face are both stored on the same side

void FaceCL::LinkTetra(const TetraCL* tp)
{
    int offset=0;

    if (tp->GetLevel()==GetLevel())
    {
        if (_Neighbors[0]) offset=1;
    }
    else
    {
        Assert(tp->GetLevel() == GetLevel()+1, DROPSErrCL("FaceCL::LinkTetra: Illegal level of green tetra"), DebugRefineEasyC);
        // tetra is stored on the same side as the parent
        offset= _Neighbors[0]==tp->GetParent() ? 2 : 3;
    }

    Assert(!_Neighbors[offset], DROPSErrCL("FaceCL::LinkTetra: Tetra linked already!"), DebugRefineEasyC);
    _Neighbors[offset]= tp;
}


void FaceCL::UnlinkTetra(const TetraCL* tp)
{
    int i;

    // find tetra in neighbor list
    for ( i=0; i<4; ++i )
        if (_Neighbors[i] == tp) break;

    Assert(i<4, DROPSErrCL("FaceCL::UnlinkTetra: No such tetra."), DebugRefineEasyC);
//    Assert(!(i<2 && IsOnNextLevel()), DROPSErrCL("FaceCL::UnlinkTetra: First unlink tetras on next level!"),
    if (i == 0)
    {
        i= 1;
        _Neighbors[0]= _Neighbors[1];
        _Neighbors[2]= _Neighbors[3];
        _Neighbors[3]= 0;
    }
    _Neighbors[i]= 0;
}

const TetraCL* FaceCL::GetNeighInTriang(const TetraCL* tp, Uint TriLevel) const
{
    Assert( tp->IsInTriang(TriLevel) && IsInTriang(TriLevel),
            DROPSErrCL("FaceCL::GetNeighInTriang: Face or Tetra not in triangulation!"), DebugRefineEasyC);
    Assert( !IsOnBoundary(),
            DROPSErrCL("FaceCL::GetNeighInTriang: Face of tetra lies on boundary!"), DebugRefineEasyC);

    const Uint oppSide= _Neighbors[ tp->GetLevel()==GetLevel() ? 0 : 2 ]==tp;

    return _Neighbors[oppSide]->IsInTriang(TriLevel) ?
           _Neighbors[oppSide] : _Neighbors[oppSide+2];
}


Point3DCL GetBaryCenter(const FaceCL& f)
{
    const TetraCL* const tp= f.GetSomeTetra();
    const Uint face= f.GetFaceNumInTetra(tp);

    return ( tp->GetVertex( VertOfFace(face, 0))->GetCoord()
           + tp->GetVertex( VertOfFace(face, 1))->GetCoord()
           + tp->GetVertex( VertOfFace(face, 2))->GetCoord() )/3.0;
}

// T e t r a C L

Point3DCL GetBaryCenter(const TetraCL& t)
{
    return 0.25*( t.GetVertex(0)->GetCoord() + t.GetVertex(1)->GetCoord()
                + t.GetVertex(2)->GetCoord() + t.GetVertex(3)->GetCoord() );
}

Point3DCL GetBaryCenter(const TetraCL& t, Uint face)
{
    return ( t.GetVertex(VertOfFace(face, 0))->GetCoord()
            +t.GetVertex(VertOfFace(face, 1))->GetCoord()
            +t.GetVertex(VertOfFace(face, 2))->GetCoord() )/3.;
}

Point3DCL GetWorldCoord(const TetraCL& t, const SVectorCL<3>& c)
{
    return (1. -c[0] -c[1] -c[2])*t.GetVertex(0)->GetCoord()
          +c[0]*t.GetVertex(1)->GetCoord()
          +c[1]*t.GetVertex(2)->GetCoord()
          +c[2]*t.GetVertex(3)->GetCoord();
}

Point3DCL GetWorldCoord(const TetraCL& t, const SVectorCL<4>& c)
{
    return c[0]*t.GetVertex(0)->GetCoord()
          +c[1]*t.GetVertex(1)->GetCoord()
          +c[2]*t.GetVertex(2)->GetCoord()
          +c[3]*t.GetVertex(3)->GetCoord();
}

Point3DCL GetWorldCoord(const TetraCL& t, Uint face, const SVectorCL<2>& c)
{
    return (1. -c[0] -c[1])*t.GetVertex(VertOfFace(face, 0))->GetCoord()
          +c[0]*t.GetVertex(VertOfFace(face, 1))->GetCoord()
          +c[1]*t.GetVertex(VertOfFace(face, 2))->GetCoord();
}

SVectorCL<3> FaceToTetraCoord(__UNUSED__ const TetraCL& t, Uint f, SVectorCL<2> c)
{
    SVectorCL<3> ret(0.);
    switch(f)
    {
      case 0:
        ret[0]= 1 -c[0] -c[1]; ret[1]= c[0]; ret[2]= c[1]; break;
      case 1:
        ret[1]= c[0]; ret[2]= c[1]; break;
      case 2:
        ret[0]= c[0]; ret[2]= c[1]; break;
      case 3:
        ret[0]= c[0]; ret[1]= c[1]; break;
      default: throw DROPSErrCL("FaceToTetraCoord: illegal face-number.");
    }
    Assert( (GetWorldCoord(t,f,c)-GetWorldCoord(t, ret)).norm_sq() < 1.e-15, DROPSErrCL("FaceToTetraCoord: inconsistent mapping!"), DebugNumericC);
    return ret;
}


World2BaryCoordCL::World2BaryCoordCL (const TetraCL& t)
{
    SMatrixCL<4,4>& m= qr_.GetMatrix();
    for (Uint v= 0; v < 4; ++v) {
        const Point3DCL& p= t.GetVertex( v)->GetCoord();
        for (Uint i= 0; i < 3; ++i) m( i, v)= p[i];
    }
    for (Uint j= 0; j < 4; ++j) m( 3, j)= 1.;
    qr_.prepare_solve();
}

BaryCoordCL
World2BaryCoordCL::operator() (const Point3DCL& p) const
{
    BaryCoordCL r( MakeBaryCoord( p[0], p[1], p[2], 1.));
    qr_.Solve( r);
    return r;
}


double TetraCL::GetVolume () const
{
    Point3DCL v[3];
    const Point3DCL p0( _Vertices[0]->GetCoord() );
    for (Uint i=1; i<NumVertsC; ++i)
        v[i-1] = _Vertices[i]->GetCoord() - p0;
    return std::fabs(  v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1])
                     - v[0][1] * (v[1][0]*v[2][2] - v[1][2]*v[2][0])
                     + v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]) ) / 6.0;
}

double
TetraCL::GetNormal(Uint face, Point3DCL& normal, double& dir) const
/**
dir is set to 1.0 if the normal points out of the tetra,
else it is set to -1.0.
Normal has unit length, but the length of the cross-product is returned,
which is useful for integration over that face
If the triangulation is consistently numbered, both tetras on a face will
return the same normal in "normal" (, of course "dir" will be different).
*/
{
    const VertexCL* v[3];
    for (Uint i=0; i<3; ++i)
        v[i]= GetVertex( VertOfFace(face, i) );
    const VertexCL* const opvert= GetVertex( OppVert(face) );

    cross_product(normal, v[1]->GetCoord()-v[0]->GetCoord(),
                          v[2]->GetCoord()-v[0]->GetCoord());
    dir= inner_prod(opvert->GetCoord() - v[0]->GetCoord(), normal) < 0.0 ? 1. : -1.;
    const double absdet2D= normal.norm();
    normal/= absdet2D;
    return absdet2D;
}

double
TetraCL::GetOuterNormal(Uint face, Point3DCL& normal) const
/**
Returns the length of the cross-product;
"normal" is set to the unit outward normal of face "face"
*/
{
    double dir;
    const double absdet2D= GetNormal(face, normal, dir);
    normal*= dir;
    return absdet2D;
}

void TetraCL::BuildEdges (EdgeContT& edgecont)
{
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
        VertexCL* const v0= _Vertices[VertOfEdge(edge,0)];
        VertexCL* const v1= _Vertices[VertOfEdge(edge,1)];
        if ( !(_Edges[edge]= v0->FindEdge(v1) ) )
        {
            std::vector<BndPointCL> commonBndVerts;
            if ( v0->IsOnBoundary() && v1->IsOnBoundary() )
            {
                commonBndVerts.reserve(2);
                std::set_intersection( v0->GetBndVertBegin(), v0->GetBndVertEnd(),
                                       v1->GetBndVertBegin(), v1->GetBndVertEnd(),
                                       std::back_inserter(commonBndVerts), BndPointSegLessCL() );
            }
            switch( commonBndVerts.size() )
            {
              case 0:
                edgecont.push_back( EdgeCL(v0, v1, GetLevel()) );   break;
              case 1:
                edgecont.push_back( EdgeCL(v0, v1, GetLevel(), commonBndVerts[0].GetBndIdx() ) );   break;
              case 2:
                edgecont.push_back( EdgeCL(v0, v1, GetLevel(), commonBndVerts[0].GetBndIdx(), commonBndVerts[1].GetBndIdx()) );   break;
              default:
                v0->DebugInfo( std::cerr); std::cerr << std::endl;
                v1->DebugInfo( std::cerr); std::cerr << std::endl;
                throw DROPSErrCL("TetraCL::BuildEdges: Found edge on more than two BndSegs!");
            }
            _Edges[edge]= &edgecont.back();
            _Edges[edge]->RecycleMe();
        }
    }
}

BndIdxT GetCommonBndSeg(const VertexCL* v0, const VertexCL* v1, const VertexCL* v2)
{
    if (!(v0->IsOnBoundary() && v1->IsOnBoundary() && v2->IsOnBoundary() ))
        return NoBndC;
    std::list<BndPointCL> intersec01, intersec012;
    std::set_intersection( v0->GetBndVertBegin(), v0->GetBndVertEnd(),
                           v1->GetBndVertBegin(), v1->GetBndVertEnd(),
                           std::back_inserter(intersec01), BndPointSegLessCL() );
    std::set_intersection( v2->GetBndVertBegin(), v2->GetBndVertEnd(),
                           intersec01.begin(), intersec01.end(),
                           std::back_inserter(intersec012), BndPointSegLessCL() );
    if (intersec012.empty() )
        return NoBndC;
    Assert( intersec012.size()==1,
        DROPSErrCL("GetCommonBndSeg: found more than one BndSeg connected to three different boundary vertices"), ~0);
    return intersec012.begin()->GetBndIdx();
}

void TetraCL::BuildAndLinkFaces (FaceContT& facecont)  // used by XXXBuilderCL
{
    for (Uint face=0; face<NumFacesC; ++face)
    {
        VertexCL* const vp0= GetVertMidVert(VertOfFace(face, 0));
        VertexCL* const vp1= GetVertMidVert(VertOfFace(face, 1));
        VertexCL* const vp2= GetVertMidVert(VertOfFace(face, 2));

        if ( !(_Faces[face]= vp0->FindFace(vp1,vp2) ) )
        {
            facecont.push_back( FaceCL(GetLevel(), GetCommonBndSeg(vp0, vp1, vp2) ) );
            _Faces[face]= &facecont.back();
            _Faces[face]->RecycleMe(vp0, vp1, vp2);
        }
        _Faces[face]->LinkTetra(this);
    }
}

void
TetraCL::SetFace(Uint f, FaceCL* fp)
{
    _Faces[f]= fp;
}

void
TetraCL::SetEdge(Uint e, EdgeCL* ep)
{
    _Edges[e]= ep;
}

void
TetraCL::SetChild(Uint c, TetraCL* cp)
{
    if (!_Children) _Children= new SArrayCL<TetraCL*, MaxChildrenC>;
    (*_Children)[c]= cp;
}

// member functions for r e f i n e m e n t

void TetraCL::RecycleReusables()
/**
is called, if the refinement rule has changed.
It recycles and rescues simplices, that will be reused:
- Edges of children, that are not edges of the parent, if they are used by the new rule
- Faces of children, that are not faces of the parent,  -- " --
- Children                                           ,  -- " --
*/
{
    const RefRuleCL& myRule= GetRefData();
    const RefRuleCL& newRule= DROPS::GetRefRule(GetRefMark() & 63);

    // prepare a list of the faces that must be recycled
    std::list<byte> commonEdges;
    std::list<byte> commonFaces;
    std::list<byte> commonChildren;

    typedef std::list<byte>::iterator SetIterT;

    std::set_intersection( myRule.Edges,  myRule.Edges +  myRule.EdgeNum,
                          newRule.Edges, newRule.Edges + newRule.EdgeNum,
                          std::back_inserter(commonEdges) );
    std::set_intersection( myRule.Faces,  myRule.Faces +  myRule.FaceNum,
                          newRule.Faces, newRule.Faces + newRule.FaceNum,
                          std::back_inserter(commonFaces) );
    std::set_intersection( myRule.Children,  myRule.Children +  myRule.ChildNum,
                          newRule.Children, newRule.Children + newRule.ChildNum,
                          std::back_inserter(commonChildren) );

    for (Uint ch=0; ch<myRule.ChildNum; ++ch)
    {
        const ChildDataCL childdat= GetChildData(myRule.Children[ch]);
        TetraCL* const child= (*_Children)[ch];
        // recycle and rescue common edges
        for (Uint edge= 0; edge<NumEdgesC; ++edge)
        {
            if ( IsParentEdge(childdat.Edges[edge]) ) continue;
            SetIterT it= std::lower_bound( commonEdges.begin(), commonEdges.end(), childdat.Edges[edge]);
            if (it != commonEdges.end() && *it == childdat.Edges[edge])
            {
                child->_Vertices[VertOfEdge(edge,0)]->ClearRemoveMark();
                child->_Vertices[VertOfEdge(edge,1)]->ClearRemoveMark();
                child->_Edges[edge]->ClearRemoveMark();
                child->_Edges[edge]->RecycleMe();
                commonEdges.erase(it);  // because edge is now already recycled
            }
        }
        // recycle and rescue common faces
        for (Uint face=0; face<NumFacesC; ++face)
        {
            if ( IsParentFace(childdat.Faces[face]) ) continue;
            SetIterT it= std::lower_bound(commonFaces.begin(), commonFaces.end(), childdat.Faces[face]);
            if (it != commonFaces.end() && *it == childdat.Faces[face])
            {
                VertexCL* const vp0= child->_Vertices[VertOfFace(face, 0)];
                VertexCL* const vp1= child->_Vertices[VertOfFace(face, 1)];
                VertexCL* const vp2= child->_Vertices[VertOfFace(face, 2)];

                child->_Faces[face]->ClearRemoveMark();
                child->_Faces[face]->RecycleMe(vp0,vp1,vp2);
                commonFaces.erase(it);  // because face is now already recycled
            }
        }
    }
    // recycle, rescue and unlink children
    for (Uint ch=0; ch<myRule.ChildNum; ++ch)
    {
        if ( std::binary_search(commonChildren.begin(), commonChildren.end(), myRule.Children[ch]) )
        {
            TetraCL* const child= (*_Children)[ch];

            child->SetNoRefMark();
            child->UnlinkFromFaces();
            child->RecycleMe();
        }
    }
}


void TetraCL::ClearAllRemoveMarks()
/**
is called if MarkEqRule(). It safes all sub simplices from removement.
*/
{
    const RefRuleCL& myRule= GetRefData();

    for (Uint ch=0; ch<myRule.ChildNum; ++ch)
    {
        const ChildDataCL childdat= GetChildData(myRule.Children[ch]);
        TetraCL* const child= (*_Children)[ch];

        for (VertexPIterator vertPIt(child->_Vertices.begin()); vertPIt!=child->_Vertices.end(); ++vertPIt)
            (*vertPIt)->ClearRemoveMark();
        for (Uint edge= 0; edge<NumEdgesC; ++edge)
            child->_Edges[edge]->ClearRemoveMark();
        for (Uint face= 0; face<NumFacesC; ++face)
            child->_Faces[face]->ClearRemoveMark();
    }
}


void TetraCL::RestrictMark()
{
    if ( IsUnrefined() )
    {
        if ( IsMarkedForRegRef() && IsRegular() )
            CommitRegRefMark();
        if ( GetLevel() == 0 && IsMarkedForRemovement() )
            SetNoRefMark();
    }
    else
        if ( IsRegularlyRef() )
        {
            bool keepAnyChild= false;
            for (ChildPIterator childPIt=GetChildBegin(); childPIt!=GetChildEnd(); ++childPIt)
                if ( (*childPIt)->IsMarkedForRemovement() )
                    (*childPIt)->SetNoRefMark();
                else keepAnyChild= true;
            if (!keepAnyChild)
            {
                SetNoRefMark();
                UnCommitRegRefMark();
            }
        }
        else
        { // tetra is irregularly refined
            bool setregrefmark= false;
            for (ChildPIterator chp= GetChildBegin(), end= GetChildEnd(); chp!=end; ++chp)
            {
                if (!setregrefmark) {
                    if ( (*chp)->IsMarkedForRef() )
                        setregrefmark= true;
                    else
                    {
                        for (const_EdgePIterator ep= (*chp)->GetEdgesBegin(), edgeend= (*chp)->GetEdgesEnd(); !setregrefmark && ep!=edgeend; ++ep)
                            if ( (*ep)->IsMarkedForRef() && (*ep)->GetLevel()!=GetLevel() )
                            // parent edges are ignored
                            {
                                setregrefmark= true;
                                break;
                            }
                    }
                }
                (*chp)->SetNoRefMark();
            }
            if (setregrefmark)
            {
                SetRegRefMark();
                CommitRegRefMark();
            }
            else
                SetNoRefMark();
        }
}

void TetraCL::CollectEdges (const RefRuleCL& refrule,
                            VertContT& vertcont, EdgeContT& edgecont,
                            const BoundaryCL& Bnd)
/**
The edges for new refinement are stored in the static TetraCL::ePtrs array.
First look for them in the recycle bin (maybe they were created before),
if the edge cannot be found, create it.
*/
{
    const Uint nextLevel= GetLevel()+1;

    // Collect obvious edges
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
        VertexCL* const vp0= _Vertices[VertOfEdge(edge, 0)];
        VertexCL* const vp1= _Vertices[VertOfEdge(edge, 1)];
        EdgeCL*   const ep = _Edges[edge];

        if ( ep->IsMarkedForRef() )
            if ( ep->IsRefined() )
            {
                _ePtrs[SubEdge(edge, 0)]= vp0->FindEdge(ep->GetMidVertex());
                _ePtrs[SubEdge(edge, 1)]= ep->GetMidVertex()->FindEdge(vp1);
                Assert(_ePtrs[SubEdge(edge, 0)], DROPSErrCL("CollectEdges: SubEdge 0 not found."), DebugRefineEasyC);
                Assert(_ePtrs[SubEdge(edge, 1)], DROPSErrCL("CollectEdges: SubEdge 1 not found."), DebugRefineEasyC);
            }
            else
            {
                ep->BuildSubEdges(edgecont, vertcont, Bnd);

                EdgeContT::iterator sub= edgecont.end();
                _ePtrs[SubEdge(edge, 1)]= &*(--sub);
                sub->RecycleMe();
                _ePtrs[SubEdge(edge, 0)]= &*(--sub);
                sub->RecycleMe();
            }
        else
        {
            _ePtrs[edge]= ep;
        }
    }

    // Collect inner edges
    byte const* tmp= std::lower_bound(refrule.Edges, refrule.Edges+refrule.EdgeNum, static_cast<byte>(NumObviousEdgesC) );
        // pointer to first inner edge entry
    for (Uint i= tmp - refrule.Edges; i<refrule.EdgeNum; ++i)
    {
        const Uint      edge = refrule.Edges[i];
        VertexCL* const vp0  = GetVertMidVert(VertOfEdge(edge, 0));
        VertexCL* const vp1  = GetVertMidVert(VertOfEdge(edge, 1));

        if ( !(_ePtrs[edge]=vp0->FindEdge(vp1)) )
        {
            if ( IsDiagonal(edge) )
                edgecont.push_back( EdgeCL(vp0, vp1, nextLevel) );
            else // lies on parent face
                edgecont.push_back( EdgeCL(vp0, vp1, nextLevel, _Faces[ParFaceOfEdge(edge)]->GetBndIdx() ) );
            _ePtrs[edge] = &edgecont.back();
            _ePtrs[edge]->RecycleMe();
        }
    }
}

void TetraCL::CollectFaces (const RefRuleCL& refrule, FaceContT& facecont)
/**
The faces for new refinement are stored in the static TetraCL::fPtrs array.
First look for them in the recycle bin (maybe they were created before),
if the face cannot be found, create it and link boundary, if necessary.
*/
{
    const Uint nextLevel= GetLevel()+1;

    for (Uint i=0; i<refrule.FaceNum; ++i)
    {
        const Uint face= refrule.Faces[i];

        if (IsParentFace(face))
            _fPtrs[face]= _Faces[face];
        else
        {
                  VertexCL* const vp0= GetVertMidVert(VertOfFace(face, 0));
            const VertexCL* const vp1= GetVertMidVert(VertOfFace(face, 1));
            const VertexCL* const vp2= GetVertMidVert(VertOfFace(face, 2));
            if (!(_fPtrs[face]= vp0->FindFace(vp1, vp2) ) )
            {
                if ( IsSubFace(face) )
                    facecont.push_back( FaceCL(nextLevel, GetFace(ParentFace(face))->GetBndIdx() ) );
                else
                    facecont.push_back( FaceCL(nextLevel) );

                _fPtrs[face] = &facecont.back();
                _fPtrs[face]->RecycleMe(vp0, vp1, vp2);
            }
        }
    }
}

void TetraCL::CollectAndLinkChildren (const RefRuleCL& refrule, TetraContT& tcont)
/**
The child tetras for new refinement are stored in the _Children array.
First look for them in the recycle bin (maybe they are still left from the old rule),
if the child cannot be found, create it.
*/
{
    if ( !_Children ) _Children= new SArrayCL<TetraCL*, MaxChildrenC>;
    Uint ChildNum= refrule.ChildNum;
    for (Uint ch=0; ch < ChildNum; ++ch)
    {
        const ChildDataCL childdat= GetChildData(refrule.Children[ch]);
        VertexCL* const vp0= GetVertMidVert(childdat.Vertices[0]);
        VertexCL* const vp1= GetVertMidVert(childdat.Vertices[1]);
        VertexCL* const vp2= GetVertMidVert(childdat.Vertices[2]);
        VertexCL* const vp3= GetVertMidVert(childdat.Vertices[3]);

        if (!( (*_Children)[ch]= vp0->FindTetra(vp1, vp2, vp3) ))
        {
            tcont.push_back(TetraCL(vp0, vp1, vp2, vp3, this) );
            (*_Children)[ch] = &tcont.back();
        }
        (*_Children)[ch]->LinkEdges(childdat);
        (*_Children)[ch]->LinkFaces(childdat);
    }
}

/**********************************************************************************************************
*
*    D e b u g   f u n c t i o n s   to verify the validity of a multigrid
*
**********************************************************************************************************/
void RecycleBinCL::DebugInfo(std::ostream& os) const
{
    os << "RecycleBinCL:\nrecycled Edges: ";
    for (EdgeContT::const_iterator it= _Edges.begin(), theend= _Edges.end(); it!=theend;++it)
        (*it)->DebugInfo(os);
    os << "recycled Faces: ";
    for (FaceWrapperContT::const_iterator it= _Faces.begin(), theend= _Faces.end(); it!=theend;++it)
        it->face->DebugInfo(os);
    os << "recycled Tetras: ";
    for (TetraContT::const_iterator it= _Tetras.begin(), theend= _Tetras.end(); it!=theend;++it)
        (*it)->DebugInfo(os);
    os << std::endl;
}

/// \todo: Checke die BoundarySegments!!

bool VertexCL::IsSane(std::ostream& os, const BoundaryCL& Bnd) const
{
    bool sane= true;

    // Check, if all boundary descriptions map to the same coordinates
    if (_BndVerts) {
        for (std::vector<BndPointCL>::const_iterator bIt(_BndVerts->begin());
             bIt != _BndVerts->end(); ++bIt) {
            if (dynamic_cast<const MeshBoundaryCL*>( Bnd.GetBndSeg(bIt->GetBndIdx())))
                continue; // We ignore MeshBoundaryCL as it does not have 2D-coordinates...
            if ((Bnd.GetBndSeg(bIt->GetBndIdx())->Map(bIt->GetCoord2D()) - GetCoord()).norm() > DoubleEpsC)
            {
                sane= false;
                os << "BndSegCL description " << bIt->GetBndIdx()
                   << " does not match the coordinates. ";
                os << "Mapping gives: " << Bnd.GetBndSeg(bIt->GetBndIdx())->Map(bIt->GetCoord2D()) << ' ';
            }
        }
    }
    // Check, that the refinement algorithm did not miss any RecycleBins
    if ( HasRecycleBin() )
    {
        sane= false;
        os << "Clear your RecycleBin!";
    }
    if (!sane) os << std::endl;
    return sane;
}


void VertexCL::DebugInfo(std::ostream& os) const
{
    os << "VertexCL:  Id " << _Id.GetIdent()
       << " Level " << GetLevel() << ' ';
    os << "RemoveMark " << _RemoveMark << std::endl << "  ";
    os << "Coord ";
    for (Uint i=0; i<3; ++i) os << _Coord[i] << " ";
    os << std::endl << "  ";
    if (_BndVerts)
    {
        os << "BndVerts ";
        for (std::vector<BndPointCL>::const_iterator biter(_BndVerts->begin()); biter!=_BndVerts->end(); ++biter)
            os << biter->GetBndIdx()<< "  "
               << biter->GetCoord2D() << "    ";
    }
    else os << "no BndVerts found";
    os << std::endl;
}

bool EdgeCL::IsSane(std::ostream& os) const
{
    bool sane= true;
    if (IsOnBoundary() && !(GetVertex(0)->IsOnBoundary() && GetVertex(1)->IsOnBoundary() ) )
    {
        sane= false;
        os << "One of the vertices is on no boundary even though the edge is. ";
    }
    if (sane)
        for (const BndIdxT *it= GetBndIdxBegin(), *end= GetBndIdxEnd(); it!=end; ++it)
        {
            if (!is_in_if( GetVertex(0)->GetBndVertBegin(), GetVertex(0)->GetBndVertEnd(), BndPointSegEqCL(*it) ) )
            {
               sane= false;
               os << "BndIdx " << *it << " is not among the boundaries of vertex 0. ";
            }
            if (!is_in_if( GetVertex(1)->GetBndVertBegin(), GetVertex(1)->GetBndVertEnd(), BndPointSegEqCL(*it) ) )
            {
               sane= false;
               os << "BndIdx " << *it << " is not among the boundaries of vertex 1. ";
            }
        }
    if (!sane) os << std::endl;
    return sane;
}

void EdgeCL::DebugInfo (std::ostream& os) const
{
    os << "EdgeCL: Level " << GetLevel()
       << "  Vertices " <<GetVertex(0)->GetId().GetIdent()
       << " " << GetVertex(1)->GetId().GetIdent() << "    "
       << " MarkForRef " << _MFR // << "(" << _localMFR << ")"
       << " RemoveMark " << _RemoveMark;
    os << std::endl << "  ";
    os << "Midvertex ";
    if ( IsRefined() ) os << GetMidVertex()->GetId().GetIdent();
    else               os << "not found ";
    if ( IsOnBoundary() )
    {
        os << " BndIndices ";
        for (const BndIdxT *it= GetBndIdxBegin(), *end= GetBndIdxEnd(); it!=end; ++it)
            os << *it << ' ';
    }
    os << '\n';
}

bool FaceCL::IsSane(std::ostream& os) const
{
    bool sane= true;
    // check, that both sides point to a tetra or a BndSeg in my level

    if (_Neighbors[0]==0 || (_Neighbors[1]!=0) == IsOnBoundary() )
    {
        sane= false;
        os << "A tetra/boundary is missing/superfluous. ";
    }
    // check, that linked neighbor tetras have me as face
    const TetraCL* t;
    for(Uint i=0; i<4; ++i)
    {
        if ( (t= _Neighbors[i]) != 0 )
        {
            if ( !is_in( t->GetFacesBegin(), t->GetFacesEnd(), this))
            {
                sane= false;
                os << "Found linked tetra, that lacks me as face:\n";
                t->DebugInfo(os);
            }
            if (t->GetLevel() != GetLevel() + i/2)
            {
                sane= false;
                os << "Found linked tetra on wrong level:\n";
                t->DebugInfo(os);
            }
        }
    }
    if ( IsOnBoundary() )
    {
        if (!(GetVertex(0)->IsOnBoundary() && GetVertex(1)->IsOnBoundary() && GetVertex(2)->IsOnBoundary() ) )
        {
            sane= false;
            os << "One of the vertices is on no boundary even though the face is. ";
        }
        else if (GetBndIdx() != GetCommonBndSeg(GetVertex(0), GetVertex(1), GetVertex(2)) )
        {
            sane= false;
            os << "BndIdx " << GetBndIdx() << " is not among the common boundaries of the face. ";
        }
    }
    if (!sane) os << std::endl;
    return sane;
}

void FaceCL::DebugInfo(std::ostream& os) const
{
    os << "FaceCL: Level " << GetLevel()
       << " RemoveMark " << _RemoveMark;
    if (IsOnBoundary() ) os << " is on boundary";
    os << "\nNeighbors in this Level:";
    for (Uint i=0; i<2; ++i)
    {
        if (i==1 && IsOnBoundary() )
            os << " BndIdx "<< GetBndIdx();
        else
            os << " Tetra " << _Neighbors[i]->GetId().GetIdent();
    }
    if (IsOnNextLevel() )
    {
        os << "\nNeighbors in next Level:";
        for (Uint i=2; i<4; ++i)
        {
            if (i==3 && IsOnBoundary() )
                os << " BndIdx "<< GetBndIdx();
            else if (_Neighbors[i])
                os << " Tetra " << _Neighbors[i]->GetId().GetIdent();
            else
                os << " not there ";
        }
    }
    os << std::endl;
}

bool TetraCL::IsSane(std::ostream& os) const
{
    bool sane= true;

    // Check, if the volume of my children adds up to my volume
    if (_Children)
    {
        double vol = 0.0;
        for (const_ChildPIterator tCPIt(GetChildBegin()); tCPIt!=GetChildEnd(); ++tCPIt)
            vol += (*tCPIt)->GetVolume();
        if ( std::fabs(GetVolume() - vol) > DoubleEpsC )
        {
            sane= false;
            os << "Volume of children does not sum up to parent's volume.";
        }
    }
    // Check, if the neighbor-connections are right:
    // If it is another tetra, check, if the numbering
    // across the face is consistent
    // If it is a boundary-segment, check,
    // if the three vertices belong to this boundary-segment
    for (Uint face=0; face<NumFacesC; ++face)
    {
        if ( IsNeighbor(face) )
        {
            const TetraCL* const np = GetNeighbor(face);
            Uint nface= _Faces[face]->GetFaceNumInTetra(np);
            std::vector<Uint> vertnum0;
            std::vector<Uint> vertnum1;
            vertnum0.reserve(3);
            vertnum1.reserve(3);
            for (Uint i=0; i<3; ++i)
            {
              vertnum0.push_back(_Vertices[VertOfFace(face, i)]->GetId().GetIdent());
              vertnum1.push_back(np->_Vertices[VertOfFace(nface,i)]->GetId().GetIdent());
            }
            // Since we assume that the initial triangulation is numbered consistently
            // there is no need to sort the vertices in the face, as they should be ordered
            // in exactly the same way in both tetras;
//                    std::sort(vertnum0.begin(), vertnum0.end());
//                    std::sort(vertnum1.begin(), vertnum1.end());
            if (vertnum0 != vertnum1)
            {
                sane= false;
                os << "Neighborhood across face " << face << " is screwed or the ordering induced\n"
                   << "by the tetras is not the same, which is a BadThing (TM). ";
            }
        }
    }
    // Check, whether the ordering of the vertices in each edge is induced by the ordering of
    // the vertices in the tetra
    for (Uint edge=0; edge <NumEdgesC; ++edge)
    {
        if (   _Edges[edge]->GetVertex(0) != _Vertices[VertOfEdge(edge, 0)]
            || _Edges[edge]->GetVertex(1) != _Vertices[VertOfEdge(edge, 1)])
        {
            sane= false;
            os << "Edge " << edge << " and this tetra are not numbered consistently. ";
        }
    }
    // Check, whether the vertices of opposing edges contain all four vertices of the tetra
    std::vector<const VertexCL*> tverts;
    tverts.reserve(4);
    tverts.push_back(_Vertices[0]); tverts.push_back(_Vertices[1]);
    tverts.push_back(_Vertices[2]); tverts.push_back(_Vertices[3]);
    std::sort(tverts.begin(), tverts.end());
    for (Uint edge=0; edge<NumEdgesC; ++edge)
    {
        std::vector<const VertexCL*> everts;
        everts.reserve(4);
        everts.push_back(_Edges[edge]->GetVertex(0)); everts.push_back(_Edges[edge]->GetVertex(1));
        everts.push_back(_Edges[OppEdge(edge)]->GetVertex(0));
        everts.push_back(_Edges[OppEdge(edge)]->GetVertex(1));
        std::sort(everts.begin(), everts.end());
        if (tverts!=everts)
        {
            sane= false;
            os << "Edge " << edge << " and its opposite " << OppEdge(edge)
               << " have wrong vertices with respect to tetra " << GetId().GetIdent() << ". ";
        }
    }
    if (!sane) os << std::endl;
    return sane;
}


void TetraCL::DebugInfo (std::ostream& os) const
{
    os << "TetraCL:  Id " << _Id.GetIdent()
       << " Level " << GetLevel() << "\t"
       << "RefRule " << GetRefRule() << " RefMark " << GetRefMark() << std::endl;

    os << "  Vertices: ";
    for (Uint i=0; i<NumVertsC; ++i)
        os << _Vertices[i]->GetId().GetIdent() << "    ";
    os << std::endl;
    os << "  Edges: ";
    for (Uint i=0; i<NumEdgesC; ++i)
        os << _Edges[i]->GetVertex(0)->GetId().GetIdent()
           << " " << _Edges[i]->GetVertex(1)->GetId().GetIdent() << "    ";
    os << std::endl;
    if (_Parent) os << "  Parent: " << _Parent->GetId().GetIdent() << "  ";
    else os << "  no parent found  ";
    if (_Children)
    {
        os << "Children: ";
        for (Uint i=0; i<MaxChildrenC; ++i)
            if ((*_Children)[i]) os << (*_Children)[i]->GetId().GetIdent() << "    ";
            else os << "not there    ";
    }
    else
        os << "no children found";
    os << std::endl;
    os << "  Neighbors: ";
    for (Uint i=0; i<NumFacesC; ++i)
    {
        if ( IsNeighbor(i) )
            os << " TetraCL " << GetNeighbor(i)->GetId().GetIdent()
               << "    ";
        else if ( IsBndSeg(i) )
            os << " BndSegCL " << GetBndIdx(i)
               << "    ";
        else
            os << " not there.    ";
    }
    os << std::endl;
}


BoundaryCL::~BoundaryCL()
{
    for (SegPtrCont::iterator It=Bnd_.begin(); It!=Bnd_.end(); ++It)
        delete *It;
    delete BndType_;
}

void BoundaryCL::SetPeriodicBnd( const BndTypeCont& type, match_fun match) const
{
    if (type.size()!=GetNumBndSeg())
        throw DROPSErrCL("BoundaryCL::SetPeriodicBnd: inconsistent vector size!");
    BndType_= new BndTypeCont(type);
    match_= match;
}

BoundaryCL::BndType PeriodicEdgesCL::GetBndType( const EdgeCL& e) const
{
    BoundaryCL::BndType type= BoundaryCL::OtherBnd;
    for (const BndIdxT *bndIt= e.GetBndIdxBegin(), *end= e.GetBndIdxEnd(); bndIt!=end; ++bndIt)
        type= std::max( type, mg_.GetBnd().GetBndType(*bndIt));
    return type;
}   

void PeriodicEdgesCL::Accumulate()
{
    // initialize MFR counters on all Per1 edges
    for (iterator It( list_.begin()), End(list_.end()); It!=End; ++It)
        It->first->_MFR= It->first->_localMFR;
    // compute sum in Per1 MFR counters
    for (iterator It( list_.begin()), End(list_.end()); It!=End; ++It)
        It->first->_MFR+= It->second->_localMFR;
    // copy Per1 MFR counter to Per2 MFR counter 
    for (iterator It( list_.begin()), End(list_.end()); It!=End; ++It)
        It->second->_MFR= It->first->_MFR;    
}

void PeriodicEdgesCL::Recompute( EdgeIterator begin, EdgeIterator end)
{
    typedef std::list<EdgeCL*> psetT;
    psetT s1, s2;
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (EdgeIterator it= begin; it!=end; ++it)
        if (it->IsOnBoundary())
        {
            BoundaryCL::BndType type= GetBndType( *it);
            if (type==BoundaryCL::Per1Bnd)
                s1.push_back( &*it);
            else if (type==BoundaryCL::Per2Bnd)
                s2.push_back( &*it);
        }
    // now we have s1.size() <= s2.size()
    // match objects in s1 and s2
    const BoundaryCL& bnd= mg_.GetBnd();
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; )
            if (bnd.Matching( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // store pair in list_
                list_.push_back( IdentifiedEdgesT( *it1, *it2));
                // remove it2 from s2
                s2.erase( it2++);
            }
            else it2++;
    }
    if (!s2.empty())
        throw DROPSErrCL( "PeriodicEdgesCL::Recompute: Periodic boundaries do not match!");
}

void PeriodicEdgesCL::DebugInfo( std::ostream& os)
{
    int num= 0;
    for (PerEdgeContT::iterator it= list_.begin(), end=  list_.end(); it!=end; ++it, ++num)
    {
    	it->first->DebugInfo( os);
    	os << "\t\t<-- " << num << " -->\n";
    	it->second->DebugInfo( os);
    	os << "===================================================================\n";
    }
    os << num << " identified edges found.\n\n";	
}


void PeriodicEdgesCL::Shrink()
{
    list_.clear();
}

void PeriodicEdgesCL::AccumulateMFR( int lvl)
{
    if (!mg_.GetBnd().HasPeriodicBnd()) return;
    Shrink();
    for (int i=0; i<=lvl; ++i)
        Recompute( mg_.GetEdgesBegin(i), mg_.GetEdgesEnd(i));
//std::cerr << " \n>>> After Recompute:\n"; DebugInfo( std::cerr);
    Accumulate();
//std::cerr << " \n>>> After Accumulate:\n"; DebugInfo( std::cerr);
    Shrink();
}

MultiGridCL::MultiGridCL (const MGBuilderCL& Builder)
    : _TriangVertex( *this), _TriangEdge( *this), _TriangFace( *this), _TriangTetra( *this)
{
    Builder.build(this);
    FinalizeModify();
}


void MultiGridCL::CloseGrid(Uint Level)
{
    Comment("Closing grid " << Level << "." << std::endl, DebugRefineEasyC);

    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; ++tIt)
    {
        Comment("Now closing tetra " << tIt->GetId().GetIdent() << std::endl, DebugRefineHardC);
        if ( tIt->IsRegular() && !tIt->IsMarkedForRegRef() )
            tIt->Close();
    }

    Comment("Closing grid " << Level << " done." << std::endl, DebugRefineEasyC);
}


void MultiGridCL::UnrefineGrid (Uint Level)
{
    Comment("Unrefining grid " << Level << "." << std::endl, DebugRefineEasyC);

    const Uint nextLevel(Level+1);

    std::for_each(_Vertices[nextLevel].begin(), _Vertices[nextLevel].end(), std::mem_fun_ref(&VertexCL::SetRemoveMark));
    std::for_each(_Edges[nextLevel].begin(),    _Edges[nextLevel].end(),    std::mem_fun_ref(&EdgeCL::SetRemoveMark));
    std::for_each(_Faces[nextLevel].begin(),    _Faces[nextLevel].end(),    std::mem_fun_ref(&FaceCL::SetRemoveMark));

    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; ++tIt)
    {
        Comment("inspecting children of tetra " << tIt->GetId().GetIdent() << "." << std::endl, DebugRefineHardC);

        if ( !tIt->IsUnrefined() ) {
            if ( tIt->IsMarkEqRule() )
                tIt->ClearAllRemoveMarks();
            else
            {
                std::for_each(tIt->GetChildBegin(), tIt->GetChildEnd(), std::mem_fun(&TetraCL::SetRemoveMark));
                if ( !tIt->IsMarkedForNoRef() ) tIt->RecycleReusables();
            }
        }
    }

    Comment("Now physically unlinking and removing superfluous tetras." << std::endl, DebugRefineEasyC);
    for (TetraIterator it= _Tetras[nextLevel].begin(), end= _Tetras[nextLevel].end(); it!=end; )
        if ( it->IsMarkedForRemovement() )
        {
            it->UnlinkFromFaces();
            _Tetras[nextLevel].erase(it++);
        }
        else
            ++it;

    Comment("Now adapting midvertex pointers on level " << Level << ". " << std::endl, DebugRefineEasyC);
    for (EdgeIterator eIt(_Edges[Level].begin()), eEnd(_Edges[Level].end()); eIt!=eEnd; ++eIt)
        if ( (eIt->IsRefined() && eIt->GetMidVertex()->IsMarkedForRemovement()) )
            eIt->RemoveMidVertex();
    Comment("Now removing superfluous faces." << std::endl, DebugRefineEasyC);
    _Faces[nextLevel].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    Comment("Now removing superfluous edges." << std::endl, DebugRefineEasyC);
    _Edges[nextLevel].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
    Comment("Now physically removing superfluous vertices." << std::endl, DebugRefineEasyC);
    _Vertices[nextLevel].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
    Comment("Unrefining grid " << Level << " done." << std::endl, DebugRefineEasyC);
}


void MultiGridCL::RefineGrid (Uint Level)
{
    Comment("Refining grid " << Level << std::endl, DebugRefineEasyC);

    const Uint nextLevel(Level+1);
    if ( Level==GetLastLevel() ) AppendLevel();

    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; ++tIt)
    {
        if ( tIt->IsMarkEqRule() ) continue;

        tIt->SetRefRule( tIt->GetRefMark() );
        if ( tIt->IsMarkedForNoRef() )
        {
            Comment("refining " << tIt->GetId().GetIdent() << " with rule 0." << std::endl, DebugRefineHardC);
            if ( tIt->_Children )
                { delete tIt->_Children; tIt->_Children=0; }
        }
        else
        {
            const RefRuleCL& refrule( tIt->GetRefData() );

            Comment("refining " << tIt->GetId().GetIdent() << " with rule " << tIt->GetRefRule() << "." << std::endl, DebugRefineHardC);
            tIt->CollectEdges           (refrule, _Vertices[nextLevel], _Edges[nextLevel], _Bnd);
            tIt->CollectFaces           (refrule, _Faces[nextLevel]);
            tIt->CollectAndLinkChildren (refrule, _Tetras[nextLevel]);
        }
    }
    for (Uint lvl= 0; lvl <= nextLevel; ++lvl)
        std::for_each( _Vertices[lvl].begin(), _Vertices[lvl].end(),
            std::mem_fun_ref( &VertexCL::DestroyRecycleBin));

    Comment("Refinement of grid " << Level << " done." << std::endl, DebugRefineEasyC);
    if (DROPSDebugC & DebugRefineHardC) if ( !IsSane(cdebug, Level) ) cdebug << std::endl;
}


void MultiGridCL::Refine()
{
    const int tmpLastLevel( GetLastLevel() );

    PeriodicEdgesCL perEdges( *this);
    ClearTriangCache();
    PrepareModify();
    for (int Level=tmpLastLevel; Level>=0; --Level)
    {
        RestrictMarks(Level);
        perEdges.AccumulateMFR( Level);
        CloseGrid(Level);
    }

    for (int Level=0; Level<=tmpLastLevel; ++Level)
    {
        if ( _Tetras[Level].empty() ) continue;
        if (Level) CloseGrid(Level);
        if ( Level != tmpLastLevel ) UnrefineGrid(Level);
        RefineGrid(Level);
    }

    while ( _Tetras[GetLastLevel()].empty() ) RemoveLastLevel();
    FinalizeModify();

    std::for_each( GetAllVertexBegin(), GetAllVertexEnd(),
        std::mem_fun_ref( &VertexCL::DestroyRecycleBin));
}


void MultiGridCL::Scale( double s)
{
    for (VertexIterator it= GetAllVertexBegin(), end= GetAllVertexEnd();
        it!=end; ++it)
        it->_Coord*= s;
}


class VertPtrLessCL : public std::binary_function<const VertexCL*, const VertexCL* , bool>
{
  public:
    bool operator() (const VertexCL* v0, const VertexCL* v1)
        { return v0->GetId() < v1->GetId(); }
};


void MultiGridCL::MakeConsistentNumbering()
// Applicable only before the first call to Refine()
// Rearranges the Vertexorder in Tetras and Edges, so that it is the one induced by
// the global vertex-numbering in level 0
{
    // correct vertex-order in the edges
    std::for_each (GetEdgesBegin(0), GetEdgesEnd(0), std::mem_fun_ref(&EdgeCL::SortVertices));

    for (TetraIterator sit= GetTetrasBegin(0), theend= GetTetrasEnd(0); sit!=theend; ++sit)
    {
        VertexCL* vp[NumVertsC];
        std::copy(sit->_Vertices.begin(), sit->_Vertices.end(), vp+0);
        // correct vertex-order in tetras
        std::sort( sit->_Vertices.begin(), sit->_Vertices.end(), VertPtrLessCL() );

        // sort edge-pointers according to new vertex-order
        EdgeCL* ep[NumEdgesC];
        std::copy(sit->_Edges.begin(), sit->_Edges.end(), ep+0);
        for (Uint edge=0; edge<NumEdgesC; ++edge)
        {
            const Uint v0= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), ep[edge]->GetVertex(0)) );
            const Uint v1= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), ep[edge]->GetVertex(1)) );
            sit->_Edges[EdgeByVert(v0, v1)]= ep[edge];
        }

        // sort face-pointers according to new vertex-order
        FaceCL* fp[NumFacesC];
        std::copy(sit->_Faces.begin(), sit->_Faces.end(), fp);
        for (Uint face=0; face<NumFacesC; ++face)
        {
            const Uint v0= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), vp[VertOfFace(face, 0)]) );
            const Uint v1= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), vp[VertOfFace(face, 1)]) );
            const Uint v2= std::distance( sit->GetVertBegin(),
                               std::find(sit->GetVertBegin(), sit->GetVertEnd(), vp[VertOfFace(face, 2)]) );
            sit->_Faces[FaceByVert(v0, v1, v2)]= fp[face];
        }

    }
}


class EdgeByVertLessCL : public std::binary_function<const EdgeCL*, const EdgeCL* , bool>
{
  public:
    bool operator() (const EdgeCL* e0, const EdgeCL* e1)
        { return    e0->GetVertex(1)->GetId() == e1->GetVertex(1)->GetId()
                 ?  e0->GetVertex(0)->GetId() <  e1->GetVertex(0)->GetId()
                 :  e0->GetVertex(1)->GetId() <  e1->GetVertex(1)->GetId(); }
};


class EdgeEqualCL : public std::binary_function<const EdgeCL*, const EdgeCL*, bool>
{
  public:
    bool operator() (const EdgeCL* e0, const EdgeCL* e1)
        { return    e0->GetVertex(1)->GetId() == e1->GetVertex(1)->GetId()
                 && e0->GetVertex(0)->GetId() == e1->GetVertex(0)->GetId(); }
};


bool MultiGridCL::IsSane (std::ostream& os, int Level) const
{
    bool sane=true;

    if (Level==-1)
    {
        for (int lvl=0; lvl<=static_cast<int>(GetLastLevel()); ++lvl)
            if ( !IsSane(os, lvl) )
                sane = false;
    }
    else
    {
        // Check Vertices
        for (const_VertexIterator vIt( GetVerticesBegin( Level));
             vIt != GetVerticesEnd( Level); ++vIt)
        {
            if ( int(vIt->GetLevel())!=Level )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                vIt->DebugInfo(os);
            }
            if ( !vIt->IsSane(os, _Bnd) )
            {
                sane=false;
                vIt->DebugInfo(os);
            }
        }
        // Check Edges
        for (const_EdgeIterator eIt( GetEdgesBegin( Level));
             eIt!=GetEdgesEnd( Level); ++eIt)
        {
            if ( int(eIt->GetLevel())!=Level )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                eIt->DebugInfo(os);
            }

            if ( !eIt->IsSane(os) )
            {
                sane=false;
                eIt->DebugInfo(os);
            }
        }
        // An edge connecting two vertices should be unique in its level
        // This is memory-expensive!
        std::list<const EdgeCL*> elist;
        ref_to_ptr<const EdgeCL> conv;
        std::transform( GetEdgesBegin( Level), GetEdgesEnd( Level),
            std::back_inserter( elist), conv);
        elist.sort( EdgeByVertLessCL());
        if (std::adjacent_find( elist.begin(), elist.end(), EdgeEqualCL()) != elist.end() )
        {
            sane = false;
            os << "Found an edge more than once in level " << Level << ".\n";
            (*std::adjacent_find( elist.begin(), elist.end(), EdgeEqualCL()))->DebugInfo( os);
        }
        // Check Faces
        for (const_FaceIterator It( GetFacesBegin( Level));
             It!=GetFacesEnd( Level); ++It)
        {
            if ( It->GetLevel() != static_cast<Uint>(Level) )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                It->DebugInfo(os);
            }
            if ( !It->IsSane(os) )
            {
                sane = false;
                It->DebugInfo(os);
            }
        }
        // Check Tetras
        for (const_TetraIterator tIt( GetTetrasBegin( Level));
             tIt!=GetTetrasEnd( Level); ++tIt)
        {
            if ( tIt->GetLevel() != static_cast<Uint>(Level) )
            {
                sane=false;
                os <<"Wrong Level (should be "<<Level<<") for\n";
                tIt->DebugInfo(os);
            }
            if ( !tIt->IsSane(os) )
            {
                sane = false;
                tIt->DebugInfo(os);
            }
        }
    }
    return sane;
}


void MultiGridCL::SizeInfo(std::ostream& os)
{
    os << GetVertices().size() << " Verts, "
       << GetEdges().size()    << " Edges, "
       << GetFaces().size()    << " Faces, "
       << GetTetras().size()   << " Tetras"
       << std::endl;
}

void MultiGridCL::ElemInfo(std::ostream& os, int Level)
{
    double hmax= -1, hmin= 1e99,
           rmax= -1, rmin= 1e99;
    DROPS_FOR_TRIANG_TETRA( (*this), Level, It) {
        double loc_max= -1, loc_min= 1e99;
        for (Uint i=0; i<3; ++i)
        {
            Point3DCL pi= It->GetVertex(i)->GetCoord();
            for (Uint j=i+1; j<4; ++j)
            {
                const double h= (It->GetVertex(j)->GetCoord() - pi).norm();
                if (h < loc_min) loc_min= h;
                if (h > loc_max) loc_max= h;
            }
        }
        if (loc_min < hmin) hmin= loc_min;
        if (loc_max > hmax) hmax= loc_max;
        const double ratio= loc_max/loc_min;
        if (ratio < rmin) rmin= ratio;
        if (ratio > rmax) rmax= ratio;
    }
    os << hmin << " <= h <= " << hmax << '\t'
       << rmin << " <= h_max/h_min <= " << rmax << std::endl;
}


void
TriangFillCL<VertexCL>::fill (MultiGridCL& mg, TriangCL<VertexCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::VertexIterator it= mg.GetAllVertexBegin( lvl),
         theend= mg.GetAllVertexEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)) c.push_back( &*it);
    TriangCL<VertexCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<EdgeCL>::fill (MultiGridCL& mg, TriangCL<EdgeCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::EdgeIterator it= mg.GetAllEdgeBegin( lvl),
         theend= mg.GetAllEdgeEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)) c.push_back( &*it);
    TriangCL<EdgeCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<FaceCL>::fill (MultiGridCL& mg, TriangCL<FaceCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::FaceIterator it= mg.GetAllFaceBegin( lvl),
         theend= mg.GetAllFaceEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)) c.push_back( &*it);
    TriangCL<FaceCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<TetraCL>::fill (MultiGridCL& mg, TriangCL<TetraCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::TetraIterator it= mg.GetAllTetraBegin( lvl),
         theend= mg.GetAllTetraEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)) c.push_back( &*it);
    TriangCL<TetraCL>::LevelCont tmp= c;
    c.swap( tmp);
}


void circumcircle(const TetraCL& t, Point3DCL& c, double& r)
{
    SMatrixCL<3,3> M;
    const Point3DCL p0= t.GetVertex(0)->GetCoord();
    const double p0_sq= p0.norm_sq();

    for (int i=0; i<3; ++i)
    {
        const Point3DCL p= t.GetVertex(i+1)->GetCoord();
        const Point3DCL p2= p - p0;
        for (Uint j=0; j<3; ++j)
            M(i, j)= p2[j];
        c[i]= 0.5*(p.norm_sq()-p0_sq);
    }
    gauss_pivot(M, c);

    if (DebugRefineEasyC
        && (std::fabs( (c-p0).norm() - (c-t.GetVertex(1)->GetCoord()).norm()) > DoubleEpsC
            || std::fabs( (c-p0).norm() - (c-t.GetVertex(2)->GetCoord()).norm()) > DoubleEpsC
            || std::fabs( (c-p0).norm() - (c-t.GetVertex(3)->GetCoord()).norm()) > DoubleEpsC) )
    {
        std::cerr << "circumcircle: Didn't find circumcenter. c: " << c << '\n';
        t.DebugInfo(std::cerr);
        throw DROPSErrCL();
    }
    r= (c-p0).norm();
}

void circumcircle(const TetraCL& t, Uint face, Point3DCL& c, double& r)
{
    const Point3DCL& v0= t.GetVertex( VertOfFace(face, 0) )->GetCoord();
    const Point3DCL& p1= t.GetVertex( VertOfFace(face, 1) )->GetCoord() - v0;
    const Point3DCL& p2= t.GetVertex( VertOfFace(face, 2) )->GetCoord() - v0;
    double m[2];
    const double p1_sq= p1.norm_sq();
    const double p2_sq= p2.norm_sq();
    const double p1p2= inner_prod(p1,p2);
    const double det= p1_sq*p2_sq - p1p2*p1p2;

    m[0]= 0.5*p2_sq*(p1_sq - p1p2)/det;
    m[1]= 0.5*p1_sq*(p2_sq - p1p2)/det;
    c= v0 + m[0]*p1 + m[1]*p2;

    if (DebugRefineEasyC
        && (   std::fabs( (c-v0).norm() - (c-(p1+v0)).norm()) > DoubleEpsC
            || std::fabs( (c-v0).norm() - (c-(p2+v0)).norm()) > DoubleEpsC) )
    {
        std::cerr << "circumcircle: Didn't find circumcenter. c: " << c << " face: " << face << '\n';
        t.DebugInfo(std::cerr);
        throw DROPSErrCL();
    }
    r= (c-v0).norm();
}


void MarkAll (DROPS::MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/, It)
        It->SetRegRefMark();
}


void UnMarkAll (DROPS::MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/, It)
        It->SetRemoveMark();
}


} // end of namespace DROPS
