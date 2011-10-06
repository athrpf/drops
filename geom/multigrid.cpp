/// \file multigrid.cpp
/// \brief classes that constitute the multigrid
/// \author LNM RWTH Aachen: Sven Gross, Eva Loch, Joerg Peters, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

/// Remarks: We should use the const-qualifier to make it difficult to
///          accidentally change the multigrid structure from anywhere
///          outside of the multigrid algorithms.
///          Thus the pointer to user data structures should probably be
///          a pointer to mutable.

#ifdef _PAR
#include "parallel/parmultigrid.h"
#include "parallel/parallel.h"
#endif

#include "geom/multigrid.h"
#include "num/solver.h"
#include <set>

namespace DROPS
{

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
#ifdef _PAR
    for (size_t i=0; i<type.size(); ++i){
        if (type[i]!=OtherBnd){
            throw DROPSErrCL("No periodic boundary conditions implemented in the parallel version, yet");
        }
    }
#endif
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
//std::cout << " \n>>> After Recompute:\n"; DebugInfo( std::cout);
    Accumulate();
//std::cout << " \n>>> After Accumulate:\n"; DebugInfo( std::cout);
    Shrink();
}


MultiGridCL::MultiGridCL (const MGBuilderCL& Builder)
    : _TriangVertex( *this), _TriangEdge( *this), _TriangFace( *this), _TriangTetra( *this), _version(0)
{
    Builder.build(this);
    FinalizeModify();
#ifdef _PAR
    ParMultiGridCL::AttachTo( *this);
    ParMultiGridCL::MarkSimplicesForUnknowns();
#endif
}

void MultiGridCL::ClearTriangCache ()
{
    _TriangVertex.clear();
    _TriangEdge.clear();
    _TriangFace.clear();
    _TriangTetra.clear();

    for (std::map<int, ColorClassesCL*>::iterator it= _colors.begin(), end= _colors.end(); it != end; ++it)
        delete it->second;
    _colors.clear();
}

void MultiGridCL::CloseGrid(Uint Level)
{
    Comment("Closing grid " << Level << "." << std::endl, DebugRefineEasyC);

    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; ++tIt)
    {
#ifdef _PAR
        AllComment("Now closing tetra " << tIt->GetGID() << std::endl, DebugRefineHardC);
#else
        Comment("Now closing tetra " << tIt->GetId().GetIdent() << std::endl, DebugRefineHardC);
#endif
        if ( tIt->IsRegular() && !tIt->IsMarkedForRegRef() )
            tIt->Close();
    }
    Comment("Closing grid " << Level << " done." << std::endl, DebugRefineEasyC);
}

#ifndef _PAR
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

        if ( !tIt->IsUnrefined() ){
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
#else
void MultiGridCL::UnrefineGrid (Uint Level)
{
    Comment("Unrefining grid " << Level << "." << std::endl, DebugRefineEasyC);

    const Uint nextLevel= Level+1;
    bool killedGhost= false;

    // mark all subsimplices on level 0 for removement. All subsimplices that
    // are needed any more, will be rescued within the refinement algorithm.
    if (Level==0)
    {
        std::for_each(_Vertices[0].begin(), _Vertices[0].end(), std::mem_fun_ref(&VertexCL::SetRemoveMark));
        std::for_each(_Edges[0].begin(),    _Edges[0].end(),    std::mem_fun_ref(&EdgeCL::SetRemoveMark));
        std::for_each(_Faces[0].begin(),    _Faces[0].end(),    std::mem_fun_ref(&FaceCL::SetRemoveMark));
    }

    std::for_each(_Vertices[nextLevel].begin(), _Vertices[nextLevel].end(), std::mem_fun_ref(&VertexCL::SetRemoveMark));
    std::for_each(_Edges[nextLevel].begin(),    _Edges[nextLevel].end(),    std::mem_fun_ref(&EdgeCL::SetRemoveMark));
    std::for_each(_Faces[nextLevel].begin(),    _Faces[nextLevel].end(),    std::mem_fun_ref(&FaceCL::SetRemoveMark));

    DynamicDataInterfaceCL::XferBegin();

    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; )
    {
        if (tIt->HasGhost() )
        {
            ++tIt;
            continue;
        }


        if ( !tIt->IsUnrefined() )
        {
            if ( tIt->IsMarkEqRule() )
                tIt->ClearAllRemoveMarks();
            else
            {
                std::for_each(tIt->GetChildBegin(), tIt->GetChildEnd(), std::mem_fun(&TetraCL::SetRemoveMark));

                // if tetra is ghost and will have no children on this proc after unref, we can delete this tetra
                if ( tIt->IsGhost() && tIt->IsMarkedForNoRef())
                {
                    if (!withUnknowns_)
                        tIt->XferDelete();
                    tIt->UnlinkFromFaces();
                    if (Level==0){  //mark subs for removement (are rescued after this loop)
                        std::for_each( tIt->GetVertBegin(), tIt->GetVertEnd(), std::mem_fun( &VertexCL::SetRemoveMark) );
                        std::for_each( tIt->GetEdgesBegin(), tIt->GetEdgesEnd(), std::mem_fun( &EdgeCL::SetRemoveMark) );
                        std::for_each( tIt->GetFacesBegin(), tIt->GetFacesEnd(), std::mem_fun( &FaceCL::SetRemoveMark) );
                    }
                    else{           // mark verts on level 0 for removement
                        for (TetraCL::const_VertexPIterator vert= tIt->GetVertBegin(), end= tIt->GetVertEnd(); vert!=end; ++vert){
                            if ((*vert)->GetLevel()==0 && !(*vert)->IsMaster()){
                                (*vert)->SetRemoveMark();
                            }
                        }
                    }
                    // remember tetra, that should be deleted
                    // !!! do not delete it now, because DDD needs still access to this tetra!!!
                    // But Prio is set to PrioKilledGhost, so HasGhost works still correct
                    ParMultiGridCL::PrioChange(&(*tIt), PrioKilledGhost);
                    toDelGhosts_.push_back(tIt++);
                    killedGhost= true;
                    continue;
                }
            }
        }
        ++tIt;
    }
    if (killedGhost)
        killedGhostTetra_=true;

    // now all removemarks are set. Now put simplices into recycle bin and clear remove marks of stilled used tetras
    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; )
    {
        if (tIt->HasGhost()){
            ++tIt;
            continue;
        }
        if ( (!tIt->IsUnrefined() || tIt->GetLevel()==0)
               && !tIt->IsMarkEqRule()
               && !tIt->IsMarkedForNoRef()
               && !tIt->IsMarkedForRemovement() )
        {
            // Maybe some subsimplices of Ghost-Tetras are marked for removement, so delete alle RemoveMarks on the Subs!
            if (tIt->GetRefRule()!=0)
                tIt->RecycleReusables();
            std::for_each( tIt->GetVertBegin(), tIt->GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
            std::for_each( tIt->GetEdgesBegin(), tIt->GetEdgesEnd(), std::mem_fun( &EdgeCL::ClearRemoveMark) );
            std::for_each( tIt->GetFacesBegin(), tIt->GetFacesEnd(), std::mem_fun( &FaceCL::ClearRemoveMark) );
        }

        ++tIt;
    }

    // rescue subs of level 0
    if (killedGhost)
    {
        for (TetraIterator tIt= _Tetras[0].begin(), tEnd= _Tetras[0].end(); tIt!=tEnd; ++tIt){
            if( tIt->IsGhost() ? !tIt->IsMarkedForNoRef() : !tIt->IsMarkedForRemovement() )
            { // rescue subs
                std::for_each( tIt->GetVertBegin(), tIt->GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
                std::for_each( tIt->GetEdgesBegin(), tIt->GetEdgesEnd(), std::mem_fun( &EdgeCL::ClearRemoveMark) );
                std::for_each( tIt->GetFacesBegin(), tIt->GetFacesEnd(), std::mem_fun( &FaceCL::ClearRemoveMark) );
            }
        }
    }

    // rescue subs that are owned by ghost tetras on next level
    ParMultiGridCL::TreatGhosts( nextLevel);
    // rescue verts special, because they can be found in different levels
    ParMultiGridCL::RescueGhostVerts( 0);

    /// \todo (of): Kann es einen Master aus hoeheren Leveln geben, der noch einen Knoten braucht, der hier im
    /// Zuge von killedGhost geloescht wird?

    // tell which simplices will be deleted
    // also if numerical data will be submitted after the refinement algorithm, no parallal information will be
    // needed on the simplices, because all datas of a tetra will be submitted. Hence this subsimplices can be
    // unsubscribe by the DDD-Sytem
    if (killedGhost)
    {
        for_each_if( _Faces[0].begin(), _Faces[0].end(),
                     std::mem_fun_ref(&FaceCL::XferDelete), std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
        for_each_if( _Edges[0].begin(), _Edges[0].end(),
                     std::mem_fun_ref(&EdgeCL::XferDelete), std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
        for_each_if( _Vertices[0].begin(), _Vertices[0].end(),
                     std::mem_fun_ref(&VertexCL::XferDelete), std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
    }

    // tetras on next level can be deleted anyway
    for_each_if( _Tetras[nextLevel].begin(), _Tetras[nextLevel].end(),
                 std::mem_fun_ref(&TetraCL::XferDelete), std::mem_fun_ref(&TetraCL::IsMarkedForRemovement) );

    // parallel infortmation about subsimplices aren't needed any more
    for_each_if( _Faces[nextLevel].begin(), _Faces[nextLevel].end(),
                 std::mem_fun_ref(&FaceCL::XferDelete), std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    for_each_if( _Edges[nextLevel].begin(), _Edges[nextLevel].end(),
                 std::mem_fun_ref(&EdgeCL::XferDelete), std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
    for_each_if( _Vertices[nextLevel].begin(), _Vertices[nextLevel].end(),
                 std::mem_fun_ref(&VertexCL::XferDelete), std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );

    // now DDD can internally delete references and information!
    DynamicDataInterfaceCL::XferEnd();

    // kill ghost tetras, that aren't need any more
    if (!withUnknowns_){
        for (std::list<TetraIterator>::iterator it=toDelGhosts_.begin(); it!=toDelGhosts_.end(); ++it)
            _Tetras[Level].erase(*it);
        toDelGhosts_.resize(0);
    }

    Comment("Now physically unlinking and removing superfluous tetras." << std::endl, DebugRefineEasyC);
    for (TetraIterator it= _Tetras[nextLevel].begin(), end= _Tetras[nextLevel].end(); it!=end; )
    {
        if ( it->IsMarkedForRemovement() )
        {
            if ( it->IsGhost() )
                throw DROPSErrCL("MultiGridCL::Unrefine: Ghost will be deleted, this is really strange!");
            it->UnlinkFromFaces();
            _Tetras[nextLevel].erase(it++);
        }
        else
            ++it;
    }

    /// \todo (of): Wieso hat Sven hier geschrieben, dass PrioVGhost-Kanten nicht die Referenz auf den MidVertex loeschen duerfen?
    Comment("Now adapting midvertex pointers on level " << Level << ". " << std::endl, DebugRefineEasyC);
    for (EdgeIterator eIt= _Edges[Level].begin(), eEnd= _Edges[Level].end(); eIt!=eEnd; ++eIt){
        if ( (eIt->IsRefined() /*&& eIt->GetHdr()->prio!=PrioVGhost*/ && eIt->GetMidVertex()->IsMarkedForRemovement()) ){
            eIt->RemoveMidVertex();
        }
    }


    // in difference to serial version, the edges, faces and vertices are deleted at other positions
    // all subsimplices will be deleted after the refinement is done. So DDD has (if needed) all the time
    // access to these unknowns.

    /// \todo (of) Muessen evtl. auch noch Vertices aus Level Level geloescht werden?

    Comment("Unrefining grid " << Level << " done." << std::endl, DebugRefineEasyC);
}
#endif

void MultiGridCL::RefineGrid (Uint Level)
{
    Comment("Refining grid " << Level << std::endl, DebugRefineEasyC);

#ifdef _PAR
    DynamicDataInterfaceCL::IdentifyBegin();    // new simplices must be identified by DDD
#endif

    const Uint nextLevel(Level+1);
    if ( Level==GetLastLevel() ) AppendLevel();

    for (TetraIterator tIt(_Tetras[Level].begin()), tEnd(_Tetras[Level].end()); tIt!=tEnd; ++tIt)
    {
        if ( tIt->IsMarkEqRule() ) continue;

        tIt->SetRefRule( tIt->GetRefMark() );
        if ( tIt->IsMarkedForNoRef() )
        {
#ifndef _PAR
            Comment("refining " << tIt->GetId().GetIdent() << " with rule 0." << std::endl, DebugRefineHardC);
#else
            AllComment("refining " << tIt->GetGID() << " with rule 0." << std::endl, DebugRefineHardC);
#endif
            if ( tIt->_Children )
                { delete tIt->_Children; tIt->_Children=0; }
        }
        else
#ifdef _PAR
            if ( !tIt->HasGhost() ) // refinement will be done on ghost tetra!
#endif
            {
                const RefRuleCL& refrule( tIt->GetRefData() );
#ifdef _PAR
                AllComment("refining " << tIt->GetGID() << " with rule " << tIt->GetRefRule() << "." << std::endl, DebugRefineHardC);
#else
                Comment("refining " << tIt->GetId().GetIdent() << " with rule " << tIt->GetRefRule() << "." << std::endl, DebugRefineHardC);
#endif
                tIt->CollectEdges           (refrule, _Vertices[nextLevel], _Edges[nextLevel], _Bnd);
                tIt->CollectFaces           (refrule, _Faces[nextLevel]);
                tIt->CollectAndLinkChildren (refrule, _Tetras[nextLevel]);
            }
    }
    for (Uint lvl= 0; lvl <= nextLevel; ++lvl)
        std::for_each( _Vertices[lvl].begin(), _Vertices[lvl].end(),
            std::mem_fun_ref( &VertexCL::DestroyRecycleBin));

    IncrementVersion();

#ifdef _PAR
    DynamicDataInterfaceCL::IdentifyEnd();
#endif

    Comment("Refinement of grid " << Level << " done." << std::endl, DebugRefineEasyC);
//    if (DROPSDebugC & DebugRefineHardC) if ( !IsSane(cdebug, Level) ) cdebug << std::endl;
}


void MultiGridCL::Refine()
{
#ifndef _PAR
    PeriodicEdgesCL perEdges( *this);
#endif
    ClearTriangCache();
    PrepareModify();
#ifdef _PAR
    killedGhostTetra_=false;
    withUnknowns_=ParMultiGridCL::UnknownsOnSimplices();
    if (!toDelGhosts_.empty())  // todo (of): als Assert schreiben!
        throw DROPSErrCL("MultiGridCL::Refine: toDelGhosts_ should be empty!");
    else
        toDelGhosts_.resize(0);
    ParMultiGridCL::AdjustLevel();      // all procs must have the same number of levels
#endif

    const int tmpLastLevel( GetLastLevel() );

    for (int Level=tmpLastLevel; Level>=0; --Level)
    {
        RestrictMarks(Level);
#ifndef _PAR
        perEdges.AccumulateMFR( Level);
#else
        // calc marks over proc boundaries
        ParMultiGridCL::CommunicateRefMarks( Level );
        ParMultiGridCL::AccumulateMFR( Level );
#endif
        CloseGrid(Level);
    }

    for (int Level=0; Level<=tmpLastLevel; ++Level)
    {
#ifndef _PAR
        if ( _Tetras[Level].empty() ) continue;
#endif
        if (Level)
            CloseGrid(Level);
        if ( Level != tmpLastLevel )
            UnrefineGrid(Level);
        RefineGrid(Level);
    }

#ifdef _PAR
    ParMultiGridCL::AdaptPrioOnSubs();

    // make killed ghost to all procs the same
    killedGhostTetra_= ProcCL::GlobalOr(killedGhostTetra_);

    // if no unknowns will be transfered after the refinement algorithm,
    // all subsimplices and killed ghosts can be deleted now
    if (!withUnknowns_)
    {
        for (Uint l=0; l<GetLastLevel(); ++l)
        {
            _Vertices[l].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
            _Edges[l].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
            _Faces[l].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
        }
        if (killedGhostTetra_){
            // todo of: Kann man das DynamicDataInterfaceCL::XferBegin und DDD_XferEnd aus der for-schleife herausziehen?
        	DynamicDataInterfaceCL::XferBegin();
            for (std::list<TetraIterator>::iterator it=toDelGhosts_.begin(); it!=toDelGhosts_.end(); ++it)
                (*it)->XferDelete();
            DynamicDataInterfaceCL::XferEnd();
            for (std::list<TetraIterator>::iterator it=toDelGhosts_.begin(); it!=toDelGhosts_.end(); ++it)
                _Tetras[(*it)->GetLevel()].erase(*it);
            toDelGhosts_.resize(0);
        }
        killedGhostTetra_= false;
    }
    while ( GetLastLevel()>0 && EmptyLevel(GetLastLevel()))
        RemoveLastLevel();

    ParMultiGridCL::AdjustLevel();
    FinalizeModify();
    ParMultiGridCL::MarkSimplicesForUnknowns();
    ClearTriangCache();

#else
    while ( _Tetras[GetLastLevel()].empty() ) RemoveLastLevel();
    FinalizeModify();
#endif

    std::for_each( GetAllVertexBegin(), GetAllVertexEnd(),
        std::mem_fun_ref( &VertexCL::DestroyRecycleBin));
}


void MultiGridCL::Scale( double s)
{
    for (VertexIterator it= GetAllVertexBegin(), end= GetAllVertexEnd();
        it!=end; ++it)
        it->_Coord*= s;
}

void MultiGridCL::Transform( Point3DCL (*mapping)(const Point3DCL&))
{
    for (VertexIterator it= GetAllVertexBegin(), end= GetAllVertexEnd();
        it!=end; ++it)
        it->_Coord= mapping(it->_Coord);
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


void SetAllEdges (TetraCL* tp, EdgeCL* e0, EdgeCL* e1, EdgeCL* e2, EdgeCL* e3, EdgeCL* e4, EdgeCL* e5)
{
    tp->SetEdge( 0, e0);
    tp->SetEdge( 1, e1);
    tp->SetEdge( 2, e2);
    tp->SetEdge( 3, e3);
    tp->SetEdge( 4, e4);
    tp->SetEdge( 5, e5);
}

void SetAllFaces (TetraCL* tp, FaceCL* f0, FaceCL* f1, FaceCL* f2, FaceCL* f3)
{
    tp->SetFace( 0, f0);
    tp->SetFace( 1, f1);
    tp->SetFace( 2, f2);
    tp->SetFace( 3, f3);
}

bool HasMultipleBndSegs (const TetraCL& t)
{
    Uint numbnd= 0;
    for (Uint i= 0; numbnd < 2 && i < NumFacesC; ++i)
        if (t.IsBndSeg( i)) ++numbnd;

    return numbnd > 1;
}

void MultiGridCL::SplitMultiBoundaryTetras()
{
    ClearTriangCache();
    PrepareModify();
    // Uint count= 0;

    // Note that new tetras can be appended to _Tetras[0] in this loop; it is assumed that
    // the pointers and iterators to the present tetras are not invalidated by appending
    // to _Tetras[0]. (This assumption must of course hold for the refinement algorithm to
    // work at all.) The new tetras never require further splitting.
    for (TetraLevelCont::iterator t= _Tetras[0].begin(), theend= _Tetras[0].end(); t != theend; ) {
        if (!HasMultipleBndSegs( *t)) {
            ++t;
            continue;
        }
        // ++count;

        // t: (v0 v1 v2 v3)

        // One new vertex: The barycenter b is in level 0 in the interior of \Omega; store it and its address.
        Point3DCL b( GetBaryCenter( *t));
        _Vertices[0].push_back( VertexCL( b, /*first level*/ 0));
        VertexCL* bp( &_Vertices[0].back());

        // Four new edges: (v0 b) (v1 b) (v2 b) (b v3); they have no midvertices and no boundary descriptions
        EdgeCL* ep[4];
        _Edges[0].push_back( EdgeCL( t->_Vertices[0], bp, /*level*/ 0));
        ep[0]= &_Edges[0].back();
        _Edges[0].push_back( EdgeCL( t->_Vertices[1], bp, /*level*/ 0));
        ep[1]= &_Edges[0].back();
        _Edges[0].push_back( EdgeCL( t->_Vertices[2], bp, /*level*/ 0));
        ep[2]= &_Edges[0].back();
        _Edges[0].push_back( EdgeCL( bp, t->_Vertices[3], /*level*/ 0));
        ep[3]= &_Edges[0].back();

        // Six new faces: (v0 v1 b) (v0 v2 b) (v0 b v3) (v1 v2 b) (v1 b v3) (v2 b v3); they have no boundary descriptions
        FaceCL* fp[6];
        for (int i= 0; i < 6; ++i) {
            _Faces[0].push_back( FaceCL( /*level*/ 0));
            fp[i]= &_Faces[0].back();
        }

        // Four new tetras: (v0 v1 v2 b) (v0 v1 b v3) (v0 v2 b v3) (v1 v2 b v3)
        TetraCL* tp[4];
        _Tetras[0].push_back( TetraCL( t->_Vertices[0], t->_Vertices[1], t->_Vertices[2], bp, /*parent*/ 0));
        tp[0]= &_Tetras[0].back();
        _Tetras[0].push_back( TetraCL( t->_Vertices[0], t->_Vertices[1], bp, t->_Vertices[3], /*parent*/ 0));
        tp[1]= &_Tetras[0].back();
        _Tetras[0].push_back( TetraCL( t->_Vertices[0], t->_Vertices[2], bp, t->_Vertices[3], /*parent*/ 0));
        tp[2]= &_Tetras[0].back();
        _Tetras[0].push_back( TetraCL( t->_Vertices[1], t->_Vertices[2], bp, t->_Vertices[3], /*parent*/ 0));
        tp[3]= &_Tetras[0].back();
        // Set the edge-pointers
        SetAllEdges( tp[0], t->_Edges[EdgeByVert(0, 1)], t->_Edges[EdgeByVert(0, 2)], t->_Edges[EdgeByVert(1, 2)], ep[0], ep[1], ep[2]);
        SetAllEdges( tp[1], t->_Edges[EdgeByVert(0, 1)], ep[0], ep[1], t->_Edges[EdgeByVert(0, 3)], t->_Edges[EdgeByVert(1, 3)], ep[3]);
        SetAllEdges( tp[2], t->_Edges[EdgeByVert(0, 2)], ep[0], ep[2], t->_Edges[EdgeByVert(0, 3)], t->_Edges[EdgeByVert(2, 3)], ep[3]);
        SetAllEdges( tp[3], t->_Edges[EdgeByVert(1, 2)], ep[1], ep[2], t->_Edges[EdgeByVert(1, 3)], t->_Edges[EdgeByVert(2, 3)], ep[3]);
        // Set the face-pointers
        SetAllFaces( tp[0], fp[3], fp[1], fp[0], t->_Faces[FaceByVert( 0, 1, 2)]);
        SetAllFaces( tp[1], fp[4], fp[2], t->_Faces[FaceByVert( 0, 1, 3)], fp[0]);
        SetAllFaces( tp[2], fp[5], fp[2], t->_Faces[FaceByVert( 0, 2, 3)], fp[1]);
        SetAllFaces( tp[3], fp[5], fp[4], t->_Faces[FaceByVert( 1, 2, 3)], fp[3]);

        // Set tetra-pointers of the new faces
        fp[0]->SetNeighbor( 0, tp[0]); fp[0]->SetNeighbor( 1, tp[1]);
        fp[1]->SetNeighbor( 0, tp[0]); fp[1]->SetNeighbor( 1, tp[2]);
        fp[2]->SetNeighbor( 0, tp[1]); fp[2]->SetNeighbor( 1, tp[2]);
        fp[3]->SetNeighbor( 0, tp[0]); fp[3]->SetNeighbor( 1, tp[3]);
        fp[4]->SetNeighbor( 0, tp[1]); fp[4]->SetNeighbor( 1, tp[3]);
        fp[5]->SetNeighbor( 0, tp[2]); fp[5]->SetNeighbor( 1, tp[3]);

        // Set tetra-pointers of the faces of t to the corresponding new tetra
        if (t->GetFace( 0)->GetNeighbor( 0) == &*t)
            t->_Faces[0]->SetNeighbor( 0, tp[3]);
        else
            t->_Faces[0]->SetNeighbor( 1, tp[3]);
        if (t->GetFace( 1)->GetNeighbor( 0) == &*t)
            t->_Faces[1]->SetNeighbor( 0, tp[2]);
        else
            t->_Faces[1]->SetNeighbor( 1, tp[2]);
        if (t->GetFace( 2)->GetNeighbor( 0) == &*t)
            t->_Faces[2]->SetNeighbor( 0, tp[1]);
        else
            t->_Faces[2]->SetNeighbor( 1, tp[1]);
        if (t->GetFace( 3)->GetNeighbor( 0) == &*t)
            t->_Faces[3]->SetNeighbor( 0, tp[0]);
        else
            t->_Faces[3]->SetNeighbor( 1, tp[0]);

        // Remove *t (now unused), increment t *before* erasing
        TetraLevelCont::iterator tmp= t;
        ++t;
        _Tetras[0].erase( tmp);
    }

    FinalizeModify();

    // std::cerr << "Split " << count << " tetras.\n";
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
        // Set all vertices, edges and faces to not needed
        for (const_VertexIterator sit( GetAllVertexBegin()); sit != GetAllVertexEnd( Level); ++sit)
            sit->SetNeeded(false);
        for (const_EdgeIterator sit( GetAllEdgeBegin()); sit != GetAllEdgeEnd( Level); ++sit)
            sit->SetNeeded(false);
        for (const_FaceIterator sit( GetAllFaceBegin()); sit != GetAllFaceEnd( Level); ++sit)
            sit->SetNeeded(false);

        // Check all levels
        for (int lvl=0; lvl<=static_cast<int>(GetLastLevel()); ++lvl)
            if ( !IsSane(os, lvl) )
                sane = false;

        // Check if all vertices, edges and faces are needed by at least one tetra
        for (const_VertexIterator sit( GetAllVertexBegin()); sit != GetAllVertexEnd( Level); ++sit)
            if (!sit->GetNeeded()){
                sane=false;
                os << "Not needed vertex:\n";
                sit->DebugInfo(os);
            }
        for (const_EdgeIterator sit( GetAllEdgeBegin()); sit != GetAllEdgeEnd( Level); ++sit)
            if (!sit->GetNeeded()){
                sane=false;
                os << "Not needed edge:\n";
                sit->DebugInfo(os);
            }
        for (const_FaceIterator sit( GetAllFaceBegin()); sit != GetAllFaceEnd( Level); ++sit)
            if (!sit->GetNeeded()){
                sane=false;
                os << "Not needed face:\n";
                sit->DebugInfo(os);
            }
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
#ifndef _PAR
    size_t numVerts= GetVertices().size(),
           numEdges= GetEdges().size(),
           numFaces= GetFaces().size(),
           numTetras= GetTetras().size(),
           numTetrasRef= numTetras - std::distance( GetTriangTetraBegin(), GetTriangTetraEnd());
    os << numVerts  << " Verts, "
       << numEdges  << " Edges, "
       << numFaces  << " Faces, "
       << numTetras << " Tetras"
       << std::endl;
#else
    int  elems[9],
        *recvbuf=0;
    if (ProcCL::IamMaster())
        recvbuf = new int[9*ProcCL::Size()];

    elems[0]= GetTriangVertex().size(); elems[1]= GetTriangEdge().size();
    elems[2]= GetTriangFace().size();   elems[3]= GetTriangTetra().size();
    elems[4] = GetVertices().size();    elems[5]=GetEdges().size();
    elems[6] = GetFaces().size();       elems[7]=GetTetras().size();
    elems[8] = elems[7] - std::distance( GetTriangTetraBegin(), GetTriangTetraEnd());

    ProcCL::Gather(elems, recvbuf, 9, ProcCL::Master());

    Uint numVerts=0, numEdges=0, numFaces=0, numTetras=0, numTetrasRef=0;
    if (ProcCL::IamMaster()){
        for (int i=0; i<ProcCL::Size(); ++i){
            numVerts  += recvbuf[i*9+4];
            numEdges  += recvbuf[i*9+5];
            numFaces  += recvbuf[i*9+6];
            numTetras += recvbuf[i*9+7];
            numTetrasRef += recvbuf[i*9+8];
        }
        os << "    On level " << GetLastLevel() << " there are:\n";
        for (int i=0; i<ProcCL::Size(); ++i){
            os << "     On Proc "<<i<<" are: "
               << recvbuf[i*9+0] << " Verts, "
               << recvbuf[i*9+1] << " Edges, "
               << recvbuf[i*9+2] << " Faces, "
               << recvbuf[i*9+3] << " Tetras"
               << '\n';
        }
        os << "  Accumulated (over all processes and levels): "
           << numVerts << " Verts, "
           << numEdges << " Edges, "
           << numFaces << " Faces, "
           << numTetras << " Tetras"
           << std::endl;
    }
    delete[] recvbuf;
#endif
    IF_MASTER
    {
        // print out memory usage.
        // before manipulating stream, remember previous precision and format flags
        const int prec= os.precision();
        const std::ios_base::fmtflags ff= os.flags();
        // print only one digit after decimal point
        os.precision(1);
        os.setf( std::ios_base::fixed);

        size_t vMem= numVerts*sizeof(VertexCL),
               eMem= numEdges*sizeof(EdgeCL),
               fMem= numFaces*sizeof(FaceCL),
               tMem= numTetras*sizeof(TetraCL) + numTetrasRef*8*sizeof(TetraCL*),
                   // also account for Children_ arrays which are allocated for all refined tetras
               Mem= vMem + eMem + fMem + tMem;
        double MemMB= double(Mem)/1024/1024;
        os << "Memory used for geometry: " << MemMB << " MB ("
           << (double(vMem)/Mem*100) << "% verts, "
           << (double(eMem)/Mem*100) << "% edges, "
           << (double(fMem)/Mem*100) << "% faces, "
           << (double(tMem)/Mem*100) << "% tetras)\n";
        // restore precision and format flags
        os.precision(prec);
        os.flags( ff);
    }
}

void MultiGridCL::ElemInfo(std::ostream& os, int Level) const
{
    double hmax= -1, hmin= 1e99,
           rmax= -1, rmin= 1e99;
    DROPS_FOR_TRIANG_CONST_TETRA( (*this), Level, It) {
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
#ifdef _PAR
    hmin = ProcCL::GlobalMin(hmin, ProcCL::Master());
    rmin = ProcCL::GlobalMin(rmin, ProcCL::Master());
    hmax = ProcCL::GlobalMax(hmax, ProcCL::Master());
    rmax = ProcCL::GlobalMax(rmax, ProcCL::Master());
#endif
    IF_MASTER
      os << hmin << " <= h <= " << hmax << '\t'
         << rmin << " <= h_max/h_min <= " << rmax << std::endl;
}

#ifdef _PAR
/// \brief Get number of distributed objects on local processor
Uint MultiGridCL::GetNumDistributedObjects() const
/** Count vertices, edges, faces and tetrahedra, that are stored on at least two
    processors. */
{
    Uint numdistVert=0, numdistEdge=0, numdistFace=0, numdistTetra=0;
    for (const_VertexIterator sit(GetVerticesBegin()), end(GetVerticesEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistVert;
    for (const_EdgeIterator sit(GetEdgesBegin()), end(GetEdgesEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistEdge;
    for (const_FaceIterator sit(GetFacesBegin()), end(GetFacesEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistFace;
    for (const_TetraIterator sit(GetTetrasBegin()), end(GetTetrasEnd()); sit!=end; ++sit)
        if (!sit->IsLocal()) ++numdistTetra;

    return numdistVert+numdistEdge+numdistFace+numdistTetra;
}

/// \brief Get number of tetrahedra of a given level
Uint MultiGridCL::GetNumTriangTetra(int Level)
{
    Uint numTetra=0;
    DROPS_FOR_TRIANG_TETRA( (*this), Level, It) {
        ++numTetra;
    }
    return numTetra;
}

/// \brief Get number of faces of a given level
Uint MultiGridCL::GetNumTriangFace(int Level)
{
    Uint numFace=0;
    DROPS_FOR_TRIANG_FACE( (*this), Level, It) {
        ++numFace;
    }
    return numFace;
}

/// \brief Get number of faces on processor boundary
Uint MultiGridCL::GetNumDistributedFaces(int Level)
{
    Uint numdistFace=0;
    DROPS_FOR_TRIANG_FACE( (*this), Level, It)
        if( It->IsOnProcBnd() )
            ++numdistFace;
    return numdistFace;
}
#endif

void
TriangFillCL<VertexCL>::fill (MultiGridCL& mg, TriangCL<VertexCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::VertexIterator it= mg.GetAllVertexBegin( lvl),
         theend= mg.GetAllVertexEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)
#ifdef _PAR
            && it->MayStoreUnk()
#endif
           )
            c.push_back( &*it);
    TriangCL<VertexCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<EdgeCL>::fill (MultiGridCL& mg, TriangCL<EdgeCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::EdgeIterator it= mg.GetAllEdgeBegin( lvl),
         theend= mg.GetAllEdgeEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)
  #ifdef _PAR
            && it->MayStoreUnk()
  #endif
           )
            c.push_back( &*it);
    TriangCL<EdgeCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<FaceCL>::fill (MultiGridCL& mg, TriangCL<FaceCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::FaceIterator it= mg.GetAllFaceBegin( lvl),
         theend= mg.GetAllFaceEnd( lvl); it != theend; ++it)
           if (it->IsInTriang( lvl)
  #ifdef _PAR
            && it->MayStoreUnk()
  #endif
           )
            c.push_back( &*it);
    TriangCL<FaceCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
TriangFillCL<TetraCL>::fill (MultiGridCL& mg, TriangCL<TetraCL>::LevelCont& c, int lvl)
{
    for (MultiGridCL::TetraIterator it= mg.GetAllTetraBegin( lvl),
         theend= mg.GetAllTetraEnd( lvl); it != theend; ++it)
        if (it->IsInTriang( lvl)
 #ifdef _PAR
            && it->IsMaster()
 #endif
           ) c.push_back( &*it);
    TriangCL<TetraCL>::LevelCont tmp= c;
    c.swap( tmp);
}

void
LocatorCL::LocateInTetra(LocationCL& loc, Uint trilevel, const Point3DCL&p, double tol)
// Assumes, that p lies in loc._Tetra and that loc._Tetra contains the barycentric
// coordinates of p therein. If these prerequisites are not met, this function might
// loop forever or lie to you. You have been warned!
// Searches p in the children of loc._Tetra up to triangulation-level trilevel
{
    const TetraCL*& t= loc._Tetra;
    SVectorCL<4>& b= loc._Coord;
    SMatrixCL<4,4> M;

    for (Uint lvl=t->GetLevel(); lvl<trilevel && !t->IsUnrefined(); ++lvl)
        {
            // Adjust relative tolerances on finer grid, so that the absolute
            // tolerances stay the same.
            tol *= 2;
            for (TetraCL::const_ChildPIterator it=t->GetChildBegin(), theend=t->GetChildEnd(); it!=theend; ++it)
            {
                MakeMatrix(**it, M);
                std::copy(p.begin(), p.end(), b.begin());
                b[3]= 1.;
                gauss_pivot(M, b);
                if ( InTetra(b, tol) )
                {
                    t= *it;
                    break;
                }
            }
        }
}

void
LocatorCL::Locate(LocationCL& loc, const MultiGridCL& MG, int trilevel, const Point3DCL& p, double tol)
/// \todo this only works for triangulations of polygonal domains, which resolve the geometry of the domain exactly (on level 0).
/// \todo this only works for FE-functions living on the finest level
{
    SVectorCL<4>& b= loc._Coord;
    SMatrixCL<4,4> M;
#ifndef _PAR
    const Uint search_level=0;
#else
    const Uint search_level=MG.GetLastLevel()-1;
#endif

    for (MultiGridCL::const_TetraIterator it= MG.GetTetrasBegin(search_level), theend= MG.GetTetrasEnd(search_level); it!=theend; ++it)
    {
        MakeMatrix(*it, M);
        std::copy(p.begin(), p.end(), b.begin());
        b[3]= 1.;
        gauss_pivot(M, b);
        if ( InTetra(b, tol) )
        {
            loc._Tetra= &*it;
            LocateInTetra(loc, MG.GetTriangTetra().StdIndex( trilevel), p, tol);
            return;
        }
    }
    loc._Tetra= 0; std::fill(b.begin(), b.end(), 0.);
}

void MarkAll (DROPS::MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, It)
        It->SetRegRefMark();
}


void UnMarkAll (DROPS::MultiGridCL& mg)
{
     DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, It)
     {
#ifdef _PAR
         if (!It->IsMaster()) std::cerr <<"Marking non-master tetra for removement!!!\n";
#endif
            It->SetRemoveMark();
     }
}

#ifdef _PAR
Ulint GetExclusiveVerts (const MultiGridCL &mg, Priority prio, int lvl)
{
    Ulint ret=0;
    for (MultiGridCL::const_TriangVertexIteratorCL it(mg.GetTriangVertexBegin(lvl)), end(mg.GetTriangVertexEnd(lvl)); it != end; ++it)
            if (it->IsExclusive(prio))
                ++ret;
    return ret;
}

Ulint GetExclusiveEdges(const MultiGridCL &mg, Priority prio, int lvl)
{
    Ulint ret=0;
    for (MultiGridCL::const_TriangEdgeIteratorCL it(mg.GetTriangEdgeBegin(lvl)), end(mg.GetTriangEdgeEnd(lvl)); it != end; ++it)
            if (it->IsExclusive(prio))
                ++ret;
    return ret;
}

std::string PrioToString(Uint prio)
{
    switch (prio)
    {
        case PrioNeutral:
            return "PrioNeutral";
        case PrioKilledGhost:
            return "PrioKilledGhost";
        case PrioVGhost:
            return "PrioVGhost";
        case PrioGhost:
            return "PrioGhost";
        case PrioMaster:
            return "PrioMaster";
        case PrioHasUnk:
            return "PrioHasUnk";
    }
    return "UnknownPrio";
}
#endif


void ColorClassesCL::compute_neighbors (MultiGridCL::const_TriangTetraIteratorCL begin,
                                        MultiGridCL::const_TriangTetraIteratorCL end,
                                        std::vector<TetraNumVecT>& neighbors)
{
    const size_t num_tetra= std::distance( begin, end);

    typedef std::tr1::unordered_map<const VertexCL*, TetraNumVecT> VertexMapT;
    VertexMapT vertexMap;
    // Collect all tetras, that have vertex v in vertexMap[v].
    for (MultiGridCL::const_TriangTetraIteratorCL sit= begin; sit != end; ++sit)
        for (int i= 0; i < 4; ++i)
            vertexMap[sit->GetVertex( i)].push_back( sit - begin);

    // For every tetra j, store all neighboring tetras in neighbors[j].
    typedef std::set<size_t> TetraNumSetT;
    std::vector<TetraNumSetT> neighborsets( num_tetra);
#   pragma omp parallel
    {
#ifndef DROPS_WIN
        size_t j;
#else
        int j;
#endif
#       pragma omp for
        for (j= 0; j < num_tetra; ++j)
            for (int i= 0; i < 4; ++i) {
                const TetraNumVecT& tetra_nums= vertexMap[(begin + j)->GetVertex( i)];
                neighborsets[j].insert( tetra_nums.begin(), tetra_nums.end());
            }
#       pragma omp for
        for (j= 0; j < num_tetra; ++j) {
            neighbors[j].resize( neighborsets[j].size());
            std::copy( neighborsets[j].begin(), neighborsets[j].end(), neighbors[j].begin());
        }
    }
}

void ColorClassesCL::fill_pointer_arrays (
    const std::list<ColorFreqT>& color_list, const std::vector<int>& color,
    MultiGridCL::const_TriangTetraIteratorCL begin, MultiGridCL::const_TriangTetraIteratorCL end)
{
    colors_.resize( color_list.size());
    for (std::list<ColorFreqT>::const_iterator it= color_list.begin(); it != color_list.end(); ++it)
        colors_[it->first].reserve( it->second);
    const size_t num_tetra= std::distance( begin, end);
    for (size_t j= 0; j < num_tetra; ++j)
        colors_[color[j]].push_back( &*(begin + j));

#ifndef DROPS_WIN
    size_t j;
#else
    int j;
#endif
    // tetra sorting for better memory access pattern
    #pragma omp parallel for
    for (j= 0; j < num_colors(); ++j)
        sort( colors_[j].begin(), colors_[j].end());
}

void ColorClassesCL::compute_color_classes (MultiGridCL::const_TriangTetraIteratorCL begin,
                                            MultiGridCL::const_TriangTetraIteratorCL end)
{
#   ifdef _PAR
        ParTimerCL timer;
#   else
        TimerCL timer;
#   endif
        timer.Start();

    const size_t num_tetra= std::distance( begin, end);

    // Build the adjacency lists (a vector of neighbors for each tetra).
    std::vector<TetraNumVecT> neighbors( num_tetra);
    compute_neighbors( begin, end, neighbors);

    // Color the tetras
    std::vector<int> color( num_tetra, -1); // Color of each tetra
    std::list<ColorFreqT> color_frequency;  // list of colors together with number of their occurrence
    std::vector<int> used_colors; // list of the colors of all neighbors (multiple occurrences of the same color or -1 are allowed)
    for (size_t j= 0; j < num_tetra; ++j) {
        for (TetraNumVecT::iterator neigh_it= neighbors[j].begin(); neigh_it != neighbors[j].end(); ++neigh_it)
            used_colors.push_back( color[*neigh_it]);
        bool color_found= false;
        std::list<ColorFreqT>::iterator it;
        for (it= color_frequency.begin(); it != color_frequency.end(); ++it)
            if (find( used_colors.begin(), used_colors.end(), it->first) == used_colors.end()) {
                color_found= true;
                break;
            }
        if (color_found) {
            color[j]= it->first;
            ++it->second;
            // Move color to the end: LRU-policy for evenly used colors.
            color_frequency.splice( color_frequency.end(), color_frequency, it);
        }
        else {
            color_frequency.push_back( std::make_pair( color_frequency.size(), 1)); // Add new color with one use
            color[j]= color_frequency.back().first;
        }
        used_colors.clear();
    }
    neighbors.clear();

    // Build arrays of pointers for the colors
    fill_pointer_arrays( color_frequency, color, begin, end);
    color.clear();

    // for (size_t j= 0; j < num_colors(); ++j)
    //     std::cout << "Color " << j << " has " << colors_[j].size() << " tetras." << std::endl;
    // std::cout << std::endl;

    timer.Stop();
    const double duration= timer.GetTime();
    std::cout << "Creation of the tetra-coloring took " << duration << " seconds, " << num_colors() << " colors used." << '\n';
}

const ColorClassesCL& MultiGridCL::GetColorClasses (int Level) const
{
    if (Level < 0)
        Level+= GetNumLevel();

    if (_colors.find( Level) == _colors.end())
        _colors[Level]= new ColorClassesCL( GetTriangTetraBegin( Level), GetTriangTetraEnd( Level));

    return *_colors[Level];
}

} // end of namespace DROPS
