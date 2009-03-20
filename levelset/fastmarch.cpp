//**************************************************************************
// File:    fastmarch.cpp                                                  *
// Content: fast marching method for reparametrization                     *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************

#include "levelset/fastmarch.h"
#include <fstream>
#include <cstring>
#ifdef _OPENMP
#  include <omp.h>
#endif


namespace DROPS
{

#ifndef _PAR
double FastMarchCL::CompValueProj( IdxT Nr, int num, const IdxT upd[3]) const
{
    double val= 1e99;

    switch (num)
    {
        case 2: // Projektion auf Edge
        {
            const Point3DCL a= Coord_[upd[1]] - Coord_[upd[0]];
            const Point3DCL b= Coord_[  Nr  ] - Coord_[upd[0]];
            const double bary= inner_prod(a,b)/a.norm_sq();
            if (bary>=0 && bary<=1)
            {
                const Point3DCL lotfuss= (1-bary)*Coord_[upd[0]] + bary*Coord_[upd[1]];
                const double y= (1-bary)*v_->Data[upd[0]] + bary*v_->Data[upd[1]];
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
        break;

        case 3: // Projektion auf Face
        {
            const Point3DCL a= Coord_[upd[1]] - Coord_[upd[0]];
            const Point3DCL b= Coord_[upd[2]] - Coord_[upd[0]];
            const Point3DCL c= Coord_[  Nr  ] - Coord_[upd[0]];
            const double bary1= inner_prod(a,c)/a.norm_sq(),
                         bary2= inner_prod(b,c)/b.norm_sq();
            if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
            {
                const Point3DCL lotfuss= (1-bary1-bary2)*Coord_[upd[0]] + bary1*Coord_[upd[1]] + bary2*Coord_[upd[2]];
                const double y= (1-bary1-bary2)*v_->Data[upd[0]] + bary1*v_->Data[upd[1]] + bary2*v_->Data[upd[2]];
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
    }

    return val;
}
#else
double FastMarchCL::CompValueProj( IdxT Nr, int num, const IdxT upd[3]) const
{
    double val= 1e99;
    if (num<2)
        return val;

    std::valarray<Point3DCL> coord(3);
    coord[0]=GlobalCoord_[upd[0]];
    coord[1]=GlobalCoord_[upd[1]];
    if (num>2)
        coord[2]=GlobalCoord_[upd[2]];
    switch (num)
    {
        case 2: // Projektion auf Edge
        {
            const Point3DCL a= coord[1] - coord[0];
            const Point3DCL b= GlobalCoord_[Nr] - coord[0];
            const double bary= inner_prod(a,b)/a.norm_sq();
            if (bary>=0 && bary<=1)
            {
                const Point3DCL lotfuss= (1-bary)*coord[0] + bary*coord[1];
                const double y= (1-bary)*GlobalV_[upd[0]] + bary*GlobalV_[upd[1]];
                val= y + (lotfuss - GlobalCoord_[Nr]).norm();
            }
        }
        break;

        case 3: // Projektion auf Face
        {
            const Point3DCL a= coord[1] - coord[0];
            const Point3DCL b= coord[2] - coord[0];
            const Point3DCL c= GlobalCoord_[Nr] - coord[0];
            const double bary1= inner_prod(a,c)/a.norm_sq(),
                         bary2= inner_prod(b,c)/b.norm_sq();
            if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
            {
                const Point3DCL lotfuss= (1-bary1-bary2)*coord[0] + bary1*coord[1] + bary2*coord[2];
                const double y= (1-bary1-bary2)*GlobalV_[upd[0]] + bary1*GlobalV_[upd[1]] + bary2*GlobalV_[upd[2]];
                val= y + (lotfuss - GlobalCoord_[Nr]).norm();
            }
        }
    }

    return val;
}
#endif


void FastMarchCL::InitZero( bool ModifyZero)
/// \param[in] ModifyZero If this flag is set, the values around the zero level of the levelset function are new computed. Otherwise
///                       the old values are kept.
{
    Comment("Init zero\n", DebugParallelNumC);
    // Knoten an der Phasengrenze als Finished markieren
    // und Distanz zur Phasengrenze bestimmen (falls ModifyZero)
    const Uint idx= v_->RowIdx->GetIdx();

    int        sign[10];
    int        num_sign[3]; // - 0 +
    IdxT       Numb[10];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    // init Coord_
    Coord_.resize( size_);

#ifdef _PAR
    Typ_.resize(size_, Far);
    tmpTyp_.resize(size_, Far);
    tmpv_.resize(size_);
#endif

    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(), end=MG_.GetTriangVertexEnd();
        it!=end; ++it)
        Coord_[it->Unknowns(idx)]= it->GetCoord();
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(), end=MG_.GetTriangEdgeEnd();
        it!=end; ++it)
        Coord_[it->Unknowns(idx)]= GetBaryCenter( *it);

    // store copy of v_.Data in Old_
    Old_.resize( size_);
    Old_= v_->Data;
    VecDescCL oldv( *v_);

    LocalP2CL<> PhiLoc;

    for (MultiGridCL::TriangTetraIteratorCL it=MG_.GetTriangTetraBegin(), end=MG_.GetTriangTetraEnd();
        it!=end; ++it)
    {
        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            Numb[v]= v<4 ? it->GetVertex(v)->Unknowns(idx)
                         : it->GetEdge(v-4)->Unknowns(idx);
            sign[v]= std::abs(Old_[Numb[v]])<1e-8 ? 0 : (Old_[Numb[v]]>0 ? 1 : -1);
            if (sign[v]==0)
                Typ_[Numb[v]]= Finished;
        }
        PhiLoc.assign( *it, oldv, NoBndDataCL<>());

        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];

            const bool intersec= (num_sign[0]*num_sign[2]!=0); // Vorzeichenwechsel

            if (!intersec) continue;

            if (!ModifyZero)
            {
                for (int vert= 0; vert<4; ++vert)
                {
                    const IdxT Nr= Numb[data.Vertices[vert]];
                    Typ_[Nr]= Finished;
                    v_->Data[Nr]= std::abs( Old_[Nr]);
                }
                continue;
            }

            // ab hier gilt intersec && ModifyZero == true

            Point3DCL Schnitt[4];
            int num= 0;
            // Berechnung der Schnittpunkte mit den Kanten des Tetra
            for (int vert= 0; vert<4; ++vert)
                if (sign[data.Vertices[vert]]==0)
                    Schnitt[num++]= Coord_[Numb[data.Vertices[vert]]];

            for (int edge= 0; edge<6 && num<4; ++edge)
            {
                const Ubyte v1= data.Vertices[ VertOfEdge( edge, 0)],
                            v2= data.Vertices[ VertOfEdge( edge, 1)];
                if (sign[v1]*sign[v2] == -1) // Vorzeichenwechsel auf edge
                {
                    const IdxT Nr1= Numb[v1],
                               Nr2= Numb[v2];
                    const double bary= InterfacePatchCL::EdgeIntersection( v1,v2, PhiLoc);
                    Schnitt[num++]= (1-bary)*Coord_[Nr1] + bary*Coord_[Nr2];
                }
            }
/*
fil << "geom {OFF " << num << " 1 0\n";
for (int i=0; i<num; ++i)
{
    for (int j=0; j<3; ++j)
        fil << Schnitt[i][j] << ' ';
    fil << '\n';
}
if (num==3)
    fil << "3 0 1 2";
else
    fil << "4 0 1 3 2";
fil << "\n}\n";
*/
            if (num<3) throw DROPSErrCL("FastMarchCL::InitZero: intersection missing");

            for (int repeat=0; repeat<num-2; ++repeat)
            { // fuer num==4 (Schnitt ABDC ist viereckig)
              // zwei Dreiecke ABC + DBC betrachten
                if (repeat) Schnitt[0]= Schnitt[3];

                const Point3DCL a= Schnitt[1] - Schnitt[0],
                                b= Schnitt[2] - Schnitt[0];

                for (int vert=0; vert<4; ++vert)
                {
                    if (sign[data.Vertices[vert]]==0) continue;

                    const IdxT Nr= Numb[data.Vertices[vert]];
                    const Point3DCL Crd= Coord_[Nr],
                                    c=   Crd - Schnitt[0];
                    double dist= std::min( c.norm(), (Crd-Schnitt[1]).norm());
                    dist= std::min( dist, (Crd-Schnitt[2]).norm());

                    const double bary1= inner_prod(a,c)/a.norm_sq(),
                                 bary2= inner_prod(b,c)/b.norm_sq();
                    if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
                    {
                       const Point3DCL lotfuss= (1-bary1-bary2)*Schnitt[0] + bary1*Schnitt[1] + bary2*Schnitt[2];
                       dist= std::min( dist, (lotfuss - Crd).norm());
                    }

                    if (Typ_[Nr] != Finished)
                    {
                        Typ_[Nr]= Finished;
                        v_->Data[Nr]= dist;
                    }
                    else{
                        v_->Data[Nr]= std::min( dist, v_->Data[Nr]);
                    }
                }
            }
        }
    }
//fil << "}\n";
}

#ifndef _PAR
void FastMarchCL::InitClose()
{
    // an Finished angrenzende Knoten mit Close markieren und dort v_ updaten
    const Uint idx= v_->RowIdx->GetIdx();

    neigh_.resize( size_);

    IdxT Numb[10];

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    for (MultiGridCL::TriangTetraIteratorCL it=MG_.GetTriangTetraBegin(), end=MG_.GetTriangTetraEnd();
        it!=end; ++it)
    {
        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            if (v<4)
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            else
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);
        }

        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            ReprTetraT t;
            bool ContainsFinished= false;
            for (int vert= 0; vert<4; ++vert)
            {
                const IdxT Nr= Numb[ data.Vertices[vert]];
                t[vert]= Nr;
                if (Typ_[Nr] == Finished)
                    ContainsFinished= true;
            }

            if (ContainsFinished)
                for (int vert= 0; vert<4; ++vert)
                    Update( t[vert]);

            // init neigh_
            for (int vert= 0; vert<4; ++vert)
                neigh_[t[vert]].push_back( t);
        }
    }
}
#else
void FastMarchCL::InitClose()
/** Init Close is only done by master-process. All DoF that have finished neighbor DoF
    are marked as close and the values are updated */
{
    Comment("Init close\n", DebugParallelNumC);
    if (!ProcCL::IamMaster())
        return;

    for (IdxT i=0; i<GlobNeigh_.size(); ++i)
    {
        bool ContainsFinished= false;
        for (Uint n=0; n<GlobNeigh_[i].size(); ++n)
        {
            // Check if a finished vert lies in tetra
            for (int vert= 0; vert<4; ++vert)
                if (GlobalTyp_[ GlobNeigh_[i][n][vert] ] == Finished)
                    ContainsFinished= true;

            // If contains finihed, update neighs
            if (ContainsFinished)
                for (int vert= 0; vert<4; ++vert)
                    Update( GlobNeigh_[i][n][vert] );
        }
    }
}
#endif

IdxT FastMarchCL::FindTrial() const
{
    double min= 1e99;
    IdxT min_idx= NoIdx;

    for (std::set<IdxT>::const_iterator it= Close_.begin(), end= Close_.end(); it!=end; ++it)
    {
#ifndef _PAR
        if (v_->Data[*it]<=min)
        {
            min= v_->Data[*it];
            min_idx= *it;
        }
#else
        if(GlobalV_[*it]<=min)
        {
            min= GlobalV_[*it];
            min_idx= *it;
        }
#endif
    }
    return min_idx;
}

#ifndef _PAR
void FastMarchCL::Update( const IdxT NrI)
{
    // Update all vertices that are not Finished
    if (Typ_[NrI] == Finished) return;

    IdxT upd[3];
    double minval= Typ_[NrI]==Close ? v_->Data[NrI] : 1e99;

    for (Uint n=0; n<neigh_[NrI].size(); ++n)
    {
        int num= 0;
        for (int j=0; j<4; ++j)
        {
            const IdxT NrJ= neigh_[NrI][n][j];
            if (Typ_[NrJ] == Finished)
            {
                upd[num++]= NrJ;
                minval= std::min( minval, v_->Data[NrJ] + (Coord_[NrJ]-Coord_[NrI]).norm());
            }
        }

        minval= std::min( minval, CompValueProj( NrI, num, upd));
    }

    v_->Data[NrI]= minval;
    if (Typ_[NrI] != Close)
    {
        Close_.insert( NrI);
        Typ_[NrI]= Close;
    }
}
#else
void FastMarchCL::Update( const IdxT NrI)
{
    // Update all vertices that are not Finished
    if (GlobalTyp_[NrI] == Finished) return;

    IdxT upd[3];
    double minval= GlobalTyp_[NrI]==Close ? GlobalV_[NrI] : 1e99;

    for (Uint n=0; n<GlobNeigh_[NrI].size(); ++n)
    {
        int num= 0;
        for (int j=0; j<4; ++j)
        {
            const IdxT NrJ= GlobNeigh_[NrI][n][j];
            Assert(NrJ<globsize_, DROPSErrCL("FastMarchCL::Update: DoF-Number to big!"), DebugParallelNumC);
            if (GlobalTyp_[NrJ] == Finished)
            {
                upd[num++]= NrJ;
                minval= std::min( minval, GlobalV_[NrJ] + (GlobalCoord_[NrJ]-GlobalCoord_[NrI]).norm());
            }
        }
        minval= std::min( minval, CompValueProj( NrI, num, upd));
    }

    GlobalV_[NrI]= minval;
    if (GlobalTyp_[NrI] != Close)
    {
        Close_.insert( NrI);
        GlobalTyp_[NrI]= Close;
    }
}
#endif

#ifndef _PAR
void FastMarchCL::Reparam( bool ModifyZero)
{
    TimerCL time;
    InitZero( ModifyZero);
    InitClose();

    IdxT next;

    while ((next= FindTrial()) != NoIdx)
    {
        Close_.erase( next);
        Typ_[next]= Finished;

        std::set<IdxT> neighVerts;
        for (Uint n=0; n<neigh_[next].size(); ++n)
        { // collect all neighboring verts in neighVerts
            for (Uint i=0; i<4; ++i)
                neighVerts.insert( neigh_[next][n][i]);
        }
        for (std::set<IdxT>::const_iterator it= neighVerts.begin(), end= neighVerts.end();
            it!=end; ++it)
        { // update all neighboring verts, mark as Close
            Update( *it);
        }
        neigh_[next].clear(); // will not be needed anymore
    }

    RestoreSigns();
    time.Stop();
    std::cout << " reparametrization took " <<time.GetTime() << " s" << std::endl;
}
#else
void FastMarchCL::Reparam( bool ModifyZero)
{
    ParTimerCL time;
    Comment("Reparametrization\n", DebugParallelNumC);

    // Create global numbering, for neighborhood
    CreateGlobNumb();

    // Init Zero-Level of levelset function
    InitZero(ModifyZero);

    // tell other procs about finished-marked DoF's
    DistributeFinished();

    // Init local neighborhood
    InitLocNeigh();

    // Collect information on master-proc
    Collect();

    if (!ProcCL::IamMaster())
        v_->Data=0.;
    else
    {
        InitClose();
        IdxT next;

        while ((next= FindTrial()) != NoIdx)
        {
            Close_.erase( next);
            GlobalTyp_[next]= Finished;
            std::set<IdxT> neighVerts;

            // collect all neighboring verts in neighVerts
            for (Uint n=0; n<GlobNeigh_[next].size(); ++n)
                for (Uint i=0; i<4; ++i)
                    neighVerts.insert( GlobNeigh_[next][n][i]);

            // update all neighboring verts, mark as Close
            for (std::set<IdxT>::const_iterator it= neighVerts.begin(), end= neighVerts.end(); it!=end; ++it)
                Update( *it);

            GlobNeigh_[next].clear();
        }
    }
    // Distribute result from master proc to other procs
    Distribute();

    // result is not accumulated so far
    ex_->Accumulate(v_->Data);

    // Assign old signs
    RestoreSigns();
    CleanUp();
    time.Stop();
    std::cout << " reparametrization took " <<time.GetTime() << " s" << std::endl;
}
#endif

void FastMarchCL::ReparamEuklid( bool ModifyZero)
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    Comment("Reparametrization\n", DebugParallelNumC);

     // Init Zero-Level of levelset function
    InitZero(ModifyZero);
#ifndef _PAR
    InitZeroSet(CoordZeroLevel_, ValueZeroLevel_);
#else
    // tell other procs about finished-marked DoF's
    DistributeFinished();

    // Distribute zerolevel
    DistributeZeroLevel();
#endif

    Comment("Calculate Euclidian distance\n", DebugParallelNumC);

#pragma omp parallel
{
#ifdef _OPENMP
#pragma omp master
    {
        std::cout << "   * Using "<<omp_get_num_threads()<<" thread(s) to compute distance ..." << std::endl;
    }
#endif

#pragma omp for
    for(int i=0; i<(int)size_; ++i)
        if (Typ_[i]!=Finished)
            v_->Data[i]=MinDist(i);

    RestoreSigns();
}
    time.Stop();
    std::cout << " reparametrization took " <<time.GetTime() << " s" << std::endl;
}

void FastMarchCL::RestoreSigns()
{ // restore signs of v_
#ifdef _OPENMP
#pragma omp master
    std::cout << "   * Using "<<omp_get_num_threads()<<" thread(s) to restore signs ..." << std::endl;
#endif

#pragma omp for schedule( static)
    for (int i=0; i<(int)Old_.size(); ++i){
        if (Old_[i]<0){
            v_->Data[i]*= -1;
        }
    }
}

void FastMarchCL::InitZeroSet(VectorCL& CoordZeroLevel, VectorCL& ValueZeroLevel)
{
    size_t NumFinished=0;
    for (IdxT i=0; i<size_; ++i)
        if (Typ_[i]==Finished){
#ifdef _PAR
            if (v_->RowIdx->GetEx().IsExclusive(i))
#endif
                ++NumFinished;
        }

    CoordZeroLevel.resize( 3*NumFinished);
    ValueZeroLevel.resize( NumFinished);
    IdxT pos1=0, pos2=0;
    for (IdxT i=0; i<size_; ++i){
        if (Typ_[i]==Finished){
#ifdef _PAR
            if (v_->RowIdx->GetEx().IsExclusive(i))
#endif
            {
                for (int j=0; j<3; ++j)
                    CoordZeroLevel[pos1++]=Coord_[i][j];
                ValueZeroLevel[pos2++]=v_->Data[i];
            }
        }
    }
}

double FastMarchCL::MinDist(IdxT nr)
{
    double min=1e99;
    for (IdxT i=0; i<ValueZeroLevel_.size(); ++i)
    {
        double dist=  (MakePoint3D( CoordZeroLevel_[3*i+0], CoordZeroLevel_[3*i+1], CoordZeroLevel_[3*i+2])-Coord_[nr]).norm()
                    + ValueZeroLevel_[i];
        min = std::min(min, dist);
    }
    Assert(min!=1e99, DROPSErrCL("FastMarchCL::MinDist: No minimal distance found"), DebugParallelNumC);
    return min;
}

// ------------------------------------------
// only parallel functions
// ------------------------------------------
#ifdef _PAR
// Init of the static variables
//--------------------------------------------
VectorBaseCL<byte> FastMarchCL::tmpTyp_      = VectorBaseCL<byte>();
VectorBaseCL<byte> FastMarchCL::Typ_         = VectorBaseCL<byte>();
VecDescCL*         FastMarchCL::v_           = 0;
VectorCL           FastMarchCL::tmpv_        = VectorCL();
IdxT               FastMarchCL::offset_      = 0;
VectorBaseCL<IdxT> FastMarchCL::allExclusive_= VectorBaseCL<IdxT>();
std::vector<IdxT>  FastMarchCL::globNumb_    = std::vector<IdxT>();
std::vector<IdxT>  FastMarchCL::locNumb_     = std::vector<IdxT>();

// Definition of Gather and Scatter-functions
//--------------------------------------------
template<typename SimplexT>
  int FastMarchCL::HandlerFinishedGather(DDD_OBJ objp, void* buf)
/** On sender-side collect typ and value of distributed DoF*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);              // pointer to simplex
    CoupMarkValST* buffer = static_cast<CoupMarkValST*>(buf);   // pointer to coupling of mark and value

    if (!sp->Unknowns.Exist() || !sp->Unknowns.Exist(v_->RowIdx->GetIdx()))
        return 1;

    Uint Nr= sp->Unknowns(v_->RowIdx->GetIdx());                // number of DoF

    buffer->mark= tmpTyp_[Nr];                              // set typ
    buffer->val=  tmpv_[Nr];                                // set value

    return 0;
}

template<typename SimplexT>
  int FastMarchCL::HandlerFinishedScatter(DDD_OBJ objp, void* buf)
/** On receiver side, update value and mark of finished DoFs*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    CoupMarkValST* buffer = static_cast<CoupMarkValST*>(buf);

    if (!sp->Unknowns.Exist() || !sp->Unknowns.Exist(v_->RowIdx->GetIdx()))
        return 1;

    Uint Nr= sp->Unknowns(v_->RowIdx->GetIdx());

    // if a finished dof is transfered
    if (buffer->mark==Finished)
    {
        // if dof is marked before transfer as finished or within transfer phase
        if (tmpTyp_[Nr]==Finished || Typ_[Nr]==Finished)
            v_->Data[Nr] = std::min(v_->Data[Nr], buffer->val);
        else
            v_->Data[Nr] = buffer->val;
        // set dof as finished
        Typ_[Nr]=Finished;
    }

    return 0;
}

template<typename SimplexT>
  int FastMarchCL::HandlerGlobDOFGather(DDD_OBJ objp, void* buf)
/** On sender side collect global number of dof on simplex (if not there send NoIdx)*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    IdxT* buffer = static_cast<IdxT*>(buf);
    if (sp->Unknowns.Exist() && sp->Unknowns.Exist(v_->RowIdx->GetIdx()))
        *buffer= globNumb_[sp->Unknowns(v_->RowIdx->GetIdx())];      // may be NoIdx if simplex is not exclusive
    return 0;
}

template<typename SimplexT>
  int FastMarchCL::HandlerGlobDOFScatter(DDD_OBJ objp, void* buf)
/** On recieved side collect global number of dof on simplex if sender has send a number*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    const IdxT* buffer = static_cast<IdxT*>(buf);
    if (*buffer!=NoIdx){
        if (sp->Unknowns.Exist() && sp->Unknowns.Exist(v_->RowIdx->GetIdx())){
            Assert(globNumb_[sp->Unknowns(v_->RowIdx->GetIdx())]==NoIdx, DROPSErrCL("FastMarchCL::HandlerGlobDOFScatter: Two exclusive simplices found!"), DebugParallelNumC);
            globNumb_[sp->Unknowns(v_->RowIdx->GetIdx())]= *buffer;
        }
    }
    return 0;
}


// Definition of the wrappers
extern "C" int HandlerFinishedGatherVertexC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerFinishedGather<VertexCL>(objp,buf);
}
extern "C" int HandlerFinishedGatherEdgeC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerFinishedGather<EdgeCL>(objp,buf);
}
extern "C" int HandlerFinishedScatterVertexC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerFinishedScatter<VertexCL>(objp,buf);
}
extern "C" int HandlerFinishedScatterEdgeC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerFinishedScatter<EdgeCL>(objp,buf);
}

extern "C" int HandlerGlobDOFGatherVertexC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerGlobDOFGather<VertexCL>(objp, buf);
}
extern "C" int HandlerGlobDOFGatherEdgeC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerGlobDOFGather<EdgeCL>(objp, buf);
}
extern "C" int HandlerGlobDOFScatterVertexC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerGlobDOFScatter<VertexCL>(objp,buf);
}
extern "C" int HandlerGlobDOFScatterEdgeC(DDD_OBJ objp, void* buf){
    return FastMarchCL::HandlerGlobDOFScatter<EdgeCL>(objp,buf);
}

/// \brief Distribute "Finished"-marks and values on these DoF
void FastMarchCL::DistributeFinished()
/// this function uses the DDD-Interfaces of master vertices and edges
{
    Comment("Distribute Finished\n", DebugParallelNumC);
    // copy data of typ
//     Assert(tmpTyp_.size()!=size_, DROPSErrCL("FastMarchCL::DistributeFinished: tmpType hast not the right length!"), DebugParallelNumC);
    tmpTyp_.resize(size_);
    tmpv_.resize(size_);
    tmpTyp_= Typ_;
    tmpv_  =v_->Data;

    DDD_IFExchange(InterfaceCL<VertexCL>::GetIF(),  sizeof(CoupMarkValST),
                   HandlerFinishedGatherVertexC,   HandlerFinishedScatterVertexC );
    DDD_IFExchange(InterfaceCL<EdgeCL>::GetIF(), sizeof(CoupMarkValST),
                   HandlerFinishedGatherEdgeC,   HandlerFinishedScatterEdgeC );
}

/// \brief Init local neighborhood
void FastMarchCL::InitLocNeigh()
{
    Comment("Init local neighborhood\n", DebugParallelNumC);
    const Uint idx= v_->RowIdx->GetIdx();
//     neigh_.resize( size_);
    IdxT Numb[10];
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    // init neigh local
    for (MultiGridCL::TriangTetraIteratorCL it=MG_.GetTriangTetraBegin(), end=MG_.GetTriangTetraEnd();
         it!=end; ++it)
    {
        // collect data on all DoF
        for (int v=0; v<10; ++v)
        {
            if (v<4)
                Numb[v]= GetGlobNum( it->GetVertex(v)->Unknowns(idx) );
            else
                Numb[v]= GetGlobNum( it->GetEdge(v-4)->Unknowns(idx) );
        }

        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            ReprTetraT t;
            for (int vert= 0; vert<4; ++vert)
            {
                const IdxT Nr= Numb[ data.Vertices[vert]];
                t[vert]= Nr;
            }
            TetraList_.push_back(t);
        }
    }
}

/// \brief Create mapping from DoF to global number
void FastMarchCL::CreateGlobNumb()
{
    Comment("Create global numbering\n", DebugParallelNumC);
    const IdxT numDoF= v_->RowIdx->NumUnknowns();
    const IdxT idx   = v_->RowIdx->GetIdx();
    IdxT numExclusiveDoF=0;
    allExclusive_.resize(ProcCL::Size());
    globNumb_.resize(0);            // if this list allready exists, delete and init new
    globNumb_.resize(numDoF,NoIdx);
    exclusive_.resize(0);           // if this list allready exists, delete and init new
    exclusive_.resize(numDoF, true);

    // If owning an exclusive DoF, put them into list, otherwise write NoIdx
    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(), end=MG_.GetTriangVertexEnd(); it!=end; ++it)
    {
        if (it->IsExclusive(PrioHasUnk)){
            globNumb_[it->Unknowns(idx)] = numExclusiveDoF++;
        }
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(), end=MG_.GetTriangEdgeEnd(); it!=end; ++it)
    {
        if (it->IsExclusive(PrioHasUnk)){
            globNumb_[it->Unknowns(idx)] = numExclusiveDoF++;
        }
    }
    locNumb_.resize(numExclusiveDoF, NoIdx);

    // Get offset
    ProcCL::Gather(numExclusiveDoF, Addr(allExclusive_), -1);
    allOffset_.resize(ProcCL::Size()+1);
    for (int i=0; i<ProcCL::Size(); ++i)
        allOffset_[i+1]  = allOffset_[i] + allExclusive_[i];
    globsize_=allOffset_[ProcCL::Size()];
    offset_=allOffset_[ProcCL::MyRank()];


    // Append offset and calculate not exclusive dof
    for (IdxT i=0; i<numDoF; ++i){
        if (globNumb_[i]!=NoIdx)
            globNumb_[i] += offset_;
        else
            exclusive_[i]=false;
    }

    // Collect global
    DDD_IFExchange(InterfaceCL<VertexCL>::GetIF(), sizeof(IdxT),
                   HandlerGlobDOFGatherVertexC,    HandlerGlobDOFScatterVertexC );
    DDD_IFExchange(InterfaceCL<EdgeCL>::GetIF(), sizeof(IdxT),
                   HandlerGlobDOFGatherEdgeC,    HandlerGlobDOFScatterEdgeC );


    // in debug mode, check if everything is right
#if DROPSDebugC&DebugParallelNumC
    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(), end=MG_.GetTriangVertexEnd(); it!=end; ++it)
        if (globNumb_[ it->Unknowns(idx) ]==NoIdx)
            it->DebugInfo(std::cerr);
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(), end=MG_.GetTriangEdgeEnd(); it!=end; ++it)
        if (globNumb_[ it->Unknowns(idx) ]==NoIdx)
            it->DebugInfo(std::cerr);

    for (IdxT i=0; i<numDoF; ++i)
        if (globNumb_[i]==NoIdx)
            throw DROPSErrCL("FastMarchCL::CreateGlobNumb: Not all local DoF numbers mapped on global DoF number");

    IdxT numExc=0;
    for (IdxT i=0; i<numDoF; ++i)
        if (IsExclusive(i))
            ++numExc;
    if (numExc!= allExclusive_[ProcCL::MyRank()])
        throw DROPSErrCL("FastMarchCL::CreateGlobNumb: number of exclusive DoF does not match!");

    IdxT maxEntry=*std::max_element(globNumb_.begin(), globNumb_.end());
    if (maxEntry>globsize_){
        std::cerr << "["<<ProcCL::MyRank()<<"] max entry is "<<maxEntry<<" and globalsize_ is "<<globsize_<<std::endl;
        throw DROPSErrCL("FastMarchCL::CreateGlobNumb: max entry in globNumb is bigger than global size");
    }
#endif
}

/// \brief Collect all information on master proc
void FastMarchCL::Collect()
{
    Comment("Collect\n", DebugParallelNumC);
    const int TypTag= 5001, CoordTag=5002 , DataTag=5003, NeighTag=5004;

    if (ProcCL::IamMaster())
    {
        IdxT allDoF=0;
        IdxT maxExclusiveDoF=0;

        for (int i=0; i<ProcCL::Size(); ++i)
        {
            allDoF          += allExclusive_[i];
            maxExclusiveDoF  = std::max(maxExclusiveDoF, allExclusive_[i]);
        }

        // allocate mem for algorithm
        GlobNeigh_.resize(allDoF);
        GlobalTyp_.resize(allDoF);
        GlobalCoord_.resize(allDoF);
        GlobalV_.resize(allDoF);

        VectorCL recvCoord(3*maxExclusiveDoF);

        ProcCL::StatusT stat;
        for (int p=0; p<ProcCL::Size(); ++p)
        {
            if (p!=ProcCL::Master()){
                // recieve types
                ProcCL::Probe(p, TypTag, stat);
                IdxT count = ProcCL::GetCount<byte>(stat);
                ProcCL::Recv(Addr(GlobalTyp_)+allOffset_[p], count, p, TypTag);

                // recieve Coord
                ProcCL::Recv(Addr(recvCoord), 3*count, p, CoordTag);
                for (IdxT i=0; i<count; ++i)
                    for (int j=0; j<3; ++j)
                       GlobalCoord_[i+allOffset_[p]][j]=recvCoord[3*i+j];

                // recieve Data
                ProcCL::Recv(Addr(GlobalV_)+allOffset_[p], count, p, DataTag);

                // recieve neighborhood
                ProcCL::Probe(p, NeighTag, stat);
                count = ProcCL::GetCount<IdxT>(stat);

                // allocate mem for recieving
                VectorBaseCL<IdxT> neighs(count);
                ProcCL::Recv(Addr(neighs), count, p, NeighTag);

                Assert(count%4==0, DROPSErrCL("FastMarchCL::Collect: Number of DoF for tetra representation is not a multiple of 4!"), DebugParallelNumC);
                for (IdxT i=0; i<count; )
                {
                    // Create Tetra
                    ReprTetraT t;
                    for (int vert=0; vert<4; ++vert)
                        t[vert]= neighs[i++];
                    // Put Tetra in neighborhood list
                    for (int vert=0; vert<4; ++vert)
                        GlobNeigh_[ t[vert] ].push_back(t);
                }
            }
            else
            {
                // Copy all data:
                IdxT j=0;
                for (IdxT i=0; i<size_; ++i)
                {
                    if (!IsExclusive(i)) continue;
                    locNumb_[j++]=i;

                    const IdxT pos=GetGlobNum(i);//allOffset_[p]+i;
                    // Typ
                    GlobalTyp_[pos]= Typ_[i];
                    // Coord
                    for (int j=0; j<3; ++j)
                        GlobalCoord_[pos][j]= Coord_[i][j];
                    // Data
                    GlobalV_[pos]= v_->Data[i];
                }
                // Copy all neigh-tetras
                for (VertexNeighT::const_iterator it(TetraList_.begin()), end(TetraList_.end()); it!=end; ++it)
                    for (int vert=0; vert<4; ++vert)
                        GlobNeigh_[ (*it)[vert] ].push_back(*it);


//                 for (IdxT i=0; i<TetraList_.size(); ++i)
//                     for (int vert=0; vert<4; ++vert)
//                         GlobNeigh_[ TetraList_[i][vert] ].push_back(TetraList_[i]);
            }
        }
    }
    else    // not master
    {
        const int me=ProcCL::MyRank();
        const IdxT myExclusive=allExclusive_[me];
        VectorBaseCL<ProcCL::RequestT> req(4);

        // Typ
        VectorBaseCL<byte> sendTypBuf(myExclusive);
        IdxT j=0;
        for (IdxT i=0; i<size_; ++i)
            if (IsExclusive(i)){
                sendTypBuf[j]= Typ_[i];
                locNumb_[j++]=i;                // rember the sequence of sending dofs
            }
        Assert(j==myExclusive, DROPSErrCL("FastMarchCL::Collect: Not the right number of types"), DebugParallelNumC);
        req[0]= ProcCL::Isend(Addr(sendTypBuf), myExclusive, ProcCL::Master(), TypTag);

        // Coord
        VectorBaseCL<double> sendCoordBuf(3*myExclusive);
        j=0;
        for (IdxT i=0; i<size_; ++i){
            if (IsExclusive(i)){
                for (int pos=0; pos<3; ++pos){
                    sendCoordBuf[3*j+pos]= Coord_[i][pos];
                }
                ++j;
            }
        }
        Assert(j==myExclusive, DROPSErrCL("FastMarchCL::Collect: Not the right number of coords"), DebugParallelNumC);
        req[1]= ProcCL::Isend(Addr(sendCoordBuf), 3*myExclusive, ProcCL::Master(), CoordTag);

        // Data
        VectorBaseCL<double> sendDataBuf(myExclusive);
        j=0;
        for (IdxT i=0; i<size_; ++i)
            if (IsExclusive(i))
                sendDataBuf[j++]= v_->Data[i];
        Assert(j==myExclusive, DROPSErrCL("FastMarchCL::Collect: Not the right number of datas"), DebugParallelNumC);
        req[2]= ProcCL::Isend(Addr(sendDataBuf), myExclusive, ProcCL::Master(), DataTag);

        // neighborhood
        IdxT numTetra=TetraList_.size();
        VectorBaseCL<IdxT> sendNeighBuf(4*numTetra);

        j=0;
        for (VertexNeighT::const_iterator it(TetraList_.begin()), end(TetraList_.end()); it!=end; ++it)
            for (int vert=0; vert<4; ++vert)
                sendNeighBuf[j++]= (*it)[vert];

//         for (IdxT i=0; i<numTetra; ++i)
//             for (int vert=0; vert<4; ++vert)
//                 sendNeighBuf[j++]= TetraList_[i][vert];
        Assert(j==4*numTetra, DROPSErrCL("FastMarchCL::Collect: Not the right number of tetras"), DebugParallelNumC);
        req[3]= ProcCL::Isend(Addr(sendNeighBuf), 4*numTetra, ProcCL::Master(), NeighTag);

        // Wait until all sends are complete, so that the buffer can deallocate
        ProcCL::WaitAll(4,Addr(req));
    }
}

void FastMarchCL::Distribute()
{
    Comment("Distribute\n", DebugParallelNumC);
    const int DataTag=5003;
    if (ProcCL::IamMaster())
    {
        VectorBaseCL<ProcCL::RequestT> req(ProcCL::Size()-1);
        int req_pos=0;
        for (int p=0; p<ProcCL::Size(); ++p)
        {
            if (p!=ProcCL::MyRank())
                req[req_pos++]= ProcCL::Isend(Addr(GlobalV_)+allOffset_[p], allExclusive_[p], p, DataTag);
        }
        // copy data
        for (IdxT i=0; i<size_; ++i)
            v_->Data[i] = GlobalV_[GetGlobNum(i)];
        ProcCL::WaitAll(ProcCL::Size()-1, Addr(req));
    }
    else
    {
        ProcCL::StatusT stat;
        ProcCL::Probe(ProcCL::Master(), DataTag, stat);
        __UNUSED__ IdxT mpi_count= ProcCL::GetCount<IdxT>(stat);
        const IdxT count=allExclusive_[ProcCL::MyRank()];
        Assert(mpi_count==count, DROPSErrCL("FastMarchCL::Distribute: Not enough data recieved"), DebugParallelNumC);
        VectorCL recvBuf(count);

        ProcCL::Recv(Addr(recvBuf), count, ProcCL::Master(), DataTag);
        for (IdxT i=0; i<count; ++i)
            v_->Data[locNumb_[i]] = recvBuf[i];
    }
}

void FastMarchCL::DistributeZeroLevel()
{
    VectorCL myCoordZeroLevel, myValueZeroLevel;
    InitZeroSet(myCoordZeroLevel, myValueZeroLevel);

    IdxT LocNumFinished=myValueZeroLevel.size(), GlobNumFinished=0;
    VectorBaseCL<IdxT> allFinished(ProcCL::Size());

    ProcCL::Gather(LocNumFinished, Addr(allFinished),-1);
    for (int p=0; p<ProcCL::Size(); ++p)
        GlobNumFinished+=allFinished[p];

    std::cout << " Euklidian FastMarching: global finished dofs " << GlobNumFinished
              << ", used memory " << (4*GlobNumFinished/1024) << " kB" << std::endl;

    CoordZeroLevel_.resize( 3*GlobNumFinished);
    ValueZeroLevel_.resize( GlobNumFinished);

    size_t pos1=0, pos2=0;
    for (int p=0; p<ProcCL::Size(); ++p)
    {
        if (p==ProcCL::MyRank())
        {
            ProcCL::Bcast(Addr(myCoordZeroLevel), 3*allFinished[p], p);
            ProcCL::Bcast(Addr(myValueZeroLevel),   allFinished[p], p);
            std::memcpy(&(CoordZeroLevel_[pos1]), Addr(myCoordZeroLevel), 3*allFinished[p]*sizeof(double));
            std::memcpy(&(ValueZeroLevel_[pos2]), Addr(myValueZeroLevel),   allFinished[p]*sizeof(double));
        }
        else
        {
            ProcCL::Bcast(&(CoordZeroLevel_[pos1]), 3*allFinished[p], p);
            ProcCL::Bcast(&(ValueZeroLevel_[pos2]),   allFinished[p], p);
        }
        pos1 += 3*allFinished[p];
        pos2 +=   allFinished[p];
    }
}

void FastMarchCL::CleanUp()
{
    Comment("Cleaning up\n", DebugParallelNumC);
    Assert(Close_.empty(), DROPSErrCL("FastMarchCL::CleanUp: Close set is not empty"), DebugParallelNumC);
    neigh_.resize(0);
    TetraList_.clear();
    GlobalV_.resize(0);
    allOffset_.resize(0);
    GlobNeigh_.resize(0);
    GlobalTyp_.resize(0);
    GlobalCoord_.resize(0);
    tmpv_.resize(0);
    Typ_.resize(0);
    tmpTyp_.resize(0);
    allExclusive_.resize(0);
    globNumb_.resize(0);
    locNumb_.resize(0);
    exclusive_.resize(0);
}
#endif      // end of parallel functions

// ===============================================================
//              variants for periodic boundaries
// ===============================================================
double FastMarchCL::CompValueProjPer( IdxT Nr, int num, const IdxT upd[3]) const
{
#ifdef _PAR
    throw DROPSErrCL("FastMarchCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif
    double val= 1e99;

    switch (num)
    {
        case 2: // Projektion auf Edge
        {
            const Point3DCL a= Coord_[upd[1]] - Coord_[upd[0]];
            const Point3DCL b= Coord_[  Nr  ] - Coord_[upd[0]];
            const double bary= inner_prod(a,b)/a.norm_sq();
            if (bary>=0 && bary<=1)
            {
                const Point3DCL lotfuss= (1-bary)*Coord_[upd[0]] + bary*Coord_[upd[1]];
                const double y= (1-bary)*v_->Data[Map(upd[0])] + bary*v_->Data[Map(upd[1])];
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
        break;

        case 3: // Projektion auf Face
        {
            const Point3DCL a= Coord_[upd[1]] - Coord_[upd[0]];
            const Point3DCL b= Coord_[upd[2]] - Coord_[upd[0]];
            const Point3DCL c= Coord_[  Nr  ] - Coord_[upd[0]];
            const double bary1= inner_prod(a,c)/a.norm_sq(),
                         bary2= inner_prod(b,c)/b.norm_sq();
            if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
            {
                const Point3DCL lotfuss= (1-bary1-bary2)*Coord_[upd[0]] + bary1*Coord_[upd[1]] + bary2*Coord_[upd[2]];
                const double y= (1-bary1-bary2)*v_->Data[Map(upd[0])] + bary1*v_->Data[Map(upd[1])] + bary2*v_->Data[Map(upd[2])];
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
    }

    return val;
}

void FastMarchCL::InitZeroPer( const BndDataCL<>& bnd, bool ModifyZero)
{
#ifdef _PAR
    throw DROPSErrCL("FastMarchCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif

    // Knoten an der Phasengrenze als Finished markieren
    // und Distanz zur Phasengrenze bestimmen (falls ModifyZero)
    const Uint idx= v_->RowIdx->GetIdx(),
               lvl= v_->GetLevel();
    int        sign[10];
    int        num_sign[3]; // - 0 +
    IdxT       Numb[10];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);
    IdxDescCL  augmIdx( P2_FE, bnd);
    augmIdx.CreateNumbering( lvl, MG_);
    const Uint augm_idx= augmIdx.GetIdx();

    // init Coord_, map_ and augmIdx
    std::vector<bool> ini( size_, false);
    IdxT k= 0;
    Coord_.resize( augmIdx.NumUnknowns());
    map_.resize( augmIdx.NumUnknowns() - size_);
    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(lvl), end=MG_.GetTriangVertexEnd(lvl);
        it!=end; ++it)
    {
        const IdxT Nr= it->Unknowns(idx);
        if (ini[Nr])
        {
            it->Unknowns(augm_idx)= size_ + k;
            map_[k++]= Nr;
        }
        else // touched for the first time
        {
            ini[Nr]= true;
            it->Unknowns(augm_idx)= Nr;
        }
        Coord_[it->Unknowns(augm_idx)]= it->GetCoord();
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(lvl), end=MG_.GetTriangEdgeEnd(lvl);
        it!=end; ++it)
    {
        const IdxT Nr= it->Unknowns(idx);
        if (ini[Nr])
        {
            it->Unknowns(augm_idx)= size_ + k;
            map_[k++]= Nr;
        }
        else // touched for the first time
        {
            ini[Nr]= true;
            it->Unknowns(augm_idx)= Nr;
        }
        Coord_[it->Unknowns(augm_idx)]= GetBaryCenter( *it);
    }
    ini.clear();

    // store copy of v_.Data in Old_
    Old_.resize( size_);
    Old_= v_->Data;

    VecDescCL oldv( *v_);
    LocalP2CL<> PhiLoc;

    neigh_.resize( size_);

    for (MultiGridCL::TriangTetraIteratorCL it=MG_.GetTriangTetraBegin(lvl), end=MG_.GetTriangTetraEnd(lvl);
        it!=end; ++it)
    {
        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            Numb[v]= v<4 ? it->GetVertex(v)->Unknowns(augm_idx)
                         : it->GetEdge(v-4)->Unknowns(augm_idx);
            const IdxT MapNr= Map(Numb[v]);
            sign[v]= std::abs(Old_[MapNr])<1e-8 ? 0 : (Old_[MapNr]>0 ? 1 : -1);
            if (sign[v]==0)
                Typ_[MapNr]= Finished;
        }
        PhiLoc.assign( *it, oldv, NoBndDataCL<>());

        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            // init num_sign, neigh_
            ReprTetraT t;
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
            {
                const Uint v= data.Vertices[vert];
                t[vert]= Numb[v];
                ++num_sign[ sign[v] + 1];
            }
            for (int vert= 0; vert<4; ++vert)
                neigh_[Map(t[vert])].push_back( t);

            const bool intersec= (num_sign[0]*num_sign[2]!=0); // Vorzeichenwechsel

            if (!intersec) continue;

            if (!ModifyZero)
            {
                for (int vert= 0; vert<4; ++vert)
                {
                    const IdxT MapNr= Map(Numb[data.Vertices[vert]]);
                    Typ_[MapNr]= Finished;
                    v_->Data[MapNr]= std::abs( Old_[MapNr]);
                }
                continue;
            }

            // ab hier gilt intersec && ModifyZero == true

            Point3DCL Schnitt[4];
            int num= 0;
            // Berechnung der Schnittpunkte mit den Kanten des Tetra
            for (int vert= 0; vert<4; ++vert)
                if (sign[data.Vertices[vert]]==0)
                    Schnitt[num++]= Coord_[Numb[data.Vertices[vert]]];

            for (int edge= 0; edge<6 && num<4; ++edge)
            {
                const Ubyte v1= data.Vertices[ VertOfEdge( edge, 0)],
                            v2= data.Vertices[ VertOfEdge( edge, 1)];
                if (sign[v1]*sign[v2] == -1) // Vorzeichenwechsel auf edge
                {
                    const IdxT Nr1= Numb[v1],
                               Nr2= Numb[v2];
                    const double bary= InterfacePatchCL::EdgeIntersection( v1,v2, PhiLoc);
                    Schnitt[num++]= (1-bary)*Coord_[Nr1] + bary*Coord_[Nr2];
                }
            }
/*
fil << "geom {OFF " << num << " 1 0\n";
for (int i=0; i<num; ++i)
{
    for (int j=0; j<3; ++j)
        fil << Schnitt[i][j] << ' ';
    fil << '\n';
}
if (num==3)
    fil << "3 0 1 2";
else
    fil << "4 0 1 3 2";
fil << "\n}\n";
*/
            if (num<3) throw DROPSErrCL("FastMarchCL::InitZero: intersection missing");

            for (int repeat=0; repeat<num-2; ++repeat)
            { // fuer num==4 (Schnitt ABDC ist viereckig)
              // zwei Dreiecke ABC + DBC betrachten
                if (repeat) Schnitt[0]= Schnitt[3];

                const Point3DCL a= Schnitt[1] - Schnitt[0],
                                b= Schnitt[2] - Schnitt[0];

                for (int vert=0; vert<4; ++vert)
                {
                    if (sign[data.Vertices[vert]]==0) continue;

                    const IdxT Nr= Numb[data.Vertices[vert]],
                               MapNr= Map(Nr);
                    const Point3DCL Crd= Coord_[Nr],
                                    c=   Crd - Schnitt[0];
                    double dist= std::min( c.norm(), (Crd-Schnitt[1]).norm());
                    dist= std::min( dist, (Crd-Schnitt[2]).norm());

                    const double bary1= inner_prod(a,c)/a.norm_sq(),
                                 bary2= inner_prod(b,c)/b.norm_sq();
                    if (bary1>=0 && bary2>=0 && bary1+bary2<=1)
                    {
                       const Point3DCL lotfuss= (1-bary1-bary2)*Schnitt[0] + bary1*Schnitt[1] + bary2*Schnitt[2];
                       dist= std::min( dist, (lotfuss - Crd).norm());
                    }

                    if (Typ_[MapNr] != Finished)
                    {
                        Typ_[MapNr]= Finished;
                        v_->Data[MapNr]= dist;
                    }
                    else
                        v_->Data[MapNr]= std::min( dist, v_->Data[MapNr]);
                }
            }
        }
    }
//fil << "}\n";
    // delete memory allocated for augmIdx
    DeleteNumbOnSimplex( augm_idx, MG_.GetAllVertexBegin(lvl), MG_.GetAllVertexEnd(lvl) );
    DeleteNumbOnSimplex( augm_idx, MG_.GetAllEdgeBegin(lvl), MG_.GetAllEdgeEnd(lvl) );
}

void FastMarchCL::InitClosePer()
{
#ifdef _PAR
    throw DROPSErrCL("FastMarchCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif
    // an Finished angrenzende Knoten mit Close markieren und dort v_ updaten
    const Uint idx= v_->RowIdx->GetIdx(),
               lvl= v_->GetLevel();

TimerCL tim;
tim.Start();
    std::set<IdxT> neighVerts;
    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(lvl), end=MG_.GetTriangVertexEnd(lvl);
        it!=end; ++it)
    {
        const IdxT Nr= it->Unknowns(idx);
        if (Typ_[Nr] == Finished)
        {
            for (Uint n=0; n<neigh_[Nr].size(); ++n)
                for (int j=0; j<4; ++j)
                    neighVerts.insert( neigh_[Nr][n][j]);
        }
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(lvl), end=MG_.GetTriangEdgeEnd(lvl);
        it!=end; ++it)
    {
        const IdxT Nr= it->Unknowns(idx);
        if (Typ_[Nr] == Finished)
        {
            for (Uint n=0; n<neigh_[Nr].size(); ++n)
                for (int j=0; j<4; ++j)
                    neighVerts.insert( neigh_[Nr][n][j]);
        }
    }

    for (std::set<IdxT>::const_iterator it= neighVerts.begin(), end= neighVerts.end();
        it!=end; ++it)
    { // update all neighboring verts, mark as Close
        UpdatePer( *it);
    }
}


void FastMarchCL::UpdatePer( const IdxT NrI)
{
#ifdef _PAR
    throw DROPSErrCL("FastMarchCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif

    const IdxT MapNrI= Map(NrI);
    // Update all vertices that are not Finished
    if (Typ_[MapNrI] == Finished) return;

    IdxT upd[3];
    double minval= Typ_[MapNrI]==Close ? v_->Data[MapNrI] : 1e99;

    for (Uint n=0; n<neigh_[MapNrI].size(); ++n)
    {
        int num= 0;
        for (int j=0; j<4; ++j)
        {
            const IdxT NrJ= neigh_[MapNrI][n][j];
            if (Typ_[Map(NrJ)] == Finished)
            {
                upd[num++]= NrJ;
                minval= std::min( minval, v_->Data[Map(NrJ)] + (Coord_[NrJ]-Coord_[NrI]).norm());
            }
        }

        minval= std::min( minval, CompValueProjPer( NrI, num, upd));
    }

    v_->Data[MapNrI]= minval;
    if (Typ_[MapNrI] != Close)
    {
        Close_.insert( MapNrI);
        Typ_[MapNrI]= Close;
    }
}


void FastMarchCL::ReparamPer( const BndDataCL<>& bnd, bool ModifyZero)
{
#ifdef _PAR
    throw DROPSErrCL("FastMarchCL: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif

    InitZeroPer( bnd, ModifyZero);
    InitClosePer();

    IdxT next;

    while ((next= FindTrial()) != NoIdx)
    {
        // remark: next < size_   =>   Map not needed for next
        Close_.erase( next);
        Typ_[next]= Finished;

        std::set<IdxT> neighVerts;
        for (Uint n=0; n<neigh_[next].size(); ++n)
        { // collect all neighboring verts in neighVerts
            for (Uint i=0; i<4; ++i)
                neighVerts.insert( neigh_[next][n][i]);
        }
        for (std::set<IdxT>::const_iterator it= neighVerts.begin(), end= neighVerts.end();
            it!=end; ++it)
        { // update all neighboring verts, mark as Close
            UpdatePer( *it);
        }
        neigh_[next].clear(); // will not be needed anymore
    }

    RestoreSigns();
}


} // end of namespace DROPS
