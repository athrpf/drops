#include "levelset/fastmarch.h"
#include <fstream>

namespace DROPS
{

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
                const double y= (1-bary)*v_.Data[upd[0]] + bary*v_.Data[upd[1]];
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
                const double y= (1-bary1-bary2)*v_.Data[upd[0]] + bary1*v_.Data[upd[1]] + bary2*v_.Data[upd[2]];
                val= y + (lotfuss - Coord_[Nr]).norm();
            }
        }
    }
    
    return val;
}

void FastMarchCL::InitZero( bool ModifyZero)
{
    // Knoten an der Phasengrenze als Finished markieren
    // und Distanz zur Phasengrenze bestimmen (falls ModifyZero)
    const IdxT num_unks= v_.RowIdx->NumUnknowns;
    const Uint idx= v_.RowIdx->GetIdx();
    int        sign[10];
    int        num_sign[3]; // - 0 + 
    IdxT       Numb[10];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    // init Coord_
    Coord_.resize( num_unks);
    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(), end=MG_.GetTriangVertexEnd();
        it!=end; ++it)
        Coord_[it->Unknowns(idx)]= it->GetCoord();
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(), end=MG_.GetTriangEdgeEnd();
        it!=end; ++it)
        Coord_[it->Unknowns(idx)]= GetBaryCenter( *it);
    
    // store copy of v_.Data in Old_
    Old_.resize( num_unks);
    Old_= v_.Data;
    
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
                    v_.Data[Nr]= std::abs( Old_[Nr]);
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
                    const double bary= Old_[Nr1]/(Old_[Nr1]-Old_[Nr2]);
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
                        v_.Data[Nr]= dist;
                    }
                    else
                        v_.Data[Nr]= std::min( dist, v_.Data[Nr]);
                }
            }
        }
    }
//fil << "}\n";    
}

void FastMarchCL::InitClose()
{
    // an Finished angrenzende Knoten mit Close markieren und dort v_ updaten
    const IdxT num_unks= v_.RowIdx->NumUnknowns;
    const Uint idx=      v_.RowIdx->GetIdx();

    neigh_.resize( num_unks);

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


IdxT FastMarchCL::FindTrial() const
{
    double min= 1e99;
    IdxT min_idx= NoIdx;
    
    for (std::set<IdxT>::const_iterator it= Close_.begin(), end= Close_.end(); it!=end; ++it)
    {
        if (v_.Data[*it]<=min)
        {
            min= v_.Data[*it];
            min_idx= *it;
        }
    }
    return min_idx;
}

void FastMarchCL::Update( const IdxT NrI)
{
    // Update all vertices that are not Finished
    if (Typ_[NrI] == Finished) return;

    IdxT upd[3];
    double minval= Typ_[NrI]==Close ? v_.Data[NrI] : 1e99;

    for (Uint n=0; n<neigh_[NrI].size(); ++n)
    {
        int num= 0;
        for (int j=0; j<4; ++j)
        {
            const IdxT NrJ= neigh_[NrI][n][j];
            if (Typ_[NrJ] == Finished)
            {
                upd[num++]= NrJ;
                minval= std::min( minval, v_.Data[NrJ] + (Coord_[NrJ]-Coord_[NrI]).norm());
            }
        }

        minval= std::min( minval, CompValueProj( NrI, num, upd));
    }

    v_.Data[NrI]= minval;
    if (Typ_[NrI] != Close)
    {
        Close_.insert( NrI);
        Typ_[NrI]= Close;
    }
}


void FastMarchCL::Reparam( bool ModifyZero)
{
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
}

void FastMarchCL::RestoreSigns()
{ // restore signs of v_
    for (IdxT i=0, N= Old_.size(); i<N; ++i)
        if (Old_[i]<0)
            v_.Data[i]*= -1;
}

} // end of namespace DROPS
