/// \file interfacePatch.cpp
/// \brief Computes 2D patches and 3D cuts of tetrahedra and interface
/// \author Sven Gross, Joerg Grande, Patrick Esser,Eva IGPM

#include "num/interfacePatch.h"

namespace DROPS
{

//*****************************************************************************
//                               InterfacePatchCL
//*****************************************************************************

const double InterfacePatchCL::approxZero_= 2.*std::numeric_limits<double>::epsilon();
const bool   InterfacePatchCL::LinearEdgeIntersection = true;
BaryCoordCL  InterfacePatchCL::AllEdgeBaryCenter_[10][10];
BaryCoordCL  InterfacePatchCL::BaryDoF_[10];

InterfacePatchCL::InterfacePatchCL()
  : RegRef_( GetRefRule( RegRefRuleC)), intersec_(0), ch_(-1)
{
    BaryDoF_[0][0]= BaryDoF_[1][1]= BaryDoF_[2][2]= BaryDoF_[3][3]= 1.;
    for (int edge=0; edge<6; ++edge)
        BaryDoF_[edge+4]= 0.5*(BaryDoF_[VertOfEdge(edge,0)] + BaryDoF_[VertOfEdge(edge,1)]);
    for (int i= 0; i < 10; ++i)
        for (int j= 0; j < 10; ++j)
            AllEdgeBaryCenter_[i][j]= BaryCenter( BaryDoF_[i], BaryDoF_[j]);
}

void InterfacePatchCL::Init( const TetraCL& t, const VecDescCL& ls, double translation)
{
    const Uint idx_ls= ls.RowIdx->GetIdx();
    ch_= -1;
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        Coord_[v]= v<4 ? t.GetVertex(v)->GetCoord() : GetBaryCenter( *t.GetEdge(v-4));
        PhiLoc_[v]= ls.Data[v<4 ? t.GetVertex(v)->Unknowns(idx_ls) : t.GetEdge(v-4)->Unknowns(idx_ls)] + translation;
        sign_[v]= Sign(PhiLoc_[v]);
    }
    barysubtetra_ = false;
}

void InterfacePatchCL::Init( const TetraCL& t, const LocalP2CL<double>& ls, double translation)
{
    ch_= -1;
    PhiLoc_= ls + translation;
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        Coord_[v]= v<4 ? t.GetVertex(v)->GetCoord() : GetBaryCenter( *t.GetEdge(v-4));
        sign_[v]= Sign(PhiLoc_[v]);
    }
    barysubtetra_ = false;
}

 SMatrixCL<3,4> GetCoordMatrix( const TetraCL& t)
{
   SMatrixCL<3,4> V;
   for (int i = 0; i<3; ++i)
       for(int j=0; j<4; ++j)
           V(i,j)= t.GetVertex(j)->GetCoord()[i];
   return V;

}

// Init for SubtetraT
void InterfacePatchCL::Init( const TetraCL& t, const SubTetraT& st, const LocalP2CL<double>& ls, double translation)
{
	st_ = st;
    ch_= -1;

    for (int v=0; v<10; ++v)
    {
    	BaryCoordCL tempBaryCoord_ = v<4 ? st[v] : BaryCenter(st[VertOfEdge(v-4,0)],st[VertOfEdge(v-4,1)]);
    	PhiLoc_[v] =ls( tempBaryCoord_) + translation;
    // collect data on all DoF
        Coord_[v]= GetCoordMatrix(t)* tempBaryCoord_;
        sign_[v]= Sign(PhiLoc_[v]);
    }
    barysubtetra_ = true;
}

// multiplication of SubTetraT with st
InterfacePatchCL::SubTetraT InterfacePatchCL::MultiplySubTetra(const InterfacePatchCL::SubTetraT & Tetrak_ )
{
 SubTetraT TetrakBary_;
 for (int i = 0; i <4; ++i)
  {
	 for (int j = 0; j <4; ++j)
	 {
		 {
			 TetrakBary_[i] += (st_[j]*Tetrak_[i][j]);
		 }
	 }
  }
 return TetrakBary_;
}

BaryCoordCL InterfacePatchCL::MultiplyBaryCoord(const BaryCoordCL& Tetrak_ )
{
 BaryCoordCL TetrakBary_;

	 for (int j = 0; j <4; ++j)
	 {
	     TetrakBary_ += (st_[j]*Tetrak_[j]);
	 }

 return TetrakBary_;
}
bool InterfacePatchCL::ComputeForChild( Uint ch)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
    ch_= ch;
    num_sign_[0]= num_sign_[1]= num_sign_[2]= 0;
    for (int vert= 0; vert<4; ++vert)
        ++num_sign_[ sign_[data.Vertices[vert]] + 1];

    intersec_= 0;
    if (num_sign_[0]*num_sign_[2]==0 && num_sign_[1]<3) {// no change of sign on child
        numchildtriangles_= 0;
        return false;
    }
    if (num_sign_[1]==4)
    {
        std::cerr << "WARNING: InterfacePatchCL: found 3-dim. zero level set, grid is too coarse!" << std::endl;
        numchildtriangles_= 0;
        return false;
    }

    // erst werden die Nullknoten in PQRS gespeichert...
    for (int vert= 0; vert<4; ++vert)
    {
        const int v= data.Vertices[vert];
        if (sign_[v]==0)
        {
            Bary_[intersec_]= BaryDoF_[v];
            PQRS_[intersec_++]= Coord_[v];
        }
    }
    // ...dann die echten Schnittpunkte auf den Kanten mit Vorzeichenwechsel
    for (int edge= 0; edge<6; ++edge)
    {
        const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                  v1= data.Vertices[ VertOfEdge( edge, 1)];
        if (sign_[v0]*sign_[v1]<0) // different sign -> 0-level intersects this edge
        {
            const double lambda= EdgeIntersection( v0,  v1, PhiLoc_);
            Bary_[intersec_]= (1-lambda)*BaryDoF_[v0] + lambda * BaryDoF_[v1];
            // bary-coords of tetra, not of subtetra!
            PQRS_[intersec_++]= (1-lambda) * Coord_[v0] + lambda * Coord_[v1];
        }
    }
    if (intersec_<3) { // Nullstellenmenge vom Mass 0!
        numchildtriangles_= 0;
        return false;
    }

    SMatrixCL<3,2> A;    // A = [ Q-P | R-P ]
    A(0,0)= PQRS_[1][0]-PQRS_[0][0];    A(0,1)= PQRS_[2][0]-PQRS_[0][0];
    A(1,0)= PQRS_[1][1]-PQRS_[0][1];    A(1,1)= PQRS_[2][1]-PQRS_[0][1];
    A(2,0)= PQRS_[1][2]-PQRS_[0][2];    A(2,1)= PQRS_[2][2]-PQRS_[0][2];
    SMatrixCL<2,2> ATA;
    ATA(0,0)=           A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    ATA(0,1)= ATA(1,0)= A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
    ATA(1,1)=           A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
    const double detATA= ATA(0,0)*ATA(1,1) - ATA(1,0)*ATA(1,0);
    sqrtDetATA_= std::sqrt( detATA);

    Point2DCL AT_i, tmp;
    for (int i=0; i<3; ++i)
    {
        // berechne B = A * (ATA)^-1 * AT
        AT_i[0]= A(i,0); AT_i[1]= A(i,1);
        Solve2x2( detATA, ATA, tmp, AT_i);
        B_[i]= A*tmp;
    }

    if (intersec_==4) // 4 intersections --> a+b != 1
    { // berechne a, b
        // Loese (Q-P)a + (R-P)b = S-P  --> lin. AGP, loese ATA * [a,b]T = AT(S-P)
        Point3DCL PS= PQRS_[3] - PQRS_[0];
        tmp[0]= A(0,0)*PS[0] + A(1,0)*PS[1] + A(2,0)*PS[2];
        tmp[1]= A(0,1)*PS[0] + A(1,1)*PS[1] + A(2,1)*PS[2];
        Solve2x2( detATA, ATA, ab_, tmp);
        //if (ab_[0]<0 || ab_[1]<0)
        //    std::cout<<"LevelsetP2CL::AccumulateBndIntegral: a or b negative"<<std::endl;
        // a,b>=0 muss erfuellt sein, da wegen edge+oppEdge==5 die Punkte P und S sich automatisch gegenueber liegen muessten...
        numchildtriangles_= 2;
    }
    else
        numchildtriangles_= 1;

    if (EqualToFace()) // interface is shared by two tetras
        sqrtDetATA_/= 2;

    // if Init for Sub TetraT has been used, coordinates must be transformed
    if (barysubtetra_ == true)
    {
        for (int k=0 ; k<intersec_; ++k)
        {
            Bary_[k] = MultiplyBaryCoord(Bary_[k]);
        }
    }
    return true; // computed patch of child;
}


Point3DCL InterfacePatchCL::GetNormal() const
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch_]);
    SMatrixCL<3,4> p1grad;
    double det; // dummy
    Point3DCL pt[4];
    SVectorCL<4> ls;
    for (int v= 0; v < 4; ++v) {
        pt[v]= Coord_ [data.Vertices[v]];
        ls[v]= PhiLoc_[data.Vertices[v]];
    }
    P1DiscCL::GetGradients( p1grad, det, pt);
    const Point3DCL n( p1grad*ls);
    return n/n.norm();
}


bool InterfacePatchCL::ComputeCutForChild( Uint ch)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
    ch_= ch;
    num_sign_[0]= num_sign_[1]= num_sign_[2]= 0;
    for (int vert= 0; vert<4; ++vert)
        ++num_sign_[ sign_[data.Vertices[vert]] + 1];

    intersec_= innersec_= 0;
    if (num_sign_[0]*num_sign_[2]==0 && num_sign_[1]<3) // no change of sign on child and no patch on a face
        return false;
    if (num_sign_[1]==4)
    {
        std::cerr << "WARNING: InterfacePatchCL: found 3-dim. zero level set, grid is too coarse!" << std::endl;
        return false;
    }

    // erst werden die Nullknoten in PQRS gespeichert...
    for (int vert= 0; vert<4; ++vert)
    {
        const int v= data.Vertices[vert];
        if (sign_[v]==0)
        {
            Bary_[intersec_]= BaryDoF_[v];
            Edge_[intersec_++]= -1;
        }
    }
    // ...dann die echten Schnittpunkte auf den Kanten mit Vorzeichenwechsel
    for (int edge= 0; edge<6; ++edge)
    {
        const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                  v1= data.Vertices[ VertOfEdge( edge, 1)];
        if (sign_[v0]*sign_[v1]<0) // different sign -> 0-level intersects this edge
        {
            const double lambda= EdgeIntersection( v0,  v1, PhiLoc_);
            Bary_[intersec_]= (1-lambda)*BaryDoF_[v0] + lambda * BaryDoF_[v1];
            Edge_[intersec_++]= edge;
            innersec_++;
        }
    }
    if (intersec_<3) return false; // Nullstellenmenge vom Mass 0!

    return true; // computed cut of child;
}

void InterfacePatchCL::DebugInfo( std::ostream& os, bool InfoForChild) const
{
    if (InfoForChild)
    {
        const ChildDataCL data= GetChildData( RegRef_.Children[ch_]);
        os << "Patch on child " << ch_ << " with " << GetNumPoints() << " intersections:\n";
        for (Uint i=0; i<GetNumPoints(); ++i)
            os << "\tP" << i << " = " << GetPoint(i) << '\n';
        os << "Signs of verts: ";
        for (int i=0; i<4; ++i)
            os << GetSign(data.Vertices[i]) << " ";
        os << std::endl;
    }
    else
    {
        os << "Signs: ";
        for (int i=0; i<10; ++i)
        {
            os << GetSign(i) << " ";
        }
        os << std::endl;
    }
}

void InterfacePatchCL::WriteGeom( std::ostream& os) const
{
    os << "geom {OFF " << intersec_ << " 1 0\n";
    for (int i=0; i<intersec_; ++i)
    {
        for (int j=0; j<3; ++j)
            os << PQRS_[i][j] << ' ';
        os << '\n';
    }
    if (IsQuadrilateral())
        os << "4 0 1 3 2";
    else
        os << "3 0 1 2";
    os << "\n}\n";
}

void InterfacePatchCL::InsertSubTetra(SubTetraT& BaryCoords, bool pos)
{
    if (pos)
        posTetras.push_back(BaryCoords);
    else
        negTetras.push_back(BaryCoords);
}

void InterfacePatchCL::ComputeSubTets()
{
    posTetras.clear();
    negTetras.clear();

    ChildDataCL data;
    SubTetraT BaryCoords;

    if (!Intersects())
    {
        for (Uint i=0; i<4; ++i)
            BaryCoords[i] = BaryDoF_[i];
        InsertSubTetra( BaryCoords, sign_[0]==1);
        return;
    }

    // Schleife ueber die Kinder
    for (Uint ch=0; ch < 8 && Intersects(); ++ch)
    {
        //std::cout << "Kind " << ch << "\n \n";
        data  = GetChildData (RegRef_.Children[ch]);

        //cuts = Child + empty set
        if (!ComputeCutForChild(ch) || (intersec_==3 && innersec_==0))
        {
            //std::cout << "cuts = Child + empty set: " << num_sign_[2] <<"\n";
            for (Uint i=0; i<4; ++i)
                BaryCoords[i]=BaryDoF_[data.Vertices[i]];
            InsertSubTetra( BaryCoords, num_sign_[2]>0);
        }
        if (intersec_==3)
            switch (innersec_)
            {
                case 1 : // cuts = Tetra + Tetra
                {
                    //std::cout << "cuts = Tetra + Tetra"<<"\n";
                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];
                    for (int i=0; i<4; ++i)
                    {
                        if (sign_[data.Vertices[i]]==0) continue;
                        BaryCoords[3]= BaryDoF_[data.Vertices[i]];
                        InsertSubTetra( BaryCoords, sign_[data.Vertices[i]]==1);
                    }
                } break;
                case 2 : // cuts = Tetra + Pyramide mit 4eckiger Grundflaeche
                {
                    // tetra
                    //std::cout << "cuts = Tetra + Pyramide"<<"\n";
                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];

                    int vertB= -1;
                    const int signB= num_sign_[0]==1 ? -1 : 1; // sign of vert B occurs only once
                    for (int i=0; vertB==-1 && i<4; ++i)
                        if (sign_[data.Vertices[i]]==signB) vertB= i;
                    BaryCoords[3]= BaryDoF_[data.Vertices[vertB]];
                    InsertSubTetra( BaryCoords, signB==1);

                    // pyramid = 2 tetras: ACDP + APQD
                    //                                     connectivity:     P-------Q
                    //                                                       | \   / |
                    //                                                       |   A   |
                    //                                                       | /   \ |
                    //                                                       C-------D
                    //A,C,D Ecken des Kindtetraeders
                    int z=0;
                    for (int i=0; i<4; ++i)
                        if (i!=vertB) BaryCoords[z++]= BaryDoF_[data.Vertices[i]];
                    BaryCoords[3]= Bary_[1];
                    InsertSubTetra( BaryCoords, signB!=1);

                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];
                    int vertD=-1;
                    for (int i=0; vertD==-1 && i<4; ++i)
                    {
                        if (sign_[data.Vertices[i]]==-signB &&
                            (i==VertOfEdge(Edge_[2],0) ||
                             i==VertOfEdge(Edge_[2],1))  )
                             vertD= i;
                    }
                    BaryCoords[3]= BaryDoF_[data.Vertices[vertD]];
                    InsertSubTetra( BaryCoords, signB!=1);
                } break;
                case 3 : // cuts = Tetra + Tetra-Stumpf
                {
                    //std::cout << "cuts = Tetra + Tetra-Stumpf\n";
                    int vertA= -1;  // cut-Tetra = APQR
                    const int signA= num_sign_[0]==1 ? -1 : 1; // sign of vert A occurs only once
                    for (int i=0; vertA==-1 && i<4; ++i)
                        if (sign_[data.Vertices[i]]==signA) vertA= i;
                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];
                    BaryCoords[3]= BaryDoF_[data.Vertices[vertA]];

                    InsertSubTetra( BaryCoords, signA==1);
                    //                                     connectivity:     R---------D
                    //                                                      /|        /|
                    //                                                     Q-|-------C |
                    //                                                      \|        \|
                    //                                                       P---------B
                    //B,C,D Ecken des Kindtetraeders
                    int vertBCD[3];
                    for (Uint i=0; i<3; ++i)
                        vertBCD[i]= (vertA==VertOfEdge(Edge_[i],0)) ? VertOfEdge(Edge_[i],1) : VertOfEdge(Edge_[i],0);

                    // QCPR
                    BaryCoords[0] = Bary_[1];
                    BaryCoords[1] = BaryDoF_[data.Vertices[vertBCD[1]]];
                    BaryCoords[2] = Bary_[0];
                    BaryCoords[3] = Bary_[2];
                    if ( signA!=1 && Edge_[1]!=-1)
                        posTetras.push_back(BaryCoords);
                    else
                        if (Edge_[1]!=-1) negTetras.push_back(BaryCoords);

                    // BCPR
                    BaryCoords[0] = BaryDoF_[data.Vertices[vertBCD[0]]];
                    if ( signA!=1 && Edge_[0]!=-1)
                        posTetras.push_back(BaryCoords);
                    else
                        if (Edge_[0]!=-1) negTetras.push_back(BaryCoords);

                    // BCDR
                    BaryCoords[2] = BaryDoF_[data.Vertices[vertBCD[2]]];
                    if ( signA!=1 && Edge_[2]!=-1)
                        posTetras.push_back(BaryCoords);
                    else
                        if (Edge_[2]!=-1) negTetras.push_back(BaryCoords);
                } break;
            }
        if (intersec_==4)
        {   // cuts = 2x 5-Flaechner mit 6 Knoten. connectivity:     R---------S
            //                                                      /|        /|
            //                                                     A-|-------B |
            //                                                      \|        \|
            //                                                       P---------Q
            //std::cout << "cuts = 2x 5-Flaechner\n";
            int vertAB[2]; // cut mit VZ==part = ABPQRS

            //erst werden die "negativen" Kinder in die Liste hinzugefuegt, dann die "positiven"
            for (int signAB = -1; signAB<=1; signAB+=2) {
                for (int i=0, k=0; i<4 && k<2; ++i)
                    if (sign_[data.Vertices[i]]==signAB) vertAB[k++]= i;
                // connectivity AP automatisch erfuellt, check for connectivity AR
                const bool AR= vertAB[0]==VertOfEdge(Edge_[2],0) || vertAB[0]==VertOfEdge(Edge_[2],1);
                // Integriere ueber Tetras ABPR, QBPR, QBSR    (bzw. mit vertauschten Rollen von Q/R)
                // ABPR    (bzw. ABPQ)
                BaryCoords[0]= BaryDoF_[data.Vertices[vertAB[0]]];
                BaryCoords[1]= BaryDoF_[data.Vertices[vertAB[1]]];
                BaryCoords[2]= Bary_[0];    BaryCoords[3]= Bary_[AR ? 2 : 1];
                InsertSubTetra( BaryCoords, signAB!=-1);
                // QBPR    (bzw. RBPQ)
                BaryCoords[0]=Bary_[AR ? 1 : 2];
                InsertSubTetra( BaryCoords, signAB!=-1);
                // QBSR    (bzw. RBSQ)
                BaryCoords[2]=Bary_[3];
                InsertSubTetra( BaryCoords, signAB!=-1);
            }
        } //intersec_==4 Ende
    } //Ende der Schleife ueber die Kinder

    // if Init for SubTetraT has bee unsed, coordinates must be transformed
     if (barysubtetra_ == true)
     {
         for (Uint k=0 ; k<posTetras.size(); ++k)
         {
             posTetras[k] = MultiplySubTetra(posTetras[k]);
         }
         for (Uint k=0 ; k<negTetras.size(); ++k)
         {
             negTetras[k] = MultiplySubTetra(negTetras[k]);
         }
     }
}

} // end of namespace DROPS

