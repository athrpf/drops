/// \file interfacePatch.cpp
/// \brief Computes 2D patches and 3D cuts of tetrahedra and interface
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Eva Loch; SC RWTH Aachen:

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

#include "num/interfacePatch.h"

namespace DROPS
{

//*****************************************************************************
//                               InterfacePatchCL
//*****************************************************************************

const double InterfacePatchCL::approxZero_= 2.*std::numeric_limits<double>::epsilon();
const bool   InterfacePatchCL::LinearEdgeIntersection = true;
BaryCoordCL   InterfacePatchCL::AllEdgeBaryCenter_[10][10];
BaryCoordCL   InterfacePatchCL::BaryDoF_[10];

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

void InterfacePatchCL::Init( const TetraCL& t, const VecDescCL& ls, const BndDataCL<>& lsetbnd, double translation)
{
    LocalP2CL<double> locphi(t, ls, lsetbnd);
    this->Init(t, locphi, translation);
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
   SMatrixCL<3,4> V(Uninitialized);
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

bool InterfacePatchCL::IntersectsChild( Uint ch) const
{
    const ChildDataCL& data= GetChildData ( ch);
    for (int i= 1; i < 4; ++i)
        if (sign_[data.Vertices[0]] != sign_[data.Vertices[i]])
            return true;
    return false;
}

int InterfacePatchCL::GetNumIntersectedSubTetras() const
{
    int intersectedSubTetra= 0;
    for (int ch=0; ch<8; ++ch) {
        if (IntersectsChild( ch))
            intersectedSubTetra++;
    }
    return intersectedSubTetra;
}

bool InterfacePatchCL::ComputeVerticesOfCut( Uint ch, bool compute_PQRS)
{
    const ChildDataCL& data= GetChildData( ch);
    ch_= ch;
    // Compute the sign of the levelset-function in the vertices of ch_.
    num_sign_[0]= num_sign_[1]= num_sign_[2]= 0;
    for (int vert= 0; vert < 4; ++vert)
        ++num_sign_[ sign_[data.Vertices[vert]] + 1];

    intersec_= innersec_= 0;
    // Return, if there is no change of sign on child and no patch on a face.
    if (num_sign_[0]*num_sign_[2]==0 && num_sign_[1]<3){
        return false;
    }
    // Warn in case of 3D-interface.
    if (num_sign_[1] == 4)
        throw DROPSErrCL("WARNING: InterfacePatchCL: found 3-dim. zero level set, grid is too coarse!");
    // first come the zero-vertices of the child ch_...
    for (int vert= 0; vert<4; ++vert)
    {
        const int v= data.Vertices[vert];
        if (sign_[v] == 0)
        {
            Bary_[intersec_]= BaryDoF_[v];
            if (compute_PQRS)
            	PQRS_[intersec_]= Coord_[v];
            Edge_[intersec_++]= -1;
        }
    }
    // ...then the real intersections on the edges with a change of sign.
    for (int edge= 0; edge<6; ++edge)
    {
        const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                  v1= data.Vertices[ VertOfEdge( edge, 1)];
        if (sign_[v0]*sign_[v1]<0) // different sign -> 0-level intersects this edge
        {
            const double lambda= EdgeIntersection( v0,  v1, PhiLoc_);
            Bary_[intersec_]= (1.-lambda)*BaryDoF_[v0] + lambda * BaryDoF_[v1];
            if (compute_PQRS)
            	PQRS_[intersec_]= (1.-lambda) * Coord_[v0] + lambda * Coord_[v1];
            Edge_[intersec_++]= edge;
            innersec_++;
        }
    }

    // zero-level of measure 0.
    if (intersec_<3){
    	return false;
	}
    // computed cut of child;
    return true;
}

void InterfacePatchCL::DebugInfo( std::ostream& os, bool InfoForChild) const
{
    if (InfoForChild)
    {
        const ChildDataCL data= GetChildData( ch_);
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

// multiplication of SubTetraT with st
InterfacePatchCL::SubTetraT InterfaceTetraCL::TransformToSubTetra(const InterfacePatchCL::SubTetraT& tetra)
{
    SubTetraT ret;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            ret[i] += st_[j]*tetra[i][j];
    return ret;
}

bool InterfaceTetraCL::ComputeCutForChild( Uint ch)
{
	return ComputeVerticesOfCut( ch);
}

void InterfaceTetraCL::InsertSubTetra(SubTetraT& BaryCoords, bool pos, Uint child)
{
    // if Init for SubTetraT has been used, coordinates must be transformed
    if (barysubtetra_ == true)
        BaryCoords= TransformToSubTetra(BaryCoords);

    if (pos) {
        posTetras.push_back(BaryCoords);
        posChildIdx.push_back(child);
    } else {
        negTetras.push_back(BaryCoords);
        negChildIdx.push_back(child);
    }
}

void InterfaceTetraCL::ComputeSubTets( Uint ch, bool clearTetras)
{
    if (clearTetras) {
        posTetras.clear();
        negTetras.clear();
    }

    //std::cout << "Kind " << ch << "\n\n";
    const ChildDataCL& data= GetChildData( ch);
    SubTetraT BaryCoords;

    //cuts = Child + empty set
    if (!ComputeCutForChild(ch) || (intersec_==3 && innersec_==0))
    {
        //std::cout << "cuts = Child + empty set: " << num_sign_[2] << "\n";
        for (Uint i=0; i<4; ++i)
            BaryCoords[i]=BaryDoF_[data.Vertices[i]];
        InsertSubTetra( BaryCoords, num_sign_[2]>0, ch);
        return;
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
                    InsertSubTetra( BaryCoords, sign_[data.Vertices[i]]==1, ch);
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
                InsertSubTetra( BaryCoords, signB==1, ch);

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
                InsertSubTetra( BaryCoords, signB!=1, ch);

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
                InsertSubTetra( BaryCoords, signB!=1, ch);
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

                InsertSubTetra( BaryCoords, signA==1, ch);
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
                    InsertSubTetra(BaryCoords, true, ch);
                else if (Edge_[1]!=-1)
                    InsertSubTetra(BaryCoords, false, ch);

                // BCPR
                BaryCoords[0] = BaryDoF_[data.Vertices[vertBCD[0]]];
                if ( signA!=1 && Edge_[0]!=-1)
                    InsertSubTetra(BaryCoords, true, ch);
                else if (Edge_[0]!=-1)
                    InsertSubTetra(BaryCoords, false, ch);

                // BCDR
                BaryCoords[2] = BaryDoF_[data.Vertices[vertBCD[2]]];
                if ( signA!=1 && Edge_[2]!=-1)
                    InsertSubTetra(BaryCoords, true, ch);
                else if (Edge_[2]!=-1)
                    InsertSubTetra(BaryCoords, false, ch);
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
        for (int signAB = -1; signAB<=1; signAB+=2) { int k = 0;
            for (int i=0; i<4 && k<2; ++i)
                if (sign_[data.Vertices[i]]==signAB) vertAB[k++]= i;
            // connectivity AP automatisch erfuellt, check for connectivity AR
            const bool AR= vertAB[0]==VertOfEdge(Edge_[2],0) || vertAB[0]==VertOfEdge(Edge_[2],1);
            // Integriere ueber Tetras ABPR, QBPR, QBSR    (bzw. mit vertauschten Rollen von Q/R)
            // ABPR    (bzw. ABPQ)
            BaryCoords[0]= BaryDoF_[data.Vertices[vertAB[0]]];
            BaryCoords[1]= BaryDoF_[data.Vertices[vertAB[1]]];
            BaryCoords[2]= Bary_[0];    BaryCoords[3]= Bary_[AR ? 2 : 1];
            InsertSubTetra( BaryCoords, signAB!=-1, ch);
            // QBPR    (bzw. RBPQ)
            BaryCoords[0]=Bary_[AR ? 1 : 2];
            InsertSubTetra( BaryCoords, signAB!=-1, ch);
            // QBSR    (bzw. RBSQ)
            BaryCoords[2]=Bary_[3];
            InsertSubTetra( BaryCoords, signAB!=-1, ch);
        }
    } //intersec_==4 Ende
}

void InterfaceTetraCL::ComputeSubTets()
{
    posTetras.clear();
    negTetras.clear();

    if (!Intersects())
    {
        SubTetraT BaryCoords;
        for (Uint i=0; i<4; ++i)
            BaryCoords[i] = BaryDoF_[i];
        InsertSubTetra( BaryCoords, sign_[0]==1, 8 /*the tetra itself as child*/);
        return;
    }

    // Schleife ueber die Kinder
    for (Uint ch= 0; ch < 8; ++ch)
        ComputeSubTets( ch, /*clearTetras=*/ false);
}

BaryCoordCL InterfaceTriangleCL::TransformToSubTetra(const BaryCoordCL& b)
{
    BaryCoordCL ret;
    for (int j = 0; j <4; ++j)
        ret += st_[j]*b[j];
    return ret;
}

bool InterfaceTriangleCL::ComputeForChild( Uint ch)
{
    const bool iscut= ComputeVerticesOfCut( ch, /*compute_PQRS*/ true);
    if (!iscut) { // no change of sign on child
        numtriangles_= 0;
        return false;
    }

    // Compute DetA_.
    SMatrixCL<3,2> A( Uninitialized); // A = [ Q-P | R-P ]
    A.col( 0, PQRS_[1] - PQRS_[0]);
    A.col( 1, PQRS_[2] - PQRS_[0]);
    QRDecompCL<3, 2> qr( A);
    DetA_= std::fabs( qr.Determinant_R());

    // Compute B = A * (ATA)^-1 * AT.
    std::memset( B_, 0, 9*sizeof( double));
    B_[0][0]= B_[1][1]= B_[2][2]= 1.;
    qr.Solve( 3, B_); // The first two components of each column are the least-squares-solution, the last is the residual.
    for (int i= 0; i < 3; ++i)
        B_[i]= A.col( 0)*B_[i][0] + A.col( 1)*B_[i][1];

    if (intersec_ == 4) // 4 intersections --> a+b != 1
    {   // Solve |(Q-P)a + (R-P)b - (S-P)| --> min; lin. least-squares problem (consistent in exact arithmetic!)

        Point3DCL PS= PQRS_[3] - PQRS_[0];
        qr.Solve( PS);
        ab_[0]= PS[0]; ab_[1]= PS[1];
        numtriangles_= 2;
    }
    else
        numtriangles_= 1;

    if (EqualToFace()) // interface is shared by two tetras
        DetA_/= 2.;

    // if Init for Sub TetraT has been used, coordinates must be transformed
    if (barysubtetra_ == true)
        for (int k=0 ; k<intersec_; ++k)
            Bary_[k] = TransformToSubTetra(Bary_[k]);
    return true; // computed patch of child;
}

Point3DCL InterfaceTriangleCL::GetNormal() const
{
    const ChildDataCL data= GetChildData( ch_);
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

} // end of namespace DROPS

