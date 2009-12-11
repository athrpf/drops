/// \file maketopo.cpp
/// \brief  computes the irregular refinement rules
/// \author LNM RWTH Aachen: Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

///  uses topo.cpp.in to compute topo.cpp

#include "geom/topo.h"
#include "geom/triang.cpp"
#include <algorithm>
#include <vector>
#include <set>
#include <utility>
#include <fstream>
#include <iostream>

using namespace DROPS;

typedef std::vector<Uint>  TetCL;
typedef std::vector<Uint>  FacCL;
typedef std::vector<TetCL> TetSeqCL;
typedef std::vector<FacCL> FacSeqCL;

///\brief Error class for generating topo.cpp
class TopoErrCL
{
  private:
    std::string ErrMesg_;
  public:
    TopoErrCL() : ErrMesg_("") {}
    TopoErrCL(const std::string& mesg) : ErrMesg_(mesg) {}
    ~TopoErrCL() {}

    std::ostream& what  (std::ostream& out) const{
        out << ErrMesg_ << std::endl;
        return out;
    }

    void handle() const{
        what(std::cout);
        std::cout.flush();
        std::abort();
    }
};

struct FaceRefRuleBaseCL
{
    Uint ChildNum;
    int  Children[4];
};

struct RawTriangCL
{
    std::vector<FacSeqCL> faces;
    std::vector<TetCL>    children;

    bool
    IsOuterFace(const FacCL&) const;
    void
    Print(std::ostream&) const;
};


FaceRefRuleBaseCL facetriang[8]=
{
    {1, {13, -1, -1, -1}},
    {2, {8,9,-1,-1}}, {2, {6,7,-1,-1}}, {3, {0,9,11,-1}},
    {2, {4,5,-1,-1}}, {3, {1,8,12,-1}}, {3, {3,6,10,-1}}, {4, {0,1,2,3}}
};

//Uint parentface[4][6]= { {1,6,2,8,9,3}, {0,5,2,7,9,3}, {0,4,1,7,8,3}, {0,4,1,5,6,2} };


FacCL GetFaceOfTetra(const TetCL& t, Uint num)
{
    FacCL f= t;
    f.erase(f.begin()+num);
    return f;
}

Uint GetEdgeOfTetra(const TetCL& t, Uint num)
{
    const Uint v0= t[VertOfEdge(num, 0)];
    const Uint v1= t[VertOfEdge(num, 1)];
    return EdgeByVert(v0, v1);
}

Uint GetFacePat(Uint refpat, Uint face)
{
    Uint facepat= 0;
    if ( RefinesEdge(refpat, EdgeOfFace(face, 0)) ) facepat |= 1;
    if ( RefinesEdge(refpat, EdgeOfFace(face, 1)) ) facepat |= 2;
    if ( RefinesEdge(refpat, EdgeOfFace(face, 2)) ) facepat |= 4;
    return facepat;
}

Uint GetTetraNumByVert(const TetCL& t)
{
    for (Uint ret=0; ret<NumAllChildrenC; ++ret)
        if (TetCL(VertOfChildAr[ret]+0, VertOfChildAr[ret]+4) == t)
            return ret;
    throw TopoErrCL("GetTetraNumByVert: not found");
}

Uint GetFaceNumByVert(const FacCL& f)
{
    for (Uint ret= 0; ret<NumAllFacesC; ++ret)
    {
        if (FacCL(VertOfFaceAr[ret]+0, VertOfFaceAr[ret]+3) == f)
            return ret;
    }
    throw TopoErrCL("GetFaceNumByVert: not found");
}

Uint GetFaceByVert(Uint v0, Uint v1, Uint v2)
{
    FacCL f(3);
    f[0]= v0; f[1]= v1; f[2]= v2;
    return GetFaceNumByVert( f );
}

bool RawTriangCL::IsOuterFace(const FacCL& f) const
{
    for (Uint face=0; face<faces.size(); ++face)
        for (Uint subface=0; subface< faces[face].size(); ++subface)
            if ( f==faces[face][subface]) return true;
    return false;
}


void RawTriangCL::Print(std::ostream& os) const
{
    os << children.size() << " Children; ";
    for (Uint i=0; i<children.size(); ++i)
        os << children[i][0] << " " << children[i][1] << " " << children[i][2] << " " << children[i][3] << "   ";
    os << '\n';
    for (Uint i=0; i<NumFacesC; ++i)
    {
        os << "Face " << i << ": " << faces[i].size() << " Subfaces; ";
        for (Uint j=0; j< faces[i].size(); ++j)
            os << faces[i][j][0] << " " << faces[i][j][1] << " " << faces[i][j][2] << "   ";
        os << "\n";
    }
    os.flush();
}


void ParseRaw(Uint rulenum, const RawRuleCL& R, RawTriangCL& r)
{
    r.faces.resize(4);
    r.children= std::vector<TetCL>(R.ChildNum, TetCL(4));
    for (int i=0; i<R.ChildNum; ++i)
        std::copy(R.Children[i]+0, R.Children[i]+4, r.children[i].begin());

    for (Uint i=0; i<NumFacesC; ++i)
    {
        int refpat= GetFacePat(rulenum, i);
        r.faces[i].resize(facetriang[refpat].ChildNum);
        for (Uint c=0; c<facetriang[refpat].ChildNum; ++c)
        {
            r.faces[i][c].resize(3);
            std::copy( VertOfFaceAr[SubFaceOfFace(i, facetriang[refpat].Children[c])]+0,
                       VertOfFaceAr[SubFaceOfFace(i, facetriang[refpat].Children[c])]+3,
                       r.faces[i][c].begin() );
        }
    }
}

/*
void
SetOuterFaces(const RawTriangCL& r, RefRuleCL& refrule, Uint& pos)
{
    refrule.FaceNum= r.faces[0].size() + r.faces[1].size() + r.faces[2].size() + r.faces[3].size();
//    Uint pos= 0;
    for (Uint face=0; face<NumFacesC; ++face)
    {
        for (Uint subface=0; subface<r.faces[face].size(); ++subface)
        {
            bool foundneighbor= false;
            refrule.Faces[pos].VertIndex= GetFaceNumByVert(r.faces[face][subface]);
            for (Uint child=0; !foundneighbor && child<r.children.size(); ++child)
                for (Uint childface=0; childface<NumFacesC; ++childface)
                    if ( GetFaceOfTetra(r.children[child], childface) == r.faces[face][subface] )
                    {
                        ++pos;
                        foundneighbor= true;
                        break;
                    }
        }
    }
}

void
SetInnerFaces(const RawTriangCL& r, RefRuleCL& refrule, Uint& pos)
{
    refrule.FaceNum+= ( r.children.size()*NumFacesC
                           -(r.faces[0].size() + r.faces[1].size() + r.faces[2].size() + r.faces[3].size()) )/2;
//    Uint pos= 0;
    for (Uint child0=0; child0<r.children.size(); ++child0)
    {
        for (Uint childface0=0; childface0<NumFacesC; ++childface0)
        {
            FacCL f= GetFaceOfTetra(r.children[child0], childface0);
            Uint fnum= GetFaceNumByVert(f);
            bool havethisface= false;
            for (Uint tmp=0; tmp < pos; ++tmp)
                if ( fnum==refrule.InnerFaces[tmp].VertIndex )
                {
                    havethisface= true;
                    break;
                }
            if ( !havethisface && !r.IsOuterFace(f) )
            {
                bool foundneighbor= false;
                refrule.Faces[pos]= fnum;
                for (Uint child1= child0+1; !foundneighbor && child1<r.children.size(); ++child1)
                    for (Uint childface1=0; childface1<NumFacesC; ++childface1)
                        if ( GetFaceOfTetra(r.children[child1], childface1) == f )
                        {
                            ++pos;
                            foundneighbor= true;
                            break;
                        }
            }
        }
    }
    for ( ; pos<MaxFacesC; ++pos)
    {
        refrule.Faces[pos].VertIndex= -1;
    }
}
*/

void SetFaces(const RawTriangCL& r, RefRuleCL& refrule)
{
    std::set<int> childFaces;

    for (Uint child=0; child<r.children.size(); ++child)
        for (Uint childface=0; childface<NumFacesC; ++childface)
        {
            const FacCL f= GetFaceOfTetra(r.children[child], childface);
            childFaces.insert( GetFaceNumByVert(f) );
        }

    refrule.FaceNum= childFaces.size();
    std::copy(childFaces.begin(), childFaces.end(), refrule.Faces);
    std::fill(refrule.Faces+refrule.FaceNum, refrule.Faces+MaxFacesC, -1);
}

void SetEdges(const RawTriangCL& r, RefRuleCL& refrule)
{
    std::set<int> childEdges;

    for (Uint i=0; i<r.children.size(); ++i)
        for (Uint j=0; j<NumEdgesC; ++j)
            childEdges.insert( GetEdgeOfTetra(r.children[i], j) );

    refrule.EdgeNum= childEdges.size();
    std::copy(childEdges.begin(), childEdges.end(), refrule.Edges);
    std::fill(refrule.Edges+refrule.EdgeNum, refrule.Edges+MaxEdgesC, -1);
}

void SetChildren(const RawTriangCL& r, RefRuleCL& refrule)
{
    refrule.ChildNum= r.children.size();
    for (Uint i=0; i<refrule.ChildNum; ++i)
        refrule.Children[i]= GetTetraNumByVert(r.children[i]);
    std::sort( refrule.Children, refrule.Children+refrule.ChildNum);
    std::fill( refrule.Children+refrule.ChildNum, refrule.Children+MaxChildrenC, -1);
}


void SetChildData(const TetCL& t, ChildDataCL& cdata)
{
    Uint VertIndex= GetTetraNumByVert(t);

    for (Uint v=0; v<NumVertsC; ++v)
        cdata.Vertices[v]= VertOfChildAr[VertIndex][v];
    for (Uint e=0; e<NumEdgesC; ++e)
        cdata.Edges[e]= EdgeByVert(t[VertOfEdge(e, 0)], t[VertOfEdge(e, 1)]);
    for (Uint f=0; f<NumFacesC; ++f)
        cdata.Faces[f]= GetFaceNumByVert( GetFaceOfTetra(t, f) );
}


void WriteChildData(const Uint VertIndex, std::ostream& os)
{
    Uint v[NumVertsC];

    os << "    ";
    for (Uint i=0; i<NumVertsC; ++i)
    {
        os << static_cast<int>(v[i]= VertOfChildAr[VertIndex][i]) << ", ";
    }
    os << "    ";
    for (Uint i=0; i<NumEdgesC; ++i)
    {
        os << static_cast<int>( EdgeByVert(v[VertOfEdge(i, 0)], v[VertOfEdge(i, 1)]) ) << ", ";
    }
    os << "    ";
    for (Uint i=0; i<NumFacesC; ++i)
    {
        os << static_cast<int>( GetFaceByVert(v[VertOfFace(i, 0)], v[VertOfFace(i, 1)], v[VertOfFace(i, 2)] ) );
        if (i != NumFacesC-1) os << ", ";
    }
}


void WriteRule(const RefRuleCL& rr, std::ostream& os)
{
    os <<   "    " << static_cast<int>(rr.EdgeNum) << ",    \t";
    for (Uint i=0; i<MaxEdgesC; ++i)
    {
        os << static_cast<int>(rr.Edges[i]) << ", ";
    }
    os << "\n    " << static_cast<int>(rr.FaceNum) << ",    \t";
    for (Uint i=0; i<MaxFacesC; ++i)
    {
        os << static_cast<int>(rr.Faces[i]) << ", ";
    }
    os << "\n    " << static_cast<int>(rr.ChildNum) << ",    \t";
    for (Uint i=0; i<MaxChildrenC; ++i)
    {
        os << static_cast<int>(rr.Children[i]);
        if (i != MaxChildrenC-1) os << ", ";
    }
}


void WriteRules(const std::vector<RefRuleCL>& refrules, std::ostream& os)
{
    os << "\n\nconst byte RefRuleAr[]= {\n";
    for (int i=0; i<64; ++i)
    {
        WriteRule(refrules[i], os);
        if (i!=63) os << ",\n\n";
    }
    os << "\n    };\n\n\n"
       << "const byte ChildDataAr[]= {\n";
    for (Uint i=0; i<NumAllChildrenC; ++i)
    {
        WriteChildData(i, os);
        if (i != NumAllChildrenC-1) os << ",\n";
    }
    os << "\n    };\n\n\n";
}


int main()
{
  try{
    std::vector<RawTriangCL> rawrules(64);
    std::vector<RefRuleCL> refrules(64);
    for(Uint i=0; i<64; ++i)
    {
        ParseRaw(i, RawRules[i], rawrules[i]);
//        rawrules[i].Print(std::cout);

        SetChildren(rawrules[i], refrules[i]);
        SetFaces(rawrules[i], refrules[i]);
        SetEdges(rawrules[i], refrules[i]);
    }
    std::string   TextLine;
    std::ifstream FileIn  ("topo.cpp.in");
    std::ofstream FileOut ("topo.cpp");
    if (!FileIn)  { std::cout << "Couldn't open file 'topo.cpp.in' for reading!" << std::endl; return -1; }
    if (!FileOut) { std::cout << "Couldn't open file 'topo.cpp' for writing!" << std::endl; return -1; }

    while (std::getline(FileIn,TextLine))
        if (TextLine.size()>=3 && TextLine.substr(0,3)==std::string("//@"))
            WriteRules(refrules, FileOut);
        else
            FileOut << TextLine << '\n';

    if (!FileIn.eof()) { std::cout << "Error reading 'topo.cpp.in'!" << std::endl; return -1; }
    if (!FileOut)      { std::cout << "Error writing 'topo.cpp'!" << std::endl; return -1; }
    FileIn.close();
    FileOut.close();
    return 0;
  }
  catch(TopoErrCL err){ err.handle(); }
}
