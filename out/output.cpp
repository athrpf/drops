/// \file output.cpp
/// \brief geometry and solution output in various formats
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen:

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

#include "out/output.h"
#ifdef _PAR
#  include "parallel/parallel.h"
#endif

namespace DROPS
{

std::ostream& GeomMGOutCL::put(std::ostream &os) const
{
    const int NumColors= 10;
    const char Color[NumColors][11] =
        {{" 1 0 0"}, {" 0 1 0"}, {" 0 0 1"}, {" 0 1 1"}, {" 1 0 1"},
         {" 1 1 0"}, {" 0.5 0 0"}, {" 0 0.5 0"}, {" 0 0 0.5"}, {" 0 0.5 0.5"}};
    Point3DCL GlobalOffset( 0.0);

#ifdef _PAR
    int me= ProcCL::MyRank(), numverts= 0;

    for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(), end= _MG->GetTriangVertexEnd();
            it!=end; ++it, ++numverts)
        GlobalOffset+= it->GetCoord();
    GlobalOffset*= _procExplode/numverts;
#endif

    os << "LIST {\n";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level); tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
#ifdef _PAR
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)
                        || tit->IsProcBnd(0) || tit->IsProcBnd(1) || tit->IsProcBnd(2) || tit->IsProcBnd(3)) )
#else
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)) )
#endif
            continue;
//        if (GetBaryCenter(*tit)[0]>0.5) continue;

        Point3DCL Offset(0.0);

        for ( TetraCL::const_VertexPIterator it=tit->GetVertBegin(), vend=tit->GetVertEnd(); it!=vend; ++it )
            Offset+= (*it)->GetCoord();
        Offset= (_explode/4)*Offset + GlobalOffset;

        os << "geom { OFF 4 4 6\n";
        for ( int i=0; i<4; i++ )
            for ( int j=0; j<3; j++ )
                os << Offset[j]+tit->GetVertex(i)->GetCoord()[j] << (j<2?" ":"\n");

#ifdef _PAR
        if (_procExplode)
            os <<   "3 1 2 3" << Color[me%10]
               << "\n3 0 2 3" << Color[me%10]
               << "\n3 0 1 3" << Color[me%10]
               << "\n3 0 1 2" << Color[me%10]
               << "\n}\n";
        else
#endif
        os <<   "3 1 2 3" << (tit->IsBndSeg(0)?Color[tit->GetBndIdx(0)%NumColors]:" 0.7 0.7 0.7")
           << "\n3 0 2 3" << (tit->IsBndSeg(1)?Color[tit->GetBndIdx(1)%NumColors]:" 0.7 0.7 0.7")
           << "\n3 0 1 3" << (tit->IsBndSeg(2)?Color[tit->GetBndIdx(2)%NumColors]:" 0.7 0.7 0.7")
           << "\n3 0 1 2" << (tit->IsBndSeg(3)?Color[tit->GetBndIdx(3)%NumColors]:" 0.7 0.7 0.7")
           << "\n}\n";
    }
    return os << "}" << std::endl;
}


std::ostream& DumpMGCL::put (std::ostream& os) const
{
    const int lastLevel ( static_cast<int>(_MG->GetLastLevel()) );

    for (int Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_VertexIterator vCIt( _MG->GetVerticesBegin( Level));
              vCIt!=_MG->GetVerticesEnd( Level); ++vCIt)
        {
            vCIt->DebugInfo(os);
        }
        os << std::endl;
    }
    os << std::endl;
    for (int Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_EdgeIterator eCIt( _MG->GetEdgesBegin( Level));
              eCIt!=_MG->GetEdgesEnd( Level); ++eCIt )
        {
            eCIt->DebugInfo(os);
        }
        os << std::endl;
    }
    os << std::endl;
    for (int Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_FaceIterator fCIt( _MG->GetFacesBegin( Level));
              fCIt!=_MG->GetFacesEnd( Level); ++fCIt )
        {
            fCIt->DebugInfo(os);
        }
        os << std::endl;
    }
    os << std::endl;
    for (int Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_TetraIterator tCIt( _MG->GetTetrasBegin( Level));
              tCIt!=_MG->GetTetrasEnd( Level); ++tCIt)
        {
            tCIt->DebugInfo(os);
        }
        os << std::endl;
    }
    return os << std::endl;
}


std::ostream& SanityMGOutCL::put (std::ostream& os) const
{
   if ( _MG->IsSane(os) ) os << "As far as I can tell the MultigridCL is sane.";
   return os;
}


TetraSectionT
intersect(const TetraCL& t, const PlaneCL& pl)
{
    std::vector<Uint> before;
    std::vector<Uint> in;
    std::vector<Uint> behind;
    before.reserve(4); in.reserve(4), behind.reserve(4);

    // Decide, where each vertex of t lies with respect to pl!
    for (Uint i=0; i<NumVertsC; ++i)
    {
        const double val= inner_prod(pl.n, t.GetVertex(i)->GetCoord()) - pl.r;
        if (val > 0.) before.push_back(i);
        else if (val < 0.) behind.push_back(i);
             else in.push_back(i);
    }

    // Now determine the cut polygon!
    const Uint bef_c= before.size();
    const Uint beh_c= behind.size();
    TetraSectionT ret;

    if (bef_c == NumVertsC || beh_c == NumVertsC) // no intersection
        return ret;
    else if (bef_c == 3 && beh_c == 1) // triangular intersection
    {
        for (Uint i=0; i<3; ++i)
        {
            const double p= CutEdge(t.GetVertex(before[i])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
            ret.push_back( ConvexComb(p, std_basis<3>(before[i]), std_basis<3>(behind[0])) );
        }
        return ret;
    }
    else if (bef_c == 1 && beh_c == 3) // as above
    {
        for (Uint i=0; i<3; ++i)
        {
            const double p= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[i])->GetCoord(), pl);
            ret.push_back( ConvexComb(p, std_basis<3>(before[0]), std_basis<3>(behind[i])) );
        }
        return ret;
    }
    else if (bef_c == 2 && beh_c == 2) // quadrilateral intersection
    {
        const double p0= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
        ret.push_back( ConvexComb(p0, std_basis<3>(before[0]), std_basis<3>(behind[0])) );
        const double p1= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[1])->GetCoord(), pl);
        ret.push_back( ConvexComb(p1, std_basis<3>(before[0]), std_basis<3>(behind[1])) );
        const double p2= CutEdge(t.GetVertex(before[1])->GetCoord(), t.GetVertex(behind[1])->GetCoord(), pl);
        ret.push_back( ConvexComb(p2, std_basis<3>(before[1]), std_basis<3>(behind[1])) );
        const double p3= CutEdge(t.GetVertex(before[1])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
        ret.push_back( ConvexComb(p3, std_basis<3>(before[1]), std_basis<3>(behind[0])) );
        return ret;
    }
    // From here on, we must deal with vertices in the cut plane; write them to ret immediatly!
    for (Uint i=0; i<in.size(); ++i)
        ret.push_back(std_basis<3>(in[i]));
    if ( (bef_c == 3 && beh_c == 0) || (bef_c == 0 && beh_c == 3) ) // one vertex in intersection
        return ret;
    else if (bef_c == 2 && beh_c == 1) // triangular intersection; cut two edges
    {
        const double p0= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
        ret.push_back( ConvexComb(p0, std_basis<3>(before[0]), std_basis<3>(behind[0])) );
        const double p1= CutEdge(t.GetVertex(before[1])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
        ret.push_back( ConvexComb(p1, std_basis<3>(before[1]), std_basis<3>(behind[0])) );
        return ret;
    }
    else if (bef_c == 1 && beh_c == 2) // as above
    {
        const double p0= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
        ret.push_back( ConvexComb(p0, std_basis<3>(before[0]), std_basis<3>(behind[0])) );
        const double p1= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[1])->GetCoord(), pl);
        ret.push_back( ConvexComb(p1, std_basis<3>(before[0]), std_basis<3>(behind[1])) );
        return ret;

    }
    else if ( (bef_c == 2 && beh_c == 0) || (bef_c == 0 && beh_c == 2) ) // edge as intersection
        return ret;
    else if (bef_c == 1 && beh_c == 1) // triangular intersection; cut one edge
    {
        const double p0= CutEdge(t.GetVertex(before[0])->GetCoord(), t.GetVertex(behind[0])->GetCoord(), pl);
        ret.push_back( ConvexComb(p0, std_basis<3>(before[0]), std_basis<3>(behind[0])) );
        return ret;
    }
    else if ( (bef_c == 1 && beh_c == 0) || (bef_c == 0 && beh_c == 1) ) // face in cut plane
        return ret;

    // We should never reach this point...
    throw DROPSErrCL("intersect: I blew it!");
}


bool
Maple3DOptionCL::WriteGlobalOptions(std::ostream& os) const
{
    bool have_previous= false;

    if ( !title.empty() )
    {
        os << "TITLE(`" << title << "`)";
        have_previous= true;
    }
    if ( !(xlabel.empty() && ylabel.empty() && zlabel.empty()) )
    {
        if (have_previous) os << ", ";
        os << "AXESLABELS(`" << xlabel << "`, `" << ylabel << "`, `" << zlabel << "`)";
        have_previous= true;
    }
    if ( !axesstyle.empty() )
    {
        if (have_previous) os << ", ";
        os << "AXESSTYLE(`" << axesstyle << "`)";
        have_previous= true;
    }
    if ( !scaling.empty() )
    {
        if (have_previous) os << ", ";
        os << "SCALING(`" << scaling << "`)";
        have_previous= true;
    }
    return have_previous;
}


bool
MapleMGOutCL::_Write_Intersection(const TetraCL& t, std::ostream& os, bool have_previous) const
{
    TetraSectionT section= intersect(t, this->_plane);
    if (section.empty()) return have_previous;
    if (have_previous) os << ", ";

    os << '[';
    for (Uint i=0; i<section.size(); ++i)
    {
        if (i!=0) os << ", ";
        const SVectorCL<3>& point= GetWorldCoord(t, section[i]);
        os << "[" << point[0] << ", " << point[1] << ", " << point[2] << ']';
    }
    os << ']';
    return true;
}


bool
MapleMGOutCL::_Write_Tetra(const TetraCL& t, std::ostream& os, bool have_previous) const
{
    if (have_previous) os << ", ";

    for (Uint i=0; i<NumFacesC; ++i)
    {
        if (i!=0) os << ", ";
        const SVectorCL<3>& p0= t.GetVertex(VertOfFace(i, 0))->GetCoord();
        const SVectorCL<3>& p1= t.GetVertex(VertOfFace(i, 1))->GetCoord();
        const SVectorCL<3>& p2= t.GetVertex(VertOfFace(i, 2))->GetCoord();
        os << "[[" << p0[0] << ", " << p0[1] << ", " << p0[2] << "], "
           << "[" << p1[0] << ", " << p1[1] << ", " << p1[2] << "], "
           << "[" << p2[0] << ", " << p2[1] << ", " << p2[2] << "]]";
    }
    return true;
}


void
MapleMGOutCL::_Intersect_Plot(std::ostream& os) const
{
    os << "POLYGONS(";
    bool have_previous= false;
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level); tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)) )
            continue;

        have_previous= _cut ? _Write_Intersection(*tit, os, have_previous)
                            : _Write_Tetra(*tit, os, have_previous);
    }
}


std::ostream&
MapleMGOutCL::put(std::ostream& os) const
{
    os << _options.name << ":= PLOT3D(";
    if ( _options.WriteGlobalOptions(os) ) os << ", ";
    _Intersect_Plot(os);
    return os << ")):" << std::endl;
}


/// \brief Write finite element function, stored in \a v, in a file, named \a filename
void WriteFEToFile( const VecDescCL& v, MultiGridCL& mg, std::string filename, bool binary, const VecDescCL* lsetp)
{
    if (!v.RowIdx->IsExtended()) {
#ifdef _PAR
        ProcCL::AppendProcNum( filename);
#endif
        std::ofstream file( filename.c_str());
        if (!file) throw DROPSErrCL("WriteFEToFile: Cannot open file "+filename+" for writing");
        v.Write( file, binary);
    }
    else { // extended FE
        IdxDescCL p1( P1_FE);
        p1.CreateNumbering( v.RowIdx->TriangLevel(), mg, *v.RowIdx);
        VecDescCL vpos(&p1), vneg(&p1);
        P1XtoP1 ( *v.RowIdx, v.Data, p1, vpos.Data, vneg.Data, *lsetp, mg);
        WriteFEToFile(vneg, mg, filename + "Neg");
        WriteFEToFile(vpos, mg, filename + "Pos");
        p1.DeleteNumbering(mg);
    }
}

/// Read a serialized finite element function from a file
/// \pre CreateNumbering of v.RowIdx must have been called before
void ReadFEFromFile( VecDescCL& v, MultiGridCL& mg, std::string filename, bool binary, const VecDescCL* lsetp)
{
    if (!v.RowIdx->IsExtended()) {

        std::cout << "Read FE "<<filename<<std::endl;
#ifdef _PAR
        ProcCL::AppendProcNum( filename);
#endif
        std::ifstream file( filename.c_str());
        if (!file) throw DROPSErrCL("ReadFEFromFile: Cannot open file "+filename);
        v.Read( file, binary);
    }
    else { // extended FE
        IdxDescCL p1( P1_FE);
        p1.CreateNumbering( v.RowIdx->TriangLevel(), mg, *v.RowIdx);
        VecDescCL vpos(&p1), vneg(&p1);
        ReadFEFromFile(vneg, mg, filename + "Neg");
        ReadFEFromFile(vpos, mg, filename + "Pos");
        P1toP1X ( *v.RowIdx, v.Data, p1, vpos.Data, vneg.Data, *lsetp, mg);
        p1.DeleteNumbering(mg);
    }
}


} // end of namespace DROPS
