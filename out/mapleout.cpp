namespace DROPS
{

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


} // end of namespace DROPS
