//**************************************************************************
// File:    output.cpp                                                       *
// Content: geometry and solution output in various formats                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - Nov, 26 2001                                           *
//**************************************************************************

#ifndef _OUTPUT_CPP_
#define _OUTPUT_CPP_

#include "out/output.h"

namespace DROPS
{

std::ostream& GeomMGOutCL::put(std::ostream &os) const
{
    const char Color[][13] =
        {" 1 0 0", " 0 1 0", " 0 0 1", " 0 1 1", " 1 0 1", " 1 1 0",
         " 0.5 0 0", " 0 0.5 0", " 0 0 0.5", " 0 0.5 0.5"};

    os << "LIST {\n";
    for ( MultiGridCL::const_TriangTetraIteratorCL tit=_MG->GetTriangTetraBegin(_level); tit!=_MG->GetTriangTetraEnd(_level); ++tit )
    {
        if ( _onlyBnd && !(tit->IsBndSeg(0) || tit->IsBndSeg(1) || tit->IsBndSeg(2) || tit->IsBndSeg(3)) )
            continue;
//        if (GetBaryCenter(*tit)[0]>0.5) continue;

        Point3DCL Offset(0.0);

        for ( TetraCL::const_VertexPIterator it=tit->GetVertBegin(), vend=tit->GetVertEnd(); it!=vend; ++it )
            Offset+= (*it)->GetCoord();
        os << "geom { OFF 4 4 6\n";
        for ( int i=0; i<4; i++ )
            for ( int j=0; j<3; j++ )
                os << _explode*Offset[j]/4+tit->GetVertex(i)->GetCoord()[j] << (j<2?" ":"\n");
        os <<   "3 1 2 3" << (tit->IsBndSeg(0)?Color[tit->GetBndIdx(0)%10]:" 0.7 0.7 0.7")
           << "\n3 0 2 3" << (tit->IsBndSeg(1)?Color[tit->GetBndIdx(1)%10]:" 0.7 0.7 0.7")
           << "\n3 0 1 3" << (tit->IsBndSeg(2)?Color[tit->GetBndIdx(2)%10]:" 0.7 0.7 0.7")
           << "\n3 0 1 2" << (tit->IsBndSeg(3)?Color[tit->GetBndIdx(3)%10]:" 0.7 0.7 0.7")
           << "\n}\n";
    }
    return os << "}" << std::endl;
}


std::ostream& DumpMGCL::put (std::ostream& os) const
{
    int lastLevel ( static_cast<int>(_MG->GetLastLevel()) );
    int Level;

    for (Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_VertexIterator vCIt(_MG->GetVertices()[Level].begin());
              vCIt!=_MG->GetVertices()[Level].end(); ++vCIt)
        {
            vCIt->DebugInfo(os);
        }
        os << std::endl;
    }
    os << std::endl;
    for (Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_EdgeIterator eCIt(_MG->GetEdges()[Level].begin());
              eCIt!=_MG->GetEdges()[Level].end(); ++eCIt )
        {
            eCIt->DebugInfo(os);
        }
        os << std::endl;
    }
    os << std::endl;
    for (Level=0; Level<=lastLevel; ++Level)
    {
        for ( MultiGridCL::const_TetraIterator tCIt(_MG->GetTetras()[Level].begin());
              tCIt!=_MG->GetTetras()[Level].end(); ++tCIt)
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


} // end of namespace DROPS


#include "out/mapleout.cpp"


#endif
