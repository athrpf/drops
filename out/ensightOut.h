//**************************************************************************
// File:    ensightOut.h                                                   *
// Content: solution output in Ensight6 Case format                        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_ENSIGHTOUT_H
#define DROPS_ENSIGHTOUT_H

#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

namespace DROPS
{

class EnsightP2SolOutCL
// output of solution in Ensight6 Case format
{
  private:
    const MultiGridCL* _MG;
    const IdxDescCL*   _idx;
    char               _decDigits;
    Uint               _timestep, _numsteps;
    double             _lasttime;
    std::ostringstream _timestr;
    std::ofstream      _case;
    std::string        _geom;
     
    void AppendTimecode( std::string&) const;
    void putTime( double time);
  
  public:
    // idx must refer to a numbering of the verts and edges of a certain
    // triangulation, i.e. use LevelsetP2CL::CreateNumbering before
    // constructing an EnsightSolOutCL-object.
    EnsightP2SolOutCL( const MultiGridCL& mg, const IdxDescCL* idx)
      : _MG( &mg), _idx( idx), _decDigits( 0), _timestep( 0), _numsteps( 0)
    {}
    ~EnsightP2SolOutCL() { if (_case.is_open()) CaseEnd(); }
      
    // call CaseBegin() before any other member function
    void CaseBegin ( const char casefileName[], Uint numsteps= 0);
    void DescribeGeom  ( const char geoName[], std::string fileName, bool timedep= false);
    void DescribeScalar( const char varName[], std::string fileName, bool timedep= false);
    void DescribeVector( const char varName[], std::string fileName, bool timedep= false);
    void putGeom      ( std::string, double t= -1);
    template<class DiscScalT>
    void putScalar    ( std::string, const DiscScalT&, double t= -1);
    template<class DiscVecT>
    void putVector    ( std::string, const DiscVecT&, double t= -1);
    // call putCaseEnd() after finishing all other output
    void CaseEnd   ();
};

//=====================================================
//              template definitions
//=====================================================

void EnsightP2SolOutCL::CaseBegin( const char casefileName[], Uint numsteps)
{
    _case.open( casefileName);
    _case << "FORMAT\ntype: ensight\n\n";
    
    _numsteps= numsteps;
    _decDigits= 1;
    while( numsteps>9)
    { ++_decDigits; numsteps/=10; }
}

void EnsightP2SolOutCL::DescribeGeom( const char geomName[], std::string fileName, bool timedep)
{
    _case << "GEOMETRY\nmodel:\t\t\t";
    if (timedep)
    {
        _case << '1';
        fileName+= std::string( _decDigits, '*');
    }
    _case << "\t\t\t" << fileName << "\n\nVARIABLE\n";
    _geom= geomName;
}

void EnsightP2SolOutCL::DescribeScalar( const char varName[], std::string fileName, bool timedep)
{
    _case << "scalar per node:\t";
    if (timedep)
    {
        _case << '1';
        fileName+= std::string( _decDigits, '*');
    }
    _case << '\t' << varName << '\t' << fileName << '\n';
}

void EnsightP2SolOutCL::DescribeVector( const char varName[], std::string fileName, bool timedep)
{
    _case << "vector per node:\t";
    if (timedep)
    {
        _case << '1';
        fileName+= std::string( _decDigits, '*');
    }
    _case << '\t' << varName << '\t' << fileName << '\n';
}

void EnsightP2SolOutCL::CaseEnd()
{
    if (_timestr.str().empty()) return;
    
    _case << "\nTIME\ntime set:\t\t1\nnumber of steps:\t" << _timestep
          << "\nfilename start number:\t0\nfilename increment:\t1\ntime values:\t\t";
    _case << _timestr.str() << "\n\n";

    _case.close();
}


void EnsightP2SolOutCL::putGeom( std::string fileName, double t)
{
    const IdxDescCL* idxDesc= _idx;
    const Uint lvl= _idx->TriangLevel,
               idx= _idx->Idx;
    
    if ( t!=-1)
    {
        AppendTimecode( fileName);
        putTime( t);
    }    
    
    std::ofstream os( fileName.c_str());
    os.flags(std::ios_base::scientific);
    os.precision(5);
    os.width(12);

    os << "DROPS geometry file: " << "\nformat: Ensight6 Case format\n";
    os << "node id given\nelement id off\n";
    os << "coordinates\n" << std::setw(8) << idxDesc->NumUnknowns << '\n';
    
    for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
        end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        os << std::setw(8) << it->Unknowns(idx)+1 << ' ' << it->GetCoord() << '\n';
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
        end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        os << std::setw(8) << it->Unknowns(idx)+1 << ' ' << GetBaryCenter( *it) << '\n';
    }

    os << "part 1\n" << _geom << "\n";
    // Ausgabe auf virtuell reg. verfeinertem Gitter, da Ensight zur Berechnung
    // der Isoflaechen anscheinend nur die Werte in den Vertices beruecksichtigt...
    os << "tetra4\n" 
       << std::setw(8) << 8*std::distance(_MG->GetTriangTetraBegin(lvl),_MG->GetTriangTetraEnd(lvl)) << '\n'; 
    for (MultiGridCL::const_TriangTetraIteratorCL it= _MG->GetTriangTetraBegin(lvl),
        end= _MG->GetTriangTetraEnd(lvl); it!=end; ++it)
    {
        RefRuleCL RegRef= GetRefRule( RegRefRuleC);
        for (int ch=0; ch<8; ++ch)
        {
            ChildDataCL data= GetChildData( RegRef.Children[ch]);
            for (int vert= 0; vert<4; ++vert)
            {
                int v= data.Vertices[vert];
                os << std::setw(8);
                if (v<4) 
                    os << it->GetVertex(v)->Unknowns(idx)+1;
                else
                    os << it->GetEdge(v-4)->Unknowns(idx)+1;
            }
            os << '\n';
        }
    }
    
/*  // alte Version mit P2-Ausgabe auf tatsaechlichem Gitter
    os << "tetra10\n" 
       << std::setw(8) << std::distance(_MG->GetTriangTetraBegin(lvl),_MG->GetTriangTetraEnd(lvl)) << '\n'; 
    for (MultiGridCL::const_TriangTetraIteratorCL it= _MG->GetTriangTetraBegin(lvl),
        end= _MG->GetTriangTetraEnd(lvl); it!=end; ++it)
    {
        for (int i=0; i<4; ++i)
            os << std::setw(8) << it->GetVertex(i)->Unknowns(idx)+1;
        // the edge order of Ensight is a little bit strange...
        os << std::setw(8) << it->GetEdge(0)->Unknowns(idx)+1;
        os << std::setw(8) << it->GetEdge(2)->Unknowns(idx)+1;
        os << std::setw(8) << it->GetEdge(1)->Unknowns(idx)+1;
        os << std::setw(8) << it->GetEdge(3)->Unknowns(idx)+1;
        os << std::setw(8) << it->GetEdge(4)->Unknowns(idx)+1;
        os << std::setw(8) << it->GetEdge(5)->Unknowns(idx)+1;
        os << '\n';
    }
*/
}

template<class DiscScalT>
void EnsightP2SolOutCL::putScalar( std::string fileName, const DiscScalT& v, double t)
{
    const Uint lvl= _idx->TriangLevel;
    int cnt= 0;
    
    if ( t!=-1)
    {
        AppendTimecode( fileName);
        putTime( t);
    }    

    std::ofstream os( fileName.c_str());
    os.flags(std::ios_base::scientific);
    os.precision(5);
    os.width(12);

    os << "DROPS data file, scalar variable:\n";
    for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
        end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        os << std::setw(12) << v.val( *it);
        if ( (++cnt)==6)
        { // Ensight expects six real numbers per line
            cnt= 0;
            os << '\n';
        }
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
        end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        os << std::setw(12) << v.val( *it);
        if ( (++cnt)==6)
        { // Ensight expects six real numbers per line
            cnt= 0;
            os << '\n';
        }
    }
    os << '\n';
}

template<class DiscVecT>
void EnsightP2SolOutCL::putVector( std::string fileName, const DiscVecT& v, double t)
{
    const Uint lvl= _idx->TriangLevel;
    int cnt= 0;
    
    if ( t!=-1)
    {
        AppendTimecode( fileName);
        putTime( t);
    }    

    std::ofstream os( fileName.c_str());
    os.flags(std::ios_base::scientific);
    os.precision(5);
    os.width(12);

    os << "DROPS data file, vector variable:\n";
    for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
        end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        for (int i=0; i<3; ++i)
            os << std::setw(12) << v.val( *it)[i];
        if ( (cnt+=3)==6)
        { // Ensight expects six real numbers per line
            cnt= 0;
            os << '\n';
        }
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
        end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        for (int i=0; i<3; ++i)
            os << std::setw(12) << v.val( *it)[i];
        if ( (cnt+=3)==6)
        { // Ensight expects six real numbers per line
            cnt= 0;
            os << '\n';
        }
    }
    os << '\n';
}

void EnsightP2SolOutCL::AppendTimecode( std::string& str) const
{
    char format[]= "%0Xi",
         postfix[8];
    format[2]= '0' + char(_decDigits);
    sprintf( postfix, format, _timestep);
    str+= postfix;
}

void EnsightP2SolOutCL::putTime( double t)
{
    if (t!=_lasttime)
    {
        _timestr << t << ' ';
        _lasttime= t;
        if (++_timestep%10==0)
            _timestr << "\n\t\t\t";
    }
}

} // end of namespace DROPS

#endif
