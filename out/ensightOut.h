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
#include "geom/multigrid.h"
#include "misc/problem.h"

namespace DROPS
{

union showInt
{
    int i;
    char s[sizeof(int)];
};

union showFloat
{
    float f;
    char s[sizeof(float)];
};


class EnsightP2SolOutCL
// output of solution in Ensight6 Case format
{
  private:
    const MultiGridCL* _MG;
    const IdxDescCL*   _idx;
    char               _decDigits;
    Uint               _timestep, _numsteps;
    double             _lasttime;
    std::ostringstream _descstr,    // stores description info
                       _timestr;    // stores time info
    std::ofstream      _case;
    std::string        _geom;
    const bool         binary_;

    void AppendTimecode( std::string&) const;
    void putTime( double time);
    void CheckFile( const std::ofstream&) const;

  public:
    // idx must refer to a numbering of the verts and edges of a certain
    // triangulation, i.e. use LevelsetP2CL::CreateNumbering before
    // constructing an EnsightSolOutCL-object.
    EnsightP2SolOutCL( const MultiGridCL& mg, const IdxDescCL* idx, bool binary=true)
      : _MG( &mg), _idx( idx), _decDigits( 0), _timestep( -1u), _numsteps( 0), _lasttime(-1), binary_( binary)
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
    void Commit    ();    // rewrites case file
    // call CaseEnd() after finishing all other output
    void CaseEnd   ();
};


class ReadEnsightP2SolCL
// read solution from Ensight6 Case format
{
  private:
    const MultiGridCL* _MG;
    const bool         binary_;

    void CheckFile( const std::ifstream&) const;

  public:
    ReadEnsightP2SolCL( const MultiGridCL& mg, bool binary=true)
      : _MG(&mg), binary_(binary) {}

    template<class BndT>
    void ReadScalar( const std::string&, VecDescCL&, const BndT&) const;
    template<class BndT>
    void ReadVector( const std::string&, VecDescCL&, const BndT&) const;
};



//=====================================================
//              template definitions
//=====================================================

void EnsightP2SolOutCL::CaseBegin( const char casefileName[], Uint numsteps)
{
    _case.open( casefileName);
    CheckFile( _case);
    _descstr << "FORMAT\ntype: ensight\n\n";

    _numsteps= numsteps;
    _decDigits= 1;
    while( numsteps>9)
    { ++_decDigits; numsteps/=10; }
}

void EnsightP2SolOutCL::DescribeGeom( const char geomName[], std::string fileName, bool timedep)
{
    _descstr << "GEOMETRY\nmodel:\t\t\t";
    if (timedep)
    {
        _descstr << '1';
        fileName+= std::string( _decDigits, '*');
    }
    _descstr << "\t\t\t" << fileName << "\n\nVARIABLE\n";
    _geom= geomName;
}

void EnsightP2SolOutCL::DescribeScalar( const char varName[], std::string fileName, bool timedep)
{
    _descstr << "scalar per node:\t";
    if (timedep)
    {
        _descstr << '1';
        fileName+= std::string( _decDigits, '*');
    }
    _descstr << '\t' << varName << '\t' << fileName << std::endl;
}

void EnsightP2SolOutCL::DescribeVector( const char varName[], std::string fileName, bool timedep)
{
    _descstr << "vector per node:\t";
    if (timedep)
    {
        _descstr << '1';
        fileName+= std::string( _decDigits, '*');
    }
    _descstr << '\t' << varName << '\t' << fileName << std::endl;
}

void EnsightP2SolOutCL::AppendTimecode( std::string& str) const
{
    char format[]= "%0Xi",
         postfix[8];
    format[2]= '0' + char(_decDigits);
    std::sprintf( postfix, format, _timestep);
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

        Commit();    // rewrite case file
    }
}

void EnsightP2SolOutCL::CheckFile( const std::ofstream& os) const
{
    if (!os) throw DROPSErrCL( "EnsightP2SolOutCL: error while opening file!");
}

void EnsightP2SolOutCL::Commit()
{
    // rewrite case file
    _case.seekp( 0, std::ios_base::beg);  // rewind to the very beginning
    _case << _descstr.str();
    if (!_timestr.str().empty())
    {
        _case << "\nTIME\ntime set:\t\t1\nnumber of steps:\t" << _timestep+1
              << "\nfilename start number:\t0\nfilename increment:\t1\ntime values:\t\t";
        _case << _timestr.str() << "\n\n";
    }
    _case.flush();
}

void EnsightP2SolOutCL::CaseEnd()
{
    Commit();
    _case.close();
}


void EnsightP2SolOutCL::putGeom( std::string fileName, double t)
{
    const IdxDescCL* idxDesc= _idx;
    const Uint lvl= _idx->TriangLevel,
               idx= _idx->GetIdx();
    if ( t!=-1)
    {
        putTime( t);
        AppendTimecode( fileName);
    }

    std::ofstream os( fileName.c_str());
    CheckFile( os);

    if(binary_)
    {
        char buffer[80];
        std::strcpy(buffer,"C Binary");     //writing of all necessary information: Binary header
        os.write(buffer,80);
        std::strcpy(buffer,"DROPS geometry file: ");         //description line 1
        os.write(buffer,80);
        std::strcpy(buffer,"format: Ensight6 Case format");  //descripiton line 2
        os.write(buffer,80);
        std::strcpy(buffer,"node id given");                 //node id status
        os.write(buffer,80);
        std::strcpy(buffer,"element id off");                //element id status
        os.write(buffer,80);
        std::strcpy(buffer,"coordinates");                   //coordinates line
        os.write(buffer,80);

        showInt sInt;                  //unions for converting ASCII int to binary
        showFloat sFlo;
        sInt.i= (int)idxDesc->NumUnknowns;            //write number of nodes
        os.write( sInt.s, sizeof(int));
        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),     //write node ids
                 end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            sInt.i=(int)it->Unknowns(idx)+1;
            os.write( sInt.s, sizeof(int));
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),          //write ids of additional nodes
                 end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            sInt.i=(int)it->Unknowns(idx)+1;
            os.write( sInt.s, sizeof(int));
        }

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),     //write coordinates of nodes
                 end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            Point3DCL c= it->GetCoord();
            for(int i=0; i<3; ++i)
            {
                sFlo.f=c[i];
                os.write( sFlo.s, sizeof(float));
            }
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),           //write coordinates of additional nodes
                 end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            Point3DCL c= GetBaryCenter( *it);
            for(int i=0; i<3; ++i)
            {
                sFlo.f=c[i];
                os.write( sFlo.s, sizeof(float));
            }
        }

        std::strcpy(buffer,"part 1");                          //part no. line
        os.write(buffer,80);
        int b=0;
        char buff[20];
        for(b=0;b<=20;b++)                                //part description line
            buff[b]=_geom[b];
        std::strcpy(buffer,buff);
        os.write(buffer,80);
        std::strcpy(buffer,"tetra4");                          // element 1 tetrahedron with 4 nodes each
        os.write(buffer,80);

        sInt.i =  8*std::distance(_MG->GetTriangTetraBegin(lvl),_MG->GetTriangTetraEnd(lvl));     //number of tetrahedra
        os.write( sInt.s, sizeof(int));

        RefRuleCL RegRef= GetRefRule( RegRefRuleC);
        for (MultiGridCL::const_TriangTetraIteratorCL it= _MG->GetTriangTetraBegin(lvl),
                      end= _MG->GetTriangTetraEnd(lvl); it!=end; ++it)
        {
            for (int ch=0; ch<8; ++ch)
            {
                ChildDataCL data= GetChildData( RegRef.Children[ch]);

                for (int vert= 0; vert<4; ++vert)
                {
                    int v= data.Vertices[vert];
                    if (v<4)
                    {
                        sInt.i= it->GetVertex(v)->Unknowns(idx)+1;
                        os.write( sInt.s, sizeof(int));
                    }
                    else
                    {
                        sInt.i=it->GetEdge(v-4)->Unknowns(idx)+1;
                        os.write( sInt.s, sizeof(int));
                    }
                }
            }
        }
     }
     else // hier startet die normale ASCII-Ausgabe
     {
         os.flags(std::ios_base::scientific);
         os.precision(5);
         os.width(12);

         os << "DROPS geometry file: " << "\nformat: Ensight6 Case format\n";
         os << "node id given\nelement id off\n";
         os << "coordinates\n" << std::setw(8) << idxDesc->NumUnknowns << '\n';

         for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
         {
             os << std::setw(8) << it->Unknowns(idx)+1;
             Point3DCL c= it->GetCoord();
             for (int i=0; i<3; ++i)
                 os << std::setw(12) << c[i];
             os << '\n';
         }
         for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
             end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
         {
             os << std::setw(8) << it->Unknowns(idx)+1;
             Point3DCL c= GetBaryCenter( *it);
             for (int i=0; i<3; ++i)
                 os << std::setw(12) << c[i];
             os << '\n';
         }

         os << "part 1\n" << _geom << "\n";
         // Ausgabe auf virtuell reg. verfeinertem Gitter, da Ensight zur Berechnung
         // der Isoflaechen anscheinend nur die Werte in den Vertices beruecksichtigt...
         os << "tetra4\n"
            << std::setw(8) << 8*std::distance(_MG->GetTriangTetraBegin(lvl),_MG->GetTriangTetraEnd(lvl)) << '\n';

         RefRuleCL RegRef= GetRefRule( RegRefRuleC);
         for (MultiGridCL::const_TriangTetraIteratorCL it= _MG->GetTriangTetraBegin(lvl),
             end= _MG->GetTriangTetraEnd(lvl); it!=end; ++it)
         {
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
    char buffer[80];
    showFloat sFlo;

    if ( t!=-1)
    {
        putTime( t);
        AppendTimecode( fileName);
    }

    std::ofstream os( fileName.c_str());
    CheckFile( os);

    if(binary_)
    {
        std::strcpy(buffer,"DROPS data file, scalar variable:");
        os.write(buffer,80);

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
           end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            sFlo.f=v.val(*it);
            os.write(sFlo.s,sizeof(float));
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
           end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            sFlo.f=v.val(*it,0.5);
            os.write(sFlo.s,sizeof(float));
        }
    }
    else //ASCII-Ausgabe
    {
        int cnt=0;
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
            os << std::setw(12) << v.val( *it, 0.5);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }
        os << '\n';
    }
}

template<class DiscVecT>
void EnsightP2SolOutCL::putVector( std::string fileName, const DiscVecT& v, double t)
{
    const Uint lvl= _idx->TriangLevel;
    char buffer[80];
    showFloat sFlo;

    if ( t!=-1)
    {
        putTime( t);
        AppendTimecode( fileName);
    }

    std::ofstream os( fileName.c_str());
    CheckFile( os);
    os.flags(std::ios_base::scientific);

    if(binary_)
    {
        std::strcpy(buffer,"DROPS data file, vector variable:");
        os.write( buffer, 80);

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            for (int i=0; i<3; ++i)
            {
                sFlo.f=v.val( *it)[i];
                os.write(sFlo.s,sizeof(float));
            }
        }
        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            for (int i=0; i<3; ++i)
            {
                 sFlo.f=v.val( *it)[i];
                 os.write(sFlo.s,sizeof(float));
            }
        }
    }
    else
    { // ASCII
        int cnt=0;
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
}

// ========== ReadEnsightP2SolCL ==========

template <class BndT>
void ReadEnsightP2SolCL::ReadScalar( const std::string& file, VecDescCL& v, const BndT& bnd) const
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
    std::ifstream is( file.c_str());
    CheckFile( is);

    if (binary_)
    {
        showFloat fl;
        char buffer[80];
        is.read( buffer, 80);           //ignore first 80 characters

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= (double)fl.f;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= fl.f;
        }
    }
    else
    { // ASCII
        char buf[256];
        double d= 0;

        is.getline( buf, 256); // ignore first line

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is >> d;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= d;
        }
        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is >> d;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= d;
        }
    }

    CheckFile( is);
}

template <class BndT>
void ReadEnsightP2SolCL::ReadVector( const std::string& file, VecDescCL& v, const BndT& bnd) const
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
    std::ifstream is( file.c_str());
    CheckFile( is);

    double d0= 0, d1= 0, d2= 0;

    if(binary_)
    {
        showFloat fl;
        //std::cout<<"READVECTOR: "<<file.c_str()<<"\n";
        char buffer[80];
        is.read( buffer, 80);       //ignore first 80 characters

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            d0=fl.f;
            is.read( fl.s, sizeof(float));
            d1=fl.f;
            is.read( fl.s, sizeof(float));
            d2=fl.f;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
             end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            d0=fl.f;
            is.read( fl.s, sizeof(float));
            d1=fl.f;
            is.read( fl.s, sizeof(float));
            d2=fl.f;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }
    }
    else
    { // ASCII
        char buf[256];

        std::ifstream is( file.c_str());
        CheckFile( is);
        is.getline( buf, 256); // ignore first line

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is >> d0 >> d1 >> d2;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;

            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is >> d0 >> d1 >> d2;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;

            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }
    }

    CheckFile( is);
}

void ReadEnsightP2SolCL::CheckFile( const std::ifstream& is) const
{
    if (!is) throw DROPSErrCL( "ReadEnsightP2SolCL: error while reading from file!");
}


} // end of namespace DROPS

#endif
