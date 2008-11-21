//**************************************************************************
// File:    ensightOut.cpp                                                 *
// Content: solution output in Ensight6 Case format                        *
// Author:  Joerg Grande, Sven Gross, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************

#include "out/ensightOut.h"

namespace DROPS{

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
    const Uint lvl= _idx->TriangLevel(),
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
        std::memset(buffer, 0, 80);
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
        sInt.i= (int)idxDesc->NumUnknowns();            //write number of nodes
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

        std::strncpy( buffer, _geom.c_str(), 20);
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
         os << "coordinates\n" << std::setw(8) << idxDesc->NumUnknowns() << '\n';

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

Ensight6OutCL::Ensight6OutCL (std::string casefileName, Uint numsteps, bool binary)
    : decDigits_( 1), timestep_( -1), numsteps_( numsteps), time_( -1.), casefile_( casefileName), binary_( binary)
{
    while( numsteps>9) {
        ++decDigits_;
        numsteps/=10;
    }
}

Ensight6OutCL::~Ensight6OutCL ()
{
    for (std::map<std::string,Ensight6VariableCL*>::iterator it= vars_.begin(); it != vars_.end(); ++it)
        delete it->second;
}

void
Ensight6OutCL::CheckFile (const std::ofstream& os) const
{
    if (!os) throw DROPSErrCL( "EnsightP2SolOutCL: error while opening file!");
}

void
Ensight6OutCL::OpenFile (std::ofstream& os, std::string varName)
{
    std::string filename( vars_[varName]->fileName());
    if (vars_[varName]->Timedep())
         AppendTimecode( filename);
    os.open( filename.c_str());
    CheckFile( os);
}

void
Ensight6OutCL::AppendTimecode(std::string& str) const
{
    std::ostringstream postfix;
    postfix.width( decDigits_);
    postfix.fill ('0');
    postfix << timestep_;
    str+= postfix.str();
}

bool
Ensight6OutCL::putTime( double t)
{
    if (t > time_) {
        time_= t;
        if (timedep_) {
            timestr_ << t << ' ';
            if (++timestep_%10==0)
                timestr_ << "\n\t\t\t";
        }
        CommitCaseFile();    // rewrite case file
        return true;
    }
    return false;
}

void
Ensight6OutCL::CommitCaseFile ()
{
    // rewrite case file
    std::ofstream caseout( casefile_.c_str());
    CheckFile( caseout);
    caseout << "FORMAT\ntype: ensight\n\n"
            << geomdesc_.str()
            << "\n\nVARIABLE\n"
            << vardesc_.str() << '\n';
            
    if (timedep_) {
        caseout << "\nTIME\ntime set:\t\t1\nnumber of steps:\t" << timestep_ + 1
                << "\nfilename start number:\t0\nfilename increment:\t1\ntime values:\t\t";
        caseout << timestr_.str() << "\n\n";
    }
}

void
Ensight6OutCL::DescribeGeom (std::string varName)
{
    Ensight6VariableCL* geo= vars_[varName];
    std::string filename=  geo->fileName();

    geomdesc_ << "GEOMETRY\nmodel:\t\t\t";
    if (geo->Timedep())
    {
        timedep_= true;
        geomdesc_ << '1';
        filename+= std::string( decDigits_, '*');
    }
    geomdesc_ << "\t\t\t" << filename;
}

void
Ensight6OutCL::putGeom (MultiGridCL& mg, int lvl, std::string geoName)
{
    std::ofstream os;
    OpenFile( os, geoName);

    IdxDescCL p2idx( P2_FE);                              // Create a temporary Idx
    p2idx.CreateNumbering( lvl, mg);
    const size_t idx= p2idx.GetIdx();

    if (binary_)
    {
        char buffer[80];
        std::memset(buffer, 0, 80);
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

        sInt.i= (int)p2idx.NumUnknowns();                  //write number of nodes
        os.write( sInt.s, sizeof(int));

        DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {     //write node ids
            sInt.i= it->Unknowns( idx) + 1;
            os.write( sInt.s, sizeof(int));
        }
        DROPS_FOR_TRIANG_EDGE( mg, lvl, it) {      //write ids of additional nodes
            sInt.i= it->Unknowns( idx) + 1;
            os.write( sInt.s, sizeof(int));
        }

        DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {    //write coordinates of nodes
            Point3DCL c= it->GetCoord();
            for(int i=0; i<3; ++i) {
                sFlo.f=c[i];
                os.write( sFlo.s, sizeof(float));
            }
        }
        DROPS_FOR_TRIANG_EDGE( mg, lvl, it) {      //write coordinates of additional nodes
            Point3DCL c= GetBaryCenter( *it);
            for(int i=0; i<3; ++i) {
                sFlo.f=c[i];
                os.write( sFlo.s, sizeof(float));
            }
        }

        std::strcpy(buffer,"part 1");                          //part no. line
        os.write(buffer,80);

        std::strncpy( buffer, geoName.c_str(), 20);
        os.write(buffer,80);
        std::strcpy(buffer,"tetra4");                          // element 1 tetrahedron with 4 nodes each
        os.write(buffer,80);

        sInt.i =  8*std::distance( mg.GetTriangTetraBegin(lvl), mg.GetTriangTetraEnd(lvl));     //number of tetrahedra
        os.write( sInt.s, sizeof(int));

        RefRuleCL RegRef= GetRefRule( RegRefRuleC);
        DROPS_FOR_TRIANG_TETRA( mg, lvl, it) {
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
        os << "coordinates\n" << std::setw(8) << p2idx.NumUnknowns() << '\n';

        DROPS_FOR_TRIANG_VERTEX( mg, lvl, it) {
            os << std::setw(8) << it->Unknowns(idx)+1;
            Point3DCL c= it->GetCoord();
            for (int i=0; i<3; ++i)
                os << std::setw(12) << c[i];
            os << '\n';
        }
        DROPS_FOR_TRIANG_EDGE( mg, lvl, it) {
            os << std::setw(8) << it->Unknowns(idx)+1;
            Point3DCL c= GetBaryCenter( *it);
            for (int i=0; i<3; ++i)
                os << std::setw(12) << c[i];
            os << '\n';
        }

        os << "part 1\n" << geoName << "\n";
        // Ausgabe auf virtuell reg. verfeinertem Gitter, da Ensight zur Berechnung
        // der Isoflaechen anscheinend nur die Werte in den Vertices beruecksichtigt...
        os << "tetra4\n"
           << std::setw(8) << 8*std::distance(mg.GetTriangTetraBegin(lvl), mg.GetTriangTetraEnd(lvl)) << '\n';

        RefRuleCL RegRef= GetRefRule( RegRefRuleC);
        DROPS_FOR_TRIANG_TETRA( mg, lvl, it) {
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
    p2idx.DeleteNumbering( mg);
}

void Ensight6OutCL::DescribeVariable (std::string varName, bool isscalar)
{
    Ensight6VariableCL* var= vars_[varName];
    std::string filename=  var->fileName();

    vardesc_ << (isscalar ? "scalar" : "vector") << " per node:\t";
    if (var->Timedep())
    {
        timedep_= true;
        vardesc_ << '1';
        filename+= std::string( decDigits_, '*');
    }
    vardesc_ << '\t' << varName << '\t' << filename << std::endl;
}

void
Ensight6OutCL::Write (double t)
{
    if (!putTime( t)) return;
    for( std::map<std::string,Ensight6VariableCL*>::iterator it= vars_.begin(); it != vars_.end(); ++it) {
        if ( t== 0. || it->second->Timedep()) it->second->put( *this);
    }
}

void
Ensight6OutCL::Register (Ensight6VariableCL& var)
{
    vars_[var.varName()]= &var;
    var.Describe( *this);
}


void ReadEnsightP2SolCL::CheckFile( const std::ifstream& is) const
{
    if (!is) throw DROPSErrCL( "ReadEnsightP2SolCL: error while reading from file!");
}

} // end of namespace DROPS
