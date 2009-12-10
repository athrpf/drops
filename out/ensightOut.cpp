/// \file ensightOut.cpp
/// \brief solution output in Ensight6 Case format
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "out/ensightOut.h"

namespace DROPS{

Ensight6OutCL::Ensight6OutCL (std::string casefileName, Uint numsteps, bool binary, __UNUSED__ bool masterout)
    : decDigits_( 1), timestep_( -1), numsteps_( numsteps), time_( -1.),
      casefile_( casefileName), binary_( binary), timedep_( false), tag_( 2001), procDigits_( 1),
      nodes_( 0), masterout_( masterout), IamMaster_( true)
{
    while( numsteps>9) {
        ++decDigits_;
        numsteps/=10;
    }

#ifdef _PAR
    int procs= ProcCL::Size();
    while (procs > 9) {
        ++procDigits_;
        procs/=10;
    }
    IamMaster_= ProcCL::IamMaster();
#endif
}

Ensight6OutCL::~Ensight6OutCL ()
{
    for (std::map<std::string,Ensight6VariableCL*>::iterator it= vars_.begin(); it != vars_.end(); ++it)
        delete it->second;
}

void
Ensight6OutCL::CheckFile (const std::ofstream& os) const
{
    if (!os) throw DROPSErrCL( "Ensight6OutCL: error while opening file!");
}

void
Ensight6OutCL::OpenFile (std::ofstream& os, std::string varName, bool appendproccode)
{
    std::string filename( vars_[varName]->fileName());
    if (vars_[varName]->Timedep())
         AppendTimecode( filename);
    if (appendproccode)
        AppendProccode( filename);
    os.open( filename.c_str());
    CheckFile( os);
}

void
Ensight6OutCL::AppendTimecode (std::string& str) const
{
    std::ostringstream postfix;
    postfix.width( decDigits_);
    postfix.fill ('0');
    postfix << timestep_;
    str+= postfix.str();
}

void
Ensight6OutCL::AppendProccode ( __UNUSED__ std::string &str) const
{
#ifndef _PAR
    return;
#else
    std::ostringstream postfix;
    postfix.width( procDigits_);
    postfix.fill ('0');
    postfix << ProcCL::MyRank();
    str+= postfix.str();
#endif
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
    if (!IamMaster_)
        return;

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

#ifndef _PAR
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

#else

/// \brief Write part of the geometry into file (filename)
/** A master process writes all coordinates of the GIDS into the file.
    If not masterout this works like this:
    Every process writes a geometry-part-file filename(proc-id). In this file
    all tetras are listed. To gain a complete ensight geometry file, paste all
    geometry-files together or use the program "pardrops2ens"
    If masterout is set:
    All processor send the connectivity part to a master processor. This proc
    writes these information into a valid ensight format and no post-processing
    has to be done.*/
void Ensight6OutCL::putGeom (MultiGridCL& mg, int lvl, std::string geoName)
{
    // create two list to buffer received coords and points
    std::queue <showInt> bufferedIDs;
    std::queue <showFloat> bufferedCoors;

    const int  me    = ProcCL::MyRank(),
    master= ProcCL::Master();

    // Collect coordinates and gids of the global nodes
    //---------------------------------------------------
    // number of vertices and edges written by this proc (= number of vals)
    const IdxT locNumUnknowns = GetExclusiveVerts(mg, PrioHasUnk, lvl) + GetExclusiveEdges(mg, PrioHasUnk, lvl);
    const IdxT globNumUnknowns = ProcCL::GlobalSum(locNumUnknowns);                         // global number of vertices
    DDD_GID *myGID   = new DDD_GID[locNumUnknowns];                                         // array of gids of exclusive vertices, owned by thiss proc
    double *myCoord = new double[3*locNumUnknowns];                                         // and the corresponding coords

    nodes_ = locNumUnknowns;
    Uint pos_GID=0, pos_Coord=0;                                                            // actual position in the above arrays
    for (MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),      // iterate over all
         end= mg.GetTriangVertexEnd(lvl); it!=end; ++it)                                  // vertices of given level
    {
        if (it->IsExclusive(PrioHasUnk))                                                    // just collect information, if the vert is exclusive
        {
            if (pos_GID>=locNumUnknowns || pos_Coord>=3*locNumUnknowns)                     // error checking
                std::cout << "["<<ProcCL::MyRank()<<"] putGeom (Vertex): "
                        << "not enough memory allocated (pos_GID " <<pos_GID<<"/"
                        <<locNumUnknowns<< ")!" << std::endl;

            Point3DCL c=it->GetCoord();                                                     // get coordinates
            myGID[pos_GID] = it->GetGID();                                                  // put gid into the array
            for (int i=0; i<3; ++i)                                                         // put coordinates into array
                myCoord[pos_Coord+i] = c[i];
            pos_GID   += 1;                                                                 // increase actual pos in arrays
            pos_Coord += 3;
        }
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(lvl),          // the same es for vertices
         end= mg.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        if (it->IsExclusive(PrioHasUnk))
        {
            if (pos_GID>locNumUnknowns || pos_Coord>3*locNumUnknowns)
                std::cout << "["<<ProcCL::MyRank()<<"] putGeom (Edge): "
                        << "not enough memory allocated (pos_GID " <<pos_GID<<"/"
                        <<locNumUnknowns<< ")!" << std::endl;

            Point3DCL c= GetBaryCenter( *it);
            myGID[pos_GID] = it->GetGID();
            for (int i=0; i<3; ++i)
                myCoord[pos_Coord+i] = c[i];
            pos_GID   += 1;
            pos_Coord += 3;
        }
    }

    std::ofstream* os= 0;

    showInt sInt;                  //unions for converting ASCII int to binary
    char buffer[80];               // is used in binary-mode
    std::memset(buffer, 0, 80);

    // Write out the collected information about vertices and coordinates
    //----------------------------------------------------------------------
    if (me!=master)                                                                     // all procs send their information to the master
    {                                                                                   // split message. Should be no performance problem
        ProcCL::RequestT req[2];                                                        // request for Isend
        req[0]= ProcCL::Isend(myGID, locNumUnknowns, master, tag_);                     // send the gids
        req[1]= ProcCL::Isend(myCoord, 3*locNumUnknowns, master, tag_+1);               // send the coordinates
        ProcCL::WaitAll(2,req);                                                         // Wait until all messages are send
        delete[] myGID;                                                                 // the memory is not used any more, so free it
        delete[] myCoord;
    }
    else                                                                                // the master writes out all information
    {
        os = new std::ofstream;
        OpenFile( *os, geoName, /*append proc-code*/ !masterout_);
        os->flags(std::ios_base::scientific);                                           // format as demanded by ensight
        os->precision(5);
        os->width(12);

        // binary output
        if(binary_) {

            std::strcpy(buffer,"C Binary");     //writing of all necessary information: Binary header
            os->write(buffer,80);
            std::strcpy(buffer,"DROPS geometry file: ");
            os->write(buffer,80);
            std::strcpy(buffer,"format: Ensight6 Case format");
            os->write(buffer,80);
            std::strcpy(buffer,"node id given");
            os->write(buffer,80);
            std::strcpy(buffer,"element id off");
            os->write(buffer,80);
            std::strcpy(buffer,"coordinates");
            os->write(buffer,80);

            sInt.i= (int) globNumUnknowns;            //write number of nodes
            os->write( sInt.s, sizeof(int));
        }
        // ascii output
        else
        {
            (*os) << "DROPS geometry file: " << "\nformat: Ensight6 Case format\n"
                  << "node id given\nelement id off\n"
                  << "coordinates\n" << std::setw(8) << globNumUnknowns << '\n';
        } // end ascii

        unsigned int recvUnknowns = 0 ;

        for (int p=0; p<ProcCL::Size(); ++p)                                            // print the gids with coords of all procs into the file
        {                                                                               // it is important that the information of proc 0 are print first, then of proc 1 and so on ( because of the ordering of the printed vals)
            IdxT numUnk;                                                                // number of exclusive verts of the other proc
            DDD_GID *Gids=0;                                                            // receive-buffer for gids
            double *coord=0;                                                            // receive-biffer for coords
            if (p!=me)  // ==> p!=master
            {                                                                           // receive gids and coords
                ProcCL::StatusT stat;
                ProcCL::Probe(p, tag_, stat);
                numUnk = ProcCL::GetCount<DDD_GID>(stat);

                Gids  = new DDD_GID[numUnk];
                coord = new double[3*numUnk];
                ProcCL::Recv(Gids, numUnk, p, tag_);
                ProcCL::Recv(coord, 3*numUnk, p, tag_+1);
            }
            else        // p==master==me
            {                                                                           // reciev-pointer == local gid/coord pointer
                numUnk = locNumUnknowns;
                Gids  = myGID;
                coord = myCoord;
            }

            if(binary_) {
                for (size_t i=0; i<numUnk; ++i) {
                        sInt.i= (int) Gids[i];
                        bufferedIDs.push(sInt);
                }
            }

            for (size_t i=0; i<numUnk; ++i)
            {
                if(binary_) {
                    showFloat sFlo;
                    for (int j=0; j<3; ++j) {
                        sFlo.f= coord[3*i+j];
                        bufferedCoors.push(sFlo);
                    }
                }
                else {
                    (*os) << std::setw(8) << Gids[i];                                          // print gid into file
                    for (int j=0; j<3; ++j)                                                 // print coord into file
                        (*os) << std::setw(12) << coord[3*i+j];
                    (*os) << std::endl;
                }
            }
            if (Gids!=0)                                                                // free memory (master frees now local array!)
            {
                delete[] Gids;
                Gids=0;
            }
            if (coord!=0)
            {
                delete[] coord;
                coord=0;
            }

            recvUnknowns += numUnk;
        }

        if(binary_ && recvUnknowns == globNumUnknowns )  {

            while(bufferedIDs.size()>0) {
                showInt tmp = bufferedIDs.front();
                os->write( tmp.s, sizeof(int));
                bufferedIDs.pop();
            }

            while(bufferedCoors.size()>0) {
                showFloat tmp = bufferedCoors.front();
                os->write( tmp.s, sizeof(float));
                //std::cout<<tmp.s <<"\n";
                bufferedCoors.pop();
            }

        }
    }

    if (masterout_)                 // in this mode the master writes out all tetrahedrons
    {
        // at first collect all information in an array tetras of length numTetras
        Uint numTetras=std::distance(mg.GetTriangTetraBegin(lvl),mg.GetTriangTetraEnd(lvl));
        Uint numAllTetra= ProcCL::GlobalSum(numTetras, master);
        Uint numMaxTetra= ProcCL::GlobalMax(numTetras, master);
        DDD_GID *tetras= new DDD_GID[8*numTetras*4];          // regrefine*numTetra*vertices
        Uint pos=0;

        for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl),    // for all tetras
            end= mg.GetTriangTetraEnd(lvl); it!=end; ++it)
        {
            RefRuleCL RegRef= GetRefRule( RegRefRuleC);                                     // get regular refine rule
            for (int ch=0; ch<8; ++ch)                                                      // iterate over all children
            {
                ChildDataCL data= GetChildData( RegRef.Children[ch]);                       // datas of children
                for (int vert= 0; vert<4; ++vert)
                {
                    int v= data.Vertices[vert];                                             // number of corresponding vert

                    if (v<4)                                                                // number correspond to original vertex
                        tetras[pos++] = it->GetVertex(v)->GetGID();
                    else                                                                    // number correspond to an edge
                        tetras[pos++] = it->GetEdge(v-4)->GetGID();
                }
            }
        }

        if (pos!=8*numTetras*4){
            throw DROPSErrCL("Ensight6OutCL:putGeom: Number of tetra information does not match");
        }

        if (me!=master)
        {
            // If not master, send collected information to master
            ProcCL::RequestT req;
            req = ProcCL::Isend(tetras, 8*numTetras*4, master, tag_);
            ProcCL::Wait(req);         // before deleting tetras, wait for finishing the transmition
            delete[] tetras;
        }
        else
        {
            // If master, write all information to file
            if(binary_) {
                    std::strcpy(buffer,"part 1");
                    os->write(buffer,80);
                    std::strcpy(buffer,"Messzelle");
                    os->write(buffer,80);
                    std::strcpy(buffer,"tetra4");
                    os->write(buffer,80);
                    sInt.i= (int) (8*numAllTetra);
                    os->write( sInt.s, sizeof(int));
            }
            else {
                (*os) << "part 1\nMesszelle\ntetra4\n"
                      << std::setw(8) << (8*numAllTetra)<<std::endl;
            }
            // this is the master proc. He has to write all information!
            DDD_GID *Gids= new DDD_GID[8*numMaxTetra*4];  // storage for all gids of one proc
            Uint numGids;                                 // number of all gids

            // for each process write out his tetras
            for (int p=0; p<ProcCL::Size(); ++p)
            {
                if (p!=me)
                {
                    ProcCL::StatusT stat;
                    ProcCL::Probe(p, tag_, stat);
                    numGids = (Uint)ProcCL::GetCount<DDD_GID>(stat);
                    ProcCL::Recv(Gids, numGids, p, tag_);
                }
                else
                {
                    std::swap(Gids,tetras);
                    numGids=8*numTetras*4;
                }

                pos=0;
                for (Uint i=0; i<numGids; i+=4)
                {
                    if(binary_){
                        for (int v=0; v<4; ++v) {
                            sInt.i= (int) Gids[i+v];
                            os->write( sInt.s, sizeof(int));
                        }

                    }
                    else{
                        for (int v=0; v<4; ++v)
                            (*os) << std::setw(8) << Gids[i+v];
                            (*os) << std::endl;
                    }
                }
                if (me==p)
                    std::swap(Gids,tetras);
            }
            delete[] Gids;
            delete[] tetras;
        }
    }

    if (os){
        os->close();
        delete os;
    }

    // Each proc writes out its own file
    if (!masterout_)
    {
        if (binary_)
            throw DROPSErrCL("Ensight6OutCL:putGeom: Non-Masterout and binary has not been implemented yet, sorry!");

        // Write out the parts (each proc writes out all his tetrahedrons of given level in a single file)
        //---------------------------------------------------------------------------------------------
        std::ofstream osp;
        OpenFile( osp, geoName, /*append proc-code*/true);
        osp.flags(std::ios_base::scientific);                                                // format demanded by ensight
        osp.precision(5);
        osp.width(12);

        osp << "part "<<(me+1)<<"\n" << geoName << "\n";                                       // part number is procid+1
        osp << "tetra4\n" << std::setw(8)
                << 8*std::distance(mg.GetTriangTetraBegin(lvl),mg.GetTriangTetraEnd(lvl))
                << '\n';
        // Ausgabe auf virtuell reg. verfeinertem Gitter, da Ensight zur Berechnung
        // der Isoflaechen anscheinend nur die Werte in den Vertices beruecksichtigt...
        for (MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(lvl),    // for all tetras
            end= mg.GetTriangTetraEnd(lvl); it!=end; ++it)
        {
            RefRuleCL RegRef= GetRefRule( RegRefRuleC);                                     // get regular refine rule
            for (int ch=0; ch<8; ++ch)                                                      // iterate over all children
            {
                ChildDataCL data= GetChildData( RegRef.Children[ch]);                       // datas of children
                for (int vert= 0; vert<4; ++vert)
                {
                    int v= data.Vertices[vert];                                             // number of corresponding vert
                    osp << std::setw(8);
                    if (v<4)                                                                // number correspond to original vertex
                        osp << it->GetVertex(v)->GetGID();
                    else                                                                    // number correspond to an edge
                        osp << it->GetEdge(v-4)->GetGID();
                }
                osp << '\n';                                                                 // go on to the next tetra
            }
        }
        osp.close();
    }
}
#endif

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
        if (!it->second->is_geom()) continue;
        if ( t== 0. || it->second->Timedep()) {
            // std::cerr << "Writing " << it->second->varName() << '\n';
            it->second->put( *this);
        }
    }
    for( std::map<std::string,Ensight6VariableCL*>::iterator it= vars_.begin(); it != vars_.end(); ++it) {
        if (it->second->is_geom()) continue;
        if ( t== 0. || it->second->Timedep()) {
            // std::cerr << "Writing " << it->second->varName() << '\n';
            it->second->put( *this);
        }
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

#ifdef _PAR           // parallel implementations

/// \brief Put procnumber behind a filename
void ReadEnsightP2SolCL::AppendProccode( std::string &str) const
{
    char format[]= ".%0Xi", postfix[8];
    format[3]= '0' + char(_procDigits);
    std::sprintf( postfix, format, ProcCL::MyRank());
    str+= postfix;
}

#endif          // end of parallel implementations

} // end of namespace DROPS
