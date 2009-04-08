//**************************************************************************
// File:    TestExchangePar.cpp                                            *
// Content: Testprogram für den Austausch von Werten                       *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, RZ RWTH Aachen                               *
// Version: 0.1                                                            *
// Begin:   March, 08th 2006                                               *
// mod:     0.2 Test of ExchangeBlockCL (29th May 2006)                    *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestExchangePar.cpp
/// \brief Testing the parallel exchange of numerical data for accumulations

 // include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include "parallel/exchange.h"
#include "partests/params.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "num/spmat.h"
#include "levelset/levelset.h"

// Standard-Header-Files für Ausgaben
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;
/****************************************************************************
* G L O B A L E   V A R I A B L E N                                         *
****************************************************************************/
DROPS::NoBndDataCL<double> Bnd;                     // Keine Randdaten
const int REF=0, MIG=1;                             // Unterscheidung, welche Datei geschrieben wird
DROPS::ParamParExchangeCL C;                        // Parameter aus dem Parameterfile
DROPS::TimeStoreCL Times(10);                       // Klasse zum Speichern der Zeiten
enum TimePart{                                      // Aufzählung, welche Zeiten gemessen werden
    Refine,
    Ex_Create,
    Ex_Create_Map,
    Ex_Create_Acc,
    Ex_Acc,
    Ex_Acc_newInterface,
    Ex_Acc_new,
    Ex_Acc_new_newInterface,
    CheckMG
};
const char line[] ="-----------------------------------------------------------------------------";

/****************************************************************************
    * S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(Refine,                  "Refinement");
    Times.SetDescriber(Ex_Create,               "100 Create Exchange list, pure");
    Times.SetDescriber(Ex_Create_Map,           "100 Create Exchange list, with index mapping");
    Times.SetDescriber(Ex_Create_Acc,           "100 Create Exchange list, with accumulated indices");
    Times.SetDescriber(Ex_Acc,                  "100 inner products with ExchangeCL (old interface)");
    Times.SetDescriber(Ex_Acc_newInterface,     "100 inner products with ExchangeCL (new interface)");
    Times.SetDescriber(Ex_Acc_new,              "100 inner products with new (accurater) ExchangeCL (old interface)");
    Times.SetDescriber(Ex_Acc_new_newInterface, "100 inner products with new (accurater) ExchangeCL (new interface)");
    Times.SetDescriber(CheckMG,                 "Checking MG");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
*****************************************************************************
*   Checkt, ob die parallele Verteilung und die MultiGrid-Struktur gesund   *
*   zu sein scheint.                                                        *
****************************************************************************/
bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg, ostream& os)
{
    DROPS::ParTimerCL time;
    bool pmg_sane = pmg.IsSane(os),
    mg_sane  = pmg.GetMG().IsSane(os);
    time.Stop(); Times.AddTime(CheckMG,time.GetMaxTime());
    return DROPS::ProcCL::Check(pmg_sane && mg_sane);
}

void PrintExchange(const DROPS::ExchangeCL& ex)
{
    const int me =DROPS::ProcCL::MyRank();
    static int ExNum=0;
    char filename[30];
    sprintf(filename, "output/%i_Exchange_%i.mg",me,ExNum++);
    ofstream file(filename);
    ex.DebugInfo(file);
    file.close();
    if (me==0)
        cout << "    --> Schreibe in Datei: " << filename << "  \n";
}

void DoRefinement(DROPS::ParMultiGridCL &pmg)
{
    DROPS::ParTimerCL time;
    if (DROPS::ProcCL::IamMaster())
        std::cout <<line<<std::endl<< " * Verteile Multigrid und verfeinere " <<C.refall<< " Mal regulär" << std::endl;
    DROPS::MultiGridCL &mg =pmg.GetMG();
    DROPS::LoadBalHandlerCL lb(mg);
    lb.DoInitDistribution(DROPS::ProcCL::Master());
    ::Times.IncCounter(lb.GetMovedMultiNodes());
    switch (C.refineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
    }

    for (int ref=0; ref<C.refall; ++ref)
    {
        DROPS::MarkAll(mg);
        time.Reset();
        pmg.Refine();
        time.Stop(); ::Times.AddTime(Refine,time.GetMaxTime());
        lb.DoMigration();
        ::Times.IncCounter(lb.GetMovedMultiNodes());
    }
    if (C.checkMG)
        CheckParMultiGrid(pmg,std::cout);
    time.Stop(); Times.AddTime(Refine, time.GetMaxTime());
}

namespace DROPS // für die Strategy und Hilfsfunktionen
{
/****************************************************************************
* S E T  N U M                                                              *
*****************************************************************************
*   Diese Funktion iteriert über alle Vertices, Edges und Tetras, auf denen *
*   Unbekannte gesetzt werden können und gibt Simplex zum ersten Index den  *
*   Wert GID+i und de zweiten Index den Wert 1/100(GID+i)+0.1, damit man    *
*   nach dem Akkumulieren die Werte checken kann.                           *
****************************************************************************/
void SetNum(MultiGridCL& mg, VecDescCL& x1, VecDescCL& x2, VecDescCL& x3)
{
    // Setze die GID als Werte auf die Vertices:
    for (MultiGridCL::TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()), end(mg.GetTriangVertexEnd()) ; sit!=end; ++sit)
    {
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x1.RowIdx->GetIdx()) )
            for (Uint i=0; i<x1.RowIdx->NumUnknownsVertex(); ++i)
                x1.Data[ sit->Unknowns(x1.RowIdx->GetIdx())+i ] = (double)(sit->GetGID()+i);
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x2.RowIdx->GetIdx()) )
            for (Uint i=0; i<x2.RowIdx->NumUnknownsVertex(); ++i)
                x2.Data[ sit->Unknowns(x2.RowIdx->GetIdx())+i ] = 1./100.*(sit->GetGID()+i)+0.1;
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x3.RowIdx->GetIdx()) ) {
            const IdxT dof= sit->Unknowns(x3.RowIdx->GetIdx());
            for (Uint i=0; i<x3.RowIdx->NumUnknownsVertex(); ++i) {
                x3.Data[ dof+i ] = (sit->GetGID()+i)+0.11;
                if (x3.RowIdx->IsExtended(dof))
                     x3.Data[ x3.RowIdx->GetXidx()[dof]+i ] = (sit->GetGID()+i)+0.55;
            }
        }
    }

    // Setze die GID als Werte auf die Edges:
    for (MultiGridCL::TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()), end(mg.GetTriangEdgeEnd()) ; sit!=end; ++sit)
    {
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x1.RowIdx->GetIdx()) )
            for (Uint i=0; i<x1.RowIdx->NumUnknownsEdge(); ++i)
                x1.Data[ sit->Unknowns(x1.RowIdx->GetIdx())+i ] = (double)(sit->GetGID()+i);
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x2.RowIdx->GetIdx()) )
            for (Uint i=0; i<x2.RowIdx->NumUnknownsEdge(); ++i)
                x2.Data[ sit->Unknowns(x2.RowIdx->GetIdx())+i ] = 1./100.*(sit->GetGID()+i)+0.1;
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x3.RowIdx->GetIdx()) ) {
            const IdxT dof= sit->Unknowns(x3.RowIdx->GetIdx());
            for (Uint i=0; i<x3.RowIdx->NumUnknownsEdge(); ++i) {
                x3.Data[ dof+i ] = (sit->GetGID()+i)+0.11;
                if (x3.RowIdx->IsExtended(dof))
                     x3.Data[ x3.RowIdx->GetXidx()[dof]+i ] = (sit->GetGID()+i)+0.55;
            }
        }
    }

    // Setze die GID als Werte auf die Tetras:
    for (MultiGridCL::TriangTetraIteratorCL sit(mg.GetTriangTetraBegin()), end(mg.GetTriangTetraEnd()) ; sit!=end; ++sit)
    {
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x1.RowIdx->GetIdx()) )
            for (Uint i=0; i<x1.RowIdx->NumUnknownsTetra(); ++i)
                x1.Data[ sit->Unknowns(x1.RowIdx->GetIdx())+i ] = (double)(sit->GetGID()+i);
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x2.RowIdx->GetIdx()) )
            for (Uint i=0; i<x2.RowIdx->NumUnknownsTetra(); ++i)
                x2.Data[ sit->Unknowns(x2.RowIdx->GetIdx())+i ] = 1./100.*(sit->GetGID()+i)+0.1;
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x3.RowIdx->GetIdx()) ) {
            const IdxT dof= sit->Unknowns(x3.RowIdx->GetIdx());
            for (Uint i=0; i<x3.RowIdx->NumUnknownsTetra(); ++i) {
                x3.Data[ dof+i ] = (sit->GetGID()+i)+0.11;
                if (x3.RowIdx->IsExtended(dof))
                     x3.Data[ x3.RowIdx->GetXidx()[dof]+i ] = (sit->GetGID()+i)+0.55;
            }
        }
    }
}


inline bool CheckValue(double ist, double soll){
    if (std::fabs(soll-ist)<DoubleEpsC*std::fabs(soll))
        return true;
    return false;
}

/****************************************************************************
* C H E C K  N U M                                                          *
*****************************************************************************
*   Diese Funktion schaut nach dem Akkumulieren nach, ob alle Werte auf den *
*   Simplices richtig ist.                                                  *
****************************************************************************/
bool CheckNum(MultiGridCL& mg, VecDescCL& x1, VecDescCL& x2, VecDescCL& x3)
{
    bool check = true;
//  const int me =ProcCL::MyRank();
    // Check vertices
    for (MultiGridCL::TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()), end(mg.GetTriangVertexEnd()) ; sit!=end; ++sit)
    {
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x1.RowIdx->GetIdx()) )
        {
            bool local_check = true;
            const DDD_HDR vert= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof= sit->Unknowns(x1.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(vert); *procList!=-1; procList+=2)
                if (procList[1]==PrioHasUnk)
                    ++num_shared;
            for (Uint i=0; i<x1.RowIdx->NumUnknownsVertex(); ++i)
                local_check = local_check && CheckValue( x1.Data[dof+i], num_shared*sit->GetGID()+i);
            check = check && local_check;
            if (!local_check)
            {
                std::cout << "["<<ProcCL::MyRank()<<"] Index "<<x2.RowIdx->GetIdx()<<": Bei Knoten " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << x1.Data[dof]
                        << " sollte sein: " << (num_shared*(double)(sit->GetGID()))
                        << ", Knoten hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x2.RowIdx->GetIdx()) )
        {
            bool local_check = true;
            const DDD_HDR vert= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof=  sit->Unknowns(x2.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(vert); *procList!=-1; procList+=2)
                if (procList[1]==PrioHasUnk)
                    ++num_shared;
            for (Uint i=0; i<x2.RowIdx->NumUnknownsVertex(); ++i)
                local_check = local_check && CheckValue( x2.Data[dof+i], num_shared*(1./100.*(sit->GetGID()+i)+0.1 ));
            check = check && local_check;
            if (!local_check)
            {
                std::cout << "["<<ProcCL::MyRank()<<"] Index "<<x2.RowIdx->GetIdx()<<": Bei Knoten " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << x2.Data[ dof]
                        << " sollte sein: " << (num_shared*(1./100.*(sit->GetGID())+0.1))
                        << ", Knoten hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }

        if ( sit->Unknowns.Exist() && x3.RowIdx && sit->Unknowns.Exist(x3.RowIdx->GetIdx()) )
        {
            bool local_check = true, local_xcheck= true;
            const DDD_HDR vert= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof= sit->Unknowns(x3.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(vert); *procList!=-1; procList+=2)
                if (procList[1]==PrioHasUnk)
                    ++num_shared;
            for (Uint i=0; i<x3.RowIdx->NumUnknownsVertex(); ++i)
                local_check = local_check && CheckValue( x3.Data[dof+i], num_shared*((sit->GetGID()+i)+0.11));
            if (x3.RowIdx->IsExtended(dof)) {
                const IdxT xdof= x3.RowIdx->GetXidx()[dof];
                for (Uint i=0; i<x3.RowIdx->NumUnknownsVertex(); ++i)
                    local_xcheck = local_xcheck && CheckValue( x3.Data[xdof+i], num_shared*((sit->GetGID()+i)+0.55 ));
            }
            check = check && local_check && local_xcheck;
            if (!local_check)
            {
                std::cout << "["<<ProcCL::MyRank()<<"] Index "<<x3.RowIdx->GetIdx()<<": Bei Knoten " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << (x3.Data[ dof])
                        << " sollte sein: " << (num_shared*((sit->GetGID())+0.11))
                        << ", Knoten hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
            if (!local_xcheck)
            {
                std::cout << "["<<ProcCL::MyRank()<<"] Index "<<x3.RowIdx->GetIdx()<<": Bei Knoten " << sit->GetGID() << " und extended DoF ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << (x3.Data[ x3.RowIdx->GetXidx()[dof]])
                        << " sollte sein: " << (num_shared*((sit->GetGID())+0.55))
                        << ", Knoten hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }
    }
    // check edges
    for (MultiGridCL::TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()), end(mg.GetTriangEdgeEnd()) ; sit!=end; ++sit)
    {
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x1.RowIdx->GetIdx()) )
        {
            bool local_check = true;
            const DDD_HDR edge= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof= sit->Unknowns(x1.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(edge); *procList!=-1; procList+=2)
                if (procList[1]==PrioHasUnk)
                    ++num_shared;
            for (Uint i=0; i<x1.RowIdx->NumUnknownsEdge(); ++i)
                local_check = local_check && CheckValue( x1.Data[dof+i], num_shared*(double)(sit->GetGID()+i));
            check = check && local_check;
            if (!local_check)
            {
                std::cout << " Index 1: Bei Ecke " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << (x1.Data[dof])
                        << " sollte sein: " << (num_shared*(double)(sit->GetGID()))
                        << ", Ecke hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x2.RowIdx->GetIdx()) )
        {
            bool local_check = true;
            const DDD_HDR edge= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof= sit->Unknowns(x2.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(edge); *procList!=-1; procList+=2)
                if (procList[1]==PrioHasUnk)
                    ++num_shared;
            for (Uint i=0; i<x2.RowIdx->NumUnknownsEdge(); ++i)
                local_check = local_check && CheckValue( x2.Data[dof+i], num_shared*(1./100.*(sit->GetGID()+i)+0.1));
            check = check && local_check;
            if (!local_check)
            {
                std::cout << " Index 2:Bei Ecke " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << (x2.Data[dof])
                        << " sollte sein: " << (num_shared*(1./100.*(sit->GetGID())+0.1))
                        << ", Ecke hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }
    }
    // check tetras
    for (MultiGridCL::TriangTetraIteratorCL sit(mg.GetTriangTetraBegin()), end(mg.GetTriangTetraEnd()) ; sit!=end; ++sit)
    {
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x1.RowIdx->GetIdx()) )
        {
            bool local_check = true;
            const DDD_HDR tetra= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof= sit->Unknowns(x1.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(tetra); *procList!=-1; procList+=2)
                if (procList[1]==PrioMaster)
                    ++num_shared;
            for (Uint i=0; i<x1.RowIdx->NumUnknownsTetra(); ++i)
                local_check = local_check && CheckValue( x1.Data[dof+i], num_shared*(double)(sit->GetGID()+i));
            check = check && local_check;
            if (!local_check)
            {
                std::cout << " Index 1: Bei Tetra " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << (x1.Data[dof])
                        << " sollte sein: " << (num_shared*(double)(sit->GetGID()))
                        << ", Tetra hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }
        if ( sit->Unknowns.Exist() && sit->Unknowns.Exist(x2.RowIdx->GetIdx()) )
        {
            bool local_check = true;
            const DDD_HDR tetra= sit->GetHdr();
            int num_shared = 0;
            const IdxT dof= sit->Unknowns(x2.RowIdx->GetIdx());
            for (int *procList = DDD_InfoProcList(tetra); *procList!=-1; procList+=2)
                if (procList[1]==PrioMaster)
                    ++num_shared;
            for (Uint i=0; i<x3.RowIdx->NumUnknownsTetra(); ++i)
                local_check = local_check && CheckValue( x2.Data[dof+i], num_shared*(1./100.*(sit->GetGID()+i)+0.1));
            check = check && local_check;
            if (!local_check)
            {
                std::cout << " Index 2: Bei Tetra " << sit->GetGID() << " ist ein Fehler passiert!\n"
                        << "  Habe den Wert: " << (x2.Data[dof])
                        << " sollte sein: " << (num_shared*(1./100.*(sit->GetGID()) +0.1))
                        << ", Tetra hat auf "<<num_shared<<" Prozessoren einen Wert"<< std::endl;
            }
        }
    }
    return ProcCL::Check(check);
}

/****************************************************************************
* C H E C K  I D X  M A P P I N G                                           *
*****************************************************************************
*   Diese Funktion überprüft, ob die ExchangeCL die Indices über Prozessor- *
*   ränder richtig berechnen kann. Dazu werden einem Prozessor (\p proc)    *
*   die verteilten Indices zugesendet. Dieser checkt dann, ob er die        *
*   gleichen berechnet hat.                                                 *
*****************************************************************************
*   Hilfsfunktionen                                                         *
*   - bool IsDist(Simplex,proc): returns if the Simplex has distributed     *
*     datas.                                                                *
*   -Uint IsDistNum(Simplex): returns number of procs, that has unknowns    *
*    on it.                                                                 *
****************************************************************************/
template<typename SimplexT>
  bool IsDist(const SimplexT& sim, int proc)
{
    if (sim.IsLocal())
        return false;
    int *procs = sim.GetProcList();

    if (sim.GetPrio()!=DROPS::PrioHasUnk)
        return false;

    for (int i=0; procs[i]!=-1; i+=2)
        if (procs[i] == proc && procs[i+1] == DROPS::PrioHasUnk)
            return true;

    return false;
}

template<typename SimplexT>
 DROPS::Uint DistNum(const SimplexT& sim)
{
    DROPS::Uint distnum=0;
    for (int* procs = sim.GetProcList(); *procs!=-1; procs+=2)
        if ( procs[1] == DROPS::PrioHasUnk)
            ++distnum;
    return distnum;
}

/// \brief Checks the mapping of one processor \a proc
bool CheckIdxMapping(const DROPS::MultiGridCL& mg, DROPS::MLIdxDescCL* idxDesc, const DROPS::ExchangeCL& ex, const int proc)
{
    DROPS::Uint idx=idxDesc->GetIdx();
    const int tag=1051;

    if (DROPS::ProcCL::MyRank()!=proc)
    {
        // Einsammeln der Information (GID, Sysnummer)
        std::vector<DROPS::Uint> GIDs;
        std::vector<DROPS::Ulint> Sysnums;
        for (DROPS::MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(),
             end= mg.GetTriangVertexEnd(); it!=end; ++it)
        {
            if (!IsDist(*it, proc))
                continue;
            GIDs.push_back(it->GetGID());
            Sysnums.push_back( it->Unknowns(idx) );
        }
        for (DROPS::MultiGridCL::const_TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(),
             end= mg.GetTriangEdgeEnd(); it!=end; ++it)
        {
            if (!IsDist(*it,proc))
                continue;
            GIDs.push_back(it->GetGID());
            Sysnums.push_back(it->Unknowns(idx) );
        }
        int count = GIDs.size();
        DROPS::Uint *SendBuf_GID = new DROPS::Uint[count];
        DROPS::Ulint *SendBuf_Sys= new DROPS::Ulint[count];
        for (int i=0; i<count; ++i)
        {
            SendBuf_GID[i]=GIDs[i];
            SendBuf_Sys[i]=Sysnums[i];
        }
        /// \todo (of) ISend verwenden
        ProcCL::Send(SendBuf_GID, count, proc, tag);
        ProcCL::Send(SendBuf_Sys, count, proc, tag+1);
        delete[] SendBuf_GID;
        delete[] SendBuf_Sys;
    }
    else
    {   // I am the testing proc
        for (int p=0; p<DROPS::ProcCL::Size(); ++p){
            if (p==proc)
                continue;
            std::map<DROPS::Uint, DROPS::Ulint> GID_Sys_Map;

            ProcCL::StatusT stat;
            ProcCL::Probe(p, tag, stat);
            int count = ProcCL::GetCount<Uint>(stat);

            DROPS::Uint  *RecvBuf_GID = new DROPS::Uint[count];
            DROPS::Ulint *RecvBuf_Sys = new DROPS::Ulint[count];
            ProcCL::Recv(RecvBuf_GID, count, p, tag);
            ProcCL::Recv(RecvBuf_Sys, count, p, tag+1);

            for (int i=0; i<count; ++i)
                GID_Sys_Map[RecvBuf_GID[i]] = RecvBuf_Sys[i];

            for (DROPS::MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(),
                 end= mg.GetTriangVertexEnd(); it!=end; ++it)
            {
                if (it->Unknowns.Exist() && it->Unknowns.Exist(idx)){
                    DROPS::Ulint sysnum=it->Unknowns(idx);
                    if (DistNum(*it)!=ex.GetNumProcs(sysnum)){
                        std::cout << "["<<ProcCL::MyRank()<<"] DistNum="<<DistNum(*it)<<", ExCL="<<ex.GetProcs(sysnum).size()+1<<std::endl;
                        throw DROPS::DROPSErrCL("Falsche ProcDimension!");
                    }
                    if (!IsDist(*it, p))
                        continue;
                    const DROPS::IdxT extSysnum=GID_Sys_Map[it->GetGID()];
                    for (Uint i=0; i<idxDesc->NumUnknownsVertex(); ++i)
                        if (extSysnum+i != ex.GetExternalIdxFromProc(sysnum+i,p)){
                            std::cout <<"["<<ProcCL::MyRank()<<"] Falsche Zuordnung (i="<<i<<"), Vert "<<it->GetGID()
                                    <<" orig Sysnum= "<<(extSysnum+i)<<", calcSysnum="<<ex.GetExternalIdxFromProc(sysnum+i,p)
                                    <<std::endl;
                            throw DROPS::DROPSErrCL("Falsche Zuordnung!");
                        }
                }
            }
            for (DROPS::MultiGridCL::const_TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(),
                 end= mg.GetTriangEdgeEnd(); it!=end; ++it)
            {
                if (it->Unknowns.Exist() && it->Unknowns.Exist(idx)){
                    DROPS::Ulint sysnum=it->Unknowns(idx);
                    if (DistNum(*it)!=ex.GetNumProcs(sysnum)){
                        std::cout << "["<<ProcCL::MyRank()<<"] DistNum="<<DistNum(*it)<<", ExCL="<<ex.GetProcs(sysnum).size()+1<<std::endl;
                        throw DROPS::DROPSErrCL("Falsche ProcDimension!");
                    }

                    if (!IsDist(*it, p))
                        continue;
                    const DROPS::IdxT extSysnum=GID_Sys_Map[it->GetGID()];
                    for (Uint i=0; i<idxDesc->NumUnknownsEdge(); ++i)
                        if (extSysnum+i != ex.GetExternalIdxFromProc(sysnum+i,p))
                            throw DROPS::DROPSErrCL("Falsche Zuordnung!");
                }
            }
        }
    }
    return true;
}

/// \brief Create lists for ExchangeCL and do (if wished) time meassurements
void CreateExchangeCL(MultiGridCL& mg, MLIdxDescCL* idx1, MLIdxDescCL* idx2, MLIdxDescCL* xidx)
{
    if (C.tests){
        ParTimerCL time;
        if (ProcCL::IamMaster())
            std::cout << "   - Timemeasurement for creating Exchange lists (and create numbering) for Index 1..."<<'\n';
        // pure
        time.Reset();
        for (int i=0; i<C.tests; ++i){
            time.Start();
            idx1->GetFinest().CreateNumbering(mg.GetLastLevel(), mg);
            time.Stop();
        }
        ::Times.AddTime(Ex_Create, time.GetMaxTime()*100./(double)C.tests);
    }

    if (ProcCL::IamMaster())
        std::cout <<  "   - Kreiere einzelne Exchange-Listen ... " << std::endl;

    if (ProcCL::IamMaster())
        std::cout <<  "     + 1. Exchange-Listen ... " << std::endl;
    idx1->GetFinest().CreateNumbering(mg.GetLastLevel(), mg);
    if (ProcCL::IamMaster())
        std::cout <<  "     + 2. Exchange-Listen ... " << std::endl;
    idx1->GetFinest().CreateNumbering(mg.GetLastLevel(), mg);

    if (C.printEx){
        PrintExchange(idx1->GetEx());
        PrintExchange(idx2->GetEx());
        PrintExchange(xidx->GetEx());
    }

    if (C.printMsgSize)
    {
        if (ProcCL::IamMaster()) std::cout << "  Größe der Nachrichten zwischen den Prozessoren\n  Für Index 1:\n";
        idx1->GetEx().SizeInfo(std::cout);
        if (ProcCL::IamMaster()) std::cout << "  Für Index 2:\n";
        idx2->GetEx().SizeInfo(std::cout);
        if (ProcCL::IamMaster()) std::cout << "  Für XFEM-Index:\n";
        xidx->GetEx().SizeInfo(std::cout);
    }
}

/// \brief Check if ExchangeCL send and accumulate correct
bool TestCorrectnessExchangeCL(MultiGridCL& mg, MLIdxDescCL* idx1, MLIdxDescCL* idx2, MLIdxDescCL* xidx)
{
    const ExchangeCL& ex1= idx1->GetEx();
    const ExchangeCL& ex2= idx2->GetEx();
    const ExchangeCL& exXfem= xidx->GetEx();
    if (ProcCL::IamMaster())
        std::cout << "   - Checke auf Richtigkeit (2 Mal) ... " << '\n'
                  << "     + Setze die GID als Werte auf die Knoten, Edges und Tetras" << '\n'
                  << "     + Tausche die Werte aus" << '\n'
                  << "     + Checke ausgetauschten Werte ..."<<std::endl;

    bool check= true;
    for (int i=0; i<2; ++i)
    {
        VecDescCL y_acc1; y_acc1.SetIdx(idx1);
        VecDescCL y_acc2; y_acc2.SetIdx(idx2);
        VecDescCL y_acc3; y_acc3.SetIdx(xidx);

        SetNum(mg,y_acc1,y_acc2,y_acc3);

        ex1.Accumulate(y_acc1.Data);
        ex2.Accumulate(y_acc2.Data);
        exXfem.Accumulate(y_acc3.Data);

        check = check && CheckNum(mg,y_acc1,y_acc2,y_acc3);

        if (!check)
            throw DROPSErrCL("ExchangeCL has not exchanged values correct");

    }
    if (check && ProcCL::IamMaster())
        std::cout << "    --> OK !" << std::endl;

    return check;
}

/// \brief Test if Sysnums on other procs can be computed correct
bool TestSysnumComputation(MultiGridCL& mg, MLIdxDescCL* idx1, MLIdxDescCL* idx2, MLIdxDescCL* xidx)
{
    const ExchangeCL& ex1= idx1->GetEx();
    const ExchangeCL& ex2= idx2->GetEx();
    const ExchangeCL& exXfem= xidx->GetEx();
    if (ProcCL::IamMaster())
        std::cout << "   - Checke die Berechnung von Sysnums auf anderen Prozessoren:\n";

    bool MappingCheck=true;
    for (int p=0; p<ProcCL::Size(); ++p)
        MappingCheck = MappingCheck &&  ProcCL::Check(CheckIdxMapping(mg, idx1, ex1, p));

    if (!MappingCheck)
        throw DROPSErrCL("Sysnumcomputation for Index 1 is not correct!");

    for (int p=0; p<ProcCL::Size(); ++p)
        MappingCheck = MappingCheck &&  ProcCL::Check(CheckIdxMapping(mg, idx2, ex2, p));

    if (!MappingCheck)
        throw DROPSErrCL("Sysnumcomputation for Index 2 is not correct!");

    for (int p=0; p<ProcCL::Size(); ++p)
        MappingCheck = MappingCheck &&  ProcCL::Check(CheckIdxMapping(mg, xidx, exXfem, p));

    if (!MappingCheck)
        throw DROPSErrCL("Sysnumcomputation for XFEM-Index is not correct!");
    else if (ProcCL::IamMaster())
        std::cout << "    --> OK !" << std::endl;

    return MappingCheck;
}

///\brief Test if ExchangeCL and ExchangeBlockCL are the same and give the same results
bool TestEquality(const ExchangeCL& ex1, const ExchangeCL& ex2, const ExchangeBlockCL& ExBlock)
{
    if (ProcCL::IamMaster())
        std::cout << "   - Checke for equality of the classes ... \n";

    bool isequal1 = ex1.IsEqual(ExBlock.GetEx(0));
    bool isequal2 = ex2.IsEqual(ExBlock.GetEx(1));

    if (ProcCL::Check(isequal1 && isequal2)){
        if (ProcCL::IamMaster())
            std::cout << "    --> OK !" << std::endl;
    }
    else
        throw DROPSErrCL("ExchangeCL differs from ExchangeBlockCL");

    return true;
}

/// \brief Create lists for an ExchangeBlockCL
void CreateExchangeBlockCL(ExchangeBlockCL& ExBlock, MLIdxDescCL* idx1, MLIdxDescCL* idx2)
{
    if (ProcCL::IamMaster())
        std::cout << "   - Kreiere eine Block-Exchange Klasse ... \n";

    ExBlock.AttachTo(idx1->GetFinest());
    ExBlock.AttachTo(idx2->GetFinest());
}

/// \brief Check if ExchangeBlockCL send and accumulate correct
bool TestCorrectnessExchangeBlockCL(MultiGridCL& mg, const ExchangeBlockCL& ExBlock,
                               MLIdxDescCL* idx1, MLIdxDescCL* idx2)
{

    if (ProcCL::IamMaster())
        std::cout << "   - Checke auf Richtigkeit (2 Mal) ... " << '\n'
                  << "     + Setze die GID als Werte auf die Knoten, Edges und Tetras" << '\n'
                  << "     + Tausche die Werte aus" << '\n'
                  << "     + Checke ausgetauschten Werte ..."<<std::endl;

    bool check= true;
    for (int i=0; i<2; ++i)
    {
        // Setzen der Werte
        VecDescCL y_acc1; y_acc1.SetIdx(idx1);
        VecDescCL y_acc2; y_acc2.SetIdx(idx2);
        VecDescCL y_acc3; y_acc3.SetIdx(idx1);
        SetNum(mg,y_acc1,y_acc2,y_acc3);

        // Einen großen Vektor für die ExchangeBlockCL erzeugen
        VectorCL y_big(y_acc1.Data.size() + y_acc2.Data.size());
        y_big[std::slice(0,                  y_acc1.Data.size(), 1)] = y_acc1.Data;
        y_big[std::slice(y_acc1.Data.size(), y_acc2.Data.size(), 1)] = y_acc2.Data;

        // Austausch der Werte
        ExBlock.Accumulate(y_big);
        ExBlock.GetEx(0).Accumulate(y_acc3.Data);

        // Checken der Werte
        y_acc1.Data = y_big[std::slice(0,                  y_acc1.Data.size(), 1)];
        y_acc2.Data = y_big[std::slice(y_acc1.Data.size(), y_acc2.Data.size(), 1)];
        check = check && CheckNum(mg,y_acc1,y_acc2,y_acc3);

        if (!check)
            throw DROPSErrCL("ExchangeCL has not exchanged values correct");
        return true;

    }
    if (check && ProcCL::IamMaster())
        std::cout << "    --> OK !" << std::endl;

    return check;
}

/// \brief Test whether all accumulation strategies and inner products give the same result
void TestInnerProducts(const ExchangeCL& ex1, const ExchangeCL& ex2, const ExchangeBlockCL& exBlock)
{
    const IdxT n1=ex1.GetNum(), n2=ex2.GetNum(), n=exBlock.GetNum();
    // Test vectors
    const VectorCL x1(0.1, n1), x2(0.1, n2), xB(0.1, n);
    // Accumulated vectors
    VectorCL x1_acc_old(x1), x1_acc_new(x1),
             x2_acc_old(x2), x2_acc_new(x2),
             x12_acc_old(n), x12_acc_new(n),
             xB_acc(xB);

    double prod1_old, prod1_new, prod2_old, prod2_new, prod_old, prod_new, prodB;
    double norm1_new, norm2_new, norm_new, normB;
    VectorCL dist(3), global_dist(3);   // distances between accumulated vectors (old-new, old-block, new-block)

    // Compute inner products (on distributed vectors)
    prod1_old= ex1.ParDotAcc(x1_acc_old, x1);
    prod1_new= ex1.ParDot(x1, false, x1, false, true, &x1_acc_new);
    prod2_old= ex2.ParDotAcc(x2_acc_old, x2);
    prod2_new= ex2.ParDot(x2, false, x2, false, true, &x2_acc_new);
    prodB    = exBlock.ParDot(xB, false, xB, false, true, &xB_acc);
    prod_old = prod1_old + prod2_old;
    prod_new = prod1_new + prod2_new;

    // Compute norms (on distributed vectors)
    norm1_new= ex1.Norm_sq(x1, false, true);
    norm2_new= ex2.Norm_sq(x2, false, true);
    norm_new = norm1_new+norm2_new;
    normB    = exBlock.Norm_sq(xB, false, true);

    // Compound small vectors to a big one
    x12_acc_old[std::slice(0,  n1, 1)] = x1_acc_old;
    x12_acc_old[std::slice(n1, n2, 1)] = x2_acc_old;
    x12_acc_new[std::slice(0,  n1, 1)] = x1_acc_new;
    x12_acc_new[std::slice(n1, n2, 1)] = x2_acc_new;

    dist[0]= norm_sq( VectorCL( x1_acc_old-x1_acc_new) + VectorCL( x2_acc_old-x2_acc_new));
    dist[1]= norm_sq( VectorCL( x12_acc_old-xB_acc));
    dist[2]= norm_sq( VectorCL( x12_acc_new-xB_acc));

    ProcCL::GlobalSum(Addr(dist), Addr(global_dist), 3, 0);

    if (ProcCL::IamMaster()){
        std::sqrt(global_dist);
        std::cout << " - Global sum of the local squared distance of the accumulated vectors:\n"
                  << "   * old interface vs. new interface: " << global_dist[0] << "\n"
                  << "   * old interface vs. blocked:       " << global_dist[1] << "\n"
                  << "   * new interface vs. blocked:       " << global_dist[2] << "\n"
                  << std::endl;
        std::cout << " - Differences of inner products (on distributed vectors):\n"
                  << "   * old interface vs. new interface: " << (prod_old-prod_new) << "\n"
                  << "   * old interface vs. blocked:       " << (prod_old-prodB)    << "\n"
                  << "   * new interface vs. blocked:       " << (prod_new-prodB)    << "\n"
                  << std::endl;
        std::cout << " - Differences of norm (on distributed vectors):\n"
                  << "   * new interface vs. blocked:       " << (norm_new-normB)    << "\n"
                  << std::endl;
    }

    // Compute inner products (on one accumulated and one distributed vectors)
    prod1_new= ex1.ParDot(x1, false, x1_acc_new, true, true);
    prod2_new= ex2.ParDot(x2, false, x2_acc_new, true, true);
    prodB    = exBlock.ParDot(xB, false, xB_acc, true, true);
    prod_new = prod1_new + prod2_new;

    // Compute norms (on accumulated vectors)
    norm1_new= ex1.Norm_sq(x1_acc_new, true, true);
    norm2_new= ex2.Norm_sq(x2_acc_new, true, true);
    norm_new = norm1_new+norm2_new;
    normB    = exBlock.Norm_sq(xB_acc, true, true);

    if (ProcCL::IamMaster()){
        std::cout << " - Differences of inner products (on one distributed and one accumulated vector):\n"
                  << "   * new interface vs. blocked:       " << (prod_new-prodB)    << "\n"
                  << std::endl;
        std::cout << " - Differences of norm (on accumulated vectors):\n"
                  << "   * new interface vs. blocked:       " << (norm_new-normB)    << "\n"
                  << std::endl;
    }
}

/// \brief Meassure time used for inner products
void MakeTimeMeassurments( const ExchangeCL& ex, MLIdxDescCL* idx)
{
    ParTimerCL time;
    const IdxT n= idx->NumUnknowns();
    VectorCL x(0.01, n), x_acc(x), y(x), tmp1(n), tmp2(n);
    VecDescCL y_acc(idx);

    if (ProcCL::IamMaster())
        std::cout << "  - Time measurements "<<'\n'
                  << "    + ParDotAcc (old and new interface)"<<'\n'
                  << "    + AccParDot (old and new interface)"<<std::endl;
    time.Reset();
    for (int i=0; i<C.tests; ++i){
        ex.ParDotAcc(x_acc,x);
    }
    time.Stop(); Times.AddTime(Ex_Acc,time.GetMaxTime()*100./(double)C.tests);

    time.Reset();
    for (int i=0; i<C.tests; ++i){
        ex.ParDot(x,false,x, false, false, &x_acc);
    }
    time.Stop(); Times.AddTime(Ex_Acc_newInterface,time.GetMaxTime()*100./(double)C.tests);

    x_acc=x;
    time.Reset();
    for (int i=0; i<C.tests; ++i)
        ex.AccParDot(x,x_acc, tmp1, tmp2);
    time.Stop(); Times.AddTime(Ex_Acc_new,time.GetMaxTime()*100./(double)C.tests);

    time.Reset();
    for (int i=0; i<C.tests; ++i)
        ex.ParDot(x,false, x_acc, false, true, &tmp1, &tmp2);
    time.Stop(); Times.AddTime(Ex_Acc_new_newInterface,time.GetMaxTime()*100./(double)C.tests);

}

double DistanceFct( const DROPS::Point3DCL& p)
{
    return p[0]-2.39;
}

/****************************************************************************
* S T R A T E G Y                                                           *
****************************************************************************/
/// \brief Test ExchangeCL and ExchangeBlockCL
/** This function tests the creating and functionality of ExchangeCL and
    ExchangeBlockCL. If wished, also time measurements are made. We are using
    P1-, P2- and P1X-FEs to test the ExchangeCL.
    /// \todo BlockExchangeCL testen
*/
void Strategy(ParMultiGridCL &pmg)
{
    ParTimerCL time;

    MultiGridCL& mg= pmg.GetMG();
    MLIdxDescCL idx1, idx2;
    VecDescCL x1, x2;

    // Teile dem Index mit, wieviele Unbekannte auf den Verts, Edges und Tetras liegen
    idx1.SetFE( P1_FE);
    idx2.SetFE( P2_FE);

    MLIdxDescCL xfemidx( P1X_FE, 1, BndCondCL(0), 0, 0.1);

    LevelsetP2CL lset( mg);
    lset.CreateNumbering( mg.GetLastLevel(), &lset.idx);
    /// \todo paralleles UpdateXNumbering()
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( DistanceFct);

    if (ProcCL::IamMaster())
        std::cout << " * Numeriere die Unbekannten, teile diese dem PMG mit, Schreibe Infos ... " << std::endl;

    // erzeuge Nummerierung zu diesen Indices
    idx1.CreateNumbering( mg.GetLastLevel(), mg, Bnd);
    idx2.CreateNumbering( mg.GetLastLevel(), mg, Bnd);
    xfemidx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi);

    // Setzte Indices auf die Vector-Describer
//     x1.SetIdx(&idx1);  pmg.AttachTo( 0, &Bnd);
//     x2.SetIdx(&idx2);  pmg.AttachTo( 1, &Bnd);
    VecDescCL xfem( &xfemidx);

    Ulint sizeinfo[3], global_info[3];
    sizeinfo[0] = idx1.NumUnknowns(); sizeinfo[1] = idx2.NumUnknowns(), sizeinfo[2]= xfemidx.NumUnknowns();
    ProcCL::GlobalSum(sizeinfo,global_info, 3, ProcCL::Master());
    if (ProcCL::IamMaster())
        std::cout << "   - Index1: " <<global_info[0]<<", Index2: " <<global_info[1]<<", XfemIdx: " <<global_info[2]<< " Unbekannte (akkumuliert)\n";

    if (ProcCL::IamMaster())
        std::cout << line << std::endl << " * Teste ExchangeCL"<< std::endl;

    ExchangeCL ExXfem;
    CreateExchangeCL(mg, &idx1, &idx2, &xfemidx);
    TestCorrectnessExchangeCL(mg, &idx1, &idx2, &xfemidx);
    TestSysnumComputation(mg, &idx1, &idx2, &xfemidx);

    if (ProcCL::IamMaster())
        std::cout << line << std::endl << " * Teste ExchangeBlockCL" << std::endl;

    ExchangeBlockCL ExBlock;
    CreateExchangeBlockCL(ExBlock, &idx1, &idx2);
    TestCorrectnessExchangeBlockCL(mg, ExBlock, &idx1, &idx2);
    TestEquality( idx1.GetEx(), idx2.GetEx(), ExBlock);

    if (ProcCL::IamMaster())
        std::cout << line << std::endl << " * Checke Paralleles InnereProduct ... " << std::endl;

    TestInnerProducts( idx1.GetEx(), idx2.GetEx(), ExBlock);

    if (C.timeMeas)
    {
        if (ProcCL::IamMaster())
            std::cout << line << std::endl << " * Mache eine Zeitmessungen (" <<C.tests<<" Laeufe) für Innere Produkte ... " << std::endl;
        MakeTimeMeassurments( idx1.GetEx(), &idx1);
    }
}
} // end of namespace DROPS

int main (int argc, char** argv)
{
    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
    try
    {
        SetDescriber();
        //DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE/*|XFER_SHOW_MSGSALL*/);

        if (argc!=2){
            std::cout << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param){
            std::cout << "error while opening parameter file\n"; return 1;
        }
        param >> C;
        param.close();
        if (DROPS::ProcCL::IamMaster())
            std::cout << C << std::endl;

        DROPS::ParTimerCL time, alltime;

        // Initialisierung der parallelen Strukturen.
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();
        DROPS::MultiGridCL   *mg;

        DROPS::Point3DCL orig(0.);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]=C.dx; e2[1]=C.dy; e3[2]= C.dz;

        if (DROPS::ProcCL::IamMaster())
            std::cout << line << std::endl << " * Erstelle das initiale Gitter ... \n";

        if (DROPS::ProcCL::IamMaster())
        {
            DROPS::BrickBuilderCL brick(orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
            mg = new DROPS::MultiGridCL(brick);
        }
        else
        {
            DROPS::EmptyBrickBuilderCL emptyBrick( orig, e1, e2, e3, C.basicref_x);
            mg = new DROPS::MultiGridCL(emptyBrick);
        }

        // Weise dem parallelen MultiGrid das MultiGrid zu
        pmg.AttachTo(*mg);
        DoRefinement(pmg);

        const double balratio=pmg.GetBalance();

        if (DROPS::ProcCL::IamMaster())
            std::cout << "     Der Balancierungsquotient lautet: "<<balratio<<std::endl<< line << std::endl;

        if (C.printMG)
            PrintMG(pmg);

        DROPS::Strategy(pmg);

        if (DROPS::ProcCL::IamMaster())
            std::cout << line << std::endl;

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        Times.Print(cout);
        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}
